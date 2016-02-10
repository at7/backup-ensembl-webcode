=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EnsEMBL::Web::Component::Variation::PairwiseLD;

use strict;
use HTML::Entities qw(encode_entities);
use POSIX qw(floor);
use base qw(EnsEMBL::Web::Component::Variation);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self = shift;
  my $object = $self->object;
  my $variant = $object->Obj;
  my $variant_name = $variant->name;

  my $hub    = $self->hub;
  
  return $self->_info('A unique location can not be determined for this Variation', $object->not_unique_location) if $object->not_unique_location;  

  my $url = $self->ajax_url('results', { focus_variant_name => $variant_name, second_variant_name => undef });  
  my $id  = $self->id;  
  my $second_variant_name = '';
  return sprintf('
    <h2>Pairwise linkage disequilibrium data by population</h2>
    <div class="navbar print_hide" style="padding-left:5px">
      <input type="hidden" class="panel_type" value="Content" />
      <form class="update_panel" action="#">
        <label for="variant">Focus variant: %s</label><br>
        <label for="variant">Enter the name for the second variant:</label>
        <input type="text" name="second_variant_name" id="variant" value="%s" size="30"/>
        <input type="hidden" name="focus_variant_name" value="%s" />
        <input type="hidden" name="panel_id" value="%s" />
        <input type="hidden" name="url" value="%s" />
        <input type="hidden" name="element" value=".results" />
        <input class="fbutton" type="submit" value="Compute" />
        <small>(e.g. rs678)</small>
      </form>
    </div>
    <div class="results">%s</div>
  ', $variant_name, $second_variant_name, $variant_name, $id, $url, $self->content_results);

}


sub format_parent {
  my ($self, $parent_data) = @_;
  return ($parent_data && $parent_data->{'Name'}) ? $parent_data->{'Name'} : '-';
}


sub get_table_headings {
  return [
    { key => 'Population', title => 'Population', sort => 'html', align => 'left' },    
    { key => 'Description', title => 'Description', sort => 'string', align => 'left' },
    { key => 'Variant1', title => 'Focus Variant', sort => 'string' },
    { key => 'Variant2', title => 'Variant 2', sort => 'string' },
    { key => 'r2', title => 'r<sup>2</sup>', sort => 'numeric', align => 'center' },
    { key => 'd_prime', title => q{D'}, sort => 'numeric', align => 'center' },
  ];
}

sub content_results {
  my $self         = shift;
  my $object       = $self->object;
  my $variant      = $object->Obj;
  my $hub          = $self->hub;
  my $species_defs = $hub->species_defs;


  my @colour_gradient = ('ffffff', $hub->colourmap->build_linear_gradient(41, 'mistyrose', 'pink', 'indianred2', 'red'));

  my $focus_variant_name  = $variant->name;
  my $second_variant_name = $hub->param('second_variant_name');

  return unless $second_variant_name;
  $second_variant_name =~ s/^\W+//;
  $second_variant_name =~ s/\s+$//;

  # set path information for LD calculations
  $Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor::BINARY_FILE = $species_defs->ENSEMBL_CALC_GENOTYPES_FILE;
  $Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor::TMP_PATH = $species_defs->ENSEMBL_TMP_TMP;

  my $ldfca = $variant->adaptor->db->get_LDFeatureContainerAdaptor;
  my $va = $variant->adaptor->db->get_VariationAdaptor;
  my $pa = $variant->adaptor->db->get_PopulationAdaptor;

  my $second_variant = $va->fetch_by_name($second_variant_name);

  if (!$second_variant) {
    return qq{<div>Could not fetch variant object for variant $second_variant_name</div>};
  }

  my $source = $variant->source_name;
  my $max_distance = $hub->param('max_distance') || 50000;
  my $min_r2 = defined $hub->param('min_r2') ? $hub->param('min_r2') : 0.8;
  my $min_d_prime = defined $hub->param('min_d_prime') ? $hub->param('min_d_prime') : 0.8;
  my $min_p_log = $hub->param('min_p_log');
  my $only_phenotypes = $hub->param('only_phenotypes') eq 'yes';
  my %mappings = %{$object->variation_feature_mapping}; # determine correct SNP location
  my ($vf, $loc);

  if (keys %mappings == 1) {
    ($loc) = values %mappings; 
  } else {
    $loc = $mappings{$hub->param('vf')};
  }
  # get the VF that matches the selected location
  foreach (@{$object->get_variation_features}) {
    if ($_->seq_region_start == $loc->{'start'}) {
      $vf = $_;
      last;
    }
  }

  my $vfs2 = $second_variant->get_all_VariationFeatures;
  my @vfs = ($vf, @$vfs2);
  my @ld_populations = @{$pa->fetch_all_LD_Populations};
  my $rows = [];

  foreach my $ld_population (@ld_populations) {
    my $description = $ld_population->description;
    $description ||= '-';

    if (length $description > 30) {
      my $full_desc = $self->strip_HTML($description);
      while ($description =~ m/^.{30}.*?(\s|\,|\.)/g) {
        $description = sprintf '%s... <span class="_ht ht small" title="%s">(more)</span>', substr($description, 0, (pos $description) - 1), $full_desc;
        last;
      }
    }
    my $pop_name  = $ld_population->name;
    my $pop_dbSNP = $ld_population->get_all_synonyms('dbSNP');

    my $pop_label = $pop_name;
    if ($pop_label =~ /^.+\:.+$/ and $pop_label !~ /(http|https):/) {
      my @composed_name = split(':', $pop_label);
      $composed_name[$#composed_name] = '<b>'.$composed_name[$#composed_name].'</b>';
      $pop_label = join(':',@composed_name);
    }

    # Population external links
    my $pop_url;
    if ($pop_name =~ /^1000GENOMES/) {
      $pop_url = $self->hub->get_ExtURL_link($pop_label, '1KG_POP', $pop_name);
    }
    else {
      $pop_url = $pop_dbSNP ? $self->hub->get_ExtURL_link($pop_label, 'DBSNPPOP', $pop_dbSNP->[0]) : $pop_label;
    }

    my @ld_values = @{$ldfca->fetch_by_VariationFeatures(\@vfs, $ld_population)->get_all_ld_values(1)};  
    foreach my $hash (@ld_values) {
      my $variation1 = $hash->{variation_name1};
      my $variation2 = $hash->{variation_name2};
      next unless ($variation1 eq $focus_variant_name || $variation2 eq $focus_variant_name);
      next unless ($variation1 eq $second_variant_name || $variation2 eq $second_variant_name);
      if ($variation1 ne $focus_variant_name) {
        ($variation1, $variation2) = ($variation2, $variation1);
      }
      my $r2 = $hash->{r2};
      my $d_prime = $hash->{d_prime};
      my $population_id = $hash->{population_id};

      my $row = {
        Variant1 => $variation1, 
        Variant2 => $variation2,
        Population => $pop_url,
        Description => $description, 
        r2 => {
          value => $r2,
          style => "background-color:#".($r2 eq '-' ? 'ffffff' : $colour_gradient[floor($r2*40)]),
        },
        d_prime => {
          value => $d_prime,
          style => "background-color:#".($d_prime eq '-' ? 'ffffff' : $colour_gradient[floor($d_prime*40)]),
        },
      };
      push @$rows, $row;
    } 
  }
 
  my $columns = $self->get_table_headings;
  my $table = $self->new_table($columns, $rows, { data_table => 1, download_table => 1, sorting => [ 'd_prime dec' ] });
  my $html = '<div style="margin:0px 0px 50px;padding:0px">'.$table->render.'</div>';

 return qq{<div class="js_panel">$html</div>};

}


sub get_all_populations {
  my $self   = shift;
  my $sample = shift;

  my @pop_names = map { $_->name } @{$sample->get_all_Populations };
  
  return (scalar @pop_names > 0) ? join('; ',sort(@pop_names)) : '-';
}

sub warning_message {
  my $self   = shift;
  my $sample = shift;
  
  return $self->_warning('Not found', qq{No genotype associated with this variant was found for the sample name '<b>$sample</b>'!});
}



1;
