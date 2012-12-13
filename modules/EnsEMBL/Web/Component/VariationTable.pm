package EnsEMBL::Web::Component::VariationTable;

use strict;
use Bio::EnsEMBL::Variation::Utils::Constants;

use base qw(EnsEMBL::Web::Component::Variation);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self = shift;
  my $hub = $self->hub;
  my $object_type = $self->hub->type;
  my $consequence_type = $hub->param('sub_table');
  my $table_id = $hub->param('table_id');
  my $icontext = $hub->param('context') || 100;
  my $sum_type = $hub->param('summary_type') || 'tree';
  my $gene_object = $self->configure($icontext, $consequence_type);
  my ($count, $msg, $html);

  my @transcripts      = sort { $a->stable_id cmp $b->stable_id } @{$gene_object->get_all_transcripts};

  if ($object_type eq 'Transcript'){
    my @temp;
    foreach (@transcripts){ 
      push (@temp, $_) if $_->stable_id eq $hub->param('t');
    }
    @transcripts = @temp;
  }
 
  $count += scalar @{$_->__data->{'transformed'}{'gene_snps'}} for @transcripts;

  if ($icontext) {
    if ($icontext eq 'FULL') {
      $msg = "<p>The <b>full</b> intronic sequence around this $object_type is used.";
    } else {
      $msg = "<p>Currently <b>$icontext"."bp</b> of intronic sequence is included either side of the exons.";
    }
    $msg .= qq( To extend or reduce the intronic sequence, use the "<b>Configure this page - Intron Context</b>" link on the left.</p>);

  }
  $msg .= qq(<p>Note: From release 68, Ensembl uses Sequence Ontology (SO) terms to describe consequences. <a href="/info/docs/variation/predicted_data.html#consequence_type_table">More information about this table</a>.</p>);

  if ($consequence_type || $count < 25) {
    $consequence_type ||= 'ALL';

    my $table_rows = $self->variation_table($consequence_type, \@transcripts);
    my $table      = $table_rows ? $self->make_table($table_rows, $consequence_type) : undef;

    $html = $self->render_content($table, $consequence_type, $table_id);
  } else {
    $html  = $self->_hint('snp_table', 'Configuring the page', $msg);
    $html .= $self->render_content($sum_type eq 'tree' ? $self->tree(\@transcripts, $gene_object) : $self->stats_table(\@transcripts, $gene_object)->render); # no sub-table selected, just show stats
  }
  
  return $html;
}

sub make_table {
  my ($self, $table_rows, $consequence_type) = @_;
  my $hub      = $self->hub;
  my $glossary = EnsEMBL::Web::DBSQL::WebsiteAdaptor->new($hub)->fetch_glossary_lookup;

  # Using explicit wdiths speeds things up and makes layout more predictable
  # u = 1unit, where unit is calculated so that total width is 100%
  my $columns = [
    { key => 'ID',       width => '12u', sort => 'html'                                                                                         },
    { key => 'chr' ,     width => '10u', sort => 'position', label => 'Chr: bp'                                                                 },
    { key => 'Alleles',  width => '16u', sort => 'string',   label => "Alle\fles",  align => 'center'                                           },
    { key => 'class',    width => '11u', sort => 'string',   label => 'Class',      align => 'center'                                           },
    { key => "Source",   width => '8u',  sort => 'string', label => "Sour\fce",                                                                                       },
    { key => 'status',   width => '6u',  sort => 'string',   label => "Val\fi\fda\ftion", align => 'center', help => $self->strip_HTML($glossary->{'Validation status'}) },
    { key => 'snptype',  width => '12u', sort => 'string',   label => 'Type',                                                                   },
    { key => 'aachange', width => '6u',  sort => 'string',   label => 'AA',         align => 'center', help => 'Amino Acid'                     },
    { key => 'aacoord',  width => '6u',  sort => 'position', label => "AA co\ford",   align => 'center', help => "Amino Acid Co-ordinate"         },
  ];
  
  # submitter data for LRGs
  splice @$columns, 5, 0, { key => 'Submitters', width => '10u', sort => 'string', align => 'center', export_options => { split_newline => 2 } } if $self->isa('EnsEMBL::Web::Component::LRG::VariationTable');

  # HGVS
  splice @$columns, 3, 0, { key => 'HGVS', width => '10u', sort => 'string', title => 'HGVS name(s)', align => 'center', export_options => { split_newline => 2 } } if $hub->param('hgvs') eq 'on';

  # add GMAF, SIFT and PolyPhen for human
  if ($hub->species eq 'Homo_sapiens') {
    push @$columns, (
      { key => 'sift',     sort => 'position_html', width => '6u', label => "SI\vFT",     align => 'center', help => $self->strip_HTML($glossary->{'SIFT'})     },
      { key => 'polyphen', sort => 'position_html', width => '6u', label => "Poly\fPhen", align => 'center', help => $self->strip_HTML($glossary->{'PolyPhen'}) },
    );

    splice @$columns, 3, 0, { key => 'gmaf', sort => 'numeric', width => '6u', label => "Glo\fbal MAF", align => 'center', help => $self->strip_HTML($glossary->{'Global MAF'}) };
  }
 
  if ($self->hub->type ne 'Transcript'){
   push @$columns, { key => 'Transcript', sort => 'string', width => '11u' };
  }

  return $self->new_table($columns, $table_rows, { data_table => 1, sorting => [ 'chr asc' ], exportable => 1, id => "${consequence_type}_table", class => 'cellwrap_inside fast_fixed_table' });
} 

sub render_content {
  my ($self, $table, $consequence_type, $table_id) = @_;
  my $stable_id = $self->object->stable_id;
  my $html;

  if ($consequence_type) {
    my $consequence_label = ucfirst($self->hub->param('table_title') || $table_id || $consequence_type);
    $consequence_label =~ s/_/ /g;
    $consequence_label =~ s/children/\(with children\)/;
    #$consequence_label .='s';
    $html = $self->toggleable_table("$consequence_label consequences", $table_id, $table, 1, qq{<span style="float:right"><a href="#$self->{'id'}_top">[back to top]</a></span>});
  } else {
    my $hub = $self->hub;
    my $current = $hub->param('summary_type') || 'tree';
    my $switched = $current eq 'tree' ? 'table' : 'tree';
    my $url = $hub->url({summary_type => $switched});
    
    $html = qq{
      <a id="$self->{'id'}_top"></a>
      <span style="float:right;">
        <a href="$url">Switch to $switched view <img src="/i/16/reload.png" height="12px"/></a>
        
      </span>
      <h2>Summary of variation consequences in $stable_id</h2>
    } . $table;
  }
  
  return $html;
}


sub stats_table {
  my ($self, $transcripts, $gene_object) = @_;
  
  my $hub         = $self->hub;
  
  my $columns = [
    { key => 'count', title => 'Number of variants', sort => 'numeric_hidden', width => '20%', align => 'right'  },   
    { key => 'view',  title => '',                   sort => 'none',           width => '5%',  align => 'center' },   
    { key => 'key',   title => '',                   sort => 'none',           width => '2%',  align => 'center' },
    { key => 'type',  title => 'Type',               sort => 'numeric_hidden', width => '20%'                    },   
    { key => 'desc',  title => 'Description',        sort => 'none',           width => '53%'                    },
  ];
  
  my (%counts, %total_counts, %ranks, %descriptions, %labels, %colours);
  my $total_counts;
  
  # colour stuff
  my $species_defs = $hub->species_defs;
  my $var_styles   = $species_defs->colour('variation');
  my $colourmap    = $hub->colourmap;
  
  my @all_cons = grep $_->feature_class =~ /transcript/i, values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  
  foreach my $con(@all_cons) {
    next if $con->SO_accession =~ /x/i;
    
    my $term = $con->SO_term;
    
    $labels{$term}       = $con->label;
    $descriptions{$term} = $con->description.' <span class="small">('.$hub->get_ExtURL_link($con->SO_accession, 'SEQUENCE_ONTOLOGY', $con->SO_accession).')</span' unless $descriptions{$term};
    
    $colours{$term} = $colourmap->hex_by_name($var_styles->{lc($con->SO_term)}->{default});
    print STDERR $con->SO_term."\t".$colours{$term}."\n";
    $ranks{$term} = $con->rank if $con->rank < $ranks{$term} || !defined($ranks{$term});
  }

  if (!exists($gene_object->__data->{'conscounts'})) {
# Generate the data the hard way - from all the vfs and tvs
    my %counts_hash;
    my %total_counts_hash;

    foreach my $tr (@$transcripts) { 
      my $tr_stable_id = $tr->stable_id;
      my $tvs          = $tr->__data->{'transformed'}{'snps'} || {};
      my $gene_snps    = $tr->__data->{'transformed'}{'gene_snps'};
      my $tr_start     = $tr->__data->{'transformed'}{'start'};
      my $tr_end       = $tr->__data->{'transformed'}{'end'};
      my $extent       = $tr->__data->{'transformed'}{'extent'};
      
      foreach (@$gene_snps) {
        my ($snp, $chr, $start, $end) = @$_;
        my $vf_id = $snp->dbID;
        my $tv    = $tvs->{$vf_id};
        
        if ($tv && $end >= $tr_start - $extent && $start <= $tr_end + $extent) {
          foreach my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
            foreach my $con (@{$tva->get_all_OverlapConsequences}) {
              my $key  = join '_', $tr_stable_id, $vf_id, $tva->variation_feature_seq;
              my $term = $con->SO_term;
              
              $counts_hash{$term}{$key} = 1 if $con;
              $total_counts_hash{$key}++;
            }
          }
        }
      }
    }
    
    foreach my $con (keys %descriptions) {
      $counts{$con} = scalar keys %{$counts_hash{$con}};
    }
    $total_counts = scalar keys %total_counts_hash;

  } else {
# Use the results of the TV count queries
    %counts = %{$gene_object->__data->{'conscounts'}};

    $total_counts = $counts{'ALL'};
  }
  
  my $warning_text = qq{<span style="color:red;">(WARNING: table may not load for this number of variants!)</span>};
  my @rows;
  
  foreach my $con (keys %descriptions) {
    my $colour_block = sprintf('<div style="background-color: %s; width: 10px;">&nbsp;</div>', $colours{$con});
    
    if ($counts{$con}) {
      my $count = $counts{$con};
      my $warning = $count > 10000 ? $warning_text : '';
      
      push @rows, {
        type  => qq{<span class="hidden">$ranks{$con}</span>$labels{$con}},
        desc  => $descriptions{$con}.' '.$warning,
        count => $count,
        view  => $self->ajax_add($self->ajax_url(undef, { sub_table => $con, update_panel => 1 }), $con),
        key   => $colour_block,
      };
    } else {
      push @rows, {
        type  => qq{<span class="hidden">$ranks{$con}</span>$labels{$con}},
        desc  => $descriptions{$con},
        count => 0,
        view  => '-',
        key   => $colour_block,
      };
    }
  }
  
  # add the row for ALL variations if there are any
  if (my $total = $total_counts) {
    my $hidden_span = qq{<span class="hidden">-</span>}; # create a hidden span to add so that ALL is always last in the table
    my $warning     = $total > 10000 ? $warning_text : '';
    
    push @rows, {
      type  => $hidden_span . 'ALL',
      view  => $self->ajax_add($self->ajax_url(undef, { sub_table => 'ALL', update_panel => 1 }), 'ALL'),
      desc  => "All variations $warning",
      count => $hidden_span . $total,
    };
  }

  return $self->new_table($columns, \@rows, { data_table => 'no_col_toggle', sorting => [ 'type asc' ], exportable => 0 });
}

sub tree {
  my ($self, $transcripts, $gene_object) = @_;
  
  my $hub         = $self->hub;
  
  # define top-level SO term
  my $top_SO_term = 'feature_variant';
  
  # get counts
  my %counts;
  
  if (!exists($gene_object->__data->{'conscounts'})) {
    # Generate the data the hard way - from all the vfs and tvs
    my %counts_hash;
    my %total_counts_hash;

    foreach my $tr (@$transcripts) { 
      my $tr_stable_id = $tr->stable_id;
      my $tvs          = $tr->__data->{'transformed'}{'snps'} || {};
      my $gene_snps    = $tr->__data->{'transformed'}{'gene_snps'};
      my $tr_start     = $tr->__data->{'transformed'}{'start'};
      my $tr_end       = $tr->__data->{'transformed'}{'end'};
      my $extent       = $tr->__data->{'transformed'}{'extent'};
      
      foreach (@$gene_snps) {
        my ($snp, $chr, $start, $end) = @$_;
        my $vf_id = $snp->dbID;
        my $tv    = $tvs->{$vf_id};
        
        if ($tv && $end >= $tr_start - $extent && $start <= $tr_end + $extent) {
          foreach my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
            foreach my $con (@{$tva->get_all_OverlapConsequences}) {
              my $key  = join '_', $tr_stable_id, $vf_id, $tva->variation_feature_seq;
              my $term = $con->SO_term;
              
              $counts_hash{$term}{$key} = 1 if $con;
            }
          }
        }
      }
    }
    
    foreach my $con (keys %counts_hash) {
      $counts{$con} = scalar keys %{$counts_hash{$con}};
    }
  } else {
    # Use the results of the TV count queries
    %counts = %{$gene_object->__data->{'conscounts'}};
  }
  
  # get SO tree
  my $tree = $self->get_SO_tree($top_SO_term);
  
  # add counts to tree
  $self->add_counts_to_tree($tree, \%counts);
  
  # add colors
  my $species_defs = $hub->species_defs;
  my $var_styles   = $species_defs->colour('variation');
  my $colourmap    = $hub->colourmap;
  
  $self->add_colours_to_tree($tree, $var_styles, $colourmap);
  
  my $html = '<ul class="tree">';
  
  $html .= $self->tree_html($tree, 1);
  
  $html .= '</ul>';
  
  return $html;
}


sub variation_table {
  my ($self, $consequence_type, $transcripts, $slice) = @_;
  my $hub         = $self->hub;
  my $show_scores = $hub->param('show_scores');
  my (@rows, $base_trans_url, $url_transcript_prefix, %handles);
  
  # create some URLs - quicker than calling the url method for every variation
  my $base_url = $hub->url({
    type   => 'Variation',
    action => 'Mappings',
    vf     => undef,
    v      => undef,
    source => undef,
  });
  
  if ($self->isa('EnsEMBL::Web::Component::LRG::VariationTable')) {
    my $gene_stable_id     = $transcripts->[0] && $transcripts->[0]->gene ? $transcripts->[0]->gene->stable_id : undef;
    $url_transcript_prefix = 'lrgt';
    
    $base_trans_url = $hub->url({
      type    => 'LRG',
      action  => 'Summary',
      lrg     => $gene_stable_id,
      __clear => 1
    });
    
    my $vfa = $hub->get_adaptor('get_VariationFeatureAdaptor', 'variation');
    
    my @var_ids =
      map {$_->{_variation_id}}
      map {$_->[0]}
      map {@{$_->__data->{transformed}{gene_snps}}}
      @$transcripts;
    
    %handles = %{$vfa->_get_all_subsnp_handles_from_variation_ids(\@var_ids)};
    
  } else {
    $url_transcript_prefix = 't';
    
    $base_trans_url = $hub->url({
      type   => 'Transcript',
      action => 'Summary',
      t      => undef,
    });
  }
  
  foreach my $transcript (@$transcripts) {
    my %snps = %{$transcript->__data->{'transformed'}{'snps'} || {}};
   
    next unless %snps;
    
    my $transcript_stable_id = $transcript->stable_id;
    my $gene_snps            = $transcript->__data->{'transformed'}{'gene_snps'} || [];
    my $tr_start             = $transcript->__data->{'transformed'}{'start'};
    my $tr_end               = $transcript->__data->{'transformed'}{'end'};
    my $extent               = $transcript->__data->{'transformed'}{'extent'};
    my $gene                 = $transcript->gene;

    foreach (@$gene_snps) {
      my ($snp, $chr, $start, $end) = @$_;
      my $raw_id               = $snp->dbID;
      my $transcript_variation = $snps{$raw_id};
      
      next unless $transcript_variation;
      
      foreach my $tva (@{$transcript_variation->get_all_alternate_TranscriptVariationAlleles}) {
        my $skip = 1;
        
        if ($consequence_type eq 'ALL') {
          $skip = 0;
        } elsif ($tva) {
          foreach my $con (map {$_->SO_term} @{$tva->get_all_OverlapConsequences}) {
            if (grep {$con eq $_} split(/\,/, $consequence_type)) {
              $skip = 0;
              last;
            }
          }
        }
        
        next if $skip;
        
        if ($tva && $end >= $tr_start - $extent && $start <= $tr_end + $extent) {
          #my $var                  = $snp->variation;
          my $validation           = $snp->get_all_validation_states || [];
          my $variation_name       = $snp->variation_name;
          my $var_class            = $snp->var_class;
          my $translation_start    = $transcript_variation->translation_start;
          my $source               = $snp->source;
          my ($aachange, $aacoord) = $translation_start ? ($tva->pep_allele_string, $translation_start) : ('-', '-');
          my $url                  = "$base_url;v=$variation_name;vf=$raw_id;source=$source";
          my $trans_url            = "$base_trans_url;$url_transcript_prefix=$transcript_stable_id";
          my $vf_allele            = $tva->variation_feature_seq;
          my $allele_string        = $snp->allele_string;
             $allele_string        = $self->trim_large_allele_string($allele_string, 'allele_' . $tva->dbID, 20) if length $allele_string > 20; # Check allele string size (for display issues)
             $allele_string        =~ s/$vf_allele/<b>$vf_allele<\/b>/g if $allele_string =~ /\/.+\//; # highlight variant allele in allele string
          
          # sort out consequence type string
          # Avoid duplicated Ensembl terms
          my %term   = map { $hub->get_ExtURL_link($_->label, 'SEQUENCE_ONTOLOGY', $_->SO_accession) => 1 } @{$tva->get_all_OverlapConsequences || []};
          my $type   = join ',<br />', keys %term;
             $type ||= '-';
          
          my $sift = $self->render_sift_polyphen($tva->sift_prediction,     $tva->sift_score);
          my $poly = $self->render_sift_polyphen($tva->polyphen_prediction, $tva->polyphen_score);
          
          # Adds LSDB/LRG sources
          if ($self->isa('EnsEMBL::Web::Component::LRG::VariationTable')) {
            my $var = $snp->variation;
            
            my $syn_sources = $var->get_all_synonym_sources;
            
            foreach my $s_source (@$syn_sources) {
              next if $s_source !~ /LSDB|LRG/;
              
              my ($synonym) = $var->get_all_synonyms($s_source);
                 $source   .= ', ' . $hub->get_ExtURL_link($s_source, $s_source, $synonym);
            }
          }
          
          my $gmaf   = $snp->minor_allele_frequency; # global maf
             $gmaf   = sprintf '%.3f <span class="small">(%s)</span>', $gmaf, $snp->minor_allele if defined $gmaf;

          my $status = join '', map {qq{<img height="20px" width="20px" title="$_" src="/i/96/val_$_.gif"/><span class="hidden export">$_,</span>}} @$validation; # validation
          
          my $row = {
            ID         => qq{<a href="$url">$variation_name</a>},
            class      => $var_class,
            Alleles    => $allele_string,
            Ambiguity  => $snp->ambig_code,
            gmaf       => $gmaf   || '-',
            status     => $status || '-',
            chr        => "$chr:$start" . ($start == $end ? '' : "-$end"),
            Source     => $source,
            Submitters => %handles && defined($handles{$snp->{_variation_id}}) ? join(", ", @{$handles{$snp->{_variation_id}}}) : undef,
            snptype    => $type,
            Transcript => qq{<a href="$trans_url">$transcript_stable_id</a>},
            aachange   => $aachange,
            aacoord    => $aacoord,
            sift       => $sift,
            polyphen   => $poly,
            HGVS       => $hub->param('hgvs') eq 'on' ? ($self->get_hgvs($tva) || '-') : undef,
          };
          
          push @rows, $row;
        }
      }
    }
  }

  return \@rows;
}

sub create_so_term_subsets {
  my ($self) = @_;

  my @all_cons = grep $_->feature_class =~ /Bio::EnsEMBL::(Feature|Transcript)/i, values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

  my %so_term_subsets;
  
  foreach my $con (@all_cons) {
    next if $con->SO_accession =~ /x/i;
    
    my $term = $con->SO_term;
    
    if (!exists($so_term_subsets{$term})) {
      $so_term_subsets{$term} = [];
    }
    push @{$so_term_subsets{$term}}, $con->SO_term;
  }

  return \%so_term_subsets;
}

sub configure {
  my ($self, $context, $consequence) = @_;
  my $object = $self->object;
  my $object_type = $self->hub->type;
  my $extent = $context eq 'FULL' ? 5000 : $context;
  my $gene_object;
  my $transcript_object;

  if ($object->isa('EnsEMBL::Web::Object::Gene')){ #|| $object->isa('EnsEMBL::Web::Object::LRG')){
    $gene_object = $object;
  } elsif ($object->isa('EnsEMBL::Web::Object::LRG')){
    my @genes       = @{$object->Obj->get_all_Genes('LRG_import')||[]};
    my $gene = $genes[0];  
    my $factory = $self->builder->create_factory('Gene');  
    $factory->createObjects($gene);
    $gene_object = $factory->object;
  } else {
    $transcript_object = $object;
    $gene_object = $self->hub->core_objects->{'gene'};
  }
  
  $gene_object->get_gene_slices(
    undef,
    [ 'context',     'normal', '100%'  ],
    [ 'gene',        'normal', '33%'   ],
    [ 'transcripts', 'munged', $extent ]
  );

  # map the selected consequence type to SO terms
  my %cons = %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  my %selected_so;
  
  $selected_so{$_} = 1 for split /\,/, $consequence;
  
  my @so_terms = keys %selected_so;

  $gene_object->store_TransformedTranscripts;      ## Stores in $transcript_object->__data->{'transformed'}{'exons'|'coding_start'|'coding_end'}
 
  my $transcript_slice = $gene_object->__data->{'slices'}{'transcripts'}[1];
  my (undef, $snps)    = $gene_object->getVariationsOnSlice($transcript_slice, $gene_object->__data->{'slices'}{'transcripts'}[2], undef,  scalar(@so_terms) ? \@so_terms : undef);


  my $vf_objs = [ map $_->[2], @$snps];

  # For stats table (no $consquence) without a set intron context ($context). Also don't try for a single transcript (because its slower)
  if (!$consequence && !$transcript_object && $context eq 'FULL') {
    my $so_term_subsets = $self->create_so_term_subsets;
    $gene_object->store_ConsequenceCounts($so_term_subsets, $vf_objs);
  }

  # If doing subtable or can't calculate consequence counts
  if ($consequence || !exists($gene_object->__data->{'conscounts'})) {
    $gene_object->store_TransformedSNPS(\@so_terms,$vf_objs); ## Stores in $transcript_object->__data->{'transformed'}{'snps'}
  }

  ## Map SNPs for the last SNP display  
  my @gene_snps = map {[
    $_->[2], $transcript_slice->seq_region_name,
    $transcript_slice->strand > 0 ?
      ( $transcript_slice->start + $_->[2]->start - 1, $transcript_slice->start + $_->[2]->end   - 1 ) :
      ( $transcript_slice->end   - $_->[2]->end   + 1, $transcript_slice->end   - $_->[2]->start + 1 )
  ]} @$snps;

  foreach (@{$gene_object->get_all_transcripts}) {
    next if $object_type eq 'Transcript' && $_->stable_id ne $self->hub->param('t'); 
    $_->__data->{'transformed'}{'extent'}    = $extent;
    $_->__data->{'transformed'}{'gene_snps'} = \@gene_snps;
  }

  return $gene_object;
}

sub get_hgvs {
  my ($self, $tva) = @_;
  my $hgvs_c = $tva->hgvs_coding;
  my $hgvs_p = $tva->hgvs_protein;
  my $hgvs;

  if ($hgvs_c) {
    if (length $hgvs_c > 35) {
      my $display_hgvs_c  = substr($hgvs_c, 0, 35) . '...';
         $display_hgvs_c .= $self->trim_large_string($hgvs_c, 'hgvs_c_' . $tva->dbID);

      $hgvs_c = $display_hgvs_c;
    }

    $hgvs .= $hgvs_c;
  }

  if ($hgvs_p) {
    if (length $hgvs_p > 35) {
      my $display_hgvs_p  = substr($hgvs_p, 0, 35) . '...';
         $display_hgvs_p .= $self->trim_large_string($hgvs_p, 'hgvs_p_'. $tva->dbID);

      $hgvs_p = $display_hgvs_p;
    }

    $hgvs .= "<br />$hgvs_p";
  }

  return $hgvs;
}

sub get_SO_tree {
  my ($self, $top_SO_term) = @_;
  
  my $oa = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;
  
  my ($top_SO_obj) = @{$oa->fetch_all_by_name($top_SO_term)};
  
  return $top_SO_obj;
}

sub add_counts_to_tree {
  my ($self, $term_obj, $counts) = @_;
  
  $self->add_counts_to_tree($_, $counts) for @{$term_obj->children};
  
  my $count = 0;
  
  $term_obj->{this_count} = defined($counts->{$term_obj->name}) ? $counts->{$term_obj->name} : 0;
  $count += $term_obj->{this_count};
  $count += (defined $_->{count} ? $_->{count} : 0) for @{$term_obj->children};
  
  push @{$term_obj->{term_list}}, ($term_obj->name, map {@{$_->{term_list}}} @{$term_obj->children});
  
  $term_obj->{count} = $count;
}

sub add_colours_to_tree {
  my ($self, $term_obj, $var_styles, $colourmap) = @_;
  $term_obj->{colour} = $colourmap->hex_by_name($var_styles->{lc($term_obj->name)}->{default}) if defined $var_styles->{lc($term_obj->name)};
  $self->add_colours_to_tree($_, $var_styles, $colourmap) for @{$term_obj->children};
}

sub tree_html {
  my ($self, $term_obj, $last) = @_;
  
  my $con = $term_obj->name;
  $con =~ s/\_/ /g;
  $con = "\u$con";
  
  $con = 'All variants' if $con eq 'Feature variant';
  
  my %include_cons = map {$_->SO_term => 1} values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  
  my @children = grep {$_->{count}} @{$term_obj->children};
 
  # don't go further in these cases only
  #undef @children if $term_obj->name =~ /feature.+ation/i;
  my $side_padding = @children ? '7px' : '5px';
  my $html = sprintf(
    '<li%s>%s<span title="%s">%s%s',
    
    # last css class
    ($last ? ' class="last"' : (@children ? ' class="parent top_level"' : '')),
    
    # toggle bit
    (@children ? '<a href="#" class="toggle open" rel="'.$term_obj->name.'"/>' : '<img src="/i/leaf.gif">'),
    
    # consequence definition etc
    (split /\"/, $term_obj->definition)[1],
    
    # colour block
    
    (defined($term_obj->{colour}) ? sprintf('<div style="background-color: %s; color: %s; width: 10px; display: inline; margin: 0 3px 0 0; padding: 0 %s;"> </div>', $term_obj->{colour}, $term_obj->{colour}, $side_padding) : ''),
    
    # name and link
    $con,
  );
  
  # this term only
  if($term_obj->{this_count} || $term_obj->{count}) {
    
    my $sub_table = @children ? join(",", grep {defined($include_cons{$_})} @{$term_obj->{term_list}}) : $term_obj->name;
    
    my $link = $self->ajax_add(
      $self->ajax_url(
        undef,
        {
          sub_table => $sub_table,
          update_panel => 1,
          table_title => $term_obj->name,
        }
      ),
      $term_obj->name,
    );
    
    $html .= sprintf(
      ' | %s (%i)',
      $link,
      $term_obj->{count},
    );
    
  }
  
  $html .= '<span class="small"> | '.$self->hub->get_ExtURL_link($term_obj->accession, 'SEQUENCE_ONTOLOGY', $term_obj->accession).'</span>' unless $con eq 'All variants';
  
  my $warning_text = $term_obj->{count} > 10000 ? qq{<span style="color:red;">(WARNING: table may not load for this number of variants!)</span>} : '';

  $html .= "</span> $warning_text";
 
  
  # iterate  
  if(@children) {
    
    @children = sort {
      ($b->name =~ /stream/) <=> ($a->name =~ /stream/) ||
      scalar @{$a->children} <=> scalar @{$b->children}
    } @children;
    
    $html .= '<div class="'.$term_obj->name.'" ><ul style="line-height: 24px;" class="toggleable">';
    
    for(my $i=0; $i<@children; $i++) {
      $html .= $self->tree_html($children[$i], ($i + 1 == @children ? 1 : 0));
    }
    
    $html .= '</ul></div>';
  }
  
  $html .= '</li>';
  
  return $html;
}


1;
