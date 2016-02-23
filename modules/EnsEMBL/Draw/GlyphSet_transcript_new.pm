=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package EnsEMBL::Draw::GlyphSet_transcript_new;

### Parent module for various glyphsets that draw transcripts
### (styles include exons as blocks, joined by angled lines across introns)

use strict;

use List::Util qw(min max);
use List::MoreUtils qw(natatime);

use base qw(EnsEMBL::Draw::GlyphSet_transcript_new_base);

sub features { return []; }

## Let us define all the renderers here...
## ... these are just all wrappers - the parameter is 1 to draw labels
## ... 0 otherwise...

sub render_normal                  { $_[0]->render_transcripts(1);           }
sub render_transcript              { $_[0]->render_transcripts(1);           }
sub render_transcript_label        { $_[0]->render_transcripts(1);           }
sub render_transcript_label_coding { $_[0]->render_transcripts(1);           }
sub render_transcript_gencode_basic{ $_[0]->render_transcripts(1);           }
sub render_transcript_nolabel      { $_[0]->render_transcripts(0);           }
sub render_collapsed_label         { $_[0]->render_collapsed(1);             }
sub render_collapsed_nolabel       { $_[0]->render_collapsed(0);             }
sub render_gene_label              { $_[0]->render_genes(1);                 }
sub render_gene_nolabel            { $_[0]->render_genes(0);                 }
sub render_as_transcript_label     { $_[0]->render_alignslice_transcript(1); }
sub render_as_transcript_nolabel   { $_[0]->render_alignslice_transcript(0); }
sub render_as_collapsed_label      { $_[0]->render_alignslice_collapsed(1);  }
sub render_as_collapsed_nolabel    { $_[0]->render_alignslice_collapsed(0);  }

sub calculate_collapsed_joins {
  my ($self,$gene_stable_id) = @_;
 
  my $hub = $self->{'config'}->hub;
  my $ga = $hub->get_adaptor('get_GeneAdaptor',
                             $self->my_config('db'),$self->species);
  my $gene = $ga->fetch_by_stable_id($gene_stable_id);
 
  my $previous_species = $self->my_config('previous_species');
  my $next_species     = $self->my_config('next_species');
  my $previous_target  = $self->my_config('previous_target');
  my $next_target      = $self->my_config('next_target');
  my $join_types       = $self->get_parameter('join_types');
  my $alt_alleles     = $gene->get_all_alt_alleles;
  my $seq_region_name = $gene->slice->seq_region_name;
  my ($target, @gene_tags);
  
  my @joins;

  if ($previous_species) {
    for ($self->get_gene_joins($gene, $previous_species, $join_types)) {
      $target = $previous_target ? ":$seq_region_name:$previous_target" : '';
      push @joins,{
        key => "$gene_stable_id:$_->[0]$target",
        colour => $_->[1],
        legend => $_->[2]
      };          
    }
    
    push @gene_tags, map { join '=', $_->stable_id, $gene_stable_id } @{$self->filter_by_target($alt_alleles, $previous_target)};
  }

  if ($next_species) {
    for ($self->get_gene_joins($gene, $next_species, $join_types)) {
      $target = $next_target ? ":$next_target:$seq_region_name" : '';
      push @joins,{
        key => "$_->[0]:$gene_stable_id$target",
        colour => $_->[1],
        legend => $_->[2]
      };
    }
    
    push @gene_tags, map { join '=', $gene_stable_id, $_->stable_id } @{$self->filter_by_target($alt_alleles, $next_target)};
  }
  my $alt_alleles_col  = $self->my_colour('alt_alleles_join');
  for (@gene_tags) {
    push @joins,{
      key => $_,
      colour => $alt_alleles_col,
      legend => 'Alternative alleles'
    };
  }
  return \@joins;
}

sub calculate_expanded_joins {
  my ($self,$gene_stable_id) = @_;

  my $hub = $self->{'config'}->hub;
  my $ga = $hub->get_adaptor('get_GeneAdaptor',
                             $self->my_config('db'),$self->species);
  my $gene = $ga->fetch_by_stable_id($gene_stable_id);

  my $previous_species = $self->my_config('previous_species');
  my $next_species     = $self->my_config('next_species');
  my $previous_target  = $self->my_config('previous_target');
  my $next_target      = $self->my_config('next_target');
  my $join_types       = $self->get_parameter('join_types');
  my $seq_region_name = $gene->slice->seq_region_name;
  my $alt_alleles = $gene->get_all_alt_alleles;
  my $alltrans    = $gene->get_all_Transcripts; # vega stuff to link alt-alleles on longest transcript
  my @s_alltrans  = sort { $a->length <=> $b->length } @$alltrans;
  my $long_trans  = pop @s_alltrans;
  my @transcripts;
  my $alt_alleles_col  = $self->my_colour('alt_alleles_join');
 
  my (@joins,%tjoins); 
  my $tsid = $long_trans->stable_id;
  
  foreach my $gene (@$alt_alleles) {
    my $vtranscripts = $gene->get_all_Transcripts;
    my @sorted_trans = sort { $a->length <=> $b->length } @$vtranscripts;
    push @transcripts, (pop @sorted_trans);
  }
  
  if ($previous_species) {
    my ($peptide_id, $homologues, $homologue_genes) = $self->get_gene_joins($gene, $previous_species, $join_types, 'ENSEMBLGENE');
    
    if ($peptide_id) {
      foreach my $h (@$homologues) {
        push @{$tjoins{$peptide_id}},{
          key => "$h->[0]:$peptide_id",
          colour => $h->[1],
          legend => $h->[2],
        };
      }
      foreach my $h (@$homologue_genes) {
        push @{$tjoins{$peptide_id}},{
          key => "$gene_stable_id:$h->[0]",
          colour => $h->[1],
          legend => $h->[2],
        };
      }
    }
  
    my $alts = $self->filter_by_target(\@transcripts,$previous_target); 
    foreach my $t (@$alts) {
      push @joins,{
        key => join('=',$t->stable_id,$tsid),
        colour => $alt_alleles_col,
        legend => 'Alternative alleles'
      };
    }
  }
  
  if ($next_species) {
    my ($peptide_id, $homologues, $homologue_genes) = $self->get_gene_joins($gene, $next_species, $join_types, 'ENSEMBLGENE');
    
    if ($peptide_id) {
      foreach my $h (@$homologues) {
        push @{$tjoins{$peptide_id}},{
          key => "$peptide_id:$h->[0]",
          colour => $h->[1],
        };
      }
      foreach my $h (@$homologue_genes) {
        push @{$tjoins{$peptide_id}},{
          key => "$h->[0]:$gene_stable_id",
          colour => $h->[1],
        };
      }
    }
   
    my $alts = $self->filter_by_target(\@transcripts,$next_target);
    foreach my $t (@$alts) {
      push @joins,{
        key => join('=',$t->stable_id,$tsid),
        colour => $alt_alleles_col,
        legend => 'Alternative alleles'
      };
    }
  }
  $tjoins{$tsid} = \@joins;
  return \%tjoins;
}

sub _get_data {
  my ($self) = @_;

  my $hub = $self->{'config'}->hub;
  return $hub->get_query('GlyphSet::Transcript')->go($self,{
    species => $self->species,
    pattern => $self->my_config('colour_key'),
    shortlabels => $self->get_parameter('opt_shortlabels'),
    label_key => $self->my_config('label_key'),
    slice => $self->{'container'},
    logic_names => $self->my_config('logic_names'),
    db => $self->my_config('db'),
  });
}
  
sub render_collapsed {
  my ($self, $labels) = @_;

  return $self->render_text('transcript', 'collapsed') if $self->{'text_export'};
  
  my $config           = $self->{'config'};
  my $container        = $self->{'container'}{'ref'} || $self->{'container'};
  my $length           = $container->length;
  my $strand           = $self->strand;
  my $selected_db      = $self->core('db');
  my $selected_gene    = $self->my_config('g') || $self->core('g');
  my $db               = $self->my_config('db');
  my $show_labels      = $self->my_config('show_labels');
  my $link             = $self->get_parameter('compara') ? $self->my_config('join') : 0;
  my $this_db          = ($self->core('db') eq $self->my_config('db'));
 
  my $ggdraw = $self->_get_data;
  foreach my $g (@$ggdraw) {
    my @exons;
    foreach my $t (@{$g->{'transcripts'}||[]}) {
      push @exons,@{$t->{'exons'}||[]};
    }
    $g->{'exons'} = \@exons;
    delete $g->{'transcripts'};
    my $gene_stable_id = $g->{'stable_id'};
    if($this_db and $gene_stable_id eq $selected_gene) {
      $g->{'highlight'} = 'highlight2';
    }
    if($link and $gene_stable_id) {
      $g->{'joins'} = $self->calculate_collapsed_joins($gene_stable_id);
    }
  }
 
  my $draw_labels = ($labels and $show_labels ne 'off');
  $self->mr_bump($ggdraw,$draw_labels,$length);
  $self->draw_collapsed_genes($length,$draw_labels,$strand,$ggdraw);

  if($config->get_option('opt_empty_tracks') != 0 && !@$ggdraw) {
    $self->no_track_on_strand;
  }
}

sub render_transcripts {
  my ($self, $labels) = @_;

  return $self->render_text('transcript') if $self->{'text_export'};
  
  my $config            = $self->{'config'};
  my $container         = $self->{'container'}{'ref'} || $self->{'container'};
  my $length            = $container->length;
  my $strand            = $self->strand;
  my $show_labels       = $self->my_config('show_labels');
  my $link              = $self->get_parameter('compara') ? $self->my_config('join') : 0;
  my $target            = $self->get_parameter('single_Transcript');
  my $target_gene       = $self->get_parameter('single_Gene');
  my $selected_gene    = $self->my_config('g') || $self->core('g');
  my $selected_trans    = $self->my_config('t') || $self->core('pt');
  my $this_db           = ($self->core('db') eq $self->my_config('db'));
  
  my ($genes, $highlights, $transcripts, $exons) = $self->features;
  my $ggdraw = $self->_get_data;
  if($link) {
    foreach my $g (@$ggdraw) {
      next unless $g->{'stable_id'};
      my $tjoins = $self->calculate_expanded_joins($g->{'stable_id'});
      foreach my $t (@{$g->{'transcripts'}}) {
        next unless $tjoins->{$t->{'stable_id'}};
        my @joins = @{$tjoins->{$t->{'stable_id'}}};
        if($t->{'translation_stable_id'}) {
          push @joins,@{$tjoins->{$t->{'translation_stable_id'}}||[]};
        }
        $t->{'joins'} = \@joins;
      }
    }
  }
  foreach my $g (@$ggdraw) {
    foreach my $t (@{$g->{'transcripts'}}) {
      if(!defined $target and $this_db) {
        if($t->{'stable_id'} eq $selected_trans) {
          $t->{'highlight'} = 'highlight2';
        } elsif($g->{'stable_id'} eq $selected_gene) {
          $t->{'highlight'} = 'highlight1';
        }
      }
    }
  }
  my @ttdraw;
  foreach my $g (@$ggdraw) {
    foreach my $t (@{$g->{'transcripts'}}) {
      next if $target and $t->{'stable_id'} ne $target;
      next unless @{$t->{'exons'}};
      push @ttdraw,$t;
    } 
  }
  my $draw_labels = ($labels and $show_labels ne 'off');
  $self->mr_bump(\@ttdraw,$draw_labels,$length);
  $self->draw_expanded_transcripts(\@ttdraw,$length,$strand,$draw_labels);
  if($config->get_option('opt_empty_tracks') != 0 && !@$genes) {
    $self->no_track_on_strand;
  }
}

sub render_alignslice_transcript {
  my ($self, $labels) = @_;

  return $self->render_text('transcript') if $self->{'text_export'};

  my $config            = $self->{'config'};
  my $container         = $self->{'container'}{'ref'} || $self->{'container'};
  my $length            = $container->length;
  my $selected_db       = $self->core('db');
  my $strand            = $self->strand;
  my $show_labels       = $self->my_config('show_labels');
  my $target            = $self->get_parameter('single_Transcript');
  my $target_gene       = $self->get_parameter('single_Gene');
  
  my ($genes, $highlights, $transcripts) = $self->features;
  my $ggdraw = $self->_get_data;
  my @ttdraw;
  foreach my $g (@$ggdraw) {
    foreach my $t (@{$g->{'transcripts'}}) {
      next if $target and $t->{'stable_id'} ne $target;
      next unless @{$t->{'exons'}};
      $t->{'start'} = min(grep {$_} map { $_->{'start'} } @{$t->{'exons'}});
      $t->{'end'} = max(grep {$_} map { $_->{'end'} } @{$t->{'exons'}});
      push @ttdraw,$t;
    } 
  }
  
  my $draw_labels = ($labels and $show_labels ne 'off');
  $self->mr_bump(\@ttdraw,$draw_labels,$length);
  $self->draw_expanded_transcripts(\@ttdraw,$length,$strand,$draw_labels,$target);
  if($config->get_option('opt_empty_tracks') != 0 && !@$genes) {
    $self->no_track_on_strand;
  }
}

sub render_alignslice_collapsed {
  my ($self, $labels) = @_;
  
  return $self->render_text('transcript') if $self->{'text_export'};

  my $config            = $self->{'config'};
  my $container         = $self->{'container'}{'ref'} || $self->{'container'};
  my $length            = $container->length;
  my $selected_db       = $self->core('db');
  my $selected_gene     = $self->core('g');
  my $pix_per_bp        = $self->scalex;
  my $strand            = $self->strand;
  my $strand_flag       = $self->my_config('strand');
  my $db                = $self->my_config('db');
  my $show_labels       = $self->my_config('show_labels');
  my $y                 = 0;
  my $h                 = 8;
  
  my ($genes, $highlights) = $self->features;
  my $ggdraw = $self->_get_data;

  foreach my $g (@$ggdraw) {
    my @exons;
    # Collapse exons
    foreach my $t (@{$g->{'transcripts'}||[]}) {
      push @exons,@{$t->{'exons'}||[]};
    }
    $g->{'exons'} = \@exons;
    $g->{'start'} = min(grep { $_ } map { $_->{'start'} } @exons);
    $g->{'end'} = max(grep { $_ } map { $_->{'end'} } @exons);
  }

  my $draw_labels = ($labels && $show_labels ne 'off');
  $self->mr_bump($ggdraw,$draw_labels,$length);
  $self->draw_collapsed_genes($length,$draw_labels,$strand,$ggdraw);

  if($config->get_option('opt_empty_tracks') != 0 && !@$genes) {
    $self->no_track_on_strand;
  }
}

sub render_genes {
  my $self = shift;

  return $self->render_text('gene') if $self->{'text_export'};
  
  my $config           = $self->{'config'};
  my $container        = $self->{'container'}{'ref'} || $self->{'container'};
  my $length           = $container->length;
  my $pix_per_bp       = $self->scalex;
  my $strand           = $self->strand;
  my $selected_gene    = $self->my_config('g') || $self->core('g');
  my $label_threshold  = $self->my_config('label_threshold') || 50e3;
  my $navigation       = $self->my_config('navigation') || 'on';
  my $link             = $self->get_parameter('compara') ? $self->my_config('join') : 0;
  my $this_db          = ($self->core('db') eq $self->my_config('db'));
 
  my $ggdraw = $self->_get_data; 
  foreach my $g (@$ggdraw) {
    my $gene_stable_id = $g->{'stable_id'};
    if($this_db and $gene_stable_id eq $selected_gene) {
      $g->{'highlight'} = 'highlight2';
    }
    if($link) {
      $g->{'joins'} = $self->calculate_collapsed_joins($gene_stable_id);
    }
  }

  my $show_navigation = $navigation eq 'on';
  delete $_->{'href'} unless $show_navigation;
 
  my $draw_labels = $self->get_parameter('opt_gene_labels');
  $draw_labels = 1 unless defined $draw_labels;
  $draw_labels = shift if @_;
  $draw_labels = 0 if $label_threshold * 1001 < $length;

  $self->draw_rect_genes($ggdraw,$length,$draw_labels,$strand);

  if($config->get_option('opt_empty_tracks') != 0 && !@$ggdraw) {
    $self->no_track_on_strand;
  }
}

# render_text will need to be reimplemented in the manner of the above
# renderers when we restore it. The old support methods have gone. Use
# history to recover it.

# Get homologous gene ids for given gene
sub get_gene_joins {
  my ($self, $gene, $species, $join_types, $source) = @_;
  
  my $config     = $self->{'config'};
  my $compara_db = $config->hub->database('compara');
  return unless $compara_db;
  
  my $ma = $compara_db->get_GeneMemberAdaptor;
  return unless $ma;
  
  my $qy_member = $ma->fetch_by_stable_id($gene->stable_id);
  return unless defined $qy_member;
  
  my $method = $config->get_parameter('force_homologue') || $species eq $config->{'species'} ? $config->get_parameter('homologue') : undef;
  my $func   = $source ? 'get_homologous_peptide_ids_from_gene' : 'get_homologous_gene_ids';
  
  return $self->$func($species, $join_types, $compara_db->get_HomologyAdaptor, $qy_member, $method ? [ $method ] : undef);
}
  
sub get_homologous_gene_ids {
  my ($self, $species, $join_types, $homology_adaptor, $qy_member, $method) = @_;
  my @homologues;
  
  foreach my $homology (@{$homology_adaptor->fetch_all_by_Member($qy_member, -TARGET_SPECIES => [$species], -METHOD_LINK_TYPE => $method)}) {
    my $colour_key = $join_types->{$homology->description};
    
    next if $colour_key eq 'hidden';
    
    my $colour = $self->my_colour($colour_key . '_join');
    my $label  = $self->my_colour($colour_key . '_join', 'text');
    
    my $tg_member = $homology->get_all_Members()->[1];
    push @homologues, [ $tg_member->gene_member->stable_id, $colour, $label ];
  }
  
  return @homologues;
}

# Get homologous protein ids for given gene
sub get_homologous_peptide_ids_from_gene {
  my ($self, $species, $join_types, $homology_adaptor, $qy_member, $method) = @_;
  my ($stable_id, @homologues, @homologue_genes);
  
  foreach my $homology (@{$homology_adaptor->fetch_all_by_Member($qy_member, -TARGET_SPECIES => [$species], -METHOD_LINK_TYPE => $method)}) {
    my $colour_key = $join_types->{$homology->description};
    
    next if $colour_key eq 'hidden';
    
    my $colour = $self->my_colour($colour_key . '_join');
    my $label  = $self->my_colour($colour_key . '_join', 'text');
    
    $stable_id    = $homology->get_all_Members()->[0]->stable_id;
    my $tg_member = $homology->get_all_Members()->[1];
    push @homologues,      [ $tg_member->stable_id,              $colour, $label ];
    push @homologue_genes, [ $tg_member->gene_member->stable_id, $colour         ];
  }
  
  return ($stable_id, \@homologues, \@homologue_genes);
}

sub filter_by_target {
  my ($self, $alt_alleles, $target) = @_;
  
  $alt_alleles = [ grep $_->slice->seq_region_name eq $target, @$alt_alleles ] if $target;
  
  return $alt_alleles;
}

#============================================================================#
#
# Helper functions....
# 
#============================================================================#

# Generate title tag which will be used to render z-menu
sub title {
  my ($self, $transcript, $gene) = @_;
  
  my $title = 'Transcript: ' . $transcript->stable_id;
  $title .= '; Gene: ' . $gene->stable_id if $gene->stable_id;
  $title .= '; Location: ' . $transcript->seq_region_name . ':' . $transcript->seq_region_start . '-' . $transcript->seq_region_end;
  
  return $title
}

sub max_label_rows { return $_[0]->my_config('max_label_rows') || 2; }

1;