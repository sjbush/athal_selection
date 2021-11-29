use strict;
use warnings;

# REQUIREMENTS
my $main_code		 = 'code.txt'; # manually created
my $mito_code  		 = 'mito_code.txt'; # according to http://www.mun.ca/biology/scarr/MtDNA_code.html
my $cds_seq    	     = 'Arabidopsis_thaliana.TAIR10.cds.all.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-38/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz
my $gene_info  	     = 'Arabidopsis_thaliana.TAIR10.gene_info_and_go_terms.tsv'; # from Ensembl BioMart: Gene stable ID | Gene description | Chromosome/scaffold name | Gene start (bp) | Gene end (bp) | Strand | Gene name | Gene type | GO term accession | GO term name | GO term evidence code | GO domain
my $orthology  	     = 'Arabidopsis_thaliana.TAIR10.orthology_relationships_and_dnds.tsv'; # from Ensembl BioMart: Gene stable ID | Arabidopsis lyrata gene stable ID | Arabidopsis lyrata chromosome/scaffold name | Arabidopsis lyrata chromosome/scaffold start (bp) | Arabidopsis lyrata chromosome/scaffold end (bp) | Arabidopsis lyrata homology type | %id. target Arabidopsis lyrata gene identical to query gene | %id. query gene identical to target Arabidopsis lyrata gene | Arabidopsis lyrata Whole-genome alignment coverage | dN with Arabidopsis lyrata | dS with Arabidopsis lyrata | Arabidopsis lyrata orthology confidence [0 low, 1 high]
my $cds_coords		 = 'Arabidopsis_thaliana.TAIR10.cds_coords.tsv'; # from 0.obtain_TAIR10_CDS_coords.pl
my $genome			 = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
my $lyrata_cds 	     = 'Arabidopsis_lyrata.v.1.0.cds.all.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-38/plants/fasta/arabidopsis_lyrata/cds/Arabidopsis_lyrata.v.1.0.cds.all.fa.gz
my $lyrata_strands   = 'Arabidopsis_lyrata.v.1.0.strands.tsv'; # from Ensembl BioMart: Gene stable ID | Strand
my $polymorphisms_19 = 'polymorphisms_19_accs.txt'; # reformatted data originally from Gan, et al. DOI: 10.1038/nature10414
my $polymorphisms_80 = 'polymorphisms_80_accs.txt'; # reformatted data originally from Cao, et al. DOI: 10.1038/ng.91
my %trans; my %mito_trans;
LOAD_MAIN_CODE($main_code);
LOAD_MITO_CODE($mito_code);

# TEMPORARY
my $temp_dir = 'CDS_alignment_temp_dir';
if (!(-d($temp_dir))) { mkdir $temp_dir or die $!; }
my $ctlFileName = "yn00.ctl";
my $tempNuc     = "temp.nuc";

# OUTPUT
my $out_file = 'Arabidopsis_thaliana.TAIR10.sequence_evolution.tsv';
open(OUT,'>',$out_file) or die $!;
print OUT "Gene name\tGene ID\tGene description\tGene type\tGene location (chr:start-end:strand)\tGene length (bp)\tNo. of protein-coding transcripts\tRepresentative transcript ID (transcript with longest CDS)\tCDS length (bp)\tExon IDs\tExon coords\tTotal exon length\t% of gene that is exonic\tCDS coords\t% of gene that is CDS\tGO terms (F)\tGO terms (P)\tGO terms (C)\tOrthologous A. lyrata gene ID\tA. lyrata gene location (chr:start-end:strand)\tA. lyrata gene length (bp)\tNo. of protein-coding transcripts in A. lyrata orthologue\tRepresentative transcript ID (transcript with longest CDS)\tCDS length (bp)\tHomology type\tReciprocal % gene identity (Atha to Alyr, Alyr to Atha)\tdN (calculated using whole gene sequence, obtained via Ensembl)\tdS (calculated using whole gene sequence, obtained via Ensembl)\tdN/dS (calculated using whole gene sequence, obtained via Ensembl)\tIs this a high confidence gene-level dN/dS estimate? (criteria: a one-to-one orthology relationship with A. lyrata of >= 75% reciprocal identity, and dS > 0.02, dS < 2 and dN < 2)\t% identity in alignment of longest A. thaliana CDS with longest A. lyrata CDS (EMBOSS needle)\t% gaps in alignment of longest A. thaliana CDS with longest A. lyrata CDS (EMBOSS needle)\tNo. of completely aligned codons\tdN (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75%)\tdS (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75%)\tdN/dS (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75%)\tIs this a high confidence CDS-level dN/dS estimate? (criteria: as gene-level high-confidence criteria, plus CDS > 150bp, >= 75% CDS identity with A. lyrata, dS > 0.02, dS < 2 and dN < 2)\tNo. of non-synonymous sites\tNo. of synonymous sites\tNo. of non-synonymous substitutions\tNo. of synonymous substitutions\tNI (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75% [polymorphism data from 19 accessions])\tNI (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75% [polymorphism data from 80 accessions])\tDOS (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75%) [polymorphism data from 19 accessions]\tDOS (calculated using completely aligned codons in CDS, by PAML yn00, providing CDS > 150bp and identity vs. A. lyrata >= 75%) [polymorphism data from 80 accessions]\n";

# STORE EXON AND CDS COORDS
my %exon_and_cds_coords = ();
open(IN,$cds_coords) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $transcript_id = $line[0]; my $transcript_coords = $line[1]; my $exon_ids = $line[5]; my $exon_coords = $line[6]; my $total_exon_length = $line[7]; my $pc_of_transcript_that_is_exonic = $line[8]; my $cds_coords = $line[9]; my $cds_length = $line[10]; my $pc_of_transcript_that_is_cds = $line[11];
	  $exon_and_cds_coords{$transcript_id}{exon_ids}    = $exon_ids;
	  $exon_and_cds_coords{$transcript_id}{cds_coords}  = $cds_coords;
	  $exon_and_cds_coords{$transcript_id}{cds_length}  = $cds_length;
	  $exon_and_cds_coords{$transcript_id}{exon_coords} = $exon_coords;
	  $exon_and_cds_coords{$transcript_id}{total_exon_length} = $total_exon_length;
	  $exon_and_cds_coords{$transcript_id}{pc_of_transcript_that_is_cds} = $pc_of_transcript_that_is_cds;
	  $exon_and_cds_coords{$transcript_id}{pc_of_transcript_that_is_exonic} = $pc_of_transcript_that_is_exonic;
	}
close(IN) or die $!;

# STORE A. THALIANA GENOME
my $chr = ''; my %col_seq = ();
open(IN,$genome) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $chr = $1;
		  if ($chr =~ /Mt/) { $chr = 'M'; }
		  if ($chr =~ /Pt/) { $chr = 'C'; }
		}
	  next if ($line =~ /^\>/);
	  $col_seq{$chr} .= $line;
	}
close(IN) or die $!;

# STORE A. THALIANA CDS
my %cds = ();
my $transcript_id = ''; $chr = ''; my $cds_start = ''; my $cds_end = ''; my $strand = ''; my $gene_id = '';
open(IN,$cds_seq) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) cds .*?\:.*?\:(.*?):(\d+):(\d+):(.*?) gene\:(.*?) .*?$/) # most of these headers read "cds chromosome:TAIR" but NOT all
		{ $transcript_id = $1; $chr = $2; $cds_start = $3; $cds_end = $4; $strand = $5; $gene_id = $6; }
	  next if ($line =~ /^\>/);
	  $cds{$gene_id}{$transcript_id}{cds_coords} = "$chr:$cds_start-$cds_end:$strand"; # IMPORTANT: these are NOT the CDS coords but the transcript coords
	  $cds{$gene_id}{$transcript_id}{seq} 	    .= $line;
	}
close(IN) or die $!;

# STORE A. LYRATA CDS
my %lyrata_cds = ();
$transcript_id = ''; $chr = ''; $cds_start = ''; $cds_end = ''; $strand = ''; $gene_id = '';
open(IN,$lyrata_cds) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) cds .*?\:.*?\:(.*?):(\d+):(\d+):(.*?) gene\:(.*?) .*?$/) # most of these headers read "cds chromosome|scaffold:v.1.0" but NOT all
		{ $transcript_id = $1; $chr = $2; $cds_start = $3; $cds_end = $4; $strand = $5; $gene_id = $6; }
	  next if ($line =~ /^\>/);
	  $lyrata_cds{$gene_id}{$transcript_id}{cds_coords} = "$chr:$cds_start-$cds_end:$strand"; # IMPORTANT: these are NOT the CDS coords but the transcript coords
	  $lyrata_cds{$gene_id}{$transcript_id}{seq} 	   .= $line;
	}
close(IN) or die $!;

# STORE A. LYRATA STRANDS
my %lyrata_strands = ();
open(IN,$lyrata_strands) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $transcript_id = $line[1]; my $strand = $line[2];
	  $lyrata_strands{$gene_id} = $strand;
	}
close(IN) or die $!;

# STORE A. THALIANA GENE INFO
my %is_gene_name_duplicated = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[6];
	  $is_gene_name_duplicated{$gene_name}{$gene_id}++ unless ($gene_name eq '');
	}
close(IN) or die $!;
my %gene_info = (); my %go_cats = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_desc = $line[1]; my $chr = $line[2]; my $gene_start = $line[3]; my $gene_end = $line[4]; my $strand = $line[5]; my $gene_name = $line[6]; my $gene_type = $line[7];
	  my $go_term = $line[8]; my $go_term_name = $line[9]; my $evidence_code = $line[10]; my $go_domain = $line[11];
	  my $num_times_this_gene_name_appears = 0;
	  if (exists($is_gene_name_duplicated{$gene_name})) { $num_times_this_gene_name_appears = scalar keys %{$is_gene_name_duplicated{$gene_name}}; }
	  if ($num_times_this_gene_name_appears > 1)
		{ $gene_name = "$gene_name/$gene_id";
		}
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  my $gene_len = ($gene_end-$gene_start)+1;
	  $gene_info{$gene_name}{loc} 		= "$chr:$gene_start-$gene_end:$strand";
	  $gene_info{$gene_name}{chr}		= $chr;
	  $gene_info{$gene_name}{gene_id}   = $gene_id;
	  $gene_info{$gene_name}{gene_len}  = $gene_len;
	  $gene_info{$gene_name}{gene_desc} = $gene_desc;
	  $gene_info{$gene_name}{gene_type} = $gene_type;
	  if ((defined($line[8])) and (defined($line[9])) and (defined($line[10])) and (defined($line[11])))
		{ next if (($evidence_code eq 'NAS') or ($evidence_code eq 'ND'));
		  $go_cats{$gene_name}{$go_domain}{"$go_term_name ($go_term)"}++;
		}
	}
close(IN) or die $!;

# GET SNP DATA FOR THE 19 GENOMES, USED TO CALCULATE NI AND TAJIMA'S D
my %polymorphisms_19_for_ni = ();
open(IN,$polymorphisms_19) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  my $chr = $line[3]; my $position = $line[4];
	  
	  # STORE SNP DATA FOR CALCULATING NI
	  my $before_codon = $line[8]; my $after_codon = $line[9];
	  my $before_aa = $line[10]; my $after_aa = $line[11]; # these are based on strand, irrespective of how the codon is reported! e.g. SNP at 2:9846086 is an I->V, supposedly for AAT->AAC (actually it's ATT->GTT as the gene, AT2G23140, is on - strand)
	  my $res1 = $line[14]; my $res2 = $line[15]; my $res3 = $line[16];
	  my @res = ("$res1","$res2","$res3");
	  if (!(exists($polymorphisms_19_for_ni{$gene_id}{$position}))) # if data DOES NOT exist for the same SNP on multiple transcripts...
		{ $polymorphisms_19_for_ni{$gene_id}{$position}{BeforeCodon} = $before_codon;
		  $polymorphisms_19_for_ni{$gene_id}{$position}{AfterCodon}  = $after_codon;
		  $polymorphisms_19_for_ni{$gene_id}{$position}{BeforeAA} 	 = $before_aa;
		  $polymorphisms_19_for_ni{$gene_id}{$position}{AfterAA}  	 = $after_aa;
		  foreach my $res (@res)
			{ push(@{$polymorphisms_19_for_ni{$gene_id}{$position}{BasesForThisCodon}},$res); }
		}
	}
close(IN) or die $!;

# GET SNP DATA FOR THE 80 GENOMES, USED TO CALCULATE NI AND TAJIMA'S D
my %polymorphisms_80_for_ni = ();
open(IN,$polymorphisms_80) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if ($line eq '');
	  my $gene_id = $line[21];
	  my $chr = $line[0]; my $position = $line[1];
	  
	  # STORE SNP DATA FOR CALCULATING NI
	  my $before_codon = $line[7]; my $after_codon = $line[8];
	  my $before_aa = $line[9]; my $after_aa = $line[10];
	  my $res1 = $line[4]; my $res2 = $line[5]; my $res3 = $line[6];
	  my @res = ("$res1","$res2","$res3");
	  $polymorphisms_80_for_ni{$gene_id}{$position}{BeforeCodon} = $before_codon;
	  $polymorphisms_80_for_ni{$gene_id}{$position}{AfterCodon}  = $after_codon;
	  $polymorphisms_80_for_ni{$gene_id}{$position}{BeforeAA} 	 = $before_aa;
	  $polymorphisms_80_for_ni{$gene_id}{$position}{AfterAA}  	 = $after_aa;
	  foreach my $res (@res)
		{ push(@{$polymorphisms_80_for_ni{$gene_id}{$position}{BasesForThisCodon}},$res); }
	}
close(IN) or die $!;

# STORE A. THALIANA TO A. LYRATA ORTHOLOGY RELATIONSHIPS
my %orthology_relationships_and_dnds = ();
open(IN,$orthology) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  next if (!(defined($line[1])));
	  my $a_lyrata_gene_id  = $line[1];
	  my $a_lyrata_coords   = "$line[2]:$line[3]-$line[4]:$lyrata_strands{$a_lyrata_gene_id}";
	  my $a_lyrata_gene_len = ($line[4]-$line[3])+1;
	  my $homology_type     = $line[5];
	  my $lyr_to_thal_pc    = sprintf("%.1f",$line[6]);
	  my $thal_to_lyr_pc    = sprintf("%.1f",$line[7]);
	  my $whole_genome_aln  = 'NA'; my $dn = 'NA'; my $ds = 'NA';
	  if (defined($line[8]))  { $whole_genome_aln = $line[8];  }
	  if (defined($line[9]))  { $dn 			  = sprintf("%.4f",$line[9]);  }
	  if (defined($line[10])) { $ds 			  = sprintf("%.4f",$line[10]); }
	  $orthology_relationships_and_dnds{$gene_id}{a_lyrata_gene_len} = $a_lyrata_gene_len;
	  $orthology_relationships_and_dnds{$gene_id}{whole_genome_aln}  = $whole_genome_aln;
	  $orthology_relationships_and_dnds{$gene_id}{a_lyrata_gene_id}  = $a_lyrata_gene_id;
	  $orthology_relationships_and_dnds{$gene_id}{a_lyrata_coords}   = $a_lyrata_coords;
	  $orthology_relationships_and_dnds{$gene_id}{thal_to_lyr_pc}    = $thal_to_lyr_pc;
	  $orthology_relationships_and_dnds{$gene_id}{lyr_to_thal_pc}    = $lyr_to_thal_pc;
	  $orthology_relationships_and_dnds{$gene_id}{homology_type}     = $homology_type;
	  $orthology_relationships_and_dnds{$gene_id}{pc_identity}       = "atha to alyr: $thal_to_lyr_pc%, alyr to atha: $lyr_to_thal_pc%";
	  $orthology_relationships_and_dnds{$gene_id}{dn} 				 = $dn;
	  $orthology_relationships_and_dnds{$gene_id}{ds} 				 = $ds;
	  my $dnds = 'NA';
	  if (($dn =~ /\d/) && ($ds =~ /\d/) && ($ds > 0)) { $dnds = sprintf("%.4f",($dn/$ds)); }
	  $orthology_relationships_and_dnds{$gene_id}{dnds} = $dnds;
	}
close(IN) or die $!;

# IDENTIFY HIGH-CONFIDENCE THALIANA-LYRATA PAIRS, THEN ALIGN THE CDS AND CALCULATE dN/dS, NI AND DOS
my @gene_names = ();
while((my $gene_name,my $irrel)=each(%gene_info))
	{ push(@gene_names,$gene_name); }
my @sorted_gene_names = sort {"\L$a" cmp "\L$b"} @gene_names;
my $gene_names_seen = 0; my $gene_names_total = @sorted_gene_names;
foreach my $gene_name (@sorted_gene_names)
	{ $gene_names_seen++;
	  print "$gene_names_seen of $gene_names_total...\n";
	  my $chr		= $gene_info{$gene_name}{chr};
	  my $loc		= $gene_info{$gene_name}{loc};
	  my $gene_id   = $gene_info{$gene_name}{gene_id};
	  my $gene_len  = $gene_info{$gene_name}{gene_len};
	  my $gene_desc = $gene_info{$gene_name}{gene_desc};
	  my $gene_type = $gene_info{$gene_name}{gene_type};
	  my $go_f_line = ''; my $go_p_line = ''; my $go_c_line = '';
	  if (exists($go_cats{$gene_name}{'molecular_function'}))
		{ my @go_terms = ();
		  while((my $go_term,my $irrel)=each(%{$go_cats{$gene_name}{'molecular_function'}}))
			{ push(@go_terms,$go_term); }
		  my @sorted_go_terms = sort {$a cmp $b} @go_terms;
		  foreach my $go_term (@sorted_go_terms)
			{ $go_f_line .= "$go_term | "; }
		}
	  if (exists($go_cats{$gene_name}{'biological_process'}))
		{ my @go_terms = ();
		  while((my $go_term,my $irrel)=each(%{$go_cats{$gene_name}{'biological_process'}}))
			{ push(@go_terms,$go_term); }
		  my @sorted_go_terms = sort {$a cmp $b} @go_terms;
		  foreach my $go_term (@sorted_go_terms)
			{ $go_p_line .= "$go_term | "; }
		}
	  if (exists($go_cats{$gene_name}{'cellular_component'}))
		{ my @go_terms = ();
		  while((my $go_term,my $irrel)=each(%{$go_cats{$gene_name}{'cellular_component'}}))
			{ push(@go_terms,$go_term); }
		  my @sorted_go_terms = sort {$a cmp $b} @go_terms;
		  foreach my $go_term (@sorted_go_terms)
			{ $go_c_line .= "$go_term | "; }
		}
	  $go_f_line =~ s/ \| $//; $go_p_line =~ s/ \| $//; $go_c_line =~ s/ \| $//;
	  if ($go_f_line eq '') { $go_f_line = 'unknown or unavailable'; }
	  if ($go_p_line eq '') { $go_p_line = 'unknown or unavailable'; }
	  if ($go_c_line eq '') { $go_c_line = 'unknown or unavailable'; }
	  my $no_of_transcripts_with_cds = 0;
	  if (exists($cds{$gene_id}))
		{ $no_of_transcripts_with_cds = scalar keys %{$cds{$gene_id}}; }
	  my @cds_lengths = ();
	  while((my $transcript_id,my $irrel)=each(%{$cds{$gene_id}}))
		{ my $seq 		 = $cds{$gene_id}{$transcript_id}{seq};
		  my $cds_coords = $cds{$gene_id}{$transcript_id}{cds_coords};
		  my $cds_length = length($seq);
		  push(@cds_lengths,[$cds_length,$transcript_id,$cds_coords]);
		}
	  my $transcript_id_with_longest_cds = 'NA'; my $length_of_longest_cds = 'NA'; my $coords_of_longest_cds = 'NA';
	  if ($#cds_lengths != -1)
		{ my @sorted_cds_lengths = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @cds_lengths;
		  $transcript_id_with_longest_cds = $sorted_cds_lengths[0][1];
		  $length_of_longest_cds 		  = $sorted_cds_lengths[0][0];
		  $coords_of_longest_cds 		  = $sorted_cds_lengths[0][2];
		}
	  my $a_lyrata_gene_id = 'NA'; my $a_lyrata_coords = 'NA'; my $a_lyrata_gene_len = 'NA'; my $no_of_lyrata_transcripts_with_cds = 'NA'; my $lyrata_transcript_id_with_longest_cds = 'NA'; my $length_of_longest_lyrata_cds = 'NA'; my $coords_of_longest_lyrata_cds = 'NA';
	  my $homology_type = 'NA'; my $pc_identity = 'NA'; my $whole_genome_aln = 'NA'; my $dn = 'NA'; my $ds = 'NA'; my $dnds = 'NA'; my $dn_cds = 'NA'; my $ds_cds = 'NA'; my $dnds_cds = 'NA'; my $n_sites_cds = 'NA'; my $s_sites_cds = 'NA'; my $n_subs_cds = 'NA'; my $s_subs_cds = 'NA';
	  my $pc_identity_in_aln = 'NA'; my $pc_gaps_in_aln = 'NA'; my $thaliana_aln = 'NA'; my $lyrata_aln = 'NA';
	  my $number_of_aligned_codons = 'NA'; my $thaliana_pep_from_alignment = 'NA'; my $lyrata_pep_from_alignment = 'NA'; my $thaliana_codons_from_alignment = ''; my $lyrata_codons_from_alignment = '';
	  my $is_high_confidence_gene_level = ''; my $is_high_confidence_cds_level = '';
	  my $ni_cds_19 = 'NA'; my $dos_cds_19 = 'NA'; my $tajimas_d_cds_19 = 'NA';
	  my $ni_cds_80 = 'NA'; my $dos_cds_80 = 'NA'; my $tajimas_d_cds_80 = 'NA';
	  my $exon_ids    					  = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{exon_ids})) 						  { $exon_ids 						 = $exon_and_cds_coords{$transcript_id_with_longest_cds}{exon_ids}; 					   }
	  my $cds_coords  					  = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{cds_coords})) 					  { $cds_coords 					 = $exon_and_cds_coords{$transcript_id_with_longest_cds}{cds_coords}; 					   }
	  my $cds_length					  = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{cds_length})) 					  { $cds_length						 = $exon_and_cds_coords{$transcript_id_with_longest_cds}{cds_length}; 					   }
	  my $exon_coords 					  = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{exon_coords})) 					  { $exon_coords 					 = $exon_and_cds_coords{$transcript_id_with_longest_cds}{exon_coords}; 					   }
	  my $total_exon_length				  = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{total_exon_length})) 			  { $total_exon_length 				 = $exon_and_cds_coords{$transcript_id_with_longest_cds}{total_exon_length}; 			   }
	  my $pc_of_transcript_that_is_cds    = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{pc_of_transcript_that_is_cds})) 	  { $pc_of_transcript_that_is_cds    = $exon_and_cds_coords{$transcript_id_with_longest_cds}{pc_of_transcript_that_is_cds};    }
	  my $pc_of_transcript_that_is_exonic = 'NA'; if (exists($exon_and_cds_coords{$transcript_id_with_longest_cds}{pc_of_transcript_that_is_exonic})) { $pc_of_transcript_that_is_exonic = $exon_and_cds_coords{$transcript_id_with_longest_cds}{pc_of_transcript_that_is_exonic}; }
	  
	  # FOR HIGH-CONFIDENCE ORTHOLOGOUS PAIRS, CALCULATE dN/dS AND NI USING THE CDS
	  if (exists($orthology_relationships_and_dnds{$gene_id}))
		{ $a_lyrata_gene_len = $orthology_relationships_and_dnds{$gene_id}{a_lyrata_gene_len};
		  $a_lyrata_gene_id  = $orthology_relationships_and_dnds{$gene_id}{a_lyrata_gene_id};
		  $a_lyrata_coords   = $orthology_relationships_and_dnds{$gene_id}{a_lyrata_coords};
		  $homology_type     = $orthology_relationships_and_dnds{$gene_id}{homology_type};
		  $pc_identity       = $orthology_relationships_and_dnds{$gene_id}{pc_identity};
		  $whole_genome_aln  = $orthology_relationships_and_dnds{$gene_id}{whole_genome_aln};
		  $dn 			     = $orthology_relationships_and_dnds{$gene_id}{dn};
		  $ds 				 = $orthology_relationships_and_dnds{$gene_id}{ds};
		  $dnds 			 = $orthology_relationships_and_dnds{$gene_id}{dnds};
		  my $thal_to_lyr_pc = $orthology_relationships_and_dnds{$gene_id}{thal_to_lyr_pc};
		  my $lyr_to_thal_pc = $orthology_relationships_and_dnds{$gene_id}{lyr_to_thal_pc};
		  
		  $no_of_lyrata_transcripts_with_cds = 0;
		  if (exists($lyrata_cds{$a_lyrata_gene_id}))
			{ $no_of_lyrata_transcripts_with_cds = scalar keys %{$lyrata_cds{$a_lyrata_gene_id}}; }
		  my @lyrata_cds_lengths = ();
		  while((my $lyrata_transcript_id,my $irrel)=each(%{$lyrata_cds{$a_lyrata_gene_id}}))
			{ my $lyrata_seq 		= $lyrata_cds{$a_lyrata_gene_id}{$lyrata_transcript_id}{seq};
			  my $lyrata_cds_coords = $lyrata_cds{$a_lyrata_gene_id}{$lyrata_transcript_id}{cds_coords};
			  my $lyrata_cds_length = length($lyrata_seq);
			  push(@lyrata_cds_lengths,[$lyrata_cds_length,$lyrata_transcript_id,$lyrata_cds_coords]);
			}
		  if ($#lyrata_cds_lengths != -1)
			{ my @sorted_lyrata_cds_lengths = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @lyrata_cds_lengths;
			  $lyrata_transcript_id_with_longest_cds = $sorted_lyrata_cds_lengths[0][1];
			  $length_of_longest_lyrata_cds 		 = $sorted_lyrata_cds_lengths[0][0];
			  $coords_of_longest_lyrata_cds 		 = $sorted_lyrata_cds_lengths[0][2];
			}
		  
		  # IDENTIFY HIGH-CONFIDENCE THALIANA/LYRATA GENE PAIRS AND ALIGN THEIR CDS WITH EMBOSS NEEDLE
		  if (($dnds =~ /\d/) and ($ds > 0.02) and ($ds < 2) and ($dn < 2) and ($homology_type eq 'ortholog_one2one') and ($thal_to_lyr_pc >= 75) and ($lyr_to_thal_pc >= 75))
			{ $is_high_confidence_gene_level = '*';
			  next if ( (!(exists($cds{$gene_id}{$transcript_id_with_longest_cds}{seq}))) or (!(exists($lyrata_cds{$a_lyrata_gene_id}{$lyrata_transcript_id_with_longest_cds}{seq}))) );
			  my $thaliana_seq = $cds{$gene_id}{$transcript_id_with_longest_cds}{seq};
			  my $lyrata_seq   = $lyrata_cds{$a_lyrata_gene_id}{$lyrata_transcript_id_with_longest_cds}{seq};
			  my $temp_file1   = "$temp_dir/thaliana.$gene_id.fa";
			  my $temp_file2   = "$temp_dir/lyrata.$a_lyrata_gene_id.fa";
			  my $temp_file3   = "$temp_dir/$gene_id.$a_lyrata_gene_id.txt";
			  open (TEMP,'>',$temp_file1) or die $!;
			  print TEMP "$thaliana_seq\n";
			  close(TEMP) or die $!;
			  open (TEMP,'>',$temp_file2) or die $!;
			  print TEMP "$lyrata_seq\n";
			  close(TEMP) or die $!;
			  my $pid;
			  eval
				{ local $SIG{ALRM} = sub { die "alarm\n"; };
				  alarm 1000;
				  if ($pid)
					{ while(1)
						{ TIMEOUT($pid) and next;
						  sleep 1;
						}
					  waitpid($pid,0);
					}
				  else
					{ (system("needle $temp_file1 $temp_file2 -gapopen 10.0 -gapextend 0.5 -outfile $temp_file3") == 0) or die $!;
					}
				  alarm 0;
				};
			  if ($@)
				{ if ($@ ne "alarm\n") { next; }
				}
			  unlink $temp_file1 or die $!; unlink $temp_file2 or die $!;
			  my $type = 'thaliana';
			  my %alns = ();
			  open(IN,$temp_file3) or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  if ($line =~ /^.*?Identity\:.*?\((.*?)%\)$/)
					{ $pc_identity_in_aln = $1; }
				  if ($line =~ /^.*?Gaps\:.*?\((.*?)%\)$/)
					{ $pc_gaps_in_aln = $1; }
				  if ($line =~ /^\s+\d/)
					{ $alns{$type} .= $line;
					  if 	($type eq 'thaliana') { $type = 'lyrata';   }
					  elsif ($type eq 'lyrata')   { $type = 'thaliana'; }
					}
				}
			  close(IN) or die $!;
			  $thaliana_aln = $alns{'thaliana'};
			  $lyrata_aln   = $alns{'lyrata'};
			  $thaliana_aln =~ s/\d+//g; $lyrata_aln =~ s/\d+//g;
			  $thaliana_aln =~ s/\s+//g; $lyrata_aln =~ s/\s+//g;
			  unlink $temp_file3 or die $!;
			  
			  # CALCULATE SEVERAL SIGNATURES OF SELECTION (RECALL THAT TAJIMA'S D HAS BEEN CALCULATED ABOVE)
			  # (1) dN/dS ACROSS THE CDS, CALCULATED FROM THE SET OF USABLE CODONS IN EACH ALIGNMENT USING PAML yn00
			  my @thaliana_aln = split(//,$thaliana_aln);
			  my @lyrata_aln   = split(//,$lyrata_aln);
			  my @aln1 = @thaliana_aln; my @aln2 = @lyrata_aln;
			  my $translated_orf_1 = ''; my $usable_codons_1 = '';
			  my $translated_orf_2 = ''; my $usable_codons_2 = '';
			  my $n_subs = 0; my $s_subs = 0;
			  for(my $x=0;$x<@aln1;$x++)
				{ my $nt1_1 = $aln1[$x]; my $nt1_2 = $aln2[$x];
				  next if ( (!(defined($aln1[$x+1]))) || (!(defined($aln1[$x+2]))) || (!(defined($aln2[$x+1]))) || (!(defined($aln2[$x+2]))) );
				  my $nt2_1 = $aln1[$x+1]; my $nt3_1 = $aln1[$x+2];
				  my $nt2_2 = $aln2[$x+1]; my $nt3_2 = $aln2[$x+2];
				  if ($nt1_1 !~ /^[ATCGN]$/) { $nt1_1 = 'N'; } # replace all ambiguity characters (R, Y, S, W, K, M, B, D, H, V) with N
				  if ($nt2_1 !~ /^[ATCGN]$/) { $nt2_1 = 'N'; }
				  if ($nt3_1 !~ /^[ATCGN]$/) { $nt3_1 = 'N'; }
				  if ($nt1_2 !~ /^[ATCGN]$/) { $nt1_2 = 'N'; }
				  if ($nt2_2 !~ /^[ATCGN]$/) { $nt2_2 = 'N'; }
				  if ($nt3_2 !~ /^[ATCGN]$/) { $nt3_2 = 'N'; }
				  my $codon_1 = "$nt1_1"."$nt2_1"."$nt3_1";
				  my $codon_2 = "$nt1_2"."$nt2_2"."$nt3_2";
				  next if (($codon_1 =~ /-/) || ($codon_2 =~ /-/));
				  my $aa_1 = 'X'; if ($chr eq 'Mt') { $aa_1 = $mito_trans{$codon_1} unless (!(exists($mito_trans{$codon_1}))); } else { $aa_1 = $trans{$codon_1} unless (!(exists($trans{$codon_1}))); }
				  my $aa_2 = 'X'; if ($chr eq 'Mt') { $aa_2 = $mito_trans{$codon_2} unless (!(exists($mito_trans{$codon_2}))); } else { $aa_2 = $trans{$codon_2} unless (!(exists($trans{$codon_2}))); }
				  if (($aa_1 ne 'X') and ($aa_2 ne 'X') and ($aa_1 ne '*') and ($aa_2 ne '*'))
					{ $translated_orf_1 .= $aa_1;
					  $translated_orf_2 .= $aa_2;
					  $usable_codons_1  .= $codon_1;
					  $usable_codons_2  .= $codon_2;
					  if ($codon_1 ne $codon_2)
						{ if ($aa_1 eq $aa_2)
							{ $s_subs++; }
						  elsif ($aa_1 ne $aa_2)
							{ $n_subs++; }
						}
					}
				  $x = $x+2;
				}
			  $thaliana_pep_from_alignment = $translated_orf_1; $thaliana_codons_from_alignment = $usable_codons_1;
			  $lyrata_pep_from_alignment   = $translated_orf_2; $lyrata_codons_from_alignment   = $usable_codons_2;
			  $number_of_aligned_codons = length($thaliana_pep_from_alignment);
			  my $sequence_length = length($thaliana_codons_from_alignment);
			  my $controlFile = <<ENDFILE;
								  seqfile = temp.nuc * sequence data file name
								  outfile = yn           * main result file
								  verbose = 0  * 1: detailed output (list sequences), 0: concise output

									icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

								weighting = 0  * weighting pathways between codons (0/1)?
							   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
							*       ndata = 1


							* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
							* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
							* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
							* 10: blepharisma nu.
							* These codes correspond to transl_table 1 to 11 of GENEBANK.

ENDFILE

			  chdir $temp_dir or die $!;
			  open(CTLFILE,">",$ctlFileName); # see also http://www.csi.uoregon.edu/projects/genetics/ntdiffs/ (mldiffs)
			  print CTLFILE $controlFile;
			  close(CTLFILE);
			  open(TEMP,">","$tempNuc") or die $!;
			  print TEMP "  2    $sequence_length\nTha       $thaliana_codons_from_alignment\nLyr       $lyrata_codons_from_alignment"; # imperative, for sequential Phylip format, that there is no last newline
			  close(TEMP) or die $!;
			  my $pid2;
			  eval
				{ local $SIG{ALRM} = sub { die "alarm\n" };
				  alarm 100;
				  if ($pid)
					{ while(1)
						{ TIMEOUT($pid2) and last;
						  sleep 1;
						}
					  waitpid($pid2,0);
					}
				  else
					{ (system("yn00") == 0) or die $!;
					}
				  alarm 0;
				};
			  my $yn_line = ''; my $yn_line_num = 0;
			  open(YN,'yn') or die $!;
			  while(<YN>)
				{ my $line = $_; chomp($line);
				  if ($line eq '(equal weighting of pathways)')
					{ $yn_line_num = $. ; }
				  if ((defined($yn_line_num)) and ($. == $yn_line_num+4))
					{ $yn_line = $line;
					}
				}
			  close(YN) or die $!;
			  my $dnds_line = '';
			  open(RES,'rst') or die $!;
			  while(<RES>)
				{ my $line = $_; chomp($line);
				  $dnds_line = $line;
				  next;
				}
			  close(RES) or die $!;
			  $dnds_line =~ s/\s+/ /g;
			  my $dn = 'NA'; my $ds = 'NA'; my $dnds = 'NA'; my $n_sites = 'NA'; my $s_sites = 'NA';
			  if (($sequence_length > 150) and ($pc_identity_in_aln >= 75)) # FILTER: retain only high-confidence CDS alignments; >150bp with >=75% identity
				{ my @dnds_line = split(/ /,$dnds_line);
				  if ( ((defined($dnds_line[5]))) || (defined($dnds_line[8])) )
					{ $dn = $dnds_line[5];
					  $ds = $dnds_line[8];
					  if (($dn !~ /\-/) || ($ds !~ /\-/))
						{ if ($ds =~ /\d+/)
							{ if ($ds != 0)
								{ $dnds = sprintf("%.4f",($dn/$ds));
								}
							}
						}
					}
				  $yn_line =~ s/\s+/\t/g;
				  my @yn_line = split(/\t/,$yn_line);
				  $s_sites = $yn_line[3];
				  $n_sites = $yn_line[4];
				  if (($dn =~ /\d/) and ($ds > 0.02) and ($ds < 2) and ($dn < 2))
					{ $is_high_confidence_cds_level = '*';
					}
				}
			  chdir ".." or die $!;
			  $dn_cds = $dn; $ds_cds = $ds; $dnds_cds = $dnds; $n_sites_cds = $n_sites; $s_sites_cds = $s_sites; $n_subs_cds = $n_subs; $s_subs_cds = $s_subs;
			  
			  # (2) CALCULATE NI FOR THE ALIGNED CDS, AFTER DETERMINING THE NUMBER OF SYN AND NON-SYN POLYMORPHISMS WITHIN THE SET OF CDS COORDS
			  my %cds_bases = ();
			  if ($cds_coords ne 'NA')
				{ my @cds_coords = split(/\, /,$cds_coords);
				  foreach my $coords (@cds_coords)
					{ if ($coords =~ /^.*?\:(\d+)\-(\d+)$/)
						{ my $start = $1; my $end = $2;
						  for(my $x=$start;$x<=$end;$x++)
							{ $cds_bases{$x}++;
							}
						}
					}
				}
			  
			  for(my $z=0;$z<=1;$z++)
				{ my %polymorphisms_for_ni = ();
				  if 	($z == 0) { %polymorphisms_for_ni = %polymorphisms_19_for_ni; }
				  elsif ($z == 1) { %polymorphisms_for_ni = %polymorphisms_80_for_ni; }
	  
				  my $s_poly = 0; my $n_poly = 0;
				  if (exists($polymorphisms_for_ni{$gene_id}))
					{ while((my $position,my $irrel)=each(%{$polymorphisms_for_ni{$gene_id}}))
						{ if (exists($cds_bases{$position}))
							{ my @pos = @{$polymorphisms_for_ni{$gene_id}{$position}{BasesForThisCodon}};
							  my $codon_from_this_seq = '';
							  foreach my $pos (@pos)
								{ my $base = substr($col_seq{$chr},$pos-1,1);
								  $codon_from_this_seq .= $base;
								}
							  my $before_codon = $polymorphisms_for_ni{$gene_id}{$position}{BeforeCodon};
							  my $after_codon  = $polymorphisms_for_ni{$gene_id}{$position}{AfterCodon};
							  next if (($codon_from_this_seq ne $before_codon) and ($codon_from_this_seq ne $after_codon));
							  my $before_aa = $polymorphisms_for_ni{$gene_id}{$position}{BeforeAA};
							  my $after_aa  = $polymorphisms_for_ni{$gene_id}{$position}{AfterAA};
							  if 	($before_aa ne $after_aa) { $n_poly++; }
							  elsif ($before_aa eq $after_aa) { $s_poly++; }
							}
						}
					}
				  if ($z == 0)
					{ $ni_cds_19  = ( ((2*$s_subs)+1)*((2*$s_poly)+1) ) / ( ((2*$n_subs)+1)*((2*$n_poly)+1) );
					  my $lni_cds = LOG10($ni_cds_19);
					  $ni_cds_19 = sprintf("%.4f",$lni_cds);
					}
				  elsif ($z == 1)
					{ $ni_cds_80  = ( ((2*$s_subs)+1)*((2*$s_poly)+1) ) / ( ((2*$n_subs)+1)*((2*$n_poly)+1) );
					  my $lni_cds = LOG10($ni_cds_80);
					  $ni_cds_80 = sprintf("%.4f",$lni_cds);
					}
			  
				  # (3) CALCULATE DOS FOR THE ALIGNED CDS
				  if ($z == 0)
					{ if (($n_poly+$s_poly) == 0)
						{ $dos_cds_19 = ($n_subs/($n_subs+$s_subs)); }
					  else
						{ $dos_cds_19 = ($n_subs/($n_subs+$s_subs)) - ($n_poly/($n_poly+$s_poly)); } # look in Stoletzki and Eyre-Walker for justification as to why this is NOT calculated (unlike Haldane, below) using information per site
					  $dos_cds_19 = sprintf("%.4f",$dos_cds_19);
					}
				  elsif ($z == 1)
					{ if (($n_poly+$s_poly) == 0)
						{ $dos_cds_80 = ($n_subs/($n_subs+$s_subs)); }
					  else
						{ $dos_cds_80 = ($n_subs/($n_subs+$s_subs)) - ($n_poly/($n_poly+$s_poly)); } # look in Stoletzki and Eyre-Walker for justification as to why this is NOT calculated (unlike Haldane, below) using information per site
					  $dos_cds_80 = sprintf("%.4f",$dos_cds_80);
					}
				}
			}
		}
	
	  print OUT "$gene_name\t$gene_id\t$gene_desc\t$gene_type\t$loc\t$gene_len\t$no_of_transcripts_with_cds\t$transcript_id_with_longest_cds\t$length_of_longest_cds\t$exon_ids\t$exon_coords\t$total_exon_length\t$pc_of_transcript_that_is_exonic\t$cds_coords\t$pc_of_transcript_that_is_cds\t$go_f_line\t$go_p_line\t$go_c_line\t$a_lyrata_gene_id\t$a_lyrata_coords\t$a_lyrata_gene_len\t$no_of_lyrata_transcripts_with_cds\t$lyrata_transcript_id_with_longest_cds\t$length_of_longest_lyrata_cds\t$homology_type\t$pc_identity\t$dn\t$ds\t$dnds\t$is_high_confidence_gene_level\t$pc_identity_in_aln\t$pc_gaps_in_aln\t$number_of_aligned_codons\t$dn_cds\t$ds_cds\t$dnds_cds\t$is_high_confidence_cds_level\t$n_sites_cds\t$s_sites_cds\t$n_subs_cds\t$s_subs_cds\t$ni_cds_19\t$ni_cds_80\t$dos_cds_19\t$dos_cds_80\n";
	}
close(OUT) or die $!;
exit 1;

sub CONVERT_TO_COORDS
	{ my $param = shift;
	  my %hash = %$param;
	  my $coords = '';
	  my @bases = ();
	  while((my $base,my $irrel)=each(%hash))
		{ push(@bases,$base); }
	  my @sorted_bases = sort {$a <=> $b} @bases;
	  @bases = @sorted_bases;
	  my $start_base = '';
	  my %bases = ();
	  for(my $x=0;$x<@bases;$x++)
		{ my $base = $bases[$x];
		  if ($x == 0)
			{ $start_base = $base; }
		  push(@{$bases{$start_base}},$base);
		  if ((defined($bases[$x+1])) && ($bases[$x+1] != $base+1))
			{ push(@{$bases{$start_base}},$base);
			  $start_base = $bases[$x+1];
			}
		}
	  my @start_bases = ();
	  while((my $start_base,my $irrel)=each(%bases))
		{ push(@start_bases,$start_base); }
	  my @sorted_start_bases = sort {$a <=> $b} @start_bases;
	  foreach my $start_base (@sorted_start_bases)
		{ my @arr = @{$bases{$start_base}};
		  my @sorted_arr = sort {$a <=> $b} @arr;
		  my $end_base   = $sorted_arr[$#sorted_arr];
		  if ($start_base != $end_base)
			{ $coords .= "$start_base-$end_base,"; }
		  else
			{ $coords .= "$start_base,"; }
		}
	  $coords =~ s/\,$//;
	  return $coords;
	}
	
sub LOG10
	{ my $n = shift;
	  return log($n)/log(10);
	}
	
sub LOAD_MAIN_CODE
	{ open(MAIN_CODE,$main_code) or die "Unable to open $main_code\n";
      while ($_= <MAIN_CODE>)
		{ /^\s*([ACGTN]{3})\s+([A-Z\*\1-8])/ or die "Error in $main_code in line $_\n";
          $trans{$1} = $2;
        }
      close(MAIN_CODE) or die $!;
	}
	
sub LOAD_MITO_CODE
	{ open(MITO_CODE,$mito_code) or die "Unable to open $mito_code\n";
      while ($_= <MITO_CODE>)
		{ /^\s*([ACGTN]{3})\s+([A-Z\*\1-8])/ or die "Error in $mito_code at line $_\n";
          $mito_trans{$1} = $2;
        }
      close(MITO_CODE) or die $!;
	}