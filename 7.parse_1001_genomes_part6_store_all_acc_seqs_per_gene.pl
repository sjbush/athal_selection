use strict;
use warnings;

# REQUIREMENTS
my $in_dir1		  = '1001_genomes_as_gene_seqs_aligned_to_col0'; # from 6.parse_1001_genomes_part5_confirm_gene_seq_aligns_back_to_col0.pl
my $in_dir2   	  = '1001_genomes_as_CDS_seqs'; # from 5.parse_1001_genomes_part4_obtain_cds_seq_per_gene.pl
my $gene_info  	  = 'Arabidopsis_thaliana.TAIR10.gene_info_and_go_terms.tsv'; # from Ensembl BioMart: Gene stable ID | Gene description | Chromosome/scaffold name | Gene start (bp) | Gene end (bp) | Strand | Gene name | Gene type | GO term accession | GO term name | GO term evidence code | GO domain
my $orthology  	  = 'Arabidopsis_thaliana.TAIR10.orthology_relationships_and_dnds.tsv'; # from Ensembl BioMart: Gene stable ID | Arabidopsis lyrata gene stable ID | Arabidopsis lyrata chromosome/scaffold name | Arabidopsis lyrata chromosome/scaffold start (bp) | Arabidopsis lyrata chromosome/scaffold end (bp) | Arabidopsis lyrata homology type | %id. target Arabidopsis lyrata gene identical to query gene | %id. query gene identical to target Arabidopsis lyrata gene | Arabidopsis lyrata Whole-genome alignment coverage | dN with Arabidopsis lyrata | dS with Arabidopsis lyrata | Arabidopsis lyrata orthology confidence [0 low, 1 high]
my $lyrata_cds 	  = 'Arabidopsis_lyrata.v.1.0.cds.all.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-38/plants/fasta/arabidopsis_lyrata/cds/Arabidopsis_lyrata.v.1.0.cds.all.fa.gz
my $lyrata_genome = 'Arabidopsis_lyrata.v.1.0.dna.toplevel.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/arabidopsis_lyrata/dna/Arabidopsis_lyrata.v.1.0.dna.toplevel.fa.gz

# TEMPORARY
my $temp_dir = 'thal_to_lyr_alignment_temp_dir';
if (!(-d($temp_dir))) { mkdir $temp_dir or die $!; }

# PARAMETERS
my $min_pc_identity = 75;

# OUTPUT
my $out_dir = 'genes_from_1001_genomes_and_lyrata_for_msa';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_file = 'coordinates_of_genes_in_each_accession.tsv';
open(OUT,'>',$out_file) or die $!;
print OUT "Gene ID\tAccession\tGene coordinates in accession\t% identity of accession CDS to the CDS of Col-0 (after pairwise alignment with EMBOSS needle)\t% identity of accession CDS to the CDS of A. lyrata (after pairwise alignment with EMBOSS needle)\n";

# STORE A. LYRATA GENOME
my $header = ''; my %genome = ();
open(IN,$lyrata_genome) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $header = $1;	}
	  next if ($line =~ /^\>/);
	  $genome{$header} .= $line;
	}
close(IN) or die $!;

# STORE A. LYRATA CDS
my %lyrata_cds = ();
my $lyrata_transcript_id = '';
open(IN,$lyrata_cds) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) cds .*?\:.*?\:.*?:\d+:\d+:.*? gene\:.*? .*?$/) # a somewhat convoluted regex because most of these headers read "cds chromosome|scaffold:v.1.0" but NOT all
		{ $lyrata_transcript_id = $1; }
	  next if ($line =~ /^\>/);
	  $lyrata_cds{$lyrata_transcript_id} .= $line;
	}
close(IN) or die $!;

# STORE THE STRANDS OF EACH A. THALIANA GENE
my %thaliana_strand = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $strand = $line[5];
	  $thaliana_strand{$gene_id} = $strand;
	}
close(IN) or die $!;

# STORE COORDINATES OF A. LYRATA GENES, BUT ONLY IF ENSEMBL PROVIDES A 'HIGH-CONFIDENCE' WHOLE-GENE dN/dS THALIANA/LYRATA ESTIMATE.
# higher-confidence criteria: the dN/dS estimate is derived using a one-to-one orthologue with >=75% reciprocal identity, and has dS < 2, dS > 0.02, and dN < 2
my %lyrata_orthologues = ();
open(IN,$orthology) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  next if (!(defined($line[1])));
	  my $a_lyrata_gene_id  = $line[1];
	  my $homology_type     = $line[5];
	  my $lyr_to_thal_pc    = sprintf("%.1f",$line[6]);
	  my $thal_to_lyr_pc    = sprintf("%.1f",$line[7]);
	  my $dn = 'NA'; my $ds = 'NA';
	  if (defined($line[9]))  { $dn = sprintf("%.4f",$line[9]);  }
	  if (defined($line[10])) { $ds = sprintf("%.4f",$line[10]); }
	  my $dnds = 'NA';
	  if (($dn =~ /\d/) && ($ds =~ /\d/) && ($ds > 0)) { $dnds = sprintf("%.4f",($dn/$ds)); }	  
	  if (($dnds =~ /\d/) and ($ds > 0.02) and ($ds < 2) and ($dn < 2) and ($homology_type eq 'ortholog_one2one') and ($thal_to_lyr_pc >= 75) and ($lyr_to_thal_pc >= 75))
		{ $lyrata_orthologues{$gene_id} = $a_lyrata_gene_id;
		}
	}
close(IN) or die $!;

# STORE % IDENTITY OF EACH GENE TO COL-0 (WE WILL USE THIS TO FILTER OUT GENES THAT AREN'T REALLY PRESENT IN A GIVEN ACCESSION)
opendir(DIR,$in_dir1) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
my %pc_identities = ();
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "storing % identity of each gene: $files_seen of $files_total...\n";
	  my $acc = '';
	  if ($file =~ /^(.*?)\.txt$/) { $acc = $1; }
	  open(IN,"$in_dir1/$file") or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $gene_id = $line[0]; my $pc_identity = $line[1];
		  $pc_identities{$gene_id}{$acc} = $pc_identity;
		}
	  close(IN) or die $!;
	}

# STORE GENE SEQUENCES, PROVIDING THEY ARE ABOVE A THRESHOLD % IDENTITY WITH COL-0
opendir(DIR,$in_dir2) or die $!;
@files = readdir(DIR);
closedir(DIR) or die $!;
@sorted_files = sort {$a cmp $b} @files;
$files_seen = 0; $files_total = @sorted_files; $files_total = $files_total-2;
my %coords_per_gene_per_acc = (); my %seqs_per_gene_per_acc = ();
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "reading sequences: $files_seen of $files_total...\n";
	  my $acc = '';
	  if ($file =~ /^(.*?)\.fa$/) { $acc = $1; }
	  my $gene_id = '';
	  open(IN,"$in_dir2/$file") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  if ($line =~ /^\>(.*?)\|(.*?)$/)
			{ $gene_id = $1; my $coords = $2;
			  $coords_per_gene_per_acc{$gene_id}{$acc} = $coords;
			}
		  else
			{ $seqs_per_gene_per_acc{$gene_id}{$acc} .= $line;
			}
		}
	  close(IN) or die $!;
	}
	
# OUTPUT DATA IN A FORMAT SUITABLE FOR MULTIPLE SEQUENCE ALIGNMENT WITH MAFFT
my @gene_ids = ();
while((my $gene_id,my $irrel)=each(%seqs_per_gene_per_acc))
	{ push(@gene_ids,$gene_id); }
my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
my $genes_seen = 0; my $genes_total = @sorted_gene_ids;
foreach my $gene_id (@sorted_gene_ids)
	{ $genes_seen++;
	  my @accs = ();
	  while((my $acc,my $irrel)=each(%{$seqs_per_gene_per_acc{$gene_id}}))
		{ push(@accs,$acc); }
	  my @sorted_accs = sort {$a cmp $b} @accs;
	  
	  next if (!(exists($lyrata_orthologues{$gene_id})));
	  
	  # A. lyrata CDS sequence
	  my $lyrata_gene_id = $lyrata_orthologues{$gene_id};
	  my $lyrata_seq = $lyrata_cds{$lyrata_gene_id};
	  
	  # test whether the A. lyrata (CDS) sequence aligns meaningfully to the set of A. thaliana (CDS) sequences - we already know it does in Col-0 but this might not be the case for all of them
	  my @possible_thaliana_seqs = ();
	  foreach my $acc (@sorted_accs)
		{ my $seq = $seqs_per_gene_per_acc{$gene_id}{$acc};
		  if (exists($pc_identities{$gene_id}{$acc}))
			{ my $pc_identity = $pc_identities{$gene_id}{$acc};
			  next if ($pc_identity < $min_pc_identity); # FILTER: this sequence does not share enough identity with this gene in Col-0; this might be because it is not really the right gene
			  push(@possible_thaliana_seqs,[$seq,$acc]);
			}
		}
	  my $temp_file1 = "$temp_dir/lyrata.$gene_id.fa";
	  open (TEMP,'>',$temp_file1) or die $!;
	  print TEMP "$lyrata_seq\n";
	  close(TEMP) or die $!;
	  my %trimmed_thalianas = (); my %pc_identities_relative_to_lyr = ();
	  my $number_of_thalianas_seen = 0; my $number_of_thalianas_total = @possible_thaliana_seqs;
	  my $there_is_a_usable_lyrata = 0;
	  for(my $x=0;$x<@possible_thaliana_seqs;$x++)
		{ $number_of_thalianas_seen++;
		  my $seq = $possible_thaliana_seqs[$x][0]; my $acc = $possible_thaliana_seqs[$x][1];
		  next if (!(defined($thaliana_strand{$gene_id})));
		  my $strand = $thaliana_strand{$gene_id};
		  if ($strand == -1)
			{ $seq =~ tr/[ATCG]/[TAGC]/;
			  $seq = reverse($seq);
			}
		  my $temp_file2 = "$temp_dir/thaliana.$gene_id.fa";
		  my $temp_file3 = "$temp_dir/$gene_id.alignment";
		  open (TEMP,'>',$temp_file2) or die $!;
		  print TEMP "$seq\n";
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
		  unlink $temp_file2 or die $!;
		  my $type = 'lyrata'; my %alns = ();
		  my $pc_identity_to_lyr = 0;
		  open(IN,$temp_file3) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  if ($line =~ /^.*?Identity\:.*?\((.*?)%\)$/)
				{ $pc_identity_to_lyr = $1; }
			  if ($line =~ /^\s+\d/)
				{ $alns{$type} .= $line;
				  if 	($type eq 'thaliana') { $type = 'lyrata';   }
				  elsif ($type eq 'lyrata')   { $type = 'thaliana'; }
				}
			}
		  close(IN) or die $!;
		  unlink $temp_file3 or die $!;
		  next if ($pc_identity_to_lyr < $min_pc_identity);

		  $pc_identities_relative_to_lyr{$gene_id}{$acc} = $pc_identity_to_lyr;
		  print "output: genes $genes_seen of $genes_total... A. thaliana accession $number_of_thalianas_seen of $number_of_thalianas_total\n";
		  $there_is_a_usable_lyrata++;
		}
	  unlink $temp_file1 or die $!;
	  
	  next if ($there_is_a_usable_lyrata == 0);
	  
	  if (!(-d("$out_dir/$gene_id"))) { mkdir "$out_dir/$gene_id" or die $!; }	  
	  
	  open(FA,'>',"$out_dir/$gene_id/Lyrata.fa") or die $!;
	  print FA ">Lyrata\n$lyrata_seq\n";
	  close(FA) or die $!;
	  
	  foreach my $acc (@sorted_accs)
		{ if (exists($pc_identities{$gene_id}{$acc}))
			{ my $seq    			 = $seqs_per_gene_per_acc{$gene_id}{$acc};
			  my $coords 			 = $coords_per_gene_per_acc{$gene_id}{$acc};
			  my $pc_identity 		 = $pc_identities{$gene_id}{$acc};
			  my $pc_identity_to_lyr = $pc_identities_relative_to_lyr{$gene_id}{$acc};
			  next if ($pc_identity < $min_pc_identity); # FILTER: this sequence does not share enough identity with this gene in Col-0; this might be because it is not really the right gene
			  next if (!(defined($thaliana_strand{$gene_id})));
			  my $strand = $thaliana_strand{$gene_id};
			  if ($strand == -1)
				{ $seq =~ tr/[ATCG]/[TAGC]/;
				  $seq = reverse($seq);
				}
			  print OUT "$gene_id\t$acc\t$coords\t$pc_identity\t$pc_identity_to_lyr\n";
			  open(FA,'>',"$out_dir/$gene_id/$acc.fa") or die $!;
			  print FA ">$acc\n$seq\n";
			  close(FA) or die $!;
			}
		}
	}

close(OUT) or die $!;
exit 1;