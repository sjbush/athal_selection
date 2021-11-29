use strict;
use warnings;

# REQUIREMENTS
my $in_dir    = '1001_genomes_as_fa'; # from 2.parse_1001_genomes_part1_make_alignments_from_vcf.pl
my $gene_info = 'Arabidopsis_thaliana.TAIR10.gene_info_and_go_terms.tsv'; # from Ensembl BioMart: Gene stable ID | Gene description | Chromosome/scaffold name | Gene start (bp) | Gene end (bp) | Strand | Gene name | Gene type | GO term accession | GO term name | GO term evidence code | GO domain

# OUTPUT
my $out_dir = '1001_genomes_as_gene_seqs';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

# STORE COORDINATES OF A. THALIANA GENES
my %gene_coords = (); my %tair10_gene_starts = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $chr = $line[2]; my $gene_start = $line[3]; my $gene_end = $line[4]; my $strand = $line[5];
	  $gene_coords{$chr}{$gene_id} = "$chr:$gene_start-$gene_end:$strand";
	  $tair10_gene_starts{$chr}{$gene_start} = $gene_id;
	}
close(IN) or die $!;

# STORE THE TAIR10-TO-ACCESSION ALIGNMENTS. USING THE TAIR10 GENE COORDINATES, WE WILL EXTRACT THE CORRESPONDING SEQUENCE FROM EACH ACCESSION.
opendir(DIR,$in_dir) or die $!;
my @subdirs = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_subdirs = sort {$a cmp $b} @subdirs;
my $subdirs_seen = 0; my $subdirs_total = @sorted_subdirs; $subdirs_total = $subdirs_total-2;
foreach my $subdir (@sorted_subdirs)
	{ next if (($subdir eq '.') or ($subdir eq '..'));
	  $subdirs_seen++;
	  next if (!(-d("$in_dir/$subdir")));
	  next if (-e("$out_dir/$subdir.fa")); # CHECKPOINT: we've seen this already
	  
	  # store all whole-chromosome accession-to-TAIR10 alignments, so that the gene sequences from each accession can be extracted (by reference to the gene location in TAIR10)
	  opendir(DIR,"$in_dir/$subdir") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  my @sorted_files = sort {$a cmp $b} @files;
	  my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
	  my %seqs = ();
	  foreach my $file (@sorted_files)
		{ next if (($file eq '.') or ($file eq '..'));
		  $files_seen++;
		  print "storing data for accession $subdirs_seen of $subdirs_total, chrs $files_seen of $files_total...\n";
		  my $acc = ''; my $chr = '';
		  if ($file =~ /^(.*?)\-(.*?)\.fa$/) { $acc = $1; $chr = $2; }
		  my $tair10_or_acc = ''; my $start_pos = '';
		  open(IN,"$in_dir/$subdir/$file") or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  if ($line =~ /^\>(.*?) \| (\d+)$/)
				{ $tair10_or_acc = $1; $start_pos = $2; }
			  next if ($line =~ /^\>/);
			  $seqs{$chr}{$tair10_or_acc}{seq} .= $line;
			  $seqs{$chr}{$tair10_or_acc}{pos}  = $start_pos;
			}
		  close(IN) or die $!;
		}
	  
	  # store the gene sequences from each accession
	  my %gene_seqs_per_acc = ();
	  my $chrs_seen = 0; my $chrs_total = scalar keys %seqs;
	  while((my $chr,my $irrel)=each(%seqs))
		{ $chrs_seen++;
		  my $tair10_chr_seq  		 = $seqs{$chr}{TAIR10}{seq}; # don't forget this has gaps in!
		  my $acc_chr_seq     		 = $seqs{$chr}{$subdir}{seq};
		  my $tair10_alignment_start = $seqs{$chr}{TAIR10}{pos};
		  my $acc_alignment_start    = $seqs{$chr}{$subdir}{pos};
		  my @tair10_chr_seq 		 = split(//,$tair10_chr_seq);
		  my @acc_chr_seq    	     = split(//,$acc_chr_seq);
		  
		  # identify the start positions, in the array, of each TAIR10 gene
		  my %actual_gene_starts = ();
		  my $actual_pos_in_tair10 = $tair10_alignment_start-1; my $actual_pos_in_acc = $acc_alignment_start-1; # or, for some reason, $actual_pos_in_acc = $acc_alignment_start-2
		  for(my $x=0;$x<@tair10_chr_seq;$x++)
			{ my $tair10_base = $tair10_chr_seq[$x];
			  my $acc_base    = $acc_chr_seq[$x];
			  if ($tair10_base ne '-')
				{ $actual_pos_in_tair10++; }
			  if ($acc_base ne '-')
				{ $actual_pos_in_acc++; }
			  if (exists($tair10_gene_starts{$chr}{$actual_pos_in_tair10}))
				{ my $gene_id = $tair10_gene_starts{$chr}{$actual_pos_in_tair10};
				  $actual_gene_starts{$gene_id}{array_location} 	  = $x;
				  $actual_gene_starts{$gene_id}{actual_pos_in_acc}    = $actual_pos_in_acc;
				  $actual_gene_starts{$gene_id}{actual_pos_in_tair10} = $actual_pos_in_tair10;
				}
			}
		  
		  # extract the TAIR10 sequence for each gene, and then - by reference to the alignment - the corresponding accession sequence
		  my $genes_seen = 0; my $genes_total = scalar keys %{$gene_coords{$chr}};
		  while((my $gene_id,my $irrel)=each(%{$gene_coords{$chr}}))
			{ $genes_seen++;
			  print "accession $subdirs_seen of $subdirs_total. extracting data from chrs $chrs_seen of $chrs_total, genes $genes_seen of $genes_total\n";
			  my $tair10_gene_coords = $gene_coords{$chr}{$gene_id};
			  if ($tair10_gene_coords =~ /^(.*?)\:(\d+)\-(\d+)\:(.*?)$/)
				{ my $this_chr = $1; my $gene_start = $2; my $gene_end = $3; my $strand = $4;
				  my $gene_length = ($gene_end-$gene_start)+1;
				  next if (!(exists($actual_gene_starts{$gene_id}{array_location})));
				  my $search_start 		   = $actual_gene_starts{$gene_id}{array_location}; # to save time, we're going to start looking for matching sequences shortly before the start of the gene of interest
				  my $actual_pos_in_tair10 = $actual_gene_starts{$gene_id}{actual_pos_in_tair10};
				  my $actual_pos_in_acc    = $actual_gene_starts{$gene_id}{actual_pos_in_acc};
				  my $start_pos_in_acc;
				  my $tair10_gene_seq = ''; my $acc_gene_seq = '';
				  for(my $x=$search_start;$x<@tair10_chr_seq;$x++)
					{ my $tair10_base = $tair10_chr_seq[$x];
					  my $acc_base    = $acc_chr_seq[$x];
					  if ($tair10_base ne '-')
						{ $actual_pos_in_tair10++; }
					  if ($acc_base ne '-')
						{ $actual_pos_in_acc++; }
					  if ($actual_pos_in_tair10 >= $gene_start)
						{ $tair10_gene_seq .= $tair10_base unless ($tair10_base eq '-');
						  $acc_gene_seq    .= $acc_base    unless ($acc_base    eq '-');
						  if ($acc_base ne '-')
							{ $start_pos_in_acc = $actual_pos_in_acc unless (defined($start_pos_in_acc));
							}
						}
					  my $current_length = length($tair10_gene_seq);
					  last if ($current_length == $gene_length);
					}
				  my $acc_gene_coords = "$chr:$start_pos_in_acc-$actual_pos_in_acc:$strand";
				  $gene_seqs_per_acc{$gene_id}{seq}    = $acc_gene_seq;
				  $gene_seqs_per_acc{$gene_id}{coords} = $acc_gene_coords;
				}
			}
		}
	  my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%gene_seqs_per_acc))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  my $out_file = "$out_dir/$subdir.fa";
	  open(OUT,'>',$out_file) or die $!;
	  foreach my $gene_id (@sorted_gene_ids)
		{ my $acc_gene_seq    = $gene_seqs_per_acc{$gene_id}{seq};
		  my $acc_gene_coords = $gene_seqs_per_acc{$gene_id}{coords};
		  print OUT ">$gene_id|$acc_gene_coords\n$acc_gene_seq\n";
		}
	  close(OUT) or die $!;
	}

exit 1;