# PURPOSE: confirm that the sequence of each gene can indeed be aligned back to Col-0; if not, the gene is likely to be experiencing presence/absence variation in that accession

use strict;
use warnings;

# REQUIREMENTS
my $in_dir    = '1001_genomes_as_gene_seqs'; # from 4.parse_1001_genomes_part3_extract_gene_seqs_from_fa.pl
my $genome    = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
my $gene_info = 'Arabidopsis_thaliana.TAIR10.gene_info_and_go_terms.tsv'; # from Ensembl BioMart: Gene stable ID | Gene description | Chromosome/scaffold name | Gene start (bp) | Gene end (bp) | Strand | Gene name | Gene type | GO term accession | GO term name | GO term evidence code | GO domain

# OUTPUT
my $out_dir = '1001_genomes_as_gene_seqs_aligned_to_col0';
my $dir_of_sh_scripts = 'shell_scripts_for_1001_genomes_as_gene_seqs_aligned_to_col0';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
if (!(-d($dir_of_sh_scripts))) { mkdir $dir_of_sh_scripts or die $!; }
my $all_sh_commands = 'run_all_sh_scripts_for_needle.sh';
open(ALL_SH,'>',$all_sh_commands) or die $!;
print ALL_SH "#!/bin/bash\n";

# TEMPORARY
my $temp_dir = 'EMBOSS_alignment_temp_dir';
if (!(-d($temp_dir))) { mkdir $temp_dir or die $!; }

# STORE A. THALIANA GENOME
my $chr = ''; my %col_seq = ();
open(IN,$genome) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $chr = $1; }
	  next if ($line =~ /^\>/);
	  $col_seq{$chr} .= $line;
	}
close(IN) or die $!;

# STORE A. THALIANA GENE INFO
my %gene_coords = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $chr = $line[2]; my $gene_start = $line[3]; my $gene_end = $line[4]; my $strand = $line[5];
	  my $gene_len = ($gene_end-$gene_start)+1;
	  $gene_coords{$gene_id}{chr} 		 = $chr;
	  $gene_coords{$gene_id}{gene_start} = $gene_start;
	  $gene_coords{$gene_id}{gene_end}   = $gene_end;
	  $gene_coords{$gene_id}{gene_len}   = $gene_len;
	}
close(IN) or die $!;

# STORE GENE SEQUENCES FROM ALL ACCESSIONS
opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
my %seqs_per_gene_per_acc = ();
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "reading sequences: $files_seen of $files_total...\n";
	  my $acc = '';
	  if ($file =~ /^(.*?)\.fa$/) { $acc = $1; }
	  my $gene_id = '';
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  if ($line =~ /^\>(.*?)\|.*?$/)
			{ $gene_id = $1; }
		  else
			{ $seqs_per_gene_per_acc{$acc}{$gene_id} .= $line; }
		}
	  close(IN) or die $!;
	}
	
# ALIGN ALL GENE SEQUENCES BACK TO COL-0 TO DETERMINE IF ANY ARE MISSING (THAT IS, PRESENCE/ABSENCE VARIATION)
my @accs = ();
while((my $acc,my $irrel)=each(%seqs_per_gene_per_acc))
	{ push(@accs,$acc); }
my @sorted_accs = sort {$a <=> $b} @accs;
my $accs_total = @sorted_accs;
for(my $x=0;$x<@sorted_accs;$x++)
	{ my $accs_seen = $x+1;
	  my $acc = $sorted_accs[$x];
	  my $out_file = "$out_dir/$acc.txt";
	  next if (-e($out_file)); # CHECKPOINT: we've seen this before
	  my $sh_file = "$dir_of_sh_scripts/$acc.sh";
	  open(SH,'>',$sh_file) or die $!;
	  print SH "#!/bin/bash\n";
	  print SH "echo -e \"Gene ID\\t% identity relative to the sequence of this gene in Col-0\" >> $out_file\n";
	  if (!(-d("$temp_dir/$acc"))) { mkdir "$temp_dir/$acc" or die $!; }
	  my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%{$seqs_per_gene_per_acc{$acc}}))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  my $genes_seen = 0; my $genes_total = @sorted_gene_ids;
	  foreach my $gene_id (@sorted_gene_ids)
		{ $genes_seen++;
		  print "creating needle files: $accs_seen of $accs_total, $genes_seen of $genes_total...\n";
		  my $tair10_chr 		= $gene_coords{$gene_id}{chr};
		  my $tair10_gene_start = $gene_coords{$gene_id}{gene_start};
		  my $tair10_gene_len   = $gene_coords{$gene_id}{gene_len};
		  my $tair10_gene_seq   = substr($col_seq{$tair10_chr},$tair10_gene_start-1,$tair10_gene_len);
	      my $acc_seq 			= $seqs_per_gene_per_acc{$acc}{$gene_id};	      
		  my $temp_file1 = "$temp_dir/$acc/$gene_id.thaliana.fa";
		  my $temp_file2 = "$temp_dir/$acc/$gene_id.acc.fa";
		  my $temp_file3 = "$temp_dir/$acc/$gene_id.thaliana.acc.txt";
		  open (TEMP,'>',$temp_file1) or die $!;
		  print TEMP "$tair10_gene_seq\n";
		  close(TEMP) or die $!;
		  open (TEMP,'>',$temp_file2) or die $!;
		  print TEMP "$acc_seq\n";
		  close(TEMP) or die $!;
		  print SH "needle $temp_file1 $temp_file2 -gapopen 10.0 -gapextend 0.5 -outfile $temp_file3\n";
		  print SH "rm $temp_file1 $temp_file2\n";
		  print SH "pc_identity=\$(grep \"Identity:\" $temp_file3 | awk '{ print \$4 }' | perl -pe 's/[()%]//g' )\n";
		  print SH "echo -e \"$gene_id\\t\$pc_identity\" >> $out_file\n";
		  print SH "rm $temp_file3\n";
		}
	  print SH "rm -r $temp_dir/$acc\n";
	  close(SH) or die $!;
	  print ALL_SH "$sh_file\n";
	}
close(ALL_SH) or die $!;
exit 1;