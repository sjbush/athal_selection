use strict;
use warnings;

# REQUIREMENTS
my $in_dir = 'individual_genes_from_1001_genomes_and_lyrata_for_msa'; # from 8.parse_1001_genomes_part7_make_one_file_per_gene_for_msa.pl

# PARAMETERS
my $num_threads = 10;
	
# OUTPUT
my $out_dir = 'msa_of_genes_in_1001_genomes_and_lyrata';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $sh_file = 'run_msa_on_thaliana_and_lyrata_genes.sh';
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";
print SH "cd $out_dir\n";

opendir(DIR,$in_dir) or die $!;
my @genes = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_genes = sort {$a cmp $b} @genes;
my $genes_seen = 0; my $genes_total = @sorted_genes; $genes_total = $genes_total-2;
foreach my $gene_id (@sorted_genes)
	{ next if (($gene_id eq '.') or ($gene_id eq '..'));
	  next if (-e("$out_dir/$gene_id.fa")); # FILTER: we've seen this before
	  $genes_seen++;
	  print "$genes_seen of $genes_total\n";
	  if (-d("$in_dir/$gene_id"))
		{ opendir(DIR,"$in_dir/$gene_id") or die $!;
		  my @files = readdir(DIR);
		  closedir(DIR) or die $!;
		  my $file_line = '';
		  foreach my $file (@files)
			{ next if (($file eq '.') or ($file eq '..'));
			  $file_line .= "$in_dir/$gene_id/$file ";
			}
		  $file_line =~ s/ $//;

		  # make multiple alignments with MAFFT
		  if (!(-e("$in_dir/$gene_id.fa"))) { print SH "cat $file_line > $in_dir/$gene_id.fa\n"; }
		  print SH "mafft --thread $num_threads $in_dir/$gene_id.fa > $out_dir/$gene_id.fa\n";
		}
	}
close(SH) or die $!;
exit 1;