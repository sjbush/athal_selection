use strict;
use warnings;

# REQUIREMENTS
my $in_dir 			= '1001_genomes'; # wget --recursive --no-parent https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/
my $ref_genome		= 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
my $ref_genome_edit = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.edited.fa'; # created when this script is run for the first time; see below
my $chr1 			= '1.fa'; # the following files are manually created from the multi-fasta $ref_genome - one chr per fasta
my $chr2 			= '2.fa';
my $chr3 			= '3.fa';
my $chr4 			= '4.fa';
my $chr5 			= '5.fa';
my $chrM 			= 'Mt.fa';
my $chrP 			= 'Pt.fa';

# PARAMETERS
my $num_threads = 10;
my @chrs = (qw/1 2 3 4 5 Mt Pt/);

# OUTPUT
my $out_dir = '1001_genomes_alignments';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $sh_file = 'convert_1001_vcf_to_fa.sh';
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";
if (!(-e($ref_genome_edit)))
	{ print SH "sed '/^>/ s/ .*//' $ref_genome > $ref_genome_edit\n";
	  print SH "cat $ref_genome_edit | awk '{ if (substr(\$0, 1, 1)==\">\") {filename=(substr(\$0,2) \".fa\")} print \$0 > filename }'\n"; # split the multi-fasta into a set of individual chromosomes
	}
print SH "cd $out_dir\n";

opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  last if ($files_seen > 10);
	  print "$files_seen of $files_total...\n";
	  my $acc = '';
	  if ($file =~ /^(.*?)\.vcf\.gz$/)
		{ $acc = $1;
		  if ($acc =~ /^(.*?)\_snp\_short\_indel\_with\_quality\_reference$/)
			{ $acc = $1; }
		}
	  next if ($file !~ /^.*?\.vcf\.gz$/);
	  
	  # create an output sub-directory to store all chromosome-specific accession-to-TAIR10 alignments
	  my $seen_already = 0;
	  foreach my $chr (@chrs)
		{ if (-e("$out_dir/$acc/$acc-$chr.TAIR10.alignment"))
			{ $seen_already++; }
		}
	  next if ($seen_already > 0); # CHECKPOINT: we've seen it all already
	  if (!(-d("$out_dir/$acc")))
		{ print SH "mkdir $out_dir/$acc\n"; }
	  print SH "cd $out_dir/$acc\n";
	  print SH "tabix -p vcf $in_dir/$file\n" unless (-e("$in_dir/$file.tbi")); # index the VCF for later use with VCFtools
	  
	  # convert accession VCF to fasta
	  print SH "cat $ref_genome | vcf-consensus $in_dir/$file > $out_dir/$acc.fa\n"; # apply the variants in each VCF to the TAIR10 reference genome, creating a consensus multi-fasta
	  print SH "sed '/^>/ s/ .*//' $out_dir/$acc.fa > $out_dir/$acc.edited.fa\n"; # remove everything after the first space in each multi-fasta header, thus retaining only the chromosome name, i.e. editing "1 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF" to "1"
	  print SH "mv $out_dir/$acc.edited.fa $out_dir/$acc.fa\n";
	  print SH "cat $out_dir/$acc.fa | awk '{ if (substr(\$0, 1, 1)==\">\") {filename=(\"$acc-\" substr(\$0,2) \".fa\")} print \$0 > filename }'\n"; # split the multi-fasta into a set of individual chromosomes
	  print SH "rm $out_dir/$acc.fa\n";
	  
	  foreach my $chr (@chrs)
		{ # align each accession chromosome to its corresponding TAIR10 chromosome
		  if (!(-e("$out_dir/$acc/$acc-$chr-TAIR10.alignment")))
			{ print SH "nucmer -t $num_threads -p $acc-$chr-TAIR10 $chr.fa $out_dir/$acc/$acc-$chr.fa\n"; # perform a whole genome alignment of each A. thaliana accession chromosome to its corresponding TAIR10 chromosome; output stored in a .delta file
			  print SH "dnadiff -p $acc-$chr-TAIR10 -d $out_dir/$acc/$acc-$chr-TAIR10.delta\n"; # to filter the delta file to retain only 1-to-1 alignments
			  print SH "rm $out_dir/$acc/$acc-$chr-TAIR10.delta\n";
			  print SH "show-aligns -r $out_dir/$acc/$acc-$chr-TAIR10.1delta $chr $chr > $out_dir/$acc/$acc-$chr-TAIR10.alignment\n"; # extract 1-to-1 alignments, reporting locations relative to the reference genome (TAIR10)
			  print SH "rm $out_dir/$acc/$acc-$chr-TAIR10.1coords $out_dir/$acc/$acc-$chr-TAIR10.1delta $out_dir/$acc/$acc-$chr-TAIR10.mcoords $out_dir/$acc/$acc-$chr-TAIR10.mdelta $out_dir/$acc/$acc-$chr-TAIR10.qdiff $out_dir/$acc/$acc-$chr-TAIR10.rdiff $out_dir/$acc/$acc-$chr-TAIR10.snps $out_dir/$acc/$acc-$chr-TAIR10.report $out_dir/$acc/$acc-$chr-TAIR10.unref $out_dir/$acc/$acc-$chr-TAIR10.unqry\n";
			}
		  print SH "rm $out_dir/$acc/$acc-$chr.fa\n";
		}
	}
close(SH) or die $!;
exit 1;