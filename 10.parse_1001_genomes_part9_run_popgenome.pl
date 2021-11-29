use strict;
use warnings;

# REQUIREMENTS
my $in_dir = 'msa_of_genes_in_1001_genomes_and_lyrata'; # from 9.parse_1001_genomes_part8_msa_all_thaliana_and_lyrata_genes.pl

# OUTPUT
my $out_dir = 'msa_of_genes_in_1001_genomes_and_lyrata_in_fa_format';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_file = 'popgenome_stats_for_genes_in_1001_genomes_and_lyrata.tsv';
my $commands_for_popgenome = 'run_popgenome_on_msa_of_genes_in_1001_genomes_and_lyrata.R';
open(R,'>',$commands_for_popgenome) or die $!;
print R "library(PopGenome)\n";
print R "out_df<-data.frame()\n";

opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "$files_seen of $files_total\n";
	  my $gene_id = '';
	  if ($file =~ /^(.*?)\.fa$/) { $gene_id = $1; }
	  my $in_file = "$in_dir/$gene_id.fa";
	  next if (!(-e($in_file)));
	  if (!(-d("$out_dir/$gene_id"))) { mkdir "$out_dir/$gene_id" or die $!; }
	  open(OUT,'>',"$out_dir/$gene_id/$gene_id.fa") or die $!;
	  open(IN,$in_file) or die $!;
	  while(<IN>) { print OUT "$_"; }
	  close(IN) or die $!;
	  close(OUT) or die $!;
	  print R "GENOME.class <- readData(\"$out_dir/$gene_id\")\n";
	  print R "GENOME.class <-set.outgroup(GENOME.class,c(\"Lyrata\"))\n";
	  print R "GENOME.class\@region.data\@outgroup\n";
	  print R "GENOME.class <- neutrality.stats(GENOME.class, detail=TRUE)\n";
	  print R "get.neutrality(GENOME.class)[[1]]\n";
	  print R "tajima_d<-GENOME.class\@Tajima.D[1]\n";
	  print R "fay_wu_h<-GENOME.class\@Fay.Wu.H[1]\n";
	  print R "fu_li_d<-GENOME.class\@Fu.Li.D[1]\n";
	  print R "zeng_e<-GENOME.class\@Zeng.E[1]\n";
	  print R "no_of_segregating_sites<-GENOME.class\@n.segregating.sites[1]\n";
	  print R "df.new<-data.frame(GENE_ID='$gene_id',TAJIMA_D=tajima_d,FAY_WU_H=fay_wu_h,FU_LI_D=fu_li_d,ZENG_E=zeng_e,NO_OF_SEGREGATING_SITES=no_of_segregating_sites)\n";
	  print R "out_df<-rbind(out_df,df.new)\n";
	}
print R "write.table(out_df,'$out_file',quote=FALSE,row.names=FALSE,sep='\\t')\n";
close(R) or die $!;
exit 1;