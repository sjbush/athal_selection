use strict;
use warnings;

# REQUIREMENTS
my $in_dir = 'genes_from_1001_genomes_and_lyrata_for_msa'; # from 7.parse_1001_genomes_part6_store_all_acc_seqs_per_gene.pl

# OUTPUT
my $out_dir = 'individual_genes_from_1001_genomes_and_lyrata_for_msa';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my $files_seen = 0; my $files_total = @files; $files_total = $files_total-2;
foreach my $file (@files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  my $acc = '';
	  if ($file =~ /^(.*?)\.fa$/)
		{ $acc = $1; }
	  my $gene_id = '';
	  my %thaliana_seqs = (); my %lyrata_seqs = ();
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  if ($line =~ /^\>(.*?)$/)
			{ $gene_id = $1; }
		  next if ($line =~ /^\>/);
		  if ($gene_id =~ /^AT\d+G\d+$/)
			{ $thaliana_seqs{$gene_id} .= $line; }
		  elsif ($gene_id =~ /^(.*?)\-Lyrata$/)
			{ $lyrata_seqs{$1} .= $line; }
		}
	  close(IN) or die $!;
	  my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%thaliana_seqs))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  my $genes_seen = 0; my $genes_total = @sorted_gene_ids;
	  foreach my $gene_id (@sorted_gene_ids)
		{ $genes_seen++;
		  print "$genes_seen of $genes_total for file $files_seen of $files_total\n";
		  if (!(-d("$out_dir/$gene_id"))) { mkdir "$out_dir/$gene_id" or die $!; }
		  if (!(-e("$out_dir/$gene_id/$acc.fa")))
			{ open(OUT,'>',"$out_dir/$gene_id/$acc.fa") or die $!;
			  print OUT ">$acc\n$thaliana_seqs{$gene_id}\n";
			  close(OUT) or die $!;
			}
		  if (!(-e("$out_dir/$gene_id/Lyrata.fa")))
			{ open(OUT,'>',"$out_dir/$gene_id/Lyrata.fa") or die $!;
			  print OUT ">Lyrata\n$lyrata_seqs{$gene_id}\n";
			  close(OUT) or die $!;
			}
		}
	}

exit 1;