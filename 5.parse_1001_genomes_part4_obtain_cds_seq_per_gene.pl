use strict;
use warnings;

# REQUIREMENTS
my $in_dir   	= '1001_genomes_as_gene_seqs'; # from 4.parse_1001_genomes_part3_extract_gene_seqs_from_fa.pl
my $coords_file = 'Arabidopsis_thaliana.TAIR10.cds_coords.tsv'; # from 1.align_TAIR10_CDS_coords.pl

# OUTPUT
my $out_dir = '1001_genomes_as_CDS_seqs';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

# STORE COORDINATES OF A. THALIANA GENE AND CDS SEQUENCES
my %coords = ();
open(IN,$coords_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[2]; my $tair10_gene_coords = $line[3]; my $tair10_cds_coords = $line[9];
	  if ($tair10_gene_coords =~ /^(.*?)\:(\d+)\-(\d+)\:(.*?)$/)
		{ my $gene_start = $2;
		  $coords{$gene_id}{cds_coords} = $tair10_cds_coords;
		  $coords{$gene_id}{gene_start} = $gene_start;
		}
	}
close(IN) or die $!;

opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
my %seqs_per_gene_per_acc = (); my %acc_gene_coords = ();
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "reading sequences: $files_seen of $files_total...\n";
	  my $acc = '';
	  if ($file =~ /^(.*?)\.fa$/) { $acc = $1; }
	  my $gene_id = ''; my $acc_gene_coords = '';
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  if ($line =~ /^\>(.*?)\|(.*?)$/)
			{ $gene_id = $1; $acc_gene_coords = $2;
			  $acc_gene_coords{$gene_id}{$acc} = $acc_gene_coords;
			}
		  else
			{ $seqs_per_gene_per_acc{$acc}{$gene_id} .= $line; }
		}
	  close(IN) or die $!;
	}
	
my @accs = ();
while((my $acc,my $irrel)=each(%seqs_per_gene_per_acc))
	{ push(@accs,$acc); }
my @sorted_accs = sort {$a <=> $b} @accs;
my $accs_total = @sorted_accs;
for(my $x=0;$x<@sorted_accs;$x++)
	{ my $accs_seen = $x+1;
	  my $acc = $sorted_accs[$x];
	  my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%{$seqs_per_gene_per_acc{$acc}}))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  my $genes_seen = 0; my $genes_total = @sorted_gene_ids;
	  my %cds_seqs = ();
	  foreach my $gene_id (@sorted_gene_ids)
		{ $genes_seen++;
		  print "accession $accs_seen of $accs_total, genes $genes_seen of $genes_total\n";
		  my $gene_seq 	  = $seqs_per_gene_per_acc{$acc}{$gene_id};
		  my @gene_seq	  = split(//,$gene_seq);
		  my $seq_length  = length($gene_seq);
		  my $cds_coords  = $coords{$gene_id}{cds_coords};
		  my $gene_start  = $coords{$gene_id}{gene_start};
		  my @cds_coords = split(/\, /,$cds_coords);
		  my %cds_positions = ();
		  foreach my $coords (@cds_coords)
			{ if ($coords =~ /^.*?\:(\d+)\-(\d+)$/)
				{ my $start = $1; my $end = $2;
				  for(my $pos=$start;$pos<=$end;$pos++)
					{ $cds_positions{$pos}++; }
				}
			}
		  my $cds_seq = '';
		  for(my $x=0;$x<@gene_seq;$x++)
			{ my $nt = $gene_seq[$x];
			  my $actual_pos = $gene_start+$x;
			  if (exists($cds_positions{$actual_pos}))
				{ $cds_seq .= $nt; }
			}
		  $cds_seqs{$gene_id} = $cds_seq;
		}
	  my $out_file = "$out_dir/$acc.fa";
	  open(OUT,'>',$out_file) or die $!;
	  foreach my $gene_id (@sorted_gene_ids)
		{ next if ( (!(defined($cds_seqs{$gene_id}))) or (!(defined($acc_gene_coords{$gene_id}{$acc}))) );
		  my $acc_cds_seq     = $cds_seqs{$gene_id};
		  my $acc_gene_coords = $acc_gene_coords{$gene_id}{$acc};
		  print OUT ">$gene_id|$acc_gene_coords\n$acc_cds_seq\n";
		}
	  close(OUT) or die $!;
	}

exit 1;