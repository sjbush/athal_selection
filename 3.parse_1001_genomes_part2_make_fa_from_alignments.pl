use strict;
use warnings;

# REQUIREMENTS
my $in_dir = '1001_genomes_alignments'; # from 2.parse_1001_genomes_part1_make_alignments_from_vcf.pl

# OUTPUT
my $out_dir = '1001_genomes_as_fa';
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

opendir(DIR,$in_dir) or die $!;
my @subdirs = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_subdirs = sort {$a cmp $b} @subdirs;
my $subdirs_seen = 0; my $subdirs_total = @sorted_subdirs; $subdirs_total = $subdirs_total-2;
foreach my $subdir (@sorted_subdirs)
	{ next if (($subdir eq '.') or ($subdir eq '..'));
	  $subdirs_seen++;
	  next if (!(-d("$in_dir/$subdir")));
	  opendir(DIR,"$in_dir/$subdir") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  my @sorted_files = sort {$a cmp $b} @files;
	  my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
	  foreach my $file (@sorted_files)
		{ next if (($file eq '.') or ($file eq '..'));
		  $files_seen++;
		  print "$subdirs_seen of $subdirs_total, $files_seen of $files_total...\n";
		  my $acc = ''; my $chr = '';
		  if ($file =~ /^(.*?)\-(.*?)\-TAIR10\.alignment$/)
			{ $acc = $1; $chr = $2; }
		  next if ($file !~ /^.*?\-.*?\-TAIR10\.alignment$/);
		  next if (-e("$out_dir/$acc/$acc-$chr.fa")); # CHECKPOINT: we have seen this already
		  my $begin = 0; my $number_in_alignment_pair = 1;
		  my $first_start_pos_for_tair10 = 0; my $first_start_pos_for_acc = 0;
		  my %alignments = ();
		  open(IN,"$in_dir/$subdir/$file") or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  if (($line =~ /^\-- BEGIN alignment \[(.*?)\].*?$/) and ($begin == 0))
				{ $begin++; }
			  elsif (($line =~ /^\-- END.*?$/) and ($begin == 1))
				{ $begin = 0; $first_start_pos_for_tair10 = 0; $first_start_pos_for_acc = 0; }
			  if (($begin > 0) and ($line =~ /^(\d+)\s+(.*?)$/))
				{ my $start_pos = $1; my $seq = $2;
				  if (($first_start_pos_for_tair10 == 0) and ($number_in_alignment_pair == 1))
					{ $first_start_pos_for_tair10 = $start_pos; }
				  elsif (($first_start_pos_for_acc == 0) and ($number_in_alignment_pair == 2))
					{ $first_start_pos_for_acc = $start_pos; }
				  if ($number_in_alignment_pair == 1)
					{ $alignments{$first_start_pos_for_tair10}{tair10} .= uc($seq); }
				  elsif ($number_in_alignment_pair == 2)
					{ $alignments{$first_start_pos_for_acc}{$acc} .= uc($seq); }
				  $number_in_alignment_pair++;
				}
			  elsif (($begin > 0) and ($line !~ /^(\d+)\s+(.*?)$/))
				{ $number_in_alignment_pair = 1; }
			}
		  close(IN) or die $!;
		  
		  # output alignments as fasta
		  if (!(-d("$out_dir/$acc"))) { mkdir "$out_dir/$acc" or die $!; }
		  open(OUT,'>',"$out_dir/$acc/$acc-$chr.fa") or die $!;
		  my @start_pos = ();
		  while((my $start_pos,my $irrel)=each(%alignments))
			{ push(@start_pos,$start_pos); }
		  my @sorted_start_pos = sort {$a <=> $b} @start_pos;
		  foreach my $start_pos (@sorted_start_pos)
			{ my $tair10_seq = $alignments{$start_pos}{tair10};
			  my $acc_seq    = $alignments{$start_pos}{$acc};
			  $tair10_seq =~ s/\./\-/g; $acc_seq =~ s/\./\-/g;
			  print OUT ">TAIR10 | $start_pos\n$tair10_seq\n";
			  print OUT ">$acc | $start_pos\n$acc_seq\n";
			}
		  close(OUT) or die $!;
		}
	}

exit 1;