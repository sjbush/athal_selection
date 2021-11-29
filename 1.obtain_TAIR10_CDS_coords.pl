# PURPOSE: produce a table of CDS coords for each protein-coding gene in A. thaliana TAIR10

use strict;
use warnings;

# REQUIREMENTS
my $main_code  = 'code.txt'; # manually created
my $mito_code  = 'mito_code.txt'; # according to http://www.mun.ca/biology/scarr/MtDNA_code.html
my $utr_5	   = 'Arabidopsis_thaliana.TAIR10.5utr.fa'; # 5' UTR sequences from Ensembl BioMart: Gene stable ID | Transcript stable ID | Exon stable ID | 5' UTR start | 5' UTR end
my $utr_3	   = 'Arabidopsis_thaliana.TAIR10.3utr.fa'; # 3' UTR sequences from Ensembl BioMart: Gene stable ID | Transcript stable ID | Exon stable ID | 3' UTR start | 3' UTR end
my $exon_seqs  = 'Arabidopsis_thaliana.TAIR10.exon_sequences.fa'; # exon sequences from Ensembl BioMart: Gene stable ID | Transcript stable ID | Chromosome/scaffold name | Transcript start (bp) | Transcript end (bp) | Exon stable ID | Exon region start (bp) | Exon region end (bp) | Strand | Exon rank in transcript | CDS start | CDS end
my $exon_rank  = 'Arabidopsis_thaliana.TAIR10.exon_ranks.tsv'; # from Ensembl BioMart: Gene stable ID | Transcript stable ID | Exon stable ID | Exon rank in transcript
my $gene_info  = 'Arabidopsis_thaliana.TAIR10.gene_info_and_go_terms.tsv'; # from Ensembl BioMart: Gene stable ID | Gene description | Chromosome/scaffold name | Gene start (bp) | Gene end (bp) | Strand | Gene name | Gene type | GO term accession | GO term name | GO term evidence code | GO domain
my $ref_genome = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'; # wget ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
my $transcript_types = 'Arabidopsis_thaliana.TAIR10.transcript_types_and_chrs.tsv'; # from Ensembl BioMart: Transcript stable ID | Transcript type | Chromosome/scaffold name
my %trans; my %mito_trans;
LOAD_MAIN_CODE($main_code);
LOAD_MITO_CODE($mito_code);

# OUTPUT
my $out_file = 'Arabidopsis_thaliana.TAIR10.cds_coords.tsv';
open(OUT,'>',$out_file) or die $!;
print OUT "Transcript ID\tTranscript coords\tGene ID\tGene coords\tGene length (bp)\tExon IDs\tExon coords\tTotal exon length (bp)\t% of gene that is exonic\tCDS coords\tCDS length\t% of gene that is CDS\n";
open(FAIL,'>','failed_cds.txt') or die $!; # if any exons were missing from $exon_seqs, the CDS-finding step will fail. We include this check because $exon_seqs is quite a large file which we plan to get from BioMart - which can time out large downloads without warning, producing truncated output

# STORE A. THALIANA GENOME. WE WILL USE THIS TO TEST THAT THE CDS COORDINATES INFERRED BY THIS SCRIPT DO INDEED PRODUCE A LEGITIMATE TRANSLATED SEQUENCE.
my $header = ''; my %genome = ();
open(IN,$ref_genome) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $header = $1;	}
	  next if ($line =~ /^\>/);
	  $genome{$header} .= $line;
	}
close(IN) or die $!;

# STORE A. THALIANA GENE LENGTHS
my %gene_lengths = (); my %gene_coords = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $chr = $line[2]; my $gene_start = $line[3]; my $gene_end = $line[4]; my $gene_type = $line[7];
	  my $gene_len = ($gene_end-$gene_start)+1;
	  $gene_lengths{$gene_id} = $gene_len;
	  $gene_coords{$gene_id} = "$chr:$gene_start-$gene_end";
	}
close(IN) or die $!;

# STORE TRANSCRIPT TYPES AND CHRS
my %transcript_types_and_chrs = ();
open(IN,$transcript_types) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $transcript_id = $line[0]; my $transcript_type = $line[1]; my $chr = $line[2];
	  $transcript_types_and_chrs{$transcript_id}{type} = $transcript_type;
	  $transcript_types_and_chrs{$transcript_id}{chr} = $chr;
	}
close(IN) or die $!;

# STORE EXON RANKS PER TRANSCRIPT
my %exon_ranks_per_transcript = ();
open(IN,$exon_rank) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $transcript_id = $line[1]; my $exon_id = $line[2]; my $exon_rank_in_this_transcript = $line[3];
	  $exon_ranks_per_transcript{$exon_id}{$transcript_id} = $exon_rank_in_this_transcript;
	}
close(IN) or die $!;

# STORE 5' UTR COORDS
my %utr_coords = ();
my $gene_id = ''; my $transcript_id = ''; my $exon_ids = ''; my $utr_5_starts = ''; my $utr_5_ends = '';
open(IN,$utr_5) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line eq 'Sequence unavailable');
	  if ($line =~ /^\>(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)$/)
		{ $gene_id = $1; $transcript_id = $2; $exon_ids = $3; $utr_5_starts = $4; $utr_5_ends = $5; }
	  $utr_coords{$transcript_id}{start_5} = $utr_5_starts;
	  $utr_coords{$transcript_id}{end_5}   = $utr_5_ends;
	}
close(IN) or die $!;

# STORE 3' UTR COORDS
$gene_id = ''; $transcript_id = ''; $exon_ids = ''; my $utr_3_starts = ''; my $utr_3_ends = '';
open(IN,$utr_3) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line eq 'Sequence unavailable');
	  if ($line =~ /^\>(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)$/)
		{ $gene_id = $1; $transcript_id = $2; $exon_ids = $3; $utr_3_starts = $4; $utr_3_ends = $5; }
	  $utr_coords{$transcript_id}{start_3} = $utr_3_starts;
	  $utr_coords{$transcript_id}{end_3}   = $utr_3_ends;
	}
close(IN) or die $!;

# DETERMINE CDS COORDS
my %exon_seqs = (); my %transcript_coords = (); my %transcript_strand = (); my %genes_per_transcript = ();
$gene_id = ''; my $transcript_ids = ''; my $chr = ''; my $transcript_start = ''; my $transcript_end = ''; my $exon_id = ''; my $exon_start = ''; my $exon_end = ''; my $strand = ''; my $exon_ranks = ''; my $cds_start = ''; my $cds_end = '';
my $eof;
open(IN,$exon_seqs) or die $!;
while(<IN>) { $eof = $.; }
close(IN) or die $!;
open(IN,$exon_seqs) or die $!;
while(<IN>)
	{ my $pc = sprintf("%.2f",(($./$eof)*100));  print "determining CDS coords: $pc%\n";
	  my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)$/)
		{ $gene_id = $1; $transcript_ids = $2; $chr = $3; $transcript_start = $4; $transcript_end = $5; $exon_id = $6; $exon_start = $7; $exon_end = $8; $strand = $9; $exon_ranks = $10; $cds_start = $11; $cds_end = $12; }
	  elsif ($line =~ /^\>(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)\|(.*?)$/)
		{ $gene_id = $1; $transcript_ids = $2; $chr = $3; $transcript_start = $4; $transcript_end = $5; $exon_id = $6; $exon_start = $7; $exon_end = $8; $strand = $9; $exon_ranks = $10; $cds_start = '';  $cds_end = '';  }
	  next if ($line =~ /^\>/);
	  if ($strand == -1)
		{ my $exon_start_strand_swapped = $exon_start;
		  my $exon_end_strand_swapped   = $exon_end;
		  $exon_start = $exon_start_strand_swapped;
		  $exon_end   = $exon_end_strand_swapped;
		}
	
	  # A WORD ON EXON RANKS:
	  # ALL EXON IDS HAVE THE FORMAT "GENE_ID.TRANSCRIPT_ID.EXON_RANK", BUT THIS CAN BE MISLEADING. EXON IDS DENOTE UNIQUE SEQUENCES, BUT SOME EXONS ARE PRESENT IN MULTIPLE TRANSCRIPTS. IN THIS CASE, THE EXON ID WILL INCLUDE THE RANK FROM ONLY *ONE* OF THE TRANSCRIPTS IT IS WITHIN. FOR INSTANCE, THE FIFTH EXON OF AT4G08535.1 IS NAMED "AT4G08535.1.exon5", BUT THIS SEQUENCE (NOTING ITS NAME CONTAINS 'exon5') IS ALSO THE *SIXTH* EXON OF AT4G08535.2.
	  # IT'S ALSO POSSIBLE FOR $#transcript_ids != $#exon_ranks. THIS WILL HAPPEN IF AN EXON IS, E.G., RANK 2 IN 2 DIFFERENT TRANSCRIPTS, BUT RANK 3 IN A THIRD TRANSCRIPT: THREE TRANSCRIPT IDs IN TOTAL BUT ONLY 2 RANKS. YOU CAN SEE THIS IF PARSING THE $exon_seqs FILE FOR "AT1G01210.3.exon2" - THIS IS THE SECOND EXON OF TRANSCRIPTS AT1G01210.1 AND AT1G01210.3 BUT THE THIRD EXON OF TRANSCRIPT AT1G01210.2. TO RESOLVE THIS ISSUE, WE WILL REFER TO A SEPARATE TABLE.
	  my @transcript_ids = split(/\;/,$transcript_ids);
	  my @exon_ranks = ();
	  foreach my $transcript_id (@transcript_ids)
		{ if (exists($exon_ranks_per_transcript{$exon_id}{$transcript_id}))
			{ my $exon_rank = $exon_ranks_per_transcript{$exon_id}{$transcript_id};
			  push(@exon_ranks,$exon_rank);
			}
		}
	  next if ($#exon_ranks == -1); # FILTER: unable to proceed because data is unavailable
	  
	  my %exon_ranks = map {$_ => 1} @exon_ranks;
	  for(my $x=0;$x<@transcript_ids;$x++)
		{ my $transcript_id = $transcript_ids[$x];
		  $genes_per_transcript{$transcript_id} = $gene_id;
		}
	  for(my $x=0;$x<@transcript_ids;$x++)
		{ my $transcript_id = $transcript_ids[$x];
		  my $exon_rank		= $exon_ranks[$x];
		  push(@{$exon_seqs{$transcript_id}},[$exon_rank,$chr,$strand,$exon_id,$exon_start,$exon_end,$transcript_id,$transcript_start,$transcript_end,$cds_start,$cds_end]);
		  $transcript_strand{$transcript_id} = $strand;
		  $transcript_coords{$transcript_id}{transcript_start} = $transcript_start;
		  $transcript_coords{$transcript_id}{transcript_end}   = $transcript_end;
		  $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{exon_id} = $exon_id;
		  $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{exon_coords} = "$chr:$exon_start-$exon_end";
		  my $cds_start = $exon_start; my $cds_end = $exon_end;
		  my $skip = 0;
		  if (exists($utr_coords{$transcript_id}))
			{ my $utr_5_start = ''; if (exists($utr_coords{$transcript_id}{start_5})) { $utr_5_start = $utr_coords{$transcript_id}{start_5}; }
			  my $utr_3_start = ''; if (exists($utr_coords{$transcript_id}{start_3})) { $utr_3_start = $utr_coords{$transcript_id}{start_3}; }
			  my $utr_5_end   = ''; if (exists($utr_coords{$transcript_id}{end_5}))   { $utr_5_end   = $utr_coords{$transcript_id}{end_5};   }
			  my $utr_3_end   = ''; if (exists($utr_coords{$transcript_id}{end_3}))   { $utr_3_end   = $utr_coords{$transcript_id}{end_3};   }
			  my @utr_5_start = split(/\;/,$utr_5_start);
			  my @utr_3_start = split(/\;/,$utr_3_start);
			  my @utr_5_end   = split(/\;/,$utr_5_end);
			  my @utr_3_end   = split(/\;/,$utr_3_end);
			  if ($strand == -1)
				{ my @utr_5_start_strand_swapped = @utr_3_start;
				  my @utr_3_start_strand_swapped = @utr_5_start;
				  my @utr_5_end_strand_swapped   = @utr_3_end;
				  my @utr_3_end_strand_swapped   = @utr_5_end;
				  @utr_5_start = @utr_5_start_strand_swapped;
				  @utr_3_start = @utr_3_start_strand_swapped;
				  @utr_5_end   = @utr_5_end_strand_swapped;
				  @utr_3_end   = @utr_3_end_strand_swapped;
				}
			  for(my $x=0;$x<@utr_5_start;$x++)
				{ my $utr_start = $utr_5_start[$x];
				  my $utr_end   = $utr_5_end[$x];
				  if (($utr_start == $exon_start) and ($utr_end == $exon_end))
					{ $skip++; }
				  elsif (($utr_start >= $exon_start) and ($utr_end <= $exon_end))
					{ $cds_start = $utr_end+1; }
				}
			  for(my $x=0;$x<@utr_3_start;$x++)
				{ my $utr_start = $utr_3_start[$x];
				  my $utr_end   = $utr_3_end[$x];
				  if (($utr_start == $exon_start) and ($utr_end == $exon_end))
					{ $skip++; }
				  elsif (($utr_start >= $exon_start) and ($utr_end <= $exon_end))
					{ $cds_end = $utr_start-1; }
				}
			}
		  $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{cds_coords} = "$chr:$cds_start-$cds_end" unless ($skip > 0);
		}
	}
close(IN) or die $!;

my @transcript_ids = ();
while((my $transcript_id,my $irrel)=each(%exon_seqs))
	{ push(@transcript_ids,$transcript_id); }
my @sorted_transcript_ids = sort {"\L$a" cmp "\L$b"} @transcript_ids;
foreach my $transcript_id (@sorted_transcript_ids)
	{ next if (!(exists($transcript_types_and_chrs{$transcript_id}{chr}))); # NOTE: sometimes $chr is undefined. This is because some of the more transcripts, such as ENSRNA049492369-T1, aren't listed in the $transcript_types file.
	  my $chr = $transcript_types_and_chrs{$transcript_id}{chr};
	  my @exon_ranks = ();
	  while((my $exon_rank,my $irrel)=each(%{$transcript_coords{$transcript_id}{exon_ranks}}))
		{ push(@exon_ranks,$exon_rank); }
	  my @sorted_exon_ranks = sort {$a <=> $b} @exon_ranks;
	  my $strand = $transcript_strand{$transcript_id};
	  if ($strand == -1)
		{ @sorted_exon_ranks = reverse(@sorted_exon_ranks); }	  
	  my $exons_in_this_transcript = ''; my $exon_coords_in_this_transcript = '';
	  my $total_exon_length = 0;
	  my @transcript_positions = ();
	  foreach my $exon_rank (@sorted_exon_ranks)
		{ my $exon_id 	  = $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{exon_id};
		  my $exon_coords = $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{exon_coords};
		  if ($exon_coords =~ /^(.*?)\:(\d+)\-(\d+)$/)
			{ my $exon_start = $2; my $exon_end = $3;
			  my $exon_length = ($exon_end-$exon_start)+1;
			  push(@transcript_positions,$exon_start);
			  push(@transcript_positions,$exon_end);
			  $total_exon_length += $exon_length;
			}
		  $exons_in_this_transcript 	  .= "$exon_id, ";
		  $exon_coords_in_this_transcript .= "$exon_coords, ";
		}
	  $exons_in_this_transcript 	  =~ s/\, $//;
	  $exon_coords_in_this_transcript =~ s/\, $//;
	  my $cds_coords_in_this_transcript = ''; my $total_cds_length = 0; my $cds_seq = '';
	  foreach my $exon_rank (@sorted_exon_ranks)
		{ if (exists($transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{cds_coords}))
			{ my $cds_coords = $transcript_coords{$transcript_id}{exon_ranks}{$exon_rank}{cds_coords};
			  $cds_coords_in_this_transcript .= "$cds_coords, ";
			  if ($cds_coords =~ /^(.*?)\:(\d+)\-(\d+)$/)
				{ my $cds_start = $2; my $cds_end = $3;
				  my $cds_length = ($cds_end-$cds_start)+1;
				  $total_cds_length += $cds_length;
				  my $seq = substr($genome{$chr},$cds_start-1,$cds_length);
				  $cds_seq .= $seq;
				}
			}
		}
	  $cds_coords_in_this_transcript =~ s/\, $//;
	  my @sorted_transcript_positions = sort {$a <=> $b} @transcript_positions;
	  my $transcript_start = $transcript_positions[0]; my $transcript_end = $transcript_positions[$#transcript_positions];
	  my $gene_id = $genes_per_transcript{$transcript_id};
	  my $gene_length = $gene_lengths{$gene_id};
	  my $gene_coords = $gene_coords{$gene_id};
	  my $pc_of_transcript_that_is_exonic = sprintf("%.2f",(($total_exon_length/$gene_length)*100));
	  my $pc_of_transcript_that_is_cds = sprintf("%.2f",(($total_cds_length/$gene_length)*100));
	  
	  # FOR PROTEIN-CODING TRANSCRIPTS, WE'LL TEST THAT THE TRANSLATED CDS SEQUENCE IS LEGITIMATE BEFORE REPORTING OUTPUT
	  my $transcript_type = $transcript_types_and_chrs{$transcript_id}{type};
	  my $failure = 0;
	  if (($transcript_type eq 'protein_coding') and ($cds_seq ne ''))
		{ if ($strand == -1)
			{ $cds_seq =~ tr/[ATCG]/[TAGC]/;
			  $cds_seq = reverse($cds_seq);
			}
		  my @cds = split(//,$cds_seq);
		  my $translated_cds = '';
		  for(my $x=0;$x<@cds;$x++)
			{ my $nt1 = $cds[$x];
			  next if ( (!(defined($cds[$x+1]))) || (!(defined($cds[$x+2]))) );
			  my $nt2 = $cds[$x+1]; my $nt3 = $cds[$x+2];
			  my $codon = "$nt1"."$nt2"."$nt3";
			  my $aa = 'X';
			  if ($chr eq 'M')
				{ $aa = $mito_trans{$codon} unless (!(exists($mito_trans{$codon}))); }
			  else
				{ $aa =      $trans{$codon} unless (!(exists(     $trans{$codon}))); }
			  $translated_cds .= $aa;
			  $x = $x+2;
			}
		  if ($translated_cds =~ /^(.*?)\*$/) { $translated_cds = $1; } # an allowance is made for a terminator * so this won't fail the subsequent test for internal stops
		  if ($translated_cds =~ /\*/)
			{ $failure++; } # print "error: internal stop in translated CDS for $transcript_id ($cds_coords_in_this_transcript): $translated_cds\n"; # these errors are due to MISSING EXONS in exon_sequences.fa
		}
	  
	  print OUT "$transcript_id\t$transcript_start-$transcript_end\t$gene_id\t$gene_coords\t$gene_length\t$exons_in_this_transcript\t$exon_coords_in_this_transcript\t$total_exon_length\t$pc_of_transcript_that_is_exonic\t$cds_coords_in_this_transcript\t$total_cds_length\t$pc_of_transcript_that_is_cds\n";
	  if ($failure > 0) { print FAIL "$transcript_id --> $exons_in_this_transcript\n"; }
	}

close(OUT) or die $!; close(FAIL) or die $!;
exit 1;

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