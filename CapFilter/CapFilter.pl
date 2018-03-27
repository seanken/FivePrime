#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use File::Basename;

my ($command, $results, $fasta, $seq_file, $soft, $barlen,$gzip, $gzip_out, $output_file, $bam_file, $genome_file, $peak_file, $bed, $peak_cutoff, $peak_header_line, @peak_ref_order, $sorted, $help);

$command = shift;

if (!defined($command)) {
	usage();
}

$output_file = '';
$help = 0;
if ($command eq 'seq') {

	$gzip_out = 0;
	$gzip = 0;
	$fasta = 0;

	my $results = GetOptions('fasta' => \$fasta,
									 'gzip' => \$gzip,
									 'gzip_out' => \$gzip_out,
									 'out=s' => \$output_file,
									 'help' => \$help);

	seq_usage() if $help;

	$seq_file = shift;

	if (!defined($seq_file)) {
		print "Error: no FASTQ file supplized. Exiting.\n";
		seq_usage();
	}
	$gzip = 1 if $seq_file =~ m/\.gz$/i;

	process_seq_file();

} elsif ($command eq 'peak') {

	#open(DEBUG, ">debug.txt") || die("Error opening file '$output_file' for writing: $!\n");
	$sorted = 0;
	$bed = 0;
	$soft=12;
	$barlen=9;
	$peak_cutoff = 50;
	print "It begins!";
	my $results = GetOptions('bed'   => \$bed,
									 'cutoff=i' => \$peak_cutoff,
									 'out=s' => \$output_file,
									 'sorted' => \$sorted,
									 'help' => \$help,
									 'soft=i' => \$soft,
									 'barcodelen=i' => \$barlen);

	print STDERR $peak_cutoff;
	peak_usage() if $help;

	$peak_file = shift;
	$bam_file = shift;
	$genome_file = shift;

	my (%ref_peaks, %genome, @bam_headers, $current_total);
	print STDERR "LOAD GENOME\n";

	%genome = parse_genome();

	print STDERR "Load peaks\n";

	%ref_peaks = parse_peaks();

	$current_total = 0;

	# Data used for tracking sorted files
	my (@ref_names, $current_ref, $current_ref_index, $current_index);
	@ref_names = keys(%ref_peaks);
	@ref_names = sort { $a cmp $b } @ref_names;
	($current_ref, $current_ref_index, $current_index) = ($ref_names[0], 0, 0);

	# TEST
	#my $test_peak = $ref_peaks{$current_ref}->[$current_index];
	#print "$test_peak->{'line'}\n$current_ref:$test_peak->{'start'}..$test_peak->{'end'}\n";

	@bam_headers = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');
	open(BAM, "samtools view $bam_file |") || die("Error opening BAM file '$bam_file': $!\n");
	print STDERR "Reading in BAM\n";
	while (my $line = <BAM>) {
		chomp $line;

		$current_total++;
		#exit 0 if $current_total>10000;
		print STDERR "\rProcessing $bam_file: $current_total" if $current_total % 1000 == 0 || eof(BAM);
		print STDERR "\n" if eof(BAM);

		my (%hit, @values, @extra, $start, $end, $strand, @overlaps);

		@values = split("\t", $line);

		for (my $x = 0; $x < scalar(@values); $x++) {
			if ($x >= scalar(@bam_headers)) {
				push @extra, $values[$x];
				next;
			}
			$hit{$bam_headers[$x]} = $values[$x];
		}

		$hit{'extra'} = \@extra if scalar(@extra);

		if (!exists($hit{'qname'}) || !exists($hit{'rname'}) || !exists($hit{'pos'}) || !exists($hit{'flag'}) || !exists($hit{'cigar'})) {
			print "Error in parsing SAM alignment - exiting.\nLast SAM data line read: $line\n";
			exit(0);
		}

		$start = $hit{'pos'};
		$end = $hit{'pos'} - 1;
		while ($hit{'cigar'} =~ m/([0-9]+)[NM]/g) {
			$end += $1;
		}
		$strand = $hit{'flag'} & 16;
		$strand = ($strand)? "-" : "+";

		next if !exists($ref_peaks{$hit{'rname'}}) || !exists($genome{$hit{'rname'}}) || $hit{'rname'} eq '*';

		if ($sorted) {
			# Update current ref if it appears we have gone past this sorted ref
			while(($current_ref ne $hit{'rname'}) && $current_ref_index != -1) {
				# We are still parsing through the end of the reads on this chromosome - don't keep updating.
				if ($hit{'rname'} lt $current_ref) {
					$current_index = -1;
					last;
				}

				# Reset to new position
				$current_ref_index++;
				$current_index = 0;

				# There are no more references to go through now that we have finished - let's stop parsing this data.
				if ($current_ref_index >= scalar(@ref_names)) {
					$current_ref_index = -1;
					last;
				}

				$current_ref = $ref_names[$current_ref_index];
			}

			next if $current_ref_index == -1;

			# Update current peak position if it appears we have gone beyond the current peak position
			while ($current_index != -1) {
				if (!exists($ref_peaks{$current_ref})) {
					$current_index = -1;
					last;
				}

				my $current_peak = $ref_peaks{$current_ref}->[$current_index];

				# We are behind on peak position - move forward to only compare the closest peaks
				if ($current_peak->{'end'} < $start) {
					$current_index++;
					if ($current_index >= scalar(@{$ref_peaks{$current_ref}})) {
						$current_index = -1;
						last;
					}
				} else {
					last;
				}
			}

			next if $current_index == -1;

			my $start_pos = 0;
			while (1) {
				last if $current_ref ne $hit{'rname'};

				my ($peak, $overlap, $edge);

				$peak = $ref_peaks{$current_ref}->[$current_index + $start_pos];
				($overlap, $edge) = overlaps($start, $end, $peak->{'start'}, $peak->{'end'});

				if ($overlap) {
					push @overlaps, [$edge, $peak, $current_index + $start_pos];
				}

				$start_pos++;

				last if $peak->{'start'} > $end || ($current_index + $start_pos) >= scalar(@{$ref_peaks{$current_ref}});
			}
		} else {
			@overlaps = find_overlap({'start' => $start, 'end' => $end}, $ref_peaks{$hit{'rname'}});
		}
	
		##Commented out!
		#if (!($hit{'qname'} =~ m/[ACGTN]\)$/)) {
		#	print "Error: SAM 'qname' field is not properly formatted - exiting.\nLast SAM qname data: $hit{'qname'}\n";
		#	exit(0);
		#}

		# The read base should be encoded in the second to last part of the 'qname' field value
		my ($genome_base, $read_base);
		##Commented is original version, bellow is our version
		#$read_base = substr($hit{'qname'}, -2, 1);
		$read_base= ($strand eq '+') ? substr($hit{'seq'},$barlen,1) : substr($hit{'seq'},-$barlen-1,1) ;
		
		if ($strand eq '-') {$read_base =~ tr/acgtACGT/tgcaTGCA/;}
		my $fordebug = substr($genome{$hit{'rname'}}, $start-4, $end-$start+8);
		
		#print DEBUG $hit{'seq'};
		#print DEBUG "\n";
		#print DEBUG $fordebug;
		#print DEBUG "\n";
		#print DEBUG $hit{'cigar'};
		#print DEBUG "\n";
		#print DEBUG $read_base;
		#print DEBUG "\n";
		
		foreach my $overlap (@overlaps) {
			my ($edge, $peak, $index);

			($edge, $peak, $index) = @{$overlap};

			# Make sure that the read strand and peak strand are identical (in cases where peak strand is noted)
			next if  ($peak->{'strand'} eq '+' && $strand eq '-') || 
						($peak->{'strand'} eq '-' && $strand eq '+');

			# Make sure the 5' base of the read overlaps the peak - assumes the first sequenced base of the read is the most 5'
			next if ($start < $peak->{'start'} && $peak->{'strand'} eq '+' && $strand eq '+') ||
					  ($end > $peak->{'end'} && $peak->{'strand'} eq '-' && $strand eq '-') ;

			$peak->{'total_reads'}++;

			# It is assumed that the base that should match the genome is *not* encoded in the read, 
			# and therefore all 'substr' calculations compensate for this 1-base offset
			#
			#commented is orginal! Also, added 2+ to elsif below
			#my $cap_pos = ($strand eq '+') ? $start - 1 : $end + 1;
			my $cap_pos = ($strand eq '+') ? $start - $soft +$barlen : $end + $soft -$barlen;
			if ($strand eq '+' && $start > 3) { #1) {
				$genome_base = substr($genome{$hit{'rname'}}, $cap_pos - 1, 1);
				$genome_base = uc($genome_base);
				if ($read_base eq 'G' && $genome_base ne $read_base) {
					$peak->{'total_unencoded'}++;
				}
			} elsif ($strand eq '-' && $end < 2+ length($genome{$hit{'rname'}})) {
				$genome_base = substr($genome{$hit{'rname'}}, $cap_pos - 1, 1);
				$genome_base =~ tr/acgtACGT/TGCATGCA/;
				if ($read_base eq 'G' && $genome_base ne $read_base) {
					$peak->{'total_unencoded'}++;
				}
			}
			#print DEBUG $genome_base;
			#print DEBUG "hi\n";
		}
	}
	close(BAM);
	#close(DEBUG);
	my $out_fh = *STDOUT;

	if ($output_file ne '') {
		open(OUT, ">$output_file") || die("Error opening file '$output_file' for writing: $!\n");
		$out_fh = *OUT;
	}

	if (!$bed) {
		print $out_fh $peak_header_line.",%-Capped\n";
	}

	# Use this for adding feature name for BED formatted files
	my $peak_id_num = 1;

	# Print out the final list of peaks that meet our filter criteria
	foreach my $ref (@peak_ref_order) {
		my $peaks = $ref_peaks{$ref};
		foreach my $peak (@{$peaks}) {
			next if !$peak->{'total_reads'};

			my $percent = $peak->{'total_unencoded'} / $peak->{'total_reads'} * 100;

			if ($percent >= $peak_cutoff) {
				my $to_print;
				if (!$bed) {
					print $out_fh "$peak->{'line'},$percent\n";
				} else {
					my $score = int($percent * 10 + 0.5);
					my @values = split("\t", $peak->{'line'});
					if (scalar(@values) < 5) {
						if (scalar(@values) < 4) {
							push @values, "peak_$peak_id_num";
							$peak_id_num++;
						}
						push @values, $score;
					} else {
						$values[4] = $score;
					}
					print $out_fh join("\t", @values)."\n";
				}
			}
		}
	}
	if ($output_file ne '') {
		close($out_fh);
	}
} else {
	print "Command '$command' not supported. Exiting.\n";
	exit(0);
}


sub parse_genome {
	my (%genome, $seq, $id);

	($seq, $id) = ('', '');

	open(DAT, $genome_file) || die("Error opening genome file '$genome_file': $!\n");
	while (my $line = <DAT>) {
		chomp $line;

		$seq .= $line if eof(DAT);

		if ($line =~ m/^>([^\s]+).*$/ || eof(DAT)) {
			if ($id ne '') {
				$genome{$id} = $seq;
			}
			$id = $1;
			$seq = '';
		} else {
			$seq .= $line;
		}
	}
	close(DAT);

	return %genome
}

sub parse_peaks {
	my ($header_line, @headers, %ref_peaks, %found_refs);

	open(DAT, $peak_file) || die("Error opening peak file '$peak_file': $!\n");
	if (!$bed) {
		$header_line = <DAT>;
		chomp $header_line;
		$peak_header_line = $header_line;
		@headers = split(",", $header_line);
	} else {
		@headers = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts');
	}
	while (my $line = <DAT>) {
		chomp $line;

		my (%peak, @values, @extra, $ref, $start, $end, $strand);

		if (!$bed) {
			@values = split(",", $line);
		} else {
			@values = split("\t", $line);
		}

		for (my $x = 0; $x < scalar(@values); $x++) {
			if ($x >= scalar(@headers)) {
				push @extra, "$values[$x]";
				next;
			}
			$peak{$headers[$x]} = $values[$x];
		}
		$peak{'extra'} = \@extra if scalar(@extra);

		if (!$bed) {
			if (!defined($peak{'Start'}) || !defined($peak{'End'}) || !defined($peak{'Strand'}) || !defined($peak{'Chromosome'})) {
				print "Error - peak file format not recognized.\n";
				exit(0);
			}
			$ref = $peak{'Chromosome'};
			$start = $peak{'Start'};
			$end = $peak{'End'};
			$strand = $peak{'Strand'};
			$strand = '.' if $peak{'Strand'} ne '+' && $peak{'Strand'} ne '-';
		} else {
			if (!defined($peak{'chromStart'}) || !defined($peak{'chromEnd'}) || !defined($peak{'chrom'})) {
				print "Error - peak file format not recognized.\n";
				exit(0);
			}
			$ref = $peak{'chrom'};
			$start = $peak{'chromStart'};
			$end = $peak{'chromEnd'};
			$strand = $peak{'strand'};
			$strand = '.' if !exists($peak{'strand'}) || ($peak{'strand'} ne '+' && $peak{'strand'} ne '-');
		}

		$ref_peaks{$ref} = [] if !exists($ref_peaks{$ref});

		push @peak_ref_order, $ref if !exists($found_refs{$ref});
		$found_refs{$ref} = 1;

		push @{$ref_peaks{$ref}}, {'start' => $start, 'end' => $end, 'strand' => $strand, 'total_reads' => 0, 'total_unencoded' => 0, 'line' => $line};
	}
	close(DAT);

	while (my ($ref, $peaks) = each %ref_peaks) {
		@{$peaks} = sort { $a->{'start'} <=> $b->{'start'} } @{$peaks};
	}

	return %ref_peaks;
}

sub process_seq_file {
	my ($out_fh, $in_fh);

	if ($output_file ne '') {
		if (!$gzip_out) {
			open($out_fh, "> $output_file") || die("Error opening file '$output_file' for writing: $!\n");
		} else {
			open($out_fh, "| gzip - > $output_file") || die("Error opening file '$output_file' for writing: $!\n");
		}
	} else {
		$out_fh = *STDOUT;
	}

	if ($gzip) {
		open($in_fh, "gzip -dc $seq_file |") || die("Error opening file '$seq_file': $!\n");
	} else {
		open($in_fh, "$seq_file") || die("Error opening file '$seq_file': $!\n");
	}

	while (!eof($in_fh)) {
		my ($seq_id, $seq, $score_id, $score);

		$seq_id = <$in_fh>;
		chomp $seq_id;
		$seq = <$in_fh>;
		chomp $seq;

		if (!$fasta) {
			$score_id = <$in_fh>;
			chomp $score_id;
			$score = <$in_fh>;
			chomp $score;
		}

		$seq_id =~ s/ /\_/g;

		my $first_base = substr($seq, 0, 1, '');
		print $out_fh "$seq_id($first_base)\n$seq\n";

		if (!$fasta) {
			substr($score, 0, 1, '');
			print $out_fh "$score_id\n$score\n";
		}
	}
	close($in_fh);

	close($out_fh) if $output_file ne '';
}

sub find_overlap {
	my ($hit, $genes, $longest, @overlaps, $total_genes, $start, $min, $max, $mid, $cont, $no_match_found, $within);

	($hit, $genes, $longest, $within) = @_;

	$longest = 40000 if !defined($longest);

	$within = 0 if !defined($within);

	$total_genes = scalar(@{$genes});
	$start = $hit->{'end'} + $within;

	($min, $max, $mid, $cont, $no_match_found) = (0, $total_genes, int($total_genes / 2 + .5), 1, 0);

	$mid = 0 if $mid >= $total_genes;

	do {
		my $comp = $genes->[$mid]->{'start'} - $within;
		if ($mid < $total_genes) {
			if ($start > ($genes->[$mid]->{'start'} - $within)) {
				$min = $mid;
				if ($max - $min <= 1) {
					$mid++;
					$cont = 0;
				} else {
					$mid += int(($max - $mid) / 2 + .5);
				}
			} elsif ($start == ($genes->[$mid]->{'start'} - $within)) {
				$cont = 0;
			} else {
				$max = $mid;
				if ($max - $min <= 1) {
					$cont = 0;
				} else {
					$mid -= int(($mid - $min) / 2 + .5);
				}
			}
		} else {
			$cont = 0;
		}
	} while ($cont);

	if ($mid < $total_genes) {
		my $gene_adjusted_start = $genes->[$mid]->{'start'} - $within;
		$gene_adjusted_start = 1 if $gene_adjusted_start < 1;
		while (($genes->[$mid]->{'start'} - $within) == $start && $mid < $total_genes - 1) {
			$mid++;
		}
	}

	$mid = $total_genes - 1 if $mid >= $total_genes;

	my $farthest = ($hit->{'start'} - ($longest + $within) > 0)? $hit->{'start'} - ($longest  + $within) : 1;

	my $does_overlap = 0;

	if (!defined($genes->[$mid]->{'start'})) {
		while (my ($key, $val) = each %{$genes->[$mid]}) {
			print "$key=$val\n";
		}
		die("WOOPS! Critical error...exiting.\n");
	}

	my $edge;
	while (($farthest) <= $genes->[$mid]->{'start'} && $mid >= 0) {
		($does_overlap, $edge) = overlaps($hit->{'start'}, $hit->{'end'}, $genes->[$mid]->{'start'}, $genes->[$mid]->{'end'}, $within);
		push @overlaps, [$edge, $genes->[$mid], $mid] if $does_overlap;

		$mid--;
	}

	return @overlaps;
}

sub overlaps {
	my ($s1, $e1, $s2, $e2, $within) = @_;

	$within = 0 if !defined($within);

	# Increase Overlap area
	$s1 -= $within;
	$e1 += $within;

	$s2 -= $within;
	$e2 += $within;

	my $overlap = 0;

	my $edge = 0;

	if ($s1 >= $s2 && $s1 <= $e2) {
		$overlap = 1;
	} elsif ($e1 >= $s2 && $e1 <= $e2) {
		$overlap = 1;
	} elsif ($s1 < $s2 && $e1 > $e2) {
		$overlap = 1;
	}

	$edge = 1 if $overlap && ($s1 < $s2 || $e1 > $e2);

	return ($overlap, $edge);
}

sub seq_usage {
	my $script = basename($0);

	print "usage: $script seq (--help --gzip --out <output_file> --fasta) <seq_file>\n\n";
	print "Options:\n";
	print "     --help                       Print this help message.\n";
	print "     --gzip                       FASTQ file is in gzip format.\n";
	print "     --gzip_out                   Output is compressed using gzip format. Requires 'output_file'(default: Not compressed).\n";
	print "     --out  <output_file>         Output file processed data is written to (default: STDOUT).\n";
	print "     --fasta                      Input file is in FASTA format. (default: FASTQ).\n";
	print "\n\nProgram Description\n\n";
	print "The 'seq' command in CapFilter processes a FASTQ/FASTA sequence file by removing the first\n";
	print "base of a sequenced read, and placing that in the header field of a FASTQ/FASTA alignment.\n";
	print "FASTQ files will also have the first letter of the score removed as well. This modified\n";
	print "header will then be parsed from the BAM file produced during alignment to generate a cap\n";
	print "filtered peak file for the end user using the 'peak' command found in the CapFilter software.\n";

	exit(0);
}

sub peak_usage {
	my $script = basename($0);

	print "usage: $script peak (--help --bed --cutoff --out <output_file>) <peak_file> <bam_file> <genome_file>\n\n";
	print "Options:\n";
	print "     --help                       Print this help message.\n";
	print "     --bed                        Peak file is in BED format.\n";
	print "     --sorted                     Input BAM file is sorted by coordinate.  Speeds up processing.\n";
	print "     --cutoff                     Confidence cutoff for a peak based on percent of reads with\n"; 
	print "                                  an unencoded 'G'. (default: 50).\n";
	print "     --out  <output_file>         Output file processed data is written to (default: STDOUT).\n";
	print "\n\nProgram Description\n\n";
	print "The 'peak' command in CapFilter takes as input a peak file, BAM alignment file produced by aligning\n";
	print "the modified FASTQ/FASTA file from the 'seq' command, and genome sequence file to identify the\n";
	print "percent of reads within a peak that start with an unencoded 'G' (based on the genome refence\n";
	print "sequence). It will automatically filter out those peaks below the cutoff.\n\n";

	exit(0);
}

sub usage {
	my $script = basename($0);

	print "usage: $script <seq|peak> (command arguments)\n\n";
	print "Command: seq -  Removes the first base of each sequenced read and modifies\n";
	print "                the header to keep track of this base.\n\n";
	print "Command: peak - Tracks the base identified by the 'seq' command to filter\n";
	print "                out identified peaks with low confidence.\n\n";
	print "For commmand arguments / help data for individual commands type '$script <command> --help'.\n";

	exit(0);
}
