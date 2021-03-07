#!/usr/bin/perl
#(c) matteo.ramazzotti@unifi.it

if (!$ARGV[1]) {
print STDERR "USAGE: FASpcr.pl [file] forward reverse diff len tol\n\n"; 
print STDERR " file: a fastq file (must be omitted with piped data)\n"; 
print STDERR "forward: 5'-sequence-3'\n"; 
print STDERR "reverse: 5'-sequence-3'\n"; 
print STDERR "   diff: the total number of mismatched allowed primers\n";
print STDERR "    len: an expected length (tol must be specified, too)\n";
print STDERR "    tol: a percent of acceptable length (len) variation\n";
exit;
}
$debug = 0;
if (-e $ARGV[0]) {
	$file = shift @ARGV;
	open($in,$file) if ($file !~ /\.gz$/);
	open($in,"zcat $file|") if ($file =~ /\.gz$/);
} else {
	$in = *STDIN;
}
$fwd = $ARGV[0];
$rev = $ARGV[1];
$mis = $ARGV[2];
if ($ARGV[3]) {
	$ARGV[4] = 0 if (!$ARGV[4]);
	$minlen = int($ARGV[3]-$ARGV[3]/100*$ARGV[4]);
	$maxlen = int($ARGV[3]+$ARGV[3]/100*$ARGV[4]);
}
#&degen;
#print "1. $fwd - $rev\n";
$rev = reverse($rev);
$rev =~ tr/ATGC/TACG/;
#print "2. $fwd - $rev\n";
#$rev =~ s/(\w)/$comp{$1}/gi;
#print "3. $fwd - $rev\n";
#exit;

$lenf = length($fwd);
$lenr = length($rev);

print STDERR "Searching for amplicons\n";
print STDERR " - starting with $fwd ($lenf bp)\n - ending with $rev ($lenr bp)\n";
print STDERR " - with no more than $mis mismatches in each primer\n";
print STDERR " - with a size between $minlen and $maxlen bp\n\n" if ($ARGV[3]);

$cnt = 0;
$tot = 0;
while($line = <$in>) {
	chomp $line;
	if ($line =~ /^>/) {
		$tot++;
		if ($seq) {
			$ampli = ampli($seq);
			print STDERR "---$ampli\n" if ($ampli && $debug);
			print STDERR "---NONE\n" if (!$ampli && $debug);
			print $name,"\n",$ampli,"\n" if ($ampli);
			$cnt++ if ($ampli);
			print STDERR "\r$cnt/$tot" if ($cnt%100 == 0);
		}
		#print STDERR $line,"\n";
		#<STDIN>;
		$seq = '';
		$name = $line;
	} else {
		$seq .= $line
	}
}
$cnt++ if ($ampli);
$ampli = ampli($seq);
print $name,"\n",$ampli,"\n" if ($ampli);
$cnt++ if ($ampli);

print STDERR "\n\n$cnt / $tot amplicon/s formed\n";

sub ampli {
	my $ampli = '';
	my $edit;
	$pos = 0;
	$totf = 0;
	while (1) {
		$part = substr($seq,$pos,$lenf);
		#hamming distance betwen two strings
		$edit = ($part ^ $fwd) =~ tr/\001-\255//;
		if ($edit <= $mis) {
			$posf = $pos;
			$totf++;
			print "$fwd\n$part\nE:$edit\nPOS:$pos\n" if ($debug);
			return(0) if ($totf > 1);
		}
		$pos++;
		last if ($pos == length($seq));
	}
	$pos = 0;
	$totr = 0;
	while (1) {
		$part = substr($seq,$pos,$lenr);
		#hamming distance betwen two strings
		$edit = ($part ^ $rev) =~ tr/\001-\255//;
		#sleep(1);
		if ($edit <= $mis) {
			$posr = $pos if ($edit <= $mis);
			$totr++ if ($edit <= $mis);
			print "$rev\n$part\nE:$edit\nPOS:$pos\n" if ($debug);
			return(0) if ($totr > 1);
		}
		$pos++;
		last if ($pos == length($seq));
	}
	$ampli = substr($seq,$posf,$posr-$posf+length($rev));
	if ($debug) {
		print "\nSEGMENT $posf - $posr, size ",($posr-$posf),"\nAMPLI LEN:",length($ampli),"\n$ampli\n";
		sleep(3);
	}
	return(0) if ($ARGV[3] && (length($ampli) < $minlen || length($ampli) > $maxlen));
	return($ampli);
}

sub degen {
(%deg) = (
"A" => "A",
"T" => "T",
"C" => "C",
"G" => "G",
"R" => "[AG]",
"Y" => "[CT]",
"S" => "[CG]",
"W" => "[AT]",
"K" => "[GT]",
"M" => "[AC]",
"B" => "[CGT]",
"D" => "[AGT]",
"H" => "[ACT]",
"V" => "[ACG]",
"N" => "[ACGT]",
);
foreach $k (keys %deg) {
	$deg{lc($k)} = lc($deg{$k});
}
(%comp) = (
"A" => "T",
"T" => "A",
"C" => "G",
"G" => "C",
"R" => "Y",
"Y" => "R",
"S" => "W",
"W" => "S",
"K" => "M",
"M" => "K",
"B" => "V",
"D" => "H",
"H" => "D",
"V" => "B",
"N" => "N",
);
}
