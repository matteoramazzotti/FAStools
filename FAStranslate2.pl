#!/usr/bin/perl
if (!$ARGV[0] || $ARGV[0] eq 'help') {
	print "\nUSAGE:FAStransalte.pl [file] [trunc] [perf] [noInc] [-2|3]\n\n";
	print "     file : read DNA from a file (from a pipe if fiel is not specified)\n";
	print "    trunc : truncate proteins after a stop codon (default not)\n";
	print "     perf : do not print proteins with unknown amino acids\n";
	print "    noInc : do not print proteins that appear incomplete (seq%3 != 0)\n";
	print "   -[2|3] : translate form frame 2 or 3 rather than 1\n\n";
	print " ERROR tolerance: by default \"wrong\" proteins are not written unless errors are explicitly allowed by specifying:\n\n";
	print " allowAll : all errors are ignored\n";
	print " allowInc : allow even if the length of CDS is not a multiple of 3\n";
	print " allowMno : allow even if the protein does not start with a Met\n";
	print " allowSno : allow even if the protein does not end with a stop codon\n";
	print "  allowIS : allow even if the protain contains Internal Stop codons\n";
	print " allowAmb : allow even if ambiguous bases are present (the amino acids will be X)\n";
	print "\n";
	exit;
}
&gen_code;
if ($ARGV[0] && -e $ARGV[0]) {
	open($in,"$ARGV[0]") if ($ARGV[0] !~ /\.gz$/);
	open($in,"zcat $ARGV[0]|") if ($ARGV[0] =~ /\.gz$/);
	shift @ARGV;
} else {
	$in = *STDIN;
}
%opts = map {$_ => 1} @ARGV;
$cnt = 0;
while($line = <$in>) {
	chomp $line;
	next if ($line !~ /\w/);
	if ($line =~ /^>/) {
		if ($seq) {
			$exit = do_translate($seq);
			$cnt++ if ($exit);
			print STDERR "\r$cnt" if ($cnt % 100 == 0);
		}
		$name = $line;
		$seq = '';
#		push(@names,$name);
	}
	else {
		$seq .= uc($line);
	}
}
$exit = do_translate($seq);
$cnt++ if ($exit);
print STDERR "\r$cnt" if ($cnt % 100 == 0);

print "--- DEBUG MODE ---\npress [Enter] at every stop\n",join "\n", keys %opts,"\n"  if ($opts{'debug'}); 
<STDIN> if ($opts{'debug'}); 

#trunc: truncate after a stop
#perf: do not print proteins with unknown amino acids
#noInc: do not print proteins that appear incomplete (seq%3 != 0)
#-2: translate from frame 2
#-3: translate from frame 3
$nonprint = 0;
sub do_translate {
	($prot,$err,$herr,$aerr) = translate($seq,$name);
	%err = %$herr;
	$prot =~ s/\*$// if (!$opts{'keepStop'}); #default remove last *
	$prot =~ s/\*.*// if ($opts{'trunc'}); #option to trim at the first internal stop if any

	#by default any error causes translation not to be written, unless errors are explicitly switched off
	#The error types are XXX = 
	#  Inc (length(cod)%3 !=0)
	#  Mno (no Met at start)
	#  Sno (no Stop at last)
	#  IS (internal stop codons)
	#  Amb (ambiguous codon/s)
	#Any error type can be disabled (so, the protein is written) at command line by setting an "allowXXX"
	#All errors can be ignored with allowAll

	$no = 0;
	if (!$opts{"allowAll"}) {
		foreach $k (keys %err) {
			$no++ if ($err{$k} && !$opts{"allow".$k});
			#print STDERR "$k $err{$k}\n";
		}
	}
	if ($opts{'debug'}) {
		foreach $k (keys %err) {
			print STDERR "ERROR $k -> $err{$k}\n";
		}
		if (($opts{'onlyErr'} && $err) || !$opts{'onlyErr'}) {
			print STDERR "$n: $err\n";
			print STDERR $seq{$n},"\n" if (!$opts{'noCod'} || !$opts{'noSeq'});
			print STDERR $prot,"\n" if (!$opts{'noProt'} || !$opts{'noSeq'});
			<STDIN>;
		}
	}
	$nonprint++ if ($no > 0);
	print LOG "$n\n" if ($no > 0 && $opts{'log'});
	return(0) if ($no > 0);
	print "$name\n$prot\n" if (!$no && !$opts{'debug'});
	return(1) if (!$no);
}
print STDERR "\n\n",$nonprint, "proteins not written due to errors\n";
close LOG;

sub translate {
	my $seq = shift;
	my $name = shift;
	my $error = '';
	my %error = ();
	my $prot = '';
	$seq =~ s/^.// if $opts{'-2'};
	$seq =~ s/^..// if $opts{'-3'};
	if (length($seq)%3 != 0) {
#		print STDERR "$name: seems incomplete...\n";
		$error .= "Inc";
		$error{'Inc'} = 1;
#		return 'no';
	}
	my $cnt = 0;
	my $stop = 0;
	my $amb = 0;
	while($seq =~ /(...)/g) {
		$cnt++;
		$cod = $1; #updated 16/4/200
		$prot .= $res{$cod} if ($res{$cod});
#		print STDERR "$name: unknown codon $1\n" if (!$res{$1});
#		return 'no' if (!$res{$1});
		$prot .= 'X' if (!$res{$cod});
		$amb++ if (!$res{$cod} || $cod =~ /^ATGC/); #updated 16/4/200
		$stop++ if ($res{$cod} eq '*');
	}
	if ($amb > 0) {
		$error{'Amb'} = $amb;
		$error .= $amb."Amb";;
	}
	if ($stop > 1 && $prot =~ /\*$/) {
#		print STDERR "$name: internal stop codons\n";
		$error{'IS'} = $stop;
		$error .= $stop."IS";;
	}
	if ($prot !~ /^M/) {
#		print STDERR "$name: no Met at start\n";
		$error{'Mno'} = 1;
		$error .= "Mno";
	}
#	return 'no' if ($prot !~ /^M/ && $opts{'-noM'});
	if ($prot !~ /\*$/) {
#		print STDERR "$name: no stop at end\n";
		$error{'Sno'} = 1;
		$error .= "Sno";
	}
#	return 'no' if ($prot !~ /\*$/ && $opts{'-noEnd'});
#	return 'no' if ($stop > 1 && $prot =~ /\*$/ && $opts{'-noInt'});
	return($prot,$error,\%error);
}

sub gen_code {
(%res) = (
'GGG' => 'G',
'GGA' => 'G',
'GGT' => 'G',
'GGC' => 'G',
'GGM' => 'G',
'GGR' => 'G',
'GGW' => 'G',
'GGS' => 'G',
'GGY' => 'G',
'GGK' => 'G',
'GGV' => 'G',
'GGH' => 'G',
'GGD' => 'G',
'GGB' => 'G',
'GGN' => 'G',

'GAG' => 'E',
'GAA' => 'E',
'GAR' => 'E',

'GAT' => 'D',
'GAC' => 'D',
'GAY' => 'D',

'GTG' => 'V',
'GTA' => 'V',
'GTT' => 'V',
'GTC' => 'V',
'GTM' => 'V',
'GTR' => 'V',
'GTW' => 'V',
'GTS' => 'V',
'GTY' => 'V',
'GTK' => 'V',
'GTV' => 'V',
'GTH' => 'V',
'GTD' => 'V',
'GTB' => 'V',
'GTN' => 'V',

'GCG' => 'A',
'GCA' => 'A',
'GCT' => 'A',
'GCC' => 'A',
'GCM' => 'A',
'GCR' => 'A',
'GCW' => 'A',
'GCS' => 'A',
'GCY' => 'A',
'GCK' => 'A',
'GCV' => 'A',
'GCH' => 'A',
'GCD' => 'A',
'GCB' => 'A',
'GCN' => 'A',

'AGG' => 'R',
'AGA' => 'R',
'AGR' => 'R',

'AGT' => 'S',
'AGC' => 'S',
'AGY' => 'S',

'AAG' => 'K',
'AAA' => 'K',
'AAR' => 'K',

'AAT' => 'N',
'AAC' => 'N',
'AAY' => 'N',

'ATA' => 'I',
'ATT' => 'I',
'ATC' => 'I',
'AT[ATC]' => 'I',
'ATM' => 'A',
'ATW' => 'A',
'ATY' => 'A',
'ATH' => 'A',

'ATG' => 'M',

'ACG' => 'T',
'ACA' => 'T',
'ACT' => 'T',
'ACC' => 'T',
'ACM' => 'T',
'ACR' => 'T',
'ACW' => 'T',
'ACS' => 'T',
'ACY' => 'T',
'ACK' => 'T',
'ACV' => 'T',
'ACH' => 'T',
'ACD' => 'T',
'ACB' => 'T',
'ACN' => 'T',

'TGA' => '*',

'TGT' => 'C',
'TGC' => 'C',
'TGY' => 'C',

'TGG' => 'W',

'TAG' => '*',
'TAA' => '*',
'TAR' => '*',

'TAT' => 'Y',
'TAC' => 'Y',
'TAY' => 'Y',

'TTG' => 'L',
'TTA' => 'L',
'TTR' => 'L',

'TTT' => 'F',
'TTC' => 'F',
'TTY' => 'F',

'TCG' => 'S',
'TCA' => 'S',
'TCT' => 'S',
'TCC' => 'S',
'TCM' => 'S',
'TCR' => 'S',
'TCW' => 'S',
'TCS' => 'S',
'TCS' => 'S',
'TCY' => 'S',
'TCK' => 'S',
'TCV' => 'S',
'TCH' => 'S',
'TCD' => 'S',
'TCB' => 'S',
'TCN' => 'S',

'CGG' => 'R',
'CGA' => 'R',
'CGT' => 'R',
'CGC' => 'R',
'CGM' => 'S',
'CGR' => 'S',
'CGW' => 'S',
'CGS' => 'S',
'CGS' => 'S',
'CGY' => 'S',
'CGK' => 'S',
'CGV' => 'S',
'CGH' => 'S',
'CGD' => 'S',
'CGB' => 'S',
'CGN' => 'S',

'CAG' => 'Q',
'CAA' => 'Q',
'CAR' => 'Q',

'CAT' => 'H',
'CAC' => 'H',
'CAY' => 'H',

'CTG' => 'L',
'CTA' => 'L',
'CTT' => 'L',
'CTC' => 'L',
'CTM' => 'L',
'CTR' => 'L',
'CTW' => 'L',
'CTS' => 'L',
'CTS' => 'L',
'CTY' => 'L',
'CTK' => 'L',
'CTV' => 'L',
'CTH' => 'L',
'CTD' => 'L',
'CTB' => 'L',
'CTN' => 'L',

'CCG' => 'P',
'CCA' => 'P',
'CCT' => 'P',
'CCC' => 'P',
'CCM' => 'P',
'CCR' => 'P',
'CCW' => 'P',
'CCS' => 'P',
'CCS' => 'P',
'CCY' => 'P',
'CCK' => 'P',
'CCV' => 'P',
'CCH' => 'P',
'CCD' => 'P',
'CCB' => 'P',
'CCN' => 'P',
);
}
