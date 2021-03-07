#!/usr/bin/perl
#(c) matteo.ramazzotti@unifi.it

if (!$ARGV[0]) {
	print "USAGE: FASfreqmat.pl alnfile c|f|string[vh] outfile\n\n";
	print "      c: count\n";
	print "      f: frequency\n";
	print " string: output variants as string, not as matrix\n";
	print "      v: output string in rows\n";
	print "      h: output string in columns\n\n";
	exit; 
}

if ($ARGV[1] eq $ARGV[$#ARGV]) {
	print STDERR "ERROR: input and output file are the same!\n\n";
	exit;
}

open(IN,$ARGV[0]) or die "cannot open $ARGV[0]";
while ($line = <IN>) {
	chomp $line;
	if ($line =~ />/) {
		$name = $line;
		push(@names,$name);
	} else {
		$seq{$name} .= $line;
	}
}
close IN;
$totseq = scalar keys %seq;

$outfile = pop @ARGV;
open(OUT,">$outfile");

%cnt = ();
foreach $pos (0..length($seq{$names[0]})-1) {
	$col = '';
	$gaps = 0;
	foreach $n (@names){
		$a = substr($seq{$n},$pos,1);
		$col .=	$a;
	}
	foreach $aa (split (//,$col)) {
		$all{$aa} = 1;
		$cnt{$aa."@".$pos}++;
	}
}
@head = split(//,'-ACDEFGHIKLMNPQRSTVWY');

#both stringv or stringh start here
goto STRING if ($ARGV[1] =~ /string/);

print OUT "Pos\tRef\t",join "\t",@head;
print OUT "\n";
@ind = split (//,$seq{$names[0]});
foreach $pos (0..length($seq{$names[0]})-1) {
	print OUT $pos+1,"\t",$ind[$pos];
	foreach $aa (@head) {
		$val = $cnt{$aa."@".$pos};
		$val = $val/$totseq if($ARGV[1] && $ARGV[1] eq 'f');
		print OUT "\t";
		print OUT $val if ($val);
		print OUT "0" if (!$val);
	}
	print OUT "\n";
}

STRING:
@ind = split (//,$seq{$names[0]});
$stringv = ''; 
$mutmax = 0;
foreach $pos (0..$#ind) {
	$stringv .= ($pos+1)."\t".$ind[$pos];
	$muttot = 0;
	foreach $aa (@head) {
		if ($cnt{$aa."@".$pos} && $aa ne $ind[$pos]) {
			$stringv .= "\t".$aa;
			$muttot++;
		}
	}
	$mutmax = $muttot if ($muttot > $mutmax);
	$stringv .= "\n";
}
print OUT $stringv;
close OUT;

if ($ARGV[1] =~ /stringh/) {
	print STDERR "Reformtting in H ";
	open(OUT,">tmp");	
	foreach $c (1..$mutmax+2) {
		$col = `cut -f$c $outfile`;
		$col =~ s/\n/\t/g;
		print OUT $col,"\n";
	}
	`mv tmp $outfile`;
	print STDERR "done.\n";
}

if ($ARGV[2]) {
	open(IN,$outfile);
	@in = <IN>;
	close IN;
	foreach $l (0..$#in) {
		chomp $in[$l];
		@{$line{$l}} = split (/\t/,$in[$l]);
	}
	$s = 0;
	$e = $ARGV[2]-1;
	$parts = int(scalar @{$line{0}}/$ARGV[2]);
	open(OUT,">tmp");
	foreach $l (0..(scalar keys %line)-1) {
		print OUT join "\t",@{$line{$l}}[$s..$e],"\n";
	}
	print OUT "\n";
	foreach $p (1..$parts) {
		foreach $c ($s..$ARGV[2]) {
			$s += $ARGV[2];
			$e += $ARGV[2];
			foreach $l (0..(scalar keys %line)-1) {
				print OUT join "\t",@{$line{$l}}[$s..$e],"\n";
			}
		}
		print OUT "\n";
	}
	`mv tmp $outfile`;
}
