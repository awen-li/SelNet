#!/usr/bin/perl
#
# subset number of individuals
#

use Math::Random qw(random_binomial);

$N2 = 100;

foreach $in (@ARGV){
	$out = "sub$N2"."$in";
	open (IN, $in) or die "failed to read\n";
	open (OUT, "> $out") or die "failed to write\n";
	$hdr = <IN>;
	print OUT "$hdr";
	while (<IN>){
		chomp;
		@A = split(" ",$_);
		@B = ();
		foreach $snp (@A){
			@cnt = split(",",$snp);
			$p = $cnt[0]/($cnt[1] + $cnt[0]);
			$x = random_binomial(1, $N2, $p);
			$y = $N2 - $x;
			push (@B,"$x,$y");
		}
		$b = join(" ",@B);
		print OUT "$b\n";
	}
	close(IN);
	close(OUT);
}

