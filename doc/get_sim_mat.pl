my %diff;
my @patlist;
my $listfile = $ARGV[0];
open FP, "$listfile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my $pat1 = $a[0];
    @patlist = (@patlist, $pat1);
    open FP1, "$pat1\_sim";
    while(<FP1>) {
	chomp();
	my @b = split("\t");
	my $pat2 = $b[1];
	$diff{$pat1}{$pat2} = $b[2];
	#print "$pat1\t$pat2\n";
    }
    close FP1;
}
close FP;


open OUT, ">sim_mat";
my %diff_sym;
foreach my $i (0..$#patlist) {
    my $pat1 = $patlist[$i];
    foreach my $j (0..$#patlist) {
	my $pat2 = $patlist[$j];
	if($diff{$pat1}{$pat2} > 0) {
	    if($diff{$pat2}{$pat1} > 0) {
		$diff_sym{$pat1}{$pat2} = ($diff{$pat1}{$pat2} + $diff{$pat2}{$pat1}) / 2;
	    } else {
		$diff_sym{$pat1}{$pat2} = $diff{$pat1}{$pat2};
	    }
	} else {
	    if($diff{$pat2}{$pat1} > 0) {
		$diff_sym{$pat1}{$pat2} = $diff{$pat2}{$pat1};
	    } else {
		$diff_sym{$pat1}{$pat2} = "NA";
	    }
	} 
	    
	print OUT "$pat1\t$pat2\t$diff_sym{$pat1}{$pat2}\n";
    }
}
close OUT;
