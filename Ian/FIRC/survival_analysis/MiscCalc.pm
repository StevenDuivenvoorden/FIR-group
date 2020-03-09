package Own::Misc::MiscCalc;

use strict;
use Exporter;

use Sort::Fields;

our $VERSION = 1.0;

our @ISA = qw(Exporter);
our @EXPORT = qw(&distcalc &lininterp &medianpluserror &remainder &sign &valinvec &vecmax &vecmin &veclength &vecmean &vecstdv &vecsum &weightedperc);
our %EXPORT_TAGS = (
	Vectors => [ qw(valinvec vecmax vecmin veclength vecmean vecstdv vecsum) ],
);


#subroutine that calculates the 2-D distance between two points 
sub distcalc {

	my ($x1,$y1,$x2,$y2) = @_;
	my $separation = sqrt(($x1-$x2)**2+($y1-$y2)**2);

	return $separation;

} # distcalc


#subroutine that linearly interpolates between the values of flux belonging to the wavelengths
# nearest to the wavelength specified in the input
sub lininterp {

	my ($lowerx, $upperx, $lowery, $uppery, $var) = @_;

	return $lowery+($uppery-$lowery)/($upperx-$lowerx)*($var-$lowerx);

} #lininterp


#subroutine that computes the remainder resulting from the division of $n by $m
sub remainder {

	my ($n, $m) = @_;
	return ($n/$m-int($n/$m))*$m;
	
} #remainder


#subroutine which implements the fortran function sign(); it is called as &sign($x, $y), the practical effect
# being that &sign($x, $y) has the same absolute value as $x, but it has the same sign as $y; thus the sign of
# $y is transferred to $x. (The case $y = 0 is a little special - it gives &sign($x, $y) always a positive sign.)
sub sign {

	my ($var1, $var2) =  @_;
	
	if ($var2 >= 0){return abs($var1);} else {return -1*abs($var1)}

} #sign


#subroutine that returns the median, the error on the median [calculated as (16%-perc - 84%-perc)/sqrt(N)]
# and an estimate of the spread in the population through an upper & lower percentile; the
# subroutine needs PDL to be installed & is called by supplying a data piddle &  percentiles, e.g.:
# ($median, $dmedian, $lowerperc, $upperperc) = &medianpluserror($datapiddle, .25, .75);
sub medianpluserror {

	eval "require PDL";
	if ($@){ die("\nCannot locate package PDL. Bailing out...\n"); } else { PDL->import(); }

	my ($indata, $lowerperc, $upperperc) = @_;
	my ($sortedata, $median, $medianerror, $lowerpercval, $upperpercval, $perc16, $perc84);

	$sortedata = qsort($indata);
	$median = index($sortedata, floor(.5*(nelem($indata)-1))) + (.5*(nelem($indata)-1)-floor(.5*(nelem($indata)-1)))*
		(index($sortedata, ceil(.5*(nelem($indata)-1)))-index($sortedata, floor(.5*(nelem($indata)-1))));

	$lowerpercval = index($sortedata, floor($lowerperc*(nelem($indata)-1))) + ($lowerperc*(nelem($indata)-1)-floor($lowerperc*(nelem($indata)-1)))*
		(index($sortedata, ceil($lowerperc*(nelem($indata)-1)))-index($sortedata, floor($lowerperc*(nelem($indata)-1))));
	$upperpercval = index($sortedata, floor($upperperc*(nelem($indata)-1))) + ($upperperc*(nelem($indata)-1)-floor($upperperc*(nelem($indata)-1)))*
		(index($sortedata, ceil($upperperc*(nelem($indata)-1)))-index($sortedata, floor($upperperc*(nelem($indata)-1))));

	$perc16 = index($sortedata, floor(.16*(nelem($indata)-1))) + (.16*(nelem($indata)-1)-floor(.16*(nelem($indata)-1)))*
	(index($sortedata, ceil(.16*(nelem($indata)-1)))-index($sortedata, floor(.16*(nelem($indata)-1))));
	$perc84 = index($sortedata, floor(.84*(nelem($indata)-1))) + (.84*(nelem($indata)-1)-floor(.84*(nelem($indata)-1)))*
		(index($sortedata, ceil(.84*(nelem($indata)-1)))-index($sortedata, floor(.84*(nelem($indata)-1))));

	$medianerror = ($perc16-$perc84)/2/sqrt(nelem($indata));

	return ($median, $medianerror, $lowerpercval, $upperpercval);

} #medianpluserror


#subtroutine that searches an array for a specific scalar or string and returns the number
# of occurrences 
sub valinvec {

	my ($target, $list) = @_;
	
	my $tracker = 0;
	for (my $i=0; $i<@$list; $i++){
		my $val = $list->[$i];
		if (($val == $target)||($val eq $target)){ $tracker++; }
	}

	return $tracker;

} #valinvec


#subroutine that computes the minimum of an array or individually passed numbers
sub vecmax {

	my @numbers=@_;
	my @sorted_numbers=sort {$b <=> $a} @numbers;

	return $sorted_numbers[0];

} #vecmax


#subroutine that computes the minimum of an array or individually passed numbers
sub vecmin {

	my @numbers=@_;
	my @sorted_numbers=sort {$a <=> $b} @numbers;

	return $sorted_numbers[0];

} #vecmin


#subroutine that calculates the length of a vector
sub veclength {

	my @invector=@_;
	my $sumsq=0;
	
	foreach my $elem (@invector) { $sumsq += $elem**2; }
	$sumsq=sqrt($sumsq);

	return($sumsq);

} #veclength


#subroutine that calculates the mean of the entries of a vector of arbitrary length
sub vecmean {

	my @invector=@_;
	my $veclength=@invector;
	
	my $sum=0;
	my $mean=0;

	my $step=0;
	while ($step<$veclength){
			
		$sum=$sum+$invector[$step];
		$step=$step+1;
		
	} #counter
		
	$mean=$sum/$veclength;

	return $mean;

} #vecmean


#subroutine that calculates the stdv of the entries of a vector of arbitrary length
# (formula of the stdv according to p.13 of the 'Statistics Manual' by Crow, Davis & Maxfield)
sub vecstdv {

	my @invector=@_;
	my $veclength=@invector;
	
	my $sum=0;
	my $sumsq=0;
	my $sigma=0;

	my $step=0;
	while ($step<$veclength){
			
		$sum=$sum+$invector[$step];
		$sumsq=$sumsq+($invector[$step])**2;
		$step=$step+1;
		
	} #counter
		
	$sigma=sqrt(abs((($veclength*$sumsq)-$sum**2)/($veclength*($veclength-1))));

	return $sigma;

} #vecstdv


#subroutine that adds up all array elements in an array of arbitrary length
sub vecsum {

	my @invector=@_;
	my $sum=0;
	
	foreach my $elem (@invector) { $sum += $elem; }
	
	return($sum);

} #vecsum


#subroutine that computes a weighted percentile; accepts a piddle of values a piddle of weights
# and the percentile of interest
sub weightedperc {

	eval "require PDL";
	if ($@){ die("\nCannot locate package PDL. Bailing out...\n"); } else { PDL->import(); }

	my ($indata, $inwghts, $perc) = @_;

	#merge data into lines to run the field-sorting module on it; sort data
	# according to increasing values of data points 
	my @lines = ();
	for (my $i=0; $i<nelem($indata); $i++){ $lines[$i] = $indata->at($i)."\t".$inwghts->at($i)."\n"; }
	undef($indata);
	undef($inwghts);

	my @newlines = fieldsort ['1n'], @lines;
	chomp(@newlines);
	undef(@lines);

	#extract new sorted data- & weight-arrays from the sorted table
	my @sorteddata = ();
	my @sortedwghts = ();
	for (my $i=0; $i<=$#newlines; $i++){ ($sorteddata[$i], $sortedwghts[$i]) = split(/\t/,$newlines[$i]); }
	undef(@newlines);

	#calculate the weighted median; define it to lie at half the sum of all
	# weights, then step through the data until this threshold is exceeded,...
	my $separval = sum(pdl(@sortedwghts))*$perc;
	undef($perc);

	my $outval = 0;
	my $wghtsum = 0;
	for (my $i=0; $i<=$#sorteddata; $i++){
	
		$wghtsum += $sortedwghts[$i];
	
		if ($wghtsum > $separval){
	
			if ($i == 0){
				$outval = $sorteddata[0];
				print STDOUT "Non-regular calculation of percentile (lower end of data range)!\n";
				last;
			}
	
			#calculate the precise value of the median by interpolating between the data
			# points
			$outval = $sorteddata[$i-1]+($sorteddata[$i]-$sorteddata[$i-1])/$sortedwghts[$i]*($separval-($wghtsum-$sortedwghts[$i]));
			
			last;
	
		}
		
		if ($i == $#sorteddata){
				$outval = $sorteddata[$i];
				print STDOUT "Non-regular calculation of percentile (upper end of data range)!\n";
		}
		
	}

	return($outval);

	undef(@sorteddata);
	undef(@sortedwghts);
	undef($separval);
	undef($wghtsum);
	undef($outval);

} #weightedperc


1;              # Modules must return true.
