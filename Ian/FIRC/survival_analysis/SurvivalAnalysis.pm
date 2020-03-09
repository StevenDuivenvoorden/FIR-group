package Own::StatisticsTools::SurvivalAnalysis;

use strict;
use Exporter;

use PDL;
use PDL::Fit::LM;
use PDL::MatrixOps;
use PDL::NiceSlice;
use PDL::GSL::INTEG;
use Sort::Fields;

use Own::Misc::MiscCalc;
use Own::Misc::Utilities;
use Own::Misc::FitLM_Funcs;


our $VERSION = 1.0;

our @ISA = qw(Exporter);
our @EXPORT = qw(&BJregress &censotransf &distribtest_left &distribtest_right &doubcense_distrib &kaplan_meier &plotdistrprep_double &plotdistrprep_left &plotdistrprep_right &sortforkaplan &survfct_median &survfct_pct);

use constant PI => (4*atan2(1,1));


#wrapper subroutine for Buckley-James regression; redirects to the pertient subroutines,
# depending on the kind of censoring present in the data
sub BJregress {

	my ($inxxarr, $inyyarr, $inflagarr, $inwghtarr, $tol, $maxiter, $verboseopt) = @_;

	my @flagarr = @{ $inflagarr };
	if ((&vecmin(@flagarr) == -1) && (&vecmax(@flagarr) == 0)){

		for (my $i=0; $i<=$#flagarr; $i++){ $flagarr[$i] += 1; }
		return (&BJregress_left($inxxarr, $inyyarr, \@flagarr, $inwghtarr, $tol, $maxiter, $verboseopt));

	} elsif ((&vecmin(@flagarr) == 0) && (&vecmax(@flagarr) == 1)){

		for (my $i=0; $i<=$#flagarr; $i++){ $flagarr[$i] = -1*($flagarr[$i]-1); }
		return (&BJregress_right($inxxarr, $inyyarr, \@flagarr, $inwghtarr, $tol, $maxiter, $verboseopt));

	} else { die("'BJregress' is not capable of handling the requested censoring scheme!\nBailing out...\n\n"); }

	undef($inxxarr);
	undef($inyyarr);
	undef($inflagarr);
	undef($inwghtarr);
	undef($tol);
	undef($maxiter);
	undef($verboseopt);
	undef(@flagarr);

} #BJregress



#subroutine that does a Buckley-James (EM algorithm) method regression according to Isobe+ '86 (p. 496ff)
# for left-censored data (upper limits)
sub BJregress_left {

	my ($inxxarr, $inyyarr, $inflagarr, $inwghtarr, $tol, $maxiter, $verboseopt) = @_;

#transform to a right-censored data set...
	my $inflect = &vecmax(@{ $inyyarr });
	my @yyarr = &censotransf(\@{ $inyyarr }, $inflect);
	undef($inyyarr);

#... & do Buckley-James regression
	my (@suboutpt) = &BJregress_right($inxxarr, \@yyarr, $inflagarr, $inwghtarr, $tol, $maxiter, $verboseopt);
	my @slope_r = @{ $suboutpt[0] };
	my $yintercept_r = $suboutpt[1];
	my $stddev = $suboutpt[2];
	my $niter = $suboutpt[3];
	undef(@suboutpt);

#transform output back to left-censoring where necessary
	my @slope = (-$slope_r[0], $slope_r[1]);
	my $yintercept = $inflect-$yintercept_r;

	return (\@slope, $yintercept, $stddev, $niter);

	undef($inxxarr);
	undef(@yyarr);
	undef($inflagarr);
	undef($inwghtarr);
	undef($tol);
	undef($maxiter);
	undef($verboseopt);
	undef($inflect);

	undef(@slope);
	undef($yintercept);
	undef($stddev);
	undef($niter);

} #BJregress_left


#subroutine that does a Buckley-James (EM algorithm) method regression according to Isobe+ '86 (p. 496ff)
# for right-censored data (lower limits)
sub BJregress_right {

	my ($inxxarr, $inyyarr, $inflagarr, $inwghtarr, $tol, $maxiter, $verboseopt) = @_;

	my $xxarr = pdl(@{ $inxxarr });
	my $yyarr= pdl(@{ $inyyarr });
	my $yyarr_rev = pdl($yyarr);
	my $flagarr = pdl(@{ $inflagarr });
	my $wghtarr = pdl(@{ $inwghtarr });

	undef($inxxarr);
	undef($inyyarr);
	undef($inflagarr);
	undef($inwghtarr);

	my $verbosedir = "BJregress_logfiles";
	if (($verboseopt eq "y")&&(!-e $verbosedir)){ system("mkdir $verbosedir"); }

#derive initial estimate of best-fitting trend line by treating censored points as detections
	my ($yf, $pf, $cf, $if) = lmfit $xxarr, $yyarr, $wghtarr, \&linefit, pdl(0, 0), { Maxiter => 300, Eps => 1e-3 };
	my $prevslope = $pf->at(0);
	my $prevyintercept = $pf->at(1);
	my ($anteprevslope, $anteprevyintercept);

	undef($yf);
	undef($pf);
	undef($cf);
	undef($if);


#iterate until convergence is achieved or the maximal number of iterations has been reached
	my @slope = (9e+99, 9e+99);
	my $yintercept = 9e+99;
	my $stddev = 9e+99;
	my $niter = 9e+99;

	for (my $iter=1; $iter<=$maxiter; $iter++){

#compute the arguments of the KM estimator (Y_i - [y-intercept] - [slope]*X_i)
		my @KMarg = list($yyarr - $prevyintercept - $prevslope*$xxarr);
#		my @KMarg = list($yyarr - $prevslope*$xxarr);
		my @ordKMarg = list(qsort(pdl(@KMarg)));

		my @ordinds = list(qsorti(pdl(@KMarg)));
		my @KMarg_ordpos = ();
		@KMarg_ordpos[@ordinds] = list(sequence(nelem($xxarr))); #position of @KMarg-values in array @ordKMarg

#define largest residual to be uncensored for computational convenience if necessary
		my @censeflags = list($flagarr);
		$censeflags[maximum_ind(pdl(@KMarg))] = 1;
		my @ordflags = @censeflags[@ordinds];

#derive the KM-distribution function
		my (@suboutpt) = &kaplan_meier_right(\@KMarg, \@censeflags);
		my @survxx = @{ $suboutpt[0] };
		my @survyy = @{ $suboutpt[1] };
		undef(@suboutpt);

		#some 'surgery' in case KMarg contains triple or more frequent entries
		if ($#survxx != $#KMarg+1){
			for (my $i=0; $i<=$#KMarg; $i++){

				if (sprintf("%.10f", $survxx[$i]) != sprintf("%.10f", $ordKMarg[$i])){

					push(@survxx, 0);
					push(@survyy, 0);
					@survxx[$i+1..$#survxx] = @survxx[$i..$#survxx-1];
					$survxx[$i] = $survxx[$i-1];
					@survyy[$i+1..$#survyy] = @survyy[$i..$#survyy-1];
					$survyy[$i] = $survyy[$i-1];

				}

			}
		}

#compute new Y-values for limits
		for (my $i=0; $i<nelem($yyarr); $i++){

			if ($censeflags[$i] == 0){

				my $KMsum = 0;
				for (my $j=$KMarg_ordpos[$i]+1; $j<nelem($yyarr); $j++){

					unless ($ordflags[$j] == 0){
						my $weight = ($survyy[$j]-$survyy[$j+1])/$survyy[$KMarg_ordpos[$i]];
						$KMsum += $weight*$survxx[$j];
						undef($weight);
					}

				} #$j

				$KMsum += $prevyintercept + $prevslope*$xxarr->at($i);
				$yyarr_rev($i) .= $KMsum;

				undef($KMsum);

			}

		} #$i

#=head
		$slope[0] = sum($yyarr_rev*($xxarr-avg($xxarr)))/sum(($xxarr-avg($xxarr))**2);
		$yintercept = avg($yyarr_rev) - $slope[0]*avg($xxarr);
#=cut
=head
		my ($yf, $pf, $cf, $if) = lmfit $xxarr, $yyarr_rev, $wghtarr, \&linefit, pdl($prevslope, $prevyintercept), { Maxiter => 300, Eps => 1e-3 };
		$slope[0] = $pf->at(0);
		$yintercept = $pf->at(1);

		undef($yf);
		undef($pf);
		undef($cf);
		undef($if);
=cut

		$stddev = sqrt(sum(($yyarr->where($flagarr==1)-avg($yyarr->where($flagarr==1)) - $slope[0]*($xxarr->where($flagarr==1)-avg($xxarr->where($flagarr==1))))**2)/(nelem($flagarr->where($flagarr==1))-2));


#if requested, print line parameters and "corrected" Y-vector (where censored values are
# updated according to prob. distrib.) to file for each iteration
		if ($verboseopt eq "y"){

			my $verbosefile = $verbosedir."/iteration.".sprintf("%05d",$iter);
			my $headerstr = "#slope\tinterc.\tstddev\n#".sprintf("%5.3e",$slope[0])."\t".sprintf("%5.3e",$yintercept)."\t".sprintf("%5.3e",$stddev)."\n#\n#xvector\tyvector\tstatus\n";
			wcols "%.5f\t%.5f\t%2d", $xxarr, $yyarr_rev, $flagarr, $verbosefile, { Header => $headerstr };

			undef($headerstr);
			undef($verbosefile);

		}


#check if convergence as been achieved
		my $del_slope = abs($slope[0]-$prevslope);
		my $del_intercept = abs($yintercept-$prevyintercept);
		my $del_tot = sqrt($del_slope**2 + $del_intercept**2);
		undef($del_slope);
		undef($del_intercept);

		unless ($del_tot <= $tol){

			if ($iter == $maxiter){

				my $del_slope = abs($slope[0]-$anteprevslope);
				my $del_intercept = abs($yintercept-$anteprevyintercept);
				my $delosc_tot = sqrt($del_slope**2 + $del_intercept**2);
				undef($del_slope);
				undef($del_intercept);

				if ($delosc_tot <= $tol){

					print STDOUT "\nOscillatory behaviour detected; now calculating average...\n\n";

					$slope[0] = ($slope[0]+$prevslope)/2;
					$yintercept = ($yintercept+$prevyintercept)/2;
					$slope[1] = $stddev/sqrt(sum(($xxarr->where($flagarr==1)-avg($xxarr->where($flagarr==1)))**2));
					$niter = $iter;

				} else {

					print STDOUT "\n'BJregress' did not achieve convergence after the ".$iter." permitted iterations.\n\n";

					my @dummyslope = (-99, -99);
					return (\@dummyslope, -99, -99, $iter);
					undef(@dummyslope);
					
				}

				undef($delosc_tot);

			} else {

				$anteprevslope = $prevslope;
				$anteprevyintercept = $prevyintercept;
				$prevslope = $slope[0];
				$prevyintercept = $yintercept;

			}

		} else {

			$slope[1] = $stddev/sqrt(sum(($xxarr->where($flagarr==1)-avg($xxarr->where($flagarr==1)))**2));
			$niter = $iter;

			if ($verboseopt eq "y"){ system("mv ".$verbosedir."/iteration.".sprintf("%05d",$niter)." ".$verbosedir."/iteration.final"); }

			last;

		}

		undef(@KMarg);
		undef(@ordKMarg);
		undef(@ordinds);
		undef(@KMarg_ordpos);
		undef(@censeflags);
		undef(@survxx);
		undef(@survyy);
		undef(@ordflags);
		undef($del_tot);

	} #$iter


	undef($tol);
	undef($maxiter);
	undef($xxarr);
	undef($yyarr);
	undef($yyarr_rev);
	undef($flagarr);
	undef($wghtarr);
	undef($prevslope);
	undef($prevyintercept);
	undef($anteprevslope);
	undef($anteprevyintercept);
	undef($verbosedir);

	return (\@slope, $yintercept, $stddev, $niter);

	undef(@slope);
	undef($yintercept);
	undef($stddev);
	undef($niter);

} #BJregress_right


#convert left-censored (incl. lower limits) to right-censored data (invert at the largest measurment in left-censored data)
sub censotransf {

	my ($inarray, $inflect) = @_;
	my $array = pdl(@{ $inarray });
	my @outarray = list($inflect-$array);
	
	return(@outarray);

} #censotransf


sub distribtest_left {

	my ($indataarr1, $inflagarr1, $indataarr2, $inflagarr2, $method, $inflect) = @_;

	my @dataarr1 = &censotransf(\@{ $indataarr1 }, $inflect);
	my @flagarr1 = @{ $inflagarr1 };
	my @dataarr2 = &censotransf(\@{ $indataarr2 }, $inflect);
	my @flagarr2 = @{ $inflagarr2 };
	
	my (@suboutpt) = &distribtest_right(\@dataarr1, \@flagarr1, \@dataarr2, \@flagarr2, $method);

	my @testresult = @{ $suboutpt[0] };
	my $alpha = $suboutpt[1];

	return(\@testresult, $alpha);

	undef(@suboutpt);
	undef(@dataarr1);
	undef(@flagarr1);
	undef(@dataarr2);
	undef(@flagarr2);
	undef(@testresult);

} #$distribtest_left


sub distribtest_right {

	my ($indataarr1, $inflagarr1, $indataarr2, $inflagarr2, $method) = @_;

	my @dataarr1 = @{ $indataarr1 };
	my @flagarr1 = @{ $inflagarr1 };
	my @dataarr2 = @{ $indataarr2 };
	my @flagarr2 = @{ $inflagarr2 };

	my @globdata = (@dataarr1, @dataarr2);
	my @globflags = (@flagarr1, @flagarr2);


#sort input data in both data sets individually, as well as in the joint data sets (in 2nd priority sort
# by descending censorship flag value)
	my (@suboutpt) = &sortforkaplan(\@dataarr1, \@flagarr1);
	my @sorteddata1 = @{ $suboutpt[0] };
	my @sortedflags1 = @{ $suboutpt[1] };

	my (@suboutpt) = &sortforkaplan(\@dataarr2, \@flagarr2);
	my @sorteddata2 = @{ $suboutpt[0] };
	my @sortedflags2 = @{ $suboutpt[1] };

	my (@suboutpt) = &sortforkaplan(\@globdata, \@globflags);
	my @sortedglobdata = @{ $suboutpt[0] };
	my @sortedglobflags = @{ $suboutpt[1] };

	undef (@suboutpt);


#retrieve array of distinct measurement values for the combined data sets
	my (@suboutpt) = &gendistinctdata(\@sortedglobdata, \@sortedglobflags);
	my $uniqueglobdatapdl = pdl(@{ $suboutpt[0] });
	my $uniqueglobflagspdl = pdl(@{ $suboutpt[1] });
	my @uniqueglobdata = list($uniqueglobdatapdl->where($uniqueglobflagspdl == 1));
	my @uniqueglobflags = list($uniqueglobflagspdl->where($uniqueglobflagspdl == 1));
	undef (@suboutpt);


#compute the various distribution descriptors
	my @d1 = ();
	my @d2 = ();
	my @d = ();

	my @n1 = ();
	my @n2 = ();
	my @n = ();

	my @m1 = ();
	my @m2 = ();
	my @m = ();


	my $counter = 0;
	for(my $i=0; $i<=$#uniqueglobdata; $i++){

		$d1[$i] = 0;
		$n1[$i] = 0;
		$m1[$i] = 0;

		for (my $j=0; $j<=$#sorteddata1; $j++){
			if (($sorteddata1[$j] == $uniqueglobdata[$i])&&($sortedflags1[$j] == 1)){ $d1[$i]++; }
			if (($sorteddata1[$j] >= $uniqueglobdata[$i])||(($sorteddata1[$j] == $uniqueglobdata[$i])&&($sortedflags1[$j] == 0))){ $n1[$i]++; }
			if (($i+1) <= $#uniqueglobdata){
				if (($sorteddata1[$j] >= $uniqueglobdata[$i])&&($sorteddata1[$j] < $uniqueglobdata[$i+1])&&($sortedflags1[$j] == 0)){ $m1[$i]++; }
			} else {
				if (($sorteddata1[$j] >= $uniqueglobdata[$i])&&($sortedflags1[$j] == 0)){ $m1[$i]++; }
			}
		}

		$d2[$i] = 0;
		$n2[$i] = 0;
		$m2[$i] = 0;

		for (my $j=0; $j<=$#sorteddata2; $j++){
			if (($sorteddata2[$j] == $uniqueglobdata[$i])&&($sortedflags2[$j] == 1)){ $d2[$i]++; }
			if (($sorteddata2[$j] >= $uniqueglobdata[$i])||(($sorteddata2[$j] == $uniqueglobdata[$i])&&($sortedflags2[$j] == 0))){ $n2[$i]++; }
			if (($i+1) <= $#uniqueglobdata){
				if (($sorteddata2[$j] >= $uniqueglobdata[$i])&&($sorteddata2[$j] < $uniqueglobdata[$i+1])&&($sortedflags2[$j] == 0)){ $m2[$i]++; }
			} else {
				if (($sorteddata2[$j] >= $uniqueglobdata[$i])&&($sortedflags2[$j] == 0)){ $m2[$i]++; }
			}
		}

		$d[$i] = $d1[$i] + $d2[$i];
		$n[$i] = $n1[$i] + $n2[$i];
		$m[$i] = $m1[$i] + $m2[$i];
	
	} #$i


#define weighting factors depending on the choice of statistical test to be performed
	my @w = ();
	my @beta = ();
		if ($method eq "Gehan"){
		@w = @n;
	} elsif ($method eq "logrank"){
		@w = list(ones($#uniqueglobdata+1));
	} elsif ($method eq "PP"){

		for (my $i=0; $i<=$#uniqueglobdata; $i++){
			$w[$i] = 1;
			$beta[$i] = 1;
			for (my $j=0; $j<=$i; $j++){
				$w[$i] *= $n[$i]/($n[$i]+1);
				$beta[$i] *= ($n[$i]+1)/($n[$i]+2);
			}
		} #$i

	} else {
		die("\nUnknown statistical test requested!\n Bailing out...\n\n");
	}


#calculate the variables for the statistical testing
	my $rankstat = 0;
	my $drankstatsq = 0;
	for (my $i=0; $i<=$#uniqueglobdata; $i++){

		$rankstat += $w[$i]*($d1[$i] - ($d[$i]*$n1[$i]/$n[$i]));

		if ($method ne "PP"){

			my $numer = $d[$i]*($w[$i]**2)*($n1[$i]/$n[$i])*($n2[$i]/$n[$i])*($n[$i]-$d[$i]);
			my $denom = ($n[$i]-1);

			if (($numer == 0)&&($denom == 0)){
				$drankstatsq += 0;
			} else {
				$drankstatsq += $numer/$denom;
			}

		} else {

			my $gamma = 2*$d2[$i] + $m2[$i];
			my $sumforerr = 0;
			for (my $j=$i+1; $j<=$#uniqueglobdata; $j++){ $sumforerr += $w[$j]*$beta[$j]; }
			$drankstatsq += $w[$i]*(1-$beta[$i])*$gamma - ($beta[$i]-$w[$i])*$gamma*($w[$i]*$gamma + 2*$sumforerr);

		}

	} #$i


#compute at which significance level the distributions differ	
	my $testval = abs($rankstat/sqrt($drankstatsq));
	my @alphaarr = gslinteg_qng(\&normalcurve,-$testval,$testval,0,1e-10, {Warn => 'y'});
	my $alpha = 1 - $alphaarr[0];

	my @testresult = ($rankstat, sqrt($drankstatsq));

	return(\@testresult, $alpha);

	undef(@dataarr1);
	undef(@flagarr1);
	undef(@dataarr2);
	undef(@flagarr2);
	undef(@globdata);
	undef(@globflags);
	undef(@sorteddata1);
	undef(@sortedflags1);
	undef(@sorteddata2);
	undef(@sortedflags2);
	undef(@sortedglobdata);
	undef(@sortedglobflags);
	undef(@uniqueglobdata);
	undef(@uniqueglobflags);
	undef(@d1);
	undef(@d2);
	undef(@d);
	undef(@n1);
	undef(@n2);
	undef(@n);
	undef(@m1);
	undef(@m2);
	undef(@m);
	undef(@w);
	undef(@beta);
	undef(@alphaarr);
	undef(@testresult);

	sub normalcurve {
		my ($x) = @_;
		return (exp(-($x**2)/2)/sqrt(2*PI));
	}

} #distribtest_right


#subroutine that computes the distribution function of a doubly censored data set according to Schmitt & al., 1993, A&A, 277, 114
sub doubcense_distrib {

	my ($inacudata, $insupdata, $ininfdata, $nbins) = @_;

	my $origacudata = pdl(@{ $inacudata });
	my $origsupdata = pdl(@{ $insupdata });
	my $originfdata = pdl(@{ $ininfdata });

#replace limits lying beyond the last measurement with a fake known measurement outside the later plotting range in order
# to ensure that the algorithm returns the correct asymptotic value as data -> +/-infty
	my $nsubsubst = nelem($origsupdata->where($origsupdata<=min($origacudata)));
	my $ninfsubst = nelem($originfdata->where($originfdata>=max($origacudata)));

	my $alldata = $origacudata->append($origsupdata)->append($originfdata);
	my $datarange = max($alldata)-min($alldata);

	my $supdata = $origsupdata->where($origsupdata>min($origacudata));
	my $infdata = $originfdata->where($originfdata<max($origacudata));
	my $acudata;
	if ($nsubsubst != 0){
		$acudata = pdl(ones($nsubsubst)*(min($alldata)-$datarange*.15))->append($origacudata);
	} else {
		$acudata = $origacudata;
	}
	unless ($ninfsubst == 0){ $acudata = $acudata->append(ones($ninfsubst)*(max($alldata)+$datarange*.15)); }


#renormalize data such that it is always larger zero
	my $renormconst = 0;
	my $rangemin = min($alldata)-$datarange*.2;
	if ($rangemin < 0){
		$renormconst = $rangemin;
		$rangemin = 0;
	}
	my $rangemax = max($alldata)+$datarange*.2-$renormconst;
	my $binwidth = sprintf("%.4f", ($rangemax-$rangemin)/$nbins);
	
	$acudata -= $renormconst;
	$supdata -= $renormconst;
	$infdata -= $renormconst;


#generate histograms of limits & accurate measurements
	my ($acuabsc, $acuhisto) = hist($acudata, $rangemin, $rangemax, $binwidth);
	my ($supabsc, $suphisto);
	unless (nelem($supdata) == 0){
		($supabsc, $suphisto) = hist($supdata, $rangemin+$binwidth/2, $rangemax+$binwidth/2, $binwidth);
	} else {
		$suphisto = zeroes(nelem($acuabsc));
	}
	my ($infabsc, $infhisto);
	unless (nelem($infdata) == 0){
		($infabsc, $infhisto) = hist($infdata, $rangemin-$binwidth/2, $rangemax-$binwidth/2, $binwidth);
	} else {
		$infhisto = zeroes(nelem($acuabsc));
	}

	my @n = list($acuhisto);
	my @u = list($suphisto);
	my @l = list($infhisto);


#initialize distribution fct. & start iteration
	my $Ntot = 0;
	for (my $i=0; $i<=$#n; $i++){ $Ntot += ($n[$i]+$u[$i]+$l[$i]); }
	my $fipdl = $acuhisto/$Ntot;
	my @f = list($fipdl);

	my $m = 1;
	my @fdiff = ();
	print STDOUT "Iterating - step\n\t";
	while ($m<=200){

		print STDOUT $m."..";

		for (my $k=0; $k<=$#f; $k++){

			my $uisum = 0;
			for (my $i=$k; $i<=$#f; $i++){
				my $fjsum = 0;
				for (my $j=0; $j<=$i; $j++){ $fjsum += $f[$j]; }
				if ($f[$k] == 0){
					$uisum += 0;
				} else {
					$uisum += $u[$i]*$f[$k]/$Ntot/$fjsum;
				}
			} #$i

			my $lisum = 0;
			for (my $i=0; $i<=$k; $i++){
				my $fjsum = 0;
				for (my $j=$i; $j<=$#f; $j++){ $fjsum += $f[$j]; }
				if ($f[$k] == 0){
					$lisum += 0;
				} else {
					$lisum += $l[$i]*$f[$k]/$Ntot/$fjsum;
				}
			} #$i

			my $oldf = $f[$k];
			$f[$k] = $n[$k]/$Ntot + $uisum + $lisum;
			if (($f[$k]-$oldf) != 0){
				$fdiff[$k] = abs(($f[$k]-$oldf)/$oldf);
			} else {
				$fdiff[$k] = 0;
			}

			undef($uisum);
			undef($lisum);

		} #$k

		my $tracker = 0;
		for (my $k=0; $k<=$#fdiff; $k++){
			if ($fdiff[$k] > 1e-3){ $tracker++; }
		} #$k

#stop iterating if all values of the distribution fct. have changed by less than 1/1e+3 in the last iteration
		if ($tracker == 0){ print STDOUT "\n"; last; }

		undef($tracker);

		$m++;

	} #$m
	undef($m);


#generate output; restrict data range to that to be plotted & generate cumulative distrib. fct. (restricted to
# the same data range)
	my @cumulf = ();
	my $finfipdl = pdl(@f);
	for (my $i=0; $i<=$#f; $i++){ $cumulf[$i] = sum($finfipdl(:$i)); };
	my $cumulfpdl = pdl(@cumulf);

	my @survabsc = list($acuabsc->where($acuhisto!=0)+$renormconst);
	my @survfct = list($cumulfpdl->where($acuhisto!=0));


#compute confidence intervals (according to Turnbull, 1974, JASA, 69, 169)
	my @nforerr = ();
	my @uforerr = ();
	my @lforerr = ();
	my @cumulforerr = ();
	my $counter = 0;
	for (my $i=0; $i<nelem($cumulfpdl); $i++){

		$uforerr[$counter] += sclr($suphisto($i));
		$lforerr[$counter] += sclr($infhisto($i));

		if ($acuhisto($i)!=0){
			$nforerr[$counter] = sclr($acuhisto($i));
			$cumulforerr[$counter] = sclr($cumulfpdl($i));
			$counter++;
		}
	
	} #$i

	undef($counter);
	$cumulforerr[$#cumulforerr] = 1;
	pop(@uforerr);
	pop(@lforerr);
	$lforerr[$#lforerr-1] += ($lforerr[$#lforerr]+$nforerr[$#nforerr]);
	pop(@nforerr);
	pop(@uforerr);
	pop(@lforerr);
	pop(@cumulforerr);

	my $DD = zeroes($#cumulforerr+1, $#cumulforerr+1);
	$DD(0,0) .= &diagDDelem(1, $cumulforerr[0], $cumulforerr[1], $nforerr[0], $nforerr[1], $lforerr[0], $uforerr[0]);

	for (my $i=1; $i<$#cumulforerr; $i++){

		$DD($i,$i) .= &diagDDelem($cumulforerr[$i-1], $cumulforerr[$i], $cumulforerr[$i+1], $nforerr[$i], $nforerr[$i+1], $lforerr[$i], $uforerr[$i]);

		for (my $j=0; $j<=$#cumulforerr; $j++){
			unless ($j==$i-1){
				next;
			} else {
				$DD($i,$j) .= &offdiagDDelem($cumulforerr[$i-1], $cumulforerr[$i], $nforerr[$i]);
				$DD($j,$i) .= $DD($i,$j);
			}
		} #$j

	} #$i

	$DD($#cumulforerr,$#cumulforerr) .= &diagDDelem($cumulforerr[$#cumulforerr-1], $cumulforerr[$#cumulforerr], 0, $nforerr[$#cumulforerr], 0, $lforerr[$#cumulforerr], $uforerr[$#cumulforerr]);
	$DD($#cumulforerr,$#cumulforerr-1) .= &offdiagDDelem($cumulforerr[$#cumulforerr-1], $cumulforerr[$#cumulforerr], $nforerr[$#cumulforerr]);
	$DD($#cumulforerr-1,$#cumulforerr) .= &offdiagDDelem($cumulforerr[$#cumulforerr-1], $cumulforerr[$#cumulforerr], $nforerr[$#cumulforerr]);


	my $VV = inv(-$DD);
	my @sigmavec = ();
	my $errstretch = 2.5; #originally had a factor 2 scaling here
	for (my $i=0; $i<=$#cumulforerr; $i++){	$sigmavec[$i] = $errstretch*sqrt(sclr($VV($i,$i))); }
	push(@sigmavec,0);

	

#adjust data range, depending on values of largest/smallest lower/upper limit(s)
	if ($nsubsubst==0){
		unshift(@survabsc,min($acuhisto)+$renormconst);
		unshift(@survfct,0);
		unshift(@sigmavec,0);
	}
	if ($ninfsubst!=0){
		pop(@survabsc);
		pop(@survfct);
		pop(@sigmavec);
	}


	return(\@survabsc, \@survfct, \@sigmavec);

	undef($origacudata);
	undef($origsupdata);
	undef($originfdata);

	undef($nsubsubst);
	undef($ninfsubst);
	undef($datarange);
	undef($rangemin);
	undef($rangemax);
	undef($binwidth);
	undef($nbins);
	undef($renormconst);

	undef($alldata);
	undef($acudata);
	undef($supdata);
	undef($infdata);

	undef($acuabsc);
	undef($supabsc);
	undef($infabsc);
	undef($acuhisto);
	undef($suphisto);
	undef($infhisto);

	undef(@n);
	undef(@u);
	undef(@l);
	undef($Ntot);

	undef($fipdl);
	undef(@f);
	undef(@fdiff);
	undef($finfipdl);

	undef(@cumulf);
	undef($cumulfpdl);

	undef(@nforerr);
	undef(@uforerr);
	undef(@lforerr);
	undef(@cumulforerr);

	undef($DD);
	undef($VV);
	
	undef(@survabsc);
	undef(@survfct);
	undef(@sigmavec);

#subroutines for the computation of the elements of the matrix D, used to compute the errors on the surv. fct.
	sub diagDDelem {

		my ($Fim1, $Fi, $Fip1, $nni, $nnip1, $lli, $uui) = @_;

		#my $outelem = -$nni/($Fi-$Fim1)**2 - $nnip1/($Fip1-$Fi)**2 - $lli/(1-$Fi)**2 - $uui/$Fi**2;
		my $outelem = 0;
		unless ($nni==0){ $outelem -= $nni/($Fi-$Fim1)**2; }
		unless ($nnip1==0){ $outelem -= $nnip1/($Fip1-$Fi)**2; }
		unless ($lli==0){ $outelem -= $lli/(1-$Fi)**2; }
		unless ($uui==0){ $outelem -= $uui/$Fi**2; }
		return ($outelem);

		undef($outelem);
		undef($Fim1);
		undef($Fi);
		undef($Fip1);
		undef($nni);
		undef($nnip1);
		undef($lli);
		undef($uui);

	} #diagDDelem

	sub offdiagDDelem {

		my ($Fi, $Fip1, $nnip1) = @_;

		my $outelem = $nnip1/($Fip1-$Fi)**2;
		return ($outelem);

		undef($outelem);
		undef($Fi);
		undef($Fip1);
		undef($nnip1);

	} #offdiagDDelem

} #doubcense_distrib


#create array of distinct measurement values, break ties btw. identical censored & uncensored values by
# means of the censorship flag (i.e. censored values are taken to be larger)
sub gendistinctdata {

	my ($indatavals, $inflagvals) = @_;
	
	my @dataarr = @{ $indatavals };
	my @flagarr = @{ $inflagvals };

	my @uniquedata = ();
	my @uniqueflags = ();
	my %seen = ();

	my $censvalflag = 0;
	for (my $i=0; $i<=$#dataarr; $i++){

		if (($i>0)&&($dataarr[$i] > $dataarr[$i-1])){ $censvalflag = 0; }

		$seen{$dataarr[$i]}++;
		if (($seen{$dataarr[$i]} > 1)&&($flagarr[$i] == 1)){

			next;

		} elsif (($seen{$dataarr[$i]} > 1)&&($flagarr[$i] == 0)){

			if ($censvalflag){

				next;

			} else {

				push @uniquedata, $dataarr[$i];
				push @uniqueflags, $flagarr[$i];
				$censvalflag++;

			}	

		} else {

			push @uniquedata, $dataarr[$i];
			push @uniqueflags, $flagarr[$i];

		}

	} #$i

	return (\@uniquedata, \@uniqueflags);

	undef(@uniquedata );
	undef(@uniqueflags);
	undef(@dataarr);
	undef(@flagarr);

	undef(%seen);

} #gendistinctdata


#wrapper subroutine for computation of Kaplan-Meier estimator; redirects to the pertient subroutines,
# depending on the kind of censoring present in the data
sub kaplan_meier {

	my ($indataarr, $inflagarr) = @_;

	my @flagarr = @{ $inflagarr };
	if ((&vecmin(@flagarr) == -1) && (&vecmax(@flagarr) == 0)){

		for (my $i=0; $i<=$#flagarr; $i++){ $flagarr[$i] += 1; }
		return (&kaplan_meier_left($indataarr, \@flagarr));

	} elsif ((&vecmin(@flagarr) == 0) && (&vecmax(@flagarr) == 1)){

		for (my $i=0; $i<=$#flagarr; $i++){ $flagarr[$i] = -1*($flagarr[$i]-1); }
		return (&kaplan_meier_right($indataarr, \@flagarr));

	} elsif ((&vecmin(@flagarr) == 0) && (&vecmax(@flagarr) == 0)){
	
		for (my $i=0; $i<=$#flagarr; $i++){ $flagarr[$i] += 1; }
		return (&kaplan_meier_right($indataarr, \@flagarr));

	} else { die("\n'kaplan_meier' is not capable of handling the requested censoring scheme!\nBailing out...\n\n"); }

	undef($indataarr);
	undef($inflagarr);
	undef(@flagarr);

} #kaplan_meier


#compute Kaplan-Meier estimator for a left-censored data set (following Feigelson & Nelson, 1985a)
sub kaplan_meier_left {

	my ($indataarr, $inflagarr) = @_;

#shift data into permitted range of values > 0
	my $tempdata = pdl(@{ $indataarr });
	my $renormconst;
	if (min($tempdata) < 0){
		$renormconst = 1.1*min($tempdata);
	} else {
		$renormconst = 0;
	}

	my $inflect = max($tempdata-$renormconst);
	my @modindataarr = list($tempdata-$renormconst);


#transform to a right-censored data set & compute according Kaplan-Meier estimator
	my @dataarr = &censotransf(\@modindataarr, $inflect);
	my @flagarr = @{ $inflagarr };

	my (@suboutpt) = &kaplan_meier_right(\@dataarr, \@flagarr);
	my $survabscissae_r = pdl(@{ $suboutpt[0] });
	my @survfct_r = @{ $suboutpt[1] };
	my @dsurvfct_r = @{ $suboutpt[2] };
	my @meaninfo_r = @{ $suboutpt[3] };


#transform right-censored statistics to left-censorship where required (& shift back to original data range)
	my @survabscissae = reverse(list(-1*$survabscissae_r + $inflect + $renormconst));
	my @survfct = reverse(@survfct_r);
	my @dsurvfct = reverse(@dsurvfct_r);
	my @meaninfo = (-1*$meaninfo_r[0] + $inflect, $meaninfo_r[1]);

	return(\@survabscissae, \@survfct, \@dsurvfct, \@meaninfo);


	undef($indataarr);
	undef($inflagarr);
	undef($tempdata);
	undef($renormconst);
	undef($inflect);
	undef(@modindataarr);
	undef(@dataarr);
	undef(@flagarr);

	undef(@suboutpt);
	undef($survabscissae_r);
	undef(@survfct_r);
	undef(@dsurvfct_r);
	undef(@meaninfo_r);

	undef(@survabscissae);
	undef(@survfct);
	undef(@dsurvfct);
	undef(@meaninfo);

} #kaplan_meier_left


#compute Kaplan-Meier estimator for a right-censored data set (following Feigelson & Nelson, 1985a)
sub kaplan_meier_right {

	my ($indataarr, $inflagarr) = @_;


#shift data into permitted range of values > 0
	my $tempdata = pdl(@{ $indataarr });
	my $renormconst;
	if (min($tempdata) < 0){
		$renormconst = 1.1*min($tempdata);
	} else {
		$renormconst = 0;
	}

	my @dataarr = list($tempdata-$renormconst);
	my @flagarr = @{ $inflagarr };


#sort input data by data value & in 2nd priority by descending censorship flag value
	my (@suboutpt) = &sortforkaplan(\@dataarr, \@flagarr);
	my @sorteddata = @{ $suboutpt[0] };
	my @sortedflags = @{ $suboutpt[1] };
	undef (@suboutpt);


#retrieve array of distinct measurement values
	my (@suboutpt) = &gendistinctdata(\@sorteddata, \@sortedflags);
	my @uniquedata = @{ $suboutpt[0] };
	my @uniqueflags = @{ $suboutpt[1] };
	undef (@suboutpt);

#compute the conditional probabilities
	my @condprob = ();
	my @d = ();
	my @n = ();

	my $counter = 0;
	foreach my $elem (@uniquedata){

		$d[$counter] = 0;
		$n[$counter] = 0;
		for (my $i=0; $i<=$#sorteddata; $i++){
			if (($sorteddata[$i] == $elem)&&($sortedflags[$i] == 1)){ $d[$counter]++; }
			if (($sorteddata[$i] >= $elem)||(($sorteddata[$i] == $elem)&&($sortedflags[$i] == 0))){ $n[$counter]++; }
		} #$i

		$condprob[$counter] = (1-$d[$counter]/$n[$counter])**$uniqueflags[$counter];

		$counter++;

	} #$elem
	$counter = 0;


#calculate the Kaplan-Meier estimator (the survivor fct.) at the jumps (i.e. at the uncensored data points + epsilon)
	my @survabsc = @uniquedata;
	my @survfct = ();
	my @dsurvfct = ();

	for (my $i=0; $i<=$#condprob; $i++){

		my ($numeratortemp, $denomtemp);
		#beneath the value of the smallest data point, the survival fct. must be unity & the error zero
		if ($i == 0){
			$numeratortemp = 0;
			$denomtemp = 1;
			$survfct[$i] = 1;
			goto sigmacalc;
		} else {
			$numeratortemp = 0;
			$denomtemp = 0;
		}

		$survfct[$i] = 1;

		for (my $j=0; $j<$i; $j++){

			$survfct[$i] *= $condprob[$j];
			$numeratortemp += $d[$j]*$uniqueflags[$j]/$n[$j]/($n[$j]-$d[$j]);
			$denomtemp += $uniqueflags[$j]*log(($n[$j]-$d[$j])/$n[$j]);

			if (($survfct[$i] == 1)&&($numeratortemp == 0)){ $denomtemp = 1; }

		}

		sigmacalc:
		$dsurvfct[$i] = sqrt($numeratortemp/$denomtemp**2);

	} #$i

	if ($sortedflags[$#sortedflags] == 1){
		push(@survfct, $survfct[$#survfct]*$condprob[$#condprob]);
		push(@dsurvfct, 0);
		push(@survabsc, $uniquedata[$#uniquedata]+($uniquedata[$#uniquedata]-$uniquedata[0])/1e+6)
	}


#calculate the mean and its error (assuming that the largest measurement is uncensored if necessary)
	my $addzeroflg = 0;
	if ($survfct[$#survfct] != 0){
		$addzeroflg++;
		push(@survfct, 0);
		push(@dsurvfct, 0);
		push(@survabsc, $uniquedata[$#uniquedata]+($uniquedata[$#uniquedata]-$uniquedata[0])/1e+6);
	}

	#undo shift in data range
	for (my $i=0; $i<=$#survabsc; $i++){ $survabsc[$i] += $renormconst; }

	my $mean = $survfct[0]*($survabsc[0]-0);
	for (my $i=1; $i<=$#survabsc; $i++){ $mean += $survfct[$i]*($survabsc[$i]-$survabsc[$i-1]); }

	my $dmeansq = 0;
	for (my $i=0; $i<=$#uniquedata; $i++){
		my $sumarg = 0;
		for (my $j=$i; $j<=$#uniquedata; $j++){
				$sumarg += $survfct[$j+1]*($survabsc[$j+1]-$survabsc[$j]);
		} #$j
		my $denom = $n[$i]*($n[$i]-$d[$i]);
		my $numer = $uniqueflags[$i]*$d[$i];
		if (($sumarg == 0)&&($denom == 0)){ next; }
		$dmeansq += ($sumarg**2)*$numer/$denom;
	} #$i
	$counter = 0;

	my @meaninfo = ($mean, sqrt($dmeansq));


#if necessary undo assumption that the largest measurement is not censored
	if ($addzeroflg == 1){
		pop(@survabsc);
		pop(@survfct);
		pop(@dsurvfct);
	}
	undef($addzeroflg);

	return(\@survabsc, \@survfct, \@dsurvfct, \@meaninfo);


	undef($indataarr);
	undef($inflagarr);
	undef($tempdata);
	undef($renormconst);
	undef(@dataarr);
	undef(@flagarr);

	undef(@sorteddata );
	undef(@sortedflags);
	undef(@uniquedata );
	undef(@uniqueflags);
	undef(@condprob);

	undef(@survabsc);
	undef(@survfct);
	undef(@dsurvfct);
	undef(@meaninfo)

} #kaplan_meier_right


#subroutine which converts the output of doubcense_distrib to an easily printable format
sub plotdistrprep_double {

	my ($inabsc, $inyvals, $inerrs, $plotrange) = @_;

	my @survabsc = @{ $inabsc };
	my @survfct = @{ $inyvals };
	my @dsurvfct = @{ $inerrs };
	my @plotlims = @{ $plotrange };

	my @plotabsc = ();
	my @plotfctvals = ();
	my @ploterrabsc = ();
	my @ploterrfctvals = ();
	my @plotuerrs = ();
	my @plotlerrs = ();


	my $counter = 0;
	for (my $i=0; $i<=$#survabsc; $i++){

		if ($i==0){

			$plotabsc[$counter] = $survabsc[$i];
			$plotfctvals[$counter] = $survfct[$i];

			$counter++;

		} else {

			$plotabsc[$counter] = $survabsc[$i]-($plotlims[1]-$plotlims[0])/1e+9;
			$plotabsc[$counter+1] = $survabsc[$i]+($plotlims[1]-$plotlims[0])/1e+9;

			$plotfctvals[$counter] = $survfct[$i-1];
			$plotfctvals[$counter+1] = $survfct[$i];

			$counter += 2;

		}

	} #$i
	undef($counter);

	push(@plotabsc, $plotlims[1]);
	unshift(@plotabsc, $plotlims[0]);
	push(@plotfctvals, @plotfctvals[$#plotfctvals]);
	unshift(@plotfctvals, @plotfctvals[0]);


	for (my $i=0; $i<=$#survabsc; $i++){

		$ploterrabsc[$i] = $survabsc[$i];

		if ($i==0){
			$plotuerrs[$i] = $survfct[$i]**(exp(-1.96*$dsurvfct[$i]));
			$plotlerrs[$i] = $survfct[0]**(exp(1.96*$dsurvfct[0]));
		} elsif ($i==$#survabsc){
			$plotuerrs[$i] = $survfct[$#survfct]**(exp(-1.96*$dsurvfct[$#dsurvfct]));
			$plotlerrs[$i] = $survfct[$i-1]**(exp(1.96*$dsurvfct[$i-1]));
		} else {
			$plotuerrs[$i] = $survfct[$i]**(exp(-1.96*$dsurvfct[$i]));
			$plotlerrs[$i] = $survfct[$i-1]**(exp(1.96*$dsurvfct[$i-1]));
		}
	
	} #$i

	push(@ploterrabsc, $survabsc[$#survabsc]+(&vecmax(@survfct)-&vecmin(@survfct))/10, $plotlims[1]);
	push(@plotuerrs, $survfct[$#survfct]**(exp(-1.96*$dsurvfct[$#dsurvfct])), $survfct[$#survfct]**(exp(-1.96*$dsurvfct[$#dsurvfct])));
	push(@plotlerrs, $survfct[$#survfct]**(exp(1.96*$dsurvfct[$#dsurvfct])), $survfct[$#survfct]**(exp(1.96*$dsurvfct[$#dsurvfct])));
	unshift(@ploterrabsc, $plotlims[0], $survabsc[0]-(&vecmax(@survfct)-&vecmin(@survfct))/10);
	unshift(@plotuerrs, $plotuerrs[0], $plotuerrs[0]);
	unshift(@plotlerrs, $plotlerrs[0], $plotlerrs[0]);


	return (\@plotabsc, \@plotfctvals, \@ploterrabsc, \@plotuerrs, \@plotlerrs);

	undef($inabsc);
	undef($inyvals);
	undef($plotrange);
	
	undef(@survabsc);
	undef(@survfct);
	undef(@plotlims);

	undef(@plotabsc);
	undef(@plotfctvals);

} #plotdistrprep_double


#subroutine which converts the output of kaplan_meier_left to an easily printable format
sub plotdistrprep_left {

	my ($inabsc, $inyvals, $inerrs, $plotrange) = @_;

	my @survabsc = @{ $inabsc };
	my @survfct = @{ $inyvals };
	my @dsurvfct = @{ $inerrs };
	my @plotlims = @{ $plotrange };

	my @plotabsc = ();
	my @plotfctvals = ();
	my @ploterrabsc = ();
	my @plotuerrs = ();
	my @plotlerrs = ();


	my $counter = 0;
	for (my $i=0; $i<=$#survabsc; $i++){

		if ($i == 0){
			$plotfctvals[$counter] = $survfct[$i];
		} else {
			$plotfctvals[$counter] = $survfct[$i-1];
		}
		$plotfctvals[$counter+1] = $survfct[$i];

		$plotabsc[$counter] = $survabsc[$i]-($plotlims[1]-$plotlims[0])/1e+9;
		$plotabsc[$counter+1] = $survabsc[$i]+($plotlims[1]-$plotlims[0])/1e+9;

		$counter += 2;
	} #$i
	$counter = 0;


	my $i=$#survabsc;
	my $prevjumpabsc;
	while ($i > 0){

		if ($i == $#survabsc){

			$ploterrabsc[$counter] = $survabsc[$i] + ($plotlims[1]-$plotlims[0])/1e+9;
			$ploterrabsc[$counter+1] = $survabsc[$i];

			$plotuerrs[$counter] = 1;
			$plotuerrs[$counter+1] = 1;

			$plotlerrs[$counter] = 1;
			$plotlerrs[$counter+1] = $survfct[$i-1]**(exp(1.96*$dsurvfct[$i-1]));

			$counter += 2;
			$prevjumpabsc = $survabsc[$i];

		} elsif ($i == 1){

				$ploterrabsc[$counter] = $prevjumpabsc - ($prevjumpabsc-$survabsc[$i])/2;
				$ploterrabsc[$counter+1] = $survabsc[$i];
				$ploterrabsc[$counter+2] = $survabsc[$i] - ($plotlims[1]-$plotlims[0])/1e+9;

				$plotuerrs[$counter] = $survfct[$i+1]**(exp(-1.96*$dsurvfct[$i+1]));
				$plotlerrs[$counter] = $survfct[$i+1]**(exp(1.96*$dsurvfct[$i+1]));
			
				if ($survfct[0] == 0){
					$plotuerrs[$counter+1] = $plotuerrs[$counter];
					$plotuerrs[$counter+2] = 0;

					$plotlerrs[$counter+1] = 0;
					$plotlerrs[$counter+2] = 0;
				} else {
					$plotuerrs[$counter+1] = $plotuerrs[$counter];#$survfct[$i]**(exp(-1.96*$dsurvfct[$i]));
					$plotlerrs[$counter+1] = $survfct[$i]**(exp(1.96*$dsurvfct[$i]));

					$plotuerrs[$counter+2] = $plotuerrs[$counter+1];
					$plotlerrs[$counter+2] = $plotlerrs[$counter+1];
				}

		} else {

			if ($survfct[$i] != $survfct[$i+1]){

				$ploterrabsc[$counter] = $prevjumpabsc - ($prevjumpabsc-$survabsc[$i])/2;
				$plotuerrs[$counter] = $survfct[$i]**(exp(-1.96*$dsurvfct[$i]));
				$plotlerrs[$counter] = $survfct[$i]**(exp(1.96*$dsurvfct[$i]));

				$counter++;
				$prevjumpabsc = $survabsc[$i];

			}

		}

		$i--;

	} #$i
	undef($i);
	undef($counter);
	undef($prevjumpabsc);

	unless ($survfct[0] == 0){
		push(@ploterrabsc, $plotlims[0]);
		push(@plotuerrs, $plotuerrs[$#plotuerrs]);
		push(@plotlerrs, $plotlerrs[$#plotlerrs]);
	}

	@ploterrabsc = reverse(@ploterrabsc);
	@plotuerrs = reverse(@plotuerrs);
	@plotlerrs = reverse(@plotlerrs);


	push(@plotabsc, $plotlims[1]);
	push(@plotfctvals, @plotfctvals[$#plotfctvals]);
	unshift(@plotabsc, $plotlims[0]);
	unshift(@plotfctvals, @plotfctvals[0]);


	return (\@plotabsc, \@plotfctvals, \@ploterrabsc, \@plotuerrs, \@plotlerrs);

	undef(@survabsc);
	undef(@survfct);
	undef(@dsurvfct);

	undef(@plotabsc);
	undef(@plotfctvals);
	undef(@ploterrabsc);
	undef(@plotuerrs);
	undef(@plotlerrs);

} #plotdistrprep_left


#subroutine which converts the output of kaplan_meier_right to an easily printable format
sub plotdistrprep_right {

	my ($inabsc, $inyvals, $inerrs, $plotrange) = @_;

	my @survabsc = @{ $inabsc };
	my @survfct = @{ $inyvals };
	my @dsurvfct = @{ $inerrs };
	my @plotlims = @{ $plotrange };

	my @plotabsc = ();
	my @plotfctvals = ();
	my @ploterrabsc = ();
	my @plotuerrs = ();
	my @plotlerrs = ();


	my $counter = 0;
	for (my $i=0; $i<=$#survabsc; $i++){

		if ($i == 0){
			$plotfctvals[$counter] = $survfct[$i];
		} else {
			$plotfctvals[$counter] = $survfct[$i-1];
		}
		$plotfctvals[$counter+1] = $survfct[$i];

		$plotabsc[$counter] = $survabsc[$i]-($plotlims[1]-$plotlims[0])/1e+9;
		$plotabsc[$counter+1] = $survabsc[$i]+($plotlims[1]-$plotlims[0])/1e+9;

		$counter += 2;
	} #$i
	$counter = 0;


	my $i=1;
	my $prevjumpabsc;
	while ($i <= $#survabsc){

		if ($i == 1){

			$ploterrabsc[$counter] = $survabsc[$i] - ($plotlims[1]-$plotlims[0])/1e+9;
			$ploterrabsc[$counter+1] = $survabsc[$i];

			$plotuerrs[$counter] = 1;
			$plotuerrs[$counter+1] = 1;

			$plotlerrs[$counter] = 1;
			$plotlerrs[$counter+1] = $survfct[$i]**(exp(1.96*$dsurvfct[$i]));

			$counter += 2;
			$prevjumpabsc = $survabsc[$i];

		} elsif ($i == $#survabsc){

				$ploterrabsc[$counter] = $prevjumpabsc + ($survabsc[$i]-$prevjumpabsc)/2;
				$ploterrabsc[$counter+1] = $survabsc[$i];
				$ploterrabsc[$counter+2] = $survabsc[$i] + ($plotlims[1]-$plotlims[0])/1e+9;

				$plotuerrs[$counter] = $survfct[$i-1]**(exp(-1.96*$dsurvfct[$i-1]));
				$plotlerrs[$counter] = $survfct[$i-1]**(exp(1.96*$dsurvfct[$i-1]));
			
				if ($survfct[$#survfct] == 0){
					$plotuerrs[$counter+1] = $plotuerrs[$counter];
					$plotuerrs[$counter+2] = 0;

					$plotlerrs[$counter+1] = 0;
					$plotlerrs[$counter+2] = 0;
				} else {
					$plotuerrs[$counter+1] = $plotuerrs[$counter];$survfct[$i]**(exp(-1.96*$dsurvfct[$i]));
					$plotlerrs[$counter+1] = $survfct[$i]**(exp(1.96*$dsurvfct[$i]));

					$plotuerrs[$counter+2] = $plotuerrs[$counter+1];
					$plotlerrs[$counter+2] = $plotlerrs[$counter+1];
				}

		} else {

			if ($survfct[$i] != $survfct[$i-1]){

				$ploterrabsc[$counter] = $prevjumpabsc + ($survabsc[$i]-$prevjumpabsc)/2;
				$plotuerrs[$counter] = $survfct[$i-1]**(exp(-1.96*$dsurvfct[$i-1]));
				$plotlerrs[$counter] = $survfct[$i-1]**(exp(1.96*$dsurvfct[$i-1]));

				$counter++;
				$prevjumpabsc = $survabsc[$i];

			}

		}

		$i++;

	} #$i
	undef($i);
	undef($counter);
	undef($prevjumpabsc);


	unless ($survfct[$#survfct] == 0){
		push(@ploterrabsc, $plotlims[1]);
		push(@plotuerrs, $plotuerrs[$#plotuerrs]);
		push(@plotlerrs, $plotlerrs[$#plotlerrs]);
	}


	push(@plotabsc, $plotlims[1]);
	unshift(@plotabsc, $plotlims[0]);
	push(@plotfctvals, $plotfctvals[$#plotfctvals]);
	unshift(@plotfctvals, $plotfctvals[0]);


	return (\@plotabsc, \@plotfctvals, \@ploterrabsc, \@plotuerrs, \@plotlerrs);

	undef(@survabsc);
	undef(@survfct);
	undef(@dsurvfct);

	undef(@plotabsc);
	undef(@plotfctvals);
	undef(@ploterrabsc);
	undef(@plotuerrs);
	undef(@plotlerrs);

} #plotdistrprep_right


#routine called by kaplan_meier_right to sort input data & censorship flags
sub sortforkaplan {

	my ($vals, $flags) = @_;
	
	my @valsarr = @{ $vals };
	my @flagsarr = @{ $flags };

	my @unsortedinfo = ();
	for (my $i=0; $i<=$#valsarr; $i++){ $unsortedinfo[$i] = $valsarr[$i]."\t".$flagsarr[$i]; }
	my @sortedinfo = fieldsort ['1n','-2n'], @unsortedinfo;

	my @ordvals = ();
	my @ordflags = ();
	for (my $i=0; $i<=$#sortedinfo; $i++){ ($ordvals[$i], $ordflags[$i]) = split(/\s+/, $sortedinfo[$i]); }

	return  (\@ordvals, \@ordflags);

	undef(@valsarr);
	undef(@flagsarr);
	undef(@unsortedinfo);
	undef(@sortedinfo);
	undef(@ordvals);
	undef(@ordflags);

} #sortforkaplan


#compute median of survival function & associated error
sub survfct_median {

	my ($inxx, $inyy, $inxxd, $indsup, $indinf, $type) = @_;

	my @xx = @{ $inxx };
	my @yy = @{ $inyy };
	my @xxd = @{ $inxxd };
	my @dsup = @{ $indsup };
	my @dinf = @{ $indinf };


	my ($infxx, $supxx, $infyy, $supyy);
	my ($median, $dinfmed, $dsupmed);

	if (($type eq "double")||($type eq "left")){

		for (my $i=0; $i<=$#xx; $i++){
			if ($yy[$i]>.5){
				$supxx = $xx[$i];
				$supyy = $yy[$i];
				last;
			}
			$infxx = $xx[$i];
			$infyy = $yy[$i];
		}
		$median = &lininterp($infyy, $supyy, $infxx, $supxx, .5);

		for (my $i=0; $i<=$#xxd; $i++){
			if ($dsup[$i]>.5){
				$supxx = $xxd[$i];
				$supyy = $dsup[$i];
				last;
			}
			$infxx = $xxd[$i];
			$infyy = $dsup[$i];
		}
		$dsupmed = $median-&lininterp($infyy, $supyy, $infxx, $supxx, .5);
	
		for (my $i=0; $i<=$#xxd; $i++){
			if ($dinf[$i]>.5){
				$supxx = $xxd[$i];
				$supyy = $dinf[$i];
				last;
			}
			$infxx = $xxd[$i];
			$infyy = $dinf[$i];
		}
		$dinfmed = &lininterp($infyy, $supyy, $infxx, $supxx, .5)-$median;

	} else {

		for (my $i=0; $i<=$#xx; $i++){
			if ($yy[$i]<.5){
				$infxx = $xx[$i];
				$infyy = $yy[$i];
				last;
			}
			$supxx = $xx[$i];
			$supyy = $yy[$i];
		}
		$median = &lininterp($infyy, $supyy, $infxx, $supxx, .5);

		for (my $i=0; $i<=$#xxd; $i++){
			if ($dsup[$i]<.5){
				$infxx = $xxd[$i];
				$infyy = $dsup[$i];
				last;
			}
			$supxx = $xxd[$i];
			$supyy = $dsup[$i];
		}
		$dsupmed = &lininterp($infyy, $supyy, $infxx, $supxx, .5)-$median;
	
		for (my $i=0; $i<=$#xxd; $i++){
			if ($dinf[$i]<.5){
				$infxx = $xxd[$i];
				$infyy = $dinf[$i];
				last;
			}
			$supxx = $xxd[$i];
			$supyy = $dinf[$i];
		}
		$dinfmed = $median-&lininterp($infyy, $supyy, $infxx, $supxx, .5);
	
	}

	undef($infxx);
	undef($supxx);
	undef($infyy);
	undef($supyy);


	return($median, $dinfmed, $dsupmed);

	undef($median);
	undef($dinfmed);
	undef($dsupmed);

} #survfct_median


#compute median of survival function & associated error
sub survfct_pct {

	my ($inxx, $inyy, $inxxd, $indsup, $indinf, $type, $perc) = @_;

	my @xx = @{ $inxx };
	my @yy = @{ $inyy };
	my @xxd = @{ $inxxd };
	my @dsup = @{ $indsup };
	my @dinf = @{ $indinf };


	my ($infxx, $supxx, $infyy, $supyy);
	my ($pct, $dinfpct, $dsuppct);

	if (($type eq "double")||($type eq "left")){

		for (my $i=0; $i<=$#xx; $i++){
			if ($yy[$i]>$perc){
				$supxx = $xx[$i];
				$supyy = $yy[$i];
				last;
			}
			$infxx = $xx[$i];
			$infyy = $yy[$i];
		}
		$pct = &lininterp($infyy, $supyy, $infxx, $supxx, $perc);

		for (my $i=0; $i<=$#xxd; $i++){
			if ($dsup[$i]>$perc){
				$supxx = $xxd[$i];
				$supyy = $dsup[$i];
				last;
			}
			$infxx = $xxd[$i];
			$infyy = $dsup[$i];
		}
		$dsuppct = $pct-&lininterp($infyy, $supyy, $infxx, $supxx, $perc);
	
		for (my $i=0; $i<=$#xxd; $i++){
			if ($dinf[$i]>$perc){
				$supxx = $xxd[$i];
				$supyy = $dinf[$i];
				last;
			}
			$infxx = $xxd[$i];
			$infyy = $dinf[$i];
		}
		$dinfpct = &lininterp($infyy, $supyy, $infxx, $supxx, $perc)-$pct;

	} else {

		for (my $i=0; $i<=$#xx; $i++){
			if ($yy[$i]<$perc){
				$infxx = $xx[$i];
				$infyy = $yy[$i];
				last;
			}
			$supxx = $xx[$i];
			$supyy = $yy[$i];
		}
		$pct = &lininterp($infyy, $supyy, $infxx, $supxx, $perc);

		for (my $i=0; $i<=$#xxd; $i++){
			if ($dsup[$i]<$perc){
				$infxx = $xxd[$i];
				$infyy = $dsup[$i];
				last;
			}
			$supxx = $xxd[$i];
			$supyy = $dsup[$i];
		}
		$dsuppct = &lininterp($infyy, $supyy, $infxx, $supxx, $perc)-$pct;
	
		for (my $i=0; $i<=$#xxd; $i++){
			if ($dinf[$i]<$perc){
				$infxx = $xxd[$i];
				$infyy = $dinf[$i];
				last;
			}
			$supxx = $xxd[$i];
			$supyy = $dinf[$i];
		}
		$dinfpct = $pct-&lininterp($infyy, $supyy, $infxx, $supxx, $perc);
	
	}

	undef($infxx);
	undef($supxx);
	undef($infyy);
	undef($supyy);


	return($pct, $dinfpct, $dsuppct);

	undef($pct);
	undef($dinfpct);
	undef($dsuppct);

} #survfct_mpct


1;              # Modules must return true.
