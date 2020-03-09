package Own::Misc::FitLM_Funcs;

use strict;
use Exporter;

use PDL;

our $VERSION = 1.0;

our @ISA = qw(Exporter);
our @EXPORT = qw(&erfunc_fixmean &gaussian &gauss_asym &gauss_asymnorm &linefit &lognormschechter &onepluszfit &schechter &schechter_fixslope &schechter_varymstar);

use constant PI => (4*atan2(1,1));


#Error function (Erf(x); http://mathworld.wolfram.com/Erf.html)
sub erfunc_fixmean {

	my ($x,$par,$ym,$dyda) = @_;

	my $avq = $GLOBAL::input;

	my ($scat) = map { $par->slice("($_)") } (0);

	$ym .= .5*(1+erf(($x-$avq)/$scat/sqrt(2)));
	my (@dy) = map {$dyda -> slice(",($_)") } (0);
	
	$dy[0] .= -1/sqrt(2*PI)*exp(-1*(($x-$avq)/$scat/sqrt(2))**2)*($x-$avq)/$scat**2;

} #erfunc


#Gaussian
sub gaussian {

	my ($x,$par,$ym,$dyda) = @_;
	my ($amp, $mean, $sig) = map { $par->slice("($_)") } (0..2);
	
	$ym .= $amp*exp(-1*(($x-$mean)**2)/(2*$sig**2));
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0..2);
	
	#derivative w/ respect to normalization
	$dy[0] .= $ym/$amp;
	#derivative w/ respect to the mode of the Gaussian
	$dy[1] .= $ym*($x-$mean)/$sig**2;
	#derivative w/ respect to sigma
	$dy[2] .= $ym/$sig*((($x-$mean)/$sig)**2 - 1);
	
}


#asymmetric Gaussian
sub gauss_asym {

	my ($x,$par,$ym,$dyda) = @_;
	my ($amp, $mean, $sig, $rr) = map { $par->slice("($_)") } (0..3);
	
	my $frag1 = $amp/($sig*($rr+1))*exp(-1*(($x->where($x<=$mean)-$mean)**2)/(2*($rr*$sig)**2));
	my $frag2 = $amp/($sig*($rr+1))*exp(-1*(($x->where($x>$mean)-$mean)**2)/(2*$sig**2));
	$ym .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0..3);
	
	#derivative w/ respect to normalization
	$dy[0] .= $ym/$amp;
	#derivative w/ respect to the mode of the Gaussian
	my $frag1 = $ym->where($x<=$mean)*($x->where($x<=$mean)-$mean)/($rr*$sig)**2;
	my $frag2 = $ym->where($x>$mean)*($x->where($x>$mean)-$mean)/$sig**2;
	$dy[1] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	#derivative w/ respect to sigma
	my $frag1 = $ym->where($x<=$mean)/$sig*((($x->where($x<=$mean)-$mean)/($rr*$sig))**2 - 1);
	my $frag2 = $ym->where($x>$mean)/$sig*((($x->where($x>$mean)-$mean)/$sig)**2 - 1);
	$dy[2] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	#derivative w/ respect to the ratio of the upper an lower width
	my $frag1 = $ym->where($x<=$mean)*(((($x->where($x<=$mean)-$mean)/($rr*$sig))**2)/$rr - 1/($rr+1));
	my $frag2 = $ym->where($x>$mean)*((($x->where($x>$mean)-$mean)/$sig)**2 - .5);
	$dy[3] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	
}


#normalized asymmetric Gaussian
sub gauss_asymnorm {

	my ($x,$par,$ym,$dyda) = @_;
	my ($mean, $sig, $rr) = map { $par->slice("($_)") } (0..2);
	
	my $frag1 = sqrt(2/PI)/($sig*(sqrt($rr)**2+1))*exp(-1*(($x->where($x<=$mean)-$mean)**2)/(2*(sqrt($rr)**2*$sig)**2));
	my $frag2 = sqrt(2/PI)/($sig*(sqrt($rr)**2+1))*exp(-1*(($x->where($x>$mean)-$mean)**2)/(2*$sig**2));
	$ym .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0..2);
	
#	#derivative w/ respect to normalization
#	$dy[0] .= $ym/$amp;
	#derivative w/ respect to the mode of the Gaussian
	my $frag1 = $ym->where($x<=$mean)*($x->where($x<=$mean)-$mean)/(sqrt($rr)**2*$sig)**2;
	my $frag2 = $ym->where($x>$mean)*($x->where($x>$mean)-$mean)/$sig**2;
	$dy[0] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	#derivative w/ respect to sigma
	my $frag1 = $ym->where($x<=$mean)/$sig*((($x->where($x<=$mean)-$mean)/(sqrt($rr)**2*$sig))**2 - 1);
	my $frag2 = $ym->where($x>$mean)/$sig*((($x->where($x>$mean)-$mean)/$sig)**2 - 1);
	$dy[1] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	#derivative w/ respect to the ratio of the upper an lower width
	my $frag1 = $ym->where($x<=$mean)*(((($x->where($x<=$mean)-$mean)/(sqrt($rr)**2*$sig))**2)/sqrt($rr)**2 - 1/(sqrt($rr)**2+1));
	my $frag2 = $ym->where($x>$mean)*((($x->where($x>$mean)-$mean)/$sig)**2 - .5);
	$dy[2] .= $frag1->append($frag2);
	undef($frag1);
	undef($frag2);
	
}


#linear fct.
sub linefit {

	my ($x,$par,$ym,$dyda) = @_;

	my ($m,$b) = map { $par->slice("($_)") } (0..1);

	$ym .= $m*$x + $b;
	my (@dy) = map {$dyda -> slice(",($_)") } (0..1);
	
	$dy[0] .= $x;
	$dy[1] .= 1;

} #linefit


#f(z) = f(z=0)*(1+z)^a
sub onepluszfit {

	my ($x,$par,$ym,$dyda) = @_;

	my ($a,$b) = map { $par->slice("($_)") } (0..1);

	$ym .= $a*(1+$x)**$b;
	my (@dy) = map {$dyda -> slice(",($_)") } (0..1);
	
	$dy[0] .= $ym/$a;
	$dy[1] .= $ym*log(1+$x);

} #onepluszfit


#Schechter fct. convolved w/ lognormal distribution (as in bivariate size-luminosity fct.)
sub lognormschechter {

	my ($x1,$x2,$par,$ym,$dyda) = @_;
	#my ($rstar,$sigr,$beta) = map { $par->slice("($_)") } (0..2);
	#my $phistar = 0.002111969;
	#my $alpha = -1.1267738;
	#my $mstar = -20.759275;
	my ($rstar,$sigr,$beta,$phistar,$alpha,$mstar) = map { $par->slice("($_)") } (0..5);
	
	$ym .= .4*log(10)*$phistar/(sqrt(2*PI)*$sigr*$x1)*(10**((-0.4)*($alpha+1)*($x2-$mstar)))*exp((-1)*(10**((-0.4)*($x2-$mstar))))*exp((-1)*(log(10)**2)/(2*$sigr**2)*(log($x1/$rstar)/log(10)-0.4*$beta*($x2-$mstar))**2);
	
	#my (@dy) = map {$dyda -> slice(",($_)") } (0..2);
	my (@dy) = map {$dyda -> slice(",($_)") } (0..5);
	
	$dy[0] .= $ym*log(10)/($rstar*$sigr**2)*(log($x1/$rstar)/log(10)-0.4*$beta*($x2-$mstar));
	$dy[1] .= $ym/$sigr*(((log(10)/$sigr)**2)*(log($x1/$rstar)/log(10)-0.4*$beta*($x2-$mstar))**2 - 1);
	$dy[2] .= .4*$ym*((log(10)/$sigr)**2)*(log($x1/$rstar)/log(10)-0.4*$beta*($x2-$mstar))*($x2-$mstar);
	$dy[3] .= $ym/$phistar;
	$dy[4] .= (-.4)*log(10)*$ym*($x2-$mstar);
	$dy[5] .= .4*log(10)*$ym*(($alpha+1)-10**((-0.4)*($x2-$mstar))-$beta*log(10)/($sigr**2)*(log($x1/$rstar)/log(10)-0.4*$beta*($x2-$mstar)));

} #lognormschechter


#Schechter function
sub schechter {

	my ($x,$par,$ym,$dyda) = @_;
	my ($phistar,$mstar,$alpha) = map { $par->slice("($_)") } (0..2);
	
	#$ym .= log(.4*log(10)*$phistar*(10**((-0.4)*($x-$mstar)*($alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar)))))/log(10);
	$ym .= .4*log(10)*$phistar*(10**((-0.4)*($x-$mstar)*($alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar))));
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0..2);
	
	#$dy[0] .= 1/($phistar*log(10));
	#$dy[1] .= (-0.4)*(10**((-0.4)*($x-$mstar)) - ($alpha+1));
	#$dy[2] .= (-0.4)*($x-$mstar);
	
	$dy[0] .= $ym/$phistar;
	$dy[1] .= 0.4*log(10)*$ym*(($alpha+1) - 10**((-0.4)*($x-$mstar)));
	$dy[2] .= (-0.4)*log(10)*$ym*($x-$mstar);
	
} #schechter


#Schechter function w/ fixed faint end slope
sub schechter_fixslope {

	my ($x,$par,$ym,$dyda) = @_;
	my ($phistar,$mstar) = map { $par->slice("($_)") } (0..1);

	#$ym .= log(.4*log(10)*$phistar*(10**((-0.4)*($x-$mstar)*($alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar)))))/log(10);
	$ym .= .4*log(10)*$phistar*(10**((-0.4)*($x-$mstar)*($main::alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar))));
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0..1);
	
	#$dy[0] .= 1/($phistar*log(10));
	#$dy[1] .= (-0.4)*(10**((-0.4)*($x-$mstar)) - ($alpha+1));
	
	$dy[0] .= $ym/$phistar;
	$dy[1] .= 0.4*log(10)*$ym*(($main::alpha+1) - 10**((-0.4)*($x-$mstar)));
	
} #schechter_fixslope


#Schechter function w/ fixed faint end slope & normalization
sub schechter_varymstar {

	my ($x,$par,$ym,$dyda) = @_;
	my ($mstar) = map { $par->slice("($_)") } (0);

	#$ym .= log(.4*log(10)*$phistar*(10**((-0.4)*($x-$mstar)*($alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar)))))/log(10);
	$ym .= .4*log(10)*$main::phistar*(10**((-0.4)*($x-$mstar)*($main::alpha+1)))*exp((-1)*(10**((-0.4)*($x-$mstar))));
	
	my (@dy) = map {$dyda -> slice(",($_)") } (0);
	
	#$dy[0] .= 1/($phistar*log(10));
	#$dy[1] .= (-0.4)*(10**((-0.4)*($x-$mstar)) - ($alpha+1));
	
	$dy[0] .= 0.4*log(10)*$ym*(($main::alpha+1) - 10**((-0.4)*($x-$mstar)));
	
} #schechter_varymstar


1;              # Modules must return true.
