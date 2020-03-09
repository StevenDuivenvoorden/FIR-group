package Own::Misc::Utilities;



use strict;

use Exporter;



our $VERSION = 1.0;



our @ISA = qw(Exporter);

our @EXPORT = qw(&buildhash &trim &round);





#subroutine that constructs & returns a hash from two input arrays (one with the hash's keys,

# the other one with the hash's values); to be called as: %newhash = &buildhash( \@keys, \@vals )

sub buildhash{



	my ($hashkeys, $hashvals) = @_;

	my %soughthash = ();



	for (my $j=0; $j < @$hashkeys; $j++){

		my $key = $hashkeys->[$j];

		$soughthash{"$key"} = $hashvals->[$j];

	}

	

	return %soughthash;



}





#subroutine that rounds numbers

sub round {



	my ($number) = shift;

	return int($number + .5 * ($number <=>0));



}





#Perl trim function to remove whitespace from the start and end of an array

sub trim {



	my (@vector) = @_;

	my $j = 0;

	

	while ($j <= $#vector){

		$vector[$j] =~ s/^\s+//;

		$vector[$j] =~ s/\s+$//;

		$j += 1;

	}

	

	return @vector;



} #trim



1;              # Modules must return true.
