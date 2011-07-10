#!/usr/bin/perl

require "digamma.pl";

for(my $x = 0.01; $x < 2; $x += 0.01) {
    print "$x,".digamma($x).",".exp(digamma($x))."\n";
}
