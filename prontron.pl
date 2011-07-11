#!/usr/bin/perl

# TODO, fix random bias
# TODO, shuffle

use strict;
use utf8;
use List::Util qw(sum min max shuffle);
use Getopt::Long;
use FindBin qw($Bin);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
binmode STDERR, ":utf8";

require "$Bin/prontron-algorithm.pl";

my ($FMIN, $FMAX, $EMIN, $EMAX) = (1,1,0,5);
my $WORD = 0;
GetOptions(
    "fmin=i" => \$FMIN,
    "fmax=i" => \$FMAX,
    "emin=i" => \$EMIN,
    "emax=i" => \$EMAX,
    "word" => \$WORD   # use word units instead of characters
);
my %config = ( "fmin" => $FMIN, "fmax" => $FMAX, "emin" => $EMIN, "emax" => $EMAX, "word" => $WORD );

if(@ARGV != 2) {
    print STDERR "Usage: disc-decode.pl PFILE WFILE < INPUT > OUTPUT\n";
    exit 1;
}

open FILE, "<:utf8", $ARGV[0] or die "$ARGV[0]: $!\n";
my (%cands);
while(<FILE>) {
    chomp; my ($f,$e,$p) = split(/\t/);
    $cands{$f} = [] if not $cands{$f}; push @{$cands{$f}}, [$e,$p];
}
close FILE;
open FILE, "<:utf8", $ARGV[1] or die "$ARGV[1]: $!\n";
my %weights = map { chomp; my ($f,$w) = split(/\t/); $f => $w } <FILE>;
close FILE;

# perform ITERS iterations
while(<STDIN>) {
    chomp;
    my ($feats, $seq, $score, $str) = decode($_,\%weights,\%cands,\%config);
    print "$str\n";
}
