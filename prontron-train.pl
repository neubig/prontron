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
my $ITERS = 50;
my $WORD = 0;
my $INAROW = 5;
my $RECHECK = 10;
GetOptions(
    "fmin=i" => \$FMIN,
    "fmax=i" => \$FMAX,
    "emin=i" => \$EMIN,
    "emax=i" => \$EMAX,
    "iters=i" => \$ITERS,
    "word" => \$WORD,      # use word units instead of characters
    "inarow" => \$INAROW,  # skip training examples we've gotten right
                           # this many times
    "recheck" => \$RECHECK # re-check skipped examples in this many times
);
my %config = ( "fmin" => $FMIN, "fmax" => $FMAX, "emin" => $EMIN, "emax" => $EMAX, "iters" => $ITERS, "word" => $WORD );

if(@ARGV != 4) {
    print STDERR "Usage: disc-train.pl FFILE EFILE PFILE WFILE\n";
    exit 1;
}

open FILE, "<:utf8", $ARGV[0] or die "$ARGV[0]: $!\n";
my @fcorp = map { chomp; $_ } <FILE>;
close FILE;
open FILE, "<:utf8", $ARGV[1] or die "$ARGV[1]: $!\n";
my @ecorp = map { chomp; $_ } <FILE>;
close FILE;
open FILE, "<:utf8", $ARGV[2] or die "$ARGV[2]: $!\n";
my (%fprobs,%probs);
while(<FILE>) {
    chomp; my ($f,$e,$p) = split(/\t/);
    $probs{"$f\t$e"} = $p; 
    $fprobs{$f} = [] if not $fprobs{$f}; push @{$fprobs{$f}}, [$e,$p];
}
close FILE;

# perform ITERS iterations
my (%weights, @consec);
my $sents = @fcorp;
foreach my $iter ( 1 .. $ITERS ) {
    my ($correct,$tot);
    # start the iteration
    print STDERR "Starting iteration $iter\n";
    my %corr;

    # for each sentence
    my @order = shuffle(0 .. $#fcorp);
    foreach my $sent (@order) {

        # skip if we've done well on this many times
        $consec[$sent] = $INAROW-1 if($consec[$sent] == $INAROW + $RECHECK);
        if($consec[$sent] >= $INAROW) { $consec[$sent]++; next; }
        $tot++;

        # do decoding
        my ($goodfeats, $goodseq, $goodscore, $goodstr) = forcedecode($fcorp[$sent],$ecorp[$sent],\%weights,\%probs,\%config);
        my ($testfeats, $testseq, $testscore, $teststr) = decode($fcorp[$sent],\%weights,\%fprobs,\%config);

        if($goodstr ne $teststr) {
            # if($iter == $ITERS) { print "WRONG: @$goodseq --> @$testseq\n"; }
            hashpluseq($goodfeats, $testfeats, -1);
            hashpluseq(\%weights, $goodfeats, 1);
            $consec[$sent] = 0;
            if($weights{"N0 <s>"}) {
                print STDERR " good: $goodstr, @$goodseq ($goodscore)\n"; 
                print STDERR " test: $teststr, @$testseq ($testscore)\n"; 
                exit;
            }
            # print STDERR " weights: ".join(" ||| ", map { "$_=>".$weights{$_} } sort keys %weights)."\n";
        } else {
            # print STDERR " correct\n";
            $correct++;
            $consec[$sent]++;
        }

        if($sent % 1000 == 0) {
            print STDERR ".";
        }

    }

    print STDERR "\nIteration $iter finished: $correct/$tot (".($correct/$tot*100)."%), writing weights\n";

    open FILE, ">:utf8", $ARGV[3] or die "$ARGV[3]: $!\n";
    for(sort keys %weights) {
        print FILE "$_\t$weights{$_}\n" if $weights{$_};
    }
    close FILE;

}
