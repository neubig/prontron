#!/usr/bin/perl

use strict;
use utf8;
use List::Util qw(sum min max shuffle);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
binmode STDERR, ":utf8";

my %beg = (
"か" => [ "が" ], "き" => [ "ぎ" ], "く" => [ "ぐ" ], "け" => [ "げ" ], "こ" => [ "ご" ], "さ" => [ "ざ" ], "し" => [ "じ" ], "す" => [ "ず" ], "せ" => [ "ぜ" ], "そ" => [ "ぞ" ], "た" => [ "だ" ], "つ" => [ "づ" ], "て" => [ "で" ], "と" => [ "ど" ], "ち" => [ "じ", "ぢ" ], "は" => [ "ば", "ぱ" ], "ひ" => [ "び", "ぴ" ], "ふ" => [ "ぶ", "ぷ" ], "へ" => [ "べ", "ぺ" ], "ほ" => [ "ぼ", "ぽ" ] 
);

my %end = ( "つ" => [ "っ" ], "く" => [ "っ" ] );

my %prons;
while(<STDIN>) {
    chomp;
    my @arr = split(/\t/);
    my %myp = ( $arr[1] => 1 );
    for(keys %myp) {
        my @myarr = split(/ /);
        if($beg{$myarr[0]}) {
            for(@{$beg{$myarr[0]}}) {
                $myarr[0] = $_;
                $myp{"@myarr"}++;
            }
        }
    }
    for(keys %myp) {
        my @myarr = split(/ /);
        if(@myarr > 1 and $end{$myarr[-1]}) {
            for(@{$end{$myarr[-1]}}) {
                $myarr[-1] = $_;
                $myp{"@myarr"}++;
            }
        }
    }
    $prons{$arr[0]} = {} if not $prons{$arr[0]};
    for(keys %myp) {
        $prons{$arr[0]}->{$_}++;
    }
}

foreach my $k (sort keys %prons) {
    foreach my $v (sort keys %{$prons{$k}}) {
        print "$k\t$v\n";
    }
}
