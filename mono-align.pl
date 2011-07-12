#!/usr/bin/perl

# A program that performs monolingual alignment by jointly maximizing the
# probability of two conditional models
#  (similar to Percy Liang et al. "Alignment by Agreement")

$| = 1;

use strict;
use utf8;
use List::Util qw(sum min max shuffle);
use Getopt::Long;
use FindBin qw($Bin);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
binmode STDERR, ":utf8";

# require "$Bin/digamma.pl";

my ($FMIN, $FMAX, $EMIN, $EMAX) = (1,1,0,5);
my $ITERS = 10;
my $CUT = 1e-2;
my $WORD = 0;
GetOptions(
    "fmin=i" => \$FMIN, # minimum length of the input
    "fmax=i" => \$FMAX, # maximum length of the input
    "emin=i" => \$EMIN, # minimum length of the output
    "emax=i" => \$EMAX, # maximum length of the output
    "iters=i" => \$ITERS, # maximum number of iterations
    "cut=f" => \$CUT,  # all pairs that do not have a maximum posterior
                       # probability of at least this much will be trimmed
    "word" => \$WORD   # use word units instead of characters
);

if(@ARGV != 3) {
    print STDERR "Usage: mono-align.pl [OPTIONS] FFILE EFILE CANDOUT\n";
    exit 1;
}

# load the corpora and save the names
my %ids = ( "\t" => 0 );
my (%fids,%eids);
open FILE, "<:utf8", $ARGV[0] or die "$ARGV[0]: $!\n";
my @fcorp = map { 
    chomp; 
    # map repeated characters
    if($WORD) { s/(\S* \S*) 々 々/$1 $1/g; s/(.) 々/$1 $1/g; }
    else      { s/(\S*\S*)々々/$1$1/g; s/(.)々/$1$1/g; }
    my @arr = map { 
        $ids{$_} = keys %ids if not $ids{$_}; 
        $fids{$_}++; $ids{$_} 
    } ( $WORD ? split(/ /, $_) : split(//, $_) ); 
    \@arr 
} <FILE>;
close FILE;
open FILE, "<:utf8", $ARGV[1] or die "$ARGV[1]: $!\n";
my @ecorp = map { chomp; my @arr = map { $ids{$_} = keys %ids if not $ids{$_}; $eids{$_}++; $ids{$_} } ( $WORD ? split(/ /, $_) : split(//, $_) ); \@arr } <FILE>;
close FILE;
@fcorp == @ecorp or die "Corpus sizes don't match\n";
my @names;
while(my ($k,$v) = each(%ids)) { $names[$v] = $k; }
my ($fzero,$ezero) = ( -log(1+keys(%fids)), -log(1+keys(%eids)) );

# subroutines
sub printstr {
    my $_ = shift;
    my $str = join(" ",map {$names[$_]} unpack("L".length($_)/4 , $_));
    $str =~ s/ *\t */\t/g;
    return $str;
}
sub logsumexp {
    my $max = max(@_);
    return log(sum(map { exp($_ - $max) } @_))+$max;
}
sub splitid {
    my $_ = shift;
    my @arr = unpack("L".length($_)/4 , $_);
    for($_ = 0; $_ < @arr; $_++) { last if not $arr[$_]; }
    return ( pack("L$_",@arr[0 .. $_-1]), pack("L".($#arr-$_),@arr[$_+1 .. $#arr]) );
}

# perform ITERS iterations
my (%probs, $lastlik, %maxpost);
foreach my $iter ( 1 .. $ITERS ) {

    # start the iteration
    print STDERR "Starting iteration $iter\n";
    my ($lik, %counts);
    %maxpost = ();

    # for each sentence
    foreach my $sent (0 .. $#fcorp) {
        my @fs = @{$fcorp[$sent]}; my @es = @{$ecorp[$sent]};
        # print "e=".join(" ",map{$names[$_]}@es).", f=".join(" ",map{$names[$_]}@fs)."\n";
        my $fl = @fs; my $el = @es; my $elp1 = $el + 1;
        my @alph; $alph[($fl+1)*($el+1)] = 0;
        $alph[$_] = -1e99 for(1 .. ($fl+1)*($el+1));
        my @beta; $beta[($fl+1)*($el+1)] = 0;
        $beta[$_] = -1e99 for(0 .. ($fl+1)*($el+1) - 2);

        # calculate forward probabilities for every valid span (s,t,u,v)
        foreach my $t ($FMIN .. $fl) {
            foreach my $v ($EMIN .. $el) {
                next if not $t+$v;
                my @psum; # store the log probs to add
                foreach my $s (max(0,$t-$FMAX) .. $t-$FMIN) {
                    foreach my $u (max(0,$v-$EMAX) .. $v-$EMIN) {
                        my $len = $t-$s+$v-$u; next if not $len;
                        my $id = pack("L".($len+1),@fs[$s .. $t-1],0,@es[$u .. $v-1]);
                        $probs{$id} = $fzero*($t-$s+1) + $ezero*($v-$u+1) if not exists $probs{$id};
                        push @psum, $alph[$s*$elp1+$u] + $probs{$id};
                        # print STDERR "($s,$t,$u,$v): $alph[$s*$elp1+$u] + $probs{$id}\n";
                    }
                }
                $alph[$t*$elp1+$v] = logsumexp(@psum);
            }
        }
        my $tot = $alph[$fl*$elp1+$el];
        $lik += $tot;

        # calculate the backward probabilities and counts
        for(my $s = $fl-$FMIN; $s >= 0; $s--) {
            for(my $u = $el-$EMIN; $u >= 0; $u--) {
                next if not $fl-$s+$el-$u;
                my @psum; # store the log probs to add
                foreach my $t ($s+$FMIN .. min($fl,$s+$FMAX)) {
                    foreach my $v ($u+$EMIN .. min($el,$u+$EMAX)) {
                        my $len = $t-$s+$v-$u; next if not $len;
                        my $id = pack("L".($len+1),@fs[$s .. $t-1],0,@es[$u .. $v-1]);
                        push @psum, $beta[$t*$elp1+$v] + $probs{$id};
                        my $mylik = $alph[$s*$elp1+$u] + $beta[$t*$elp1+$v] + $probs{$id};
                        # print STDERR "($s,$t,$u,$v) ".$alph[$s*$elp1+$u]." + ".$beta[$t*$elp1+$v]." + ".$probs{$id}."\n" if $iter > 1;
                        die "bad likelihood $mylik (/$tot)" if $mylik-$tot > 0.01;
                        my $post = exp($mylik-$tot);
                        $counts{$id} += $post; $maxpost{$id} = max($maxpost{$id},$post);
                        # print STDERR "counts{".printstr($id)."} @ ($s,$t,$u,$v) += exp($mylik-$tot) -> $counts{$id}\n";
                    }
                }
                $beta[$s*$elp1+$u] = logsumexp(@psum);
            }
        }

        if($sent % 1000 == 0) {
            print STDERR ".";
        }
    }
    
    print STDERR "\nIter $iter finished, likelihood=$lik, normalizing...\n";

    # count bidirectional probs
    my (%fcounts,%ecounts);
    while(my ($k,$v) = each(%counts)) {
        my ($fk,$ek) = splitid($k);
        # print STDERR "k=".printstr($k)." f=".printstr($fk)." e=".printstr($ek)." v=$v\n";
        $fcounts{$fk} += $v; $ecounts{$ek} += $v;
    } 

    # calculate the denominator
    my $denom;
    while(my ($k,$v) = each(%counts)) {
        my ($fk,$ek) = splitid($k);
        next if (not $fcounts{$fk}) or (not $ecounts{$ek});
        my $newv = $v/$fcounts{$fk}*$v/$ecounts{$ek};
        # print STDERR "recalcs ".printstr($k).": $newv = $v/$fcounts{$fk}*$v/$ecounts{$ek}\n";;
        $counts{$k} = $newv;
        $denom += $newv;
    }
    
    # calculate the new probabilities
    while(my ($k,$v) = each(%counts)) {
        $probs{$k} = log(max($counts{$k}/$denom,1e-100));
        die "bad prob ".printstr($k)." = $probs{$k}" if($probs{$k} > 0);
        # print STDERR "probs probs{".printstr($k)."} = log($counts{$k}/$denom) = ".log(max($counts{$k}/$denom,1e-100))."\n";
    }

    # if likelihoods didn't change very much
    last if($iter > 2 && ($lastlik-$lik)/$lastlik < 0.01);
    $lastlik = $lik; 

}

open FILE, ">:utf8", $ARGV[2] or die "$ARGV[2]: $!\n";
foreach my $id (sort keys %probs) {
    if($maxpost{$id} > $CUT) {
        # print FILE printstr($id)."\t$probs{$id}\t$maxpost{$id}\n";
        print FILE printstr($id)."\n"; # "\t$probs{$id}\t$maxpost{$id}\n";
    }
}
close;
