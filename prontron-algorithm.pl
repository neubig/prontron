
use strict;
use utf8;

my %phones = (
    "あ" => "a", "い" => "i", "う" => "u", "え" => "e", "お" => "o", "か" => "k a", "き" => "k i", "く" => "k u", "け" => "k e", "こ" => "k o", "が" => "g a", "ぎ" => "g i", "ぐ" => "g u", "げ" => "g e", "ご" => "g o", "さ" => "s a", "し" => "s i", "す" => "s u", "せ" => "s e", "そ" => "s o", "ざ" => "z a", "じ" => "z i", "ず" => "z u", "ぜ" => "z e", "ぞ" => "z o", "た" => "t a", "ち" => "t i", "つ" => "t u", "て" => "t e", "と" => "t o", "だ" => "d a", "ぢ" => "d i", "づ" => "d u", "で" => "d e", "ど" => "d o", "な" => "n a", "に" => "n i", "ぬ" => "n u", "ね" => "n e", "の" => "n o", "は" => "h a", "ひ" => "h i", "ふ" => "h u", "へ" => "h e", "ほ" => "h o", "ば" => "b a", "び" => "b i", "ぶ" => "b u", "べ" => "b e", "ぼ" => "b o", "ぱ" => "p a", "ぴ" => "p i", "ぷ" => "p u", "ぺ" => "p e", "ぽ" => "p o", "ま" => "m a", "み" => "m i", "む" => "m u", "め" => "m e", "も" => "m o", "や" => "y a", "ゆ" => "y u", "よ" => "y o", "ら" => "r a", "り" => "r i", "る" => "r u", "れ" => "r e", "ろ" => "r o", "っ" => "x", "ゃ" => "xy a", "ゅ" => "xy u", "ょ" => "xy o", "ぁ" => "xa", "ぃ" => "xi", "ぅ" => "xu", "ぇ" => "xe", "ぉ" => "xo", "わ" => "w a"
);

sub phoneticmap {
    $_ = shift;
    # print "size: ".%phones."\n";
    return join(" ", map { $phones{$_} ? $phones{$_} : $_ } split(/ /));
}

# transform a symbol into its representation
sub makerep {
    my $str = shift;
    my @ret;
    $_ = $str; s/ /_/g; s/\t/|/g; push @ret, $_; # joint
    my ($f,$e) = split(/\|/); push @ret, $f; push @ret, $e; # f and e
    $e =~ s/_/ /g; push @ret, $e; # e words
    # print "phoneticmap($e) -> ".phoneticmap($e)."\n";
    push @ret, phoneticmap($e);
    return join("\t",@ret);
}
sub hashproduct {
    my ($sparse,$dense) = @_;
    my $ret = 0;
    while(my ($k,$v) = each(%$sparse)) {
        $ret += $v*$dense->{$k};
    }
    return $ret;
}
sub hashpluseq {
    my ($left,$right,$mult) = @_;
    return if not $right;
    while(my ($k,$v) = each(%$right)) {
        $left->{$k} += $v*$mult if $v;
    }
}

# the length over n-gram representations
my @nlen = (2,2,2,2,2);
my $init = join("\t", map { join(" ", map {"<s>"} ( 2 .. $_ )) } @nlen);
# print STDERR "initial state: $init\n";

sub makefeat {
    my ($prev,$next) = @_;
    # # print STDERR "MAKEFEAT:\n prev: $prev\n next: $next\n";
    my @pa = split(/\t/,$prev); push @pa, "" while(@pa < @nlen);
    my @na = split(/\t/,$next); push @na, "" while(@na < @nlen);
    # die "bad arrays @pa, @na" if($#pa != $#na);
    my (%feat,@ncont);
    foreach my $i (0 .. $#pa) {
        my @arr = split(/ /,$na[$i]);
        $feat{"L$i"} += @arr if @arr;
        $feat{"L$i-".@arr}++;
        @arr = (split(/ /,$pa[$i]),@arr);
        my $l = $nlen[$i];
        foreach my $j (0 .. $#arr) {
            my $id = "N$i";
            foreach my $k ($j .. min($#arr,$j+$l-1)) {
                $id .= " ".$arr[$k];
                $feat{$id}++ if $k >= $l-1;
            }
        }
        push @ncont, join(" ",@arr[max(@arr-$l+1,0) .. $#arr]);
        # print STDERR " pa[$i] == $pa[$i] --- arr=@arr, ncont=@ncont\n";
    }
    # # print STDERR " feat: ".join(" ||| ", map { "$_=>".$feat{$_} } sort keys %feat)."\n ncont=".join(" ||| ",@ncont)."\n";
    return (\%feat, join("\t",@ncont));
}

sub forcedecode {
    die "forcedecode requires 5 arguments" if (@_ != 5);
    my ($f,$e,$weights,$cands,$config) = @_;
    my ($FMIN, $FMAX) = ($config->{"fmin"}, $config->{"fmax"});
    my ($EMIN, $EMAX) = ($config->{"emin"}, $config->{"emax"});

    # forward pass for the viterbi forced decoding
    my @fs = ($config->{"word"} ? split(/ /,$f) : split(//,$f));
    my @es = ($config->{"word"} ? split(/ /,$e) : split(//,$e));
    my $fl = @fs; my $el = @es; my $elp1 = ($el+1);
    # make the stacks, [ score, prev, string, feat ]
    my @stacks = ( { $init => [ 0, 0, 0, 0 ] } ); $stacks[($fl+1)*$elp1] = 0;
    foreach my $t ($FMIN .. $fl) {
        foreach my $v ($EMIN .. $el) {
            next if not $t+$v;
            my %nstack;
            foreach my $s (max(0,$t-$FMAX) .. $t-$FMIN) {
                my $fid = join(" ",@fs[$s .. $t-1]);
                foreach my $u (max(0,$v-$EMAX) .. $v-$EMIN) { 
                    my $eid = join(" ",@es[$u .. $v-1]);
                    my $id = "$fid\t$eid";
                    my $pstack = $stacks[$s*$elp1+$u];
                    next if (not exists $cands->{$id}) or (not $pstack);
                    my $rep = makerep($id);
                    while(my ($k,$v) = each(%$pstack)) {
                        my ($nfeat,$ncontext) = makefeat($k,$rep);
                        my $score = hashproduct($nfeat,$weights) + $v->[0];
                        if((not $nstack{$ncontext}) or # if previous doesn't exist
                            ($nstack{$ncontext}->[0] < $score) or # or this is better
                            (($nstack{$ncontext}->[0] == $score) and (rand() < 0.5))) { # or equal at p=0.5
                            $nstack{$ncontext} = [ $score, $v, $id, $nfeat ]; # add to the stack
                        }
                    }
                }
            }
            $stacks[$t*$elp1+$v] = \%nstack;
        }
    }
    my $last = $stacks[$fl*$elp1+$el];
    
    # go backward through to get the forced match
    my (%feats, @seq, $score, $str);
    if(!$last || (%$last == 0)) {  print STDERR " NO VITERBI MATCH at forced for @fs, @es\n"; next; }
    else {
        my $curr = [-1e99,0,0,0];
        while(my($k,$v) = each(%$last)) {
            my ($nfeat,$ncontext) = makefeat($k,$init);
            my $score = hashproduct($nfeat,$weights) + $v->[0];
            if(($curr->[0] < $score) or (($curr->[0] == $score) and (rand() < 0.5))) {
                $curr = [ $score, $v, 0, $nfeat ]; # add to the stack
            }
        }
        $score = $curr->[0];
        while($curr) {
            # print STDERR " ".join(" ",@$curr)."\n";
            # print STDERR " feat: ".join(" ||| ", map { "$_=>".$curr->[3]->{$_} } sort keys %{$curr->[3]})."\n" if $curr->[3];
            hashpluseq(\%feats,$curr->[3], 1);
            unshift @seq, $curr->[2] if $curr->[2];
            $curr = $curr->[1];
        }
        $str = join(" ",map { my $str = $_; $str =~ s/.*\t//g; $_ ? split(/ /,$str) : () } @seq);
        # print STDERR " feat: ".join(" ||| ", map { "$_=>".$feats{$_} } sort keys %feats)."\n";
    }

    return (\%feats, \@seq, $score, $str);

}

sub decode {
    die "decode requires 4 arguments" if (@_ != 4);
    my ($f,$weights,$fcands,$config) = @_;
    my ($FMIN, $FMAX) = ($config->{"fmin"}, $config->{"fmax"});

    # get the viterbi non-forced decoding
    # make the stacks, [ score, prev, string, feat ]
    my @fs = ($config->{"word"} ? split(/ /,$f) : split(//,$f));
    my $fl = @fs;
    my @stacks = ( { $init => [ 0, 0, 0, 0 ] } ); $stacks[$fl+1] = 0;

    # forward pass for the viterbi non-forced decoding
    foreach my $t ($FMIN .. $fl) {
        my %nstack = ( $t == 0 ? %{$stacks[0]} : () );
        foreach my $s (max(0,$t-$FMAX) .. $t-$FMIN) {
            my $pstack = $stacks[$s];
            my $fid = join(" ",@fs[$s .. $t-1]);
            if ((not $fcands->{$fid}) and ($t-$s==1)) {
                print STDERR "WARNING: Unknown character $fid\n";
                $fcands->{$fid} = [ [ $fid, 0 ] ];
            }
            next if (not exists $fcands->{$fid}) or (not $pstack);
            foreach my $eid (map { $_->[0] } @{$fcands->{$fid}}) {
                my $id = "$fid\t$eid";
                # # print STDERR "($s,$t) cands{$id}=$cands{$id}\n";
                my $rep = makerep($id);
                while(my ($k,$v) = each(%$pstack)) {
                    my ($nfeat,$ncontext) = makefeat($k,$rep);
                    my $score = hashproduct($nfeat,$weights) + $v->[0];
                    if((not $nstack{$ncontext}) or # if previous doesn't exist
                        ($nstack{$ncontext}->[0] < $score) or # or this is better
                        (($nstack{$ncontext}->[0] == $score) and (rand() < 0.5))) { # or equal at p=0.5
                        $nstack{$ncontext} = [ $score, $v, $id, $nfeat ]; # add to the stack
                    }
                }
            }
        }
        $stacks[$t] = \%nstack;
    }
    my $last = $stacks[$fl];

    # go backward through to get the unforced match
    my (%feats, @seq, $score, $str);
    if(!$last or (%$last == 0)) {  print STDERR " NO VITERBI MATCH at unforced for @fs\n"; next; }
    else {
        my $curr = [-1e99,0,0,0];
        while(my($k,$v) = each(%$last)) {
            my ($nfeat,$ncontext) = makefeat($k,$init);
            my $score = hashproduct($nfeat,$weights) + $v->[0];
            if(($curr->[0] < $score) or (($curr->[0] == $score) and (rand() < 0.5))) {
                $curr = [ $score, $v, 0, $nfeat ]; # add to the stack
            }
        }
        $score = $curr->[0];
        while($curr) {
            # print STDERR " ".join(" ",@$curr)."\n";
            # print STDERR " feat: ".join(" ||| ", map { "$_=>".$curr->[3]->{$_} } sort keys %{$curr->[3]})."\n" if $curr->[3];
            hashpluseq(\%feats,$curr->[3], 1);
            unshift @seq, $curr->[2] if $curr->[2];
            $curr = $curr->[1];
        }
        $str = join(" ",map { my $str = $_; $str =~ s/.*\t//g; split(/ /,$str) } @seq);
        # print STDERR " feat: ".join(" ||| ", map { "$_=>".$feats{$_} } sort keys %feats)."\n";
    }

    return (\%feats, \@seq, $score, $str);

}

1
