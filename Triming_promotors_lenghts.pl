#!/usr/bin/perl -w

#This script take FASTA sequences as input in one row (.prepared file input) and take 1 Kb upstream

my $file = shift @ARGV;

open IN1, $file or die "Can't read $file: $!\n";

while( my $seqs = <IN1> ){

my @seqs = split(/\n/, $seqs);

foreach $seq (@seqs) {

if ($seq =~ /^>/) {         



$seqend = substr($seq,0,30);
chomp $seqend;
print $seqend,"\n";
} else {
$seqend2 = substr($seq,-1000);
chomp $seqend2;
print $seqend2,"\n";

}
}
}

