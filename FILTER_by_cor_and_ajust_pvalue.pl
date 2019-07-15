#!/usr/bin/perl -w

##This script is filtering a correlation pariwisse file  based on its third (cor) and fifth col (ajust_pvalue). 
# perl FILTER_by_cor_and_ajust_pvule.pl 0.5 0.01 input_file | wc -l

 
use strict;

unless( @ARGV == 3 ){

   die "\n\nThe idea is filter the input file based on 2 cut-offs values;  1). cor (third column)  and 2). Ajusted pvalue (fifth column)\n\nenter correlation cut-off, space, pajust cuttof, space  and then input file name\n\n";
}

my $cutoff1 = shift @ARGV;
my $cutoff2 = shift @ARGV;
my $file =  shift @ARGV ;



open IN , $file or die "Can't open $file for reading: $!\n";

while(my $line = <IN> ){

   chomp( $line );

   my @line = split(/\,/, $line);
     

if ($line[0] eq 'id1' ) {
         next ;
     }

   

 else {  print $line[0],",", $line[1],",", $line[2], ",", $line[4],"\n" if ($line[2]) >= ( $cutoff1 )&&($cutoff2)>=($line[4]);
    }
 
                 }
