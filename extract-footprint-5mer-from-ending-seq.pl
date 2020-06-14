#!/usr/bin/perl

###########################################################
# extrac "footprint" (subsequence with length >= 5) from the ending sequence
# the input ending sequence file contains three columns
# IGHV-gene-name functionality-info ending-seq
# 
# usage:
#  [.pl] ending-seq-file |sort -k 3,3nr -k 1,1
# to sort the list by seq length (from longest to shortest). 
# or
#  [.pl] ending-seq-file |sort -k 1,1i
# to sort the list of alphabet pattern
#
# to get fp from a subset of ending seq
# (e.g., F(functional) without D)  cat ending-seq-file|grep F |grep -v ORF |grep -v D > tmp; [.pl] tmp |sort -k 3,3nr -k 1,1
# (e.g. those with location) grep -v "\/OR" ending-seq-file | grep -v D | grep -v "V(" | grep -v NL > tmp; [.pl] tmp |sort -k 3,3nr -k 1,1
###########################################################

if($#ARGV !=0){
 print "[.pl] ending-seq.txt(e.g.,fimr-all-ending.txt)\n";
 exit;
}

$f=$ARGV[0];
if(-e $f){ if($f=~ m/\.gz/){ open(IN,"gunzip -c $f | "); }else{ open(IN,"< $f"); } }

$header=0; # if there is a header, change the value to 1
$delimit='\s+'; # if it's comma-delimited, change to ',' (note: double quotation doesn't work)

$line=0;
while($_=<IN>){
 if($line >= $header){
 chop; 
 @a=split(/$delimit/, $_);

 $end_tmp=$a[$#a];
 $Lend_tmp=length($end_tmp);
 if($end_tmp =~ m/n/){}	# doesn't use the ending seq with "n" in it
 else{

  for($L=5; $L<=$Lend_tmp; $L++){
   $i=0;
   while( length($tmp =substr($end_tmp, $i,$L))== $L){
    $count{$ tmp} ++;
    $gene{$tmp} .= $a[0]."-".$a[1].",";
    $i++;
   }
  }

 }
} $line++; }



###########################################################
# summarize the sub-seq (footprint) result
###########################################################

for $k (keys %count){
 $lfp=length($k);
 print "$k $count{$k} $lfp $gene{$k}\n";
}
