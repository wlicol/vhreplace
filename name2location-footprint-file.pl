#!/usr/bin/perl

###########################################################
# simplify the IGHV gene annotation in the footprint file.
# instead of listing full IGHV gene name (gene-location*allele),
# only the location information is kept (if the location info is available).
#
# input footprint file has this format:
# (footprint-seq) (num-time) (footprint-length) (full-list-of-IGHV-genes-with-this-fp)
# caaaagata 5 9 IGHV3-9*01-F,IGHV3-9*02-F,IGHV3-9*03-F,IGHV3-43*01-F,IGHV3-43*02-F,
# ... 
# output will look like this (note the second column value is changed from 5 to 2,
# because only two locations show this fp):
# caaaagata 2 9 3-9,3-43
# ...
###########################################################

if($#ARGV !=1){
 print "[.pl] foot-print-file(e.g.fimr-functional-noD-fp.txt) 0/1(without-with-gene-fam-number-in-output)\n";
 exit;
}

$f=$ARGV[0];
if(-e $f){ if($f=~ m/\.gz/){ open(IN,"gunzip -c $f | "); }else{ open(IN,"< $f"); } }
$fl=$ARGV[1];

###########################################################
# e.g., a[3] = IGHV7-4-1*05-F, IGHV2-70D*04-F
# b[0]= IGHV7-4-1*05-F => IGHV7-4.1*05-F => c[0]=7 and d[0]=4.1
# b[1]= IGHV2-70D*04-F =>  c[0]=2 and d[0]=70D (do not consider). similarly do no consider NL1,...
# at the moment, atypical c[0] is kept (e.g. IGHV4/OR15, IGHV(III), etc.)

while($_=<IN>){
 my %flL=();
 chop;
 @a=split(/\s+/, $_);
 $a[$#a] =~ s/,$//;
 @b=split(/\,/, $a[$#a]);

 for($i=0; $i<=$#b;$i++){
  $b[$i] =~ s/\-/\./g;
  $b[$i] =~ s/\./\-/;
  @c=split(/\-/, $b[$i]);

  $c[0] =~ s/IGHV//;
  @d=split(/\*/, $c[1]);

  if($d[0] =~ m/[a-zA-Z]/){}
  else{
   $flL{ $d[0]} = $c[0];
  }

  # print "$i $b[$i] $c[0] $d[0]\n";
 }

print "$a[0] ";
$nl=0; $listL="";

if($fl==0){ for $k (sort {$a <=> $b} keys %flL){ $nl++; $listL .="$k,"; } }
elsif($fl==1){ for $k (sort {$a <=> $b} keys %flL){ $nl++; $listL .="$flL{$k}-$k,"; } }
$listL =~ s/,$//;
print " $nl $a[2] $listL\n";
}

