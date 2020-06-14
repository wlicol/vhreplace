#!/usr/bin/perl

###########################################################
# use either intermediate file or canonical intermediate file as 
# input # determine whether fp-covering regions are disjoint 
# (indication of multiple replacements)
###########################################################

if($#ARGV !=0){
 print "[.pl] (canonical)intermediate-file\n";
 exit;
}

$file=$ARGV[0];
if(-e $file){ if($file=~ m/\.gz/){ open(IN,"gunzip -c $file | "); }else{ open(IN,"< $file"); } }
else{ print "couldn't find file: $file\n"; exit; } 

# if no header, change HEADER value to 0
$header=1;
# if the delimitor is comma, change it to ',';
$delimit='\s+';


$line=0;
while($_=<IN>){
 if($line < $header){}
 else{
  chop;
  @a=split(/$delimit/, $_);
  if($a[$#a] =~ m/fp/){

 #--------- for each position in the n1/n2 seq, link to a list of genes whose fp cover it
 $n12seq=$a[2];
 $n12seqL=length($n12seq);

 # 1-init
 for($i=1; $i <=  $n12seqL; $i++){ $linkedG[$i]=""; }

print "$line $a[0] $a[2] $a[$#a] ";
# print "\n";
 #--------- going through all fp's
 @b=split(/\//, $a[$#a]);

 # example of a fp info: fp=aggga;pos=2,8,[2];vg=7
 for($i=0; $i<=$#b; $i++){
  
  # $1 is fp string, $2 is start position(s); $3 is vg position(s)
  # ighv genes can be represented as 3-7, or 7 
  if( $b[$i] =~ /fp=(\w+);pos=(.+)\[\d+\];vg=(.+)/){
  $ha1=$1;
  $vgl=$3;
  $fpL=length($ha1);
  $ha2 = $2;
  $ha2 =~ s/\,$//;
  # print "\ntest => $b[$i]. 1=$ha1 ($fpL) 2=$ha2 3=$vgl\n";



  if($ha2 =~ m/,/){
   # more complicated. need to go through each matching position
   @startA= split(/,/, $ha2);
   for($ii=0; $ii<=$#startA; $ii++){
    $start= $startA[$ii];
    $end= $start + $fpL-1;
    for($j=$start; $j<=$end; $j++){ $linkedG[$j] .= $vgl.",";}
   }

  }else{
   # straightforward 
   $start= $ha2;
   $end = $start +$fpL-1;
   for($j=$start; $j<=$end; $j++){ $linkedG[$j] .= $vgl.",";}
  }


  } # if fp is consistent with the regex pattern


  
  #---------- for each fp, scan all positions, add VG position info to $linkedG

 }

 for($k=1; $k<=$n12seqL; $k++){
  $linkedG[$k] =~ s/,$//;
  #--- get a non-redundant list
  @b=split(/,/, $linkedG[$k]);
  %count=(); for($kk=0; $kk<=$#b; $kk++){ 
   if($b[$kk] =~ m/\-/){ 
    @bb=split(/\-/, $b[$kk]);
    $count{ $bb[1]}++;
   }else{
    $count{ $b[$kk]}++; 
   }
  }
  $newlist=""; for $kkk (sort {$a <=> $b} keys %count){ $newlist .= $kkk.","; } $newlist=~ s/,$//;
  $linkedG[$k] = $newlist;
  # print "\t$k $linkedG[$k]\n";
 }

 # check matching pattern
 $up=$down=0;
 $past= -1;
 for($k=1; $k<=$n12seqL; $k++){
  if(length($linkedG[$k]) >=1){
   if($past ==0){ $up++;}
   print "*";
   $past=1;
  }else{
   if($past ==1){ $down++;}
   print "_";
   $past=0;
  }
 }
 $both=$up+$down;
 print " up=$up down=$down both=$both ";
 if($up > 1 | $down > 1){ print " MUL ";}
 print "\n";

  } #only if there is a fp-match
} $line++;}

