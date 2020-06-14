#!/usr/bin/perl

###########################################################
# this is one of the few program to parse the intermediate result
# which may look like this:
# (line number) (IGHV_gene_location_where_n1/n2_from) (footprint-matching-information)
#  16  43 ggcggaaaaagtattgtag      fp=aaaag;pos=7,[1];vg=9,43/fp=aagta;pos=9,[1];vg=34.1/
# ...
# 26590  48 tgcggttaca
# 26602  53 tccctggag     fp=ctggag;pos=4,[1];vg=75/fp=ctgga;pos=4,[1];vg=75/fp=tggag;pos=5,[1];vg=75/
# ...
# more explanation on the last column:
# footprint aaaag is present in this n1 sequence, starting from position 7 (to position 11=7+length(aaaag)-1). 
# this fp appear only once ([1]). this fp can be traced to the ending seq of IGHV genes with locations 9, 43.
# next matched footprint is aagta,...
#
# this program will remove any "initial" IGHV locations (e.g., 43) which is larger or equal to the current IGHV location (=43).
# resulting in a new line
# 16  43 ggcggaaaaagtattgtag      fp=aaaag;pos=7,[1];vg=9,43/fp=aagta;pos=9,[1];vg=34.1/
# 26590  48 tgcggttaca
# 26602  53 tccctggag
# ...
# those fp remain as said to be "canonical". non-canonical fp's are removed.
###########################################################


if($#ARGV !=0){
 print "[.pl] intermediate-file\n";
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
 if($line <  $header){
	print "$_";
 }else{ 
 chop;
 @a=split(/$delimit/, $_);

# there are several situations:
# the final IGHV can be either 1-69 or 69
# the final IGHV can have multiple genes: 3-30,3-30.5
 if($a[$#a] =~ m/fp/){
  $finalL="";

  # -------------- process final IGHV location
  if($a[1] =~ m/,/){
   @c=split(/,/, $a[1]);
   for($i=0; $i<=$#c; $i++){
    if($c[$i] =~ m/\-/){ @b=split(/\-/, $c[$i]); $finalL .= $b[1].","; }
    else{ $finalL .= $c[$i].","; }
   }
   $finalL =~ s/,$//;
  }else{ 
   if($a[1] =~ m/\-/){ @b=split(/\-/, $a[1]); $finalL=$b[1]; }
   else{ $finalL=$a[1]; }
  }

print "$a[0]  $a[1] $a[2] ";
  # -------------- with comma, (e.g.), 4,59 -> make it easy to pick the largest location number (59)
  if($finalL =~ m/,/){
   @b=split(/,/, $finalL);
   $tmp=0; for($i=0; $i<=$#b; $i++){ if($b[$i] > $tmp){ $tmp=$b[$i];} } $finalL=$tmp;
  }

  # -------------- go through each footprint
  @b=split(/\//, $a[$#a]);
  # go through each fp
  for($i=0; $i<=$#b; $i++){
   @c=split(/vg=/, $b[$i]);
   # second part, c[1], contains locations (with "-" or not, with "," or not)
   @d=split(/,/, $c[1]);

   $survL="";
   for($j=0; $j<=$#d; $j++){
    if($d[$j] =~ m/\-/){
     @e=split(/\-/, $d[$j]);
     if($e[1] < $finalL){ $survL .= $d[$j].","; }
    }else{
     if($d[$j] < $finalL){ $survL .= $d[$j].",";}
    }
   }
if( length($survL)>0){
 $survL =~ s/,$//;
 print "$c[0]vg=$survL\/";
}
  }
print "\n";

 }else{
  # don't make any changes
  print "$_\n";
 }

} $line++; }

