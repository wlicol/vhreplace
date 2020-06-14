#!/usr/bin/perl

###########################################################
# similar to Zhixin Zhang's program (Front. Immunol. (2014), 5:40
# given a list of footprints (substring  (L>=5) of ending sequences: 
# the portion of IGHV sequences after RSS motif), and a list of n1 n2 
# sequences, determine the  fp-matching status. the result is written 
# to an intermediate file (with information on how many times a 
# fp-matching event occurs  and where the matching occurs).
#
# another program will parse the  intermediate file to get statistics easily.
#
# Usage: 
# [.pl (processed) sequence_file footprint_file n1/n2
# where
# processed_sequence_file: it contains n1 or/and n2 sequence information
# footprint_file: a list of 5(or longer)-mers that are subsequences of the
#        IGHV ending sequences
# n1/n2: choice of analyzing n1 or n2 sequences
#
# format of the processed_sequence_file:
# ---------------------------------------------------------------------------------------
# sequence_id v-gene_and_allele n1-region n2-region v-domain_functionality subset_missing
# AT-01-0001-H1 Homsap_IGHV1-69*01_F_or_Homsap_IGHV1-69D*01_F  NA ttccta productive NA
# ...
#
# format for footprint_file
# (fp_seq) (num_location) (fp_length) (location_list)
# ---------------------------------
# caaaagata  2 9 9,43
# caaggagac  1 9 10
# ...
###########################################################

if($#ARGV !=3){
 print "[.pl] (processed)sequence_file footprint_file n1/n2 0/1(without/with-gene-fam-num)\n";
 exit;
}

if( ($ARGV[2] eq "n1") | ($ARGV[2] eq "n2") ){
 $fn12 = $ARGV[2];
}else{
 print "$ARGV[2] needs to be either n1 or n2\n";
 exit;
}

$fl=$ARGV[3];
if($fl !=0 & $fl !=1){
 print "$fl should be 0 or 1\n";
 exit;
}

#----------------------------------------------------------------------------------------
# open two input files
#----------------------------------------------------------------------------------------

# standardized n1/n2 seq file
# e.g. joy-CLL26k.txt
$file=$ARGV[0];
if(-e $file){ if($file=~ m/\.gz/){ open(IN1,"gunzip -c $file | "); }else{ open(IN1,"< $file"); } }
else{ print "couldn't find file: $file\n"; exit; } 

# annotated footprint file
# e.g. self-generated-footprints-dec5.txt
$f2 =$ARGV[1];
if(-e $f2){ if($f2=~ m/\.gz/){ open(IN2,"gunzip -c $f2 | "); }else{ open(IN2,"< $f2"); } }
else{ print "couldn't find file: $f2\n"; exit; } 

#----------------------------------------------------------------------------------------
# save the inputs to arrays
# "guessing" whether it's space/tab/comma delimtied [by counting the number of space/tab/comma] 
# "guessing" whether it has header [by searching n1, n2 name]
# if this guessing does not work, then the code will be rewritten
#----------------------------------------------------------------------------------------

# file1
$stext=""; while($_=<IN1>){ $stext .= $_; } close(IN1);

$nspace= ($stext =~ tr/\s+//);
$ntab= ($stext =~ tr/\t//);
$ncomma= ($stext =~ tr/,//);

#sequence delimit type
if($nspace < $ntab){
 if($ntab < $ncomma){ $sDL=","; }
 else{ $sDL= "\\t";}
}else{ $sDL = "\\s+"; }
# header or not, assuming the head contains n1/N1 or n2/N2
@stextA=split(/\n/, $stext);
if( ($stextA[0] =~ m/[nN]1/) | ($stextA[0] =~ m/[nN]2/) ){ $StartLine=1;} else{ $StartLine=0;}

# file2
$fptext=""; while($_=<IN2>){ $fptext .= $_; } close(IN2);
$nspace2= ($fptext =~ tr/\s+//);
$ntab2= ($fptext =~ tr/\t//);
$ncomma2= ($fptext =~ tr/,//);
if($nspace2 < $ntab2){
 if($ntab2 < $ncomma2){ $sDL2=","; }
 else{ $sDL2= "\\t";}
}else{ $sDL2 = "\\s+"; }

# header or not, assuming the head contains fp or finger(print)
@fptextA=split(/\n/, $fptext);
if( ($fptextA[0] =~ m/fp/) | ($fptextA[0] =~ m/finger/) ){ $StartLine2=1;} else{ $StartLine2=0;}


#----------------------------------------------------------------------------------------
# go through text from fp files
#----------------------------------------------------------------------------------------

$fptmp="";
for($i=$StartLine2; $i<=$#fptextA; $i++){
 @a=split($sDL2, $fptextA[$i]);
 $fptmp .= $a[0]." ";
 # clean all space. only keep comma-delimited array
 $a[$#a] =~ s/\s+//g;
 $vgfp{ $a[0]}=$a[$#a];
 # change this to a[1]
 # $uniq{ $a[0]}=$a[2];
 $uniq{ $a[0]}=$a[1];
}

$fptmp=~ s/\s+$//;
$fptmp=~ s/^\s+//;
@fp=split(/\s+/, $fptmp);

#----------------------------------------------------------------------------------------
# go through text from sequence (n1,n2) files
# Note: empty n1 or n2 is indicated by "NA"
# Note: in out footprint definition, the length should be 5 or longer
#----------------------------------------------------------------------------------------

# columns are fixed:
# 0. sequenceID. 1.vgeneName 2.n1 3.n2 4.functionality 5.subset/other

# print headers
print "row V-gene-location $fn12-sequence details\n";

for($i=$StartLine; $i <= $#stextA; $i++){
 @a=split($sDL, $stextA[$i]);

 if($fn12 eq "n1"){
  $n1=$a[2];
 }
 elsif($fn12 eq "n2"){
  $n1=$a[3];
 }

 $n1L=length($n1);
 
#----------- only if n1 is not NA, and 5-mer or longer
if( ($n1 ne "NA") & (length($n1) >=5) ){

 #----------- create a list of non-redundant locations for the sequence gene list
 my %gL=();	# create a list of non-redundant location from gene list
 @b=split(/\_or\_/, $a[1]);
 for($j=0; $j<=$#b; $j++){
  $b[$j] =~ s/Homsap_//;
  $b[$j] =~ s/IGHV//;
  $b[$j] =~ s/\_F//;
  $b[$j] =~ s/\-/\./g;
  $b[$j] =~ s/\./\-/;

  # skip names NL, and remove D
  if($b[$j] =~ m/NL/){}
  else{
   $b[$j] =~ s/D//;
   @c=split(/\-/, $b[$j]);
   @d=split(/\*/, $c[1]);
   $gL{ $d[0] } = $c[0];
  }
 }

  $listL="";
  if($fl==0){ for $k (sort {$a <=> $b} keys %gL){ $listL .="$k,"; } }
  elsif($fl==1){ for $k (sort {$a <=> $b} keys %gL){ $listL .="$gL{$k}-$k,"; } }
  $listL =~ s/,$//;
  
print "$i  $listL $n1\t";

 #----------- search fp
 for($ll=0; $ll<=$#fp; $ll++){
  $fptmp=$fp[$ll];
  $fpLtmp=length($fptmp);
  if($n1L >= $fpLtmp){	# if n1-length is shorter, of course no match

  $matched_fp_pos="";
  $n1tmp=$n1;$npos=0;
  while($n1tmp =~ m/$fptmp/g){ # record all matching positions
   $mpos=pos($n1tmp);
   $npos++;
   $mstart = $mpos-$fpLtmp+1;
   $matched_fp_pos .= $mstart.",";
  }

  if(length($matched_fp_pos)>0){
   $matched_fp_pos =~ s/,$//;
   # ------------ explain what's written here
   # fp= footprint_seq; pos_matching=2,7,[2]; vg=reverse-map-ighv-location-of-that-fp;
   # print "fp=$fptmp;pos=$matched_fp_pos,[$npos];vg=$vgfp{$fptmp};uniq=$uniq{$fptmp}/";
   # don't really need to print uniq{}
   print "fp=$fptmp;pos=$matched_fp_pos,[$npos];vg=$vgfp{$fptmp}/";
  }

  } # end of checking n1 length
 } # end of looping through all fp's


print "\n";


} # non-empty non-NA or short  n1

} # end of looping through sequence file



