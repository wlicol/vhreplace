#!/usr/bin/perl

###########################################################
# the input file (the name to be read from the command line)
# should contain a sequence-name line start by
# >
# then followed by (potentially) multiple lines of DNA sequence.
# the ending sequence is defined by those after the RSS (tactgtg).
#
# Usage: [.pl] IGHV-sequence-file
# 
# further filtering by IGHV names (e.g.)
#     [.pl] IGHV-sequence-file |sort -t "-" -k 2,2n
# to remove names which contain "(", then sorted by their chromosome location.
###########################################################

if($#ARGV !=0){
 print "[.pl] ighv_allP.txt\n";
 exit;
}

$f=$ARGV[0];
if(-e $f){ if($f=~ m/\.gz/){ open(IN,"gunzip -c $f | "); }else{ open(IN,"< $f"); } }


###########################################################
# list of pseudo IGHV genes "not in frame" (136 of them)
###########################################################
@PnotIF= qw(IGHV7-NL1*01 IGHV7-NL1*02 IGHV7-NL1*03 IGHV7-NL1*04 IGHV7-NL1*05 
IGHV1/OR16-1*01 IGHV(II)-1-1*01 IGHV1/OR16-2*01 IGHV(III)-2-1*01 IGHV(III)-2-1*02 
IGHV1/OR16-3*01 IGHV1/OR16-4*01 IGHV1/OR16-4*02 IGHV(III)-5-1*01 IGHV(III)-5-2*01 
IGHV(III)-5-2*02 IGHV1/OR15-6*01 IGHV1/OR15-6*02 IGHV3-6*01 IGHV3-6*02 IGHV3/OR16-6*01 
IGHV3/OR15-7*04 IGHV3/OR16-7*01 IGHV3/OR16-7*02 IGHV3/OR16-7*03 IGHV3-11*02 IGHV3/OR16-11*01 
IGHV(III)-11-1*01 IGHV1-12*01 IGHV1-12*02 IGHV(III)-13-1*01 IGHV(III)-13-1*02 
IGHV1-14*01 IGHV1-14*02 IGHV(II)-15-1*01 IGHV(III)-16-1*01 IGHV(III)-16-1*02 IGHV1-17*01 
IGHV1-17*02 IGHV(II)-22-1*01 IGHV(II)-22-1*02 IGHV(II)-22-1D*01 IGHV(III)-22-2*01 
IGHV(III)-22-2D*01 IGHV(III)-25-1*01 IGHV(III)-25-1*02 IGHV(II)-26-2*01 IGHV(III)-26-1*01 
IGHV(III)-26-1*02 IGHV7-27*01 IGHV(II)-28-1*01 IGHV(II)-28-1*02 IGHV(II)-28-1*03 
IGHV(II)-30-1*01 IGHV(II)-30-1*02 IGHV(II)-30-21*01 IGHV(II)-30-32*01 IGHV(II)-30-41*01 
IGHV(II)-30-51*01 IGHV(II)-30-51*02 IGHV(II)-31-1*01 IGHV(II)-33-1*01 IGHV3-36*01 
IGHV3-36*02 IGHV3-36*03 IGHV3-37*01 IGHV3-37*02 IGHV(III)-38-1*01 IGHV(III)-38-1*02 
IGHV(III)-38-1D*01 IGHV7-40*01 IGHV7-40*04 IGHV7-40D*01 IGHV(II)-40-1*01 IGHV3-41*01 
IGHV3-42*01 IGHV3-42*02 IGHV3-42*03 IGHV3-42D*01 IGHV3-42D*04 IGHV(II)-43-1*01 IGHV(II)-43-1*02 
IGHV(II)-43-1D*01 IGHV(II)-43-1D*02 IGHV(II)-44-2*01 IGHV(III)-44*01 IGHV(III)-44*02 
IGHV(III)-44D*01 IGHV(IV)-44-1*01 IGHV(IV)-44-1*02 IGHV(II)-46-1*01 IGHV3-47*03 
IGHV(III)-47-1*01 IGHV(II)-49-1*01 IGHV3-50*01 IGHV(II)-51-2*01 IGHV(II)-53-1*01 
IGHV(II)-53-1*02 IGHV3-54*03 IGHV7-56*01 IGHV7-56*02 IGHV3-57*01 IGHV3-57*02 IGHV3-57*03 
IGHV3-60*01 IGHV3-60*02 IGHV3-60*03 IGHV(II)-60-1*01 IGHV(II)-60-1*02 IGHV3-62*02 
IGHV(II)-62-1*01 IGHV(II)-62-1*02 IGHV(II)-62-1*03 IGHV3-65*01 IGHV3-65*02 IGHV3-65*03 
IGHV(II)-65-1*01 IGHV1-67*01 IGHV1-67*02 IGHV(II)-67-1*01 IGHV(III)-67-2*01 IGHV(III)-67-3*01 
IGHV(III)-67-4*01 IGHV(II)-74-1*01 IGHV3-75*01 IGHV3-76*01 IGHV3-76*02 IGHV(III)-76-1*01 
IGHV5-78*02 IGHV(II)-78-1*01 IGHV(II)-78-1*02 IGHV3-79*01 IGHV3-79*02 IGHV4-80*01 IGHV4-80*02 IGHV(III)-82*01);

foreach $item (@PnotIF){ $nif{ $item}=1; }

###########################################################
# go through the sequence file. save the result in
# %seq: $seq{ $ighv-gene-name}
###########################################################

$line=0;
while($_=<IN>){
  chop;

#### ighv gene information line. should contain ">"
 if($_ =~ m/>/){

if($line >0){
  # summarize the previous result
  $seq{ $name} = $seq_tmp;
}

  $_=~ s/>//;

  # IMGT (http://www.imgt.org) FAST header format
  # The FASTA header contains 15 fields separated by '|':
  #1. IMGT/LIGM-DB accession number(s)
  #2. IMGT gene and allele name
  #3. species      (Homo sapiens)
  #4. IMGT allele functionality (P)
  #5. exon(s), region name(s), or extracted label(s) (V-REGION)
  #6. start and end positions in the IMGT/LIGM-DB accession number(s)
  #7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
  #8. codon start, or 'NR' (not relevant) for non coding labels
  #9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
  #10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
  #11. +n, -n, and/or nS: num of added, deleted, and/or substituted nt to correct seq errors, or 'not corrected' if non corrected seq errors
  #12. number of amino acids (AA): this field indicates that the sequence is in amino acids
  #13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
  #14. partial (if it is)
  #15. reverse complementary (if it is)
  @a=split(/\|/, $_);

  # not-in-frame pseudo is label as "(P")
  if($nif{ $a[1]} ==1){
   $name = $a[1]." (P)";
  }else{
   $name = $a[1]." ".$a[3];
  }
  $seq_tmp="";
#### sequence itself  
 }else{ $seq_tmp .=$_; }

$line++; }

# miss the last one
$seq{ $name} = $seq_tmp;


###########################################################
# go though all sequences
###########################################################
$RSS="tactgtg";
for $k (keys %seq){
 ($beforeRSS, $afterRSS)=split(/$RSS/, $seq{$k});
 if(length( $afterRSS) > 0){
  print "$k $afterRSS\n";
 }
 
}


