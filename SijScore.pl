#!usr/bin/perl -w -s
use strict;
use warnings;


BEGIN{
	push(@INC,"G:\\Sij");
}

my %scoring_matrix;
my @aminoacid= qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V X);
my $TotalfrequecyAmbiguity=0;
my $totalsubstitutions=0;
my $TotalQijExcSimilar=0;


for(my $i=0;$i<20;$i++){
	for(my $j=0;$j<20;$j++){
		 $scoring_matrix{$aminoacid[$i].$aminoacid[$j]}=0;
	}
}


my @Blocks;

open(ST,"<Block.mblk") or die($!);

foreach(<ST>){
	chomp($_);
	push(@Blocks,$_);
}


my $blocklength=@Blocks;
my $coloumnlength=length($Blocks[0]);

print "Coloumn length : $coloumnlength\n";

my @Qij=0;

my $divideQij=$coloumnlength*$blocklength*(($blocklength-1)/2);


#amino acid hash

my %amino;

foreach(@aminoacid){
	$amino{$_}=0;
}

print $coloumnlength." ".$blocklength."\n";

#Calculating Frequencies of each Amino acid substitutions Coloumn wise 

for(my $col=0;$col<$coloumnlength-2;$col++){
	for(my $row1=0;$row1<$blocklength-1;$row1++){	#reference row
		my $aminorow1=substr($Blocks[$row1],$col,1);
		
		for(my $row2=0;$row2<$blocklength-1;$row2++){ #traversing row
			my $aminorow2=substr($Blocks[$row2],$col,1);
			
			#print $aminorow1." ".$aminorow2." ".$row1." - ".$row2."\n";
				$scoring_matrix{$aminorow1.$aminorow2}+=1;	
				
				if($aminorow1 eq 'X' || $aminorow2 eq 'X'){
					$TotalfrequecyAmbiguity-=1;
				}else{
					$totalsubstitutions+=1;
				}
		}	
	}
	#last;
}

print "Total substitutions : ".$totalsubstitutions."\n";
print "Total substitutions : ".$divideQij."\n";

#print each frequencies after it

#Calculating Qij for each amino acids 

my %fuij;	#Hash that maintains summation Qi of all Qij1 Qij2 Qij3....

foreach(@aminoacid){
	$fuij{$_}=0;
}

foreach(@aminoacid){
	my $aminoacidatfirst=$_;

	foreach(@aminoacid){
		my $aminoacidatsecond=$_;
		$fuij{$aminoacidatfirst}+=$scoring_matrix{$aminoacidatfirst.$aminoacidatsecond};
		#$totalsubstitutions+=$scoring_matrix{$aminoacidatfirst.$aminoacidatsecond};
	}
}


foreach(keys %fuij){
	print $_." --- > ".$fuij{$_}."\n";
}



print $scoring_matrix{$aminoacid[1].'F'};

#printing the matrix


print "\nFij Matrix\n";

print "\n \tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\n";
print "A\t".$scoring_matrix{"AA"}."\t".$scoring_matrix{"AR"}."\t".$scoring_matrix{"AN"}."\t".$scoring_matrix{"AD"}."\t".$scoring_matrix{"AC"}."\t".$scoring_matrix{"AQ"}."\t".$scoring_matrix{"AE"}."\t".$scoring_matrix{"AG"}."\t".$scoring_matrix{"AH"}."\t".$scoring_matrix{"AI"}."\t".$scoring_matrix{"AL"}."\t".$scoring_matrix{"AK"}."\t".$scoring_matrix{"AM"}."\t".$scoring_matrix{"AF"}."\t".$scoring_matrix{"AP"}."\t".$scoring_matrix{"AS"}."\t".$scoring_matrix{"AT"}."\t".$scoring_matrix{"AW"}."\t".$scoring_matrix{"AY"}."\t".$scoring_matrix{"AV"}."\n";
print "R\t".$scoring_matrix{"RA"}."\t".$scoring_matrix{"RR"}."\t".$scoring_matrix{"RN"}."\t".$scoring_matrix{"RD"}."\t".$scoring_matrix{"RC"}."\t".$scoring_matrix{"RQ"}."\t".$scoring_matrix{"RE"}."\t".$scoring_matrix{"RG"}."\t".$scoring_matrix{"RH"}."\t".$scoring_matrix{"RI"}."\t".$scoring_matrix{"RL"}."\t".$scoring_matrix{"RK"}."\t".$scoring_matrix{"RM"}."\t".$scoring_matrix{"RF"}."\t".$scoring_matrix{"RP"}."\t".$scoring_matrix{"RS"}."\t".$scoring_matrix{"RT"}."\t".$scoring_matrix{"RW"}."\t".$scoring_matrix{"RY"}."\t".$scoring_matrix{"RV"}."\n";
print "N\t".$scoring_matrix{"NA"}."\t".$scoring_matrix{"NR"}."\t".$scoring_matrix{"NN"}."\t".$scoring_matrix{"ND"}."\t".$scoring_matrix{"NC"}."\t".$scoring_matrix{"NQ"}."\t".$scoring_matrix{"NE"}."\t".$scoring_matrix{"NG"}."\t".$scoring_matrix{"NH"}."\t".$scoring_matrix{"NI"}."\t".$scoring_matrix{"NL"}."\t".$scoring_matrix{"NK"}."\t".$scoring_matrix{"NM"}."\t".$scoring_matrix{"NF"}."\t".$scoring_matrix{"NP"}."\t".$scoring_matrix{"NS"}."\t".$scoring_matrix{"NT"}."\t".$scoring_matrix{"NW"}."\t".$scoring_matrix{"NY"}."\t".$scoring_matrix{"NV"}."\n";
print "D\t".$scoring_matrix{"DA"}."\t".$scoring_matrix{"DR"}."\t".$scoring_matrix{"DN"}."\t".$scoring_matrix{"DD"}."\t".$scoring_matrix{"DC"}."\t".$scoring_matrix{"DQ"}."\t".$scoring_matrix{"DE"}."\t".$scoring_matrix{"DG"}."\t".$scoring_matrix{"DH"}."\t".$scoring_matrix{"DI"}."\t".$scoring_matrix{"DL"}."\t".$scoring_matrix{"DK"}."\t".$scoring_matrix{"DM"}."\t".$scoring_matrix{"DF"}."\t".$scoring_matrix{"DP"}."\t".$scoring_matrix{"DS"}."\t".$scoring_matrix{"DT"}."\t".$scoring_matrix{"DW"}."\t".$scoring_matrix{"DY"}."\t".$scoring_matrix{"DV"}."\n";
print "C\t".$scoring_matrix{"CA"}."\t".$scoring_matrix{"CR"}."\t".$scoring_matrix{"CN"}."\t".$scoring_matrix{"CD"}."\t".$scoring_matrix{"CC"}."\t".$scoring_matrix{"CQ"}."\t".$scoring_matrix{"CE"}."\t".$scoring_matrix{"CG"}."\t".$scoring_matrix{"CH"}."\t".$scoring_matrix{"CI"}."\t".$scoring_matrix{"CL"}."\t".$scoring_matrix{"CK"}."\t".$scoring_matrix{"CM"}."\t".$scoring_matrix{"CF"}."\t".$scoring_matrix{"CP"}."\t".$scoring_matrix{"CS"}."\t".$scoring_matrix{"CT"}."\t".$scoring_matrix{"CW"}."\t".$scoring_matrix{"CY"}."\t".$scoring_matrix{"CV"}."\n";
print "Q\t".$scoring_matrix{"QA"}."\t".$scoring_matrix{"QR"}."\t".$scoring_matrix{"QN"}."\t".$scoring_matrix{"QD"}."\t".$scoring_matrix{"QC"}."\t".$scoring_matrix{"QQ"}."\t".$scoring_matrix{"QE"}."\t".$scoring_matrix{"QG"}."\t".$scoring_matrix{"QH"}."\t".$scoring_matrix{"QI"}."\t".$scoring_matrix{"QL"}."\t".$scoring_matrix{"QK"}."\t".$scoring_matrix{"QM"}."\t".$scoring_matrix{"QF"}."\t".$scoring_matrix{"QP"}."\t".$scoring_matrix{"QS"}."\t".$scoring_matrix{"QT"}."\t".$scoring_matrix{"QW"}."\t".$scoring_matrix{"QY"}."\t".$scoring_matrix{"QV"}."\n";
print "E\t".$scoring_matrix{"EA"}."\t".$scoring_matrix{"ER"}."\t".$scoring_matrix{"EN"}."\t".$scoring_matrix{"ED"}."\t".$scoring_matrix{"EC"}."\t".$scoring_matrix{"EQ"}."\t".$scoring_matrix{"EE"}."\t".$scoring_matrix{"EG"}."\t".$scoring_matrix{"EH"}."\t".$scoring_matrix{"EI"}."\t".$scoring_matrix{"EL"}."\t".$scoring_matrix{"EK"}."\t".$scoring_matrix{"EM"}."\t".$scoring_matrix{"EF"}."\t".$scoring_matrix{"EP"}."\t".$scoring_matrix{"ES"}."\t".$scoring_matrix{"ET"}."\t".$scoring_matrix{"EW"}."\t".$scoring_matrix{"EY"}."\t".$scoring_matrix{"EV"}."\n";
print "G\t".$scoring_matrix{"GA"}."\t".$scoring_matrix{"GR"}."\t".$scoring_matrix{"GN"}."\t".$scoring_matrix{"GD"}."\t".$scoring_matrix{"GC"}."\t".$scoring_matrix{"GQ"}."\t".$scoring_matrix{"GE"}."\t".$scoring_matrix{"GG"}."\t".$scoring_matrix{"GH"}."\t".$scoring_matrix{"GI"}."\t".$scoring_matrix{"GL"}."\t".$scoring_matrix{"GK"}."\t".$scoring_matrix{"GM"}."\t".$scoring_matrix{"GF"}."\t".$scoring_matrix{"GP"}."\t".$scoring_matrix{"GS"}."\t".$scoring_matrix{"GT"}."\t".$scoring_matrix{"GW"}."\t".$scoring_matrix{"GY"}."\t".$scoring_matrix{"GV"}."\n";
print "H\t".$scoring_matrix{"HA"}."\t".$scoring_matrix{"HR"}."\t".$scoring_matrix{"HN"}."\t".$scoring_matrix{"HH"}."\t".$scoring_matrix{"HC"}."\t".$scoring_matrix{"HQ"}."\t".$scoring_matrix{"HE"}."\t".$scoring_matrix{"HG"}."\t".$scoring_matrix{"HH"}."\t".$scoring_matrix{"HI"}."\t".$scoring_matrix{"HL"}."\t".$scoring_matrix{"HK"}."\t".$scoring_matrix{"HM"}."\t".$scoring_matrix{"HF"}."\t".$scoring_matrix{"HP"}."\t".$scoring_matrix{"HS"}."\t".$scoring_matrix{"HT"}."\t".$scoring_matrix{"HW"}."\t".$scoring_matrix{"HY"}."\t".$scoring_matrix{"HV"}."\n";
print "I\t".$scoring_matrix{"IA"}."\t".$scoring_matrix{"IR"}."\t".$scoring_matrix{"IN"}."\t".$scoring_matrix{"ID"}."\t".$scoring_matrix{"IC"}."\t".$scoring_matrix{"IQ"}."\t".$scoring_matrix{"IE"}."\t".$scoring_matrix{"IG"}."\t".$scoring_matrix{"IH"}."\t".$scoring_matrix{"II"}."\t".$scoring_matrix{"IL"}."\t".$scoring_matrix{"IK"}."\t".$scoring_matrix{"IM"}."\t".$scoring_matrix{"IF"}."\t".$scoring_matrix{"IP"}."\t".$scoring_matrix{"IS"}."\t".$scoring_matrix{"IT"}."\t".$scoring_matrix{"IW"}."\t".$scoring_matrix{"IY"}."\t".$scoring_matrix{"IV"}."\n";
print "L\t".$scoring_matrix{"LA"}."\t".$scoring_matrix{"LR"}."\t".$scoring_matrix{"LN"}."\t".$scoring_matrix{"LD"}."\t".$scoring_matrix{"LC"}."\t".$scoring_matrix{"LQ"}."\t".$scoring_matrix{"LE"}."\t".$scoring_matrix{"LG"}."\t".$scoring_matrix{"LH"}."\t".$scoring_matrix{"LI"}."\t".$scoring_matrix{"LL"}."\t".$scoring_matrix{"LK"}."\t".$scoring_matrix{"LM"}."\t".$scoring_matrix{"LF"}."\t".$scoring_matrix{"LP"}."\t".$scoring_matrix{"LS"}."\t".$scoring_matrix{"LT"}."\t".$scoring_matrix{"LW"}."\t".$scoring_matrix{"LY"}."\t".$scoring_matrix{"LV"}."\n";
print "K\t".$scoring_matrix{"KA"}."\t".$scoring_matrix{"KR"}."\t".$scoring_matrix{"KN"}."\t".$scoring_matrix{"KD"}."\t".$scoring_matrix{"KC"}."\t".$scoring_matrix{"KQ"}."\t".$scoring_matrix{"KE"}."\t".$scoring_matrix{"KG"}."\t".$scoring_matrix{"KH"}."\t".$scoring_matrix{"KI"}."\t".$scoring_matrix{"KL"}."\t".$scoring_matrix{"KK"}."\t".$scoring_matrix{"KM"}."\t".$scoring_matrix{"KF"}."\t".$scoring_matrix{"KP"}."\t".$scoring_matrix{"KS"}."\t".$scoring_matrix{"KT"}."\t".$scoring_matrix{"KW"}."\t".$scoring_matrix{"KY"}."\t".$scoring_matrix{"KV"}."\n";
print "M\t".$scoring_matrix{"MA"}."\t".$scoring_matrix{"MR"}."\t".$scoring_matrix{"MN"}."\t".$scoring_matrix{"MD"}."\t".$scoring_matrix{"MC"}."\t".$scoring_matrix{"MQ"}."\t".$scoring_matrix{"ME"}."\t".$scoring_matrix{"MG"}."\t".$scoring_matrix{"MH"}."\t".$scoring_matrix{"MI"}."\t".$scoring_matrix{"ML"}."\t".$scoring_matrix{"MK"}."\t".$scoring_matrix{"MM"}."\t".$scoring_matrix{"MF"}."\t".$scoring_matrix{"MP"}."\t".$scoring_matrix{"MS"}."\t".$scoring_matrix{"MT"}."\t".$scoring_matrix{"MW"}."\t".$scoring_matrix{"MY"}."\t".$scoring_matrix{"MV"}."\n";
print "F\t".$scoring_matrix{"FA"}."\t".$scoring_matrix{"FR"}."\t".$scoring_matrix{"FN"}."\t".$scoring_matrix{"FD"}."\t".$scoring_matrix{"FC"}."\t".$scoring_matrix{"FQ"}."\t".$scoring_matrix{"FE"}."\t".$scoring_matrix{"FG"}."\t".$scoring_matrix{"FH"}."\t".$scoring_matrix{"FI"}."\t".$scoring_matrix{"FL"}."\t".$scoring_matrix{"FK"}."\t".$scoring_matrix{"FM"}."\t".$scoring_matrix{"FF"}."\t".$scoring_matrix{"FP"}."\t".$scoring_matrix{"FS"}."\t".$scoring_matrix{"FT"}."\t".$scoring_matrix{"FW"}."\t".$scoring_matrix{"FY"}."\t".$scoring_matrix{"FV"}."\n";
print "P\t".$scoring_matrix{"PA"}."\t".$scoring_matrix{"PR"}."\t".$scoring_matrix{"PN"}."\t".$scoring_matrix{"PD"}."\t".$scoring_matrix{"PC"}."\t".$scoring_matrix{"PQ"}."\t".$scoring_matrix{"PE"}."\t".$scoring_matrix{"PG"}."\t".$scoring_matrix{"PH"}."\t".$scoring_matrix{"PI"}."\t".$scoring_matrix{"PL"}."\t".$scoring_matrix{"PK"}."\t".$scoring_matrix{"PM"}."\t".$scoring_matrix{"PF"}."\t".$scoring_matrix{"PP"}."\t".$scoring_matrix{"PS"}."\t".$scoring_matrix{"PT"}."\t".$scoring_matrix{"PW"}."\t".$scoring_matrix{"PY"}."\t".$scoring_matrix{"PV"}."\n";
print "S\t".$scoring_matrix{"SA"}."\t".$scoring_matrix{"SR"}."\t".$scoring_matrix{"SN"}."\t".$scoring_matrix{"SD"}."\t".$scoring_matrix{"SC"}."\t".$scoring_matrix{"SQ"}."\t".$scoring_matrix{"SE"}."\t".$scoring_matrix{"SG"}."\t".$scoring_matrix{"SH"}."\t".$scoring_matrix{"SI"}."\t".$scoring_matrix{"SL"}."\t".$scoring_matrix{"SK"}."\t".$scoring_matrix{"SM"}."\t".$scoring_matrix{"SF"}."\t".$scoring_matrix{"SP"}."\t".$scoring_matrix{"SS"}."\t".$scoring_matrix{"ST"}."\t".$scoring_matrix{"SW"}."\t".$scoring_matrix{"SY"}."\t".$scoring_matrix{"SV"}."\n";
print "T\t".$scoring_matrix{"TA"}."\t".$scoring_matrix{"TR"}."\t".$scoring_matrix{"TN"}."\t".$scoring_matrix{"TD"}."\t".$scoring_matrix{"TC"}."\t".$scoring_matrix{"TQ"}."\t".$scoring_matrix{"TE"}."\t".$scoring_matrix{"TG"}."\t".$scoring_matrix{"TH"}."\t".$scoring_matrix{"TI"}."\t".$scoring_matrix{"TL"}."\t".$scoring_matrix{"TK"}."\t".$scoring_matrix{"TM"}."\t".$scoring_matrix{"TF"}."\t".$scoring_matrix{"TP"}."\t".$scoring_matrix{"TS"}."\t".$scoring_matrix{"TT"}."\t".$scoring_matrix{"TW"}."\t".$scoring_matrix{"TY"}."\t".$scoring_matrix{"TV"}."\n";
print "W\t".$scoring_matrix{"WA"}."\t".$scoring_matrix{"WR"}."\t".$scoring_matrix{"WN"}."\t".$scoring_matrix{"WD"}."\t".$scoring_matrix{"WC"}."\t".$scoring_matrix{"WQ"}."\t".$scoring_matrix{"WE"}."\t".$scoring_matrix{"WG"}."\t".$scoring_matrix{"WH"}."\t".$scoring_matrix{"WI"}."\t".$scoring_matrix{"WL"}."\t".$scoring_matrix{"WK"}."\t".$scoring_matrix{"WM"}."\t".$scoring_matrix{"WF"}."\t".$scoring_matrix{"WP"}."\t".$scoring_matrix{"WS"}."\t".$scoring_matrix{"WT"}."\t".$scoring_matrix{"WW"}."\t".$scoring_matrix{"WY"}."\t".$scoring_matrix{"WV"}."\n";
print "Y\t".$scoring_matrix{"YA"}."\t".$scoring_matrix{"YR"}."\t".$scoring_matrix{"YN"}."\t".$scoring_matrix{"YD"}."\t".$scoring_matrix{"YC"}."\t".$scoring_matrix{"YQ"}."\t".$scoring_matrix{"YE"}."\t".$scoring_matrix{"YG"}."\t".$scoring_matrix{"YH"}."\t".$scoring_matrix{"YI"}."\t".$scoring_matrix{"YL"}."\t".$scoring_matrix{"YK"}."\t".$scoring_matrix{"YM"}."\t".$scoring_matrix{"YF"}."\t".$scoring_matrix{"YP"}."\t".$scoring_matrix{"YS"}."\t".$scoring_matrix{"YT"}."\t".$scoring_matrix{"YW"}."\t".$scoring_matrix{"YY"}."\t".$scoring_matrix{"YV"}."\n";
print "V\t".$scoring_matrix{"VA"}."\t".$scoring_matrix{"VR"}."\t".$scoring_matrix{"VN"}."\t".$scoring_matrix{"VD"}."\t".$scoring_matrix{"VC"}."\t".$scoring_matrix{"VQ"}."\t".$scoring_matrix{"VE"}."\t".$scoring_matrix{"VG"}."\t".$scoring_matrix{"VH"}."\t".$scoring_matrix{"VI"}."\t".$scoring_matrix{"VL"}."\t".$scoring_matrix{"VK"}."\t".$scoring_matrix{"VM"}."\t".$scoring_matrix{"VF"}."\t".$scoring_matrix{"VP"}."\t".$scoring_matrix{"VS"}."\t".$scoring_matrix{"VT"}."\t".$scoring_matrix{"VW"}."\t".$scoring_matrix{"VY"}."\t".$scoring_matrix{"VV"}."\n";




#
#Entering Quij Frequecy of occurence values
#
#$quij_matrix{$Hashkey}+=$scoring_matrix{$Hashkey}/$fuij{$_};
#


my %quij_matrix;

foreach(keys %scoring_matrix){
	my $char=substr($_,0,1);
	#my $divide=$fuij{substr($_,0,1)};
	if($totalsubstitutions==0){
		$quij_matrix{$_}=0;
	}else{
		$quij_matrix{$_}=$scoring_matrix{$_}/(($totalsubstitutions-$TotalfrequecyAmbiguity)/2);	#ponder on it divide by zero can be the case
		#$quij_matrix{$_}=$scoring_matrix{$_}/($divideQij-$TotalfrequecyAmbiguity)/2;
	}
	
}


print "\nQij Matrix\n";


print "\n \tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\n";
print "A\t".$quij_matrix{"AA"}."\t".$quij_matrix{"AR"}."\t".$quij_matrix{"AN"}."\t".$quij_matrix{"AD"}."\t".$quij_matrix{"AC"}."\t".$quij_matrix{"AQ"}."\t".$quij_matrix{"AE"}."\t".$quij_matrix{"AG"}."\t".$quij_matrix{"AH"}."\t".$quij_matrix{"AI"}."\t".$quij_matrix{"AL"}."\t".$quij_matrix{"AK"}."\t".$quij_matrix{"AM"}."\t".$quij_matrix{"AF"}."\t".$quij_matrix{"AP"}."\t".$quij_matrix{"AS"}."\t".$quij_matrix{"AT"}."\t".$quij_matrix{"AW"}."\t".$quij_matrix{"AY"}."\t".$quij_matrix{"AV"}."\n";
print "R\t".$quij_matrix{"RA"}."\t".$quij_matrix{"RR"}."\t".$quij_matrix{"RN"}."\t".$quij_matrix{"RD"}."\t".$quij_matrix{"RC"}."\t".$quij_matrix{"RQ"}."\t".$quij_matrix{"RE"}."\t".$quij_matrix{"RG"}."\t".$quij_matrix{"RH"}."\t".$quij_matrix{"RI"}."\t".$quij_matrix{"RL"}."\t".$quij_matrix{"RK"}."\t".$quij_matrix{"RM"}."\t".$quij_matrix{"RF"}."\t".$quij_matrix{"RP"}."\t".$quij_matrix{"RS"}."\t".$quij_matrix{"RT"}."\t".$quij_matrix{"RW"}."\t".$quij_matrix{"RY"}."\t".$quij_matrix{"RV"}."\n";
print "N\t".$quij_matrix{"NA"}."\t".$quij_matrix{"NR"}."\t".$quij_matrix{"NN"}."\t".$quij_matrix{"ND"}."\t".$quij_matrix{"NC"}."\t".$quij_matrix{"NQ"}."\t".$quij_matrix{"NE"}."\t".$quij_matrix{"NG"}."\t".$quij_matrix{"NH"}."\t".$quij_matrix{"NI"}."\t".$quij_matrix{"NL"}."\t".$quij_matrix{"NK"}."\t".$quij_matrix{"NM"}."\t".$quij_matrix{"NF"}."\t".$quij_matrix{"NP"}."\t".$quij_matrix{"NS"}."\t".$quij_matrix{"NT"}."\t".$quij_matrix{"NW"}."\t".$quij_matrix{"NY"}."\t".$quij_matrix{"NV"}."\n";
print "D\t".$quij_matrix{"DA"}."\t".$quij_matrix{"DR"}."\t".$quij_matrix{"DN"}."\t".$quij_matrix{"DD"}."\t".$quij_matrix{"DC"}."\t".$quij_matrix{"DQ"}."\t".$quij_matrix{"DE"}."\t".$quij_matrix{"DG"}."\t".$quij_matrix{"DH"}."\t".$quij_matrix{"DI"}."\t".$quij_matrix{"DL"}."\t".$quij_matrix{"DK"}."\t".$quij_matrix{"DM"}."\t".$quij_matrix{"DF"}."\t".$quij_matrix{"DP"}."\t".$quij_matrix{"DS"}."\t".$quij_matrix{"DT"}."\t".$quij_matrix{"DW"}."\t".$quij_matrix{"DY"}."\t".$quij_matrix{"DV"}."\n";
print "C\t".$quij_matrix{"CA"}."\t".$quij_matrix{"CR"}."\t".$quij_matrix{"CN"}."\t".$quij_matrix{"CD"}."\t".$quij_matrix{"CC"}."\t".$quij_matrix{"CQ"}."\t".$quij_matrix{"CE"}."\t".$quij_matrix{"CG"}."\t".$quij_matrix{"CH"}."\t".$quij_matrix{"CI"}."\t".$quij_matrix{"CL"}."\t".$quij_matrix{"CK"}."\t".$quij_matrix{"CM"}."\t".$quij_matrix{"CF"}."\t".$quij_matrix{"CP"}."\t".$quij_matrix{"CS"}."\t".$quij_matrix{"CT"}."\t".$quij_matrix{"CW"}."\t".$quij_matrix{"CY"}."\t".$quij_matrix{"CV"}."\n";
print "Q\t".$quij_matrix{"QA"}."\t".$quij_matrix{"QR"}."\t".$quij_matrix{"QN"}."\t".$quij_matrix{"QD"}."\t".$quij_matrix{"QC"}."\t".$quij_matrix{"QQ"}."\t".$quij_matrix{"QE"}."\t".$quij_matrix{"QG"}."\t".$quij_matrix{"QH"}."\t".$quij_matrix{"QI"}."\t".$quij_matrix{"QL"}."\t".$quij_matrix{"QK"}."\t".$quij_matrix{"QM"}."\t".$quij_matrix{"QF"}."\t".$quij_matrix{"QP"}."\t".$quij_matrix{"QS"}."\t".$quij_matrix{"QT"}."\t".$quij_matrix{"QW"}."\t".$quij_matrix{"QY"}."\t".$quij_matrix{"QV"}."\n";
print "E\t".$quij_matrix{"EA"}."\t".$quij_matrix{"ER"}."\t".$quij_matrix{"EN"}."\t".$quij_matrix{"ED"}."\t".$quij_matrix{"EC"}."\t".$quij_matrix{"EQ"}."\t".$quij_matrix{"EE"}."\t".$quij_matrix{"EG"}."\t".$quij_matrix{"EH"}."\t".$quij_matrix{"EI"}."\t".$quij_matrix{"EL"}."\t".$quij_matrix{"EK"}."\t".$quij_matrix{"EM"}."\t".$quij_matrix{"EF"}."\t".$quij_matrix{"EP"}."\t".$quij_matrix{"ES"}."\t".$quij_matrix{"ET"}."\t".$quij_matrix{"EW"}."\t".$quij_matrix{"EY"}."\t".$quij_matrix{"EV"}."\n";
print "G\t".$quij_matrix{"GA"}."\t".$quij_matrix{"GR"}."\t".$quij_matrix{"GN"}."\t".$quij_matrix{"GD"}."\t".$quij_matrix{"GC"}."\t".$quij_matrix{"GQ"}."\t".$quij_matrix{"GE"}."\t".$quij_matrix{"GG"}."\t".$quij_matrix{"GH"}."\t".$quij_matrix{"GI"}."\t".$quij_matrix{"GL"}."\t".$quij_matrix{"GK"}."\t".$quij_matrix{"GM"}."\t".$quij_matrix{"GF"}."\t".$quij_matrix{"GP"}."\t".$quij_matrix{"GS"}."\t".$quij_matrix{"GT"}."\t".$quij_matrix{"GW"}."\t".$quij_matrix{"GY"}."\t".$quij_matrix{"GV"}."\n";
print "H\t".$quij_matrix{"HA"}."\t".$quij_matrix{"HR"}."\t".$quij_matrix{"HN"}."\t".$quij_matrix{"HH"}."\t".$quij_matrix{"HC"}."\t".$quij_matrix{"HQ"}."\t".$quij_matrix{"HE"}."\t".$quij_matrix{"HG"}."\t".$quij_matrix{"HH"}."\t".$quij_matrix{"HI"}."\t".$quij_matrix{"HL"}."\t".$quij_matrix{"HK"}."\t".$quij_matrix{"HM"}."\t".$quij_matrix{"HF"}."\t".$quij_matrix{"HP"}."\t".$quij_matrix{"HS"}."\t".$quij_matrix{"HT"}."\t".$quij_matrix{"HW"}."\t".$quij_matrix{"HY"}."\t".$quij_matrix{"HV"}."\n";
print "I\t".$quij_matrix{"IA"}."\t".$quij_matrix{"IR"}."\t".$quij_matrix{"IN"}."\t".$quij_matrix{"ID"}."\t".$quij_matrix{"IC"}."\t".$quij_matrix{"IQ"}."\t".$quij_matrix{"IE"}."\t".$quij_matrix{"IG"}."\t".$quij_matrix{"IH"}."\t".$quij_matrix{"II"}."\t".$quij_matrix{"IL"}."\t".$quij_matrix{"IK"}."\t".$quij_matrix{"IM"}."\t".$quij_matrix{"IF"}."\t".$quij_matrix{"IP"}."\t".$quij_matrix{"IS"}."\t".$quij_matrix{"IT"}."\t".$quij_matrix{"IW"}."\t".$quij_matrix{"IY"}."\t".$quij_matrix{"IV"}."\n";
print "L\t".$quij_matrix{"LA"}."\t".$quij_matrix{"LR"}."\t".$quij_matrix{"LN"}."\t".$quij_matrix{"LD"}."\t".$quij_matrix{"LC"}."\t".$quij_matrix{"LQ"}."\t".$quij_matrix{"LE"}."\t".$quij_matrix{"LG"}."\t".$quij_matrix{"LH"}."\t".$quij_matrix{"LI"}."\t".$quij_matrix{"LL"}."\t".$quij_matrix{"LK"}."\t".$quij_matrix{"LM"}."\t".$quij_matrix{"LF"}."\t".$quij_matrix{"LP"}."\t".$quij_matrix{"LS"}."\t".$quij_matrix{"LT"}."\t".$quij_matrix{"LW"}."\t".$quij_matrix{"LY"}."\t".$quij_matrix{"LV"}."\n";
print "K\t".$quij_matrix{"KA"}."\t".$quij_matrix{"KR"}."\t".$quij_matrix{"KN"}."\t".$quij_matrix{"KD"}."\t".$quij_matrix{"KC"}."\t".$quij_matrix{"KQ"}."\t".$quij_matrix{"KE"}."\t".$quij_matrix{"KG"}."\t".$quij_matrix{"KH"}."\t".$quij_matrix{"KI"}."\t".$quij_matrix{"KL"}."\t".$quij_matrix{"KK"}."\t".$quij_matrix{"KM"}."\t".$quij_matrix{"KF"}."\t".$quij_matrix{"KP"}."\t".$quij_matrix{"KS"}."\t".$quij_matrix{"KT"}."\t".$quij_matrix{"KW"}."\t".$quij_matrix{"KY"}."\t".$quij_matrix{"KV"}."\n";
print "M\t".$quij_matrix{"MA"}."\t".$quij_matrix{"MR"}."\t".$quij_matrix{"MN"}."\t".$quij_matrix{"MD"}."\t".$quij_matrix{"MC"}."\t".$quij_matrix{"MQ"}."\t".$quij_matrix{"ME"}."\t".$quij_matrix{"MG"}."\t".$quij_matrix{"MH"}."\t".$quij_matrix{"MI"}."\t".$quij_matrix{"ML"}."\t".$quij_matrix{"MK"}."\t".$quij_matrix{"MM"}."\t".$quij_matrix{"MF"}."\t".$quij_matrix{"MP"}."\t".$quij_matrix{"MS"}."\t".$quij_matrix{"MT"}."\t".$quij_matrix{"MW"}."\t".$quij_matrix{"MY"}."\t".$quij_matrix{"MV"}."\n";
print "F\t".$quij_matrix{"FA"}."\t".$quij_matrix{"FR"}."\t".$quij_matrix{"FN"}."\t".$quij_matrix{"FD"}."\t".$quij_matrix{"FC"}."\t".$quij_matrix{"FQ"}."\t".$quij_matrix{"FE"}."\t".$quij_matrix{"FG"}."\t".$quij_matrix{"FH"}."\t".$quij_matrix{"FI"}."\t".$quij_matrix{"FL"}."\t".$quij_matrix{"FK"}."\t".$quij_matrix{"FM"}."\t".$quij_matrix{"FF"}."\t".$quij_matrix{"FP"}."\t".$quij_matrix{"FS"}."\t".$quij_matrix{"FT"}."\t".$quij_matrix{"FW"}."\t".$quij_matrix{"FY"}."\t".$quij_matrix{"FV"}."\n";
print "P\t".$quij_matrix{"PA"}."\t".$quij_matrix{"PR"}."\t".$quij_matrix{"PN"}."\t".$quij_matrix{"PD"}."\t".$quij_matrix{"PC"}."\t".$quij_matrix{"PQ"}."\t".$quij_matrix{"PE"}."\t".$quij_matrix{"PG"}."\t".$quij_matrix{"PH"}."\t".$quij_matrix{"PI"}."\t".$quij_matrix{"PL"}."\t".$quij_matrix{"PK"}."\t".$quij_matrix{"PM"}."\t".$quij_matrix{"PF"}."\t".$quij_matrix{"PP"}."\t".$quij_matrix{"PS"}."\t".$quij_matrix{"PT"}."\t".$quij_matrix{"PW"}."\t".$quij_matrix{"PY"}."\t".$quij_matrix{"PV"}."\n";
print "S\t".$quij_matrix{"SA"}."\t".$quij_matrix{"SR"}."\t".$quij_matrix{"SN"}."\t".$quij_matrix{"SD"}."\t".$quij_matrix{"SC"}."\t".$quij_matrix{"SQ"}."\t".$quij_matrix{"SE"}."\t".$quij_matrix{"SG"}."\t".$quij_matrix{"SH"}."\t".$quij_matrix{"SI"}."\t".$quij_matrix{"SL"}."\t".$quij_matrix{"SK"}."\t".$quij_matrix{"SM"}."\t".$quij_matrix{"SF"}."\t".$quij_matrix{"SP"}."\t".$quij_matrix{"SS"}."\t".$quij_matrix{"ST"}."\t".$quij_matrix{"SW"}."\t".$quij_matrix{"SY"}."\t".$quij_matrix{"SV"}."\n";
print "T\t".$quij_matrix{"TA"}."\t".$quij_matrix{"TR"}."\t".$quij_matrix{"TN"}."\t".$quij_matrix{"TD"}."\t".$quij_matrix{"TC"}."\t".$quij_matrix{"TQ"}."\t".$quij_matrix{"TE"}."\t".$quij_matrix{"TG"}."\t".$quij_matrix{"TH"}."\t".$quij_matrix{"TI"}."\t".$quij_matrix{"TL"}."\t".$quij_matrix{"TK"}."\t".$quij_matrix{"TM"}."\t".$quij_matrix{"TF"}."\t".$quij_matrix{"TP"}."\t".$quij_matrix{"TS"}."\t".$quij_matrix{"TT"}."\t".$quij_matrix{"TW"}."\t".$quij_matrix{"TY"}."\t".$quij_matrix{"TV"}."\n";
print "W\t".$quij_matrix{"WA"}."\t".$quij_matrix{"WR"}."\t".$quij_matrix{"WN"}."\t".$quij_matrix{"WD"}."\t".$quij_matrix{"WC"}."\t".$quij_matrix{"WQ"}."\t".$quij_matrix{"WE"}."\t".$quij_matrix{"WG"}."\t".$quij_matrix{"WH"}."\t".$quij_matrix{"WI"}."\t".$quij_matrix{"WL"}."\t".$quij_matrix{"WK"}."\t".$quij_matrix{"WM"}."\t".$quij_matrix{"WF"}."\t".$quij_matrix{"WP"}."\t".$quij_matrix{"WS"}."\t".$quij_matrix{"WT"}."\t".$quij_matrix{"WW"}."\t".$quij_matrix{"WY"}."\t".$quij_matrix{"WV"}."\n";
print "Y\t".$quij_matrix{"YA"}."\t".$quij_matrix{"YR"}."\t".$quij_matrix{"YN"}."\t".$quij_matrix{"YD"}."\t".$quij_matrix{"YC"}."\t".$quij_matrix{"YQ"}."\t".$quij_matrix{"YE"}."\t".$quij_matrix{"YG"}."\t".$quij_matrix{"YH"}."\t".$quij_matrix{"YI"}."\t".$quij_matrix{"YL"}."\t".$quij_matrix{"YK"}."\t".$quij_matrix{"YM"}."\t".$quij_matrix{"YF"}."\t".$quij_matrix{"YP"}."\t".$quij_matrix{"YS"}."\t".$quij_matrix{"YT"}."\t".$quij_matrix{"YW"}."\t".$quij_matrix{"YY"}."\t".$quij_matrix{"YV"}."\n";
print "V\t".$quij_matrix{"VA"}."\t".$quij_matrix{"VR"}."\t".$quij_matrix{"VN"}."\t".$quij_matrix{"VD"}."\t".$quij_matrix{"VC"}."\t".$quij_matrix{"VQ"}."\t".$quij_matrix{"VE"}."\t".$quij_matrix{"VG"}."\t".$quij_matrix{"VH"}."\t".$quij_matrix{"VI"}."\t".$quij_matrix{"VL"}."\t".$quij_matrix{"VK"}."\t".$quij_matrix{"VM"}."\t".$quij_matrix{"VF"}."\t".$quij_matrix{"VP"}."\t".$quij_matrix{"VS"}."\t".$quij_matrix{"VT"}."\t".$quij_matrix{"VW"}."\t".$quij_matrix{"VY"}."\t".$quij_matrix{"VV"}."\n";

#handle Pij Probability of occurrence for expected probability


my %Pij;

#Summation of All Qij's

foreach(keys %quij_matrix){
	my $firstsubstitute=substr($_,0,1);
	my $secondsubstitute=substr($_,1,1);
	
	if($firstsubstitute eq $secondsubstitute){
		
	}else{
		$TotalQijExcSimilar+=$quij_matrix{$_};
	}

}

print "\nTotal Qij Except i==i and j==j : ".$TotalQijExcSimilar."\n\n";


foreach(keys %quij_matrix){
	my $firstchar2=substr($_,0,1);
	my $secondchar2=substr($_,1,1);

		if($TotalQijExcSimilar==0){
			#$Pij{$firstchar2}=$quij_matrix{$firstchar2.$secondchar2};
		}else{
			$Pij{$firstchar2}=$quij_matrix{$firstchar2.$firstchar2}+$quij_matrix{$firstchar2.$secondchar2}/($TotalQijExcSimilar/2);
		}
}


#
#Pij Ends here
#


#Eij Calculations

my %Eij;

foreach(keys %quij_matrix){
	my $firstchar3=substr($_,0,1);
	my $secondchar3=substr($_,1,1);
	
	if($firstchar3 eq $secondchar3){
		$Eij{$_}=$Pij{$firstchar3}*$Pij{$secondchar3};
	}else{
		$Eij{$_}=2*($Pij{$firstchar3}+$Pij{$secondchar3});
	}
}

print $Eij{"AR"}." ... ".$Eij{"RA"}."\n";

#
#Eij Ends here
#


###Sij implementation



my %Sij;

foreach(keys %Eij){
	if($Eij{$_}==0 || $quij_matrix{$_}==0){		#handles divide by zero exception
		$Sij{$_}=0;
	}else{
		$Sij{$_}=int(log(($quij_matrix{$_}/$Eij{$_})/log(2)));
	}
}


print "\nSij Matrix\n";


print "\n \tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\n";
print "A\t".$Sij{"AA"}."\t".$Sij{"AR"}."\t".$Sij{"AN"}."\t".$Sij{"AD"}."\t".$Sij{"AC"}."\t".$Sij{"AQ"}."\t".$Sij{"AE"}."\t".$Sij{"AG"}."\t".$Sij{"AH"}."\t".$Sij{"AI"}."\t".$Sij{"AL"}."\t".$Sij{"AK"}."\t".$Sij{"AM"}."\t".$Sij{"AF"}."\t".$Sij{"AP"}."\t".$Sij{"AS"}."\t".$Sij{"AT"}."\t".$Sij{"AW"}."\t".$Sij{"AY"}."\t".$Sij{"AV"}."\n";
print "R\t".$Sij{"RA"}."\t".$Sij{"RR"}."\t".$Sij{"RN"}."\t".$Sij{"RD"}."\t".$Sij{"RC"}."\t".$Sij{"RQ"}."\t".$Sij{"RE"}."\t".$Sij{"RG"}."\t".$Sij{"RH"}."\t".$Sij{"RI"}."\t".$Sij{"RL"}."\t".$Sij{"RK"}."\t".$Sij{"RM"}."\t".$Sij{"RF"}."\t".$Sij{"RP"}."\t".$Sij{"RS"}."\t".$Sij{"RT"}."\t".$Sij{"RW"}."\t".$Sij{"RY"}."\t".$Sij{"RV"}."\n";
print "N\t".$Sij{"NA"}."\t".$Sij{"NR"}."\t".$Sij{"NN"}."\t".$Sij{"ND"}."\t".$Sij{"NC"}."\t".$Sij{"NQ"}."\t".$Sij{"NE"}."\t".$Sij{"NG"}."\t".$Sij{"NH"}."\t".$Sij{"NI"}."\t".$Sij{"NL"}."\t".$Sij{"NK"}."\t".$Sij{"NM"}."\t".$Sij{"NF"}."\t".$Sij{"NP"}."\t".$Sij{"NS"}."\t".$Sij{"NT"}."\t".$Sij{"NW"}."\t".$Sij{"NY"}."\t".$Sij{"NV"}."\n";
print "D\t".$Sij{"DA"}."\t".$Sij{"DR"}."\t".$Sij{"DN"}."\t".$Sij{"DD"}."\t".$Sij{"DC"}."\t".$Sij{"DQ"}."\t".$Sij{"DE"}."\t".$Sij{"DG"}."\t".$Sij{"DH"}."\t".$Sij{"DI"}."\t".$Sij{"DL"}."\t".$Sij{"DK"}."\t".$Sij{"DM"}."\t".$Sij{"DF"}."\t".$Sij{"DP"}."\t".$Sij{"DS"}."\t".$Sij{"DT"}."\t".$Sij{"DW"}."\t".$Sij{"DY"}."\t".$Sij{"DV"}."\n";
print "C\t".$Sij{"CA"}."\t".$Sij{"CR"}."\t".$Sij{"CN"}."\t".$Sij{"CD"}."\t".$Sij{"CC"}."\t".$Sij{"CQ"}."\t".$Sij{"CE"}."\t".$Sij{"CG"}."\t".$Sij{"CH"}."\t".$Sij{"CI"}."\t".$Sij{"CL"}."\t".$Sij{"CK"}."\t".$Sij{"CM"}."\t".$Sij{"CF"}."\t".$Sij{"CP"}."\t".$Sij{"CS"}."\t".$Sij{"CT"}."\t".$Sij{"CW"}."\t".$Sij{"CY"}."\t".$Sij{"CV"}."\n";
print "Q\t".$Sij{"QA"}."\t".$Sij{"QR"}."\t".$Sij{"QN"}."\t".$Sij{"QD"}."\t".$Sij{"QC"}."\t".$Sij{"QQ"}."\t".$Sij{"QE"}."\t".$Sij{"QG"}."\t".$Sij{"QH"}."\t".$Sij{"QI"}."\t".$Sij{"QL"}."\t".$Sij{"QK"}."\t".$Sij{"QM"}."\t".$Sij{"QF"}."\t".$Sij{"QP"}."\t".$Sij{"QS"}."\t".$Sij{"QT"}."\t".$Sij{"QW"}."\t".$Sij{"QY"}."\t".$Sij{"QV"}."\n";
print "E\t".$Sij{"EA"}."\t".$Sij{"ER"}."\t".$Sij{"EN"}."\t".$Sij{"ED"}."\t".$Sij{"EC"}."\t".$Sij{"EQ"}."\t".$Sij{"EE"}."\t".$Sij{"EG"}."\t".$Sij{"EH"}."\t".$Sij{"EI"}."\t".$Sij{"EL"}."\t".$Sij{"EK"}."\t".$Sij{"EM"}."\t".$Sij{"EF"}."\t".$Sij{"EP"}."\t".$Sij{"ES"}."\t".$Sij{"ET"}."\t".$Sij{"EW"}."\t".$Sij{"EY"}."\t".$Sij{"EV"}."\n";
print "G\t".$Sij{"GA"}."\t".$Sij{"GR"}."\t".$Sij{"GN"}."\t".$Sij{"GD"}."\t".$Sij{"GC"}."\t".$Sij{"GQ"}."\t".$Sij{"GE"}."\t".$Sij{"GG"}."\t".$Sij{"GH"}."\t".$Sij{"GI"}."\t".$Sij{"GL"}."\t".$Sij{"GK"}."\t".$Sij{"GM"}."\t".$Sij{"GF"}."\t".$Sij{"GP"}."\t".$Sij{"GS"}."\t".$Sij{"GT"}."\t".$Sij{"GW"}."\t".$Sij{"GY"}."\t".$Sij{"GV"}."\n";
print "H\t".$Sij{"HA"}."\t".$Sij{"HR"}."\t".$Sij{"HN"}."\t".$Sij{"HH"}."\t".$Sij{"HC"}."\t".$Sij{"HQ"}."\t".$Sij{"HE"}."\t".$Sij{"HG"}."\t".$Sij{"HH"}."\t".$Sij{"HI"}."\t".$Sij{"HL"}."\t".$Sij{"HK"}."\t".$Sij{"HM"}."\t".$Sij{"HF"}."\t".$Sij{"HP"}."\t".$Sij{"HS"}."\t".$Sij{"HT"}."\t".$Sij{"HW"}."\t".$Sij{"HY"}."\t".$Sij{"HV"}."\n";
print "I\t".$Sij{"IA"}."\t".$Sij{"IR"}."\t".$Sij{"IN"}."\t".$Sij{"ID"}."\t".$Sij{"IC"}."\t".$Sij{"IQ"}."\t".$Sij{"IE"}."\t".$Sij{"IG"}."\t".$Sij{"IH"}."\t".$Sij{"II"}."\t".$Sij{"IL"}."\t".$Sij{"IK"}."\t".$Sij{"IM"}."\t".$Sij{"IF"}."\t".$Sij{"IP"}."\t".$Sij{"IS"}."\t".$Sij{"IT"}."\t".$Sij{"IW"}."\t".$Sij{"IY"}."\t".$Sij{"IV"}."\n";
print "L\t".$Sij{"LA"}."\t".$Sij{"LR"}."\t".$Sij{"LN"}."\t".$Sij{"LD"}."\t".$Sij{"LC"}."\t".$Sij{"LQ"}."\t".$Sij{"LE"}."\t".$Sij{"LG"}."\t".$Sij{"LH"}."\t".$Sij{"LI"}."\t".$Sij{"LL"}."\t".$Sij{"LK"}."\t".$Sij{"LM"}."\t".$Sij{"LF"}."\t".$Sij{"LP"}."\t".$Sij{"LS"}."\t".$Sij{"LT"}."\t".$Sij{"LW"}."\t".$Sij{"LY"}."\t".$Sij{"LV"}."\n";
print "K\t".$Sij{"KA"}."\t".$Sij{"KR"}."\t".$Sij{"KN"}."\t".$Sij{"KD"}."\t".$Sij{"KC"}."\t".$Sij{"KQ"}."\t".$Sij{"KE"}."\t".$Sij{"KG"}."\t".$Sij{"KH"}."\t".$Sij{"KI"}."\t".$Sij{"KL"}."\t".$Sij{"KK"}."\t".$Sij{"KM"}."\t".$Sij{"KF"}."\t".$Sij{"KP"}."\t".$Sij{"KS"}."\t".$Sij{"KT"}."\t".$Sij{"KW"}."\t".$Sij{"KY"}."\t".$Sij{"KV"}."\n";
print "M\t".$Sij{"MA"}."\t".$Sij{"MR"}."\t".$Sij{"MN"}."\t".$Sij{"MD"}."\t".$Sij{"MC"}."\t".$Sij{"MQ"}."\t".$Sij{"ME"}."\t".$Sij{"MG"}."\t".$Sij{"MH"}."\t".$Sij{"MI"}."\t".$Sij{"ML"}."\t".$Sij{"MK"}."\t".$Sij{"MM"}."\t".$Sij{"MF"}."\t".$Sij{"MP"}."\t".$Sij{"MS"}."\t".$Sij{"MT"}."\t".$Sij{"MW"}."\t".$Sij{"MY"}."\t".$Sij{"MV"}."\n";
print "F\t".$Sij{"FA"}."\t".$Sij{"FR"}."\t".$Sij{"FN"}."\t".$Sij{"FD"}."\t".$Sij{"FC"}."\t".$Sij{"FQ"}."\t".$Sij{"FE"}."\t".$Sij{"FG"}."\t".$Sij{"FH"}."\t".$Sij{"FI"}."\t".$Sij{"FL"}."\t".$Sij{"FK"}."\t".$Sij{"FM"}."\t".$Sij{"FF"}."\t".$Sij{"FP"}."\t".$Sij{"FS"}."\t".$Sij{"FT"}."\t".$Sij{"FW"}."\t".$Sij{"FY"}."\t".$Sij{"FV"}."\n";
print "P\t".$Sij{"PA"}."\t".$Sij{"PR"}."\t".$Sij{"PN"}."\t".$Sij{"PD"}."\t".$Sij{"PC"}."\t".$Sij{"PQ"}."\t".$Sij{"PE"}."\t".$Sij{"PG"}."\t".$Sij{"PH"}."\t".$Sij{"PI"}."\t".$Sij{"PL"}."\t".$Sij{"PK"}."\t".$Sij{"PM"}."\t".$Sij{"PF"}."\t".$Sij{"PP"}."\t".$Sij{"PS"}."\t".$Sij{"PT"}."\t".$Sij{"PW"}."\t".$Sij{"PY"}."\t".$Sij{"PV"}."\n";
print "S\t".$Sij{"SA"}."\t".$Sij{"SR"}."\t".$Sij{"SN"}."\t".$Sij{"SD"}."\t".$Sij{"SC"}."\t".$Sij{"SQ"}."\t".$Sij{"SE"}."\t".$Sij{"SG"}."\t".$Sij{"SH"}."\t".$Sij{"SI"}."\t".$Sij{"SL"}."\t".$Sij{"SK"}."\t".$Sij{"SM"}."\t".$Sij{"SF"}."\t".$Sij{"SP"}."\t".$Sij{"SS"}."\t".$Sij{"ST"}."\t".$Sij{"SW"}."\t".$Sij{"SY"}."\t".$Sij{"SV"}."\n";
print "T\t".$Sij{"TA"}."\t".$Sij{"TR"}."\t".$Sij{"TN"}."\t".$Sij{"TD"}."\t".$Sij{"TC"}."\t".$Sij{"TQ"}."\t".$Sij{"TE"}."\t".$Sij{"TG"}."\t".$Sij{"TH"}."\t".$Sij{"TI"}."\t".$Sij{"TL"}."\t".$Sij{"TK"}."\t".$Sij{"TM"}."\t".$Sij{"TF"}."\t".$Sij{"TP"}."\t".$Sij{"TS"}."\t".$Sij{"TT"}."\t".$Sij{"TW"}."\t".$Sij{"TY"}."\t".$Sij{"TV"}."\n";
print "W\t".$Sij{"WA"}."\t".$Sij{"WR"}."\t".$Sij{"WN"}."\t".$Sij{"WD"}."\t".$Sij{"WC"}."\t".$Sij{"WQ"}."\t".$Sij{"WE"}."\t".$Sij{"WG"}."\t".$Sij{"WH"}."\t".$Sij{"WI"}."\t".$Sij{"WL"}."\t".$Sij{"WK"}."\t".$Sij{"WM"}."\t".$Sij{"WF"}."\t".$Sij{"WP"}."\t".$Sij{"WS"}."\t".$Sij{"WT"}."\t".$Sij{"WW"}."\t".$Sij{"WY"}."\t".$Sij{"WV"}."\n";
print "Y\t".$Sij{"YA"}."\t".$Sij{"YR"}."\t".$Sij{"YN"}."\t".$Sij{"YD"}."\t".$Sij{"YC"}."\t".$Sij{"YQ"}."\t".$Sij{"YE"}."\t".$Sij{"YG"}."\t".$Sij{"YH"}."\t".$Sij{"YI"}."\t".$Sij{"YL"}."\t".$Sij{"YK"}."\t".$Sij{"YM"}."\t".$Sij{"YF"}."\t".$Sij{"YP"}."\t".$Sij{"YS"}."\t".$Sij{"YT"}."\t".$Sij{"YW"}."\t".$Sij{"YY"}."\t".$Sij{"YV"}."\n";
print "V\t".$Sij{"VA"}."\t".$Sij{"VR"}."\t".$Sij{"VN"}."\t".$Sij{"VD"}."\t".$Sij{"VC"}."\t".$Sij{"VQ"}."\t".$Sij{"VE"}."\t".$Sij{"VG"}."\t".$Sij{"VH"}."\t".$Sij{"VI"}."\t".$Sij{"VL"}."\t".$Sij{"VK"}."\t".$Sij{"VM"}."\t".$Sij{"VF"}."\t".$Sij{"VP"}."\t".$Sij{"VS"}."\t".$Sij{"VT"}."\t".$Sij{"VW"}."\t".$Sij{"VY"}."\t".$Sij{"VV"}."\n";


my $Entropy=0;

foreach(keys %quij_matrix){
	$Entropy+=$quij_matrix{$_}*$Sij{$_};
}

print "@ Entropy :  ".$Entropy."\n";




=head

my $Expect=0;
foreach(keys %quij_matrix){
	my $firstchar4=substr($_,0,1);
	my $secondchar4=substr($_,1,1);
	$Expect+=$Pij{$firstchar4}*$Pij{$secondchar4}*$Sij{$_};
}

print "@ ExpectValue :  ".$Expect."\n";

Error Entry points 235 300

=cut


print "END\n\n";

