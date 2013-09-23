#!c:\perl\bin\perl.exe -w -s
use strict;
use warnings;
use BlockSeqGetter;
use BLASTmake;
require MakeBlastDB;
use DownSpeciesFasta;
use TaxParentGet;
use FastaToBlockGet;
use BlockParser;
use WWW::Mechanize;
use fastatoblocks_proxy;

BEGIN{
	push(@INC,"c:\\users\\vscripts\\desktop");
	
}


open(BLOCK,">C:\\Users\\Anand\\Desktop\\Build\\block.mblk") or die($!);

#Download The Given Species Swiss-Prot sequences in Fasta format

#DownSpeciesFasta::DownSpeciesFasta("Entamoeba","histolytica");

#Download the Swiss-prot sequences in Family level

#MakeBlastDB::MakeBlastDB(5759);

#Run the Blast Query given org and family database

#BLASTmake::BLASTmake("DBFASTA","SpeciesFASTA.txt");

#Retrive the blast sequence from blast output

#BlockSeqGetter::BlockSeqGetter("BlockInput.txt");

#Submit BLOCKs to server

#FastaToBlockGet::FastaToBlockGet('anandsprodigy@gmail.com',"http://127.0.0.0:8080/","12.txt","blockprepared.txt");

#BLOCKS file cleaning 

#BlockParser::BlockParser("BlockOut.blk");

#Calculating BLOSUM matrix

#system("perl SijScore.pl");

my $multiplefasta=">1\\nARGGSRRDDWRAEEERRAADDARRGGSRRARRGGGGSSRRGGGARRSGGGS\\n>2\\nARARGGRDDWRAEEERRAADDARRGGSRRARRGGGGSSRRGGGARRSGGGS";


my @blocks = blocksretrieve($multiplefasta,"2");


foreach(@blocks){
	print BLOCK $_;
}



close(BLOCK);


system("perl C:\\Users\\Anand\\Desktop\\Build\\SijScore.pl");