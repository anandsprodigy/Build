#!c:\perl\bin\perl.exe -w -s
use strict;
use warnings;
use BlockSeqGetter;
use BLASTmake;
require MakeBlastDB;
use DownSpeciesFasta;
use TaxParentGet;
use FastaToBlockGet;

BEGIN{
	push(@INC,"C:\\Users\\Anand\\Desktop\\Build");
}

FastaToBlockGet::FastaToBlockGet();