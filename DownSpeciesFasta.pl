#!c:\perl\bin\perl.exe -w
use strict;
use warnings;
use DownSpeciesFasta;

BEGIN{
	push(@INC,"C:\\Users\\Anand\\Desktop\\Build");
}

DownSpeciesFasta::DownSpeciesFasta("Entamoeba","histolytica");