#!c:\perl\bin\perl.exe -w
use strict;
use BLASTmake;
require MakeBlastDB;

our $DBN;

BLASTmake::BLASTmake("DBFASTA","SpeciesFASTA.txt");

