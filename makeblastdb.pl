#!c:\perl\bin\perl.exe -w
use strict;
use warnings;
use MakeBlastDB;
use TaxParentGet;

MakeBlastDB::MakeBlastDB(TaxParentGet::GetParentTaxID(33208));