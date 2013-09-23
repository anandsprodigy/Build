package MakeBlastDB;

use WWW::Mechanize;
use Cwd;


BEGIN{
    my $CWDir=Cwd::_win32_cwd();
    push(@INC,$CWDir);
}

sub MakeBlastDB{

my ($TaxonomyID)=@_;

my $CWDir=Cwd::_win32_cwd();

open(AD,">".$CWDir."/DBFASTA") or die($!);

print "Downloading FASTA files of $TaxonomyID from Swiss-Prot...\n\n..";

my $link="http://www.uniprot.org/uniprot/?query=(taxonomy%3a".$TaxonomyID.")+AND+reviewed%3ayes&sort=score&format=fasta";
		
	my $mech = WWW::Mechanize->new();
	$mech->get($link);
	my @content = $mech->content;
		
	foreach(@content){
            print AD $_;
	}
        
close(AD);

our $DBN="DBFASTA";

print "Making BLAST database...\n\n";

system("formatdb.exe -p T -i ".$CWDir."/DBFASTA");

print "Blast database has been created...\n\n";

}

1;
