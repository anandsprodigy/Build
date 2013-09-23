package DownSpeciesFasta;

use WWW::Mechanize;
use Cwd;

sub DownSpeciesFasta{
    
our $CWDir=Cwd::_win32_cwd();

open(AB,">$CWDir\\SpeciesFASTA.txt") or die($!);

my ($Genus,$Species)=@_;

print "Downloading species FASTA sequences...\n";

my $link="http://www.uniprot.org/uniprot/?query=(organism%3a$Genus+$Species)+AND+reviewed%3ayes&sort=score&format=fasta";
		
	my $mech = WWW::Mechanize->new();
	$mech->get($link);
	my @content = $mech->content;
		
	foreach(@content){
            print AB $_;
	}
        
close(AB);

print "The Species Fasta proteome has been downloaded\n";

}

1;
