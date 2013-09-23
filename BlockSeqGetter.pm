package BlockSeqGetter;

use WWW::Mechanize;
use Cwd;

sub BlockSeqGetter{
    my ($BlockInput)=@_;
    
    my $CWDir=Cwd::_win32_cwd();
    
    print "Retrieving FASTA sequence for BLOCKS server...\n\n";
    
    open(AH,"<".$CWDir."/".$BlockInput)or die($!);
    open(AI,">".$CWDir."/BlockInputFasta.fasta")or die($!);
    
    foreach(<AH>){
    
   	my $link="http://www.uniprot.org/uniprot/?query=".$_."&format=fasta";

		
	my $mech = WWW::Mechanize->new();
	$mech->get($link);
	my @content = $mech->content();
	#print $content[0];	#to print FASTA block of seq
	foreach(@content){
           print AI $_;
	}
        
    }
    
    close(AH);
    close(AI);
    
    print "BLOCKS input Fasta file created successfully...\n\n";
    
}


1;