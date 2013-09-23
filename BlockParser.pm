package BlockParser;

sub BlockParser{

my ($Blockfile)=shift;

open(SA,"<$Blockfile") or die($!);
open(SB,">ParsedBlock.txt") or die($!);

print "Cleaning BLOCK file...\n\n";


while(<SA>){
	if($_ =~ m/.*\)\s(\D+)\s+(\d+)/){		#takes considerations about % clustering
		if($2>0){
			print SB $1."\n";
			print SB $1."\n";
		}
		
		
	}
	
}

print "Blocks parsed successfully...\n\n";

close(SA);
close(SB);

}

1;