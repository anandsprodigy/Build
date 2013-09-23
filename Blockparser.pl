package BlockParser;

sub BlockParser{

my ($Blockfile)=shift;

open(SA,"<BlockOut.blk") or die($!);
open(SB,">ParsedBlock") or die($!);

while(<SA>){
	if($_ =~ m/.*\)\s(\D+)\s+(\d+)/){		#takes considerations about % clustering
		if($2>50){
			print SB $1."\n";
		}
	}
	
}

print "Blocks parsed successfully...\n\n";

close(SA);
close(SB);

}

1;