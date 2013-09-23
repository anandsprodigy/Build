package BLASTmake;
use Cwd;


sub BLASTmake{
    
    
    my ($DBN,$QueryFile)=@_;
    
    print "Running BLAST process...\n\n";
    
    system("blastall -p blastp -M BLOSUM50 -e 0.01 -d ".$DBN." -i ".$QueryFile." -o BlastResult.txt -m 8");
    
    print "BLAST has been completed...\n\n";
    
    our $CWDir=Cwd::_win32_cwd();
    
   open(AE,"<".$CWDir."/BlastResult.txt") or die($!);
   open(AF,">".$CWDir."/FilteredResult.txt") or die($!);
   open(AG,">".$CWDir."/BlockInput.txt") or die($!);
    
    foreach(<AE>){
        #print $_;
#sp|P38505|CALBP_ENTHI	sp|P38505|CALBP_ENTHI	100.00	134	0	0	1	134	1	134	5e-094	 259
        #if ($_ =~ m/.*\|(.*)\|.*\t.*\|(.*)\|.*\t(\d+.\d+)\t.*/) {
        if ($_ =~ m/\w+\|(\S+)\|\S+\t\w+\|(\S+)\|\S+\t(\d+.\d+)\t.*/) {
        	#S or s can be problomatic
            if ($3 > 50.00) {
                print AF $_;
                print AG $1."\n";
            }
        }
        
    }
    
    print "BlastResult has been filtered and Block Input Uniprot ID has been written...\n\n";
    
    close(AE);
    close(AF);
    close(AG);
}

1;
