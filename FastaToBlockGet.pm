package FastaToBlockGet;

use WWW::Mechanize;
use LWP::UserAgent;
use LWP::Simple;
use LWP::UserAgent;
my $ua = LWP::UserAgent->new();


sub FastaToBlockGet{

#CWD

my ($Email,$Proxy,$BlockIn,$BlockOut)=(shift,shift,shift,shift);

my $filecontent;

open(AS,"<$BlockIn") or die($!);
open(AT,">$BlockOut") or die($!);


while(<AS>){
	$filecontent.=$_;
}

print "Sending the FASTA file to server,wait while BLOCK server processes result...\n\n";


my $url="http://blocks.fhcrc.org/blocks/make_blocks.html";
my $Mechanize=WWW::Mechanize->new();

if($Proxy eq "http://127.0.0.0:8080/"){
	print "No proxy server used...\n\n";
}else{
	$Mechanize->proxy(['http', 'ftp'], 'http://127.0.0.0:8080/');
}

$Mechanize->get($url);
$Mechanize->set_fields(address=>$Email,desc =>'fffs',sequences =>$filecontent);
my $result=$Mechanize->submit();


my @content=$result->content();

foreach(@content){
	#if($_ =~ m/.*bm_format.pl\?(.*)\">here.*/){
	if($_ =~ m/.*bm_format.pl\?(\d+)\">here.*/){
			$SubmitID=$1;
	}
}

print "Block submit ID : \t".$SubmitID."\n\n";

#my $BlockGet=$Mechanize->follow_link( text => 'here');

#my @Block=$Mechanize->text();

#foreach(@Block){
#	print AT $_;
#}

print "Block retrieved successfully and results were sent to $Email...\n\n";


chomp($SubmitID);
chomp($SubmitID);

print $SubmitID."\n\n\n";

my $mblk= "http://blocks.fhcrc.org/blocks-bin/catfile.sh?../tmp/bm/$SubmitID/$SubmitID.gblks";
my $string="http://blocks.fhcrc.org/blocks-bin/catfile.sh?../tmp/bm/$SubmitID/$SubmitID.mblks";

$Mechanize->get($mblk);

my $filename="c:\\simple.html";

my @content1 = $Mechanize->content(format => "gblk");

print $Mechanize->content(format => "gblk");

foreach(@content1){
	if($_ =~ m/.*\)\s(\D+)\s+(\d+)/){		#takes considerations about % clustering
		print "$_\n";
	}
	
	print $_;
}



close(AS);
close(AT);
}

1;