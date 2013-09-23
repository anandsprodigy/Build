require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
use WWW::Mechanize;
use LWP::UserAgent;
my $mech;
my $ua = LWP::UserAgent->new();


sub blocksretrieve
{
   my $multiplefasta = shift;
   my $proxy=shift;
   if($proxy eq "1")
   {
      $ua->proxy(['http'],"");
      $mech = WWW::Mechanize->new();
   }
   else
   {
    $mech = WWW::Mechanize->new(noproxy => 1);
   }
   my $id;
$mech->get("http://blocks.fhcrc.org/blocks/make_blocks.html");
$mech->set_fields(desc =>'gmail',sequences =>"$multiplefasta");
my $result=$mech->submit();

my @content=$result->content();

my $SubmitID;

foreach(@content){
	#if($_ =~ m/.*bm_format.pl\?(.*)\">here.*/){
	if($_ =~ m/.*bm_format.pl\?(\d+)\">here.*/){
			$SubmitID=$1;
	}
}

print "ID ".$SubmitID;

my $args ="../tmp/bm/$SubmitID/$SubmitID.mblks";
my $resp = $ua->get("http://blocks.fhcrc.org/blocks-bin/catfile.sh?../tmp/bm/9115/9115.mblks");
my @google=%{$resp};
$google[3]=~s/(<PRE>\n|<\/PRE>\n)//g;
my @returnarray = split("//\n",$google[3]);
return @returnarray;
}


1;