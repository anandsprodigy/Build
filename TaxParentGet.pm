 package TaxParentGet;
  
    use WWW::Mechanize;
    
    BEGIN{
        $ENV{"BLAST"}="C:\\Program Files\\NCBI\\bin\\blastall.exe";
        our $TaxFile="C:\\Program Files\\NCBI\\Tax\\nodes.dmp";
    }
    
    sub GetParentTaxID{
    
        my ($TaxID)=@_;
    
        open(TAX,"<",$TaxFile) or die($!);

        foreach(<TAX>){
            
            if ($_ =~ m/(\d+)\t\|\t(\d+)\t\|\t(\w+).*/) {
                
                
                if ($1 == $TaxID && $3 eq 'species' || $1 == $TaxID && $3 eq 'subspecies') {
                    $TaxID=$2;
                    GetParentTaxID($TaxID);
                    last;
                }elsif($1 == $TaxID && $3 eq 'genus' || $1 == $TaxID && $3 eq 'subgenus'){
                    $TaxID=$2;
                    GetParentTaxID($TaxID);
                    last;
                }elsif($1 == $TaxID && $3 eq 'family' || $1 == $TaxID && $3 eq 'subfamily' || $1 == $TaxID && $3 eq 'superfamily'){
                    $TaxID=$2;
                    GetParentTaxID($TaxID);
                    last;
                }elsif($1 == $TaxID && $3 eq 'order' ||  $1 == $TaxID && $3 eq 'suborder' || $1 == $TaxID && $3 eq 'superorder' || $1 == $TaxID && $3 eq 'infraorder'){
                    $TaxID=$2;
                    GetParentTaxID($TaxID);
                    last;
                }elsif($1 == $TaxID && $3 eq 'class' || $1 == $TaxID && $3 eq 'subclass' || $1 == $TaxID && $3 eq 'superclass' || $1 == $TaxID && $3 eq 'infraclass'){
                    return $TaxID;
                }elsif(($1 == $TaxID && $3 eq 'phylum') || ($1 == $TaxID && $3 eq 'subphylum') || ($1 == $TaxID && $3 eq 'superphylum')){
                    return $TaxID;
                }elsif($1 == $TaxID && $3 eq 'no rank' || $1 == $TaxID && $3 eq 'varietas'){
                    return $TaxID;
                }elsif($1 == $TaxID && $3 eq 'superkingdom'|| $1 == $TaxID && $3 eq 'kingdom' ){
                    return $TaxID;
                }
                
            }
            
        }


    close(TAX);
    }
    
    1;