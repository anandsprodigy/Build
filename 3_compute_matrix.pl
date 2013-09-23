#!/usr/bin/perl
#*****************

###############################################################################
# Copyright (C) 2011 Université Bordeaux 2                                    #
# 									      #
# Contributors : Claire Lemaitre                 			      #
# 									      #
# Contact : claire.lemaitre@inria.fr                              	      #
# 									      #
# This file is part of matrix-builder.					      #
# 									      #
# matrix-builder is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by	      #
#  the Free Software Foundation, either version 3 of the License, or	      #
# (at your option) any later version.					      #
# 									      #
# matrix-builder is distributed in the hope that it will be useful,	      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of	      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		      #
# GNU General Public License for more details.				      #
# 									      #
# You should have received a copy of the GNU General Public License	      #
# along with matrix-builder.  If not, see <http://www.gnu.org/licenses/>.     #
###############################################################################


#launch the blosum program with the block file and rename and format the output matrix file

use strict;
use warnings;
use Getopt::Std;


sub Usage
{
    print STDOUT "
-----------------------------------------------------------------------
Script to run the program blosum on a block database
-----------------------------------------------------------------------
Usage : $0
Obligatory
   -i  :  input block database file
   -o  :  output matrix file
Facultative
   -c  :  value of the clustering parameter (default='n' for no clustering)
   -s  :  scale of the matrix (default=0 for letting blosum decide)

Example :
$0 -i blocks.blks -c 60 -o MOLLI60

Details :
 - The program blosum outputs 3 files + STDOUT :  blosumX.iij, blosumX.sij, blosumX.qij (X the clustering coefficient). These files are moved in the directory tmp which is deleted at the end (if you want to keep them uncomment the corresponding line)
 - Scale of the matrix (-scale) :
    - 0: let blosum decide based on entropy
    - n: 1/n bits
 - warning : you may have to modify the path to your blosum executable inside the script or put a copy in the current dir
-----------------------------------------------------------------------
		\n";
    exit(0);	
}


MAIN: 
{

    ## You can change the path of the blosum program here :
    my $blosum_program="./blosum";
 
   
    ## Getting the parameters 
    #########################
    my %opts;
    getopts('i:o:c:s:', \%opts);

    my $blk_file=$opts{i};
    &Usage if not defined $blk_file;

    my $output=$opts{o};
    &Usage if not defined $output;

    my $cluster=$opts{c};
    $cluster="n" if not defined $cluster;

    my $scale=$opts{s};
    $scale=0 if not defined $scale;


    ## Temporary dir
    #####################
    my $tmp_dir="tmp/";
    if(! -d $tmp_dir){
	mkdir($tmp_dir);
    }
    my $delete_tmp=1;

    
    ## run blosum
    #############
    my $output_file=$tmp_dir."matrix.out";
    my $blosum_command=$blosum_program." ".$blk_file." 0 9999 ".$cluster." ".$scale." > ".$output_file;
    #print "$blosum_command\n";
    `$blosum_command`;

    ## move the blosum output files
    #########################
    my $source="blosum".$cluster.".*ij";
    my $dest=$tmp_dir;
    `mv $source $dest`;


    ## formatting the final matrix
    ##############################
    my $matrix_iij=$tmp_dir."blosum".$cluster.".iij";
    
    ## Get the triangular matrix
    my %matrix=read_half_matrix($matrix_iij);
    my @aa_order=read_aa_order($matrix_iij);

    ## get the minimal value
    my $mini=mini_value(%matrix);

    ## print the comments + the full matrix
    print_comments($matrix_iij,$output);
    print_full_matrix($output,\%matrix,$mini,@aa_order);

    
    if($delete_tmp){
	`rm -rf $tmp_dir`;
    }

}



sub read_aa_order {
    my $matrix_file=$_[0];
    open(R,"<$matrix_file") or die "cannot open file $matrix_file";

    my @aa_names;
    while(<R>){
	if(! /^\#/){
	    if(/\D\s+\D\s+\D\s+\D/){
		@aa_names=split();
	    }
	}
    }
    close(R);
    return @aa_names;
}

sub read_half_matrix {

    my $matrix_file=$_[0];
    open(R,"<$matrix_file") or die "cannot open file $matrix_file";

    my %matrix=();
    my @aa_names;
    my $i;
    while(<R>){
	if(! /^\#/){
	    if(/\D\s+\D\s+\D\s+\D/){
		@aa_names=split();
		$i=0;
	    }
	    else{
		if(/\d/){
		    my @ligne=split();
		    for(my $j=0;$j<=$#ligne;$j++) {
			$matrix{$aa_names[$i]}{$aa_names[$j]}=$ligne[$j];
		    }
		    $i=$i+1;
		}
	    }
	}
    }
    close(R);

    return %matrix;
}

sub mini_value {

    my (%matrix)=@_;
    my @all=();
    foreach my $aa1 (keys(%matrix)) {
	foreach my $aa2 (keys(%matrix)) {
	    #print "ok\n";
	    push(@all,$matrix{$aa1}{$aa2}) if defined($matrix{$aa1}{$aa2});
	}
    }
    my @all_sorted=sort {$a<=>$b} (@all);
    return $all_sorted[0];
}


sub print_comments {

    my ($matrix_file,$output_file)=@_;
    open(R,"<$matrix_file") or die "cannot open file $matrix_file";
    open(OUTPUT,">$output_file") or die "cannot open file $output_file";
    while(my $ligne=<R>){
	if($ligne =~ /^\#/){
	    print OUTPUT $ligne;
	}
    }
    close(R);
    close(OUTPUT);
}


sub print_full_matrix {

    my ($output_file,$matrix,$mini_value,@aa_order)=@_;
    my $sep="  ";
    push(@aa_order,"*");
    open(OUTPUT,">>$output_file") or die "cannot open file $output_file";
    print OUTPUT "   ".join($sep,@aa_order)."\n";
    $$matrix{"*"}{"*"}=1;
    my $value;
    my $toprint;
    foreach my $aa1 (@aa_order) {
	print OUTPUT $aa1;
	foreach my $aa2 (@aa_order) {
	    if (defined($$matrix{$aa1}{$aa2})) {
		$value=$$matrix{$aa1}{$aa2};
	    }
	    else{
		if (defined($$matrix{$aa2}{$aa1})) {
		    $value=$$matrix{$aa2}{$aa1};
		}
		else{
		    $value=$mini_value;
		}
	    }
	    
	    if(length($value)==1){
		$toprint="  ".$value;
	    }
	    else{
		$toprint=" ".$value;
	    }
	    print OUTPUT $toprint;
	}

	print OUTPUT "\n";
    }
    close(OUTPUT)

}

