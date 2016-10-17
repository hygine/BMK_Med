#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use Getopt::Long;
use Data::Dumper qw(Dumper);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Math::BigFloat;
my $Title = "exp_id_trans_drow_pheatmap";
my $version="1.0.0";

my ( $in, $symbol, $FC,$FDR,$islog,$scale);
GetOptions(
	"in:s"		=> \$in,
	"symbolfile:s"		=> \$symbol,
	"FC:s"		=>\$FC,
	"FDR:s"		=>\$FDR,
	"scale:s"		=>\$scale,
	"islog"		=>\$islog,
	"help|h"		=> \&USAGE,
) or &USAGE;
&USAGE unless ( $in );
$in = &ABSOLUTE_DIR($in);
$symbol = &ABSOLUTE_DIR($symbol);


##################################
###### symbol ID hash ######
##################################
my %symbol_ID_hash;
#my %para;
#my %sample;

open IN, "$symbol" || die;
while (<IN>) {
	chomp;
	next if (/^$/);
	next if (/^#/);
	next if (/^X./);
	my @tmp = split /\t/, $_;
	my $gene_id = $tmp[0];
	my $symbol_id = $tmp[1];
	$symbol_ID_hash{ $gene_id } = $symbol_id;
}
close IN ;
#print Dumper \%symbol_ID_hash;



open IN, "$in" || die;
open OUT, ">$in.trans.$FDR.$FC.txt" || die;
my $head = <IN>;chomp($head);
my @tmp_array = split /\t+/, $head;chomp(@tmp_array);
my $num = @tmp_array;
my $head_1 = join "\t",@tmp_array[0..$num-4],"gene_symbol";
my $head_2 = join "\t",@tmp_array[$num-3..$num-1];
print OUT "$head_1\n";
my $final_num = 0;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @tmp = split /\t/, $_;
	my $tmp_1 = join "\t",@tmp[0..$num-4];
	my $tmp_2 = join "\t",@tmp[$num-3..$num-1];
	my $nosymbol_num = 0;
	if($symbol_ID_hash{$tmp[0]}){
		$tmp_1 = join "\t",$tmp_1,$symbol_ID_hash{$tmp[0]};
	}else{
		$nosymbol_num ++;
		$tmp_1 = join "\t",$tmp_1,"--$nosymbol_num--";
	}
	#print "aaaaaaa\t$tmp_1\n";
	my $fdr = $tmp[$num-3];
	$fdr = $fdr+0;
	my $fc=abs($tmp[$num-2]) ;
	my $logfc =  &log2($FC);
	
	#print "bbbbbb\t$fdr\t$FDR\t$fc\t$logfc\n";
	if( ($fdr <= $FDR) && ( $fc >= &log2($FC) ) ){
		print OUT "$tmp_1\n";
		$final_num ++;
	}
}
print "selected Gene number: $final_num\n";
close IN;
close OUT;

open HH,">$in.trans.R" or die "$!";
if($islog){
	print HH <<RRRRR;
library(pheatmap)
data=read.table(file = "$in.trans.$FDR.$FC.txt",header=T,sep="\\t",comment.char = "&")
n = ncol(data)
mattr1=as.matrix(data[,-c(1,n)])
rownames(mattr1) = data[,n]
#mattr1 = log2(mattr1)
mattr2=as.matrix(mattr1)
rownames(mattr2) = data[,1]
mattr1_log=log2(mattr1+1)
mattr2_log=log2(mattr2+1)
#pdf(file = "$in.geneID.log.$FDR.$FC.pdf");
pheatmap(mattr2_log,color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, scale = "$scale",
         show_rownames = T,
         show_colnames = T,filename = "$in.geneID.log.$FDR.$FC.pdf",
         main = "DEG Cluster Heatmap")
#dev.off()

#pdf(file = "$in.genesymbol.log.$FDR.$FC.pdf");
pheatmap(mattr1_log,color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, scale = "$scale",
         show_rownames = T,
         show_colnames = T,filename = "$in.genesymbol.log.$FDR.$FC.pdf",
         main = "DEG Cluster Heatmap")
#dev.off()

RRRRR
}else{
	print HH <<RRRRR;
library(pheatmap)
data=read.table(file = "$in.trans.$FDR.$FC.txt",header=T,sep="\\t",comment.char = "&")
n = ncol(data)
mattr1=as.matrix(data[,-c(1,n)])
rownames(mattr1) = data[,n]
#mattr1 = log2(mattr1)
mattr2=as.matrix(mattr1)
rownames(mattr2) = data[,1]
#pdf(file = "$in.geneID.$FDR.$FC.pdf");
pheatmap(mattr2,color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, scale = "$scale",
         show_rownames = T,
         show_colnames = T,filename = "$in.geneID.$FDR.$FC.pdf",
         main = "DEG Cluster Heatmap")
#dev.off()

#pdf(file = "$in.genesymbol.$FDR.$FC.pdf");
pheatmap(mattr1,color=colorRampPalette(c("green","black","red"))(100),
         border_color = NA, scale = "$scale",
         show_rownames = T,
         show_colnames = T,filename = "$in.genesymbol.$FDR.$FC.pdf",
         main = "DEG Cluster Heatmap")
#dev.off()

RRRRR

}
close HH;

`Rscript $in.trans.R`;




sub log2 {
	my $n = shift;
	return log($n)/log(2);
}

sub ABSOLUTE_DIR {
    my $cur_dir = `pwd`;
    chomp($cur_dir);
    my ($in) = @_;
    my $return = "";
    if ( -f $in ) {
        my $dir  = dirname($in);
        my $file = basename($in);
        chdir $dir;
        $dir = `pwd`;
        chomp $dir;
        $return = "$dir/$file";
    }
    elsif ( -d $in ) {
        chdir $in;
        $return = `pwd`;
        chomp $return;
    }
    else {
        warn "Warning just for file and dir\n";
        exit;
    }
    chdir $cur_dir;
    return $return;
}

sub MKDIR {    # &MKDIR($out_dir);
    my ($dir) = @_;
    mkdir( $dir) if ( !-d $dir );
}

sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------
   Program: $Script
   Version: 1.0.0
   Contact: linhj  <linhj\@biomarker.com.cn>
     Usage:
            --in              input  file              [forced ]
            --symbolfile      gene symbol file
            --FC              Fold change cutoff
            --FDR             FDR cutoff
            --scale           if the values should be centered and scaled in either the row 
                              direction or the column direction, or none.[forced: none | row | column]
            --islog           set the number to log2. [T/F]

            --h               help documents

   Example:
            perl $Script -in L02_vs_L01.DEG_final.xls -symbolfile id_trans.xls.xls -FDR 0.001 -FC 120 -scale none -islog
            perl $Script -in L02_vs_L01.DEG_final.xls -symbolfile id_trans.xls.xls -FDR 0.001 -FC 120 -scale row


----------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
