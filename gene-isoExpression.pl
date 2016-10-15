#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cuffnorm_dir,$final_gtf,$outdir);

GetOptions(
				"help|?" =>\&USAGE,
				"c:s"=>\$cuffnorm_dir,
				"f:s"=>\$final_gtf,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($cuffnorm_dir and $final_gtf and $outdir);
#$final_track_dir=&ABSOLUTE_DIR($final_track_dir);
$cuffnorm_dir=&ABSOLUTE_DIR($cuffnorm_dir);
`mkdir $outdir` unless (-d $outdir);
$outdir=&ABSOLUTE_DIR($outdir);
# ------------------------------------------------------------------
# get saved genes and transcripts track list
# ------------------------------------------------------------------
#KB742382.1      XLOC_000001     ENSAPLG00000016114      5801    32652   +       TCONS_00000003  ENSAPLT00000016844      17476   32485   2005
#KB742382.1      XLOC_000003     Duck_newGene_1  121709  128447  +       TCONS_00000010  Duck_newGene_1.2        121709  128447  1092
#my  $raw_track_file=`ls $final_track_dir/Compare.gtf`;
my  $raw_track_file=`ls $final_gtf `;
#my $new_track_file=`ls $final_track_dir/*.newGene_filtered_final_tracking.list`;
chomp($raw_track_file);
#chomp($new_track_file);
my %track;
my %gene_track;
#my %gene_length;
#my %gene_locus_start;
#my %gene_locus_end;

open (TRACK,"$raw_track_file") or die $!;

while (<TRACK>) {
	chomp;
	my @col=split /\t/;
	#next if ($col[7]=~/CUFF\./);
	$col[8]=~/gene_id\s\"([^\"]+)\";\stranscript_id\s\"([^\"]+)\";/;
	my $trans_id=$2;
    my $gene_id=$1;
	$track{$trans_id}{'gene_id'}=$1;
	$track{$trans_id}{'stand'}=$col[6];
	$track{$trans_id}{'trans_id'}=$2;

	$gene_track{$gene_id}{'gene_id'}=$1;
	$gene_track{$gene_id}{'stand'}=$col[6];
	$gene_track{$gene_id}{'trans_id'}=$2;

	#$track{$trans_id}{'length'}=$col[10];
	my $len=$col[4]-$col[3]+1; ###20160728 linhj
	$track{$trans_id}{'locus'}="$col[0]:$col[3]-$col[4]";
	if (exists $track{$trans_id}{'length'}){
		$track{$trans_id}{'length'}+=$len;
        $len = $track{$trans_id}{'length'};
	}else{
		$track{$trans_id}{'length'}=$len;
	}
#	print $track{$trans_id}{'trans_id'};
	#if(exists $track{$gene_id}{'start'}){
    #    push @{$track{$gene_id}{'start'}},$col[3];
    #}else{
		#$track{$gene_id}{'start'}=$col[3];
        push @{$gene_track{$gene_id}{'start'}},$col[3];
	#}
    #print "aaaaaa\t$track{$gene_id}{'start'}\n";
	#if(exists $track{$gene_id}{'end'}){
    #    push @{$track{$gene_id}{'end'}},$col[4];
    #}else{
		#$track{$gene_id}{'end'}=$col[4];
        push @{$gene_track{$gene_id}{'end'}},$col[4];
	#}
	$gene_track{$gene_id}{'chr'}=$col[0];


    #if (exists $track{$gene_id}{'length'}){
	#	push @{$track{$gene_id}{'length'}},$len;
	#}else{
	#	#$track{$gene_id}{'length'}=$len;
        push @{$gene_track{$gene_id}{'length'}},$len;
	#}

}

close TRACK;

# ------------------------------------------------------------------
# abstract data from isoforms.fpkm_tracking and isoforms.count_tracking
# ------------------------------------------------------------------
#tracking_id     class_code      nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage        T1_FPKM T1_conf_lo      T1_conf_hi      T1_status       T2_FPKM
#TCONS_00000003  =       ENSAPLT00000016844      XLOC_000001     ENSAPLG00000016114      TSS3    KB742382.1:5800-32652   2005    -       0.743331        0       3.85708 OK
my $fpkm_file="$cuffnorm_dir/isoforms.fpkm_tracking";
my $gene_fpkm_file="$cuffnorm_dir/genes.fpkm_tracking";
my $count_file="$cuffnorm_dir/isoforms.count_tracking";
my $gene_count_file="$cuffnorm_dir/genes.count_tracking";
my @sample;
my @gene_sample;

open (FPKM,$fpkm_file) or die $!;
chomp(my $head=<FPKM>);
my @attributes_fpkm=(split /\t/,$head);
for (my $i=9;$i<@attributes_fpkm;$i+=4) {
	my ($sample)=$attributes_fpkm[$i]=~/(\S+)_FPKM$/;
	push @sample,$sample;
}

open (GENE_FPKM,$gene_fpkm_file) or die $!;
chomp(my $head2=<GENE_FPKM>);
my @gene_attributes_fpkm=(split /\t/,$head);
for (my $i=9;$i<@gene_attributes_fpkm;$i+=4) {
	my ($sample2)=$gene_attributes_fpkm[$i]=~/(\S+)_FPKM$/;
	push @gene_sample,$sample2;
}

while (<FPKM>) {
	if (exists $track{(split /\t/)[0]}) {

		my @col=split /\t/;

		for (my $i=0;$i<(@col-9)/4;$i++) {
			$track{$col[0]}{$sample[$i]."_FPKM"}=$col[9+4*$i];
		}
	}
}

while (<GENE_FPKM>) {
	if (exists $gene_track{(split /\t/)[0]}) {

		my @col2=split /\t/;

		for (my $i=0;$i<(@col2-9)/4;$i++) {
			$gene_track{$col2[0]}{$gene_sample[$i]."_FPKM"}=$col2[9+4*$i];
		}
	}
}


close FPKM;
close GENE_FPKM;

open (COUNT,$count_file) or dir $!;
<COUNT>;

while (<COUNT>) {
	if (exists $track{(split /\t/)[0]}) {
		chomp;
		my @col=split /\t/;

		for (my $i=0;$i<(@col-1)/5;$i++) {
			$track{$col[0]}{$sample[$i]."_count"}=$col[1+5*$i];
		}
	}
}

close COUNT;

open (GENE_COUNT,$gene_count_file) or dir $!;
<GENE_COUNT>;

while (<GENE_COUNT>) {
	if (exists $gene_track{(split /\t/)[0]}) {
		chomp;
		my @col2=split /\t/;

		for (my $i=0;$i<(@col2-1)/5;$i++) {
			$gene_track{$col2[0]}{$gene_sample[$i]."_count"}=$col2[1+5*$i];
		}
	}
}

close GENE_COUNT;


# ------------------------------------------------------------------
# gene count fpkm output 
# ------------------------------------------------------------------
my $gene_fpkm_counts;

if (@gene_sample>0) {
	for (my $i=0;$i<@gene_sample;$i++) {
		$gene_fpkm_counts.="$gene_sample[$i]\_FPKM\t$gene_sample[$i]\_count\t";
	}
	$gene_fpkm_counts =~ s/\t$//;
}

open (LOG,">$outdir/err2.log") or die $!;
open (OUT,">$outdir/AllSample.genes_expression.xls") or die $!;
#print OUT "#transcript_id\tlength\tstrand\tgene_id\tlocus\t$gene_fpkm_counts\n";
print OUT "#gene_id\tlength\tstrand\tlocus\t$gene_fpkm_counts\n";

#for my $i (sort {$track{$a}{'gene_id'} cmp $track{$b}{'gene_id'}} keys %track) {
#for my $i (sort {($gene_track{$a}=~ /(\d+)$/)[0] <=> ($gene_track{$b}=~ /(\d+)$/)[0]} keys %gene_track) {
for my $i (sort {($a=~ /(\d+)$/)[0] <=> ($b=~ /(\d+)$/)[0]} keys %gene_track) {
#	my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
    my $start = min(@{$gene_track{$i}{'start'}});
    my $end = max(@{$gene_track{$i}{'end'}});
    my $locus = $gene_track{$i}{'chr'}.":".$start."-".$end;
	my $basic_inf=join "\t",($gene_track{$i}{'gene_id'},max(@{$gene_track{$i}{'length'}}),$gene_track{$i}{'stand'},$locus);
	my $print=$basic_inf;

	for (my $s=0;$s<@gene_sample;$s++) {
		if (!defined $gene_track{$i}{$gene_sample[$s]."_FPKM"} || !defined $gene_track{$i}{$gene_sample[$s]."_count"}) {
			delete $gene_track{$i};
			print LOG "WARNING: $i in sample $gene_sample[$s] can't get expression enrichment value.\n";
			next;
		}
		my $sample_inf=join "\t",($gene_track{$i}{$gene_sample[$s]."_FPKM"},$gene_track{$i}{$gene_sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
	}

	#print OUT "$print\n" if (exists $track{$i});
#	print OUT "$print\n" if ((exists $gene_track{$i}) && (min(@{$gene_track{$i}{'length'}}) >= 5)); ###2016.0725  edit by linhj
    print OUT "$print\n" if ((min(@{$gene_track{$i}{'length'}}) >= 1)); ###2016.0725  edit by linhj
}

close OUT;
close LOG;

for (my $s=0;$s<@gene_sample;$s++) {

	open (EXP,">$outdir/$gene_sample[$s].geneExpression.xls") or die $!;
#	print EXP "#transcript_id\tlength\tstrand\tgene_id\tlocus\t$sample[$s]\_FPKM\t$sample[$s]\_count\n";
	print EXP "#gene_id\tlength\tstrand\tlocus\t$gene_sample[$s]\_FPKM\t$gene_sample[$s]\_count\n";

#	for my $i (sort {$gene_track{$a}{'gene_id'} cmp $gene_track{$b}{'gene_id'}} keys %gene_track) {
#    for my $i (sort {$gene_track{$a} cmp $gene_track{$b}} keys %gene_track) {
    for my $i (sort {($a=~ /\_(\d+)/)[0] <=> ($b=~ /\_(\d+)/)[0]} keys %gene_track) {
        #print "a$i\a" ;
#		my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
        my $start = min(@{$gene_track{$i}{'start'}});
        my $end = max(@{$gene_track{$i}{'end'}});
        my $locus = $gene_track{$i}{'chr'}.":".$start."-".$end;

		my $basic_inf=join "\t",($gene_track{$i}{'gene_id'},max(@{$gene_track{$i}{'length'}}),$gene_track{$i}{'stand'},$locus);
		my $print=$basic_inf;

		my $sample_inf=join "\t",($gene_track{$i}{$gene_sample[$s]."_FPKM"},$gene_track{$i}{$gene_sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
		#print EXP "$print\n";
		print EXP "$print\n"  if (min(@{$gene_track{$i}{'length'}}) >= 1) ; ###2016.0725  edit by linhj
	}

	close EXP;
}



# ------------------------------------------------------------------
# trans count fpkm output 
# ------------------------------------------------------------------
my $fpkm_counts;

if (@sample>0) {
	for (my $i=0;$i<@sample;$i++) {
		$fpkm_counts.="$sample[$i]\_FPKM\t$sample[$i]\_count\t";
	}

	$fpkm_counts=~s/\t$//;
}

open (LOG,">$outdir/err.log") or die $!;
open (OUT,">$outdir/AllSample.isoforms_expression.xls") or die $!;
#print OUT "#transcript_id\tlength\tstrand\tgene_id\tlocus\t$fpkm_counts\n";
print OUT "#transcript_id\tlength\tstrand\tlocus\t$fpkm_counts\n";

#for my $i (sort {$track{$a}{'gene_id'} cmp $track{$b}{'gene_id'}} keys %track) {

#for my $i (sort {($track{$a}{'gene_id'} =~ /\_(\d+)/g)[0] <=> ($track{$b}{'gene_id'} =~ /\_(\d+)/g)[0]} keys %track) {
for my $i (sort {($track{$a}{'trans_id'} =~ /\_(\d+)/g)[0] <=> ($track{$b}{'trans_id'} =~ /\_(\d+)/g)[0]} keys %track) {
#	my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
	my $basic_inf=join "\t",($track{$i}{'trans_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'locus'});  #################
	my $print=$basic_inf;

	for (my $s=0;$s<@sample;$s++) {
		if (!defined $track{$i}{$sample[$s]."_FPKM"} || !defined $track{$i}{$sample[$s]."_count"}) {
			delete $track{$i};
			print LOG "WARNING: $i in sample $sample[$s] can't get expression enrichment value.\n";
			next;
		}
		my $sample_inf=join "\t",($track{$i}{$sample[$s]."_FPKM"},$track{$i}{$sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
	}

	#print OUT "$print\n" if (exists $track{$i});
	print OUT "$print\n" if ((exists $track{$i}) && ($track{$i}{'length'} >= 5)); ###2016.0725  edit by linhj
}

close OUT;
close LOG;

for (my $s=0;$s<@sample;$s++) {

	open (EXP,">$outdir/$sample[$s].isoExpression.xls") or die $!;
#	print EXP "#transcript_id\tlength\tstrand\tgene_id\tlocus\t$sample[$s]\_FPKM\t$sample[$s]\_count\n";
	print EXP "#transcript_id\tlength\tstrand\tlocus\t$sample[$s]\_FPKM\t$sample[$s]\_count\n";

#	for my $i (sort {$track{$a}{'gene_id'} cmp $track{$b}{'gene_id'}} keys %track) {
#   for my $i (sort {int(($track{$a}{'gene_id'} =~ /(\d+)$/)[0]) <=> int(($track{$b}{'gene_id'} =~ /(\d+)$/)[0])} keys %track) {
    for my $i (sort {($track{$a}{'trans_id'} =~ /(\d+)$/)[0] <=> ($track{$b}{'trans_id'} =~ /(\d+)$/)[0]} keys %track) {
#		my $basic_inf=join "\t",($track{$i}{'iso_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'gene_id'},$track{$i}{'locus'});
		my $basic_inf=join "\t",($track{$i}{'trans_id'},$track{$i}{'length'},$track{$i}{'stand'},$track{$i}{'locus'});  ######################
		my $print=$basic_inf;

		my $sample_inf=join "\t",($track{$i}{$sample[$s]."_FPKM"},$track{$i}{$sample[$s]."_count"});
		$print=$print."\t".$sample_inf;
		#print EXP "$print\n";
		print EXP "$print\n"  if ($track{$i}{'length'} >= 1) ; ###2016.0725  edit by linhj
	}

	close EXP;
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:	$0
     Version:	$version
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2014.04.15
      Modify:	linhj <linhj\@biomarker.com.cn> 
      Update:	2016.09.22
 Description:	This program is used to abstact genes & transcripts expression from Cuffnorm or Cuffdiff results, enrichment per sample.
       Usage:
		Options:
		-f <str>	input gtf file,forced

		-c <str>	input directory,Cuffdiff/Cuffnorm directory,forced

		-o <str>	output directory,forced

		-h		help

        Example:
		perl $0 -f Compare/Compare.gtf -c Tophat_Cufflinks/Cuffnorm/ -o gene-isoExpression/

USAGE
	print $usage;
	exit;
}
