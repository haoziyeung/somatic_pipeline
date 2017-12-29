use strict;
use warnings;
use utf8;
use lib '/home/yanghao/perl5/lib/perl5/';
use Getopt::Long qw/GetOptions/;
use feature 'say';
use Data::Dumper;
use File::Basename;
my $base = dirname $0;

my %cgi1;
open CGI1,"sed 1d $base/cgi_biomarkers_per_variant.tsv|";

while(<CGI1>){
	chomp;
	my @a = split /\t/;
	next unless $a[16];
	next if $a[16] =~ /del|ins|_/;
	$a[16] =~ /(chr.*):g.(\d+)(\S)>(\S)/;
	my $k = join '-',($1,$2,$3,$4);
	$cgi1{$k}{'Association'} .= "$a[3];";
	$cgi1{$k}{'Biomarker'} .= "$a[4];";
	$cgi1{$k}{'Drug'} .= "$a[6];";
	$cgi1{$k}{'Drug family'} .= "$a[7];";
	$cgi1{$k}{'Drug full name'} .= "$a[8];";
	$cgi1{$k}{'Drug status'} .= "$a[9];";
	$cgi1{$k}{'Evidence level'} .= "$a[10];";
	$cgi1{$k}{'Source'} .= "$a[14];";
	$cgi1{$k}{'Primary Tumor type'} .= "$a[23];";
}

say "Association\tBiomarker\tDrug\tDrug family\tDrug full name\tDrug status\tEvidence level\tSource\tPrimary Tumor type\tChr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.ensGene\tGene.ensGene\tGeneDetail.ensGene\tExonicFunc.ensGene\tAAChange.ensGene\tcytoBand\tesp6500siv2_all\t1000g2015aug_all\tavsnp150\tgnomAD_genome_ALL\tgnomAD_genome_AFR\tgnomAD_genome_AMR\tgnomAD_genome_ASJ\tgnomAD_genome_EAS\tgnomAD_genome_FIN\tgnomAD_genome_NFE\tgnomAD_genome_OTH\tgnomAD_exome_ALL\tgnomAD_exome_AFR\tgnomAD_exome_AMR\tgnomAD_exome_ASJ\tgnomAD_exome_EAS\tgnomAD_exome_FIN\tgnomAD_exome_NFE\tgnomAD_exome_OTH\tgnomAD_exome_SAS\tcosmic70\tFilter\tGT\tAF\tAD\tREF_F1R2\tREF_F2R1\tALT_F1R2\tALT_F2R1\tGT\tAF\tAD\tREF_F1R2\tREF_F2R1\tALT_F1R2\tALT_F2R1";

my $mul = shift;
open F,"sed 1d $mul|";
while(<F>){
	chomp;
	my @a = split /\t/;
	next if $a[-5] =~ /alt_allele_in_normal/;

	my $k = "chr$a[0]-$a[1]-$a[3]-$a[4]";
	my $med;
	if(exists $cgi1{$k}){
		$med = join "\t",($cgi1{$k}{'Association'},$cgi1{$k}{'Biomarker'},$cgi1{$k}{'Drug'},$cgi1{$k}{'Drug family'},$cgi1{$k}{'Drug full name'},$cgi1{$k}{'Drug status'},$cgi1{$k}{'Evidence level'},$cgi1{$k}{'Source'},$cgi1{$k}{'Primary Tumor type'});
	}else{
		$med = join "\t",('.')x9;
	}

	my @k = split /:/,$a[-3];
	my @tumor = split /:/,$a[-2];
	my @normal = split /:/,$a[-1];

	my %t_h;
	@t_h{@k} = @tumor;

	my %n_h;
	@n_h{@k} = @normal;

	#print Dumper \%t_h;
	say join "\t",(
		$med,@a[0..35,102,-5],
		$t_h{'GT'},$t_h{'AF'},$t_h{'AD'},$t_h{'REF_F1R2'},$t_h{'REF_F2R1'},$t_h{'ALT_F1R2'},$t_h{'ALT_F2R1'},
		$n_h{'GT'},$n_h{'AF'},$n_h{'AD'},$n_h{'REF_F1R2'},$n_h{'REF_F2R1'},$n_h{'ALT_F1R2'},$n_h{'ALT_F2R1'},
	);
}
