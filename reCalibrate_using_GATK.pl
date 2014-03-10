#!/usr/bin/perl
#v2.0
#Vishal Koparde	

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use Vutil;

sub usage;
sub getFileSize;

our $java="/usr/global/blp/jdk1.7.0_45/bin/java";
our $gatkBaseDir="/usr/global/blp/GenomeAnalysisTK-2.8-1-g932cd3a/";
our $gatk_60g="$java -Xmx60g -jar ${gatkBaseDir}/GenomeAnalysisTK.jar ";
our $reference="/data/refdb/genomes/Homo_sapiens/UCSC/hg18.fa";
our $samtools="/data2/illumina/tools/samtools-0.1.18/samtools";

usage() if (scalar @ARGV==0);

my ($bamFileFolder,$force,$help);

my $result = GetOptions ( "bamFileFolder=s" => \$bamFileFolder,
                          "force|f" => \$force,
			  "help|h" => \$help
			  );

usage() if ($help);

my @bamFileList = &getBamFiles($bamFileFolder);
my $reAlignFolder="${bamFileFolder}/reAligned";
my $reCalibrateFolder="${bamFileFolder}/reCalibrated";
my $unifiedGenotyperFolder="${bamFileFolder}/unifiedGenotyper";
my $lastjobid;
my @jobids;


foreach my $bamFile (@bamFileList) {
	print "$bamFile\n";

	my $cmd1="$samtools index $bamFile";
	system($cmd1) unless(-e "${bamFile}.bai");

	my ($bamFileName,$bamFilePath,$bamFileExt);
	($bamFileName,$bamFilePath,$bamFileExt) = fileparse($bamFile,qr"\..[^.]*$");

	$lastjobid=-1;

	$lastjobid = &reAlign("${bamFileFolder}/${bamFile}",$reAlignFolder);
	push @jobids,$lastjobid unless $lastjobid == -1;

	$lastjobid = &reCalibrate("${reAlignFolder}/${bamFileName}.realigned.bam",$reCalibrateFolder,$lastjobid);
	push @jobids,$lastjobid unless $lastjobid == -1;;

}

sleep(10);
#wait for all jobs to finish
foreach my $jobid (@jobids) {
	while (`/usr/global/sge/bin/lx24-amd64/qstat|grep $jobid|wc -l` == 1) {
		sleep(10);
	}
}

#check for all reCalibrated bam files
foreach my $bamFile (@bamFileList) {
	my ($bamFileName,$bamFilePath,$bamFileExt);
	($bamFileName,$bamFilePath,$bamFileExt) = fileparse($bamFile,qr"\..[^.]*$");
	my $reCalBamFile="${reCalibrateFolder}/${bamFileName}.realigned.recal.bam";
	my $reCalBaiFile="${reCalibrateFolder}/${bamFileName}.realigned.recal.bai";
	Vutil::fileCheck($reCalBamFile,"Something went wrong with reCalibration!");
	Vutil::fileCheck($reCalBaiFile,"Something went wrong with reCalibration!");
}

#UnifiedGenoTyper
$lastjobid=&unifiedGenotyper($reCalibrateFolder,$unifiedGenotyperFolder);
Vutil::fileCheck("${unifiedGenotyperFolder}/raw.vcf","Something went wrong with UnifiedGenotyper!");
Vutil::fileCheck("${unifiedGenotyperFolder}/raw_SNP.vcf","Something went wrong with UnifiedGenotyper!");
Vutil::fileCheck("${unifiedGenotyperFolder}/raw_INDEL.vcf","Something went wrong with UnifiedGenotyper!");

#VariantRecalibrator
$lastjobid=&variantRecalibrator($unifiedGenotyperFolder);
Vutil::fileCheck("${unifiedGenotyperFolder}/raw_SNP.recal.vcf","Something went wrong with VariantRecalibrator!");

#FilterIndels
$lastjobid=&filterIndels($unifiedGenotyperFolder);
Vutil::fileCheck("${unifiedGenotyperFolder}/raw_INDEL.filtered.vcf","Something went wrong while filtering INDELs!");

#combineVariants
$lastjobid=&combineVariants($unifiedGenotyperFolder);


exit;


sub combineVariants {
#                     _     _          __     __         _             _       
#  ___ ___  _ __ ___ | |__ (_)_ __   __\ \   / /_ _ _ __(_) __ _ _ __ | |_ ___ 
# / __/ _ \| '_ ` _ \| '_ \| | '_ \ / _ \ \ / / _` | '__| |/ _` | '_ \| __/ __|
#| (_| (_) | | | | | | |_) | | | | |  __/\ V / (_| | |  | | (_| | | | | |_\__ \
# \___\___/|_| |_| |_|_.__/|_|_| |_|\___| \_/ \__,_|_|  |_|\__,_|_| |_|\__|___/
#        
my ($ugDir) = @_;
my $snpvcf = "${ugDir}/raw_SNP.recal.vcf";
my $indelvcf = "${ugDir}/raw_INDEL.filtered.vcf";
my $outvcf = "${ugDir}/analysisReady.vcf";

my $cv_cmd1="$gatk_60g -T CombineVariants";
$cv_cmd1.=" -R $reference";
$cv_cmd1.=" --variant $snpvcf";
$cv_cmd1.=" --variant $indelvcf";
$cv_cmd1.=" -genotypeMergeOptions UNIQUIFY";
$cv_cmd1.=" -o $outvcf";
my $cv_job1=new Qsub(name=>"CV1",outFile=>"CV1.tmp.out",wd=>$ugDir,cmd=>$cv_cmd1);
$cv_job1->submit();
$cv_job1->waitForCompletion();
return $cv_job1->{jobid};
}




sub filterIndels {
#  __ _ _ _           ___           _      _     
# / _(_) | |_ ___ _ _|_ _|_ __   __| | ___| |___ 
#| |_| | | __/ _ \ '__| || '_ \ / _` |/ _ \ / __|
#|  _| | | ||  __/ |  | || | | | (_| |  __/ \__ \
#|_| |_|_|\__\___|_| |___|_| |_|\__,_|\___|_|___/
#                                                                                           
my ($ugDir) = @_;
my $inputvcf = "${ugDir}/raw_INDEL.vcf";
my $filteredvcf = "${ugDir}/raw_INDEL.filtered.vcf";
return -1 if ( (-e $filteredvcf) and (getFileSize($filteredvcf) > 0) and not defined $force);

my $fi_cmd1="$gatk_60g -T VariantFiltration";
$fi_cmd1.=" -R $reference";
$fi_cmd1.=" --variant $inputvcf";
$fi_cmd1.=" --filterExpression \"QD < 2.0\"";
$fi_cmd1.=" --filterExpression \"ReadPosRankSum < -20.0\"";
$fi_cmd1.=" --filterExpression \"InbreedingCoeff < -0.8\"";
$fi_cmd1.=" --filterExpression \"FS > 200.0\"";
$fi_cmd1.=" --filterName QDFilter";
$fi_cmd1.=" --filterName ReadPosFilter";
$fi_cmd1.=" --filterName InbreedingFilter";
$fi_cmd1.=" --filterName FSFilter";
my $fi_job1=new Qsub(name=>"FI1",outFile=>"FI1.tmp.out",wd=>$ugDir,cmd=>$fi_cmd1);
$fi_job1->submit();
$fi_job1->waitForCompletion();
return $fi_job1->{jobid};
}


sub variantRecalibrator {
#                 _             _   ____                _ _ _               _             
#__   ____ _ _ __(_) __ _ _ __ | |_|  _ \ ___  ___ __ _| (_) |__  _ __ __ _| |_ ___  _ __ 
#\ \ / / _` | '__| |/ _` | '_ \| __| |_) / _ \/ __/ _` | | | '_ \| '__/ _` | __/ _ \| '__|
# \ V / (_| | |  | | (_| | | | | |_|  _ <  __/ (_| (_| | | | |_) | | | (_| | || (_) | |   
#  \_/ \__,_|_|  |_|\__,_|_| |_|\__|_| \_\___|\___\__,_|_|_|_.__/|_|  \__,_|\__\___/|_|   
#                                                                                         
my ($ugDir) = @_;
my $inputvcf = "${ugDir}/raw_SNP.vcf";
my $recalFile = "${ugDir}/raw_SNP.recal";
my $tranchesFile = "${ugDir}/raw_SNP.tranches";
my $recalvcf = "${ugDir}/raw_SNP.recal.vcf";
return -1 if ( (-e $recalFile) and (getFileSize($recalFile) > 0) and (-e $tranchesFile) and (getFileSize($tranchesFile) > 0) and not defined $force);

my $vr_cmd1="$gatk_60g -T VariantRecalibrator";
$vr_cmd1.=" -R $reference";
$vr_cmd1.=" -input $inputvcf";
$vr_cmd1.=" --maxGaussians 6";
$vr_cmd1.=" -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 /data/refdb/GATK_bundle/ftp.broadinstitute.org/bundle/2.8/hg18/hapmap_3.3.hg18.vcf";
$vr_cmd1.=" -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 /data/refdb/GATK_bundle/ftp.broadinstitute.org/bundle/2.8/hg18/1000G_omni2.5.hg18.vcf";
$vr_cmd1.=" -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=6.0 /data/refdb/GATK_bundle/ftp.broadinstitute.org/bundle/2.8/hg18/dbsnp_138.hg18.vcf";
$vr_cmd1.=" -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff";
$vr_cmd1.=" -mode SNP";
$vr_cmd1.=" -recalFile $recalFile";
$vr_cmd1.=" -tranchesFile $tranchesFile";
$vr_cmd1.=" -rscriptFile raw_SNP.plots.R";
my $vr_job1=new Qsub(name=>"VR1",outFile=>"VR1.tmp.out",wd=>$ugDir,cmd=>$vr_cmd1);
$vr_job1->submit();

my $vr_cmd2="$gatk_60g -T ApplyRecalibration";
$vr_cmd2.=" -R $reference";
$vr_cmd2.=" -input raw_SNP.vcf";
$vr_cmd2.=" -recalFile $recalFile";
$vr_cmd2.=" -tranchesFile $tranchesFile";
$vr_cmd2.=" -o $recalvcf";
my $vr_job2=new Qsub(name=>"VR2",outFile=>"VR2.tmp.out",wd=>$ugDir,cmd=>$vr_cmd2,waitfor=>($vr_job1->{jobid}));
$vr_job2->submit();

$vr_job2->waitForCompletion();
return $vr_job2->{jobid}; 
}


sub unifiedGenotyper {
#             _  __ _          _  ____                  _                         
# _   _ _ __ (_)/ _(_) ___  __| |/ ___| ___ _ __   ___ | |_ _   _ _ __   ___ _ __ 
#| | | | '_ \| | |_| |/ _ \/ _` | |  _ / _ \ '_ \ / _ \| __| | | | '_ \ / _ \ '__|
#| |_| | | | | |  _| |  __/ (_| | |_| |  __/ | | | (_) | |_| |_| | |_) |  __/ |   
# \__,_|_| |_|_|_| |_|\___|\__,_|\____|\___|_| |_|\___/ \__|\__, | .__/ \___|_|   
#                                                           |___/|_|              
my ($recalDir,$ugDir) = @_;
my $outvcf="${ugDir}/raw.vcf";
my $outvcf1="${ugDir}/raw_SNP.vcf";
my $outvcf2="${ugDir}/raw_INDEL.vcf";

return -1 if ( (-e $outvcf) and (getFileSize($outvcf) > 0) and (-e $outvcf1) and (getFileSize($outvcf1) > 0) and (-e $outvcf2) and (getFileSize($outvcf2) > 0) and not defined $force);

mkdir($ugDir) unless(-d $ugDir);

my @bamFileList = &getBamFiles($recalDir);
my $ug_cmd1="$gatk_60g -T UnifiedGenotyper";
$ug_cmd1.=" -R $reference";
$ug_cmd1.=" --genotype_likelihoods_model BOTH";
$ug_cmd1.=" --output_mode EMIT_VARIANTS_ONLY";
$ug_cmd1.=" --min_base_quality_score 25";
$ug_cmd1.=" -dcov 500";
$ug_cmd1.=" -nt 20";
foreach my $bamFile (@bamFileList) {
	$ug_cmd1.=" -I ${recalDir}/${bamFile}";
}
$ug_cmd1.= " -o ${ugDir}/raw.vcf";
my $ug_job1=new Qsub(name=>"UG1",outFile=>"UG1.tmp.out",wd=>$ugDir,cmd=>$ug_cmd1,nproc=>"20");
$ug_job1->submit();

my $ug_cmd2="$gatk_60g -T SelectVariants";
$ug_cmd2.= " -R $reference";
$ug_cmd2.= " --variant ${ugDir}/raw.vcf";
$ug_cmd2.= " -o ${ugDir}/raw_SNP.vcf";
$ug_cmd2.= " -selectType SNP";
my $ug_job2=new Qsub(name=>"UG2",outFile=>"UG2.tmp.out",wd=>$ugDir,cmd=>$ug_cmd2,waitfor=>($ug_job1->{jobid}));
$ug_job2->submit();

my $ug_cmd3="$gatk_60g -T SelectVariants";
$ug_cmd3.= " -R $reference";
$ug_cmd3.= " --variant ${ugDir}/raw.vcf";
$ug_cmd3.= " -o ${ugDir}/raw_INDEL.vcf";
$ug_cmd3.= " -selectType INDEL";
my $ug_job3=new Qsub(name=>"UG3",outFile=>"UG3.tmp.out",wd=>$ugDir,cmd=>$ug_cmd3,waitfor=>($ug_job1->{jobid}));
$ug_job3->submit();

$ug_job2->waitForCompletion();
$ug_job3->waitForCompletion();
return $ug_job3->{jobid};
}



sub getBamFiles{
#            _   ____                  _____ _ _           
#  __ _  ___| |_| __ )  __ _ _ __ ___ |  ___(_) | ___  ___ 
# / _` |/ _ \ __|  _ \ / _` | '_ ` _ \| |_  | | |/ _ \/ __|
#| (_| |  __/ |_| |_) | (_| | | | | | |  _| | | |  __/\__ \
# \__, |\___|\__|____/ \__,_|_| |_| |_|_|   |_|_|\___||___/
# |___/
my ($bamFilePath) = @_;
my @bamFileList;
opendir (DIR,$bamFilePath);
while (my $file=readdir(DIR)) {
	next if $file eq '.' or $file eq '..';
	next unless ($file =~ m/\.bam$/);
	push @bamFileList,$file;
}
closedir (DIR);
return @bamFileList;
}




sub reAlign {
#             _    _ _             
# _ __ ___   / \  | (_) __ _ _ __  
#| '__/ _ \ / _ \ | | |/ _` | '_ \ 
#| | |  __// ___ \| | | (_| | | | |
#|_|  \___/_/   \_\_|_|\__, |_| |_|
#                      |___/       
my ($bamFile,$outDirPath) = @_;

my ($bamFileName,$bamFilePath,$bamFileExt);
($bamFileName,$bamFilePath,$bamFileExt) = fileparse($bamFile,qr"\..[^.]*$");

my $intervalsFile="${outDirPath}/${bamFileName}.intervals";
my $outBamFile="${outDirPath}/${bamFileName}.realigned.bam";
my $outBaiFile="${outDirPath}/${bamFileName}.realigned.bai";

return -1 if ((-e $outBamFile) and (getFileSize($outBamFile) > 0) and (-e $outBaiFile) and (getFileSize($outBaiFile) > 0) and (not defined $force));

mkdir($outDirPath) unless(-d $outDirPath);

my $reAlign_cmd1="$gatk_60g -T RealignerTargetCreator";
$reAlign_cmd1.=" -R $reference";
$reAlign_cmd1.=" -I $bamFile";
$reAlign_cmd1.=" --out $intervalsFile";
$reAlign_cmd1.=" -nt 4";
#print "$reAlign_cmd1\n";
#system($reAlign_cmd1);
my $jname="A1.${bamFileName}";
my $jout="A1.${bamFileName}.out";
my $reAlign_job1=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$reAlign_cmd1,nproc=>"4");
$reAlign_job1->submit();

my $reAlign_cmd2="$gatk_60g -T IndelRealigner";
$reAlign_cmd2.=" -R $reference";
$reAlign_cmd2.=" -I $bamFile";
$reAlign_cmd2.=" -targetIntervals $intervalsFile";
$reAlign_cmd2.=" -o $outBamFile";
#print "$reAlign_cmd2\n";
#system($reAlign_cmd2);
$jname="A2.${bamFileName}";
$jout="A2.${bamFileName}.out";
my $reAlign_job2=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$reAlign_cmd2,waitfor=>($reAlign_job1->{jobid}));
$reAlign_job2->submit();

#$jname="A3.${bamFileName}";
#$jout="A3.${bamFileName}.out";
#my $reAlign_cmd3="$samtools index ${outBamFile}";
#my $reAlign_job3=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$reAlign_cmd3,waitfor=>($reAlign_job2->{jobid}));
#$reAlign_job3->submit();

return $reAlign_job2->{jobid};

}




sub reCalibrate {
#           ____      _ _ _               _       
# _ __ ___ / ___|__ _| (_) |__  _ __ __ _| |_ ___ 
#| '__/ _ \ |   / _` | | | '_ \| '__/ _` | __/ _ \
#| | |  __/ |__| (_| | | | |_) | | | (_| | ||  __/
#|_|  \___|\____\__,_|_|_|_.__/|_|  \__,_|\__\___|
#       
my ($bamFile,$outDirPath,$lastjobid) = @_; 

my ($bamFileName,$bamFilePath,$bamFileExt);
($bamFileName,$bamFilePath,$bamFileExt) = fileparse($bamFile,qr"\..[^.]*$");

my $recalCSVFile="${outDirPath}/${bamFileName}.recal.csv";
my $recalBamFile="${outDirPath}/${bamFileName}.recal.bam";
my $recalBaiFile="${outDirPath}/${bamFileName}.recal.bai";

return -1 if ((-e $recalBamFile) and (getFileSize($recalBamFile) > 0) and (-e $recalBaiFile) and (getFileSize($recalBaiFile) > 0) and (not defined $force));

mkdir($outDirPath) unless(-d $outDirPath);

my $recal_cmd1="$gatk_60g -T BaseRecalibrator";
$recal_cmd1.=" -R $reference";
$recal_cmd1.=" -knownSites /data/refdb/GATK_bundle/ftp.broadinstitute.org/bundle/2.8/hg18/dbsnp_138.hg18.vcf";
$recal_cmd1.=" -I $bamFile";
#$recal_cmd1.=" -cov QualityScoreCovariate";
#$recal_cmd1.=" -cov CycleCovariate";
#$recal_cmd1.=" -cov ContextCovariate";
$recal_cmd1.=" -o $recalCSVFile";
#print "$recal_cmd1\n";
#system($recal_cmd1);
my $jname="C1.${bamFileName}";
my $jout="C1.${bamFileName}.out";
my $reCalibrate_job1;
if ($lastjobid!=-1) {
$reCalibrate_job1=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$recal_cmd1,waitfor=>$lastjobid);
} else {
$reCalibrate_job1=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$recal_cmd1);
}
$reCalibrate_job1->submit();

my $recal_cmd2="$gatk_60g -T PrintReads";
$recal_cmd2.=" -R $reference";
$recal_cmd2.=" -I $bamFile";
$recal_cmd2.=" -o $recalBamFile";
$recal_cmd2.=" -BQSR $recalCSVFile";
#print "$recal_cmd2\n";
#system($recal_cmd2);
$jname="C2.${bamFileName}";
$jout="C2.${bamFileName}.out";
my $reCalibrate_job2=new Qsub(name=>$jname,wd=>$outDirPath,outFile=>$jout,cmd=>$recal_cmd2,waitfor=>($reCalibrate_job1->{jobid}));
$reCalibrate_job2->submit();

return $reCalibrate_job2->{jobid};
}



sub usage()
{
print <<EOF;
This script takes a folder containing bam files and
1. reAligns
2. reCalibrates
3. calls variants .... all using GATK

Author : Vishal N Koparde, Ph. D.
Created : 140115
Modified : 140105
Version : 2.0

options:
--bamFileFolder		folder with all the bam files
--force			overwrite all files
--help
EOF
exit 1;
}


sub getFileSize {
	my ($fileName) = @_;
	my $fileSize = -s $fileName;
	return $fileSize;
}
