#!/usr/bin/perl 

#use strict;
#use warnings;

##### I. Set arguments ######################################################################
use Cwd 'getcwd';
$curdir = getcwd;

use File::Which;

use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
$evaluex = 10; # -e
$evaluey = 10; # (= -e)
$pcd = 0; # -D (Use 16S nucleotide substitution rate. Poisson Corrected Discance is default.)  
$savematch = "N"; 
$savetep = "N"; 
$subdirbmatch = (); # --outblastlist
$subdirtep = (); # --out2devalue
$outdir = "MPASS_Output"; # -outdir
$logfile = (); # -L
$threadnum = 1; # -t
$indistlist = "DISTANCElist"; # -d
$outmtx = "MPASSmatrix\.txt"; # -o
$help = ""; # -h
$version = ""; # -v 
GetOptions(
	'e=i' => \$evaluex,
	'D' => \$pcd,
	'outblastlist=s' => \$subdirbmatch,
	'out2devalue=s' => \$subdirtep,
	'outdir=s' => \$outdir,
	'L=s' => \$logfile,
	't=i' => \$threadnum,
	'd=s' => \$indistlist,
	'o=s' => \$outmtx,
	'h' => \$help,
	'v' => \$version,
) || die "Invalid options";

if ($help == 1) {
	print "MPASS\.pl\n";
	print "version 2022\.2\.2\n\n";
	print "Usage\: perl MPASS\.pl \[options\]\n\n";
	print "Options\:\n";
	print "\t\-d Filename for the genomic distances\(Default\: DISTANCElist\) \[String\]\n";
	print "\t\-o Filename of distance matrix\(Default\: MPASSmatrix\.txt\) \[String\]\n";
	print "\t\-e E\-value threshold for BLASTP\(Default\: 10\) \[Real\]\n";
	print "\t\-t Number of threads for metaSPAdes and BLASTP\(Default\: 1\) \[Integer\]\n\n";
	print "\t\-\-outdir Directory for the distance matrix and list, and log files\(Default\: MPASS_Output\) \[String\]\n";
	print "\t\-\-outblastlist Name of subdirectory for BLASTP data \[String\]\n";
	print "\t\-\-out2devalue Name of subdirectory for BLASTP\-based besthit pairs \[String\]\n\n";
	print "\t\-D Use 16S nucleotide substitution rate as the genomic distance\(Default\: Poisson corrected distance\)\n\n";
	print "\t\-h Help\n";
	print "\t\-v Version\n\n\n";
	print "Requirements\:\n";
	print "\tPaired-end fastq files\n\n";
	print "\t\[Softwares \(Commands and URL for installation\)\]\n";
	print "\tSeqtk \(conda install \-c bioconda seqtk\)\n";
	print "\tmetaSPAdes\.py \(conda install \-c bioconda spades\)\n";
	print "\tMetaGeneMark \(http\:\/\/exon\.gatech\.edu\/license_download\.cgi\)\n";
	print "\tBLAST \(conda install \-c bioconda blast\-legacy\)\n\n";
	print "\tSetting PATH for the MetaGeneMark is requied\.\n";
	print "\tConfiguration file of MetaGeneMark \(MetaGeneMark_v1\.mod\) is required in the current working directory\.\n\n\n";
	print "References and URL\:\n";
	print "\tseqtk\: https\:\/\/github\.com\/lh3\/seqtk\n";
	print "\tmetaSPAdes\: Nurk S\, Meleshko D\, Korobeynikov A\, Pevzner PA\.\, Genome res\.\, 2017\;27\(5\)\:824\-34\.\n";
	print "\tMetaGeneMark\: Zhu W\, Lomsadze A\, Borodovsky M\.\, Nucleic acids res\.\, 2010\;38\(12\)\:e132\.\n";
	print "\tBLAST\: Altschul SF\, Gish W\, Miller W\, Myers EW\, Lipman DJ\.\, J mol\. biol\.\, 1990\;215\(3\)\:403\-10\.\n\n\n";
	exit;
}
if ($version == 1) {
	print "MPASS\.pl\nversion 2022\.2\.2\n\n";
	exit;
}

##### II. Checking errors ###################################################################
##### A. Softwares and config files #########################################################
$checktools = 1;
$pathlen = 0;

$toolpath = which("seqtk");
$pathlen = length($toolpath);
$checktools = $checktools * $pathlen;

$toolpath = which("metaspades\.py");
$pathlen = length($toolpath);
$checktools = $checktools * $pathlen;

$toolpath = which("gmhmmp"); # MetaGeneMark
$pathlen = length($toolpath);
$checktools = $checktools * $pathlen;

$toolpath = which("blastall");
$pathlen = length($toolpath);
$checktools = $checktools * $pathlen;

$checkgmhmmpconfig = 0;
opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;
foreach $file (@dir) {
	chomp $file;
	if ($file eq "MetaGeneMark_v1\.mod") {
		$checkgmhmmpconfig = 1;
	}
}
$checktools = $checktools * $checkgmhmmpconfig;

if ($checktools == 0) {
	print "Error\: Softwares or config files are not installed\.\n";
	$help = 1;
}

@dir = ();
$file = "";

##### B. Specification of matrix file name ##################################################
if (length($outmtx) == 0) {
	$help = 1;
	print "Error\: Matrix file name is not specified\(\-o option\)\.\n";
}

##### C. Duplication and number of filenames ################################################
opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

$count = 0;
$filedup = 0;
%filecheckhash = ();
$first8char = "";
foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_[12]\.fastq$/) {
		++$count;
		$first8char = $file;
		$first8char =~ s/_[12]\.fastq$//;
		if (length($first8char > 8)) {
			$first8char = substr($first8char, 0, 8);
		}
		if(exists($filecheckhash{$first8char})) {
			$filecheckhash{$first8char} = $filecheckhash{$first8char} + 1;
			if ($filecheckhash{$first8char} > 2) {
				$filedup = 1;
			}
		} else {
			$filecheckhash{$first8char} = 1;
		}
	}
}
if($count < 6) {
	print "Error\: The number of sequence datasets is less than 3\.\n";
	$help = 1;
}
if($filedup == 1) {
	print "Error\: Duplication of filenames\(first 8 characters\)\.\n";
	$help = 1;
}
@dir = ();
$file = "";
$count = "";

##### D. Already created files (matrix, distance list, log file) ############################
# To avoid overwriting, the matrix, the distance list, and log files are created in a different folder (not the current directory).
opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file eq $outmtx) {
		print "Error\: $file aleardy exists\.\n";
		$help = 1;
	}
	if ($file eq $indistlist) {
		print "Error\: $file aleardy exists\.\n";
		$help = 1;
	}
	if ($file eq "MPASS_All_log") {
		print "Error\: $file aleardy exists\.\n";
		$help = 1;
	}
	if ($help == 1) {
		last;
	}
}
@dir = ();
$file = "";

##### E. If errors occur, show help and force termination ################################### 
if ($help == 1) {
	print "MPASS\.pl\n";
	print "version 2022\.2\.2\n\n";
	print "Usage\: perl MPASS\.pl \[options\]\n\n";
	print "Options\:\n";
	print "\t\-d Filename for the genomic distances\(Default\: DISTANCElist\) \[String\]\n";
	print "\t\-o Filename of distance matrix\(Default\: MPASSmatrix\.txt\) \[String\]\n";
	print "\t\-e E\-value threshold for BLASTP\(Default\: 10\) \[Real\]\n";
	print "\t\-t Number of threads for metaSPAdes and BLASTP\(Default\: 1\) \[Integer\]\n\n";
	print "\t\-\-outdir Directory for the distance matrix and list, and log files\(Default\: MPASS_Output\) \[String\]\n";
	print "\t\-\-outblastlist Name of subdirectory for BLASTP data \[String\]\n";
	print "\t\-\-out2devalue Name of subdirectory for BLASTP\-based besthit pairs \[String\]\n\n";
	print "\t\-D Use 16S nucleotide substitution rate as the genomic distance\(Default\: Poisson corrected distance\)\n\n";
	print "\t\-h Help\n";
	print "\t\-v Version\n\n\n";
	print "Requirements\:\n";
	print "\tPaired-end fastq files\n\n";
	print "\t\[Softwares \(Commands and URL for installation\)\]\n";
	print "\tSeqtk \(conda install \-c bioconda seqtk\)\n";
	print "\tmetaSPAdes\.py \(conda install \-c bioconda spades\)\n";
	print "\tMetaGeneMark \(http\:\/\/exon\.gatech\.edu\/license_download\.cgi\)\n";
	print "\tBLAST \(conda install \-c bioconda blast\-legacy\)\n\n";
	print "\tSetting PATH for the MetaGeneMark is requied\.\n";
	print "\tConfiguration file of MetaGeneMark \(MetaGeneMark_v1\.mod\) is required in the current working directory\.\n\n\n";
	print "References and URL\:\n";
	print "\tseqtk\: https\:\/\/github\.com\/lh3\/seqtk\n";
	print "\tmetaSPAdes\: Nurk S\, Meleshko D\, Korobeynikov A\, Pevzner PA\.\, Genome res\.\, 2017\;27\(5\)\:824\-34\.\n";
	print "\tMetaGeneMark\: Zhu W\, Lomsadze A\, Borodovsky M\.\, Nucleic acids res\.\, 2010\;38\(12\)\:e132\.\n";
	print "\tBLAST\: Altschul SF\, Gish W\, Miller W\, Myers EW\, Lipman DJ\.\, J mol\. biol\.\, 1990\;215\(3\)\:403\-10\.\n\n\n";
	exit;
}

##### III. BLAST setting and Make files and subdirectries ###################################
unless ($evaluex == 10) {
	$evaluey = $evaluex;
}
if ($subdirbmatch =~ /[a-zA-Z0-9]/) {
	$savematch = "Y";
}
if ($subdirtep =~ /[a-zA-Z0-9]/) {
	$savetep = "Y";
}

system("mkdir -m 777 MPASStmp_GPD"); #directory for the (down-sampled) fastq and fasta files.
system("cp MetaGeneMark_v1\.mod $curdir\/MPASStmp_GPD\/MetaGeneMark_v1\.mod"); #MetaGeneMark_v1.mod is required for MetaGeneMark.
system("mkdir -m 777 MPASStmp_CDM"); #directory to construction of distance matrix for the metaphylogenomic tree.
chdir("MPASStmp_CDM");
if ($savematch eq "Y") {
	system("mkdir -m 777 $subdirbmatch"); #directory for the BLAST best-matched pairs
}
if ($savetep eq "Y") {
	system("mkdir -m 777 $subdirtep"); #directory for the e-values of query-query and query-reference pairs.
}
chdir("\.\.");
if ($outdir =~ /[a-zA-Z0-9]/) {
	system("mkdir -m 777 $outdir"); #directory for the distance matrix, distance list and log files.
}

if ($logfile =~ /[a-zA-Z0-9]/) { #When the analysis is interrupted, analysis can be resumed from the MPASS-log file.
	open (CHECKLOG, $logfile);
	@checklog = <CHECKLOG>;
} else {
	@checklog = ();
}

$logcount = 0; #For checking each step whether such steps were already recorded in the log file.
open OUTLOG, ">phyl_proteome_log"; #Open new log file.

###### PATH for BLAST ####################################################################### 
$afexecute = "formatdb"; # Path for the formatdb.
$blast_execute = "blastall"; # Path for the Legacy-BLAST.


##### IV. Main ##############################################################################
##### 1. Down-sampling (seqtk) ##############################################################
# Requirement: paired-fastq files (_1.fastq, _2.fastq)
# First 8 characters of fastq-files (without "_1.fastq", "_2.fastq") are used for the OTU name in the distance matrix for metaphylogenomic tree construction.
print "\nSTEP1\: Down\-sampling\(seqtk\)\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

@dir2 = @dir;
$minseqnum = "";
$sw = 0;
foreach $file (@dir2) {
	chomp $file;
	if ($file =~ /_1\.fastq$/) {
		open (IN, $file);
		@lines = <IN>;
		$seqnum = @lines;
		$seqnum = $seqnum / 4;
		if ($sw == 0) {
			$minseqnum = $seqnum;
			$sw = 1;
		} else {
			if ($minseqnum > $seqnum) {
				$minseqnum = $seqnum;
			}
		}
		close IN;
	}
}

print "\tNumber of reads in each fastq file\: $minseqnum\n";

@dir2 = @dir;
foreach $file (@dir2) { #Down-sampling and file name shortening
	chomp $file;
	if ($file =~ /_1\.fastq$/) {
		$file2 = $file;
		$file2 =~ s/_1\.fastq$/_2\.fastq/;
		$seqtkfile1 = $file;
		if(length($seqtkfile1) > 16) {
			$seqtkfile1 = substr($seqtkfile1, 0, 8);
			$seqtkfile1 = $seqtkfile1 . "_1\.fastq";
		}
		$seqtkfile1 =~ s/_1\.fastq$/_MPASSseqtk_1\.fastq/;		
		$seqtkfile2 = $seqtkfile1;
		$seqtkfile2 =~ s/_MPASSseqtk_1\.fastq$/_MPASSseqtk_2\.fastq/;
		$sysin = "seqtk sample \-s 1 $file $minseqnum \> $seqtkfile1";
		print "\t\t$sysin\n";
		system($sysin);
		$sysin = "mv $seqtkfile1 MPASStmp_GPD";
		system($sysin);
		$sysin = "seqtk sample \-s 1 $file2 $minseqnum \> $seqtkfile2";
		print "\t\t$sysin\n";
		system($sysin);
		$sysin = "mv $seqtkfile2 MPASStmp_GPD";
		system($sysin);
	}
}
@dir = ();
@dir2 = ();
$file = "";
$file2 = "";
$sysin = "";
$sw = 0;

##### 2. Assembling (metaSPAdes) ############################################################
#Input: _MPASSseqtk_1.fastq, _MPASSseqtk_2.fastq
#Output: _contigs.fasta
print "STEP2\: Assembling\(metaSPAdes\)\n";
chdir("MPASStmp_GPD");

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_MPASSseqtk_1\.fastq$/) {
		$file2 = $file;
		$file2 =~ s/_MPASSseqtk_1\.fastq$/_MPASSseqtk_2\.fastq/;
		$spadesout = $file;
		$spadesout =~ s/_MPASSseqtk_1\.fastq$/_contigs\.fasta/;
		$sysin = "metaspades\.py \-t $threadnum \-1 $file \-2 $file2 \-o MPASStmp_mSPA";
		print "\t$sysin\n";
		$sysin = "$sysin \> \/dev\/null 2\>\&1";
		system($sysin);
		chdir("MPASStmp_mSPA");
		$sysin = "cp contigs\.fasta $curdir\/MPASStmp_GPD\/$spadesout";
		system($sysin);
		chdir("\.\.");
		$sysin = "rm \-rf MPASStmp_mSPA";
		system($sysin);
		unlink($file);
		unlink($file2);
	}
}
@dir = ();
$file = "";
$file2 = "";
$sysin = "";

##### 3. Gene prediction (MetaGeneMark) #####################################################
#Input: _contigs.fasta
#Output: _MGM_aa.fa, _MGM_nt.fa
print "STEP3\: Gene prediction\(MetaGeneMark\)\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

$mgmseqlist = "MPASStmp_MGM_list\.txt";

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_contigs\.fasta$/) {
		$profasta = $file;
		$profasta =~ s/_contigs\.fasta$/_MGM_aa\.fa/; 
		$ntfasta = $file;
		$ntfasta =~ s/_contigs\.fasta$/_MGM_nt\.fa/;
		$sysin = "gmhmmp \-a \-d \-f L \-m MetaGeneMark_v1\.mod $file \-o $mgmseqlist \-A $profasta \-D $ntfasta";
		print "\t$sysin\n";
		system($sysin);
		#unlink($file);
		unlink($mgmseqlist);
	}
}
@dir = ();
$file = "";
$sysin = "";

##### 4. Remove partial and short sequences (Custom Perl script) ############################
#4-1. Make DNA and Protein sequence lists ###################################################
#Input: _contigs.fasta
#Output: _aa (FASTA files without file extension)
#GeneID<TAB>Description for protein sequences<TAB>Description for scaffolds<TAB>Sequence

print "STEP4\: Remove partial and short sequences\(custom script\)\n";
#print "\tMake DNA and Protein sequence lists\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_MGM_aa\.fa$/) {
		$outrpss = $file;
		$outrpss =~ s/_MGM_aa\.fa$/_MGM_aa_seqlist\.txt/;
	} elsif ($file =~ /_MGM_nt\.fa$/) {
		$outrpss = $file;
		$outrpss =~ s/_MGM_nt\.fa$/_MGM_nt_seqlist\.txt/;
	}
	if ($file =~ /_MGM_[an][at]\.fa$/) {
		open (IN, $file);
		open OUTRPSS, ">$outrpss";

		$sw = 0;
		while ($line = <IN>) {
			chomp $line;
			if ($line =~ /[a-zA-Z]/) {
				if ($line =~ /^\>/) {
					@data = split(/\t/, $line);
					$gene = $data[0];
					$gene =~ s/^\>//;
					$gene =~ s/\|.*//g;
					if ($sw == 0) {
						print OUTRPSS "$gene\t$data[0]\t$data[1]\t";
						$sw = 1;
					} else {
						print OUTRPSS "\n$gene\t$data[0]\t$data[1]\t";
					}
				} else {
					print OUTRPSS $line;
				}
			}
		}
		print OUTRPSS "\n";
		close OUTRPSS;
		close IN;
		unlink($file);
	}
}
@dir = ();
$file = "";
$outrpss = "";
$sw = 0;
$gene = "";

##### 4-2. Extract nucleotide sequences containing initiation and termination codons ########

#print "\tExtract complete CDSs from _nt_seqlist\.txt\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_MGM_nt_seqlist\.txt$/) {
		open (IN, $file);
		$outrpss = $file;
		$outrpss =~ s/_MGM_nt_seqlist\.txt$/_nt_sl_compCDS\.txt/;
		open OUTRPSS, ">$outrpss";
				
		$sw1 = 0;
		$sw2 = 0;
		while ($line = <IN>) {
			chomp $line;
			@data = split(/\t/, $line);
			if ($data[3] =~ /^ATG.*/) { #ATG, GTG, ATA and TTG are defined as initiation codons.
				$sw1 = 1;
			} elsif ($data[3] =~ /^GTG.*/) {
				$sw1 = 1;
			} elsif ($data[3] =~ /^ATA.*/) {
				$sw1 = 1;
			} elsif ($data[3] =~ /^TTG.*/) {
				$sw1 = 1;
			}
			if ($data[3] =~ /TAA$/) {
				$sw2 = 1;
			} elsif ($data[3] =~ /TAG$/) {
				$sw2 = 1;
			} elsif ($data[3] =~ /TGA$/) {
				$sw2 = 1;
			}
			if ($sw1 == 1) {
				if ($sw2 == 1) {
					print OUTRPSS "$line\n";
				}
			}
			$sw1 = 0;
			$sw2 = 0;
		}
		close IN;
		close OUTRPSS;
		unlink($file);
	}
}
@dir = ();
$file = "";
$outrpss = "";
$sw1 = 0;
$sw2 = 0;

##### 4-3. Extract full length amino acid sequences #########################################

#print "\tExtract Complete CDS from _aa_seqlist\.txt\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_MGM_aa_seqlist\.txt$/) { # IN2
		$file2 = $file;
		$file2 =~ s/_MGM_aa_seqlist\.txt$/_nt_sl_compCDS\.txt/; # IN1
		open(IN1, $file2);
		open(IN2, $file);
		$outrpss = $file;
		$outrpss =~ s/_MGM_aa_seqlist\.txt$/_aa_sl_compCDS\.txt/;
		open OUTRPSS, ">$outrpss";
		%aalist = ();
		while($line = <IN2>) { # _MGM_aa_seqlist.txt
			chomp $line;
			if($line =~ /[a-zA-Z0-9]/) {
				@data = split(/\t/, $line);
				$aalist{$data[0]} = $line;
			}
		}
		while($line = <IN1>) { # _nt_sl_compCDS.txt
			chomp $line;
			@data = split(/\t/, $line);
			if(exists($aalist{$data[0]})) {
				print OUTRPSS "$aalist{$data[0]}\n";
			}
		}
		close IN1;
		close IN2;
		close OUTRPSS;
		unlink($file);
		unlink($file2);
	}
}
@dir = ();
$file = "";
$file2 = "";
$outrpss = "";

##### 4-4. Remove protein sequences shorter than 100aa ######################################

#print "\tRemove short protein sequences \(\< 100 amino acids\)\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
    chomp $file;
    if ($file =~ /_aa_sl_compCDS\.txt$/) {
        open (IN, $file);
        $outrpss = $file;
        $outrpss =~ s/_aa_sl_compCDS\.txt$/_aa_sl_compCDS100\.txt/;
        open OUTRPSS, ">$outrpss";

        while ($line = <IN>) {
            chomp $line;
            if ($line =~ /[a-zA-Z]/) {
                @data = split(/\t/, $line);
                if (length($data[3]) > 99) {
                    print OUTRPSS "$line\n";
                }
            }
        }
        close IN;
        close OUTRPSS;
        unlink($file);
    }
}
@dir = ();
$file = "";
$outrpss = "";

##### 4-5. Make protein fasta file from sequence list #######################################

#print "\tMake protein fasta files \(Full\-length\, \>99 amino acids\)\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
	chomp $file;
	if ($file =~ /_aa_sl_compCDS100\.txt$/) {
		$outrpss = $file;
		$outrpss =~ s/_aa_sl_compCDS100\.txt$/_aa_sl_compCDS_fasta\.txt/;
		open (IN, $file);
		open OUTRPSS, ">$outrpss";

		while ($line = <IN>) {
			chomp $line;
			if ($line =~ /[a-zA-Z]/) {
				@data = split(/\t/, $line);
				$gene = $data[2];
				$gene =~ s/^\>//;
				$gene = "\>" . $data[0] . "\_" . $gene;
				print OUTRPSS "$gene\n$data[3]\n"; # E.g. >gene1_NODE_1_length_aaaa_cov_bbbb
			}
		}
		close IN;
		close OUTRPSS;
		unlink($file);
	}
}
@dir = ();
$file = "";
$outrpss = "";
$gene = "";

##### 4-6. Change FASTA headers #############################################################
#Requirement: Latest version of metaSPAdes.py (22' 4/3, bioconda)

#print "\tChange FASTA headers\n";
print "\tFollowing FASTA files \(_aa\) are created\.\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

foreach $file (@dir) {
    chomp $file;
	if ($file =~ /_aa_sl_compCDS_fasta\.txt$/) {
        $outrpss = $file;
		$outrpss =~ s/_aa_sl_compCDS_fasta\.txt$/_aa/;
        open (IN, $file);
        open OUTRPSS, ">$outrpss";

        while ($line = <IN>) {
            chomp $line;
            if ($line =~ /^\>/) { # >gene_t_NODE_m_length_aaaa_cov_bbbb -> >gene_t_NODE_m_covbbbb_lengthaaaa
				@data = split(/_/, $line);
				$gene = $data[0] . "_" . $data[1] . "_" . $data[2] . "_" . $data[3] . "_" . $data[6] . $data[7] . "_" . $data[4] . $data[5] . "\n";
                print OUTRPSS $gene;
            } else {
                print OUTRPSS "$line\n";
            }
        }
		close IN;
        close OUTRPSS;
		print "\t\t$outrpss\n";
		unlink($file);
    }
}
@dir = ();
$file = "";
$outrpss = "";
$gene = "";


##### 5. Construction of distance matrix for the metaphylogenomic tree  #####################
print "STEP5\: Construction of distance matrix\(BLASTP and custom scripts\)\n";

##### 5-1. Move current directory and _aa files to $outdir ##################################
$moveto = $curdir . "\/MPASStmp_CDM";
$sysin = "mv \*_aa $moveto";
system($sysin);
chdir($moveto);

##### 5-2. Calculation of metagenomic distances  ############################################
##### 5-2-1. Formatting the BLAST Databases ##################################################
print "\tFormatting the BLAST Databases\n";

opendir DIR, ".";
@dir = readdir DIR;
closedir DIR;

@dir1 = @dir;
foreach $db_file (@dir1) {
	if ($db_file =~ /_aa$/) {
		system ($afexecute . " \-i " . $db_file . " \-p T");
	}
}

##### 5-2-2. BLAST ##########################################################################
#The end of the filenames with "_aa" is used as the query or reference metagenomes. 
print "\tBLASTP \(\.out files\)\n";
print "\t\tE\-value\: $evaluex\, Threads\: $threadnum\n";
@dir1 = @dir;
@dir2 = @dir;
foreach $query_file (@dir1) {

	if ($query_file =~ /_aa$/ ) {
		print "\t\t\tQuery DB\: $query_file\n";
		foreach $db_file (@dir2) {
			if ($db_file =~ /_aa$/) {
				
				foreach $filecheck (@checklog) {
					chomp $filecheck;
					if ($filecheck eq "$query_file\_and\_$db_file\.out") {
						$logcount = 1;
						last;
					}
				}
				if ($logcount =~ 0) {
					#print "\t\t\tInput Reference DB \t$db_file\n";
					if ($query_file eq $db_file) {
						#print "blastall \-p blastp \-i $query_file \-d $db_file \-e $evaluex \-m8 \-o $query_file_and_$db_file\.out \-a $threadnum \-b 1\n";
						system ($blast_execute . " -p" . " blastp" . " -i " . $query_file . " -d " . $db_file . " -e " . $evaluex . " -m8" . " -o " . "$query_file\_and\_$db_file\.out" . " -a " . $threadnum . " -b 1 \> \/dev\/null 2\>\&1");
					} else {
						#print "blastall \-p blastp \-i $query_file \-d $db_file \-e $evaluey \-m8 \-o $query_file_and_$db_file\.out \-a $threadnum \-b 1\n";					}
						system ($blast_execute . " -p" . " blastp" . " -i " . $query_file . " -d " . $db_file . " -e " . $evaluey . " -m8" . " -o " . "$query_file\_and\_$db_file\.out" . " -a " . $threadnum . " -b 1 \> \/dev\/null 2\>\&1");
					}
				} else {
					$logcount = 0;
				}
				print OUTLOG "$query_file\_and\_$db_file\.out\n";
			}
		}
	}
}
closedir DIR;

##### 5-2-3. Output of Best-matched pairs ###################################################
#Open files with the ".out" extension and output the result as ".tpht" files.
opendir DIR, ".";
@dir = readdir DIR;

print "\t\tOutput best\-matched pairs \(\.tpht files\)\n";

$queryname2 = 0;
foreach $blast_result (@dir) {
	if ($blast_result =~ /\.out$/) {
		#print "\t\tBest\-Hit\:Input\t$blast_result\n";
		$outname = $blast_result;
		$outname =~ s/\.out$/\.tpht/;
		
		$logcount = 0;
		foreach $filecheck (@checklog) {
			chomp $filecheck;
			if ($filecheck eq $outname) {
				$logcount = 1;
				last;
			}
		}
		if ($logcount == 0) {
			open (BLASTRES, $blast_result);
			$output = ();
	
			while ($line = <BLASTRES>) {
				chomp $line;
				$queryname = $line;
				$queryname =~ s/\s.*//g;
				unless ($queryname2 eq $queryname) {
					$output = $output . $line . "\n";
				}	
				$queryname2 = $queryname;
			}
			close BLASTRES;
			open OUT, ">$outname";
			print OUT $output;
			close OUT;
			#print "\t\t\tOutput\t$outname\n";
			
			print OUTLOG "$outname\n";			
			unlink ($blast_result);
		}
	} elsif ($blast_result =~ /\_aa\.phr$/) {
		unlink ($blast_result); 
	} elsif ($blast_result =~ /\_aa\.pin$/) {
		unlink ($blast_result);
	} elsif ($blast_result =~ /\_aa\.psq$/) {
		unlink ($blast_result);
	}
}
closedir DIR;

##### 5-2-4. Concatenation of query and ref BLAST results ###################################
#Concatenate the query and reference ".tpht" files and output them as ".asm" files.

print "\t\tConcatenate query and reference BLAST results \(\.asm files\)\n";

opendir DIR, ".";
@dir = readdir DIR;
@dir1 = @dir;
@dir2 = @dir;
foreach $firstfile (@dir1) {
	if ($firstfile =~ /\.tpht$/ ) {
		$ffname = $firstfile;
		chomp $ffname;
		$ffname =~ s/\.tpht$//;
		@swords = split (/\_and\_/, $ffname);
		if ($swords[0] eq $swords[1]){
			#print "\t\tInput\(1st\)\: $firstfile\n";
			foreach $secondfile (@dir2) {
				if ($secondfile =~ /\.tpht$/) {
					$sfname = $secondfile;
					chomp $sfname;
					$sfname =~ s/\.tpht//;
					@swords = split (/\_and\_/, $sfname);
					unless ($swords[0] eq $swords[1]) {
						#print "\t\t\tInput\(2nd\)\: $secondfile\n";
						if ($firstfile eq $secondfile) {
							next;
						} 
						$checkname1 = $firstfile;
						$checkname2 = $secondfile;
						$checkname1 =~ s/aa\_.*//g;
						$checkname2 =~ s/aa\_.*//g;
						if ($checkname1 eq $checkname2) {
							$outF = $firstfile;
							$outR = $secondfile;
							$outF =~ s/\.tpht//g;
							$outR =~ s/\.tpht//g;
							$outfile = "1" . $outF . "\_2" . $outR . "\.asm";
							
							$logcount = 0;
							foreach $filecheck (@checklog) {
								chomp $filecheck;
								if ($filecheck eq $outfile) {
									$logcount = 1;
									last;
								}
							}
							if ($logcount == 0) {
								$output = "";
								open (FIRSTFILE, $firstfile);
								open (SECONDFILE, $secondfile);

								%hash2 = ();
								while ($list2 = <SECONDFILE>) {
									chomp $list2;
									@element2 = split(/\t/, $list2);
									$hash2{$element2[0]} = $list2;
									@element2 = ();
									$list2 = "";
								}
								while ($list1 = <FIRSTFILE>) {
									chomp $list1;
									@element1 = split (/\t/, $list1);
									if (exists($hash2{$element1[0]})) { #query \t reference1 \t reference2 \t e-value1 \t evalue2
										@element2 = split(/\t/, $hash2{$element1[0]});
										$output = "$output$element1[0]\t$element1[1]\t$element2[1]\t$element1[10]\t$element2[10]\n";
									}
								}
								close FIRSTFILE;
								close SECONDFILE;
								open OUT, ">$outfile";
								print OUT $output;
								close OUT;
							}
						
							#print "\t\t\t\tOutput\: $outfile\n";
							print OUTLOG "$outfile\n";
						}
					}
				}
			}
		}
	}
}
closedir DIR;	

opendir CHECKDIR, ".";
@checkdir = readdir CHECKDIR;
foreach $checktpht (@checkdir) {
	if ($checktpht =~ /\.tpht$/) {
		if ($savematch eq "Y") {
			system ("mv $checktpht $subdirbmatch");
		} elsif ($savematch eq "N") {
			unlink ($checktpht);
		}
	}
}
closedir CHECKDIR; 

##### 5-2-5. Estimation of metagenomic distances (without conversion to rRNA-based distances) ###
print "\tEstimate the metagenomic distances\n";

open OUT, ">$indistlist";
print OUT "FILENAME\tDistance\(not based on the rRNA\)\n";

opendir DIR, ".";
@dir = readdir DIR;

@distdata = ();
foreach $in_filename (@dir) {
	if ($in_filename =~ /\.asm$/) {
		open (INFILE, $in_filename);
		#print "\t\t$in_filename\n";

		$line_num = 0;
		$evaluea1 = 0;
		$evaluea2 = 0;
		$evalueb1 = 0;
		$evalueb2 = 0;
		$result = 0;
		$tan1 = 0;
		$tan2 = 0;
		$tan3 = 0;
				
		while ($line = <INFILE>) {
			@element = split (/\t/, $line);
			
			#Coverage of query genes ##
			#metaspades: [0]: "gene", [1]: GeneNo., [2]: "Node", [3]: NodeNo., [4]: Coverage (E.g. cov****), [5]: Length
			$queryinfo = $element[0];						
			@queryelement = split(/\_/, $queryinfo);		
			$cov = $queryelement[4];						
			$cov =~ s/^cov//;								
			$queryinfo = "";								
			@queryelement = ();								
		
		
			$evaluea1 = $element[3];
			$evalueb1 = $element[4];
			$evaluea1 =~ s/\s//g;
			$evalueb1 =~ s/\s//g;

		#Convert e-values to logarithm.
		#If e-values are "0", log "e-value" is set to "-180". (for NCBI-BLAST)
			if ($evaluea1 =~ /e/) {
				$logFa = $evaluea1;
				$logFa =~ s/e.*//;
				$logRa = $evaluea1;
				$logRa =~ s/.*e//;
				$logFa = log $logFa;
				$log10 = log 10;
				$logFa = $logFa / $log10;
				$evaluea1 = $logFa + $logRa;
			} elsif ($evaluea1 == 0) {
				$evaluea1 = -180;
			} else {
				$logFa = $evaluea1;	  
				$logFa = log $logFa;
				$log10 = log 10;
				$logFa = $logFa / $log10;
				$evaluea1 = $logFa;
			}
			if ($evalueb1 =~ /e/) {
				$logFb = $evalueb1;
				$logFb =~ s/e.*//;
				$logRb = $evalueb1;
				$logRb =~ s/.*e//;
				$logFb = log $logFb;
				$log10 = log 10;
				$logFb = $logFb / $log10;
				$evalueb1 = $logFb + $logRb;
			} elsif ($evalueb1 == 0) {
				$evalueb1 = -180;
			} else {
				$logFb = $evalueb1;	  
				$logFb = log $logFb;
				$log10 = log 10;
				$logFb = $logFb / $log10;
				$evalueb1 = $logFb;
			}
			
			#Baseline correction
			$log10 = log 10;
			$shiftx = log $evaluex;
			$shiftx = $shiftx / $log10;
			$shifty = log $evaluey;
			$shifty = $shifty / $log10;
			$evaluea1 = $evaluea1 - $shiftx;
			$evalueb1 = $evalueb1 - $shifty;
			
			#E-value x Coverage
			$evaluea1 = $evaluea1 * $cov;
			$evalueb1 = $evalueb1 * $cov;
			
			$evaluea2 = $evaluea1 + $evaluea2;
			$evalueb2 = $evalueb1 + $evalueb2;
					
			$line_num = $line_num + $cov;
			$cov = "";
		}
		
		$evalue_avga = $evaluea2 / $line_num; #Averaging the query-query e-values
		$evalue_avgb = $evalueb2 / $line_num; #Averaging the query-reference e-values
				
		#query-reference e-value / query-query e-value
		$result = $evalue_avgb / $evalue_avga;
		
		#Estimation of genomic distance
		#(m-1)/(m+1) x -1, (m = query-ref e-value / query-query e-value)
		$tan1 = $result - 1;
		$tan2 = $result + 1;
		$tan3 = $tan1 / $tan2;
		$result = $tan3 * -1;

		$outdata = "$in_filename\t$result\n";
		push (@distdata, $outdata);
		close INFILE;
		@element = ();

		if ($savetep eq "Y") {
			system ("mv $in_filename $subdirtep");
		} elsif ($savetep eq "N") {
			unlink ($in_filename);
		}
	}
}
print OUT @distdata;
close OUT;
closedir DIR;

print OUTLOG "Metagenomic distance list\: $indistlist\n";
print "\t\tMetagenomic distance list\: $indistlist\n";

##### 5-3. Make a distance matrix  ######################################################
print "\tConvert rRNA\-based distances and Make a distance matrix\n";

open (IN, $indistlist);
chomp $outmtx;
open OUT, ">$outmtx";

$e = 2.71828182845904;

# rRNA-based genomic distance ######
# rRNA_based_distance = 4.142 x e ^ (2.824 x genomic_distance)  (Poisson corrected distance)
# rRNA_based_distance = 5.001 x e ^ (2.324 x genomic_distance)	(Original substitution rate)
if ($pcd == 0) { #Poisson corrected distance
	$const1 = 4.142;
	$const2 = 2.824;
} elsif($pcd == 1) { #Original 16S substitution rate
	$const1 = 5.001;
	$const2 = 2.324;
}

$seqcount = ""; #Number of query (For PHYLIP distance matrix construction)
@qlist = ();
%qhash = ();
%hash = ();

#Estimation of rRNA-based genomic distances
while ($line = <IN>) {
	chomp $line;
	if ($line =~ /[0-9]/) {
		$line =~ s/^1//;
		$line =~ s/\.asm//;
		$line =~ s/_aa//g;
		$line =~ s/_and_.*_and_/\t/;
		@data = split(/\t/, $line);

		### rRNA-based genomic distance ######
		# rRNA_based_distance = Const1 x e ^ (Const2 x genomic_distance)
		$dist = $data[2]; #genomic distance
		$dist = $const2 * $dist;
		$dist = $e ** $dist;
		$dist = $const1 * $dist;
		$qandr = $data[0] . "\t" . $data[1];
		$hash{$qandr} = $dist;

		### Count query genomes
		unless (exists($qhash{$data[0]})) {
			$qhash{$data[0]} = 1;
		}

		@data = ();
		$qandr = "";
		$dist = "";
	}
}

$seqcount = keys(%qhash);
@qlist = keys(%qhash);
@qlist = sort {$a cmp $b} @qlist;


# Averaging the rRNA-based distances (Distances from reciprocal pairs of 2 metagenomes)
foreach $key (keys(%hash)) {
	unless ($hash{$key} =~ /\tavg$/) { 
		$qpair = $key;

		$revpair = $key;
		@data = split(/\t/, $revpair);
		$revpair = $data[1] . "\t" . $data[0];

		# Averaging the Distances
		$avg = ($hash{$qpair} + $hash{$revpair}) / 2;
		$avg = $avg . "\tavg";
		$hash{$qpair} = $avg;
		$hash{$revpair} = $avg;

		$qpair = "";
		$revpair = "";
		@data = ();
		$avg = "";
	}

	# If both query and reference are identical, their distance is "0.000000000000".
	$pair = $key;
	@data = split(/\t/, $pair);
	$qandq = $data[0] . "\t" . $data[0];
	unless (exists($hash{$qandq})) {
		$hash{$qandq} = "0\.000000000000\tavg";
	}
	
	$pair = "";
	@data = ();
	$qandq = "";
	$key = "";
}

# Make distance matrix
@qlist2 = @qlist;
print OUT "   $seqcount\n"; # Three spaces in the front of $seqcount for PHYLIP-style distance matrix.

foreach $query (@qlist) {
	# For adjusting the number of characters of query names
	if (length($query) < 8) {
		$j = 8 - length($query);
		for ($i = 0; $i < $j; ++$i) {
			$query = $query . " ";
		}
		$j = "";
	} else {
		$query = substr($query, 0, 8); # Initial 8 characters from longer name are used
	}
	print OUT "$query  "; 

	$query =~ s/\s+$//; 
	@qlist2_2 = @qlist2;
	foreach $ref (@qlist2_2) {
		unless (length($ref) < 8) {
			$ref = substr($ref, 0, 8);
		}
		$ref =~ s/\s+$//;

		if ($query eq $ref) { # If both query and refrerence are identical, print out "0.000000".
			print OUT "  0\.000000";
		} else {
			$qandr = $query . "\t" . $ref;
			$dist = $hash{$qandr};
			$dist =~ s/\tavg$//;

			# For adjusting the number of digits
			if ($dist =~ /\./) {
				$dist = $dist . "000000";
			} else {
				$dist = $dist . "\.000000";
			}
			$dist = substr($dist, 0, 8);
			print OUT "  $dist";

			$qandr = "";
			$dist = "";
		}
	}
	print OUT "\n";
	@qlist2_2 = ();
	$query = "";
}

close IN;
close OUT;
print "\t\tDistance matrix for PHYLIP\: $outmtx\n\n";
print OUTLOG "Distance matrix for PHYLIP\: $outmtx\n";
print OUTLOG "Distance matrix is created\!\n\n";
print OUTLOG "MPASS\.pl is completed\!\n";
close OUTLOG;

$outdir = $curdir . "\/" . $outdir;
opendir DIR, "."; #Move distance matrix, list of distances(not rRNA-based distances) and log files to the MPASS_Output directory.
@dir = readdir DIR;
foreach $file (@dir) {
	chomp $file;
	if ($file eq $outmtx) {
		system ("mv $file $outdir");
	} elsif ($file eq $indistlist) {
		system ("mv $file $outdir");
	} elsif ($file eq "formatdb\.log") {
		system ("mv $file $outdir");
	} elsif ($file eq "error\.log") {
		system ("mv $file $outdir");
	} elsif ($file eq "MPASS_All_log") {
		system ("mv $file $outdir");
	}
}
closedir DIR;
chdir("\.\.");

print "MPASS\.pl is completed\.\n\n";

exit;