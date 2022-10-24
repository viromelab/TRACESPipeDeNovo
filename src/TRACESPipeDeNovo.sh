#!/bin/bash
#
STUDY="data";
MIN_READ_LENGTH="30";
MIN_NC="0.9";
MIN_CONTIG_LENGTH="100";
THREADS="8";
FW_READS="";
RV_READS="";
ADAPTERS="adapters.fa";
DATABASE="VDB.fa";
REPORTS_DIR="reports";
RESULTS_DIR="results";
HUMAN_READS="0";
HUMAN_FASTA="chm13v2.0.fa.gz";
#
FALCON_PARAM_DB=" -n $THREADS -m 6:1:1:0/0 -m 13:50:1:0/0 -m 19:500:1:5/10 -g 0.85 -c 30 ";
FALCON_PARAM_CONTIGS=" -n $THREADS -m 11:1:1:0/0 -m 13:50:1:0/0 -g 0.85 -c 10 ";
NC_PARAM=" -t $THREADS -m 3:1:1:0:0:0.85/0:0:0 -m 12:20:1:1:0:0.9/0:0:0 ";
GECO3_PARAM=" -v -rm 6:1:1:0:0.8/0:0:0 -rm 13:50:1:0:0.9/0:0:0 -rm 19:500:1:10:0.9/5:10:0,9 ";
#
# =============================================================================
# FUNCTIONS -------------------------------------------------------------------
#
CHECK_FILE () {
  if [ ! -f $1 ];
    then
    echo -e "\e[31mERROR: $1 file not found!\e[0m"
    echo "For addition information, see the instructions at the web page."
    exit 1;
    fi
  }
#
CHECK_FILE_GZIPED () {
  if (file $1 | grep -q compressed ) ; then
    echo "$1 is compressed";
  else
    echo -e "\e[31mERROR: $1 file is not compressed with Gzip!\e[0m"
    echo "Before running, use: gzip $1";
    exit;
  fi
  }
#
PROGRAM_EXISTS () {
  printf "Checking $1 ... ";
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: ./TRACESPipeDeNovo.sh --install" >&2;
    exit 1;
    else
    echo -e "\e[42mSUCCESS!\e[49m";
    fi
  }
#
CHECK_PROGRAMS () {
  PROGRAM_EXISTS "trimmomatic";
  PROGRAM_EXISTS "metaviralspades.py";
  PROGRAM_EXISTS "FALCON";
  PROGRAM_EXISTS "GeCo3";
  PROGRAM_EXISTS "./AltaiR";
  }
#
SHOW_MENU () {
  echo " -------------------------------------------------------- ";
  echo "                                                          ";
  echo " TRACESPipeDeNovo.sh : TRACESPipe de-novo version v1.0    ";
  echo "                                                          ";
  echo " This is a de-novo version of TRACESPipe. It provide      ";
  echo " automatic reconstruction of viral genomes and performs   ";
  echo " extended analyses such as classification.                ";
  echo "                                                          ";
  echo " Program options ---------------------------------------- ";
  echo "                                                          ";
  echo " -h, --help                     Show this,                ";
  echo " -i, --install                  Installation (w/ conda),  ";
  echo "                                                          ";
  echo " -f, --filter-human             Filter human reads,       ";
  echo "                                                          ";
  echo " -t  <INT>, --threads <INT>     Number of threads,        ";
  echo "                                                          ";
  echo " -mr <INT>, --min-read <INT>    Minimum size of reads,    ";
  echo " -mn <INT>, --min-nc <INT>      Minimum NC to classify,   ";
  echo " -mc <INT>, --min-contig <INT>  Minimum contig length,    ";
  echo "                                                          ";
  echo " -r1 <STR>, --reads1 <STR>      FASTQ reads (forward),    ";
  echo " -r2 <STR>, --reads2 <STR>      FASTQ reads (reverse),    ";
  echo "                                                          ";
  echo " -db <STR>, --database <STR>    FASTA Viral Database,     ";
  echo "                                                          ";
  echo " -o  <STR>, --output <STR>      Output folder name.       ";
  echo "                                                          ";
  echo " Example -----------------------------------------------  ";
  echo "                                                          ";
  echo " bash TRACESPipeDeNovo.sh --reads1 reads_fw.fq.gz \\      ";
  echo "      --reads2 reads_rv.fq.gz --database VDB.mfa \\       ";
  echo "      --output output_analysis --threads 8                ";
  echo "                                                          ";
  echo " -------------------------------------------------------  ";
  }
#
# =============================================================================
# MAIN ------------------------------------------------------------------------
#
if [[ "$#" -lt 1 ]];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help|?)
      HELP=1;
      shift
    ;;
    -c|--install|--compile)
      INSTALL=1;
      shift
    ;;
    -f|--filter-human|--no-human)
      HUMAN_READS=1;
      shift
    ;;
    -t|--threads)
      THREADS="$2";
      shift 2;
    ;;
    -r1|-r|--input1|--reads|--reads1)
      FW_READS="$2";
      RUN=1;
      shift 2;
    ;;
    -r2|--input2|--reads2)
      RV_READS="$2";
      RUN=1;
      shift 2;
    ;;
    -db|--database)
      DATABASE="$2";
      shift 2;
    ;;
    -mn|--min-nc)
      MIN_NC="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -mr|--min-read)
      MIN_READ_LENGTH="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -mc|--min-contig)
      MIN_CONTIG_LENGTH="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -o|--output)
      RESULTS_DIR="$2";
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: TRACESPipeDeNovo.sh -h"
    exit 1;
    ;;
  esac
  done
#
set -- "${POSITIONAL[@]}" # restore positional parameters
#
################################################################################
#
if [[ "$HELP" -eq "1" ]];
  then
  SHOW_MENU;
  exit;
  fi
#
# =============================================================================
# INSTALLATION ----------------------------------------------------------------
#
if [[ "$INSTALL" -eq "1" ]];
  then
  #
  # INSTALL ALTAIR
  git clone https://github.com/cobioders/altair
  cd altair/src/
  cmake . ; make;
  cp AltaiR ../../
  cd ../../
  #
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  #
  conda install -c bioconda trimmomatic -y
  conda install -c bioconda bowtie2 -y
  conda install -c bioconda geco3 -y
  conda install -c bioconda samtools=1.9 -y
  conda install -c bioconda spades=3.5.15 -y
  conda install -c cobilab falcon -y
  #
  wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
  gunzip chm13v2.0.fa.gz
  bowtie2-build chm13v2.0.fa host_DB
  #
  printf ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n>PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n" > adapters.fa
  #
  exit;
  #
  fi
#
# =============================================================================
# RUN -------------------------------------------------------------------------
#
if [[ "$RUN" -eq "1" ]];
  then
  #
  # ===========================================================================
  # CHECK_VERSIONS OF PROGRAMS AND FILES AND SET CONFIGS ----------------------
  #
  echo "Checking data ...";
  CHECK_FILE "$ADAPTERS";
  CHECK_FILE "$FW_READS";
  CHECK_FILE "$RV_READS";
  CHECK_FILE "$DATABASE";
  #
  rm -fr $REPORTS_DIR/
  rm -fr $RESULTS_DIR/
  mkdir -p $REPORTS_DIR/
  mkdir -p $RESULTS_DIR/
  #
  DATE=`date`;
  echo "NEW ANALYSIS: $DATE " > $REPORTS_DIR/report_stderr.txt
  echo "NEW ANALYSIS: $DATE " > $REPORTS_DIR/report_stdout.txt
  #
  CHECK_FILE_GZIPED "$FW_READS" > $REPORTS_DIR/report_stdout.txt
  CHECK_FILE_GZIPED "$RV_READS" > $REPORTS_DIR/report_stdout.txt
  #
  # ===========================================================================
  # REMOVE HUMAN READS --------------------------------------------------------
  #
  if [[ "$HUMAN_READS" -eq "1" ]];
    then
    #
    echo "Removing host data ...";
    CHECK_FILE "$HOST_FASTA";
    #
    # https://www.aespindola.com/post/dehosting/
    # https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
    bowtie2 -p $THREADS -x host_DB -1 $FW_READS -2 $RV_READS --un-conc-gz SAMPLE_host_removed > $STUDY-data.sam 2>> $REPORTS_DIR/report_stderr.txt
    samtools view -bS $STUDY-data.sam > $STUDY-data.bam 2>> $REPORTS_DIR/report_stderr.txt
    samtools view -b -f 12 -F 256 $STUDY-data.bam > $STUDY-data-no-human.bam 2>> $REPORTS_DIR/report_stderr.txt
    samtools sort -n -o $STUDY-sorted-data.bam $STUDY-data-no-human.bam 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
    samtools fastq -1 $STUDY-reads-r1.fq -2 $STUDY-reads-r2.fq $STUDY-sorted-data.bam 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
    #
    else
    #
    echo "Uncompressing reads ...";
    zcat $FW_READS > $STUDY-reads-r1.fq;
    zcat $RV_READS > $STUDY-reads-r2.fq;
    #
  fi
  #
  # ===========================================================================
  # INFORM SIZE OF UNCOMPRESSED READS -----------------------------------------
  #
  FIL_READS_SIZE_1=`ls -lah $STUDY-reads-r1.fq | awk '{ print $5; }'`;
  FIL_READS_SIZE_2=`ls -lah $STUDY-reads-r2.fq | awk '{ print $5; }'`;
  echo "READS1: $FIL_READS_SIZE_1 - READS2: $FIL_READS_SIZE_2";
  #
  # ===========================================================================
  # TRIMMING ------------------------------------------------------------------
  #
  echo "Trimming data ...";
  CHECK_FILE "$ADAPTERS";
  CHECK_FILE "$STUDY-reads-r1.fq";
  CHECK_FILE "$STUDY-reads-r2.fq";
  #
  rm -f $STUDY-fw-pr.fq.gz $STUDY-fw-unpr.fq.gz $STUDY-rv-pr.fq.gz $STUDY-rv-unpr.fq.gz;
  #
  trimmomatic PE -threads $THREADS -phred33 $STUDY-reads-r1.fq $STUDY-reads-r2.fq $STUDY-fw-pr.fq $STUDY-fw-unpr.fq $STUDY-rv-pr.fq $STUDY-rv-unpr.fq ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MIN_READ_LENGTH 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  #
  cat $REPORTS_DIR/report_stderr.txt | grep "Input Read Pairs:";
  cat $STUDY-fw-unpr.fq $STUDY-rv-unpr.fq > $STUDY-unpr.fq
  #
  # ===========================================================================
  # DE-NOVO ASSEMBLY ----------------------------------------------------------
  #
  echo "De-novo assembling data ...";
  CHECK_FILE "$STUDY-fw-pr.fq";
  CHECK_FILE "$STUDY-rv-pr.fq";
  CHECK_FILE "$STUDY-fw-unpr.fq";
  CHECK_FILE "$STUDY-rv-unpr.fq";
  CHECK_FILE "$STUDY-unpr.fq";
  #
  metaspades.py --threads $THREADS -1 $STUDY-fw-pr.fq -2 $STUDY-rv-pr.fq -s $STUDY-unpr.fq -o $RESULTS_DIR/denovo_$STUDY 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  #
  # ===========================================================================
  # ALIGNMENT-FREE REFERENCE-BASED CLASSIFICATION OF CONTIGS ------------------
  #
  echo "Classifying contigs with reference-based approach ...";
  #
  CHECK_FILE "$RESULTS_DIR/denovo_$STUDY/scaffolds.fasta";
  CHECK_FILE "$DATABASE";
  #
  FALCON -v -t 5000 -F $FALCON_PARAM_DB -x $STUDY-metagenomics-all-scaffolds.csv $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta $DATABASE 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  cp $STUDY-metagenomics-all-scaffolds.csv $RESULTS_DIR/
  # 
  # ===========================================================================
  # IDENTIFY UNKNOWN CONTIGS --------------------------------------------------
  #
  echo "Identifying contigs with minimum size of $MIN_CONTIG_LENGTH ...";
  #
  # CALCULATE NC --------------------------------------------------------------
  #
  ./AltaiR nc $NC_PARAM $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta > $STUDY-NC-all.csv;
  cp $STUDY-NC-all.csv $RESULTS_DIR/
  #
  # CALCULATE NRC FOR EACH READ ABOVE THRESHOLD TO THE VIRAL ------------------
  #
  printf "# CONTIG_NAME\tCONTIG_NC\tBEST_REF_LENGTH\tBEST_REF_SIMILARITY\tBEST_REF_NAME\tCONTIG_BPS_USING_BEST_REF\n" > $RESULTS_DIR/FINAL-RESULTS.txt;
  #
  grep ">" $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta | tr -d '>' > $STUDY-HEADERS.txt;
  mapfile -t T_DATA < $STUDY-HEADERS.txt
  #
  for line in "${T_DATA[@]}" #
    do
    #
    PATTERN=`echo "$line" | awk '{ print $1 }'`;
    ./AltaiR filter --pattern "$PATTERN" < $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta > $STUDY-TMP.fa
    #
    NC_VALUE=`grep "$PATTERN" $RESULTS_DIR/$STUDY-NC-all.csv | awk '{ print $2 }'`;
    CONTIG_SIZE=`grep -v ">" $STUDY-TMP.fa | tr -d -c "ACGT" | wc -c`;
    # echo "CONTIG_SIZE=$CONTIG_SIZE, NC_VALUE=$NC_VALUE, MIN_NC=$MIN_NC";
    #
    TOP_VALUE="-";
    if (($(bc <<<"$CONTIG_SIZE < $MIN_CONTIG_LENGTH")));
      then  
      #  
      TOP_VALUE="contig_size_smaller_than_min";
      printf "$PATTERN\t$NC_VALUE\t$TOP_VALUE\n" >> $RESULTS_DIR/FINAL-RESULTS.txt;
      #
    elif (($(bc <<<"$NC_VALUE > $MIN_NC")));
      then
      mkdir -p $RESULTS_DIR/contigs
      #
      FALCON -v -t 1000 -F $FALCON_PARAM_CONTIGS -x $RESULTS_DIR/contigs/$STUDY-$PATTERN.txt $STUDY-TMP.fa $DATABASE 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
      TOP_VALUE=`head -n 1 $RESULTS_DIR/contigs/$STUDY-$PATTERN.txt | awk '{print $2"\t"$3"\t"$4 }'`;
      #
      # PERFORM LOOPBACK NRC:
      BEST_REF_NAME=`head -n 1 $RESULTS_DIR/contigs/$STUDY-$PATTERN.txt  | tr -d ">" | awk '{print $4 }' | tr '_' '\t' | awk '{ print $1; }'`;
      ./AltaiR filter --pattern "$BEST_REF_NAME" < $DATABASE > $STUDY-BEST-REF.fa
      REV_BPS=`GeCo3 $GECO3_PARAM -r $STUDY-BEST-REF.fa $STUDY-TMP.fa |& grep "Total bytes" | awk '{ print $8 }'`; 
      cp $STUDY-BEST-REF.fa xREF.fa
      cp $STUDY-TMP.fa xTAR.fa
      printf "$PATTERN\t$NC_VALUE\t$TOP_VALUE\t$REV_BPS\n" >> $RESULTS_DIR/FINAL-RESULTS.txt;
      #
    else
      TOP_VALUE="contig_sequence_with_nc_smaller_than_min";
      printf "$PATTERN\t$NC_VALUE\t$TOP_VALUE\n" >> $RESULTS_DIR/FINAL-RESULTS.txt;
      fi
    #
    done
  #
  # ===========================================================================
  # REFERENCE-FREE CLASSIFICATION ---------------------------------------------
  #
  echo "Reference-free classification of unknown contigs ...";
  # RUN RFSC ADAPTATION
  #
  # ===========================================================================
  # CLOSE RUN -----------------------------------------------------------------
  #
  fi
#
# Clean temporary data
echo "Cleaning temporary data ...";
rm -f $STUDY-reads-r1.fq $STUDY-reads-r2.fq;
rm -f $STUDY-fw-pr.fq $STUDY-fw-unpr.fq $STUDY-rv-pr.fq $STUDY-rv-unpr.fq;
rm -f $STUDY-metagenomics.csv $STUDY-NC-all.csv $STUDY-HEADERS.txt;
rm -f $STUDY-data.sam $STUDY-data.bam $STUDY-data-no-human.bam $STUDY-sorted-data.bam;
#
###############################################################################
#
