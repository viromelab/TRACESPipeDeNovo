#!/bin/bash
#
STUDY="datax";
MIN_READ_LENGTH="30";
MIN_SIMILARITY="20";
THREADS="8"; 
FW_READS="";
RV_READS="";
ADAPTERS="adapters.fa";
DATABASE="VDB.mfa";
REPORTS_DIR="reports";
RESULTS_DIR="results";
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
  echo " -t  <INT>, --threads <INT>     Number of threads,        ";
  echo "                                                          ";
  echo " -mr <INT>, --min-read <INT>    Minimum size of reads,    ";
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
    -si|--similarity)
      MIN_SIMILARITY="$2";
      SHOW_HELP=0;
      shift 2;
    ;;
    -mr|--min-read)
      MIN_READ_LENGTH="$2";
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
  conda install -c bioconda spades=3.5.15 -y
  conda install -c cobilab falcon -y
  #
  printf ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n>PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n" > adapters.fa
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
  # ===========================================================================
  # TRIMMING ------------------------------------------------------------------
  #
  echo "Trimming data ...";
  CHECK_FILE "$ADAPTERS";
  CHECK_FILE "$FW_READS";
  CHECK_FILE "$RV_READS";
  rm -f $STUDY-fw-pr.fq.gz $STUDY-fw-unpr.fq.gz $STUDY-rv-pr.fq.gz $STUDY-rv-unpr.fq.gz;
  trimmomatic PE -threads $THREADS -phred33 $FW_READS $RV_READS $STUDY-fw-pr.fq $STUDY-fw-unpr.fq $STUDY-rv-pr.fq $STUDY-rv-unpr.fq ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MIN_READ_LENGTH 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
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
  #metaviralspades.py --threads $THREADS -1 $STUDY-fw-pr.fq -2 $STUDY-rv-pr.fq -s $STUDY-unpr.fq -o $RESULTS_DIR/denovo_$STUDY 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  metaspades.py --threads $THREADS -1 $STUDY-fw-pr.fq -2 $STUDY-rv-pr.fq -s $STUDY-unpr.fq -o $RESULTS_DIR/denovo_$STUDY 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  #
  # ===========================================================================
  # REFERENCE-BASED CLASSIFICATION OF CONTIGS ---------------------------------
  #
  echo "Classifying contigs with reference-based approach ...";
  #CHECK_FILE "$RESULTS_DIR/denovo_$STUDY/before_rr.fasta";
  #FALCON -v -n $THREADS -t 5000 -F -m 6:1:1:0/0 -m 13:50:1:0/0 -m 19:500:1:5/10 -g 0.85 -c 10 -x $STUDY-metagenomics.csv $RESULTS_DIR/denovo_$STUDY/before_rr.fasta $DATABASE 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  CHECK_FILE "$RESULTS_DIR/denovo_$STUDY/scaffolds.fasta";
  CHECK_FILE "$DATABASE";
  FALCON -v -n $THREADS -t 5000 -F -m 6:1:1:0/0 -m 13:50:1:0/0 -m 19:500:1:5/10 -g 0.85 -c 30 -x $STUDY-metagenomics-all-scaffolds.csv $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta $DATABASE 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
  cp $STUDY-metagenomics-all-scaffolds.csv $RESULTS_DIR/
  # 
  # ===========================================================================
  # IDENTIFY UNKNOWN CONTIGS --------------------------------------------------
  #
  echo "Identifying unknown contigs ...";
  #
  # CALCULATE NC --------------------------------------------------------------
  #
  NC_PARAM=" -t $THREADS -m 11:1:1:0:0:0.85/0:0:0 -m 13:50:1:1:0:0.9/0:0:0 ";
  ./AltaiR nc $NC_PARAM $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta > $STUDY-NC-all.csv;
  cp $STUDY-NC-all.csv $RESULTS_DIR/
  #
  # CALCULATE NRC FOR EACH READ ABOVE THRESHOLD TO THE VIRAL ------------------
  #
  rm -f $RESULTS_DIR/FINAL-RESULTS.txt
  grep ">" $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta > $STUDY-HEADERS.txt;
  mapfile -t T_DATA < $1
  #
  for line in "${T_DATA[@]}" #
    do
    #
    ./AltaiR filter --pattern "$line" < $RESULTS_DIR/denovo_$STUDY/scaffolds.fasta > $STUDY-TMP.fa
    #
    NC_VALUE=`grep "$line" $STUDY-NC-all.csv | awk '{ print $2 }'`;
    TOP_VALUE="---";
    #
    if [[ "$MIN_SIMILARITY" > "$NC_VALUE" ]];
      then
      #echo "$NAME $NC ";
      mkdir -p $RESULTS_DIR/contigs
      FALCON -v -n $THREADS -t 5000 -F -m 7:1:1:0/0 -m 13:50:1:0/0 -g 0.85 -c 10 -x $RESULTS_DIR/contigs/$STUDY-$line $STUDY-TMP.fa $DATABASE 1>> $REPORTS_DIR/report_stdout.txt 2>> $REPORTS_DIR/report_stderr.txt
      TOP_VALUE=`head -n 1 $RESULTS_DIR/contigs/$STUDY-$line`;
      #
      fi
    #
    printf "$line\t$NC_VALUE\t$TOP_VALUE\n" >> $RESULTS_DIR/FINAL-RESULTS.txt;
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
rm -f $STUDY-fw-pr.fq $STUDY-fw-unpr.fq $STUDY-rv-pr.fq $STUDY-rv-unpr.fq;
rm -f $STUDY-metagenomics.csv $STUDY-NC-all.csv $STUDY-HEADERS.txt
#
###############################################################################
#
