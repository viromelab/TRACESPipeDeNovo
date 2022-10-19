#!/bin/bash
#
# ==============================================================================
#
lzma -k -f -d VDB.fa.lzma
#
ls *_forward.fq.gz > .r1_names.txt;
ls *_reverse.fq.gz > .r2_names.txt;
#
mapfile -t READS_NAMES_1 < .r1_names.txt;
mapfile -t READS_NAMES_2 < .r2_names.txt;
#
idx=0;
for reads in "${READS_NAMES_1[@]}"
  do
  #
  R1=`echo $reads | awk '{ print $1 }'`;
  R2=`echo ${READS_NAMES_2[$idx]} | awk '{ print $1 }'`;
  echo "Running $R1 and $R2 ...";
  #
  ./TRACESPipeDeNovo.sh --similarity 20 --reads1 $R1 --reads2 $R2 --database VDB.fa --threads 16 --output results_denovo_R-$R1-$R2-ready
  #
  echo "Done!";
  ((++idx));
  done

