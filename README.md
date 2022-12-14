# TRACESPipeDeNovo

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Assembly&message=De-Novo&color=green)](#)
[![Release](https://img.shields.io/static/v1.svg?label=Release&message=v1.1.0&color=orange)](#)

TRACESPipeDeNovo is the reference-free version of TRACESPipe. TRACESPipeDeNovo is oriented for unknown viral genome classification, specifically to apply directly in environment samples. It provides the de-novo reconstruction of viral genomes. The classification is reference and feature-based. The viral database is freely provided exclusively for the classification phase. Assembly is exclusively reference-free.

## 1. Installation ##

To install TRACESPipeDeNovo, please run
```
git clone https://github.com/viromelab/TRACESPipeDeNovo.git
cd TRACESPipeDeNovo/src/
lzma -d VDB.fa.lzma
chmod +x *.sh
./TRACESPipeDeNovo.sh --install
```

## 2. Running example ##

The TRACESPipeDeNovo package includes viral FASTQ reads for a quick demonstration.
To run the example, please, first, install the tool. Then, run the following command
```
./TRACESPipeDeNovo.sh --threads 8 --reads1 reads_forward.fq.gz \
--reads2 reads_forward.fq.gz --database VDB.fa --output test_viral_analysis
```

## 3. Usage ##

To see the option of TRACESPipeDeNovo, please run the following command
```
./TRACESPipeDeNovo.sh -h
```
This command will output the following content
```
 --------------------------------------------------------
                                                         
 TRACESPipeDeNovo.sh : TRACESPipe de-novo version v1.0   
                                                         
 This is a de-novo version of TRACESPipe. It provide     
 automatic reconstruction of viral genomes and performs  
 extended analyses such as classification.               
                                                         
 Program options ----------------------------------------
                                                         
 -h, --help                     Show this,               
 -i, --install                  Installation (w/ conda), 
                                                         
 -t  <INT>, --threads <INT>     Number of threads,       
                                                         
 -mr <INT>, --min-read <INT>    Minimum size of reads,   
                                                         
 -r1 <STR>, --reads1 <STR>      FASTQ reads (forward),   
 -r2 <STR>, --reads2 <STR>      FASTQ reads (reverse),   
                                                         
 -db <STR>, --database <STR>    FASTA Viral Database,    
                                                         
 -o  <STR>, --output <STR>      Output folder name.      
                                                         
 Example ----------------------------------------------- 
                                                         
 bash TRACESPipeDeNovo.sh --reads1 reads_fw.fq.gz \     
      --reads2 reads_rv.fq.gz --database VDB.mfa \      
      --output output_analysis --threads 8               
                                                         
 -------------------------------------------------------
```

## 4. Future updates ##

Inclusion of the following viruses:
```
read alignment on reconstructed sequences for coverage plots
```

## 5. Citation ##

Please, cite:
```
Pratas, D., Toppinen, M., Py??ri??, L., Hedman, K., Sajantila, A. and Perdomo, M.F., 2020. 
A hybrid pipeline for reconstruction and analysis of viral genomes at multi-organ level.
GigaScience, 9(8), p.giaa086.
```
[PDF Link](https://doi.org/10.1093/gigascience/giaa086)

## 6. License ##

GPL v3.

For more information see LICENSE file or visit
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

