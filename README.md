# TRACESPipeDeNovo

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Speed](https://img.shields.io/static/v1.svg?label=Assembly&message=De-Novo&color=green)](#)
[![Release](https://img.shields.io/static/v1.svg?label=Release&message=v2.1.0&color=orange)](#)

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
<work in process...>
```

## 4. Future updates ##

Inclusion of the following viruses:
```
read alignment on reconstructed sequences for coverage plots
```

## 5. Citation ##

Please, cite:
```
Pratas, D., Toppinen, M., Pyöriä, L., Hedman, K., Sajantila, A. and Perdomo, M.F., 2020. 
A hybrid pipeline for reconstruction and analysis of viral genomes at multi-organ level.
GigaScience, 9(8), p.giaa086.
```
[PDF Link](https://doi.org/10.1093/gigascience/giaa086)

## 6. License ##

GPL v3.

For more information see LICENSE file or visit
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

