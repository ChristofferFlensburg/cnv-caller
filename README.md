cnv-caller
==========

A cnv, snv and clonality caller.

Best way to run it
==================
is to download the self contained example,
make sure that it runs properly, and then run on your own data.

1) download and unzip https://www.dropbox.com/s/egse96gb8n1xq4t/example.zip?dl=0

2) make sure you have all the needed R-packages. (or proceed to 3 and handle the errors)

3) Run the example. (this should confirm that the package works on your computer)

4) Skip to 6), or replace the .R files with the latest version from github if you want.

5) Rerun example, switching on forceredo. (this should confirm that the latest version of the package works on your computer)

6) Set up your data in the same way as in the example.

7) Run on your data by changing the directory paths in exampleMain.R and fill in metadata in metaData.txt.


Second best way to run it is
============================

1) download all the .R files into a directory that you will be running from.

2) Set up your data:

  a) A directory with the bam files of the samples you want to analyse. Include normals.
  
  b) A directory with vcf files of the positions of SNP/SNV (and short indels) that may be interesting. Be very liberal, filtering and flagging is applied.
  
  c) A directory containing subdirectories 'bam' and 'vcf' like in 1) and 2), but with a pool of normal samples run on the same technology, and aligned/variant called with the same pipeline. Can entirely, partially or not at all overlap with the normal samples that are being analysed.
  
  d) directory with dbSNP data. hg19 can be downloaded from https://www.dropbox.com/s/nst57cbvcy9lcgg/hg19_dbSNP.zip?dl=0. Only one supported atm.
  
  e) A bed file with your capture regions. Gene names as you want them to be displayed in the analysis should go in column 4, and gc content of the region goes in column 5. gc content can be added for example with bedtools nuc, but columns may need to be reorganised.
  
3) Fill out a metadata file for the samples according to instructions in main.R.

4) Edit main.R to point at your set up data.

5) open R in the .R directory and source('main.R'). The run will take many hours if you have a big data-set, to run it in screen if possible.
