
###############################################################################
#
# The main analysis function, calling this will analyse everything and produce
# most of the plots and output you need for a paper. The current version will not
# write the paper for you.
#
# Fill out the metadata, and provide the paths to the required data as
# described below.
#
###############################################################################

#path to the metadata file: a tab separated file with the named columns as in the example file.
#path to bamfiles are relative from the directory of the metadata file.
metaDataFile = 'path/to/metaData.txt'

#vector of paths to vcf files. Should include germline and somatic SNVs (and short indels if you want).
#These just point to the locations to be analysed, the sample data isn't used.
vcfFiles = list.files('path/to/vcfDirectory', pattern='*.vcf$', full.names=T)

#the capture regions. These need to be in a specific format, essentially a tab separated .bed file with columns:
#chromosome start end geneName gcContent
#further columns are allowed but will not be used.
#the gene names in this file will appear in the output, so you probably want these to be with symbol names.
#translation can for example be done with biomart.
#gc content can be found through bedtools nuc.
captureRegionsFile = '/path/to/captureRegions.gc.bed'

#db SNP directory. hg19 can be downloaded from https://www.dropbox.com/s/nst57cbvcy9lcgg/hg19_dbSNP.zip?dl=0
#they can also be generated using the script at URL.
dbSNPdirectory = '/path/to/dbSNP'

#directory with the pool of normals. These normal samples don't have to be from the same individuals,
#but should be sequenced on the same technolgy, and aligned and variant called with the same pipeline.
#the more normals, the better.
#These samples are assumed to be diploid and not have any variants of interest.
#The directory has to have a subdirectory names 'bam' with (links to) the bamfiles and index files.
#Also a 'vcf' subdirectory with vcf files of the bams. Again, only the positions in the vcfs are important.
normalDirectory = '/path/to/normalDirectory'

#path to where you want the saved data should be stored
#these will not need to be manually accessed
#The parent directory must exist, but the last level will be created if needed
Rdirectory = '/path/to/R'

#path to where you want the plots and other output to go
#this is where you will look at your output
#The parent directory must exist, but the last level directory will be created if needed
plotDirectory = '/path/to/plots'

inputFiles = list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFile, 'normalDirectory'=normalDirectory,
  'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
runtimeSettings = list('cpus'=cpus)
settings = list('genome'=genome, 'BQoffset'=BQoffset)

analyse(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings)
