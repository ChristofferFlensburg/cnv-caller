
###############################################################################
#
# The main analysis function, calling this will analyse everything and produce
# plots of a lot of different things.
#
# Fill out the metadata file, and provide the paths to the required data as
# described below.
#
###############################################################################

#This load all the methods. Any issued with uninstalled dependencies will appear here.
source('analyse.R')
loadMethods()

#path to the metadata file: a tab separated file with the named columns as in the example file.
#path to bamfiles are relative from the directory of the metadata file.
metaDataFile = normalizePath('../metaData.txt')

#vector of paths to vcf files. Should include germline and somatic SNVs (and short indels if you want).
#These just point to the locations to be analysed, the sample data isn't used, and additional
#filtering and variant QC is done.
vcfFiles = normalizePath(list.files('../vcf', pattern='*.vcf$', full.names=T))

# The capture regions. A bed file with gene names as fourth column.
# So a tab separated chr start end gene. Use unpadded regions.
# The gene names will appear in plots and output, so you may want these
#   in a human readable format.
captureRegionsFile = normalizePath('../captureRegions/captureRegions.gc.bed')

#db SNP directory. hg19 can be downloaded from https://www.dropbox.com/s/nst57cbvcy9lcgg/hg19_dbSNP.zip?dl=0
#they can also be generated using the script at URL.
dbSNPdirectory = normalizePath('../dbSNP')

#directory with the pool of normals. These normal samples don't have to (but can) be from the same individuals,
#but has to have the same capture regions. The more normals, the better, minimum of 2 to estimate variance.
#These samples should not have a cancer content of more than about 1%.
#The directory has to have a subdirectory names 'bam' with (links to) the bamfiles and index files.
normalDirectory = normalizePath('../normals')

#path to where you want the saved R data to be stored.
#these will not need to be manually accessed.
#The parent directory must exist, but the last level will be created if needed
Rdirectory = '../R'

#path to where you want the plots and other output to go
#this is where you will look at your output
#The parent directory must exist, but the last level directory will be created if needed
plotDirectory = '../plots'

#maximum number of cpus to be used. This can also affect memory usage.
cpus=6

#The base quality phred offset. This can be read from fastqc analysis for example.
#Most technologies have 33 these days, but some have 64.
BQoffset = 33

#the reference genome. Only hg19 is supported this release.
genome = 'hg19'
fastaFile = '../../data/reference/hg19/hg19.fa'

#grouping up settings and metadata
inputFiles =
  list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFiles,
       'normalDirectory'=normalDirectory, 'fastaFile'=fastaFile,
       'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)

#All runtime information is logged to runtimeTracking.log in the Rdirectory.
#Here is an option to also output to terminal as the analysis is going.
runtimeSettings = list('cpus'=cpus, 'outputToTerminalAsWell'=T)
settings = list('genome'=genome, 'BQoffset'=BQoffset)

#If you want to redo some part of the analyses (rather than access the previously saved data
#from earlier runs) you need to change this setting. For example forceRedoEverything().
forceRedo = forceRedoNothing()

#The two parameters of the CNV calling analysis. Essentially setting how much trust to put in
#the read depth and the SNP frequencies.
#SystematicVariance = 0 is the most sensitive to read depth signal, but also more prone to false positives.
#maxCov = Inf is the most sensitive to SNP frequency signals, but also more prone to false positives.
#systematicVariance = 0.03 and maxCov = 150 is default.
parameters = list('systematicVariance'=0.03, 'maxCov'=150)

#Run the analysis.
data = analyse(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings, parameters)

#Run the VEP afterburned to annotate variants.
#Remakes some of the plots and output, adding annotation.
#VEP must be callable from the terminal with the command vep.
postAnalyseVEP(Rdirectory, plotDirectory, parameters)
