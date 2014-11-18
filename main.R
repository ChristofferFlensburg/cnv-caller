
###############################################################################
#
# The main analysis function, calling this will analyse everything and produce
# most of the plots and output you need for a paper. The current version will not
# write the paper for you.
#
# Input:
#   inputFiles:
#     A named list with entries
#     'metaData' absolute path to the tab separated meta data file.
#     'vcfFiles' vector of absolute paths to vcf files.
#
#  outputDirectories:
#    A named list with entries
#    'R' absolute path to the directory where data from the analysis will be stored.
#    'plots' absolute path to the directory where plots will be placed.
#
#  forceRedo:
#    A named logical vector indicating what steps should be forced to be redone,
#    even if saved results are already in place. forceRedo() gives the default
#    setting with all the possible steps set to FALSE.
#
###############################################################################

#path to the metadata file
metaDataFile = 'path/to/metaData.txt'

#vector of paths to vcf files. Should include germline and somatic SNVs (and short indels if you want).
#These just point to the locations to be analysed, the sample data isn't used.
vcfFiles = list.files('path/to/vcfDirectory', pattern='*.vcf$', full.names=T)

#the capture regions. These need to be in a specific format, essentially a tab separated .bed file with columns:
#chromosome start end geneName gcContent
captureRegionsFile = '/path/to/captureRegions.gc.bed'

#db SNP directory. these can be downloaded from PATH, or created from .flat files using the script SCRIPT.
dbSNPdirectory = '/path/to/dbSNP'

#directory with the pool of normals. These normal samples don't have to be from the same individuals,
#but should be sequenced on the same technolgy, and aligned and variant called with the same pipeline.
#These samples are assumed to be diploid and not have any variants of interest.

#path to where you want the saved data should be stored
#these will not need to be manually accessed
#The parent directory must exist, but the last level will be created if needed
Rdirectory = '/path/to/R'

#path to where you want the plots and other output to go
#this is where you will look at your output
#The parent directory must exist, but the last level directory will be created if needed
plotDirectory = '/path/to/plots'

inputFiles = list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFile,
  'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)

analyse(inputFiles, outputDirectories, forceRedo, runtimeSettings)


forceRedo = function() {

}
