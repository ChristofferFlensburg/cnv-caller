
#call dependencies
#source('Rscripts/functions.R')
#source('Rscripts/analysisFunctions.R')

#required packages
  require(limma)
  require(edgeR)
  require(Rsubread)
  require(parallel)

analyse = function(inputFiles, outputDirectories, forceRedo, runtimeSettings) {

}

#set default options if not specified
if ( !exists('v') ) v = 'trackProgress'
if ( !exists('cpus') ) cpus=10
if ( !exists('genome') | class(genome) != 'character' ) genome='hg19'

#sanity check input to the R script.
if ( !exists('sampleMetaDataFile') ) stop('Need to set sampleMetaDataFile.')
if ( class(sampleMetaDataFile) != 'character' ) stop('sampleMetaDataFile needs to be of class character.')
if ( !file.exists(sampleMetaDataFile) ) stop('sampleMetaDataFile ', sampleMetaDataFile, ' not found.')

if ( !exists('vcfFiles') ) stop('Need to set vcfFiles!')
if ( class(vcfFiles) != 'character' ) stop('vcfFiles needs to be of class character.')
for ( vcfFile in vcfFiles ) if ( !file.exists(vcfFile) ) stop('vcfFile ', vcfFile, ' not found.')

if ( !exists('normalDirectory') ) stop('Need to set normalDirectory!')
if ( class(normalDirectory) != 'character' ) stop('normalDirectory needs to be of class character.')
if ( !file.exists(normalDirectory) ) stop('normalDirectory ', normalDirectory, ' not found.')

if ( !exists('technology') ) {
  if ( !exists('captureRegionsFile') ) stop('Need to set captureRegionsFile.')
  if ( class(captureRegionsFile) != 'character' ) stop('captureRegionsFile needs to be of class character.')
  if ( !file.exists(captureRegionsFile) ) stop('captureRegionsFile ', captureRegionsFile, ' not found.')
}
  
if ( !exists('Rdirectory') ) stop('Need to set Rdirectory.')
if ( class(Rdirectory) != 'character' ) stop('Rdirectory needs to be of class character.')
if ( !file.exists(Rdirectory) ) {
  dirSuccess = dir.create(Rdirectory)
  if ( !dirSuccess ) stop('Failed to create the Rdirectory ', Rdirectory,
                          '. Note that parent directory must exist.')
}

if ( !exists('plotDirectory') ) stop('Need to set plotDirectory.')
if ( class(plotDirectory) != 'character' ) stop('plotDirectory needs to be of class character.')
if ( !file.exists(plotDirectory) ) {
  dirSuccess = dir.create(plotDirectory)
  if ( !dirSuccess ) stop('Failed to create the plotDirectory ', plotDirectory,
                          '. Note that parent directory must exist.')
}


#set up logfile and log start of run.
logFile = paste0(Rdirectory, '/runtimeTracking.log')
catLog = function(...) cat(..., file=logFile, append=T)
cat('Runtime tracking and QC information printed to ', logFile, '.\n', sep='')
catLog(as.character(Sys.time()),
       '\n\n\n######################################################################\n',
       'Starting run with input files:',
       '\nsampleMetaDataFile:', sampleMetaDataFile,
       '\nvcfFiles:\n')
catLog(vcfFiles, sep='\n')
catLog('Normal directory:', normalDirectory, '\n')
catLog('capture regions:', captureRegionsFile, '\n')
catLog('Plotting to', plotDirectory, '\n')
catLog('Saving R files to', Rdirectory, '\n')
catLog('Genome is', genome, '\n')
catLog('Running on at most', cpus, 'cpus.\n')


if ( !(genome %in% c('hg19', 'mm10')) ) stop('Only genomes that are supported atm are hg19 and mm10, sorry.\nNew genomes can easily be added though, please contact the authors.\n')

#REMOVE FROM PUBLIC RELEASE
if ( exists('technology') ) {
  if ( technology %in% c('Agilent', 'kinomeCapture') ) BQoffset = as.integer(33)
  if ( technology == 'Illumina' ) BQoffset = as.integer(64)
}
if ( !exists('BQoffset') ) {
  BQoffset = 33
  catLog('phred offset BQoffset not specified, so set to default 33. Set this manually if this is not correct.\n')
}


#set forceRedo parameters to false unless already specified
if ( !exists('forceRedoAgilentCaptureNames') ) forceRedoAgilentCaptureNames = F
if ( !exists('forceRedoCount') ) forceRedoCount = F
if ( !exists('forceRedoNormalCount') ) forceRedoNormalCount = F
if ( !exists('forceRedoFit') ) forceRedoFit = F
if ( !exists('forceRedoVolcanoes') ) forceRedoVolcanoes = F
if ( !exists('forceRedoDifferentRegions') ) forceRedoDifferentRegions = F
if ( !exists('forceRedoSNPs') ) forceRedoSNPs = F
if ( !exists('forceRedoVariants') ) forceRedoVariants = F
if ( !exists('forceRedoNormalSNPs') ) forceRedoNormalSNPs = F
if ( !exists('forceRedoNormalVariants') ) forceRedoNormalVariants = F
if ( !exists('forceRedoMatchFlag') ) forceRedoMatchFlag = F
if ( !exists('forceRedoScatters') ) forceRedoScatters = F
if ( !exists('forceRedoOutputSomatic') ) forceRedoOutputSomatic = F
if ( !exists('forceRedoNewVariants') ) forceRedoNewVariants = F
if ( !exists('forceRedoSNPprogression') ) forceRedoSNPprogression = F
if ( !exists('forceRedoCNV') ) forceRedoCNV = F
if ( !exists('forceRedoCNVplots') ) forceRedoCNVplots = F
if ( !exists('forceRedoSummary') ) forceRedoSummary = F
if ( !exists('forceRedoStories') ) forceRedoStories = F
if ( !exists('forceRedoRiver') ) forceRedoRiver = F

#make sure that if something is redone, then all depending steps are redone as well
if ( forceRedoAgilentCaptureNames ) {
  catLog('Redoing capture regions, so need to redo coverage counts and SNPs as well.\n')
  forceRedoCount = T
  forceRedoNormalCount = T
  forceRedoSNPs = T
}
if ( forceRedoCount ) {
  catLog('Redoing coverage counts, so need to redo linear analysis.\n')
  forceRedoFit = T
}
if ( forceRedoNormalCount ) {
  catLog('Redoing normal coverage counts, so need to redo linear analysis of coverage.\n')
  forceRedoFit = T
}
if ( forceRedoFit ) {
  catLog('Redoing linear analysis of coverage, so need to redo CNVs, volcano plots and diffent regions sheet.\n')
  forceRedoCNV = T
  forceRedoVolcanoes = T
  forceRedoDifferentRegions = T
}
if ( forceRedoSNPs ) {
  catLog('Redoing SNPs, so need to redo quality flagging of variants.\n')
  forceRedoVariants = T
}
if ( forceRedoVariants ) {
  catLog('Redoing quality flagging of variants, so need to redo flagmatching with normals.\n')
  forceRedoMatchFlag = T
}
if ( forceRedoNormalSNPs ) {
  catLog('Redoing normal SNPs, so need to redo quality flagging of normal variants.\n')
  forceRedoNormalVariants = T
}
if ( forceRedoNormalVariants ) {
  catLog('Redoing normal quality flagging of variants, so need to redo flag matching with normals.\n')
  forceRedoMatchFlag = T
}
if ( forceRedoMatchFlag ) {
  catLog('Redoing flagmatching of normals, so need to redo frequency scatters, variant sheets, frequency progressions and CNVs.\n')
  forceRedoScatters = T
  forceRedoNewVariants = T
  forceRedoSNPprogression = T
  forceRedoCNV = T
}
if ( forceRedoCNV ) {
  catLog('Redoing CNVs, so need to redo CNV plots, summary plot and clonality stories.\n')
  forceRedoCNVplots = T
  forceRedoSummary = T
  forceRedoStories = T
}
if ( forceRedoStories ) {
  catLog('Redoing stories, so need to redo river plots.\n')
  forceRedoRiver = T
}

#REMOVE FROM PUBLIC RELEASE!
if ( exists('technology') )
{
  if ( genome == 'hg19' & technology == 'Agilent' ) {
    captureRegions = importCaptureRegions(file = '/wehisan/general/academic/grp_leukemia_genomics/AGILENT_EXOME/HUMAN_v5/hg19.gc.bed', genome=genome, gcColumn=6)
    captureRegions = fixAgilentRegionNames(captureRegions, genome, cpus=cpus, forceRedo = forceRedoAgilentCaptureNames)
    cat('Using capture regions at /wehisan/general/academic/grp_leukemia_genomics/AGILENT_EXOME/HUMAN_v5/hg19.gc.bed\n')
  }
  else if ( genome == 'hg19' & technology == 'Illumina') {
    captureRegions = importCaptureRegions(file = '/wehisan/general/academic/grp_leukemia_genomics/F13TSFAPHT0461_HUMwkbX/analysis/captureRegions.gc.bed', genome=genome, gcColumn=8)
    cat('Using capture regions at /wehisan/general/academic/grp_leukemia_genomics/F13TSFAPHT0461_HUMwkbX/analysis/captureRegions.gc.bed\n')
  }
  else if (genome == 'mm10' & technology == 'Agilent') {
    captureRegions = importCaptureRegions(file = '/wehisan/general/academic/grp_leukemia_genomics/AGILENT_EXOME/MOUSE_v1/S0276129/S0276129_Regions.mm10.gc.bed', genome=genome, gcColumn=6)
    captureRegions = fixAgilentRegionNames(captureRegions, genome, cpus=cpus, forceRedo = forceRedoAgilentCaptureNames)
    cat('Using capture regions at /wehisan/general/academic/grp_leukemia_genomics/AGILENT_EXOME/MOUSE_v1/S0276129/S0276129_Regions.mm10.gc.bed\n')
  }
  else if (genome == 'hg19' & technology == 'kinomeCapture') {
    captureRegions = importCaptureRegions(file = '/wehisan/general/academic/grp_leukemia_genomics/DLBCL\ PROJECT/analysis/kinomeCaptureRegions.gc.bed', genome=genome, gcColumn=6)
    captureRegions = fixKinomeRegionNames(captureRegions, genome, cpus=cpus, forceRedo = forceRedoAgilentCaptureNames)
    cat('Using capture regions at /wehisan/general/academic/grp_leukemia_genomics/DLBCL\ PROJECT/analysis/kinomeCaptureRegions.gc.bed\n')
  }
  else {
    cat('Dont know about capture regions for genome', genome, 'with technology', technology, '.\n')
  }
  if ( technology == 'Agilent' & genome == 'hg19' )
    normalDirectory = '/wehisan/general/academic/grp_leukemia_genomics/normals/Agilent'
  else if ( technology == 'Illumina' & genome == 'hg19' )
    normalDirectory = '/wehisan/general/academic/grp_leukemia_genomics/normals/Illumina'
  else if ( technology == 'Agilent' & genome == 'mm10' )
    normalDirectory = '/wehisan/general/academic/grp_leukemia_genomics/normals/Agilent/mouse'
  else if ( technology == 'kinomeCapture' & genome == 'hg19' )
    normalDirectory = '/wehisan/general/academic/grp_leukemia_genomics/normals/kinome'
}

normalRdirectory = paste0(normalDirectory, '/R')
if ( !file.exists(normalRdirectory) ) dir.create(normalRdirectory)
externalNormalBams = list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$', full.names=T)
names(externalNormalBams) = gsub('.bam$', '', list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$'))
externalNormalVcfs = list.files(path=paste0(normalDirectory, '/vcf'), pattern = '*.vcf$', full.names=T)
catLog('Normal bamfiles are:\n')
catLog(externalNormalBams, sep='\n')
catLog('Normal vcf files are:\n')
catLog(externalNormalVcfs, sep='\n')

#MAKE THIS DEFAULT ON RELASE! remove wrapping if.
if ( !exists('captureRegions') ) {
  captureRegions = try(importCaptureRegions(captureRegionsFile, gcColumn=5))
  if ( class(captureRegions) != 'GRanges' ) {
    catLog('Failed to import capture regions, aborting.\n')
    stop('Failed to import capture regions, aborting.\n')
  }
  if ( length(captureRegions) == 0 ) {
    catLog('Empty capture regions, aborting.\n')
    stop('Empty capture regions, aborting.\n')
  }
  if ( !('region' %in% colnames(mcols(captureRegions)) & 'gc' %in% colnames(mcols(captureRegions))) ) {
    catLog('Need both region name and gc content of capture regions, aborting.\n')
    stop('Need both region name and gc content of capture regions, aborting.\n')
  }
  if ( any(is.na(captureRegions$gc)) ) {
    catLog('NA gc content value in capture regions, aborting.\n')
    stop('NA gc content value in capture regions, aborting.\n')
  }
  catLog('Imported capture regions with', length(captureRegions), 'regions and',
         length(unique(captureRegions$region)), 'unique gene names.\n')
  catLog('Mean GC content is ', round(mean(captureRegions$gc), 3), '.\n', sep='')
}

#read in metadata
sampleMetaData = try(importSampleMetaData(sampleMetaDataFile))
if ( class(sampleMetaData) != 'data.frame' ) {
  catLog('Failed to import meta data, aborting.\n')
  stop('Failed to import meta data, aborting.\n')
}
if ( !all(c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL') %in% colnames(sampleMetaData)) ) {
  missing = c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL')[!(c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL') %in% colnames(sampleMetaData))]
  catLog('Missing columns in meta data:' , missing, ', aborting.\n')
  stop('Missing columns in meta data:' , missing, ', aborting.\n')
}
if ( !all(captureRegions$NORMAL %in% c('YES', 'NO')) ) {
  catLog('Want only YES or NO in normal column, aborting.\n')
  stop('Want only YES or NO in normal column, aborting.\n')
}
bamFiles = paste0(dirname(sampleMetaDataFile), '/', sampleMetaData$BAM)
names = make.names(sampleMetaData$NAME, unique=T)
individuals = sampleMetaData$INDIVIDUAL
timePoints = sampleMetaData$TIMEPOINT
names(timePoints) = names(individuals) = names
normals = as.logical(gsub('YES', 'T', gsub('NO', 'F', sampleMetaData$NORMAL)))
samplePairs = metaToSamplePairs(names, individuals, normals, v)
timeSeries = metaToTimeSeries(names, individuals, normals, v)

catLog('#############################################################\n\n',
       as.character(Sys.time()),'\n',
       'Imported and sanity checked meta data. Looking good so far!\n',
       'metadata:\n')
catLog(colnames(sampleMetaData), '\n', sep='   ')
for ( row in 1:nrow(sampleMetaData) )
  catLog(as.matrix(sampleMetaData[row,]), '\n', sep='   ')
catLog('\ntimeSeries:\n')
for ( ts in timeSeries )
  catLog(ts, '\n', sep='   ')
catLog('\n#############################################################\n\n')

#compare coverage of samples to the pool of normals, using limma-voom.
fitS = try(runDEforSamples(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory,
  normalRdirectory, v=v, cpus=cpus, forceRedoFit=forceRedoFit, forceRedoCount=forceRedoCount,
  forceRedoNormalCount=forceRedoNormalCount))
if ( class('fitS') == 'try-error' ) {
  catLog('Error in runDEforSamples! Input was:')
  catLog('bamFiles:', bamFiles, '\nnames:', names, '\nexternalNormalBams:', externalNormalBams,
         '\nstart(captureRegions)[1:4]:', start(captureRegions)[1:4], '\nRdirectory:', Rdirectory,
         '\nplotDirectory:', plotDirectory, '\nnormalRdirectory:', normalRdirectory, '\ncpus:', cpus,
         '\nforceRedoFit:', forceRedoFit, '\nforceRedoCount', forceRedoCount,
         '\nforceRedoNormalCount', forceRedoNormalCount, '\n')
  stop('Error in runDEforSamples.')
}


#Plot volcanoes and output an excel file with top DE regions.
ret = try(makeFitPlots(fitS, plotDirectory, v=v,
  forceRedoVolcanoes=forceRedoVolcanoes, forceRedoDifferentRegions=forceRedoDifferentRegions))
if ( class(ret) == 'try-error' ) {
  catLog('Error in makeFitPlots, will continue analysis anyway.')
  warning('Error in makeFitPlots.')
}


saveFile = paste0(Rdirectory, '/allVariants.Rdata')
{
  if ( file.exists(saveFile) & !forceRedoMatchFlag & !forceRedoVariants & !forceRedoNormalVariants ) {
    catLog('Loading final version of combined variants.\n')
    load(file=saveFile)
  }
  else {
    #import, filter and QC the variants. Save to file.
    #The information about normals is used for QC, as there will be only true frequencies of 0, 0.5 and 1 in those samples.
    variants = try(getVariants(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset,
      Rdirectory, filterBoring=T, cpus=cpus, v=v, forceRedoSNPs=forceRedoSNPs,
      forceRedoVariants=forceRedoVariants))
    if ( class(variants) == 'try-error' ) {
      catLog('Error in getVariants.\n')
      stop('Error in getVariants.')
    }
    #TODO: HOW TO HANDLE DB-SNPS???????????? REQUEST THE INFORMATION IN THE VCF??
    #TODO: chromosome names?? with or without 'chr'??

    #Get variants from the external normals
    normalVariants =
      try(getVariants(externalNormalVcfs, externalNormalBams, names(externalNormalBams),
                      captureRegions, genome, BQoffset, normalRdirectory, filterBoring=F, cpus, v,
                      forceRedoSNPs=forceRedoNormalSNPs, forceRedoVariants=forceRedoNormalVariants))
    if ( class(normalVariants) == 'try-error' ) {
      catLog('Error in getVariants for normals.\n')
      stop('Error in getVariants for normals.')
    }
    
    #share variants with normals
    allVariants = try(matchFlagVariants(variants, normalVariants, individuals, normals, Rdirectory, v=v, forceRedoMatchFlag=forceRedoMatchFlag))
    if ( class(allVariants) == 'try-error' | !all(c('variants', 'normalVariants') %in% names(allVariants)) ) {
      catLog('Error in matchFlagVariants.\n')
      stop('Error in matchFlagVariants.')
    }
  }
}

variants = allVariants$variants
normalVariants = allVariants$normalVariants

#Make the scatter plots of the pairs
scatter = try(makeScatterPlots(variants, samplePairs, timePoints, plotDirectory, genome=genome, cpus=cpus, v=v, forceRedo=forceRedoScatters))
if ( class(scatter) == 'try-error' ) {
  catLog('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.\n')
  warning('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.')
}

#identify and output new variants
newVar = try(outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, v, forceRedo=forceRedoNewVariants))
if ( class(newVar) == 'try-error' ) {
  catLog('Error in outputNewVariants! Continuing anyway.\n')
  warning('Error in outputNewVariants! Continuing anyway.')
}

#output somatic variants
somatics = try(outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, v, forceRedo=forceRedoOutputSomatic))
if ( class(somatics) == 'try-error' ) {
  catLog('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.\n')
  warning('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.')
}

#do multi-sample heatmaps and frequency progression
progression = try(makeSNPprogressionPlots(variants, timeSeries, plotDirectory, cpus=cpus, v=v, forceRedo=forceRedoSNPprogression))
if ( class(progression) == 'try-error' ) {
  catLog('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.\n')
  warning('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.')
}




#call CNVs compared to the normals.
#If the same individual has a normal sample, use that sample to mark germline SNPs.
cnvs =
  try(callCNVs(variants=variants, normalVariants=normalVariants, fitS=fitS,
           names=names, individuals=individuals, normals=normals, Rdirectory=Rdirectory,
               genome=genome, cpus=cpus, v=v, forceRedoCNV=forceRedoCNV))
if ( class(cnvs) == 'try-error' ) {
  catLog('Error in callCNVs!\n')
  stop('Error in callCNVs.')
}


#make CNV plots
cnvplot = try(makeCNVplots(cnvs, plotDirectory=plotDirectory, genome, v=v, forceRedoCNVplots=forceRedoCNVplots))
if ( class(cnvplot) == 'try-error' ) {
  catLog('Error in makeCNVplots! Continuing, but these plots are kindof useful.\n')
  warning('Error in makeCNVplots! Continuing, but these plots are kindof useful.')
}

#make summary plots
summary = try(makeSummaryPlot(variants, cnvs, normals, individuals, timePoints, plotDirectory,
  genome, cpus, v, forceRedo=forceRedoSummary))
if ( class(summary) == 'try-error' ) {
  catLog('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.\n')
  warning('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.')
}


#combine SNPs and CNVs into stories of subclones.
setVariantLoss(normalVariants$variants, v=v)
stories = try(getStories(variants=variants, cnvs=cnvs, timeSeries=timeSeries, Rdirectory=Rdirectory,
  plotDirectory=plotDirectory, cpus=cpus, v=v, forceRedo=forceRedoStories))
if ( class(stories) == 'try-error' ) {
  catLog('Error in getStories!\n')
  stop('Error in getStories!')
}

river = try(makeRiverPlots(stories, variants, plotDirectory, forceRedo=forceRedoRiver))
if ( class(river) == 'try-error' ) {
  catLog('Error in makeRiverPlots!\n')
  warning('Error in makeRiverPlots!')
}
