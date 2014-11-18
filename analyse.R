
#call dependencies
#source('Rscripts/functions.R')
#source('Rscripts/analysisFunctions.R')

#required packages
require(limma)
require(edgeR)
require(Rsubread)
require(parallel)

analyse = function(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings) {

  if ( !all(c('metaDataFile', 'vcfFiles', 'normalDirectory', 'captureRegionsFile', 'dbSNPdirectory') %in% names(inputFiles)) )
    stop('inputFiles need all entries: metaDataFile, vcfFiles, normalDirectory, captureRegionsFile, dbSNPdirectory.')
  if ( !all(c('Rdirectory', 'plotDirectory') %in% names(outputDirectories)) )
    stop('outputDirectories need all entries: Rdirectory, plotDirectory.')
  if ( !all(c('BQoffset', 'genome') %in% names(settings)) )
    stop('settings need all entries: BQoffset, genome.')
  if ( !all(c('cpus') %in% names(runtimeSettings)) )
    stop('runtimeSettings need all entries: cpus.')
  
  cpus = runtimeSettings$cpus
  genome = settings$genome
  sampleMetaDataFile = inputFiles$metaDataFile
  vcfFiles = inputFiles$vcfFiles
  normalDirectory = inputFiles$normalDirectory
  dbSNPdirectory = inputFiles$dbSNPdirectory
  captureRegionsFile = inputFiles$captureRegionsFile

  Rdirectory = outputDirectories$Rdirectory
  plotDirectory = outputDirectories$plotDirectory
  BQoffset = settings$BQoffset
  
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
  assign('catLog', function(...) cat(..., file=logFile, append=T), envir = .GlobalEnv)
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

  normalRdirectory = paste0(normalDirectory, '/R')
  if ( !file.exists(normalRdirectory) ) dir.create(normalRdirectory)
  externalNormalBams = list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$', full.names=T)
  names(externalNormalBams) = gsub('.bam$', '', list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$'))
  externalNormalVcfs = list.files(path=paste0(normalDirectory, '/vcf'), pattern = '*.vcf$', full.names=T)
  catLog('Normal bamfiles are:\n')
  catLog(externalNormalBams, sep='\n')
  catLog('Normal vcf files are:\n')
  catLog(externalNormalVcfs, sep='\n')

  source('importCaptureRegions.R')
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
  
  #read in metadata
  source('importSampleMetaData.R')
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
  source('runDE.R')
  source('XRank.R')
  fitS = try(runDE(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory,
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
  source('makeFitPlots.R')
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
      source('getVariants.R')
      variants = try(getVariants(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset, dbSNPdirectory,
        Rdirectory, plotDirectory, filterBoring=T, cpus=cpus, v=v, forceRedoSNPs=forceRedoSNPs,
        forceRedoVariants=forceRedoVariants))
      if ( class(variants) == 'try-error' ) {
        catLog('Error in getVariants.\n')
        stop('Error in getVariants.')
      }
      
      #Get variants from the external normals
      normalVariants =
        try(getVariants(externalNormalVcfs, externalNormalBams, names(externalNormalBams), captureRegions,
                        genome, BQoffset, dbSNPdirectory, normalRdirectory, plotDirectory, filterBoring=F, cpus, v,
                        forceRedoSNPs=forceRedoNormalSNPs, forceRedoVariants=forceRedoNormalVariants))
      if ( class(normalVariants) == 'try-error' ) {
        catLog('Error in getVariants for normals.\n')
        stop('Error in getVariants for normals.')
      }
      
      #share variants with normals
      source('matchFlagVariants.R')
      allVariants = try(matchFlagVariants(variants, normalVariants, individuals, normals,
        Rdirectory, v=v, forceRedoMatchFlag=forceRedoMatchFlag))
      if ( class(allVariants) == 'try-error' | !all(c('variants', 'normalVariants') %in% names(allVariants)) ) {
        catLog('Error in matchFlagVariants.\n')
        stop('Error in matchFlagVariants.')
      }
    }
  }
  variants = allVariants$variants
  normalVariants = allVariants$normalVariants

  #Make the scatter plots of the pairs
  source('makeScatterPlots.R')
  scatter = try(makeScatterPlots(variants, samplePairs, timePoints, plotDirectory,
    genome=genome, cpus=cpus, v=v, forceRedo=forceRedoScatters))
  if ( class(scatter) == 'try-error' ) {
    catLog('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.')
  }
  
  #identify and output new variants
  source('outputNewVariants.R')
  newVar = try(outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, v, forceRedo=forceRedoNewVariants))
  if ( class(newVar) == 'try-error' ) {
    catLog('Error in outputNewVariants! Continuing anyway.\n')
    warning('Error in outputNewVariants! Continuing anyway.')
  }
  
  #output somatic variants
  source('outputSomaticVariants.R')
  somatics = try(outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, v, forceRedo=forceRedoOutputSomatic))
  if ( class(somatics) == 'try-error' ) {
    catLog('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.\n')
    warning('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.')
  }
  
  #do multi-sample heatmaps and frequency progression
  source('makeSNPprogressionPlots.R')
  progression = try(makeSNPprogressionPlots(variants, timeSeries, plotDirectory, cpus=cpus, v=v, forceRedo=forceRedoSNPprogression))
  if ( class(progression) == 'try-error' ) {
    catLog('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.')
  }

  #call CNVs compared to the normals.
  source('callCNVs.R')
  cnvs =
    try(callCNVs(variants=variants, normalVariants=normalVariants, fitS=fitS,
                 names=names, individuals=individuals, normals=normals, Rdirectory=Rdirectory,
                 genome=genome, cpus=cpus, v=v, forceRedoCNV=forceRedoCNV))
  if ( class(cnvs) == 'try-error' ) {
    catLog('Error in callCNVs!\n')
    stop('Error in callCNVs.')
  }


  #make CNV plots
  source('makeCNVplots.R')
  cnvplot = try(makeCNVplots(cnvs, plotDirectory=plotDirectory, genome, v=v, forceRedoCNVplots=forceRedoCNVplots))
  if ( class(cnvplot) == 'try-error' ) {
    catLog('Error in makeCNVplots! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeCNVplots! Continuing, but these plots are kindof useful.')
  }

  #make summary plots
  source('makeSummaryPlot.R')
  summary = try(makeSummaryPlot(variants, cnvs, normals, individuals, timePoints, plotDirectory,
    genome, cpus, v, forceRedo=forceRedoSummary))
  if ( class(summary) == 'try-error' ) {
    catLog('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.')
  }

  #combine SNPs and CNVs into stories of subclones.
  source('getStories.R')
  setVariantLoss(normalVariants$variants, v=v)
  stories = try(getStories(variants=variants, cnvs=cnvs, timeSeries=timeSeries, Rdirectory=Rdirectory,
    plotDirectory=plotDirectory, cpus=cpus, v=v, forceRedo=forceRedoStories))
  if ( class(stories) == 'try-error' ) {
    catLog('Error in getStories!\n')
    stop('Error in getStories!')
  }
  
  source('makeRiverPlots.R')
  river = try(makeRiverPlots(stories, variants, plotDirectory, forceRedo=forceRedoRiver))
  if ( class(river) == 'try-error' ) {
    catLog('Error in makeRiverPlots!\n')
    warning('Error in makeRiverPlots!')
  }

  catLog('Run done! :)')
  return()
}
