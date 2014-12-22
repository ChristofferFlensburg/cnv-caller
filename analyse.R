
#' Analyse exomes
#'
#' This function runs a full SNV, SNP, CNV clonality analysis in the input exome data.
#' @param inputFiles A named list of input files, containing the entries metaDataFile, vcfFiles, normalDirectory, captureRegionsFile and dbSNPdirectory
#' @param outputDirectories A named list of output directories, containing the entries Rdirectory and plotDirectory where the saved data and plots will be stored respectively.
#' @param settings A named list containing the entries genome and BQoffset. The only genome supporter atm is 'hg19', and the BQ offset is 33 for most exomes, altough some have 64. Check your fastqc files if you are not sure.
#' @param forceRedo A named list of logicals controling if existing saved data should be loaded or regenerated (overwriting the previous saved data). shortcuts to create these lists are forceRedoNothing() and forceRedoEverything().
#' @param runtimeSettings A named list containing the entries cpus and outputToTerminalAsWell. cpus is an integer controling the maximum number of cpus used in parallel, and outputToTerminalAsWell prints log data to the R session as well as to the log file.
#' @keywords analyse exomes CNV clonality
#' @export
#' @examples
#' analyse()

analyse = function(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings) {
  source('debug.R')

  if ( !all(c('Rdirectory', 'plotDirectory') %in% names(outputDirectories)) )
    stop('outputDirectories need all entries: Rdirectory, plotDirectory.')
  Rdirectory = outputDirectories$Rdirectory
  if ( !exists('Rdirectory') ) stop('Need to set Rdirectory.')
  if ( class(Rdirectory) != 'character' ) stop('Rdirectory needs to be of class character.')
  
  if ( !file.exists(Rdirectory) ) {
    dirSuccess = dir.create(Rdirectory)
    if ( !dirSuccess ) stop('Failed to create the Rdirectory ', Rdirectory,
                            '. Note that parent directory must exist.')
  }

  #set up logfile and log start of run.
  outputToTerminalAsWell = runtimeSettings$outputToTerminalAsWell
  if ( is.null(outputToTerminalAsWell) | class(outputToTerminalAsWell) != 'logical' ) outputToTerminalAsWell=F
  logFile = paste0(Rdirectory, '/runtimeTracking.log')
  assign('catLog', function(...) cat(..., file=logFile, append=T), envir = .GlobalEnv)
  if ( outputToTerminalAsWell )
    assign('catLog', function(...) {cat(..., file=logFile, append=T); cat(...)}, envir = .GlobalEnv)

  catLog('\n\n\n', as.character(Sys.time()),
         '\n######################################################################\n')


  libraryLoaded = library(limma, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load limma.\n')
    stop('Failed to load limma.')
  }
  libraryLoaded = library(edgeR, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load edgeR.\n')
    stop('Failed to load edgeR.')
  }
  libraryLoaded = library(Rsubread, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load Rsubread.\n')
    stop('Failed to load Rsubread.')
  }
  libraryLoaded = library(parallel, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load parallel.\n')
    stop('Failed to load parallel.')
  }
  libraryLoaded = library(Rsamtools, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load Rsamtools.\n')
    stop('Failed to load Rsamtools.')
  }
  libraryLoaded = library(GenomicRanges, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load GenomicRanges.\n')
    stop('Failed to load GenomicRanges.')
  }
  libraryLoaded = library(R.oo, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load R.oo.\n')
    stop('Failed to load R.oo.')
  }
  libraryLoaded = library(rtracklayer, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load rtracklayer.\n')
    stop('Failed to load rtracklayer.')
  }
  libraryLoaded = library(WriteXLS, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load WriteXLS.\n')
    stop('Failed to load WriteXLS.')
  }
  libraryLoaded = library(fields, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load fields.\n')
    stop('Failed to load fields.')
  }
  catLog('All libraries loaded successfully.\n')

  if ( !all(c('metaDataFile', 'vcfFiles', 'normalDirectory', 'captureRegionsFile', 'dbSNPdirectory') %in% names(inputFiles)) )
    stop('inputFiles need all entries: metaDataFile, vcfFiles, normalDirectory, captureRegionsFile, dbSNPdirectory.')
  if ( !all(c('BQoffset', 'genome') %in% names(settings)) )
    stop('settings need all entries: BQoffset, genome.')
  if ( !all(c('cpus') %in% names(runtimeSettings)) )
    stop('runtimeSettings need all entries: cpus.')
  
  cpus = runtimeSettings$cpus
  
  genome = settings$genome
  sampleMetaDataFile = inputFiles$metaDataFile
  vcfFiles = inputFiles$vcfFiles
  normalDirectory = inputFiles$normalDirectory
  normalRdirectory = paste0(normalDirectory, '/R')
  if ( !file.exists(normalRdirectory) ) {
    dirSuccess = dir.create(normalRdirectory)
    if ( !dirSuccess ) stop('Failed to create the normal R directory ', normalRdirectory,
                            '. Note that parent directory must exist.')
  }
  dbSNPdirectory = inputFiles$dbSNPdirectory
  captureRegionsFile = inputFiles$captureRegionsFile

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
  
    
  if ( !exists('plotDirectory') ) stop('Need to set plotDirectory.')
  if ( class(plotDirectory) != 'character' ) stop('plotDirectory needs to be of class character.')
  if ( !file.exists(plotDirectory) ) {
    dirSuccess = dir.create(plotDirectory)
    if ( !dirSuccess ) stop('Failed to create the plotDirectory ', plotDirectory,
                            '. Note that parent directory must exist.')
  }



  cat('Runtime tracking and QC information printed to ', logFile, '.\n', sep='')
  catLog('Starting run with input files:',
         '\nsampleMetaDataFile:', sampleMetaDataFile,
         '\nvcfFiles:\n')
  catLog(vcfFiles, sep='\n')
  catLog('Normal directory:', normalDirectory, '\n')
  catLog('dbSNP directory:', dbSNPdirectory, '\n')
  catLog('capture regions:', captureRegionsFile, '\n')
  catLog('Plotting to', plotDirectory, '\n')
  catLog('Saving R files to', Rdirectory, '\n')
  catLog('Genome is', genome, '\n')
  catLog('Running on at most', cpus, 'cpus.\n')


  if ( !(genome %in% c('hg19', 'mm10')) ) stop('Only genomes that are supported atm are hg19 and mm10, sorry.\nNew genomes can easily be added though, please contact the authors.\n')



  #set forceRedo parameters to false unless already specified
  forceRedoCount = forceRedo$forceRedoCount
  forceRedoNormalCount = forceRedo$forceRedoNormalCount
  forceRedoFit = forceRedo$forceRedoFit
  forceRedoVolcanoes = forceRedo$forceRedoVolcanoes
  forceRedoDifferentRegions = forceRedo$forceRedoDifferentRegions
  forceRedoSNPs = forceRedo$forceRedoSNPs
  forceRedoVariants = forceRedo$forceRedoVariants
  forceRedoNormalSNPs = forceRedo$forceRedoNormalSNPs
  forceRedoNormalVariants = forceRedo$forceRedoNormalVariants
  forceRedoMatchFlag = forceRedo$forceRedoMatchFlag
  forceRedoScatters = forceRedo$forceRedoScatters
  forceRedoOutputSomatic = forceRedo$forceRedoOutputSomatic
  forceRedoNewVariants = forceRedo$forceRedoNewVariants
  forceRedoSNPprogression = forceRedo$forceRedoSNPprogression
  forceRedoCNV = forceRedo$forceRedoCNV
  forceRedoCNVplots = forceRedo$forceRedoCNVplots
  forceRedoSummary = forceRedo$forceRedoSummary
  forceRedoStories = forceRedo$forceRedoStories
  forceRedoRiver = forceRedo$forceRedoRiver

#make sure that if something is redone, then all depending steps are redone as well
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

  externalNormalBams = list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$', full.names=T)
  names(externalNormalBams) = gsub('.bam$', '', list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$'))
  externalNormalVcfs = list.files(path=paste0(normalDirectory, '/vcf'), pattern = '*.vcf$', full.names=T)
  catLog('Normal bamfiles are:\n')
  catLog(externalNormalBams, sep='\n')
  catLog('Normal vcf files are:\n')
  catLog(externalNormalVcfs, sep='\n')

  source('importCaptureRegions.R')
  captureRegions = try(importCaptureRegions(captureRegionsFile, gcColumn=5, genome=genome))
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
  names(normals) = names
  samplePairs = metaToSamplePairs(names, individuals, normals)
  timeSeries = metaToTimeSeries(names, individuals, normals)

  if ( any(!file.exists(bamFiles)) ) {
    catLog('Missing (or misnamed) bam files:' , bamFiles[!file.exists(bamFiles)], ', aborting.\n')
    stop('Missing (or misnamed) bam files:' , bamFiles[!file.exists(bamFiles)], ', aborting.\n')
  }
  bamIndexFiles = paste0(bamFiles, '.bai')
  if ( any(!file.exists(bamIndexFiles)) ) {
    catLog('Missing (or misnamed) bam index files:' , bamIndexFiles[!file.exists(bamIndexFiles)], ', aborting.\n')
    stop('Missing (or misnamed) bam index files:' , bamIndexFiles[!file.exists(bamIndexFiles)], ', aborting.\n')
  }
  
  catLog('##################################################################################################\n\n',
         as.character(Sys.time()),'\n',
         'Imported and sanity checked meta data. Looking good so far!\n',
         'metadata:\n')
  catLog('',colnames(sampleMetaData), '\n', sep='   ')
  for ( row in 1:nrow(sampleMetaData) )
    catLog('',as.matrix(sampleMetaData[row,]), '\n', sep='   ')
  catLog('\n timeSeries:\n')
  for ( ts in timeSeries )
    catLog('',ts, '\n', sep='   ')
  catLog('\n##################################################################################################\n\n')

  #compare coverage of samples to the pool of normals, using limma-voom.
  source('runDE.R')
  source('XRank.R')
  fitS = try(runDE(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory,
    normalRdirectory, cpus=cpus, forceRedoFit=forceRedoFit, forceRedoCount=forceRedoCount,
    forceRedoNormalCount=forceRedoNormalCount))
  if ( class(fitS) == 'try-error' ) {
    catLog('Error in runDEforSamples! Input was:')
    catLog('bamFiles:', bamFiles, '\nnames:', names, '\nexternalNormalBams:', externalNormalBams,
           '\nstart(captureRegions)[1:4]:', start(captureRegions)[1:4], '\nRdirectory:', Rdirectory,
           '\nplotDirectory:', plotDirectory, '\nnormalRdirectory:', normalRdirectory, '\ncpus:', cpus,
           '\nforceRedoFit:', forceRedoFit, '\nforceRedoCount', forceRedoCount,
           '\nforceRedoNormalCount', forceRedoNormalCount, '\n')
    dumpInput(Rdirectory, list('bamFiles'=bamFiles, 'names'=names, 'externalNormalBams'=externalNormalBams,
                               'captureRegions'=captureRegions, 'Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory,
                               'normalRdirectory'=normalRdirectory, 'cpus'=cpus, 'forceRedoFit'=forceRedoFit,
                               'forceRedoCount'=forceRedoCount, 'forceRedoNormalCount'=forceRedoNormalCount))
    stop('Error in runDEforSamples.')
  }

  #Plot volcanoes and output an excel file with top DE regions.
  source('makeFitPlots.R')
  ret = try(makeFitPlots(fitS, plotDirectory, genome,
    forceRedoVolcanoes=forceRedoVolcanoes, forceRedoDifferentRegions=forceRedoDifferentRegions))
  if ( class(ret) == 'try-error' ) {
    catLog('Error in makeFitPlots, will continue analysis anyway.')
    warning('Error in makeFitPlots.')
  }


  source('getVariants.R')
  source('matchFlagVariants.R')
  saveFile = paste0(Rdirectory, '/allVariants.Rdata')
  {
    if ( file.exists(saveFile) & !forceRedoMatchFlag & !forceRedoVariants & !forceRedoNormalVariants ) {
      catLog('Loading final version of combined variants.\n')
      load(file=saveFile)
    }
    else {
      #import, filter and QC the variants. Save to file.
      #The information about normals is used for QC, as there will be only true frequencies of 0, 0.5 and 1 in those samples.
      variants = try(getVariants(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset, dbSNPdirectory,
        Rdirectory, plotDirectory, cpus=cpus, forceRedoSNPs=forceRedoSNPs,
        forceRedoVariants=forceRedoVariants))
      if ( class(variants) == 'try-error' ) {
        catLog('Error in getVariants.\n')
        stop('Error in getVariants.')
      }
      
      #Get variants from the external normals
      normalVariants =
        try(getNormalVariants(variants, externalNormalBams, names(externalNormalBams), captureRegions,
                        genome, BQoffset, Rdirectory, Rdirectory, plotDirectory, cpus=cpus,
                        forceRedoSNPs=forceRedoNormalSNPs, forceRedoVariants=forceRedoNormalVariants))
      if ( class(normalVariants) == 'try-error' ) {
        catLog('Error in getVariants for normals.\n')
        dumpInput(Rdirectory, list('variants'=variants, 'externalNormalBams'=externalNormalBams,
                                   'names'=names(externalNormalBams), 'captureRegions'=captureRegions,
                                   'genome'=genome, 'BQoffset'=BQoffset, 'dbSNPdirectory'=dbSNPdirectory,
                                   'Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory, 'cpus'=cpus,
                                   'forceRedoSNPs'=forceRedoNormalSNPs, 'forceRedoVariants'=forceRedoNormalVariants))
        stop('Error in getVariants for normals.')
      }
      
      #share variants with normals
      allVariants = try(matchFlagVariants(variants, normalVariants, individuals, normals,
        Rdirectory, forceRedoMatchFlag=forceRedoMatchFlag))
      if ( class(allVariants) == 'try-error' | !all(c('variants', 'normalVariants') %in% names(allVariants)) ) {
        catLog('Error in matchFlagVariants.\n')
        dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'individuals'=individuals, 'normals'=normals,
                                   'Rdirectory'=Rdirectory, 'forceRedoMatchFlag'=forceRedoMatchFlag))
        stop('Error in matchFlagVariants.')
      }
    }
  }
  variants = allVariants$variants
  normalVariants = allVariants$normalVariants
  a=try(setVariantLoss(normalVariants$variants))
  if ( class(a) == 'try-error' ) {
    catLog('Error in setVariantLoss(normalVariants)\n')
    stop('Error in setVariantLoss(normalVariants)!')
  }
  

  #Make the scatter plots of the pairs
  source('makeScatterPlots.R')
  scatter = try(makeScatterPlots(variants, samplePairs, timePoints, plotDirectory,
    genome=genome, cpus=cpus, forceRedo=forceRedoScatters))
  if ( class(scatter) == 'try-error' ) {
    catLog('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.')
  }
  
  #identify and output new variants
  source('outputNewVariants.R')
  newVar = try(outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, forceRedo=forceRedoNewVariants))
  if ( class(newVar) == 'try-error' ) {
    catLog('Error in outputNewVariants! Continuing anyway.\n')
    warning('Error in outputNewVariants! Continuing anyway.')
  }
  
  #output somatic variants
  source('outputSomaticVariants.R')
  somatics = try(outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, forceRedo=forceRedoOutputSomatic))
  if ( class(somatics) == 'try-error' ) {
    catLog('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.\n')
    warning('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.')
  }
  
  #do multi-sample heatmaps and frequency progression
  source('makeSNPprogressionPlots.R')
  progression = try(makeSNPprogressionPlots(variants, timeSeries, normals, plotDirectory, cpus=cpus, forceRedo=forceRedoSNPprogression))
  if ( class(progression) == 'try-error' ) {
    catLog('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.')
  }

  #call CNVs compared to the normals.
  source('callCNVs.R')
  cnvs =
    try(callCNVs(variants=variants, normalVariants=normalVariants, fitS=fitS,
                 names=names, individuals=individuals, normals=normals, Rdirectory=Rdirectory,
                 genome=genome, cpus=cpus, forceRedoCNV=forceRedoCNV))
  if ( class(cnvs) == 'try-error' ) {
    catLog('Error in callCNVs!\n')
    dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'fitS'=fitS,
                               'names'=names, 'individuals'=individuals, 'normals'=normals, 'Rdirectory'=Rdirectory,
                               'genome'=genome, 'cpus'=cpus, 'forceRedoCNV'=forceRedoCNV))
    stop('Error in callCNVs.')
  }


  #make CNV plots
  source('makeCNVplots.R')
  cnvplot = try(makeCNVplots(cnvs, plotDirectory=plotDirectory, genome, forceRedoCNVplots=forceRedoCNVplots))
  if ( class(cnvplot) == 'try-error' ) {
    catLog('Error in makeCNVplots! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeCNVplots! Continuing, but these plots are kindof useful.')
  }

  #make summary plots
  source('makeSummaryPlot.R')
  summary = try(makeSummaryPlot(variants, cnvs, normals, individuals, timePoints, plotDirectory,
    genome, cpus, forceRedo=forceRedoSummary))
  if ( class(summary) == 'try-error' ) {
    catLog('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.')
  }

  #combine SNPs and CNVs into stories of subclones.
  source('getStories.R')
  stories = try(getStories(variants=variants, normalVariants=normalVariants, cnvs=cnvs, timeSeries=timeSeries, normals=normals, Rdirectory=Rdirectory,
    plotDirectory=plotDirectory, cpus=cpus, forceRedo=forceRedoStories))
  if ( class(stories) == 'try-error' ) {
    catLog('Error in getStories!\n')
    dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'cnvs'=cnvs,
                               'timeSeries'=timeSeries, 'Rdirectory'=Rdirectory,
                               'plotDirectory'=plotDirectory, 'cpus'=cpus, 'forceRedo'=forceRedoStories))
    stop('Error in getStories!')
  }
  
  source('makeRiverPlots.R')
  river = try(makeRiverPlots(stories, variants, genome, plotDirectory, forceRedo=forceRedoRiver))
  if ( class(river) == 'try-error' ) {
    catLog('Error in makeRiverPlots!\n')
    warning('Error in makeRiverPlots!')
  }

  catLog('Run done! Have fun with the output! :)\n\n')
  
  return(list('fit'=fitS, 'variants'=variants, 'normalVariants'=normalVariants, 'cnvs'=cnvs, 'stories'=stories))
}

#' Loads saved data
#'
#' This function returns the data produced from an analyse() run.
#' @param Rdirectory A character string pointing to the Rdirectory of the analyse() run.
#' @keywords load saved data
#' @export
#' @examples
#' loadData()

loadData = function(Rdirectory) {
  saveFiles = list.files(Rdirectory, pattern = '*.Rdata', full.names=T)
  names = gsub('.Rdata$', '', basename(saveFiles))
  names(names) = saveFiles
  cat('Loading..')
  for ( file in saveFiles ) {
    cat(gsub('.Rdata$', '', basename(file)), '..', sep='')
    names[file] = load(file=file)
  }
  cat('done.\n')
  ret = lapply(names, function(name) get(name))
  names(ret) = gsub('.Rdata$', '', basename(saveFiles))

  return(ret)
}

#' Loads methods.
#'
#' This function loads the analysis functions used in analyse().
#' @keywords load methods analyse
#' @export
#' @examples
#' loadMethods()

loadMethods = function() {
  assign('catLog', function(...) cat(...), envir = .GlobalEnv)
  source('debug.R')
  source('importCaptureRegions.R')
  source('importSampleMetaData.R')
  source('runDE.R')
  source('XRank.R')
  source('makeFitPlots.R')
  source('getVariants.R')
  source('matchFlagVariants.R')
  source('makeScatterPlots.R')
  source('outputNewVariants.R')
  source('outputSomaticVariants.R')
  source('makeSNPprogressionPlots.R')
  source('callCNVs.R')
  source('makeCNVplots.R')
  source('makeSummaryPlot.R')
  source('getStories.R')
  source('makeRiverPlots.R')
}

#' returns input that uses saved data if present.
#'
#' This function returns the input 'forceRedo' for analyse(), so that saved data from previous runs is used if present.
#' @keywords forceRedo
#' @export
#' @examples
#' loadData()

forceRedoNothing = function() list(
  'forceRedoCount'=F,
  'forceRedoNormalCount'=F,
  'forceRedoFit'=F,
  'forceRedoVolcanoes'=F,
  'forceRedoDifferentRegions'=F,
  'forceRedoSNPs'=F,
  'forceRedoVariants'=F,
  'forceRedoNormalSNPs'=F,
  'forceRedoNormalVariants'=F,
  'forceRedoMatchFlag'=F,
  'forceRedoScatters'=F,
  'forceRedoOutputSomatic'=F,
  'forceRedoNewVariants'=F,
  'forceRedoSNPprogression'=F,
  'forceRedoCNV'=F,
  'forceRedoCNVplots'=F,
  'forceRedoSummary'=F,
  'forceRedoStories'=F,
  'forceRedoRiver'=F)

#' returns input that doesn't use saved data.
#'
#' This function returns the input 'forceRedo' for analyse(), so that saved data from previous runs is never used, and any present saved data is overwritten.
#' @keywords forceRedo
#' @export
#' @examples
#' loadData()

forceRedoEverything = function() list(
  'forceRedoCount'=T,
  'forceRedoNormalCount'=T,
  'forceRedoFit'=T,
  'forceRedoVolcanoes'=T,
  'forceRedoDifferentRegions'=T,
  'forceRedoSNPs'=T,
  'forceRedoVariants'=T,
  'forceRedoNormalSNPs'=T,
  'forceRedoNormalVariants'=T,
  'forceRedoMatchFlag'=T,
  'forceRedoScatters'=T,
  'forceRedoOutputSomatic'=T,
  'forceRedoNewVariants'=T,
  'forceRedoSNPprogression'=T,
  'forceRedoCNV'=T,
  'forceRedoCNVplots'=T,
  'forceRedoSummary'=T,
  'forceRedoStories'=T,
  'forceRedoRiver'=T)
