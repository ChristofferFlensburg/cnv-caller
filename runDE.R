#required packages
require(limma)
require(edgeR)
require(Rsubread)

#This function take a set of bamfiles, a set of normal bamfiles, and capture regions as input.
#The function runs differential coverage one each sample vs the pool of normals using limma-voom.
#the counts are loess and gc corrected.
runDE = function(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory, normalRdirectory, v='', cpus=1,
  forceRedoFit=F, forceRedoCount=F, forceRedoNormalCount=F) {
  catLog('Starting differential coverage analysis by sample.\n')

  fitPsaveFile = paste0(Rdirectory, '/fitP.Rdata')
  if ( file.exists(fitPsaveFile) & !forceRedoFit & !forceRedoCount & !forceRedoNormalCount ) {
    catLog('Loading saved differential coverage results.\n')
    load(file=fitPsaveFile)
    return(fitP)
  }

  catLog('Preparing capture regions for featureCounts..')
  captureAnnotation = try(captureRegionToAnnotation(captureRegions))
  if ( class('captureAnnotation') == 'try-error' ) stop('Error in captureRegionToAnnotation.')
  if ( !('GeneID' %in% colnames(captureAnnotation)) ) stop('captureAnnotation does not have a GeneID.')
  captureAnnotation$GeneID = gsub('.*:', '', captureAnnotation$GeneID)
  catLog('done.\n')
  
  fCsSaveFile = paste0(Rdirectory, '/fCs.Rdata')
  if ( !file.exists(fCsSaveFile) | forceRedoCount ) {
    catLog('Counting reads over capture regions.\n')
    fCs = try(featureCounts(bamFiles, annot.ext=captureAnnotation, useMetaFeatures=F,
      allowMultiOverlap=T, isPairedEnd=T, minMQS=10, nthreads=cpus))
    if ( class(fCs) != 'list' ) {
      catLog('Error in featureCounts.\nInput was\nbamFiles:', bamFiles,
             '\ncaptureAnnotation[1:10,]:', as.matrix(captureAnnotation[1:10,]), '\n')
      stop('Error in featureCounts.')
    }
    catLog('Got a sample count matrix of size', dim(fCs$counts), '\n')
    colnames(fCs$counts) = names
    catLog('Saving sample counts to ', fCsSaveFile, '..', sep='')
    save(fCs, file=fCsSaveFile)
    catLog('done.\n')
  }
  else {
    catLog('Loading counts from file..')
    load(fCsSaveFile)
    catLog('done.\n')
  }

  normalFCsSaveFile = paste0(normalRdirectory, '/normalFCs.Rdata')
  if ( !file.exists(normalFCsSaveFile) | forceRedoNormalCount ) {
    catLog('Counting normal reads over capture regions.\n')
    normalFCs = try(featureCounts(externalNormalBams, annot.ext=captureAnnotation, useMetaFeatures=F,
      allowMultiOverlap=T, isPairedEnd=T, minMQS=10, nthreads=cpus))
    if ( class(normalFCs) != 'list' ) {
      catLog('Error in featureCounts.\nInput was\nexternalNormalBams:', externalNormalBams,
             '\ncaptureAnnotation[1:10,]:', as.matrix(captureAnnotation[1:10,]), '\n')
      stop('Error in featureCounts of normals.')
    }
    catLog('Got a normals count matrix of size', dim(normalFCs$counts), '\n')
    colnames(normalFCs$counts) = names(externalNormalBams)
    catLog('Saving normals counts to ', normalFCsSaveFile, '..', sep='')
    save(normalFCs, file=normalFCsSaveFile)
    catLog('done.\n')
  }
  else {
    catLog('Loading normals counts from file..')
    load(normalFCsSaveFile)
    catLog('done.\n')
  }
  catLog('Merging sample and normals counts..')
  counts = cbind(fCs$counts, normalFCs$counts)
  annotation = fCs$annotation
  catLog('done.\n')

  catLog('Determining sex..')
  xes = grep('^X', annotation$Chr)
  yes = grep('^Y', annotation$Chr)
  if ( length(xes)==0 | length(yes)==0 ) {
    catLog('Capture regions not present in both chromosome X and Y, skip sex matching.\n')
    sex = rep('female', ncol(counts))
  }
  else {
    maleScore =
      colsums(10*libNorm(counts)[yes,])/sum(annotation$Length[yes]+300) -
        colsums(libNorm(counts)[xes,])/sum(annotation$Length[xes]+300)
    sex = ifelse(maleScore > 0, 'male', 'female')
    names(sex) = colnames(counts)
    catLog('done.\nSAMPLE', 'SCORE', 'SEX', '\n', sep='   ')
    for ( i in 1:length(sex) ) catLog(names[i], maleScore[i], sex[i], '\n', sep='   ')
  }

  catLog('Setting up design matrix for linear analysis..')
  group = c(names, rep('normal', length(externalNormalBams)))
  if ( sum(sex[group == 'normal'] == 'male') > 1 & sum(sex[group == 'normal'] == 'female') > 1 ) {
    design = model.matrix(~0+group+sex)
    catLog('Got enough normals of both sexes, will sex-match normals..')
  }
  else {
    design = model.matrix(~0+group)
    catLog('Not enough normals of both sexes, will not sex-match normals..')
  }
  colnames(design) = gsub('^group', '', gsub('^sex', '',colnames(design)))
  contrastList = c(lapply(names, function(name) paste0(name, '-normal')), list(levels=colnames(design)))
  contrasts = do.call(makeContrasts, contrastList)
  catLog('done.\n')
  catLog('Design matrix is\n', colnames(design), '\n', sep='   ')
  for ( row in rownames(design) ) catLog(design[row,], '\n', sep=' ')

  #MA plots
  diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
  if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
  catLog('Making MA plots of coverage in diagnostics directory..')
  normalMean = rowMeans(counts[,group=='normal'])
  MAdirectory = paste0(diagnosticPlotsDirectory, '/MAplots/')
  if ( !file.exists(MAdirectory) ) dir.create(MAdirectory)
  for ( col in colnames(counts) ) {
    catLog(col, '..', sep='')
    png(paste0(MAdirectory, col, '.png'), height=2000, width=4000, res=144)
    plotMA(counts[,col], normalMean, loess=T, span=0.5, verbose=F, main = paste0(col, ' vs normals (before loess normalisation)'))
    dev.off()
  }
  catLog('done.\n')  

  #loess normalise normals, then match the non-normal M-A dependence to the normals.
  #This is not needed for most samples, but will really improve accuracy when there is an M-A bias.
  catLog('Loess normalising counts to normals..')
  counts[,group=='normal'] = loessNormAll(counts[,group=='normal'], span=0.5)
  counts[,group=='normal'] = loessNormAll(counts[,group=='normal'], span=0.5)
  counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal'], counts[,group=='normal'], span=0.5)
  counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal'], counts[,group=='normal'], span=0.5)
  catLog('done.\n')

  #gc correcting
  catLog('Correcting gc bias..')
  genes = rownames(counts)
  regionNames = captureRegions$region
  gc = as.numeric(captureRegions$gc)
  widths = width(captureRegions)
  gc = sapply(genes, function(gene) {
    is = which(regionNames == gene)
    return(sum(gc[is]*widths[is])/sum(widths[is]))
  })
  correctionFactor = list()
  GCdirectory = paste0(diagnosticPlotsDirectory, '/GCplots/')
  if ( !file.exists(GCdirectory) ) dir.create(GCdirectory)
  for (col in colnames(counts)) {
    catLog(col, '..', sep='')
    png(paste0(GCdirectory, col, '.png'), height=2000, width=4000, res=144)
    lo = loess(cov~gc, data=data.frame(cov=log((1+counts[,col])/(annotation$Length)), gc=gc),
      weights=noneg(annotation$Length-100), span=0.3, family='symmetric', degrees=1)
    los = exp(predict(lo, gc))
    correctionFactor[[col]] = pmin(10, pmax(0.1, (1/los)/mean(1/los)))
    correctionFactor[[col]] = correctionFactor[[col]]*sum(counts[,col])/sum(counts[,col]*correctionFactor[[col]])
    plotColourScatter(gc, (counts[,col]/(annotation$Length)), ylim=c(0,2), cex=pmin(3,sqrt(annotation$Length/1000))*1.5, main=col, xlab='GC content', ylab='reads/bp', verbose=F)
    lines((0:100)/100, exp(predict(lo, (0:100)/100)), lwd=5, col=mcri('orange'))
    legend('topleft', 'loess fit', lwd=5, col=mcri('orange'))
    dev.off()
  }
  correctionFactor = do.call(cbind, correctionFactor)
  catLog('done.\n')


  catLog('Running voom..')
  jpeg(paste0(diagnosticPlotsDirectory, '/voomVariance.jpg'), height=1000, width=2000)
  voomWeights = voom(counts, design=design, plot=T)
  dev.off()
  voomCounts = voom(counts*correctionFactor, design=design)
  voomCombined = voomCounts
  voomCombined$weights = voomWeights$weights
  catLog('limma..')
  fitP = lmFit(voomCombined, design=design)
  fitP = eBayes(fitP)
  fitP = contrasts.fit(fitP, contrasts)
  fitP = eBayes(fitP)
  catLog('XRank..')
  fitP = XRank(fitP, plot=F, cpus=cpus, verbose=F)
  fitP$x = annotationToX(annotation, genome=genome)
  x1x2 = annotationToX1X2(annotation, genome=genome)
  fitP$x1 = x1x2[,1]
  fitP$x2 = x1x2[,2]
  fitP$chr = annotationToChr(annotation)
  fitP$longNames = annotation$GeneID
  catLog('done.\n')

  #remember stats from counting
  catLog('Importing stats about feature counts..')
  tots = colSums(fCs$stat[,-1])
  assigned = fCs$stat[1,-1]
  names(tots) = names(assigned) = names
  fitP$totNReads = tots
  fitP$assignedNReads = assigned
  catLog('done.\n')
  
  catLog('Saving fit..')
  save(fitP, file=fitPsaveFile)
  catLog('done.\n')
  return(fitP)
}










#helper function transforming GRanges capture regions to a format that can be inserted into featureCounts.
captureRegionToAnnotation = function(cR) {
  ret = data.frame(GeneID = cR$region, Start = start(cR), End = end(cR), Strand=strand(cR), Chr = seqnames(cR))
  return(ret)
}

#helper functions doing loess normalisation from MA plots.
loessNormAll = function(counts, ...) {
  newCounts = do.call(cbind, lapply(1:ncol(counts), function(col) loessNorm(counts[,col], counts[,-col], ...)[,1]))
  colnames(newCounts) = colnames(counts)
  rownames(newCounts) = rownames(counts)
  return(newCounts)
}
loessNormAllToReference = function(counts, reference, ...) {
  newCounts = do.call(cbind, lapply(1:ncol(counts), function(col) loessNormToReference(counts[,col], reference, ...)))
  colnames(newCounts) = colnames(counts)
  rownames(newCounts) = rownames(counts)
  return(newCounts)
}
loessNorm = function(counts1, counts2, span=0.5) {
  if ( is.matrix(counts1) ) x = 1 + rowSums(counts1)
  else x = 1 + counts1
  if ( is.matrix(counts2) ) y = 1 + rowSums(counts2)
  else y = 1 + counts2
  M = log(x/y)
  A = log(x*y)/2
  lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
  loM = predict(lo, A)
  loM = loM - sum(loM*exp(A))/sum(exp(A))

  counts1 = counts1*exp(-loM/2)
  counts2 = counts2*exp(loM/2)
  return(round(cbind(counts1, counts2)))
}
loessNormToReference = function(counts1, reference, span=0.5) {
  if ( is.matrix(counts1) ) x = 1 + rowSums(counts1)
  else x = 1 + counts1
  if ( is.matrix(reference) ) y = 1 + rowSums(reference)
  else y = 1 + reference
  M = log10(x/y)
  A = log10(x*y)/2
  lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
  loM = predict(lo, A)
  loM = loM - sum(loM*10^(A))/sum(10^(A))

  counts1 = counts1*10^(-loM)
  return(round(counts1))
}

#helper functions that extracts information from annotation objects.
annotationToX = function(annotation, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(chrLengths(genome)), 'outside')
  chr =gsub('chr','',gsub(';.*', '', as.character(annotation$Chr)))
  start = as.numeric(gsub(';.*', '', annotation$Start))
  end =  as.numeric(gsub('.*;', '', annotation$End))
  x = (end + start)/2 + prevChrL[chr]

  return(x)
}
annotationToX1X2 = function(annotation, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(chrLengths(genome)), 'outside')
  chr =gsub('chr','',gsub(';.*', '', as.character(annotation$Chr)))
  start = as.numeric(gsub(';.*', '', annotation$Start))
  end =  as.numeric(gsub('.*;', '', annotation$End))
  x1 = start + prevChrL[chr]
  x2 = end + prevChrL[chr]

  return(data.frame(x1, x2))
}
annotationToChr = function( annotation ) {
  return(gsub('chr','',gsub(';.*', '',annotation$Chr)))
}

#helper functions that returns the lengths of the chromosomes of a genome.
humanAllChrLengths = function() {
lengths = c(249250621, 106433, 547496, 243199373, 198022430, 191154276, 590426, 189789, 191469, 180915260, 171115067, 4622290, 4795371,
4610396, 4683263, 4833398, 4611984, 4928567, 159138663, 182896, 146364022, 38914, 37175, 141213431, 90085, 169874, 187035, 36148,
135534747, 135006516, 40103, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 1680828, 37498, 81310, 174588, 41001,
78077248, 4262, 59128983, 92689, 159169, 63025520, 48129895, 27682, 51304566, 155270560, 59373566, 166566, 186858, 164239, 137718,
172545, 172294, 172149, 161147, 179198, 161802, 155397, 186861, 180455, 179693, 211173, 15008, 128374, 129120, 19913, 43691, 27386,
40652, 45941, 40531, 34474, 41934, 45867, 39939, 33824, 41933, 42152, 43523, 43341, 39929, 36651, 38154, 36422, 39786, 38502, 16571)

names(lengths) = c('1','1_gl000191_random','1_gl000192_random','2','3','4','4_ctg9_hap1','4_gl000193_random','4_gl000194_random',
'5','6','6_apd_hap1','6_cox_hap2','6_dbb_hap3','6_mann_hap4','6_mcf_hap5','6_qbl_hap6','6_ssto_hap7','7','7_gl000195_random',
'8','8_gl000196_random','8_gl000197_random','9','9_gl000198_random','9_gl000199_random','9_gl000200_random','9_gl000201_random',
'10','11','11_gl000202_random','12','13','14','15','16','17','17_ctg5_hap1','17_gl000203_random','17_gl000204_random',
'17_gl000205_random','17_gl000206_random','18','18_gl000207_random','19','19_gl000208_random',
'19_gl000209_random', '20', '21', '21_gl000210_random', '22', 'X', 'Y', 'Un_gl000211', 'Un_gl000212', 'Un_gl000213', 'Un_gl000214',
'Un_gl000215', 'Un_gl000216', 'Un_gl000217', 'Un_gl000218', 'Un_gl000219', 'Un_gl000220', 'Un_gl000221', 'Un_gl000222', 'Un_gl000223',
'Un_gl000224', 'Un_gl000225', 'Un_gl000226', 'Un_gl000227', 'Un_gl000228', 'Un_gl000229', 'Un_gl000230', 'Un_gl000231', 'Un_gl000232',
'Un_gl000233', 'Un_gl000234', 'Un_gl000235', 'Un_gl000236', 'Un_gl000237', 'Un_gl000238', 'Un_gl000239', 'Un_gl000240', 'Un_gl000241',
'Un_gl000242', 'Un_gl000243',  'Un_gl000244', 'Un_gl000245', 'Un_gl000246', 'Un_gl000247', 'Un_gl000248', 'Un_gl000249', 'M')

  return(lengths)
}
humanChrLengths = function() {
  humanAllChrLengths()[as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y', 'M'))]
}
chrLengths = function(genome='hg19') {
  if ( genome == 'hg19' ) return(humanChrLengths())
  else stop('chrLengths doesnt know about genome', genome, '\n')
}
