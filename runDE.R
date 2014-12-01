#required packages
require(limma)
require(edgeR)
require(Rsubread)
require(parallel)

#This function take a set of bamfiles, a set of normal bamfiles, and capture regions as input.
#The function runs differential coverage one each sample vs the pool of normals using limma-voom.
#the counts are loess and gc corrected.
runDE = function(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory, normalRdirectory, cpus=1,
  forceRedoFit=F, forceRedoCount=F, forceRedoNormalCount=F) {
  catLog('Starting differential coverage analysis by sample.\n')

  fitPsaveFile = paste0(Rdirectory, '/fitP.Rdata')
  if ( file.exists(fitPsaveFile) & !forceRedoFit & !forceRedoCount & !forceRedoNormalCount ) {
    catLog('Loading saved differential coverage results.\n')
    load(file=fitPsaveFile)
    catLog('Loaded saved differential coverage results of dimension', dim(fitP), '\n')
    return(fitP)
  }

  catLog('Preparing capture regions for featureCounts..')
  captureAnnotation = try(captureRegionToAnnotation(captureRegions))
  if ( class('captureAnnotation') == 'try-error' ) stop('Error in captureRegionToAnnotation.')
  if ( !('GeneID' %in% colnames(captureAnnotation)) ) stop('captureAnnotation does not have a GeneID.')
  catLog('done.\n')
  
  fCsSaveFile = paste0(Rdirectory, '/fCs.Rdata')
  if ( !file.exists(fCsSaveFile) | forceRedoCount ) {
    catLog('Counting reads over capture regions.\n')
    fCs = try(featureCounts(bamFiles, annot.ext=captureAnnotation, useMetaFeatures=T,
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
    catLog('Loaded counts of dimension', dim(fCs$counts), '\n')
  }

  normalFCsSaveFile = paste0(normalRdirectory, '/normalFCs.Rdata')
  if ( !file.exists(normalFCsSaveFile) | forceRedoNormalCount ) {
    catLog('Counting normal reads over capture regions.\n')
    normalFCs = try(featureCounts(externalNormalBams, annot.ext=captureAnnotation, useMetaFeatures=T,
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
    catLog('Loaded normal counts of dimension', dim(normalFCs$counts), '\n')
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
    plotfile = paste0(MAdirectory, col, '.png')
    if ( !file.exists(plotfile) | forceRedoFit ) {
      png(plotfile, height=2000, width=4000, res=144)
      plotMA(counts[,col], normalMean, loess=T, span=0.5, verbose=F, main = paste0(col, ' vs normals (before loess normalisation)'))
      dev.off()
    }
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
  widths = noneg(width(captureRegions))+1
  gc = sapply(genes, function(gene) {
    is = which(regionNames == gene)
    if ( length(is) == 0 ) return(0)
    return(sum(gc[is]*widths[is])/sum(widths[is]))
  })
  correctionFactor = list()
  GCdirectory = paste0(diagnosticPlotsDirectory, '/GCplots/')
  if ( !file.exists(GCdirectory) ) dir.create(GCdirectory)
  for (col in colnames(counts)) {
    catLog(col, '..', sep='')
    lo = loess(cov~gc, data=data.frame(cov=log((1+counts[,col])/(annotation$Length)), gc=gc),
      weights=noneg(annotation$Length-100), span=0.3, family='symmetric', degrees=1)
    los = exp(predict(lo, gc))
    correctionFactor[[col]] = pmin(10, pmax(0.1, (1/los)/mean(1/los)))
    correctionFactor[[col]] = correctionFactor[[col]]*sum(counts[,col])/sum(counts[,col]*correctionFactor[[col]])
    plotfile = paste0(GCdirectory, col, '.png')
    if ( !file.exists(plotfile) | forceRedoFit ) {
      png(plotfile, height=2000, width=4000, res=144)
      plotColourScatter(gc, (counts[,col]/(annotation$Length)), ylim=c(0,2), cex=pmin(3,sqrt(annotation$Length/1000))*1.5, main=col, xlab='GC content', ylab='reads/bp', verbose=F)
      lines((0:100)/100, exp(predict(lo, (0:100)/100)), lwd=5, col=mcri('orange'))
      legend('topleft', 'loess fit', lwd=5, col=mcri('orange'))
      dev.off()
    }
  }
  correctionFactor = do.call(cbind, correctionFactor)
  catLog('done.\n')

  #MA plots after loess and GC correcrion
  diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
  if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
  catLog('Making MA plots of coverage after loess and GC correction in diagnostics directory..')
  normalMean = rowMeans(counts[,group=='normal'])
  MAdirectory = paste0(diagnosticPlotsDirectory, '/MAplotsAfter/')
  if ( !file.exists(MAdirectory) ) dir.create(MAdirectory)
  for ( col in colnames(counts) ) {
    catLog(col, '..', sep='')
    plotfile = paste0(MAdirectory, col, '.png')
    if ( !file.exists(plotfile) | forceRedoFit ) {
      png(plotfile, height=2000, width=4000, res=144)
      plotMA(counts[,col], normalMean, loess=T, span=0.5, verbose=F, main = paste0(col, ' vs normals (before loess normalisation)'))
      dev.off()
    }
  }
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
  
  catLog('Saving fit of dimension', dim(fitP), '..')
  save(fitP, file=fitPsaveFile)
  catLog('done.\n')

  catLog('Returning fit of dimension', dim(fitP), '\n')
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


#helper function that MA-plots two vectors
plotMA = function(x, y, col=mcri('darkblue'), libNorm = F, span=0.2,
  xlab='A = log2(x*y)/2', ylab='M = log2(x/y)', loess=F, cex=0.6, pch=16, verbose=T, ...) {
  Nx = sum(x)
  Ny = sum(y)
  if ( libNorm ) {
    x = x/Nx
    y = y/Ny
  }
  xmin = ymin = 1
  if ( libNorm ) {
    xmin = smear/Nx
    ymin = smear/Ny
  }
  x = xmin*0.2 + noneg(xmin*0.8 + x + pmax(-0.45*xmin, pmin(0.45*xmin, rnorm(length(x), 0, xmin*0.15))))
  y = ymin*0.2 + noneg(ymin*0.8 + y + pmax(-0.45*xmin, pmin(0.45*xmin, rnorm(length(y), 0, ymin*0.15))))
  A = log2(x*y)/2
  M = log2(x/y)
  plot(A, M, cex=cex, pch=pch, xlab=xlab, ylab=ylab, yaxt='n', col=col, ...)
  segments(rep(-30, 5), c(0,-1,1, log2(c(10, 0.1))), rep(30,5), c(0,-1,1, log2(c(10, 0.1))), col=rgb(.5, .5, .5, .2), lwd=2)
    points(A, M, cex=2/3*cex, pch=pch,
           col=mcri('blue', 0.4))
    points(A, M, cex=1/2*cex, pch=pch,
           col=mcri('azure', 0.1))
    points(A, M, cex=1/3*cex, pch=pch,
           col=mcri('green', 0.02))
  axis(2, at=c(log2(0.1), -1, 0,1,log2(10)), labels=c('log2(0.1)', '-1', '0', '1', 'log2(10)'), cex.axis=1)
    
  if ( loess ) {
    if ( verbose ) cat('Calculating loess fit...')
    lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
    if ( verbose ) cat('done.\n')
    As = min(A) + (0:100)/100*(max(A) - min(A))
    lines(As, predict(lo, As), col=mcri('orange'), lwd=6)
  }
}

#helper function that returns the mcri version of the provided colour(s) if available
#Otherwise returns the input. Call without argument to see available colours
mcri = function(col=0, al=1) {
  if ( col[1] == 0 ) {
    cat('Use: mcri(\'colour\'), returning an official MCRI colour.\nAvailable MCRI colours are:\n\ndarkblue\nblue\nlightblue\nazure\ngreen\norange\nviolet\ncyan\nred\nmagenta (aka rose).\n\nReturning default blue.\n')
    return(mcri('blue'))
  }
  if ( length(col) > 1 ) return(sapply(col, function(c) mcri(c, al)))
  if ( is.numeric(col) ) {
    col = (col %% 8) + 1
    if ( col == 1 ) col = 'blue'
    else if ( col == 2 ) col = 'orange'
    else if ( col == 3 ) col = 'green'
    else if ( col == 4 ) col = 'magenta'
    else if ( col == 5 ) col = 'cyan'
    else if ( col == 6 ) col = 'red'
    else if ( col == 7 ) col = 'violet'
    else if ( col == 8 ) col = 'darkblue'
    else col = 'black'
  }
  ret = 0
  if ( col == 'darkblue') ret = rgb(9/255, 47/255, 94/255, al)
  if ( col == 'blue') ret = rgb(0, 83/255, 161/255, al)
  if ( col == 'lightblue') ret = rgb(0, 165/255, 210/255, al)
  if ( col == 'azure') ret = rgb(0, 173/255, 239/255, al)
  if ( col == 'green') ret = rgb(141/255, 198/255, 63/255, al)
  if ( col == 'orange') ret = rgb(244/255, 121/255, 32/255, al)  
  if ( col == 'violet') ret = rgb(122/255, 82/255, 199/255, al)  
  if ( col == 'cyan') ret = rgb(0/255, 183/255, 198/255, al)  
  if ( col == 'red') ret = rgb(192/255, 80/255, 77/255, al)  
  if ( col == 'magenta' | col == 'rose') ret = rgb(236/255, 0/255, 140/255, al)
  if ( ret == 0 ) ret = do.call(rgb, as.list(c(col2rgb(col)/255, al)))
  return(ret)
}


#helper function for plotting
#A pimped up version of the usual boring plot.
#Features include pch=16, mcri colour scheme, correlation in title and highlighting for dense regions.
plotColourScatter = function(x, y, xlab='', ylab='', col=mcri('darkblue'), main='cor',
  add=F, cex=1,verbose=T,...) {
  if ( verbose ) cat('Correlation is ', cor(x,y), '.\n', sep='')
  if ( main == 'cor' ) main = paste('Correlation is', signif(cor(x,y), 2))
  if ( !add ) plot(x, y, cex=cex*0.6, pch=16, xlab=xlab, ylab=ylab,
                   col=col, main=main, ...)
  else points(x, y, cex=cex*0.6, pch=16, col=col, ...)
  points(x, y, cex=cex*0.4, pch=16, col=mcri('blue', 0.4))
  points(x, y, cex=cex*0.3, pch=16, col=mcri('azure', 0.1))
  points(x, y, cex=cex*0.2, pch=16, col=mcri('green', 0.02))
}

#helper function that replaces negative values with 0.
noneg = function(x) return(ifelse(x < 0, 0, x))

#helper wrappers of colSums etc that handle non-matrices.
colsums = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (colSums(mx))
}
rowsums = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (rowSums(mx))
}
colmeans = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (colMeans(mx))
}
rowmeans = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (rowMeans(mx))
}

#normalise all columns to the average column size
libNorm = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  sizes = colSums(mx)
  av = mean(colSums(mx))
  ret = t(t(mx)*av/sizes)
  return(ret)
}
