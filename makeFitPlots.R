
require(WriteXLS)


#takes a fit object for the differential coverage between samples and normals, and plots volcano plots
#also outputs differential coverage to an excel file.
makeFitPlots = function(fit, plotDirectory, genome, forceRedoVolcanoes=F, forceRedoDifferentRegions=F) {
  dirname = paste0(plotDirectory, '/volcanoes/')
  if ( !file.exists(dirname) ) dir.create(dirname)

  catLog('Plotting volcanoes to ', dirname, '..', sep='')
  for (col in colnames(fit) ) {
    volcanoFile = paste0(dirname, col, '.jpg')
    if ( !file.exists(volcanoFile) | forceRedoVolcanoes ) {
      catLog(col, '..', sep='')
      jpeg(volcanoFile, height = 10, width = 15, res=300, units='in')
      plotVolcano(fit, coef=col)
      dev.off()
    }
  }
  catLog('done!\n', sep='')


  differentRegionFile = paste0(plotDirectory, '/differentRegionsSamples.xls')
  if ( !file.exists(differentRegionFile) | forceRedoDifferentRegions ) {
    catLog('Writing different regions to ', differentRegionFile, '..', sep='')
    tops = list()
    for (col in colnames(fit) ) {
      catLog(col, '..', sep='')
      fdr = p.adjust(fit$p.value[,col], method='fdr')
      ord = order(fit$XRank[,col])
      if ( length(ord) > 65534 ) {
        ord = ord[1:65534]
        catLog('outputting top 65k DE regions only, for excel...')
      }
      
      tops[[col]] = data.frame(
            'captureRegion'=rownames(fit),
            'chr'=xToChr(fit$x, genome=genome),
            'pos'=xToPos(fit$x, genome=genome),
            'LFC'=fit$coefficients[,col],
            'width'=abs(fit$coefficients/fit$t)[,col],
            'p-value'=fit$p.value[,col],
            'FDR'=fdr,
            'correctedLFC'=fit$best.guess[,col],
            'XRank'=fit$XRank[,col]
            )[ord,]
    }
    WriteXLS('tops', differentRegionFile)
    catLog('done!\n', sep='')
  }
}

#helepr functions converting from genomic coordinates to chr+bp
xToChr = function(x, genome='hg19') {
  chrL = chrLengths(genome)
  ret = rep('', length(x))
  for ( chr in names(chrL) ) {
    ret[x > 0 & x < chrL[chr]] = chr
    x = x - chrL[chr]
  }
  return(ret)
}
xToPos = function(x, genome='hg19') {
  chr = xToChr(x, genome)
  pos = x - cumsum(chrLengths(genome))[chr] + chrLengths(genome)[chr]
  return(pos)
}
