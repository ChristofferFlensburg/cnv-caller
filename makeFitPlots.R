
require(WriteXLS)


#takes a fit object for the differential coverage between samples and normals, and plots volcano plots
#also outputs differential coverage to an excel file.
makeFitPlots = function(fit, plotDirectory, v, forceRedoVolcanoes=F, forceRedoDifferentRegions=F) {
  volcanoFile = paste0(plotDirectory, '/volcanoSamples.pdf')
  if ( !file.exists(volcanoFile) | forceRedoVolcanoes ) {
    catLog('Plotting volcanoes to ', volcanoFile, '..', sep='')
    pdf(volcanoFile, height = 10, width = 15, compress=T)
    for (col in colnames(fit) ) {
      catLog(col, '..', sep='')
      plotVolcano(fit, coef=col)
    }
    dev.off()
    catLog('done!\n', sep='')
  }
  else catLog('Volcano plots already present in ', volcanoFile, '\n', sep='')


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
            'chr'=xToChr(fit$x),
            'pos'=xToPos(fit$x),
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
