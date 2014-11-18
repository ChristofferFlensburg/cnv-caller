
#plots the CNV calls and clonalities, and the underlying coverage and frequencies that the calls are based on.
makeCNVplots = function(cnvs, plotDirectory, genome='hg19', v='', forceRedoCNVplots=F) {
  CNVplotDirectory = paste0(plotDirectory, '/CNV/')
  if ( !file.exists(CNVplotDirectory) ) dir.create(CNVplotDirectory)
  for ( name in names(cnvs) ) {
    filename = paste0(CNVplotDirectory, name, '.pdf')
    if ( !file.exists(filename) | forceRedoCNVplots ) {
      catLog('Plotting CNVs to ', filename, '.\n', sep='')
      pdf(filename, width=20, height=10, compress=T)
      plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome)
      plotCR(cnvs[[name]]$clusters, genome=genome)
      for( chr in names(chrLengths(genome)) ) {
        plotCR(cnvs[[name]]$clusters, errorBars=F, genome=genome, chr=chr, alpha=0)
        plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome, chr=chr, alpha=0.3, add=T)
        plotCR(cnvs[[name]]$clusters, errorBars=T, genome=genome, chr=chr, add=T)
      }
      dev.off()
    }
  }
}




#The plotting function for CNV calls.
plotCR = function(cR, showClonality=T, errorBars=T, chr='all', genome='hg19', alpha=1, add=F, ...) {
  showClonality = showClonality & 'subclonality' %in% names(cR)
  if ( nrow(cR) == 0 ) return()
  if ( chr != 'all' ) {
    chrL = chrLengths(genome=genome)
    xlim = c(cumsum(chrL)[chr]-chrL[chr], cumsum(chrL)[chr])
    cR = cR[cR$x1 > xlim[1] & cR$x2 < xlim[2],]
    if ( nrow(cR) == 0 ) return()
  }
  else xlim=c(min(cR$x1), max(cR$x2))
  xlimMar = xlim
  xlimMar[1] = xlim[1] - (xlim[2]-xlim[1])*0.01
  xlimMar[2] = xlim[2] + (xlim[2]-xlim[1])*0.01

  ylim = c(-2.3, 1.1)
  if ( !showClonality ) ylim = c(-1.1, 1.1)
  if ( !add ) {
    par(oma=rep(0, 4))
    par(mar=c(0, 0, 0, 0))
    plot(0, xlim=xlimMar, ylim=ylim, type='n', xlab='', ylab='', axes=F, ...)
    segments(xlim[1], c(0.1,0.6,0.6+log2(1.5)/2,1.1, -0.1,-1.1+2/3, -0.6, -1.1),
             xlim[2], c(0.1,0.6,0.6+log2(1.5)/2,1.1, -0.1,-1.1+2/3, -0.6, -1.1), lwd=3,
             col=c(rgb(0,0,0,0.3), rgb(0,0.8,0,0.3), rgb(1,0.5,0,0.3), rgb(1,0,0,0.3),
               rgb(0,0.8,0,0.3), rgb(1,0.5,0,0.3), rgb(1,0,0,0.3), rgb(0,0,0, 0.3)))
    text(xlim[2] + (xlim[2]-xlim[1])*0.005, c(0.1,0.6,0.6+log2(1.5)/2,1.1, -0.1,-1.1+2/3, -0.6, -1.1),
         c('A', 'AB', 'ABB', 'ABBB', 'AB', 'ABB', 'ABBB', 'A'), adj=0, cex=0.7,
         col=c(rgb(0,0,0,0.3), rgb(0,0.8,0,0.3), rgb(1,0.5,0,0.3), rgb(1,0,0,0.3),
           rgb(0,0.8,0,0.3), rgb(1,0.5,0,0.3), rgb(1,0,0,0.3), rgb(0,0,0, 0.3)))
    text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.032, 0.6, srt=90, 'LFC coverage', cex=0.8)
    text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.032, -0.6, srt=90, 'MAF', cex=0.8)
    text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.032, -1.8, srt=90, 'clonality', cex=0.8)
    addChromosomeLines(ylim=c(-2.3, 1.15), col=mcri('green', 0.6), lwd=1, genome=genome)
    segments(2*xlim[1]-xlim[2], c(0, -1.2), 2*xlim[2]-xlim[1], c(0, -1.2), lwd=5)
    axis(side=2, at=0.1 + (0:4)/4, labels=c(-1, -0.5, 0, 0.5, 1), line=-3)
    axis(side=2, at=-1.1 + (0:5)/5, labels=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), line=-3)
    axis(side=2, at=-2.3 + (0:5)/5, labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line=-3)
  }
  
  x = (cR$x1 + cR$x2)/2
  y = 0.6+cR$M/2
  col = rgb(pmin(1,pmax(0,-2*cR$M)),0,pmin(1,pmax(0,2*cR$M)), alpha)
  lwd = sqrt(0.05/cR$width)
  points(x, pmax(0, y), cex=lwd, col=col, pch=16)
  if ( errorBars ) segments(x, pmax(0,y-cR$width/2), x, pmax(0, y+cR$width/2), lwd=lwd, col=col)
  segments(cR$x1, pmax(0,y), cR$x2, pmax(0, y), lwd=lwd, col=col)
  tooHigh = y > 1.2
  tooLow = y < 0
  if ( any(tooHigh) ) {
    arrows(x[tooHigh], 1.12, x[tooHigh], 1.18, lwd=2, length=0.1, col=rgb(0,0,1,alpha))
    text(x[tooHigh], 1.09, round(cR$M[tooHigh],2), cex=0.8, col=rgb(0,0,1,alpha))
  }
  if ( any(tooLow) ) {
    arrows(x[tooLow], 0.11, x[tooLow], 0.05, lwd=2, length=0.1, col=rgb(1,0,0,alpha))
    text(x[tooLow], 0.14, round(cR$M[tooLow],2), cex=0.8, col=rgb(1,0,0,alpha))
  }

  xf = (cR$x1+cR$x2)[cR$cov > 0]/2
  cRf = cR[cR$cov > 0,]
  f = refUnbias(cRf$var/cRf$cov)
  if ( 'f' %in% colnames(cR) ) f = cR$f[cR$cov > 0]
  yf = -1.1 + 2*f
  ferr = 2*sqrt(cRf$var)/cRf$cov
  cex = pmax(0.2, pmin(3, sqrt(cRf$var/100)))
  #plot the average MAFs, opaqueness from how likely a non-het is.
  pHet = cRf$pHet
  if ( 'call' %in% colnames(cRf) ) pHet[cRf$call == 'AB'] = 1
  col = rgb(pmax(0,pmin(1, 4*(0.5-f))),0,0, (1-pHet)*alpha)
  points(xf, yf, cex=cex, pch=16, col=col)
  if ( errorBars ) segments(xf, pmin(-0.1, yf+ferr), xf, pmax(-1.1, yf-ferr), col=col, lwd=cex)
  segments(cRf$x1, yf, cRf$x2, yf, lwd=cex, col=col)
  #plot f=0.5, opaqueness from how likely a het is.
  col = rgb(0,0,0, pHet*alpha)
  yf = rep(-0.1, length(xf))
  points(xf, yf, cex=cex, pch=16, col=col)
  if ( errorBars ) segments(xf, pmin(-0.1, yf+ferr), xf, pmax(-1.1, yf-ferr), col=col, lwd=cex)
  segments(cRf$x1, yf, cRf$x2, yf, lwd=cex, col=col)

  if ( 'call' %in% colnames(cR) ) {
    called = which(cR$call != 'AB')
    if ( length(called) > 0 ) {
      y = rep(c(-1.05, -0.95, -0.85), length(called))[1:length(called)]
      text(x[called], y-0.05, cR$call[called], cex=0.9, col=mcri('green', alpha))
    }
  }

  if ( showClonality ) {
    subcloneMax = sort(unique(cR$subclonality + cR$subclonalityError), decreasing=T)
    subclone = sort(unique(cR$subclonality), decreasing=T)
    subcloneMin = sort(unique(cR$subclonality - cR$subclonalityError), decreasing=T)
    subcloneSigma = sapply(subclone, function(subF) mean(((cR$clonality-subF)/cR$clonalityError)[cR$subclonality == subF]))
    if ( length(subclone) == 1 ) subcloneCol = subcloneColBG = rgb(0,0,0,alpha)
    else {
      subcloneCol = c(rgb(0,0,0,alpha), randomCol(1:(length(subclone)-1), a=alpha, noBlack=T))
      subcloneColBG = c(rgb(0,0,0,alpha), randomCols(1:(length(subclone)-1),
        pmax(0.02, pmin(0.25, (0.01/pmax(subcloneSigma,1)/(subcloneMax-subcloneMin))[-1]))*alpha, noBlack=T))
    }
    col = sapply(cR$subclonality, function(sub) subcloneCol[which(subclone == sub)])
    rect(xlim[1], subcloneMin - 2.3, xlim[2], pmin(-1.3, subcloneMax -2.3),
         col=subcloneColBG, border=F)
    y = -2.3 + cR$clonality
    error = cR$clonalityError
    points(x, y, pch=16, cex=pmin(2, pmax(0.2, 0.02/error)), col=col)
    segments(x, y-error, x, pmin(-1.3, y+error), lwd=pmax(0.2, pmin(2, 0.05/error)), col=col)
    segments(cR$x1, y, cR$x2, y, lwd=pmax(0.2, pmin(2, 0.05/error)), col=col)
  }
}
