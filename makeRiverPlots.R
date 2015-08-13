
#takes the subclone evolution from the stories and plots rivers, showing which events are aprt of which subclone.
makeRiverPlots = function(stories, variants, genome, cpus=1, plotDirectory, forceRedo=forceRedoRiver) {
  if ( length(stories) == 0 ) return
  riverDirectory = paste0(plotDirectory, '/rivers/')
  if ( !file.exists(riverDirectory) ) dir.create(riverDirectory)
  for ( ts in names(stories) ) {
    file = paste0(riverDirectory, ts, '-river.pdf')
    if ( file.exists(file) & !forceRedo ) next
    catLog('Making riverplot for ', ts, '..', sep='')
    labels = lapply(stories[[ts]]$clusters$storyList, function(rows) storyToLabel(stories[[ts]]$all[rows,], variants$SNPs, genome=genome))
    pdf(file, width=15, height=10)
    cloneCols = plotRiver(stories[[ts]]$cloneTree, stories[[ts]]$clusters$cloneStories, labels)
    plotStories(stories[[ts]]$clusters$cloneStories, variants$SNPs, col=cloneCols, genome=genome)
    heatmapStories(stories[[ts]]$all, stories[[ts]]$clusters$storyList, variants$SNPs, labels=do.call(c, labels), genome=genome)
    for ( subclone in names(stories[[ts]]$clusters$storyList) ) {
      i = which(names(stories[[ts]]$clusters$storyList) == subclone)
      plotStories(stories[[ts]]$all[stories[[ts]]$clusters$storyList[[subclone]],], variants$SNPs, alpha=0.2)
      plotStories(stories[[ts]]$clusters$cloneStories[i,], variants$SNPs, add=T, col=cloneCols[i])
    }
    dev.off()
    catLog('done.\n')
    catLog('Outputting data on stories ', ts, '..', sep='')
    excelFile = paste0(riverDirectory, ts, '-river.xls')
    output = do.call(rbind, lapply(stories[[ts]]$clusters$storyList, function(sL) stories[[ts]]$all[sL,]))
    clone = do.call(c, lapply(1:length(stories[[ts]]$clusters$storyList), function(i) rep(i, length(stories[[ts]]$clusters$storyList[[i]]))))
    output$clone = clone
    chr = xToChr(output$x1, genome)
    start = xToPos(output$x1, genome)
    end = xToPos(output$x2, genome)
    label = unlist(labels)
    clonality = as.data.frame(output$stories)
    error = as.data.frame(output$errors)
    output = cbind(chr = chr, start=start, end=end, name=label, clone=output$clone, clonality=clonality, error=error)
    WriteXLS('output', excelFile)
    catLog('done.\n')
  }
}

#helper function converting internal event names to informative labels for the plots.
storyToLabel = function(stories, SNPs, genome) {
  call = stories$call
  isSNP = grepl('[0-9]', call)
  ret = ''
  SNPs = SNPs[SNPs$x %in% stories$x1,]
  SNPs = SNPs[as.character(stories$x1[isSNP]),]
  ret[isSNP] = paste0(SNPs$reference, ' -> ', gsub('^[0-9]*', '', call[isSNP]), ': ',
       SNPs$inGene, ' (', xToChr(stories$x1[isSNP], genome=genome), ')')
  dist = stories$x2-stories$x1
  distText = ifelse(dist >= 1e6, paste0(round(dist/1e6), 'Mbp '), ifelse(dist >= 1e3, paste0(round(dist/1e3), 'kbp '), paste0(dist, 'bp ')))
  ret[!isSNP & !is.na(dist)] = paste0(distText, call, ' (', xToChr(stories$x1, genome=genome), ')')[!isSNP & !is.na(dist)]
  ret[!isSNP & is.na(dist)] = call[!isSNP & is.na(dist)]
  return(ret)
}

#main plotting function
plotRiver = function(cloneTree, cloneStories, cloneLabels, normalise=T, xlim='default', ylim='default', labels=T, setPar=T) {
  names(cloneLabels) = rownames(cloneStories)
  stories = abs(cloneStories$stories)
  rownames(stories) = rownames(cloneStories)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(abs(stories[names(cloneTree[[1]]),i])))
  leadingNormal = F
  if ( any(purity==0) ) {
    leadingNormal = T
    purity[purity==0] = 1
  }
  stories[stories < cloneStories$errors*2] = 0
  if ( normalise ) stories = t(t(stories)/purity)
  if ( !leadingNormal ) {
    stories = cbind(rep(0, nrow(stories)), stories)
    colnames(stories)[1] = 'germline'
    if ( any(rownames(stories)=='germline') ) stories['germline','germline'] = 1
  }
  x = 1:ncol(stories)

  if ( setPar ) {
    par(oma=rep(0, 4))
    par(mar=c(0, 4, 0, 0))
  }
  if ( !labels & xlim[1] == 'default' ) xlim = c(1, max(x))
  if ( labels & xlim[1] == 'default' ) xlim = c(1, max(x)+ceiling(nrow(cloneStories)/2))
  if ( ylim[1] == 'default' ) ylim = c(-0.02,1)
  plot(1, type='n', xlim=xlim, ylim=ylim, xaxt='n', frame.plot=F,
       ylab='clonality', xlab = '')
  cloneCols = addSubclone(cloneTree[[1]], stories, ylims = matrix(rep(c(0,1), ncol(stories)), nrow=2), margin=0.02)$usedCols
  for (i in 1:nrow(cloneStories)) {
    clone = rownames(cloneStories)[i]
    xText = max(x) + ceiling(i/2) - 0.9
    y0 = i/2 - floor(i/2)
    if ( labels ) {
      clone = names(cloneCols)[i]
      muts = cloneLabels[[clone]]
      if ( length(muts) > 15 ) muts = c(muts[1:15], 'and more...')
      col = cloneCols[i]
      text(rep(xText, length(muts)), y0 + 0.5 - (1:length(muts))/33, muts, col=col, adj=0, cex=0.9)
    }
  }

  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=5, col=rgb(0.7, 0.7, 0.7, 0.3))
  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=2, col=rgb(0.3, 0.3, 0.3, 0.3))
  text(1:ncol(stories), -0.02, colnames(stories), srt=20, cex=0.9)

  par(oma=rep(0, 4))
  par(mar=rep(4, 4))

  return(cloneCols)
}


#plots a set of parallel disjoint clones, and recurs to each clones subclones to be plotted on top.
addSubclone = function(cT, stories, ylims, colourPool=c(), margin=0.02, preNorm=1) {
  if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
  subClones = names(cT)
  subStories = stories[names(cT),,drop=F]
  maxSize = ylims[2,]-ylims[1,]
  cloneSum = margin*maxSize+colsums(t(margin*maxSize + t(subStories)))
  unitarityViolating = any(preNorm*colsums(subStories)*0.9 > maxSize)
  cloneSum[cloneSum == 0] = 1   #this happens if all subclones are 0 at a sample. this avoids NaNs.
  norm = pmin(1, maxSize/cloneSum)
  base = ylims[1,]
  usedCols = c()
  for ( subClone in subClones ) {
    subStory = subStories[subClone,]*norm
    range = rbind(base+margin*maxSize*norm, base + margin*maxSize*norm + subStory)
    range[2,] = pmax(base+margin*maxSize*norm, range[2,])
    subrange = rbind(range[1,], range[2]-margin*maxSize*norm)
    addStream(range, col=colourPool[1], violating=unitarityViolating)
    usedCols = c(usedCols, colourPool[1])
    names(usedCols)[length(usedCols)] = subClone
    colourPool = colourPool[-1]
      if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
    if ( length(cT[[subClone]]) > 0 ) {
      out = addSubclone(cT[[subClone]], stories, range, colourPool, margin, preNorm=norm)
      usedCols = c(usedCols, out$usedCols)
      colourPool = out$colourPool
    }
    base = range[2,]
  }
  return(list('colourPool'=colourPool, usedCols=usedCols))
}

#plots the stream corresponding to a clone
addStream = function(ylims, col='grey', violating=F) {
  for ( sample in 2:ncol(ylims) ) {
    x1 = sample-1
    x2 = sample
    y1h = ylims[2, sample-1]
    y1l = ylims[1, sample-1]
    y2h = ylims[2, sample]
    y2l = ylims[1, sample]
    range = c(0,1)
    if ( y1h == y1l & y2h > y2l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) < sample) )  {
      range[1] = 1-sqrt(y2h-y2l)
    }
    if ( y2h == y2l & y1h > y1l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) > sample)) {
      range[2] = sqrt(y1h-y1l)
    }
    addStreamSegment(x1, x2, y1l, y1h, y2l, y2h, range=range, col=col,violating=violating)
  }
}

#helper third degree polynomial
third = function(x, x0, a, b, c) a + b*(x-x0) + c*(x-x0)^3

#plots a smooth stream segment.
addStreamSegment = function(x1, x2, y1low, y1high, y2low, y2high, range=c(0,1), col, parts = 100, violating=F) {
  sh = (y1high-y2high)/2
  sl = (y1low-y2low)/2
  z = (x1-x2)/2

  ah = (y1high+y2high)/2
  bh = 3*sh/(2*z)
  ch = -sh/(2*z*z*z)
  al = (y1low+y2low)/2
  bl = 3*sl/(2*z)
  cl = -sl/(2*z*z*z)
  x0 = (x1+x2)/2

  x = seq(from=x1, to=x2, length.out=parts)
  yhigh = third(x, x0, ah, bh, ch)
  ylow = third(x, x0, al, bl, cl)

  if ( range[1] == 0 & range[2] < 1 ) {
    xnorm = (x - x1)/(x2-x1)
    r = range[2]
    temp = ylow + (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)
    yhigh = yhigh - (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)
    ylow = ifelse(xnorm < r, temp, yhigh)
  }
  if ( range[1] > 0 & range[2] == 1 ) {
    xnorm = (x2 - x)/(x2-x1)
    r = 1-range[1]
    temp = ylow + (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)
    yhigh = yhigh - (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)
    ylow = ifelse(xnorm < r, temp, yhigh)
  }
  if ( !violating )
    polygon(c(x, rev(x)), c(yhigh, rev(pmin(ylow, yhigh))), col=col, border=NA)
  if ( violating )
    polygon(c(x, rev(x)), c(yhigh, rev(pmin(ylow, yhigh))), col=col, border=NA, density=10, angle=45, lwd=5)
}


plotStories = function(stories, SNPs, col='default', lty='default', add=F, alpha=1, xlab='sample', ylab='clonality', lwd='default', errorBars=T, setPar=T, legend=T, labels=T, xlim='default', genome='hg19', xSpread=0.25,...) {
  names = rownames(stories)
  clon = stories$stories
  ce = stories$errors

  Nsample = ncol(clon)
  Nmut = nrow(clon)

  if ( col[1] == 'default' | length(col) != Nmut ) {
    segcol = randomCols(row(clon)[,2:Nsample], a=alpha)
    errcol = randomCols(row(clon), a=alpha)
  }
  else if ( length(col) == Nmut ) {
    if ( all(rownames(clon) %in% names(col)))
      col = col[rownames(clon)]
    segcol = rep(col, Nsample-1)
    errcol = rep(col, Nsample)
  }

  if ( lty[1] == 'default' | length(lty) != Nmut ) {
    iMx = floor(1+(row(clon)-1)/8)
    seglty = randomLtys(iMx[,2:Nsample])
    errlty = randomLtys(iMx)
  }
  else if ( length(lty) == Nmut ) {
    seglty = rep(lty, Nsample-1)
    errlty = rep(lty, Nsample)
  }

  
  if ( !add ) {
    if ( setPar ) {
      par(oma=rep(0, 4))
      par(mar=c(0, 4, 0, 0))
    }
    if ( xlim[1] == 'default' ) xlim = c(1, Nsample*1.3)
    plot(0,0, type='n', xlim=xlim, ylim=c(-.02,1), xlab=xlab, ylab=ylab, xaxt='n', frame.plot=F, ...)
  }
  if (lwd[1] == 'default' ) lwd = 1/(sqrt(0.1^2+ce[,1:(Nsample-1)]^2 + ce[,2:Nsample]^2)/0.2)^2
  segments(col(clon)[,1:(Nsample-1)]+(row(clon)[,1:(Nsample-1)] - 0.5 - Nmut/2)/Nmut*xSpread, clon[,1:(Nsample-1)],
           col(clon)[,2:Nsample    ]+(row(clon)[,2:Nsample    ] - 0.5 - Nmut/2)/Nmut*xSpread, clon[,2:Nsample    ],
           col=segcol, lwd=lwd, lty=seglty)
  if ( errorBars ) {
    if (lwd[1] == 'default' ) lwd = 1/(sqrt(0.1^2+ce^2)/0.2)^2
    segments(col(clon)+(row(clon) - 0.5 - Nmut/2)/Nmut*xSpread, noneg(clon - ce),
             col(clon)+(row(clon) - 0.5 - Nmut/2)/Nmut*xSpread, clon + ce,
             col=errcol, lwd=lwd, lty=errlty)
  }

  if ( !add ) {
    if ( legend ) {
      lbls = storyToLabel(stories, SNPs, genome)
      legCex = pmin(1, pmax(0.5, 45/length(lbls)))
      legend('topright', lbls, lwd=2, col=errcol[1:nrow(stories)], lty=errlty[1:nrow(stories)], cex= legCex)
      
    }
    if ( labels ) text(1:Nsample, -0.02, colnames(stories$stories), srt=20, cex=0.9)
    if ( setPar ) {
      par(oma=rep(0, 4))
      par(mar=rep(4, 4))
    }
  }
}

heatmapStories = function(stories, storyList, SNPs, labels=NA, genome='hg19') {
  stories = stories[do.call(c, storyList),]
  clone = do.call(c, lapply(1:length(storyList), function(i) rep(i, length(storyList[[i]]))))
  sideCol = randomCols(clone)
  
  clon = stories$stories

  if ( is.na(labels)[1] )
    labels = storyToLabel(stories, SNPs, genome)
  rownames(clon) = labels


  if ( nrow(clon) < 1000 )
    makeHeatmap(clon, RowSideColors=sideCol, Colv=NA, label='clonality')
  else {
    catLog('Too many stories for the built-in heatmap clustering, using default row ordering.\n')
    makeHeatmap(clon, RowSideColors=sideCol, Colv=NA, Rowv=NA)
  }
}


makeHeatmap = function(mx, nCol=200, col='default', maxVal='default', minVal='default', scale='none', label='', ...) {
  if ( maxVal == 'default' ) maxVal = max(mx, na.rm=T)
  if ( minVal == 'default' ) minVal = min(mx, na.rm=T)
  if ( col[1] == 'default' )
    col = colourGradient(cols=mcri(c('black', 'grey', 'cyan', 'blue', 'red')),
      anchors=c(0, 0.02, 0.2, 0.5, 1), steps=nCol)
  ret = heatmap(mx, col=col, scale=scale, ...)
  
  barXmax = par('usr')[1]*0.88+par('usr')[2]*0.12
  barXmin = par('usr')[1]*0.92+par('usr')[2]*0.08
  barYmin = par('usr')[3]*0.8+par('usr')[4]*0.2
  barYmax = par('usr')[3]*0.2+par('usr')[4]*0.8
  lowY = barYmin + (0:(nCol-1))*(barYmax-barYmin)/nCol
  highY = barYmin + (1:nCol)*(barYmax-barYmin)/nCol
  lowX = rep(barXmin, nCol)
  highX = rep(barXmax, nCol)
  barCols = col
  rect(lowX, lowY, highX, highY, col=barCols, border=F)
  segments(rep(barXmin, 3), barYmin + c(0.002, 0.5, 0.998)*(barYmax-barYmin),
           rep(barXmin-(barXmax-barXmin)*0.1, 3), barYmin + c(0.002, 0.5, 0.998)*(barYmax-barYmin),
           lwd=2, col=barCols[c(1, round(nCol/2), nCol)])
  minDist = minVal
  maxDist = maxVal
  text(rep(barXmin-(barXmax-barXmin)*0.2, 4), barYmin + c(0.003, 0.5, 0.997, 1.07)*(barYmax-barYmin),
       c(round(c(minDist, (minDist+maxDist)/2, maxDist), 2), label), adj=c(1, 0.5))

  invisible(ret)
}


#This function takes two colours (in a format that can be handled by col2rgb, such as "red", or rgb(0.1, 0.2, 0.3))
#and two weights and returns a weighted average between the two colours.
#Mainly a helper function for colourGradient, but can potentially find uses on itself.
combineColours = combineColors = function (col1, col2, w1=0.5, w2=0.5) {
  rgb1 = col2rgb(col1)
  rgb2 = col2rgb(col2)
  if ( w1 == Inf & w2 == Inf ) {
    w1 = 0.5
    w2 = 0.5
  }
  if ( w1 == Inf ) {
    w1 = 1
    w2 = 0
  }
  if ( w2 == Inf ) {
    w1 = 0
    w2 = 1
  }
  combinedRgb = (rgb1*w1 + rgb2*w2)/(w1+w2)
  combinedColour = do.call(rgb, as.list(pmin(1, pmax(0, combinedRgb/255))))
  return(combinedColour)
}
#Takes a vector of colours (in a format that can be handled by col2rgb, such as "red", or rgb(0.1, 0.2, 0.3))
#and an optional sorted vector of anchor points between 0 and 1. Returns a vector of colours of length @steps
#that gradually goes through the colours in the vector, hitting each colour at  fraction through the vector
#set by the anchor points. Defaults to a 100-step vector from blue through white to red.
colourGradient = colorGradient = function(cols=mcri(c('red', 'orange', 'white', 'cyan', 'blue')), steps=100, anchors=(1:length(cols)-1)/(length(cols)-1) ) {
  if ( length(anchors) != length(cols) ) {
    warning(paste0('colourGradient: Length of cols and anchors has to be the same. They are ', length(cols), ' and ', length(anchors), '. Returning default colour gradient.'))
  }
  N = length(cols)
  x = (0:(steps-1))/(steps-1)
  if ( any(anchors != sort(anchors)) ) warning('colourGradient: anchors not sorted. Sorting.')
  ord = order(anchors)
  anchors = anchors[ord]
  cols = cols[ord]
  anchors = c(-Inf, anchors, Inf)
  cols = c(cols[1], cols, cols[length(cols)])
  col2i = sapply(x, function(X) which(anchors > X)[1])
  col1i = col2i-1
  col1 = cols[col1i]
  col2 = cols[col2i]
  w1 = 1/abs(x-anchors[col1i])
  w2 = 1/abs(x-anchors[col2i])
  gradient = sapply(1:length(col1), function(i) combineColours(col1[i], col2[i], w1[i], w2[i]))
  names(gradient) = x
  return(gradient)
}
