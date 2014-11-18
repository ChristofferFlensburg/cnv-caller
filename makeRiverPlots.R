
#takes the subclone evolution from the stories and plots rivers, showing which events are aprt of which subclone.
makeRiverPlots = function(stories, variants, plotDirectory, forceRedo=forceRedoRiver) {
  if ( length(stories) == 0 ) return
  riverDirectory = paste0(plotDirectory, '/rivers/')
  if ( !file.exists(riverDirectory) ) dir.create(riverDirectory)
  for ( ts in names(stories) ) {
    file = paste0(riverDirectory, ts, '-river.pdf')
    if ( file.exists(file) & !forceRedo ) next
    catLog('Making riverplot for ', ts, '..', sep='')
    labels = lapply(stories[[ts]]$clusters$storyList, function(rows) storyToLabel(stories[[ts]]$all[rows,], variants$SNPs))
    pdf(file, width=15, height=10)
    plotRiver(stories[[ts]]$cloneTree, stories[[ts]]$clusters$cloneStories, labels)
    plotStories(stories[[ts]]$clusters$cloneStories)
    for ( subclone in names(stories[[ts]]$clusters$storyList) ) {
      plotStories(stories[[ts]]$all[stories[[ts]]$clusters$storyList[[subclone]],])
    }
    dev.off()
    catLog('done.\n')
  }
}

#helper function converting internal event names to informative labels for the plots.
storyToLabel = function(stories, SNPs) {
  call = stories$call
  isSNP = grepl('[0-9]', call)
  ret = ''
  ret[isSNP] = paste0(SNPs[as.character(stories$x1), ]$reference, ' -> ', gsub('^[0-9]*', '', call), ': ',
       SNPs[as.character(stories$x1), ]$inGene, ' (', xToChr(stories$x1), ')')[isSNP]
  dist = stories$x2-stories$x1
  distText = ifelse(dist >= 1e6, paste0(round(dist/1e6), 'Mbp '), ifelse(dist >= 1e3, paste0(round(dist/1e3), 'kbp '), paste0(dist, 'bp ')))
  ret[!isSNP] = paste0(distText, call, ' (', xToChr(stories$x1), ')')[!isSNP]
  return(ret)
}

#main plotting function
plotRiver = function(cloneTree, cloneStories, cloneLabels) {
  names(cloneLabels) = rownames(cloneStories)
  stories = abs(cloneStories$stories)
  rownames(stories) = rownames(cloneStories)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(abs(stories[names(cloneTree[[1]]),i])))
  stories[stories < cloneStories$errors*2] = 0
  stories = t(t(stories)/purity)
  x = 1:ncol(stories)

  par(oma=rep(0, 4))
  par(mar=c(0, 4, 0, 0))
  plot(1, type='n', xlim=c(1, max(x)+ceiling(nrow(cloneStories)/2)), ylim=c(0,1), xaxt='n', frame.plot=F,
       ylab='clonality', xlab = '')
  cloneCols = addSubclone(cloneTree[[1]], stories, ylims = matrix(rep(c(0,1), ncol(stories)), nrow=2), margin=0.02)$usedCols
  for (i in 1:nrow(cloneStories)) {
    clone = rownames(cloneStories)[i]
    xText = max(x) + ceiling(i/2) - 0.9
    y0 = i/2 - floor(i/2)
    clone = names(cloneCols)[i]
    muts = cloneLabels[[clone]]
    if ( length(muts) > 15 ) muts = c(muts[1:15], 'and more...')
    col = cloneCols[i]
    text(rep(xText, length(muts)), y0 + 0.5 - (1:length(muts))/33, muts, col=col, adj=0, cex=0.9)
  }

  segments(1:ncol(stories), 0.02, 1:ncol(stories), 1, lwd=5, col=rgb(0.7, 0.7, 0.7, 0.3))
  segments(1:ncol(stories), 0.02, 1:ncol(stories), 1, lwd=2, col=rgb(0.3, 0.3, 0.3, 0.3))
  text(1:ncol(stories), -0.02, colnames(stories))

  par(oma=rep(0, 4))
  par(mar=rep(4, 4))
}


#plots a set of parallel disjoint clones, and recurs to each clones subclones to be plotted on top.
addSubclone = function(cT, stories, ylims, colourPool=c(), margin=0.02) {
  if ( length(colourPool) == 0 ) colourPool = mcri(c('darkblue', 'blue', 'green', 'red', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'black', 'grey'))
  subClones = names(cT)
  subStories = stories[names(cT),,drop=F]
  maxSize = ylims[2,]-ylims[1,]
  cloneSum = margin*maxSize+colsums(t(margin*maxSize + t(subStories)))
  cloneSum[cloneSum == 0] = 1   #this happens if all subclones are 0 at a sample. this avoids NaNs.
  norm = pmin(1, maxSize/cloneSum)
  base = ylims[1,]
  usedCols = c()
  for ( subClone in subClones ) {
    subStory = subStories[subClone,]*norm
    range = rbind(base+margin*maxSize*norm, base + margin*maxSize*norm + subStory)
    range[2,] = pmax(base+margin*maxSize*norm, range[2,])
    subrange = rbind(range[1,], range[2]-margin*maxSize*norm)
    addStream(range, col=colourPool[1])
    usedCols = c(usedCols, colourPool[1])
    names(usedCols)[length(usedCols)] = subClone
    colourPool = colourPool[-1]
    if ( length(cT[[subClone]]) > 0 ) {
      out = addSubclone(cT[[subClone]], stories, range, colourPool, margin)
      usedCols = c(usedCols, out$usedCols)
      colourPool = out$colourPool
    }
    base = range[2,]
  }
  return(list('colourPool'=colourPool, usedCols=usedCols))
}

#plots the stream corresponding to a clone
addStream = function(ylims, col='grey') {
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
    addStreamSegment(x1, x2, y1l, y1h, y2l, y2h, range=range, col=col)
  }
}

#helper third degree polynomial
third = function(x, x0, a, b, c) a + b*(x-x0) + c*(x-x0)^3

#plots a smooth stream segment.
addStreamSegment = function(x1, x2, y1low, y1high, y2low, y2high, range=c(0,1), col, parts = 100) {
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
  polygon(c(x, rev(x)), c(yhigh, rev(pmin(ylow, yhigh))), col=col, border=NA)
}
