


#Takes variants and plots the frequencies of any two samples from the same individual against each other.
#Goes to some effort in marking interesting SNVs.
makeScatterPlots = function(variants, samplePairs, timePoints, plotDirectory, genome='hg19', cpus=1, forceRedo=F) {
  scatterDirectory = paste0(plotDirectory, '/scatters')
  if ( !file.exists(scatterDirectory) ) dir.create(scatterDirectory)
  for ( pair in samplePairs ) {
    dir1 = paste0(scatterDirectory, '/', pair[1])
    if ( !file.exists(dir1) ) dir.create(dir1)
    dir2 = paste0(dir1, '/', pair[2])
    if ( !file.exists(dir2) | forceRedo ) {
      if ( !file.exists(dir2) ) dir.create(dir2)
      boring = variants$variants[[pair[1]]]$var == 0 & variants$variants[[pair[2]]]$var == 0
      q1 = variants$variants[[pair[1]]][!boring,]
      q2 = variants$variants[[pair[2]]][!boring,]
      ps=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, verbose=F, doPlot=F)
      psuf=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, plotFlagged=F, verbose=F, doPlot=F)
      
      outfile = paste0(dir2, '/all.png')
      catLog('Plotting to', outfile, '\n')
      png(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1)
      dev.off()
      
      outfile = paste0(dir2, '/allNamed.png')
      png(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1,
                     print=T, printRedCut=0.25)
      dev.off()

      outfile = paste0(dir2, '/allFlagged.png')
      png(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=ps,
                     main=paste0('all variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1)
      dev.off()
      
      catLog('Now by chromsome. Preparing..')
      chrs = xToChr(q1$x, genome=genome)
      catLog('done!\n  Plotting chr:')
      for ( chr in names(chrLengths(genome)) ) {
        outfile = paste0(dir2, '/chr', chr, '.png')
        png(outfile, width = 10, height=10, res=144, units='in')
        catLog(chr, '..', sep='')
        use = chrs == chr
        qualityScatter(q1[use,], q2[use,], variants$SNPs, ps=ps[use],
                       main=paste0('all variants: ', pair[1], ' vs ', pair[2], ', chr', chr),
                       xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                       ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1, print=T,
                       printRedCut=0.25, plotPosition=T, verbose=F)
        dev.off()
      }
      catLog('done!\n')
    }
  }
}


#helper function that generates fancy scatter plots of the frequencies of two samples.
#Requires two variants objects, and a SNPs object.
qualityScatter = function(q1, q2, SNPs, ps = NA, covScale=100, maxCex=1.5, minCov=10, main='', xlab='variant frequency: sample1', ylab='variant frequency: sample2', plotFlagged=T, cpus=1, verbose=T, print = F, printRedCut = 0.99, printOnlyNovel=F, plotPosition=F, genome='hg19', xlim=c(0,1), ylim=c(0,1), outputHighlighted=F, frame.plot=F, legend=T, redCut=0.75, forceCol=NA, add=F, GoI=c(), printCex=1, doPlot=T, minSomaticP=0, ...) {
  use = q1$var > 0 | q2$var > 0
  q1 = q1[use,]
  q2 = q2[use,]

  if ( minSomaticP > 0 & 'somaticP' %in% names(q1) & 'somaticP' %in% names(q2)  ) {
    use = pmax(q1$somaticP) >= minSomaticP
    q1 = q1[use,]
    q2 = q2[use,]    
  }
  
  freq1 = q1$var/q1$cov
  freq1[is.na(freq1)] = -0.02
  freq2 = q2$var/q2$cov
  freq2[is.na(freq2)] = -0.02

  q1$cov = q1$var + q1$ref
  q2$cov = q2$var + q2$ref

  temp = options()$scipen
  options(scipen = 100)
  xChar = as.character(q1$x)
  SNPs = SNPs[!duplicated(SNPs$x),]
  SNPxChar = as.character(SNPs$x)
  rownames(SNPs) = SNPxChar
  options(scipen = temp)
  keep = SNPs$x %in% q1$x
  if ( sum(keep) == 0 ) return()
  SNPs = SNPs[keep,]
  SNPs = SNPs[xChar,]
  db = q1$db
  flag1 = q1$flag
  flag2 = q2$flag

  #remove a few flags, such as single variant read in only one sample
  singleIn1 = grepl('Svr', flag1) & !grepl('Svr', flag2) & freq2 > 0
  flag1[singleIn1] = gsub('Svr', '', flag1[singleIn1])
  singleIn2 = grepl('Svr', flag2) & !grepl('Svr', flag1) & freq1 > 0
  flag2[singleIn2] = gsub('Svr', '', flag2[singleIn2])
  flag1 = gsub('Srr', '', flag1)
  flag2 = gsub('Srr', '', flag2)
  minorIn1 = grepl('Mv', flag1) & !grepl('Mv', flag2) & freq2 > 0
  flag1[minorIn1] = gsub('Mv', '', flag1[minorIn1])
  minorIn2 = grepl('Mv', flag2) & !grepl('Mv', flag1) & freq1 > 0
  flag2[minorIn2] = gsub('Mv', '', flag2[minorIn2])
  
  clean = flag1 == '' & flag2 == ''
  
  if ( verbose ) catLog(sum(clean), 'out of', nrow(q1), 'non-zero variants are not flagged!\n')

  if ( !plotFlagged ) {
    q1 = q1[clean,]
    q2 = q2[clean,]
    freq1 = q1$var/q1$cov
    freq1[is.na(freq1)] = -0.02
    freq2 = q2$var/q2$cov
    freq2[is.na(freq2)] = -0.02
    x = as.character(q1$x)
    SNPs = SNPs[clean,]
    db = q1$db
    clean = rep(T, nrow(q1))
  }

  if ( is.na(ps[1]) ) {
    doP = which(q1$cov >= max(1, minCov) & q2$cov >= max(1, minCov))
    if ( cpus == 1 )
      psCov = sapply(doP, function(i) fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value)
    else {
      psCov = unlist(mclapply(doP, function(i)
        fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value, mc.cores=cpus))
    }
    ps = rep(1, nrow(q1))
    ps[doP] = psCov
    names(ps) = rownames(q1)
  }
  dof = max(20, sum(clean & freq1+freq2 > 0 & freq1+freq2 < 2))
  if ( verbose ) catLog('MHC done with effective dof', dof, '\n')
  red = pmin(1, pmax(0, (-log10(ps)/log10(dof) - redCut)))
  if ( any(is.na(red)) ) {
    red[is.na(red)] = 0
    warning('Got NA red colour in scatter.')
  }
  if ( any(is.na(db)) ) {
    db[is.na(db)] = 0
    warning('Got NA sb in scatter.')
  }

  use = q1$cov >= minCov & q2$cov >= minCov
  if ( any(is.na(use)) ) {
    warning(paste0(sum(is.na(use)), ' NA entries in use'))
    use = use & !is.na(use)
  }
  if ( sum(use) == 0 ) invisible(ps)
  cleanOrder = which(clean&use)[order((red-0.1*db)[clean&use])]
  col = rep('black', length(red))
  col[!clean&use] = rgb(0.6 + red[!clean&use]*0.4, 0.6+red[!clean&use]*0.2, 0.6)
  col[clean&use] = rgb(red,0,ifelse(db, 0,(1-red)))[clean&use]
  if ( !is.na(forceCol) ) col = rep(forceCol[1], length(col))
  cex = pmin(maxCex, sqrt(sqrt(q1$cov*q2$cov)/covScale))
  if ( 'germline' %in% names(q1) ) pch = ifelse(clean, ifelse(db, 4, ifelse(q1$germline & !is.na(q1$germline), 3, 19)), ifelse(db, 4, 1))
  else pch = ifelse(clean, ifelse(db, 4, 19), ifelse(db, 4, 1))
  if ( !add & doPlot ) {
    plot(1, type='n', xlim=xlim, ylim=ylim, xlab = xlab, ylab = ylab, main=main, frame.plot=frame.plot, ...)
    segments(c(0,1,0,0,0, 0.5, 0), c(0,0, 0,1,0, 0, 0.5), c(0, 1, 1, 1, 1, 0.5, 1), c(1, 1, 0, 1, 1, 1, 0.5), col = rgb(0.8, 0.8, 0.8), lwd=0.3)
  }
  if ( !plotPosition & doPlot ) {
    if ( legend & !add ) {
      if ( plotFlagged )
        legend('bottomright', c('clean non-db', 'flagged non-db', 'clean germline-like non-db', 'clean db', 'flagged db', 'significantly different', 'high coverage', 'low coverage', 'protein altering', 'COSMIC Census Gene'), pch=c(19, 1, 3, 4, 4, 19, 19, 19, 1, 1), col=c('blue', 'grey', 'blue', 'black', 'grey', 'red', 'black', 'black', 'orange', 'green'), pt.cex=c(1,1,1,1,1,1,1,0.3, 1.5, 2), pt.lwd=c(1,1,1,1,1,1,1,1, 2, 4), bg='white')
      else
        legend('bottomright', c('not in dbSNP', 'in dbSNP', 'germline-like non-db', 'significantly different', 'low coverage', 'protein altering', 'COSMIC Census Gene'), pch=c(19, 4, 3, 19, 19, 1, 1), col=c('blue', 'black', 'blue', 'red', 'black', 'orange', 'green'), pt.cex=c(1,1,1,1,0.3, 1.5, 2), pt.lwd=c(1,1,1,1,1, 2, 4), bg='white')
    }
    points(freq1[!clean&use], freq2[!clean&use], cex=cex[!clean&use],
           lwd=pmin(maxCex, sqrt(sqrt(q1$cov*q2$cov)[!clean&use]/covScale)), pch=pch[!clean&use], col=col[!clean&use])
    points(freq1[cleanOrder], freq2[cleanOrder], cex=cex[cleanOrder],
           lwd=pmin(maxCex, sqrt(sqrt(q1$cov*q2$cov)[cleanOrder]/covScale)), pch=pch[cleanOrder], col=col[cleanOrder])
  }

  if ( 'severity' %in% names(q1) & 'severity' %in% names(q2) & doPlot ) {
    severity = pmin(q1$severity, q2$severity)
    severe = clean & use & severity <= 11 & !db
    points(freq1[severe], freq2[severe], cex=cex[severe]+(12-severity[severe])/10,
           lwd=(12-severity[severe])/2, pch=1, col='orange')
    if ( 'isCosmicCensus' %in% names(q1) & 'isCosmicCensus' %in% names(q2) ) {
      isCosmic = (q1$isCosmicCensus | q2$isCosmicCensus) & severe
      points(freq1[isCosmic], freq2[isCosmic], cex=cex[isCosmic]+(12-severity[isCosmic])/10+1,
             lwd=3, pch=1, col='green')
    }
  }
  
  if ( plotPosition & doPlot ) {
    segCex = pmin(1,(cex^2/7*abs((freq1+freq2)*(2-freq1-freq2))))
    segCex[freq1 < 0 | freq2 < 0] = 0
    pos = SNPs$start/chrLengths(genome)[as.character(SNPs$chr)]
    col = rgb(pmin(1,pmax(0, 1-abs(3*pos-0.5))), pmin(1,pmax(0, 1-abs(3*pos-1.5))), pmin(1,pmax(0, 1-abs(3*pos-2.5))), segCex)
    col[!clean] = rgb(0.5, 0.5, 0.5, segCex[!clean])
    lty = ifelse(use, ifelse(clean, 1, 2), 0)
    segments(freq1[!clean&use], freq2[!clean&use], pos[!clean&use], rep(1.02,sum(!clean&use)), lwd=segCex[!clean&use],
             col=col[!clean&use], lty=lty[!clean&use])
    segments(freq1[clean&use], freq2[clean&use], pos[clean&use], rep(1.02, sum(clean&use)), lwd=segCex[clean&use],
             col=col[clean&use], lty=lty[clean&use])
    col = rgb(pmin(1,pmax(0, 1-abs(3*pos-0.5))), pmin(1,pmax(0, 1-abs(3*pos-1.5))), pmin(1,pmax(0, 1-abs(3*pos-2.5))))
    col[!clean&use] = rgb(0.5, 0.5, 0.5)
    points(freq1[!clean&use], freq2[!clean&use], cex= cex[!clean&use], pch=pch[!clean&use], col = col[!clean&use])
    points(freq1[clean&use], freq2[clean&use], cex= cex[clean&use], pch=pch[clean&use], col = col[clean&use])
    points(pos[!clean&use], rep(1.02, length(pos[!clean&use])), pch=16, cex=0.5, col=col[!clean&use])
    points(pos[clean&use], rep(1.02, length(pos[clean&use])), pch=16, cex=0.5, col=col[clean&use])
    text(1.02, 1.02, 'SNPs')
  }

  if ( print & doPlot ) {
    toPrint = red > printRedCut
    if ( 'severity' %in% names(q1) & 'severity' %in% names(q2) ) {
      severity = pmin(q1$severity, q2$severity)
      severe = clean & use & severity <= 11 & !db
      if ( 'isCosmicCensus' %in% names(q1) & 'isCosmicCensus' %in% names(q2) ) {
        isCosmic = (q1$isCosmicCensus | q2$isCosmicCensus) & severe
        toPrint = toPrint | isCosmic
        col[isCosmic] = 'darkgreen'
        cex[isCosmic] = pmax(1.5, 1.5*cex[isCosmic])
      }
    }
    if ( length(GoI) > 0 )
      toPrint = toPrint | SNPs$inGene %in% GoI
    if ( printOnlyNovel ) toPrint = toPrint & freq1 < 0.1 & freq2 > 0.2
    if ( sum(toPrint) > 0 ) {
      printNames = gsub('.+:', '', SNPs[toPrint,]$inGene)
      printNames[is.na(printNames)] = 'i'
      if ( length(toPrint) > 0 )
        text(freq1[toPrint], freq2[toPrint] + 0.015*pmax(0.6, cex[toPrint]), printNames, col = col[toPrint], cex = pmax(0.6, cex[toPrint])*printCex)
      if ( 'isCosmicCensus' %in% names(q1) & 'isCosmicCensus' %in% names(q2) ) {
        isCosmic = (q1$isCosmicCensus | q2$isCosmicCensus) & severe
        if ( any(isCosmic) )
          text(freq1[isCosmic], freq2[isCosmic] + 0.015*pmax(0.6, cex[isCosmic]), printNames[isCosmic[which(toPrint)]],
               col = col[isCosmic], cex = pmax(0.6, cex[isCosmic])*printCex)
      }
      if ( verbose ) catLog('Highlighting', length(unique(printNames)), 'genes.\n')
      if ( outputHighlighted ) {
        out = data.frame('chr'=SNPs[toPrint,]$chr, 'pos'=SNPs[toPrint,]$start, 'x'=SNPs[toPrint,]$x, 'var'=q1$variant[toPrint], 'gene'=printNames,
          'f1'=freq1[toPrint], 'f2'=freq2[toPrint])
        print(out)
      }
    }
  }
  
  invisible(ps)
}
