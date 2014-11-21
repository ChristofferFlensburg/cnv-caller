


#Takes variants and plots the frequencies of any two samples from the same individual against each other.
#Goes to some effort in marking interesting SNVs.
makeScatterPlots = function(variants, samplePairs, timePoints, plotDirectory, genome='hg19', cpus=1, v='', forceRedo=F) {
  scatterDirectory = paste0(plotDirectory, '/scatters')
  if ( !file.exists(scatterDirectory) ) dir.create(scatterDirectory)
  for ( pair in samplePairs ) {
    dir1 = paste0(scatterDirectory, '/', pair[1])
    if ( !file.exists(dir1) ) dir.create(dir1)
    dir2 = paste0(dir1, '/', pair[2])
    if ( !file.exists(dir2) | forceRedo ) {
      dir.create(dir2)
      boring = variants$variants[[pair[1]]]$var == 0 & variants$variants[[pair[2]]]$var == 0
      q1 = variants$variants[[pair[1]]][!boring,]
      q2 = variants$variants[[pair[2]]][!boring,]
      ps=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, verbose=F)
      psuf=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, plotFlagged=F, verbose=F)
      
      outfile = paste0(dir2, '/all.jpg')
      catLog('Plotting to', outfile, '\n')
      jpeg(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1)
      dev.off()
      
      outfile = paste0(dir2, '/allNamed.jpg')
      jpeg(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1,
                     print=T, printRedCut=0.25)
      dev.off()

      outfile = paste0(dir2, '/allFlagged.jpg')
      jpeg(outfile, width = 10, height=10, res=300, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=ps,
                     main=paste0('all variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1], ' (', timePoints[pair[1]], ')'),
                     ylab=paste0('variant frequency for ', pair[2], ' (', timePoints[pair[2]], ')'), cpus=1)
      dev.off()
      
      catLog('Now by chromsome. Preparing..')
      chrs = xToChr(q1$x, genome=genome)
      catLog('done!\n  Plotting chr:')
      for ( chr in names(chrLengths(genome)) ) {
        outfile = paste0(dir2, '/chr', chr, '.jpg')
        jpeg(outfile, width = 10, height=10, res=300, units='in')
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
qualityScatter = function(q1, q2, SNPs, ps = NA, covScale=100, maxCex=1.5, minCov=10, main='', xlab='variant frequency: sample1', ylab='variant frequency: sample2', plotFlagged=T, cpus=1, verbose=T, print = F, printRedCut = 0.99, printOnlyNovel=F, plotPosition=F, genome='hg19', xlim=c(0,1), ylim=c(0,1), outputHighlighted=F) {
  use = q1$var > 0 | q2$var > 0
  q1 = q1[use,]
  q2 = q2[use,]
  freq1 = q1$var/q1$cov
  freq1[is.na(freq1)] = -0.02
  freq2 = q2$var/q2$cov
  freq2[is.na(freq2)] = -0.02

  x = as.character(q1$x)
  keep = SNPs$x %in% x
  if ( sum(keep) == 0 ) return()
  SNPs = SNPs[keep,]
  SNPs = SNPs[x,]
  db = SNPs$db
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
  
  if ( verbose ) cat(sum(clean), 'out of', nrow(q1), 'non-zero variants are not flagged!\n')

  if ( !plotFlagged ) {
    q1 = q1[clean,]
    q2 = q2[clean,]
    freq1 = q1$var/q1$cov
    freq1[is.na(freq1)] = -0.02
    freq2 = q2$var/q2$cov
    freq2[is.na(freq2)] = -0.02
    x = as.character(q1$x)
    SNPs = SNPs[clean,]
    db = db[clean]
    clean = rep(T, nrow(q1))
  }

  if ( is.na(ps[1]) ) {
    if ( cpus == 1 )
      ps = sapply(1:nrow(q1), function(i) fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value)
    else {
      require(parallel)
      ps = unlist(mclapply(as.list(1:nrow(q1)), function(i)
        fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value, mc.cores=cpus))
    }
    names(ps) = rownames(q1)
  }
  dof = max(20, sum(clean & freq1+freq2 > 0 & freq1+freq2 < 2))
  if ( verbose ) cat('MHC done with effective dof', dof, '\n')
  red = pmin(1, pmax(0, (-log10(ps)/log10(dof) - 0.75)))

  covScale = 100
  maxSize = 1.5
  covCut = 10
  use = q1$cov >= covCut & q2$cov >= covCut
  if ( sum(use) == 0 ) invisible(ps)
  cleanOrder = which(clean&use)[order((red-0.1*db)[clean&use])]
  col = rep('black', length(red))
  col[!clean&use] = rgb(0.6 + red[!clean&use]*0.4, 0.6+red[!clean&use]*0.2, 0.6)
  col[clean&use] = rgb(red,0,(1-red)*pmin(1, q1$somaticP+q2$somaticP))[clean&use]
  cex = pmin(maxSize, sqrt(sqrt(q1$cov*q2$cov)/covScale))
  pch = ifelse(clean, ifelse(db, 4, 19), ifelse(db, 4, 1))
  plot(1, type='n', xlim=xlim, ylim=ylim, xlab = xlab, ylab = ylab, main=main)
  segments(c(0,1,0,0,0, 0.5, 0), c(0,0, 0,1,0, 0, 0.5), c(0, 1, 1, 1, 1, 0.5, 1), c(1, 1, 0, 1, 1, 1, 0.5), col = rgb(0.8, 0.8, 0.8), lwd=0.3)
  if ( !plotPosition ) {
    if ( plotFlagged )
      legend('bottomright', c('clean non-db', 'flagged non-db', 'clean db', 'flagged db', 'significantly different', 'low coverage'), pch=c(19, 1, 4, 4, 19, 19), col=c('blue', 'grey', 'black', 'grey', 'red', 'black'), pt.cex=c(1,1,1,1,1,0.3))
    else
      legend('bottomright', c('somatic', 'db', 'significantly different', 'low coverage'), pch=c(19, 4, 19, 19), col=c('blue', 'black', 'red', 'black'), pt.cex=c(1,1,1,0.3))
    points(freq1[!clean&use], freq2[!clean&use], cex=cex[!clean&use],
           lwd=pmin(maxSize, sqrt(sqrt(q1$cov*q2$cov)[!clean&use]/covScale)), pch=pch[!clean&use], col=col[!clean&use])
    points(freq1[cleanOrder], freq2[cleanOrder], cex=cex[cleanOrder],
           lwd=pmin(maxSize, sqrt(sqrt(q1$cov*q2$cov)[cleanOrder]/covScale)), pch=pch[cleanOrder], col=col[cleanOrder])
  }
  
  if ( plotPosition ) {
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

  if ( print ) {
    toPrint = red > printRedCut
    if ( printOnlyNovel ) toPrint = toPrint & freq1 < 0.1 & freq2 > 0.2
    if ( sum(toPrint) > 0 ) {
      printNames = gsub('.+:', '', SNPs[toPrint,]$inGene)
      printNames[is.na(printNames)] = 'i'
      if ( length(toPrint) > 0 )
        text(freq1[toPrint], freq2[toPrint] + 0.015*pmax(0.6, cex[toPrint]), printNames, col = col[toPrint], cex = pmax(0.6, cex[toPrint]))
      if ( verbose ) {
        cat('Highlighting', length(unique(printNames)), 'genes.\n')
        if (any(DLBCLgenes() %in% printNames) )  cat('DLBCL gene! ', intersect(DLBCLgenes(), printNames), '\n')
        if (any(AMLgenes() %in% printNames) )  cat('AML gene! ', intersect(AMLgenes(), printNames), '\n')
      }
      if ( outputHighlighted ) {
        out = data.frame('chr'=SNPs[toPrint,]$chr, 'pos'=SNPs[toPrint,]$start, 'x'=SNPs[toPrint,]$x, 'var'=q1$variant[toPrint], 'gene'=printNames,
          'f1'=freq1[toPrint], 'f2'=freq2[toPrint])
        print(out)
      }
    }
  }
  
  invisible(ps)
}
