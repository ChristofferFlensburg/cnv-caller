require(WriteXLS)
require(fields)


#prints heatmaps and line plots of the frequency of the SNVs in the samples of the same individual
#also outputs the results to excel files.
makeSNPprogressionPlots = function(variants, timeSeries, plotDirectory, cpus=1, v='', forceRedo=F) {
  msDirectory = paste0(plotDirectory, '/multiSample')
  if ( length(timeSeries) == 0 ) return()
  if ( !file.exists(msDirectory) ) dir.create(msDirectory)
  for ( i in 1:length(timeSeries) ) {
    ts = timeSeries[[i]]
    outfile = paste0(msDirectory, '/', names(timeSeries)[i], '.pdf')
    excelFileDB = paste0(msDirectory, '/', names(timeSeries)[i], '_DB.xls')
    excelFileNotDB = paste0(msDirectory, '/', names(timeSeries)[i], '_NotDB.xls')
    if ( !file.exists(outfile) | forceRedo ) {
      catLog('Plotting SNP progression to ', outfile, '..', sep='')
      pdf(outfile, width = 15, height=10, compress=T)
      qualityProgression(variants$variants[ts], variants$SNPs, nondb=F, excelFile=excelFileDB, main='dbSNPs only', v=v)
      qualityProgression(variants$variants[ts], variants$SNPs, db=F, excelFile=excelFileNotDB, main='non-dbSNPs only', v=v)
      dev.off()
      catLog('done.\n')
    }
  }
}


#helper function that does all the work for the frequency progression plots.
qualityProgression = function(qs, SNPs, db=T, nondb=T, excelFile='', main='', v='') {
  catLog('Finding common variants..')
  if ( !nondb ) qs = lapply(qs, function(q) q[q$db,])
  if ( !db ) qs = lapply(qs, function(q) q[!q$db,])
  common = Reduce(union, lapply(qs, rownames))
  if ( length(common) == 0 ) return()
  qs = lapply(qs, function(q) {
    new = common[!(common %in% rownames(q))]
    if ( length(new) > 0 ) {
      x = as.numeric(gsub('[AGNTCagtcn+-].*$', '', new))
      newQ = q[q$x %in% x,]
      newQ = newQ[!duplicated(newQ$x),]
      rownames(newQ) = newQ$x
      newQ = newQ[as.character(x),]
      newQ$var = 0
      newQ$flag = ''
      rownames(newQ) = new
      newQ$variant = gsub('^[0-9]+', '', new)
      q = rbind(q, newQ)
      q = q[common,]
    }
    return(q)
  })
  catLog(length(common), '. Cleaning..', sep='')
  var = do.call(cbind, lapply(qs, function(q) q$var))
  keep = rowSums(var) > 0
  qs = lapply(qs, function(q) q[keep,])
  common = common[keep]
  catLog(length(common), '..', sep='')
  fs = do.call(cbind, lapply(qs, function(q) q$var/q$cov))
  rownames(fs) = common
  cov = do.call(cbind, lapply(qs, function(q) q$cov))
  rownames(cov) = common
  var = do.call(cbind, lapply(qs, function(q) q$var))
  rownames(var) = common
  flag = do.call(cbind, lapply(qs, function(q) q$flag))
  rownames(flag) = common
  f = rowSums(var)/rowSums(cov)
  fs[is.na(fs)] = -0.02
  catLog('p-values..')
  ps = matrix(pBinom(as.integer(cov), as.integer(var), rep(f, ncol(cov))), ncol=ncol(cov))
  p = apply(ps, 1, fisherTest)[2,]

  SNPs = SNPs[SNPs$x %in% qs[[1]]$x,]
  gene = paste0(gsub('.+:', '', SNPs[as.character(qs[[1]]$x),]$inGene), ' (', SNPs[as.character(qs[[1]]$x),]$chr, ')')
  if ( db ) gene = SNPs[as.character(qs[[1]]$x),]$chr
  
  catLog('colours..')
  clean = rowSums(matrix(flag %in% c('', 'Svr'), ncol=ncol(flag))) == ncol(flag)
  dof = max(20, sum(clean & rowSums(fs) > 0 & rowSums(fs) < ncol(fs)))
  importance = pmin(1, noneg(-log10(p)/log10(dof) - 0.75))
  weight = pmin(1, pmax(importance, sqrt(rowMeans(cov/300))))
  doColour = importance > 0.5 & clean
  recurringGenes = gene[doColour][duplicated(gene[doColour])]
  isRecurringGene = gene %in% recurringGenes & doColour
  hue = rep(0, length(doColour))
  hue[isRecurringGene] = as.integer(as.factor(gene[isRecurringGene]))
  hue = hue/(max(hue)+1)
  col = D3colours(weight, importance, hue)

  if ( sum(doColour) > 1 ) {
    catLog('plotting heatmap..')
    rGcol = unique(D3colours(rep(1, sum(isRecurringGene)), rep(1, sum(isRecurringGene)), hue[isRecurringGene]))
    rGcolFlag = unique(D3colours(rep(0.4, sum(isRecurringGene)), rep(0.3, sum(isRecurringGene)), hue[isRecurringGene]))
    rG = unique(gene[isRecurringGene])
    names(rGcol) = names(rGcolFlag) = rG
    clusterOrder =
      heatmap(fs[doColour,,drop=F], cexCol=1, labRow=gene[doColour],
              RowSideColors = ifelse(gene[doColour] %in% rG, ifelse(clean[doColour],
                rGcol[gene[doColour]], rGcolFlag[gene[doColour]]), ifelse(clean[doColour], rgb(0.65, 0.65, 0.65),'grey')),
              col=sapply((0:100)/100, function(heat) D3colours(1, heat, heat)), margins=c(8,15), scale='none', main=main)
    colorbar.plot(par('usr')[1]*0.9+par('usr')[2]*0.1, par('usr')[3]*0.5+par('usr')[4]*0.5,
                  strip.width = 0.05, strip.length = 0.4, 0:1000,
                  col=sapply((0:1000)/1000, function(heat) D3colours(1, heat, heat)), margins=c(10,5), horizontal=F)
    segments(par('usr')[1]*0.93+par('usr')[2]*0.07, par('usr')[3]*0.5+par('usr')[4]*0.5,
             par('usr')[1]*0.925+par('usr')[2]*0.075, par('usr')[3]*0.5+par('usr')[4]*0.5, lwd=2, D3colour(c(1,0.5,0.5)))
    if ( length(rG) > 0 ) {
      legend('right', rG, col = rGcol, lwd=10, bg='white')
    }
    fs = fs[,clusterOrder[[2]]]
  }
  else catLog('skipping heatmap (no important variants)..')

  N = ncol(fs)

  catLog('frequency progression..')
  is = order(importance+isRecurringGene)
  is = is[!qs[[1]]$db[is] | importance[is] > 0]
  plot(0, type='n', xlim=c(1, length(qs)*1.2), ylim=c(0,1), main=main)
  segments(col(fs)[is, 1:(N-1)], fs[is, 1:(N-1)], col(fs)[is, 2:N], fs[is, 2:N], lwd=weight[is]+importance[is], col=col[is])
  text(1:N, 1.02, colnames(fs), cex=0.7)

  if ( sum(doColour) > 1 ) {
      if ( length(rG) > 0 ) legend('right', rG, col = rGcol, lwd=3, bg='white')
  }
  catLog('done!\n')

  if ( excelFile != '' ) {
    catLog('Output plotted variants to', excelFile, '...')
    colnames(fs) = names(qs)
    multiSampleData = data.frame(
      'gene'=gsub('.+:', '', SNPs[as.character(qs[[1]]$x[doColour]),]$inGene),
      'chr'=SNPs[as.character(qs[[1]]$x[doColour]),]$chr,
      'start'=qs[[1]]$x[doColour],
      'end'=qs[[1]]$x[doColour],
      'reference'=qs[[1]]$reference[doColour],
      'variant'=qs[[1]]$variant[doColour],
      'f'=fs[doColour,],
      'cov'=cov[doColour,],
      'var'=var[doColour,])
    WriteXLS('multiSampleData', excelFile)
    catLog('done!\n')
  }
  
}
