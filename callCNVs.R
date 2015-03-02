


#This function takes the cancer and normal variants, and the coverage analysis as input
#It tries to find the germline het variants for each sample.
#It then clusters regions of the genome that have similar het frequencies and coverage.
#It then calls the CN of each clustered region, as well as clonality, uncertainty estimate and p-value
#for normal CN.
#returns a list of data frames with each region of the genome on a row.
callCNVs = function(variants, normalVariants, fitS, names, individuals, normals, Rdirectory, plotDirectory, genome='hg19', cpus=1, forceRedoCNV=F) {
  clustersSaveFile = paste0(Rdirectory, '/clusters.Rdata')
  if ( file.exists(clustersSaveFile) & !forceRedoCNV ) {
    catLog('Loading saved CNV results.\n')
    load(file=clustersSaveFile)
    return(clusters)
  }

  #some capture regions can name regions very far apart the same name, which will be treated as a
  #single region, messing up resolution. So remove any regions larger than 5Mbp.
  #If this causes problems, rename the capture regions to avoid equal names at long distance.
  fitS = subsetFit(fitS, rows = abs(fitS$x2 - fitS$x1) < 5e6)

  #identify cancer-normal pairs, check which cancers have normals from the same individual
  correspondingNormal = findCorrespondingNormal(names, individuals, normals)

  #run cnv on the samples, using cancer-normal where available
  clusters = lapply(names, function(name) {
    catLog('\nCalling CNVs for ', name, '.\n', sep='')
    if ( !is.na(correspondingNormal[name]) )
      return(callCancerNormalCNVs(cancerVariants=variants$variants[[name]],
                                  normalVariants=variants$variants[[correspondingNormal[name]]],
                                  moreNormalVariants=normalVariants$variants,
                                  fit = subsetFit(fitS, cols=paste0(name, '-normal')),
                                  plotDirectory, name,
                                  genome=genome, cpus=cpus))
    else
      return(callCancerNormalCNVs(cancerVariants=variants$variants[[name]],
                                  normalVariants=FALSE,
                                  moreNormalVariants=normalVariants$variants,
                                  fit = subsetFit(fitS, cols=paste0(name, '-normal')),
                                  plotDirectory, name,
                                  genome=genome, cpus=cpus))
  })

  names(clusters) = names
  save(clusters, file=clustersSaveFile)
  return(clusters)
}


#the high level function that controls the steps of the CNV calling for given sample and normal variant objects.
callCancerNormalCNVs = function(cancerVariants, normalVariants, moreNormalVariants, fit, plotDirectory, name, genome='hg19', cpus=1) {
  require(parallel)
  #estimate reference bias and variance from the selected normal hets.
  if ( class(normalVariants) == 'logical') {
    a=try(setVariantLoss(moreNormalVariants))
    if ( class(a) == 'try-error' ) {
      catLog('Error in setVariantLoss(moreNormalVariants)\n')
      stop('Error in setVariantLoss(moreNormalVariants)!')
    }
  }
  else {
    a=try(setVariantLoss(normalVariants))
    if ( class(a) == 'try-error' ) {
      catLog('Error in setVariantLoss(normalVariants)\n')
      stop('Error in setVariantLoss(normalVariants)!')
    }
  }

  #select good germline het variants from normals:
  if ( class(normalVariants) == 'logical') {
    use = try(selectGermlineHetsFromCancer(cancerVariants, moreNormalVariants, fit$sex, cpus=cpus))
    if ( class(use) == 'try-error' ) {
      catLog('Error in selectGermlineHetsFromCancer(cancerVariants, moreNormalVariants, cpus=cpus)\n')
      stop('Error in selectGermlineHetsFromCancer(cancerVariants, moreNormalVariants, cpus=cpus)!')
    }
  }
  else {
    use = try(selectGermlineHets(normalVariants, moreNormalVariants, fit$sex, cpus=cpus))
    if ( class(use) == 'try-error' ) {
      catLog('Error in selectGermlineHets(normalVariants, moreNormalVariants, cpus=cpus)\n')
      stop('Error in selectGermlineHets(normalVariants, moreNormalVariants, cpus=cpus)!')
    }
  }
  is = rownames(cancerVariants) %in% use
  cancerVariants = cancerVariants[is,]
  
  #summarise by capture region
  catLog('Summarising capture regions..')
  mC = maxCov()
  effectiveCov = round((cancerVariants$cov+2*cancerVariants$cov^2/(mC/2))/(1+cancerVariants$cov/(mC/2))^2)
  effectiveVar = round(cancerVariants$var/cancerVariants$cov*effectiveCov)
  effectiveFreqs = data.frame(var=mirrorDown(effectiveVar, cov=effectiveCov),
    cov=effectiveCov, x=cancerVariants$x)
  cancerCR = try(unifyCaptureRegions(effectiveFreqs, fit, cpus=cpus))
  if ( class(cancerCR) == 'try-error' ) {
    catLog('Error in unifyCaptureRegions(cancerVariants, fit, cpus=cpus)\n')
    stop('Error in unifyCaptureRegions(cancerVariants, fit, cpus=cpus)!')
  }
  catLog('done!\n')

  #correct width if variance is underestimated
  cancerCR = boostCRwidth(cancerCR)

  #run clustering algorithm
  cancerCluster = mergeChromosomes(cancerCR, effectiveFreqs, genome=genome, cpus=cpus)

  #run post processing, such as correcting normalisation from AB regions, calling CNVs and clonalities.
  catLog('Postprocessing...')
  post = postProcess(cancerCluster, cancerCR, effectiveFreqs, plotDirectory, name, genome, cpus=cpus)
  cancerCluster = post$clusters
  catLog('found total of', sum(cancerCluster$call != 'AB' & !grepl('\\?', cancerCluster$call)), 'CNVs..')
  cancerCR$M = cancerCR$M - post$meanM
  #plotCR(cancerCluster)

  #return clustered regions with calls, as well as the raw capture region data.
  catLog('done!\n')
  effectiveFreqs = data.frame(var=effectiveVar, cov=effectiveCov, x=cancerVariants$x)
  return(list(clusters=cancerCluster, CR=cancerCR,
              eFreqs=effectiveFreqs))
}


#helper function that selects germline het SNPs in the presence of a normal sample from the same individual.
selectGermlineHets = function(normalVariants, moreNormalVariants, sex, minCoverage = 10, cpus=1) {
  #only bother with variants that have enough coverage so that we can actually see a change in frequency
  catLog('Taking variants with minimum coverage of', minCoverage, '...')
  decentCoverage = normalVariants$cov >= minCoverage
  use = rownames(normalVariants)[decentCoverage]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(decentCoverage), 'variants.\n')
  if ( length(use) == 0 ) return(use)

  #only use dbSNPs
  catLog('Taking variants in dbSNP positions...')
  isDB = normalVariants$db
  use = use[isDB]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(isDB), 'variants.\n')
  if ( length(use) == 0 ) return(use)

  #dont use insertions or deletions, too error-prone.
  catLog('Removing indels..')
  indel = grepl('[+-]', normalVariants$variant)
  use = use[!indel]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(!indel), 'variants.\n')

  #pick het in matching normal sample
  catLog('Taking variants that are het in the normal sample..')
  normalF = normalVariants$var/normalVariants$cov
  normalHet = pBinom(normalVariants$cov, normalVariants$var, refBias(0.5)) > 0.1 & abs(normalF-0.5) < 0.25 & normalVariants$flag == ''
  use = use[normalHet]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(normalHet), 'variants.\n')
  if ( length(use) == 0 ) return(use)
  
  #again to filter out noisy variants, filter on the RIB statistic in the cancer sample
  catLog('Taking variants that have low expected rate of incorrect basecalls in the normal sample..')
  highNormalRIB = normalVariants$RIB > 0.03
  use = use[!highNormalRIB]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(!highNormalRIB), 'variants.\n')

  #be extra picky with the mapping and base quality
  catLog('Taking variants that have good mapping and base quality in the normal sample..')
  highQ = normalVariants$pbq > 0.001 & normalVariants$pmq > 0.001
  use = use[highQ]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(highQ), 'variants.\n')

  #remove variants from male X and Y, and female Y-chromsomes
  catLog('Removing variants in male X and Y, and female Y...')
  notHet = xToChr(normalVariants$x) %in% ifelse(sex=='male', c('X','Y'), 'Y')
  use = rownames(normalVariants)[!notHet]
  normalVariants = normalVariants[use,]
  catLog('done! Discarded', sum(notHet), 'variants. This is an indication of the number of false SNPs in other chromsomes.\n')
  if ( length(use) == 0 ) return(use)

  return(use)
}

#helper function that selects germline het SNPs in the absence of a normal sample from the same individual.
selectGermlineHetsFromCancer = function(cancerVariants, moreNormalVariants, sex, minCoverage = 10, cpus=1) {
  #only bother with variants that have enough coverage so that we can actually see a change in frequency
  catLog('Taking variants with minimum coverage of', minCoverage, '...')
  decentCoverage = cancerVariants$cov >= minCoverage
  use = rownames(cancerVariants)[decentCoverage]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(decentCoverage), 'variants.\n')
  
  #only use dbSNPs
  catLog('Taking variants in dbSNP positions...')
  isDB = cancerVariants$db
  use = use[isDB]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(isDB), 'variants.\n')
  
  #dont use insertions or deletions, too error-prone.
  catLog('Removing indels..')
  indel = grepl('[+-]', cancerVariants$variant)
  use = use[!indel]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(!indel), 'variants.\n')

  #take variants between 5% and 95%
  catLog('Taking unflagged cancer variants that have frequency above 5% and below 95%..')
  cancerF = cancerVariants$var/cancerVariants$cov
  normalHet = abs(cancerF-0.5) < 0.45 & cancerVariants$flag == ''
  use = use[normalHet]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(normalHet), 'variants.\n')

  #require support from at least 2 reads
  catLog('Taking cancer variants that have support from at least 2 reads..')
  cancerF = cancerVariants$var/cancerVariants$cov
  supported = cancerVariants$var > 1 & cancerVariants$var < cancerVariants$cov-1
  use = use[supported]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(supported), 'variants.\n')

  
  #again to filter out noisy variants, filter on the RIB statistic in the cancer sample
  catLog('Taking variants that have low expected rate of incorrect basecalls in the sample..')
  highCancerRIB = cancerVariants$RIB > 0.03
  use = use[!highCancerRIB]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(!highCancerRIB), 'variants.\n')

  #be extra picky with the mapping and base quality
  catLog('Taking variants that have good mapping and base quality..')
  highQ = cancerVariants$pbq > 0.001 & cancerVariants$pmq > 0.001
  use = use[highQ]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(highQ), 'variants.\n')

  #remove variants from male X and Y, and female Y-chromsomes
  catLog('Removing variants in male X and Y, and female Y...')
  notHet = xToChr(cancerVariants$x) %in% ifelse(sex=='male', c('X','Y'), 'Y')
  use = rownames(cancerVariants)[!notHet]
  cancerVariants = cancerVariants[use,]
  catLog('done! Discarded', sum(notHet), 'variants. This is an indication of the number of false SNPs in other chromsomes.\n')
  if ( length(use) == 0 ) return(use)

  return(use)
}

unifyCaptureRegions = function(eFreqs, fit, cpus=1) {
  #group SNPs by capture region
  x = eFreqs$x
  snps = lapply(1:nrow(fit), function(row) which(x < fit$x2[row] & x > fit$x1[row]))

  uniFreq = do.call(rbind, lapply(snps, function(is) c(sum(eFreqs$var[is]), sum(eFreqs$cov[is]))))
  colnames(uniFreq) = c('var', 'cov')
  
  pHet = rep(0.5, length(snps))
  pHet[uniFreq[,'cov'] > 0] = unlist(mclapply(snps[uniFreq[,'cov'] > 0], function(is) fisherTest(pBinom(eFreqs$cov[is], eFreqs$var[is], refBias(0.5)))[2], mc.cores=cpus))

  pAlt = rep(0.5, length(snps))
  pAlt[uniFreq[,'cov'] > 0] = unlist(lapply(snps[uniFreq[,'cov'] > 0], function(is)
        fisherTest(pTwoBinom(eFreqs$cov[is], eFreqs$var[is], sum(eFreqs$var[is])/sum(eFreqs$cov[is])))[2]))

  odsHet = rep(0.5, length(snps))
  odsHet[uniFreq[,'cov'] > 0] = unlist(lapply(snps[uniFreq[,'cov'] > 0], function(is) {
    fAlt = refBias(alternativeFrequency(eFreqs[is,])['f'])
    dpHet = dbinom(eFreqs$var[is], eFreqs$cov[is], refBias(0.5))
    dpAlt = twoBinom(eFreqs$var[is], eFreqs$cov[is], fAlt, refBiasMirror(fAlt))
    odsHet = dpHet/(dpHet+dpAlt)
    return(fisherTest(odsHet)[2])
  }))

  cR = data.frame(x1=fit$x1, x2=fit$x2, M=fit$coefficients[,1], width=fit$coefficients[,1]/fit$t[,1], df=fit$df.total,
    uniFreq[,1:2], pHet=pHet, pAlt=pAlt, odsHet=odsHet)
  return(cR)
}

alternativeFrequency = function(efs, plot=F) {
  f = (0:200)/200
  rbf = refBias(f)
  fps = sapply(1:nrow(efs), function(row) {
    p=pBinom(rep(efs$cov[row], 201), rep(efs$var[row], 201), rbf)
    return(p/sum(p))
    })
  MAFchi2 = apply(fps + fps[201:1,], 1, fisherTest)[1,]
  minChi2 = which(MAFchi2 == min(MAFchi2))[1]
  within1pc = which(MAFchi2 < min(MAFchi2)*1.01 & f <= 0.5)
  error = max(1/200, max(abs(f[within1pc]-f[minChi2])), 1/sqrt(sum(efs$cov)))
  if ( plot ) plot(f[1:101], MAFchi2[1:101], ylim=c(min(MAFchi2), min(MAFchi2)*2), pch=16, col=ifelse(1:101 %in% within1pc, 'red', 'black'))
  return(c('f'=f[minChi2], 'ferr'=error))
}


#merges regions in each chromosome.
mergeChromosomes = function(cR, eFreqs, genome='hg19', cpus=1, ...) {
  chrs = xToChr(cR$x1, genome=genome)
  
  catLog('Merging capture regions with same coverage and MAF: ')
  clusters = mclapply(unique(chrs), function(chr) {
    catLog(chr, '..', sep='')
    ret = mergeRegions(cR[chrs == chr,], ...)
    singleSNP = calledFromSingleSNP(ret, eFreqs)
    if ( any(singleSNP) ) ret = mergeRegions(ret, forceMerge=singleSNP, ...)
    ret$altStatErr = ret$nullStatErr = ret$altStat = ret$nullStat = ret$stat = ret$f = rep(NA, nrow(ret))
    ret$postHet = rep(1, nrow(ret))
    if ( any(ret$cov > 0) ) {
      ret$f[ret$cov > 0] = refUnbias(ret$var[ret$cov>0]/ret$cov[ret$cov>0])
      ret[ret$cov>0,] = redoHetCalculations(ret[ret$cov>0,], eFreqs, cpus=1)
      ret = mergeRegions(ret, ...)
    }
    return(ret)
  }, mc.cores=cpus)
  
  clusters = do.call(rbind, clusters)
  catLog('done!\n')
  return(clusters)
}
mergeRegions = function(cR, minScore = 0.05, plot=F, debug=F, forceMerge=NA) {
  merged = c()
  scores = c()
  if ( debug ) Nloop = 0
  cR = cR[order(cR$x1+cR$x2),]
  pairScore = sameCNV(cR)
  if ( !is.na(forceMerge[1]) & length(forceMerge) == length(pairScore) ) pairScore[forceMerge] = 1
  while(dim(cR)[1] > 1) {
    #find the pair of regions with the largest pairing probability
    best = which(pairScore == max(pairScore))[1]
    if ( debug && Nloop < 1 ) {
      print (cR)
      catLog('DEBUG: Merging pair at ', cR$x2[best], ' with score ', pairScore[best], ', ', nrow(cR), ' regions.\n', sep='')
      plotCR(cR)
      points(cR$x1[-1], pairScore, pch=4)
      segments((cR$x1[best+1]+cR$x2[best])/2, 0, (cR$x1[best+1]+cR$x2[best])/2, 1)
      a = scan()
      if ( length(a) > 0 ) Nloop = a
      else Nloop = 1
      catLog('DEBUG: looping another ', Nloop, ' times.\n', sep='')
    }
    if ( debug ) Nloop = Nloop - 1

    if ( plot & length(scores) == 0 ) {
      plotCR(cR)
    }

    #break loop if no clusters that are sufficiently similar
    if ( pairScore[best] < minScore ) break
    #otherwise save the merged regions
    merged = c(merged, best)
    scores = c(scores, pairScore[best])
    
    #merge the regions, using the provided weight
    cR$x1[best] = min(cR$x1[best],cR$x1[best+1])
    cR$x2[best] = max(cR$x2[best],cR$x2[best+1])
    cR$var[best] = cR$var[best] + cR$var[best+1]
    cR$cov[best] = cR$cov[best] + cR$cov[best+1]
    cR$M[best] = (cR$M[best]/cR$width[best]^2 + cR$M[best+1]/cR$width[best+1]^2)/(1/cR$width[best]^2+1/cR$width[best+1]^2)
    cR$width[best] = 1/sqrt(1/cR$width[best]^2 + 1/cR$width[best+1]^2)
    cR$pHet[best] = if ( cR$cov[best] + cR$cov[best+1] > 0 ) stoufferTest(c(cR$pHet[best], cR$pHet[best+1]), c(cR$cov[best], cR$cov[best+1]))[2] else 0.5
    cR = cR[-(best+1),]
    pairScore = pairScore[-best]
    if ( best > 1 ) pairScore[best-1] = sameCNV(cR[(best-1):best,])
    if ( best <= length(pairScore) )  pairScore[best] = sameCNV(cR[best:(best+1),])
  }
  
  if ( plot ) {
    plotCR(cR)
    layout(1)
  }
  return(cR)
}

calledFromSingleSNP = function(cR, eFreqs) {
  first = 1:(nrow(cR)-1)
  second = 2:nrow(cR)
  
  x = eFreqs$x
  snps = lapply(first, function(row) which(x < cR$x2[row] & x > cR$x1[row]))
  singleSNP = sapply(snps, length) == 1 | sapply(snps, function(is) max(eFreqs$x[is]) - min(eFreqs$x[is]) < fragmentLength())

  width = cR$width + 0*systematicVariance()
  meanM = (cR$M[first]/width[first]^2 + cR$M[second]/width[second]^2)/(1/width[first]^2 + 1/width[second]^2)
  MP1 = pt(-abs(cR$M[first] - meanM)/width[first], df = cR$df[first])
  MP2 = pt(-abs(cR$M[second] - meanM)/width[second], df = cR$df[second])
  pM = sapply(first, function(i) min(p.adjust(c(MP1[i], MP2[i]), method='fdr')))

  return(pM > 0.05 & singleSNP)
}

fragmentLength = function() {return(300)}

#takes clustered regions and post processes them by calling copy numbers and clonality.
postProcess = function(clusters, cRs, eFreqs, plotDirectory, name, genome='hg19', cpus=1) {
  catLog('corrected frequencies..')
  cf = correctedFrequency(clusters, eFreqs, cpus=cpus)
  clusters$f = cf[,'f']
  clusters$ferr = cf[,'ferr']
  catLog('frequency probability..')
  diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
  if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
  diagnosticMAFDirectory = paste0(diagnosticPlotsDirectory, '/MAFstat/')
  if ( !file.exists(diagnosticMAFDirectory) ) dir.create(diagnosticMAFDirectory)
  plotFile = paste0(diagnosticMAFDirectory, name, '.pdf')
  catLog('plotting MAF stats to', plotFile, '..')
  pdf(plotFile, width=15, height=7)
  clusters = redoHetCalculations(clusters, eFreqs, cpus=cpus, plot=T)
  dev.off()
  catLog('call CNVs..')
  clusters = addCall(clusters, eFreqs)
  catLog('renormalise from AB calls..')
  renorm =  normaliseCoverageToHets(clusters)
  clusters = renorm$clusters
  catLog('redo CNV calls..')
  clusters = addCall(clusters, eFreqs)
  clusters = fixFalseSNPcall(clusters, eFreqs)
  sCAN = sameCallAsNeighbour(clusters)
  if ( any(sCAN) ) {
    catLog('found', sum(sCAN),'neighbouring regions with same call: merge and redo postprocessing!\n')
    clusters = forceMerge(clusters, sCAN)
    ret = postProcess(clusters, cRs, eFreqs, plotDirectory, name, genome=genome, cpus=cpus)
    ret$meanM = ret$meanM + renorm$meanM
    return(ret)
  }
  catLog('find subclones..')
  clusters = findSubclones(clusters)
  rownames(clusters) = make.names(paste0('chr', xToChr(clusters$x1, genome=genome)), unique=T)
  clusters = clusters[order(clusters$x1),]
  return(list(clusters=clusters, meanM=renorm$meanM))
}

fixFalseSNPcall = function(clusters, eFreqs) {
  isAA = gsub('\\?', '', clusters$call) == 'AA'
  smallRegion = clusters$x2 - clusters$x1 < 1e7
  x = eFreqs$x
  snps = lapply(1:nrow(clusters), function(row) which(x < clusters$x2[row] & x > clusters$x1[row]))
  basedOnFewSNPs = sapply(snps, length) < 5
  basedOnSmallRegion = sapply(snps, function(is) length(is) > 0 && max(eFreqs$x[is]) - min(eFreqs$x[is]) < 1e5)
  inconsistentEvidence = clusters$sigma > 3
  falseCall = isAA & smallRegion & (basedOnFewSNPs | basedOnSmallRegion | inconsistentEvidence)
  if ( any(falseCall) ) {
    catLog('found', sum(falseCall), 'AA calls from false SNPs..')
    isab = sapply(which(falseCall), function(i) unlist(isAB(clusters[i,], eFreqs[snps[[i]],])))
    clusters[falseCall,]$sigma = isab[2,]
    clusters[falseCall,]$pCall = clusters[falseCall,]$pHet
    clusters[falseCall,]$clonality = 1
    clusters[falseCall,]$clonalityError = 0
    if ( 'subclonality' %in% names(clusters) ) clusters[falseCall,]$subclonality = 1
    if ( 'subclonalityError' %in% names(clusters) ) clusters[falseCall,]$subclonalityError = 0
    clusters[falseCall,]$call = ifelse(isab[2,] > 10, 'AB??', ifelse(isab[2,] > 5, 'AB?', 'AB'))
    }
  return(clusters)
}

sameCallAsNeighbour = function(clusters) {
  chr = xToChr(clusters$x1)
  first = 1:(nrow(clusters)-1)
  second = 2:nrow(clusters)
  sameCall = chr[first] == chr[second] &
  gsub('\\?','',clusters$call[first]) == gsub('\\?','',clusters$call[second]) &
  abs(clusters$clonality[first] - clusters$clonality[second]) <=
    sqrt(clusters$clonalityError[first]^2 + clusters$clonalityError[second]^2)*2
  return(sameCall)
}

forceMerge = function(clusters, toMerge) {
  chrs = xToChr(clusters$x1)
  redoChrs = unique(chrs[which(toMerge)])
  newClusters = list()
  for ( chr in redoChrs ) {
    is = which(chrs == chr)
    newClusters[[chr]] = mergeRegions(clusters[is,], forceMerge=toMerge[is][-length(is)])
  }
  clusters = clusters[!(chrs %in% redoChrs),]
  for ( chr in redoChrs ) {
    clusters = rbind(clusters, newClusters[[chr]])
  }
  clusters = clusters[order(clusters$x1),]
  return(clusters)
}

correctedFrequency = function(clusters, eFreqs, maxCov = maxCov(), cpus=1) {
  ret = do.call(rbind, mclapply(1:nrow(clusters), function(row) {
    efs = eFreqs[eFreqs$x > clusters$x1[row] & eFreqs$x < clusters$x2[row],]
    if ( nrow(efs) == 0 ) return(c('f'=NA, 'ferr'=NA))
    return(alternativeFrequency(efs))
  }, mc.cores=cpus))
  return(ret)
}

maxCov = function() {return(get('.maxCov', envir = .GlobalEnv))}

#helper function that calculates the probability that a region has 50% frequency.
redoHetCalculations = function(clusters, eFreqs, plot=F, cpus=1) {
  x = eFreqs$x
  snps = lapply(1:nrow(clusters), function(row) which(x < clusters$x2[row] & x > clusters$x1[row]))

  Z = matrix(rep(0, length(snps)*5), nrow=5)
  Zfreq = mclapply(snps[clusters$cov > 0], function(is) {
    efs = eFreqs[is,]

    fAlt = refBias(alternativeFrequency(efs)['f'])
    
    #z score for each SNP
    z = sapply(1:length(efs$cov), function(i)
      whichDistribution(y1 = twoBinom(0:floor(efs$cov[i]*refBias(0.5)), efs$cov[i], refBias(0.5), refBias(0.5)),
                        y2 = twoBinom(0:floor(efs$cov[i]*refBias(0.5)), efs$cov[i], fAlt, refBiasMirror(fAlt)),
                        x = min(floor(efs$cov[i]*refBias(0.5)), efs$var[i])+1))
    ret = sum(z[2,])/sum(z[1,])
    meanRetNull = sum(z[3,])/sum(z[5,])
    errorRetNull = sqrt(sum(z[7,]-z[3,]^2))/sum(z[5,])
    meanRetAlt = sum(z[4,])/sum(z[6,])
    errorRetAlt = sqrt(sum(z[8,]) - sum(z[4,]^2))/sum(z[6,])
    return(c(ret, meanRetNull, meanRetAlt, errorRetNull, errorRetAlt))
  }, mc.cores=cpus)
  Zfreq = do.call(cbind, Zfreq)
  Z[,clusters$cov > 0] = Zfreq

  if ( plot ) {
    plot(Z[1,], ylim=c(-0.3,0.3), type='n', xlab='region', ylab='mean log likelihood ratio')
    segments(0,0,1000,0)
    x = 1:ncol(Z)
    points(x-0.05, Z[2,], pch=16, col='blue')
    segments(x-0.05, Z[2,]-Z[4,], x-0.05, Z[2,]+Z[4,], pch=16, col='blue')
    points(x+0.05, Z[3,], pch=16, col='red')
    segments(x+0.05, Z[3,]-Z[5,], x+0.05, Z[3,]+Z[5,], pch=16, col='red')
    points(Z[1,], pch=4, cex=1.5, lwd=2)
    legend('topright', c('expected from null', 'expected from CNA', 'measured'), col=c('blue', 'red', 'black'),
           lwd = c(1,1,0.001), pt.lwd=c(1,1,2), pt.cex=c(1,1,1.5), pch=c(16,16,4), bg='white')
    plot(Z[1,], ylim=c(min(Z[1:3,]), max(Z[1:3,])), type='n', xlab='region', ylab='mean log likelihood ratio')
    segments(0,0,1000,0)
    x = 1:ncol(Z)
    points(x-0.05, Z[2,], pch=16, col='blue')
    segments(x-0.05, Z[2,]-Z[4,], x-0.05, Z[2,]+Z[4,], pch=16, col='blue')
    points(x+0.05, Z[3,], pch=16, col='red')
    segments(x+0.05, Z[3,]-Z[5,], x+0.05, Z[3,]+Z[5,], pch=16, col='red')
    points(Z[1,], pch=4, cex=1.5, lwd=2)
    legend('topright', c('expected from null', 'expected from CNA', 'measured'), col=c('blue', 'red', 'black'),
           lwd = c(1,1,0.001), pt.lwd=c(1,1,2), pt.cex=c(1,1,1.5), pch=c(16,16,4), bg='white')
  }
  
  pHet = pnorm(Z[1,], Z[2,], Z[4,])
  pHet[is.nan(pHet)] = 0.5
  pAlt = pnorm(-Z[1,], -Z[3,], Z[5,])
  pAlt[is.nan(pAlt)] = 0.5
  
  clusters$pHet = pHet
  clusters$pAlt = pAlt
  clusters$postHet = 0.5*pHet/(0.5*pHet+(1-0.5)*pAlt)  #prior belief that 50% of the regions are normal
  clusters$stat = Z[1,]
  clusters$nullStat = Z[2,]
  clusters$altStat = Z[3,]
  clusters$nullStatErr = Z[4,]
  clusters$altStatErr = Z[5,]

  return(clusters)
}

#between two input function y1 and y2 on some common set (y1 and y2 must be same length)
#if you pick x, what is the likleyhood ratio that it came from y1, not y2? Also returns weight, indicating how
#different the functions are, ie how much power x has to differentiate the two functions.
whichDistribution = function(y1, y2, x) {
  y1 = y1/sum(y1)
  y2 = y2/sum(y2)
  ord1 = order(y1)
  y1 = y1[ord1]
  y2 = y2[ord1]
  x = which(ord1==x)
  p1 = cumsum(y1)
  ord2 = order(y2)
  y1 = y1[ord2]
  y2 = y2[ord2]
  p1 = p1[ord2]
  x = which(ord2==x)
  p2 = cumsum(y2)
  
  lods = log(y1/y2)
  weight = ifelse(pmax(p1, p2) > 0.05, 1, 0)
  retStat = weight*lods
  
  meanStatNull = sum(y1*retStat)
  meanWeightNull = sum(y1*weight)
  meanStatAlt = sum(y2*retStat)
  meanWeightAlt = sum(y2*weight)
  meanStat2Null = sum(y1*retStat^2)
  meanStat2Alt = sum(y2*retStat^2)
    
  return(c('weight'=weight[x], 'stat'=retStat[x], 'meanStatNull'=meanStatNull, 'meanStatAlt'=meanStatAlt,
           'meanWeightNull'=meanWeightNull, 'meanWeightAlt'=meanWeightAlt,
           'meanStat2Null'=meanStat2Null, 'meanStat2Alt'=meanStat2Alt))
}

#this functions estimates the variance of total statistic due to the addition of this observation.
#to be combined into a single variance later on.
meanDistribution = function(y1, y2, Wnull, Walt, totNull, totAlt) {
  y1 = y1/sum(y1)
  y2 = y2/sum(y2)
  ord1 = order(y1)
  y1 = y1[ord1]
  y2 = y2[ord1]
  p1 = cumsum(y1)
  ord2 = order(y2)
  y1 = y1[ord2]
  y2 = y2[ord2]
  p1 = p1[ord2]
  p2 = cumsum(y2)
  
  weight = pmax(p1, p2)
  stat = log(y1/y2)
  
  meanStatNull = sum(y1*weight*log(y1/y2))
  meanStatAlt = sum(y2*weight*log(y1/y2))
  totStatNull = (y1*weight*log(y1/y2) + tot)/(y1*weight + W)
  varTotNull = mean(totStatNull^2) - mean(totStatNull)^2

    mean(y1^2*weight*log(y1/y2)^2) - mean(y1*weight*log(y1/y2))^2
  varRetAlt = mean(y2^2*weight*log(y1/y2)^2) - mean(y2*weight*log(y1/y2))^2
  
  return(c('weight'=weight, 'stat'=retStat, 'meanStatNull'=meanRetNull, 'meanStatAlt'=meanRetAlt,
           'varStatNull'=varRetNull, 'varStatAlt'=varRetAlt))
}


#remormalise the coverage so that the AB calls are average of 0.
normaliseCoverageToHets = function(clusters) {
  is = which(clusters$call == 'AB' & abs(clusters$M) < 0.3)
  if ( length(is) > 0 )
    meanM = sum((clusters$M/clusters$width^2)[is])/sum(1/clusters$width[is]^2)
  else meanM = 0
  catLog('Shifting overall normalisation to correct for a mean LFC of', meanM, 'in AB calls.\n')
  clusters$M = clusters$M - meanM
  return(list(clusters=clusters, meanM=meanM))
}

#calls the allelic copy number and clonality of the region.
addCall = function(clusters, eFreqs) {
  catLog('Calling CNVs in clustered regions..')
  for ( row in 1:nrow(clusters) ) {
    clusters$call[row] = '???'
    clusters$clonality[row] = 0
    clusters$clonalityError[row] = Inf
    clusters$sigma[row] = 3
    efs = eFreqs[eFreqs$x > clusters$x1[row] & eFreqs$x < clusters$x2[row],]
    isab = isAB(clusters[row,], efs, sigmaCut=clusters$sigma[row])
    clusters$call[row] = 'AB'
    clusters$clonality[row] = 1
    clusters$clonalityError[row] = 0
    clusters$sigma[row] = isab$sigma
    clusters$pCall[row] = clusters$postHet[row]
    if ( !(isab$call) ) {
      for ( tryCall in allCalls() ) {
        iscnv = isCNV(clusters[row,], efs, callTofM(tryCall)['M'], callTofM(tryCall)['f'], callPrior(tryCall),
          sigmaCut=max(3, clusters$sigma[row]))
        if ( (iscnv$call & clusters$sigma[row] > 3) |
            (iscnv$call & (iscnv$clonality > clusters$clonality[row] | (iscnv$clonality == clusters$clonality[row] & iscnv$sigma < clusters$sigma[row]))) ) {
          clusters$clonality[row] = iscnv$clonality
          clusters$clonalityError[row] = iscnv$clonalityError
          clusters$sigma[row] = iscnv$sigma
          clusters$call[row] = tryCall
          clusters$pCall[row] = iscnv$pCall
        }
      }
    }
  }
  aBitWeird = clusters$sigma > 5 | clusters$clonalityError > clusters$clonality/2 | clusters$pCall < 1e-4
  veryWeird = clusters$sigma > 10 | clusters$clonalityError > clusters$clonality | clusters$pCall < 1e-8
  clusters$call[aBitWeird] = paste0(clusters$call[aBitWeird], '?')
  clusters$call[veryWeird] = paste0(clusters$call[veryWeird], '?')
  catLog('done!\n')
  return(clusters)
}

#probability that a region is AB.
isAB = function(cluster, eFreqs, sigmaCut=3) {
  if ( sum(eFreqs$cov) == 0 ) pF = 0.5
  else if ( cluster$stat >= 0 ) pF = 1
  else {
    pF = cluster$postHet
  }
  pM = 2*pt(-abs(cluster$M)/(cluster$width+systematicVariance()), df=cluster$df)  #allow systematic effects
  pBoth = sapply(1:length(pF), function(row) fisherTest(c(pF[row], pM[row]))[2])
  sigma = abs(qnorm(pBoth/2, 0, 1))
  return(list(call=sigma < sigmaCut, sigma = sigma, clonalityError=0))
}

#the considered calls in the algorithm. (CL is complete loss, ie loss of both alleles.)
allCalls = function() {
  return(c('AB', 'A', 'AA', 'AAA', 'AAAA', 'AAB', 'AAAB', 'AAAAB', 'AAAAAB', 'AAAAAAB', 'AABB', 'CL'))
}

#returns the prior of a call. Prior is proprotional to 1 divided by the number of removed or added chromosomes.
#no CNV is given a prior 5 times as high as A and AAB.
callPrior = function(call) {
  priors = c('AB'=5, 'A'=1, 'AA'=1/2, 'AAA'=1/3, 'AAAA'=1/4, 'AAB'=1, 'AAAB'=1/2, 'AAAAB'=1/3, 'AAAAAB'=1/4, 'AAAAAAB'=1/5,
    'AABB'=1/2, 'CL'=1/2)
  priors = priors/sum(priors)
  if ( call %in% names(priors) ) return(priors[call])
  else return(0.1)
}

#helper function that returns the expected coverage and frequency of a CNV call.
callTofM = function(call) {
  nA = nchar(call) - nchar(gsub('A', '', call))
  nB = nchar(call) - nchar(gsub('B', '', call))
  f = min(nA, nB)/(nA+nB)
  if ( call %in%  c('CL', '') ) f = 0.5  #if complete loss, assume not completely clonal and a normal AB background.
  return(c(f=f, M=(nA+nB)/2))
}

#estimates how likely a certain call is for a given region.
isCNV = function(cluster, efs, M, f, prior, sigmaCut=3) {
  #set an estimate of the error on the measured MAF in the cluster
  ferr = if ( is.nan(cluster$f) ) Inf else cluster$ferr

  #set weights for coverage and frequency in determining clonality.
  linearMWidth = abs(2^(cluster$M+cluster$width) - 2^(cluster$M-cluster$width))/2 + 2^cluster$M*systematicVariance()
  if ( is.na(cluster$f) ) cf = 0 else cf = cluster$f

  #If no basis to determine clonality, the no point in continuing
  if ( (f == 0.5 | is.na(cluster$f)) & M == 1 )
    return(list(call=F, clonality=0, sigma = Inf, clonalityError=Inf, pCall=1))

  #estimate clonality as a weighted mean from the coverage and frequency
  freqClonality = (0.5-cf)/(0.5-cf+M*(cf-f))
  clonalityErrorF = ferr*(freqClonality*(M-1)+1)/(0.5-cf+M*(cf-f))
  if ( f == 0.5 | is.na(cluster$f) ) {
    freqClonality = 1
    clonalityErrorF = Inf
  }
  if ( f == 0.5 | is.na(cluster$f) ) fweight = 0
  else fweight = 1/clonalityErrorF^2
  covClonality = (2^cluster$M - 1)/(M - 1)
  clonalityErrorM = linearMWidth/abs(M-1)
  if ( M == 1 ) {
    covClonality = 1
    clonalityErrorM = Inf
  }
  if ( !is.na(clonalityErrorM) ) Mweight = 1/clonalityErrorM^2
  else Mweight = 0
  clonality = (freqClonality*fweight + covClonality*Mweight)/(fweight+Mweight)
  clonality = min(1.05, clonality)

  #propagate errors on f and M to the clonality. Give bonus uncertainty if f and M disagree.
  clonalityDifference =
    if ( fweight > 0 & Mweight > 0 ) noneg(abs(freqClonality-covClonality) - sqrt(clonalityErrorF^2 + clonalityErrorM^2))
    else 0
  if ( fweight == 0 ) {     #if information from only one source, add some extra uncertainty to the clonality
    clonalityErrorF = 0.3
    fweight = sqrt(Mweight)
  }
  if ( Mweight == 0 ) {
    clonalityErrorM = 0.3
    Mweight = sqrt(fweight)
  }
  clonalityError =
    sqrt(fweight^2*clonalityErrorF^2 + Mweight^2*clonalityErrorM^2 + fweight*Mweight*clonalityDifference^2)/(fweight + Mweight)
  #if negative clonality, add uncertainty to the clonality and force back to 0.
  if ( clonality < 0 ) {
    clonalityError = sqrt(clonalityError^2 + abs(clonality)^2)
    clonality = 0
  }
  
  #update f and M for the estimated clonality
  fClone = (clonality*(2*f*M-1) + 1)/(clonality*(2*M - 2) + 2)
  MClone = log2(1 + (M-1)*clonality)
  #calculate sigma for these f and M
  if ( cluster$cov == 0 ) pF = pCall = exp(-1)
  else {
    pCall = pF = fisherTest(pBinom(efs$cov, efs$var, refBias(fClone)))['pVal']
    pF = (pF+1e-10)
  }
  pM = 2*pt(-noneg(abs(cluster$M - MClone))/(cluster$width+systematicVariance()), df=cluster$df)  #allow systematic error
  pBoth = fisherTest(c(pF, pM))[2]
  if ( cluster$cov == 0 ) pBoth = pM
  
  #add prior
  pBoth = prior*pBoth/(prior*pBoth + (1-prior)*(1-pBoth))
  sigma = abs(qnorm(pBoth/2, 0, 1))

  call = sigma < sigmaCut & clonality > 0.1
  if ( is.na(call) ) call = F
  return(list(call=as.logical(call), clonality=pmin(1,as.numeric(clonality)),
              sigma = sigma, clonalityError=as.numeric(clonalityError), pCall=pCall))
}

#groups up CNV regions with similar clonalities to estimate the subclonal structure of the sample
findSubclones = function(cR) {
  catLog('Finding subclones..')  
  clones = order(-cR$clonality)
  clones = clones[cR$clonality[clones] < 1]
  regions = as.list(clones)
  clonality = cR$clonality[clones]
  error = cR$clonalityError[clones]
  while( length(clonality) > 1 ) {
    pairScore = sameClone(clonality, error)
    best = which(pairScore == max(pairScore))[1]

    if ( pairScore[best] < 0.05 ) break
    clonality[best] = sum((clonality/error^2)[best:(best+1)])/sum(1/error[best:(best+1)]^2)
    error[best] = 1/sqrt(sum(1/error[best:(best+1)]^2))
    regions[[best]] = unlist(c(regions[best], regions[best+1]))
    regions = regions[-(best+1)]
    clonality = clonality[-(best+1)]
    error = error[-(best+1)]
  }
  subclonality = rep(1, nrow(cR))
  subclonalityError = rep(0, nrow(cR))
  if ( length(regions) > 0 ) {
    for ( i in 1:length(regions) ) {
      subclonality[regions[[i]]] = clonality[i]
      subclonalityError[regions[[i]]] = error[i]
    }
  }
  cR = cbind(cR, subclonality, subclonalityError)
  catLog('done!\n')  
  return(cR)
}

#helper function that subsets fit objects properly, with all the extra columns.
subsetFit = function(fit, rows=NA, cols=NA) {
  if ( is.na(cols)[1] ) cols = 1:ncol(fit)
  if ( is.na(rows)[1] ) rows = 1:nrow(fit)
  fit = fit[rows,cols]
  if ('best.guess' %in% names(fit) ) fit$best.guess = fit$best.guess[rows,cols, drop=F]
  if ('posterior' %in% names(fit) ) fit$posterior = lapply(fit$posterior[cols], function(post) post[rows,,drop=F])
  if ('prior' %in% names(fit) ) fit$prior[cols]
  if ('XRank' %in% names(fit) )   fit$XRank = fit$XRank[rows, cols, drop=F]
  if ('postWidth' %in% names(fit) ) fit$postWidth = fit$postWidth[rows, cols, drop=F]
  if ('x' %in% names(fit) ) fit$x = fit$x[rows]
  if ('x1' %in% names(fit) ) fit$x1 = fit$x1[rows]
  if ('x2' %in% names(fit) ) fit$x2 = fit$x2[rows]
  if ('chr' %in% names(fit) ) fit$chr = fit$chr[rows]
  if ('longNames' %in% names(fit) ) fit$longNames = fit$longNames[rows]
  if ('sex' %in% names(fit) ) fit$sex = fit$sex[cols]

  return(fit)
}

systematicVariance = function() {return(get('.systematicVariance', envir = .GlobalEnv))}

#the posterior probability that two capture regions belong the same CNVinterval
sameCNV = function(cR) {
  if ( dim(cR)[1] < 2 ) return(1)

  #find the prior from the gap lengths
  first = 1:(nrow(cR)-1)
  second = 2:nrow(cR)
  dx = abs((cR$x1+cR$x2)[second]/2 - (cR$x1+cR$x2)[first]/2)
  prior = exp(-dx*CNVregionsPerBP())

  #find the probabilities of getting measure values if the SNP frequencies are equal
  #a fisher test would be more appropriate, but can be computationally slow. These two binomals give very similar results.
  f = (cR$var[first] + cR$var[second])/(cR$cov[first] + cR$cov[second])
  p1 = pBinom(cR$cov[first], cR$var[first], f)
  p2 = pBinom(cR$cov[second], cR$var[second], f)
  simesFP = sapply(first, function(i) min(p.adjust(c(p1[i], p2[i]), method='fdr')))
  
  #get posterior using the dx-based prior above. alternative hypothesis is a flat distribution on the frequency.
  #include the possibility of both being 50% hets, in which case a different variance (coverage) between the regions
  #can cause different mean MAF. We catch that case by grouping regions that are consistent with 50%, independently of MAF.
  if ( 'postHet' %in% names(cR) ) pBoth50 = sapply(first, function(i) min(p.adjust(c(cR$postHet[i], cR$postHet[i+1]), method='fdr')))
  else pBoth50 = sapply(first, function(i) min(p.adjust(c(cR$pHet[i], cR$pHet[i+1]), method='fdr')))
  postBoth50 = prior*pBoth50/(prior*pBoth50 + (1-prior))
  pF = pmax(simesFP, pBoth50)
  noFreq = cR$cov[first] == 0 | cR$cov[second] == 0

  #find the probability densities (P*dM) of getting measure values if the regions have the same fold change
  width = cR$width + 0*systematicVariance()
  meanM = (cR$M[first]/width[first]^2 + cR$M[second]/width[second]^2)/(1/width[first]^2 + 1/width[second]^2)
  MP1 = pt(-abs(cR$M[first] - meanM)/width[first], df = cR$df[first])
  MP2 = pt(-abs(cR$M[second] - meanM)/width[second], df = cR$df[second])
  pM = sapply(first, function(i) min(p.adjust(c(MP1[i], MP2[i]), method='fdr')))

  pMandF = sapply(first, function(i) min(p.adjust(c(pF[i], pM[i]), method='fdr')))

  #get posteriors, and take average
  post = prior*pMandF/(prior*pMandF + (1-prior))
  names(post) = cR$x1[first]

  return(post)
}

#prior of density of CNV region breakpoints. This corresponds to around 30 breakpoints in a sample.
CNVregionsPerBP = function() {return(1/1e8)}

#helper function doing the stouffer Test
stoufferTest = function(p, w) {
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c(Z = Z, p.value = p.val))
}

#helper function calculating the posterior of clonalities being the same
sameClone = function(clonality, error, prior = 0.99) {
  first = 1:(length(clonality)-1)
  second = 2:length(clonality)
  
  totalError = sqrt(error[first]^2 + error[second]^2)
  sigma = abs(clonality[first]-clonality[second])/totalError
  pval = 2*pnorm(-sigma, mean=0, sd=1)
  #alternative hypothesis is a flat distribution on difference in clonality from -0.5 to 0.5.
  post = pval*prior/(pval*prior + totalError*(1-prior))

  return(post)
}

#helper function that use neighbours to study the variance estimate
#will add a constant to width to correct for underrestimated variance
boostCRwidth = function(CR) {
  CR = CR[order(CR$x1),]
  x = (CR$x1+CR$x2)/2
  
  x = (-1000:1000)/100
  y1 = rt(100000, CR$df[1])
  y2 = rt(100000, CR$df[1])
  the = median(abs(y1-y2))
  breaks = (0:500)/5

  last = nrow(CR)
  changeBoost = function(boost) median(abs((CR$M[-1]-CR$M[-last])/sqrt((CR$width[-1]+boost)^2+(CR$width[-last]+boost)^2)))
  changes = sapply((0:50)/1000, changeBoost)
  widthBoost = ((0:50)/1000)[max(c(1, which(changes > the)))]
  CR$width = CR$width + widthBoost
  
  return(CR)
}
