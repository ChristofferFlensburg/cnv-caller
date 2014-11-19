


#This function takes the cancer and normal variants, and the coverage analysis as input
#It tries to find the germline het variants for each sample.
#It then clusters regions of the genome that have similar het frequencies and coverage.
#It then calls the CN of each clustered region, as well as clonality, uncertainty estimate and p-value
#for normal CN.
#returns a list of data frames with each region of the genome on a row.
callCNVs = function(variants, normalVariants, fitS, names, individuals, normals, Rdirectory, genome='hg19', cpus=1, v='', forceRedoCNV=F) {
  clustersSaveFile = paste0(Rdirectory, '/clusters.Rdata')
  if ( file.exists(clustersSaveFile) & !forceRedoCNV ) {
    catLog('Loading saved CNV results.\n')
    load(file=clustersSaveFile)
    return(clusters)
  }

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
                                  genome=genome, v=v, cpus=cpus))
    else
      return(callCancerNormalCNVs(cancerVariants=variants$variants[[name]],
                            normalVariants=FALSE,
                            moreNormalVariants=normalVariants$variants,
                            fit = subsetFit(fitS, cols=paste0(name, '-normal')),
                                  genome=genome, v=v, cpus=cpus))
  })

  names(clusters) = names
  save(clusters, file=clustersSaveFile)
  return(clusters)
}


#the high level function that controls the steps of the CNV calling for given sample and normal variant objects.
callCancerNormalCNVs = function(cancerVariants, normalVariants, moreNormalVariants, fit, genome='hg19', v='', cpus=1) {
  require(parallel)
  #estimate reference bias and variance from the selected normal hets.
  if ( class(normalVariants) == 'logical') setVariantLoss(moreNormalVariants, v=v)
  else setVariantLoss(normalVariants, v=v)

  #select good germline het variants from normals:
  if ( class(normalVariants) == 'logical') use = selectGermlineHetsFromCancer(cancerVariants, moreNormalVariants, v=v, cpus=cpus)
  else use = selectGermlineHets(normalVariants, moreNormalVariants, v=v, cpus=cpus)
  is = rownames(cancerVariants) %in% use
  cancerVariants = cancerVariants[is,]
  
  #summarise by capture region
  catLog('Summarising capture regions..')
  cancerCR = unifyCaptureRegions(cancerVariants, fit, cpus=cpus)
  catLog('done!\n')

  #run clustering algorithm
  cancerCluster = mergeChromosomes(cancerCR, genome=genome, v=v, cpus=cpus)

  #run post processing, such as correcting normalisation from AB regions, calling CNVs and clonalities.
  catLog('Postprocessing...')
  cancerFreqs = data.frame(var=mirrorDown(cancerVariants$var, cov=cancerVariants$cov),
    cov=cancerVariants$cov, x=cancerVariants$x)
  post = postProcess(cancerCluster, cancerCR, cancerFreqs, genome, v=v)
  cancerCluster = post$clusters
  cancerCR$M = cancerCR$M - post$meanM
  #plotCR(cancerCluster)

  #return clustered regions with calls, as well as the raw capture region data.
  catLog('done!\n')
  return(list(clusters=cancerCluster, CR=cancerCR,
              freqs=data.frame(var=cancerVariants$var, cov=cancerVariants$cov, x=cancerVariants$x)))
}


#helper function that selects germline het SNPs in the presence of a normal sample from the same individual.
selectGermlineHets = function(normalVariants, moreNormalVariants, minCoverage = 10, v='', cpus=1) {
  #only bother with variants that have enough coverage so that we can actually see a change in frequency
  if ( v == 'trackProgress' ) cat('Taking variants with minimum coverage of', minCoverage, '...')
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
  
  #pick het in matching normal sample
  catLog('Taking variants that are het in the normal sample..')
  normalF = normalVariants$var/normalVariants$cov
  normalHet = pBinom(normalVariants$cov, normalVariants$var, refBias(0.5)) > 0.1 & abs(normalF-0.5) < 0.25 & normalVariants$flag == ''
  use = use[normalHet]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(normalHet), 'variants.\n')
  if ( length(use) == 0 ) return(use)

  #keep only variants that are present in analysis of other normals
  if ( v == 'trackProgress' ) cat('Taking variants that are analysed for all normals..')
  use = Reduce(intersect, lapply(moreNormalVariants, rownames), rownames(normalVariants))
  moreNormalVariants = lapply(moreNormalVariants, function(vs) {
    is = rownames(vs) %in% use
    return(vs[is,])
  })
  normalVariants = normalVariants[use,]
  catLog('done! Got', length(use), 'variants.\n')
  if ( length(use) == 0 ) return(use)
  
  #to avoid chosing noise that happens to be around f=0.5, require at least one normal that is not het.
  catLog('Taking variants that are polymorphic..')
  norVar = do.call(cbind, lapply(moreNormalVariants, function(vs) vs[use,]$var))
  norCov = do.call(cbind, lapply(moreNormalVariants, function(vs) vs[use,]$cov))
  is0 = norVar/norCov < 0.1 #too strict with many normals?
  is1 = norVar/norCov > 0.9
  isHet = abs((norVar-norCov*refBias(0.5))/sqrt((norCov+1)*refBias(0.5))) < 3
  rownames(isHet) = rownames(is0) = rownames(is1) = use
  polymorphic = rowSums(!isHet & (is0 | is1)) > 1
  use = use[polymorphic]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(polymorphic), 'variants.\n')
  if ( length(use) == 0 ) return(use)
  
  #again to filter out noisy variants, require that all normals are conistent with ref, het or hom.
  if ( v == 'trackProgress' ) cat('Taking variants that are ref, het or hom in all normals..')
  consistent = rowSums(!(is0[use,,drop=F] | is1[use,,drop=F] | isHet[use,,drop=F])) == 0
  use = use[consistent]
  normalVariants = normalVariants[use,]
  catLog('done! Got', sum(consistent), 'variants.\n')
  if ( length(use) == 0 ) return(use)

  #again to filter out noisy variants, filter on the RIB statistic in the cancer sample
  catLog('Taking variants that have low expected rate of incorrect basecalls in the sample..')
  highCancerRIB = cancerVariants$RIB > 0.03
  use = use[!highCancerRIB]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(!highCancerRIB), 'variants.\n')

  return(use)
}

#helper function that selects germline het SNPs in the absence of a normal sample from the same individual.
selectGermlineHetsFromCancer = function(cancerVariants, moreNormalVariants, minCoverage = 10, v='', cpus=1) {
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
  
  #pick het in matching normal sample
  catLog('Taking unflagged cancer variants that have frequency above 5% and below 95%..')
  cancerF = cancerVariants$var/cancerVariants$cov
  normalHet = abs(cancerF-0.5) < 0.45 & cancerVariants$flag == ''
  use = use[normalHet]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(normalHet), 'variants.\n')

  #keep only variants that are present in analysis of other normals
  catLog('Taking variants that are analysed for all normals..')
  use = Reduce(intersect, lapply(moreNormalVariants, rownames), rownames(cancerVariants))
  moreNormalVariants = lapply(moreNormalVariants, function(vs) {
    is = rownames(vs) %in% use
    return(vs[is,])
  })
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', length(use), 'variants.\n')
  
  #to avoid chosing noise that happens to be around f=0.5, require at least one normal that is not het.
  catLog('Taking variants that are polymorphic..')
  norVar = do.call(cbind, lapply(moreNormalVariants, function(vs) vs[use,]$var))
  norCov = do.call(cbind, lapply(moreNormalVariants, function(vs) vs[use,]$cov))
  is0 = norVar/norCov < 0.03 #stricter cut here than in normal, to avoid noisy germline refs around 5-10%
  is1 = norVar/norCov > 0.97
  isHet = abs((norVar-norCov*refBias(0.5))/sqrt((norCov+1)*refBias(0.5))) < 3
  rownames(isHet) = rownames(is0) = rownames(is1) = use
  polymorphic = rowSums(!isHet & (is0 | is1)) > 1
  use = use[polymorphic]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(polymorphic), 'variants.\n')
  
  #again to filter out noisy variants, require that all normals are conistent with ref, het or hom.
  catLog('Taking variants that are ref, het or hom in all normals..')
  consistent = rowSums(!(is0[use,] | is1[use,] | isHet[use,])) == 0
  use = use[consistent]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(consistent), 'variants.\n')

  #again to filter out noisy variants, filter on the RIB statistic in the cancer sample
  catLog('Taking variants that have low expected rate of incorrect basecalls in the sample..')
  highCancerRIB = cancerVariants$RIB > 0.03
  use = use[!highCancerRIB]
  cancerVariants = cancerVariants[use,]
  catLog('done! Got', sum(!highCancerRIB), 'variants.\n')

  return(use)
}

unifyCaptureRegions = function(variants, fit, cpus=1) {
  #group SNPs by capture region
  x = variants$x
  snps = lapply(1:nrow(fit), function(row) which(x < fit$x2[row] & x > fit$x1[row]))

  #mirror down to MAFs
  variants$var = mirrorDown(variants$var, variants$cov)

  uniFreq = do.call(rbind, lapply(snps, function(is) c(sum(variants$var[is]), sum(variants$cov[is]))))
  colnames(uniFreq) = c('var', 'cov')
  
  pHet = rep(0.5, length(snps))
  pHet[uniFreq[,'cov'] > 0] = unlist(mclapply(snps[uniFreq[,'cov'] > 0], function(is) fisherTest(pBinom(variants$cov[is], variants$var[is], refBias(0.5)))[2], mc.cores=cpus)) 

  pAlt = rep(0.5, length(snps))
  pAlt[uniFreq[,'cov'] > 0] = unlist(lapply(snps[uniFreq[,'cov'] > 0], function(is)
        fisherTest(pTwoBinom(variants$cov[is], variants$var[is], sum(variants$var[is])/sum(variants$cov[is])))[2]))

  odsHet = rep(0.5, length(snps))
  odsHet[uniFreq[,'cov'] > 0] = unlist(lapply(snps[uniFreq[,'cov'] > 0], function(is) {
    fAlt = sum(variants$var[is])/sum(variants$cov[is])
    dpHet = dbinom(variants$var[is], variants$cov[is], refBias(0.5))
    dpAlt = twoBinom(variants$var[is], variants$cov[is], fAlt, refBiasMirror(fAlt))
    odsHet = dpHet/(dpHet+dpAlt)
    return(fisherTest(odsHet)[2])
  }))

  cR = data.frame(x1=fit$x1, x2=fit$x2, M=fit$coefficients[,1], width=fit$coefficients[,1]/fit$t[,1], df=fit$df.total,
    uniFreq[,1:2], pHet=pHet, pAlt=pAlt, odsHet=odsHet)
  return(cR)
}


#merges regions in each chromosome.
mergeChromosomes = function(cR, genome='hg19', v='', cpus=1, ...) {
  chrs = xToChr(cR$x1, genome=genome)
  #clusters = list()
  catLog('Merging capture regions with same coverage and MAF: ')

  clusters = mclapply(unique(chrs), function(chr) {
    catLog(chr, '..', sep='')
    return(mergeRegions(cR[chrs == chr,], ...))
  }, mc.cores=cpus)
  
  clusters = do.call(rbind, clusters)
  catLog('done!\n')
  return(clusters)
}
mergeRegions = function(cR, minScore = 0.05, plot=F, debug=F) {
  merged = c()
  scores = c()
  if ( debug ) Nloop = 0
  pairScore = sameCNV(cR) #update to only refresh the affected row/column
  while(dim(cR)[1] > 1) {
    #find the pair of regions with the largest pairing probability
    best = which(pairScore == max(pairScore))[1]
    if ( debug && Nloop < 1 ) {
      catLog('DEBUG: Merging pair at ', cR$x2[best], ' with score ', pairScore[best], ', ', nrow(cR), ' regions.\n', sep='')
      print (cR)
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
      layout(matrix(1:2, nrow=1))
      plotCR(cR)
    }

    #break loop if no clusters that are sufficiently similar
    if ( pairScore[best] < minScore ) break
    #otherwise save the merged regions
    merged = c(merged, best)
    scores = c(scores, pairScore[best])
    
    #merge the regions, using the provided weight
    cR$x2[best] = cR$x2[best+1]
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

#takes clustered regions and post processes them by calling copy numbers and clonality.
postProcess = function(clusters, cRs, freqs, genome='hg19', v='') {
  clusters$f = refUnbias(clusters$var/clusters$cov)
  clusters = redoHetCalculations(clusters, freqs)
  renorm =  normaliseCoverageToHets(clusters)
  clusters = renorm$clusters
  clusters = addCall(clusters, freqs, v=v)
  clusters = findSubclones(clusters, v=v)
  rownames(clusters) = make.names(paste0('chr', xToChr(clusters$x1, genome=genome)), unique=T)
  clusters = clusters[order(clusters$x1),]
  return(list(clusters=clusters, meanM=renorm$meanM))
}

#helper function that fixes the probability that a region has 50% frequency.
redoHetCalculations = function(clusters, freqs) {
  x = freqs$x
  snps = lapply(1:nrow(clusters), function(row) which(x < clusters$x2[row] & x > clusters$x1[row]))

  pHet = rep(0.5, length(snps))
  pHet[clusters$cov > 0] = unlist(mclapply(snps[clusters$cov > 0], function(is) fisherTest(pBinom(freqs$cov[is], freqs$var[is], refBias(0.5)))[2], mc.cores=cpus)) 

  pAlt = rep(0.5, length(snps))
  pAlt[clusters$cov > 0] = unlist(lapply(snps[clusters$cov > 0], function(is)
        fisherTest(pTwoBinom(freqs$cov[is], freqs$var[is], sum(freqs$var[is])/sum(freqs$cov[is])))[2]))

  odsHet = rep(0.5, length(snps))
  odsHet[clusters$cov > 0] = unlist(lapply(snps[clusters$cov > 0], function(is) {
    fAlt = sum(freqs$var[is])/sum(freqs$cov[is])
    dpHet = dbinom(freqs$var[is], freqs$cov[is], refBias(0.5))
    dpAlt = twoBinom(freqs$var[is], freqs$cov[is], fAlt, refBiasMirror(fAlt))
    odsHet = dpHet/(dpHet+dpAlt)
    return(fisherTest(odsHet)[2])
  }))

  clusters$pHet = pHet
  clusters$pAlt = pAlt
  clusters$odsHet = odsHet

  return(clusters)
}

#remormalise the coverage so that the AB calls are average of 0.
normaliseCoverageToHets = function(clusters) {
  is = which(clusters$pHet > 0.05 & abs(clusters$M) < 0.3)
  if ( length(is) > 0 )
    meanM = sum((clusters$M/clusters$width^2)[is])/sum(1/clusters$width[is]^2)
  else meanM = 0
  clusters$M = clusters$M - meanM
  return(list(clusters=clusters, meanM=meanM))
}

#calls the allelic copy number and clonality of the region.
addCall = function(clusters, freqs, v='') {
  catLog('Calling CNVs in clustered regions..')
  for ( row in 1:nrow(clusters) ) {
    clusters$call[row] = '???'
    clusters$clonality[row] = 0
    clusters$clonalityError[row] = Inf
    clusters$sigma[row] = 3
    fs = freqs[freqs$x > clusters$x1[row] & freqs$x < clusters$x2[row],]
    isab = isAB(clusters[row,], fs, sigmaCut=clusters$sigma[row])
    clusters$call[row] = 'AB'
    clusters$clonality[row] = 1
    clusters$clonalityError[row] = 0
    clusters$sigma[row] = isab$sigma
    clusters$pCall[row] = clusters$pHet[row]
    if ( !(isab$call) ) {
      for ( tryCall in allCalls() ) {
        iscnv = isCNV(clusters[row,], fs, log2(callTofM(tryCall)['M']), callTofM(tryCall)['f'], callPrior(tryCall),
          sigmaCut=max(3, clusters$sigma[row]))
        if ( (iscnv$call & clusters$sigma[row] > 3) | (iscnv$call & iscnv$clonality > clusters$clonality[row]) ) {
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
isAB = function(cluster, freqs, sigmaCut=3) {
  if ( sum(freqs$var)/sum(freqs$cov) > refBias(0.5) - 0.03 | sum(freqs$cov) == 0 ) pF = 1
  else {
    pF = cluster$odsHet
  }
  pM = 2*pt(-abs(cluster$M)/(cluster$width+0.05), df=cluster$df)  #allow 0.05 off from systematic effects
  pBoth = sapply(1:length(pF), function(row) fisherTest(c(pF[row], pM[row]))[2])
  sigma = abs(qnorm(pBoth/2, 0, 1))
  return(list(call=sigma < sigmaCut, sigma = sigma, clonalityError=0))
}

#the considered calls in the algorithm. (CL is complete loss, ie loss of both alleles.)
allCalls = function() {
  return(c('A', 'AA', 'AAA', 'AAAA', 'AAB', 'AAAB', 'AAAAB', 'AAAAAB', 'AAAAAAB', 'AABB', 'CL'))
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
isCNV = function(cluster, freqs, M, f, prior, sigmaCut=3) {
  #set an estimate of the error on the measured MAF in the cluster
  ferr = if ( is.nan(cluster$f) ) Inf else sqrt(cluster$var+1)/cluster$cov

  #set weights for coverage and frequency in determining clonality.
  fweight = abs(0.5-f)/ferr
  Mweight = abs(M)/(cluster$width+0.05)
  if ( is.nan(cluster$f) ) cf = 0 else cf = cluster$f

  #If no basis to determine clonality, the no point in continuing
  if ( fweight == 0 & Mweight == 0 )
    return(list(call=F, clonality=0, sigma = Inf, clonalityError=Inf, pCall=1))

  #estimate clonality as a weighted mean.
  clonality = ((cf-0.5)/(f-0.5)*fweight + cluster$M/M*Mweight)/(fweight+Mweight)
  if ( f == 0.5 ) clonality = cluster$M/M
  if ( M == 0 ) clonality = (cf-0.5)/(f-0.5)
  clonality = min(1, clonality)

  #propagate errors on f and M to the clonality. Give bonus uncertainty if f and M disagree.
  clonalityErrorF = ferr/(0.5-f)
  clonalityErrorM = (cluster$width+0.05)/M   #assume 0.05 systematic error on baseline
  clonalityDifference =
    if ( fweight > 0 & Mweight > 0 ) abs(((cf-0.5)/(f-0.5)-clonality)*fweight - (cluster$M/M-clonality)*Mweight)/(fweight+Mweight)
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
    sqrt(fweight^2*clonalityErrorF^2 + Mweight^2*clonalityErrorM^2 + fweight*Mweight*clonalityDifference^2)/
      sqrt(fweight^2 + Mweight^2)
  #if negative clonality, add uncertainty to the clonality and force back to 0.
  if ( clonality < 0 ) {
    clonalityError = sqrt(clonalityError^2 + abs(clonality)^2)
    clonality = 0
  }
  
  #update f and M for the estimated clonality
  fClone = 0.5 - (0.5-f)*clonality
  MClone = M*clonality
  #calculate sigma for these f and M
  if ( cluster$cov == 0 ) pF = pCall = exp(-1)
  else {
    pCall = pF = fisherTest(pBinom(freqs$cov, freqs$var, refBias(fClone)))['pVal']
    pF = (pF+1e-10)   #50% prior belief that the frequency is f
  }
  pM = 2*pt(-noneg(abs(cluster$M - MClone) -0.05)/cluster$width, df=cluster$df)  #allow 0.05 systematic error
  pBoth = sapply(1:length(pF), function(row) fisherTest(c(pF[row], pM[row]))[2])
  #add prior
  pBoth = prior*pBoth/(prior*pBoth + (1-prior)*(1-pBoth))
  sigma = abs(qnorm(pBoth/2, 0, 1))


  call = sigma < sigmaCut & clonality > 0.1
  if ( is.na(call) ) call = F
  return(list(call=as.logical(call), clonality=as.numeric(clonality),
              sigma = sigma, clonalityError=as.numeric(clonalityError), pCall=pCall))
}

#groups up CNV regions with similar clonalities to estimate the subclonal structure of the sample
findSubclones = function(cR, v='') {
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

  return(fit)
}


#the posterior probability that two capture regions belong the same CNVinterval
sameCNV = function(cR) {
  if ( dim(cR)[1] < 2 ) return(1)

  #find the prior from the gap lengths
  first = 1:(nrow(cR)-1)
  second = 2:nrow(cR)
  dx = abs(cR$x1[second] - cR$x2[first])
  prior = exp(-dx*CNVregionsPerBP())

  #find the probabilities of getting measure values if the SNP frequencies are equal
  f = (cR$var[first] + cR$var[second])/(cR$cov[first] + cR$cov[second])
  p1 = dbinom(cR$var[first], cR$cov[first], f)
  p2 = dbinom(cR$var[second], cR$cov[second], f)
  
  #get posterior using the dx-based prior above. alternative hypothesis is a flat distribution on the frequency.
  postF1 = prior*p1/(prior*p1 + (1-prior)/(1+cR$cov[first]))
  postF2 = prior*p2/(prior*p2 + (1-prior)/(1+cR$cov[second]))
  postF = sqrt(postF1*postF2)
  noFreq = cR$cov[first] == 0 | cR$cov[second] == 0

  #find the probability densities (P*dM) of getting measure values if the regions have the same fold change
  meanM = (cR$M[first]/cR$width[first]^2 + cR$M[second]/cR$width[second]^2)/(1/cR$width[first]^2 + 1/cR$width[second]^2)
  MP1 = pt(-abs(cR$M[first] - meanM)/cR$width[first], df = cR$df[first])
  MP2 = pt(-abs(cR$M[second] - meanM)/cR$width[second], df = cR$df[second])

  #get probability densities for alternative hypothesis: a fold change of at least 1.3 in either direction.
  meanMAup = (cR$M[first]/cR$width[first]^2 + (cR$M[second]+log2(1.3))/cR$width[second]^2)/(1/cR$width[first]^2 + 1/cR$width[second]^2)
  MP1Aup = pt(-abs(cR$M[first] - meanMAup)/cR$width[first], df = cR$df[first])
  MP2Aup = pt(-abs((cR$M[second]+log2(1.3)) - meanMAup)/cR$width[second], df = cR$df[second])
  meanMAdown = (cR$M[first]/cR$width[first]^2 + (cR$M[second]-log2(1.3))/cR$width[second]^2)/(1/cR$width[first]^2 + 1/cR$width[second]^2)
  MP1Adown = pt(-abs(cR$M[first] - meanMAdown)/cR$width[first], df = cR$df[first])
  MP2Adown = pt(-abs((cR$M[second]-log2(1.3)) - meanMAdown)/cR$width[second], df = cR$df[second])

  diffM = abs(cR$M[first] - cR$M[second]) > log2(1.3)
  MP1Aup[diffM] = MP2Aup[diffM] = 1
  MP1Adown[diffM] = MP2Adown[diffM] = 0


  #get posteriors, and take average
  postM1 = prior*MP1/(prior*MP1 + (1-prior)*(MP1Aup+MP1Adown)/2)
  postM2 = prior*MP2/(prior*MP2 + (1-prior)*(MP2Aup+MP2Adown)/2)
  postM = sqrt(postM1*postM2)

  post = sapply(first, function(i) if ( noFreq[i] ) postM[i] else fisherTest(c(postM[i], postF[i]))[2])
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
