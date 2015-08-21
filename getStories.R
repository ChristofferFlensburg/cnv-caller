
#summarises somatic SNVs and CNV calls into subclone evolution over samples
#and determines which subclones are subclones of which other subclones.
getStories = function(variants, normalVariants, cnvs, timeSeries, normals, genome, Rdirectory, plotDirectory, cpus=1, forceRedo=F) {
  setVariantLoss(normalVariants$variants)
  stories = list()
  saveFile = paste0(Rdirectory, '/stories.Rdata')
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading stories.\n')
    load(file=saveFile)
    return(stories)
  }
  if ( length(timeSeries) > 0 ) {
    for ( i in 1:length(timeSeries) ) {
      ts = timeSeries[[i]]
      name = names(timeSeries)[i]
      qs = variants$variants[ts]
      catLog('\nTracking clonal evolution in ', name, '..', sep='')
      
      #select somatic SNPs
      somaticMx = do.call(cbind, lapply(qs, function(q) q$somaticP > 0.95))
      somatic = apply(somaticMx, 1, any)
      somaticQs = lapply(qs, function(q) q[somatic,])

      #switch to effective coverage, to not overestimate the accuracy of high-coverage SNPs
      #not too low though, keep at least 100.
      mC = max(100, maxCov())
      somaticQs = lapply(somaticQs, function(q) {
        d = q$cov
        effectiveCov = round(d*(1 + d/mC)/(1 + d/mC + d^2/mC^2))
        effectiveVar = round(q$var/q$cov*effectiveCov)
        q$var = effectiveVar
        q$cov = effectiveCov
        return(q)
      })
      
      #set clonality of SNPs from frequency and local CNV, return stories
      catLog('SNVs..\n')
      snpStories = findSNPstories(somaticQs, cnvs[ts], normals[ts], filter=T)
      catLog('Keeping', nrow(snpStories), 'SNV stories.\n')
      
      #combine CNV calls over sample into stories
      catLog('Tracking clonal evolution in CNVs..')
      cnvStories = getCNVstories(cnvs[ts], normals[ts], genome, filter=T)
      catLog('Keeping', nrow(cnvStories), 'CNV stories.\n')

      #merge SNV and CNA stories.
      catLog('merge events into clones..')
      allStories = data.frame(row.names=c('germline', rownames(snpStories), rownames(cnvStories)), stringsAsFactors=F)
      allStories$x1 = c(NA, snpStories$x1, cnvStories$x1)
      allStories$x2 = c(NA, snpStories$x2, cnvStories$x2)
      allStories$call = c('germline', as.character(snpStories$call), as.character(cnvStories$call))
      allStories$stories = as.matrix(rbind(matrix(rep(1, length(qs)), nrow=1), snpStories$stories, cnvStories$stories))
      allStories$errors = as.matrix(rbind(matrix(rep(0, length(qs)), nrow=1), snpStories$errors, cnvStories$errors))
      rownames(allStories$stories) = rownames(allStories$errors) = rownames(allStories)
      
      #clusters stories into clones
      clusteredStories = storiesToCloneStories(allStories, variants$SNPs, plotProgress=F, plot=F, cpus=cpus)
      germlineCluster = which(apply(clusteredStories$cloneStories$stories + 1e-3 > 1, 1, all) &
        apply(clusteredStories$cloneStories$errors-1e-5 < 0, 1, all))
      rownames(clusteredStories$cloneStories)[germlineCluster] = clusteredStories$cloneStories$call[germlineCluster] = 'germline'
      rownames(clusteredStories$cloneStories$stories)[germlineCluster] = rownames(clusteredStories$cloneStories$errors)[germlineCluster] = 'germline'

      #add in previously filtered SNV and CNA stories if they fit with the found clones
      #accept somatic SNVs down to 0.5 somaticP for these
      somaticMx = do.call(cbind, lapply(qs, function(q) q$somaticP > 0.5))
      somatic = apply(somaticMx, 1, any)
      somaticQs = lapply(qs, function(q) q[somatic,])
      filteredSnpStories = findSNPstories(somaticQs, cnvs[ts], normals[ts], filter=F)
      filteredSnpStories = filteredSnpStories[!(rownames(filteredSnpStories) %in% snpStories),]
      filteredCnvStories = getCNVstories(cnvs[ts], normals[ts], genome, filter=F)
      filteredCnvStories = filteredCnvStories[!(rownames(filteredCnvStories) %in% cnvStories),]
      
      filteredStories = combineStories(filteredSnpStories, filteredCnvStories)
      filteredStories = filteredStories[!(rownames(filteredStories) %in% rownames(allStories)),]
      
      if ( nrow(filteredStories) > 0 )
        clusteredStories = mergeStories(clusteredStories, filteredStories)
      
      if ( any(!(unlist(clusteredStories$storyList) %in% rownames(allStories))) ) {
        addedStories = unlist(clusteredStories$storyList)[!(unlist(clusteredStories$storyList) %in% rownames(allStories))]
        allStories = combineStories(allStories, filteredStories[addedStories,])
      }
      
      #pick out the variants that behave like germline, ie present clonaly in all samples
      if ( nrow(clusteredStories$cloneStories) > 1 ) {
        germlineVariants = clusteredStories$storyList[[which(rownames(clusteredStories$cloneStories) == 'germline')]]
        germlineVariants = germlineVariants[germlineVariants != 'germline']
      }
      else germlineVariants = c()

      #combine the clustered stories into clonal evolution (figure out which is subclone of which)
      catLog('decide subclone structure..')
      cloneTree = findCloneTree(clusteredStories$cloneStories)
      
      stories[[name]] = list('all'=allStories, 'clusters'=clusteredStories, 'cloneTree'=cloneTree, 'germlineVariants'=germlineVariants)
      catLog('done!\n')
    }
  }

  #flag variants that ended up in the germline clone
  catLog('\nMarking SNV that behave like germline SNPs..')
  variants$variants = lapply(variants$variants, function(q) {
    q$germline = rep(NA, nrow(q))
    return(q)
  })
  if ( length(timeSeries) > 0 ) {
    for ( ind in names(stories) ) {
      story = stories[[ind]]
      if ( length(story$germlineVariants) > 0 ) {
        SNVs = story$germlineVariants[grep('^[0-9]', story$germlineVariants)]
        if ( length(SNVs) > 0 ) {
          samples = timeSeries[[ind]]
          variants$variants[samples] = lapply(variants$variants[samples], function(q) {
            q$germline = rownames(q) %in% SNVs
            return(q)
          })
        }
      }
    }
  }
  allVariants = list('variants'=variants, 'normalVariants'=normalVariants)
  allVariantSaveFile = paste0(Rdirectory, '/allVariants.Rdata')
  save('allVariants', file=allVariantSaveFile)
  catLog('done!\n')


  #return clustered and raw stories
  stories = list('stories'=stories, 'variants'=variants, 'normalVariants'=normalVariants)
  catLog('Saving stories..')
  save('stories', file=saveFile)
  catLog('done!\n\n\n')
  return(stories)
}

findSNPstories = function(somaticQs, cnvs, normal, filter=T) {
  if ( nrow(somaticQs[[1]]) == 0 ) return(data.frame(x1=integer(), x2=integer(), call=character(), stories=numeric(), errors=numeric(), stringsAsFactors=F))
  if ( filter ) {
    cov10 = rowMeans(do.call(cbind, lapply(somaticQs, function(q) q$cov))) >= 10
    somaticQs = lapply(somaticQs, function(q) q[cov10,])
  }
  somaticQs = findSNPclonalities(somaticQs, cnvs)
  clonality = matrix(sapply(somaticQs, function(q) q$clonality), ncol=length(somaticQs))
  clonalityError = matrix(sapply(somaticQs, function(q) q$clonalityError), ncol=length(somaticQs))
  ret = data.frame(x1=somaticQs[[1]]$x, x2=somaticQs[[1]]$x, call=rownames(somaticQs[[1]]), row.names=rownames(somaticQs[[1]]), stringsAsFactors=F)
  ret$stories = clonality
  ret$errors = clonalityError

  colnames(ret$stories) = colnames(ret$errors) = names(somaticQs)
  if ( !filter ) return(ret)
  
  allSmall = rowSums(is.na(ret$errors) | ret$stories - ret$errors*2 < 0 | ret$errors > 0.2) == ncol(ret$stories)
  ret = ret[!allSmall,,drop=F]
  uncertain = rowMeans(ret$errors) > 0.2
  ret = ret[!uncertain,,drop=F]
  indel = grepl('[-\\+]', rownames(ret))
  ret = ret[!indel,,drop=F]
  if ( any(normal) ) {
    presentInNormal = ret$stories[,normal] > ret$errors[,normal] | ret$stories[,normal] > 0.2
    if ( class(presentInNormal) == 'matrix' ) presentInNormal = apply(presentInNormal, 1, any)
    ret = ret[!presentInNormal,,drop=F]
    catLog('Filtered ', sum(allSmall), ' small, ', sum(uncertain), ' uncertain, ', sum(indel), ' indel and ', sum(presentInNormal), ' present in normal stories.\n', sep='')
  }
  else catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain and ', sum(indel), ' indel stories.\n', sep='')

  return(ret)
}

findLocalCNV = function(qs, cnvs) {
  x = qs[[1]]$x
  for ( i in 1:length(cnvs) ) {
    call = sapply(x, function(X) c(cnvs[[i]]$clusters$call[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 'AB')[1])
    clonality = sapply(x, function(X) c(cnvs[[i]]$clusters$clonality[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 1)[1])
    clonalityError = sapply(x, function(X) c(cnvs[[i]]$clusters$clonalityError[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 0)[1])
    if ( any(grepl('CL', call)) ) {
      call[grepl('CL', call)] = 'AB'
      clonality[grepl('CL', call)] = 1 - clonality[grepl('CL', call)]
    }
    qs[[i]]$CNV = call
    qs[[i]]$CNVclonality = clonality
    qs[[i]]$CNVclonalityError = clonalityError
  }
  return(qs)
}

#takes a quality variant object, that has gone through findLocalCNV, and adds a clonality and clonalityError column
findSNPclonalities = function(somaticQs, cnvs) {
  somaticQs = findLocalCNV(somaticQs, cnvs)
  somaticQs = lapply(somaticQs, function(q) {
    nA = nchar(q$CNV) - nchar(gsub('A', '', q$CNV))
    nB = nchar(q$CNV) - nchar(gsub('B', '', q$CNV))

    noCNV = q$CNV %in% c('AB', 'AB?', 'AB??')

    #calculate expected frequencies range if the SNP is on the A or B allele of the CNV, or in the AB background.
    cloneMin = ifelse(noCNV, 0, pmax(0, q$CNVclonality-q$CNVclonalityError*1.5))
    cloneMax = ifelse(noCNV, 0, pmin(1, q$CNVclonality+q$CNVclonalityError*1.5))
    fAmin = nA*cloneMin/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fAmax = nA*cloneMax/((nA+nB)*cloneMax + 2*(1-cloneMax))
    fBmin = nB*cloneMin/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fBmax = nB*cloneMax/((nA+nB)*cloneMax + 2*(1-cloneMax))
    fNmax = (1-cloneMin)/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fNmin = (1-cloneMax)/((nA+nB)*cloneMax + 2*(1-cloneMax))

    #calculate p-values for the SNV living on the different alleles.
    f = refUnbias(q$var/q$cov)
    f[q$cov==0] = 0
    closestFA = ifelse(f < fAmin, fAmin, ifelse(f > fAmax, fAmax, f))
    closestFB = ifelse(f < fBmin, fBmin, ifelse(f > fBmax, fBmax, f))
    closestFN = ifelse(f < fNmax, f, fNmax)
    pA = pBinom(q$cov, q$var, closestFA)
    pB = pBinom(q$cov, q$var, closestFB)
    pN = pBinom(q$cov, q$var, closestFN)

    clonalityA = pmin(1, ifelse(nA==0 & f==0, 0, 2*f/(nA+f*(2-nA-nB))))
    clonalityB = pmin(1, ifelse(nB==0 & f==0, 0, 2*f/(nB+f*(2-nA-nB))))
    clonalityN = pmin(1, 1 - (1 - f*2)/(f*(nA+nB-2)+1))

    #frequency error estimate: add up the poissonian width with the RIB as independent normal error sources.
    fErr = sqrt(frequencyError(q$var, q$cov, p0=0.05)^2 + q$RIB^2)
    fErr[q$cov==0] = 1
    #propagate to clonality
    clonalityHighA = 2/(2+nA/(f+fErr)-nA-nB)
    clonalityHighB = 2/(2+nB/(f+fErr)-nA-nB)
    clonalityHighN = 1 - (1 - (f+fErr)*2)/(pmin(1, f+fErr)*(nA+nB-2)+1)
    #Not allowing the lower limit to go into negatives reduces the error estimate of low frequencies
    #and is associated with a prior bias towards clonality 0 for low frequencies.
    clonalityLowA = pmax(0, 2/(2+nA/(f-fErr)-nA-nB))
    clonalityLowB = pmax(0, 2/(2+nB/(f-fErr)-nA-nB))
    clonalityLowN = pmax(0, 1 - (1 - (f-fErr)*2)/((f-fErr)*(nA+nB-2)+1))
      
    clonalityErrorA = abs(clonalityHighA - clonalityLowA)/2
    clonalityErrorB = abs(clonalityHighB - clonalityLowB)/2
    clonalityErrorN = abs(clonalityHighN - clonalityLowN)/2
    #clonalityError[f==0] = fErr[f==0]

    #if the B allele is lost, and the SNV is clonal in the N-cells, assume that the SNV
    #was present in the cells with the lost B-allele as well, giving a clonality of 1.
    lostSNV = nB==0 & abs(clonalityN+q$CNVclonality - 1) < q$CNVclonalityError
    if ( any(lostSNV) ) clonalityErrorN[lostSNV] = 1

    consistentA = pA > 0.05 & clonalityLowA <= cloneMax
    consistentB = pB > 0.05 & clonalityLowB <= cloneMax
    consistentN = pN > 0.05 & clonalityLowN <= 1-cloneMin
    homeAllele =
      ifelse(q$CNV %in% c('AB', 'AB?', 'AB??', 'CL'),
             'N',
             ifelse(consistentA & (!consistentB | clonalityB <= clonalityA) & (!consistentN | clonalityN <= clonalityA),
                    'A',
                    ifelse(consistentB & (!consistentN | clonalityN <= clonalityB),
                           'B',
                           'N')))
    clonality = ifelse(homeAllele == 'A', clonalityA, ifelse(homeAllele == 'B', clonalityB, clonalityN))
    homeAllele[f==0] = 'N'
    clonality[f==0] = 0
    clonalityError = ifelse(homeAllele == 'A', clonalityErrorA, ifelse(homeAllele == 'B', clonalityErrorB, clonalityErrorN))

    #check for calls where several different home alleles are possible, and increase
    #error estimate accordingly
    uncertainCalls = which(consistentA + consistentB + consistentN > 1)
    if ( length(uncertainCalls) > 0 ) {
      clonalityVar = sapply(uncertainCalls, function(i) var(c(clonalityA[i], clonalityB[i], clonalityN[i])[c(consistentA[i], consistentB[i], consistentN[i])]))
      clonalityError[uncertainCalls] = sqrt(clonalityError[uncertainCalls]^2 + clonalityVar)
    }

    #check if the snv can be a germline SNP, ie present both on the N and one of the A or B alleles.
    #if so, make sure that the clonality error reaches 100%.
    SNPfA = (q$CNVclonality*nA + 1-q$CNVclonality)/(2*(1-q$CNVclonality) + q$CNVclonality*(nA+nB))
    SNPfB = (q$CNVclonality*nB + 1-q$CNVclonality)/(2*(1-q$CNVclonality) + q$CNVclonality*(nA+nB))
    canBeSNPA = pBinom(q$cov, q$var, SNPfA)
    canBeSNPB = pBinom(q$cov, q$var, SNPfB)
    canBeSNP = canBeSNPA > 0.1 | canBeSNPB > 0.1
    clonalityError[canBeSNP] = pmax(clonalityError[canBeSNP], (1-clonality)[canBeSNP])


    q$homeAllele = homeAllele
    q$clonality = noneg(clonality)
    q$clonalityError = abs(clonalityError)
    return(q)
    })
  return(somaticQs)
}

frequencyError = function(var, cov, p0=0.15, reportBothEnds=F) {
  covOri = cov
  varOri = var

  converged = cov == 0
  bsolutions = ifelse(cov==0, 1, 0)
  var = var[!converged]
  cov = cov[!converged]
  b = var - round(sqrt(var))
  p = ifelse(b < 0, 0, 1-pbinom(var-1, cov, noneg(b)/cov))
  while ( any(!converged) ) {
    tooLow = p < p0
    bnew = b + sign(p0-p)
    pnew = ifelse(bnew < 0, 0, 1-pbinom(var-1, cov, noneg(bnew)/cov))
    solved = sign(p-p0) != sign(pnew-p0) 
    #handle converged cases
    if ( any(solved) ) {
      bsolutions[which(!converged)[solved]] = ((b*(pnew-p0)+bnew*(p0-p))/(pnew-p))[solved]
      converged[which(!converged)[solved]] = rep(T, sum(solved))
    }
    #update not converged cases, removing solved values
    p = pnew[!solved]
    b = bnew[!solved]
    cov = cov[!solved]
    var = var[!solved]
  }

  var = varOri
  cov = covOri
  converged = cov == 0
  tsolutions = ifelse(cov==0, 1, 0)
  var = var[!converged]
  cov = cov[!converged]
  t = var + round(sqrt(var))
  p = ifelse(t > cov, 0, pbinom(var, cov, pmin(cov, t)/cov))
  while ( any(!converged) ) {
    tooLow = p < p0
    tnew = t + sign(p-p0)
    pnew = ifelse(tnew > cov, 0, pbinom(var, cov, pmin(cov, tnew)/cov))
    solved = sign(p-p0) != sign(pnew-p0) 
    #handle converged cases
    if ( any(solved) ) {
      tsolutions[which(!converged)[solved]] = ((t*(pnew-p0)+tnew*(p0-p))/(pnew-p))[solved]
      converged[which(!converged)[solved]] = rep(T, sum(solved))
    }
    #update not converged cases, removing solved values
    p = pnew[!solved]
    t = tnew[!solved]
    cov = cov[!solved]
    var = var[!solved]
  }

  if ( reportBothEnds ) return(cbind(pmax(0,bsolutions/covOri), pmin(1,tsolutions/covOri)))
  error = pmax(tsolutions-varOri, varOri-bsolutions)/covOri
  return(error)
}

getCNVstories = function(cnvs, normal, genome, filter=T) {
  regions = splitRegions(cnvs)
  catLog(nrow(regions), ' regions..', sep='')
  events = splitEvents(cnvs, regions)
  catLog(nrow(regions), ' events.\nExtracting stories, comparing SNP shift directions..', sep='')
  stories = cnvsToStories(cnvs, events, normal, genome, filter=filter)
  return(stories)
}

splitRegions = function(cnvs) {
  regions = data.frame('x1' = c(), 'x2' = c(), stringsAsFactors=F)
  cnvX1 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x1)))
  cnvX2 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x2)))
  x1 = x2 = -Inf
  while (x1 < max(cnvX1)) {
    x1 = min(cnvX1[cnvX1 >= x2])
    x2 = min(cnvX2[cnvX2 > x1])
    regions = rbind(regions, data.frame('x1'=x1, 'x2'=x2, stringsAsFactors=F))
  }
  return(regions)
}

splitEvents = function(cnvs, regions) {
  events = data.frame('x1' = c(), 'x2' = c(), 'call' = c(), stringsAsFactors=F)
  cnvX1 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x1))
  cnvX2 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x2))
  cnvCalls = unlist(sapply(cnvs, function(cnv) cnv$cluster$call))
  for ( i in 1:nrow(regions) ) {
    calls = unique(cnvCalls[cnvX1 <= regions$x1[i] & cnvX2 >= regions$x2[i]])
    calls = calls[!(calls %in% c('AB', 'AB?', 'AB??'))]
    for ( call in calls ) events = rbind(events, data.frame('x1' = regions$x1[i], 'x2' = regions$x2[i], 'call' = call, stringsAsFactors=F))
  }
  events = events[!grepl('\\?', events$call),]
  return(events)
}

cnvsToStories = function(cnvs, events, normal, genome, filter=T) {
  stories = errors = data.frame(stringsAsFactors=F)
  if ( nrow(events) > 0 ) {
    i=1
    while ( i <= nrow(events) ) {
      call = as.character(events$call[i])
      clonalities = extractClonalities(cnvs, events[i,])
      stories = rbind(stories, noneg(clonalities$clonality))
      errors = rbind(errors, clonalities$clonalityError)
      if ( any(clonalities$clonality < 0) & any(clonalities$clonality > 0) ) {
        negativeClonalities = clonalities
        negativeClonalities$clonality = noneg(-negativeClonalities$clonality)
        negEvent = events[i,]
        negEvent$call = reverseCall(negEvent$call)
        events = rbind(events, negEvent)[order(c(1:nrow(events), i+0.5)),]
        stories = rbind(stories, noneg(negativeClonalities$clonality))
        errors = rbind(errors, negativeClonalities$clonalityError)
        i = i + 1
      }
      i = i + 1
    }
  }
  else {
    events$stories=matrix(,nrow=0, ncol=length(cnvs))
    events$errors=matrix(,nrow=0,ncol=length(cnvs))
    return(events)
  }
  colnames(errors) = colnames(stories) = names(cnvs)
  ret = events
  rownames(ret) = make.names(paste0('chr', xToChr(ret$x1, genome), '-', events$call), unique=T)
  ret$stories = as.matrix(stories)
  ret$errors = as.matrix(errors)
  ret$errors[is.na(ret$errors)] = 1

  if ( !filter ) return(ret)
  
  allSmall = rowSums(ret$stories - ret$errors*1.5 < 0.15 | ret$errors > 0.2) == ncol(ret$stories)
  ret = ret[!allSmall,,drop=F]
  uncertain = rowMeans(ret$errors) > 0.2
  ret = ret[!uncertain,,drop=F]
  notSignificant = rowSums(abs(ret$stories) < ret$errors*1.5 | is.na(ret$errors)) >= ncol(ret$stories)
  ret = ret[!notSignificant,,drop=F]
  smallRegion = ret$x2 - ret$x1 < 5e6
  ret = ret[!smallRegion,,drop=F]
  if ( any(normal) ) {
    presentInNormal = ret$stories[,normal,drop=F] > ret$errors[,normal,drop=F] | ret$stories[,normal,drop=F] > 0.2
    presentInNormal = apply(presentInNormal, 1, any)
    ret = ret[!presentInNormal,,drop=F]
    catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain, ', sum(smallRegion), ' small region and ', sum(presentInNormal), ' present in normal stories.\n', sep='')
  }
  else catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain and ', sum(smallRegion), ' small region stories.\n', sep='')
  falseSNPcalls = ret$call == 'AA' & (ret$x2 - ret$x1 < 2e6 | rowMeans(ret$errors) > 0.1 ) 
  ret = ret[!falseSNPcalls,,drop=F]
  catLog('Filtered ', sum(falseSNPcalls) , ' stories that are likely based on false SNPs.\n', sep='')
  return(ret)
}

extractClonalities = function(cnvs, event) {
  nSample = length(cnvs)
  subCR = do.call(rbind, lapply(cnvs, function(cs) mergeToOneRegion(cs$CR[cs$CR$x2 >= event$x1 & cs$CR$x1 <= event$x2,,drop=F], cs$eFreqs)))
  subFreqs = lapply(cnvs, function(cs) cs$eFreqs[cs$eFreqs$x >= event$x1 & cs$eFreqs$x <= event$x2,])
  subFreqsMirror = lapply(cnvs, function(cs) {
    ret=cs$eFreqs[cs$eFreqs$x >= event$x1 & cs$eFreqs$x <= event$x2,]
    ret$var = mirrorDown(ret$var, ret$cov)
    return(ret)
    })
  fM = callTofM(as.character(event$call))
  ret = as.data.frame(do.call(rbind, lapply(1:nSample, function(sample)
    unlist(isCNV(cluster=subCR[sample,], efs=subFreqsMirror[[sample]], M=fM[2], f=fM[1],
                 prior=callPrior(as.character(event$call)))))), stringsAsFactors=F)
  direction = rep(0, nSample)

  if ( fM[1] != 0.5 ) {
    freqX = unique(unlist(lapply(subFreqs, function(freq) freq$x)))
    if ( length(freqX) > 0 ) {
      directionScores = do.call(cbind, lapply(1:nSample, function(i) freqToDirectionProb(subFreqs[[i]], freqX, fM[1], ret$clonality[i])))
      strongestSample = which(colSums(directionScores^2) == max(colSums(directionScores^2)))
      direction = sapply(1:nSample, function(i) scalarNorm(directionScores[,strongestSample], directionScores[,i]))
    }
  }

  ret$clonality[ret$sigma > 5] = 0
  ret$clonality[ret$sigma > 5] = 0
  ret$clonalityError[ret$sigma > 5] = sqrt(ret$clonalityError[ret$sigma > 5]^2 + 1/ret$sigma[ret$sigma > 5]^2)
  ret$clonality[!is.na(direction) & direction < 0] = -ret$clonality[!is.na(direction) & direction < 0]
  if ( !any(ret$clonality > 0) & any(ret$clonality < 0) ) ret$clonality = - ret$clonality
  
  return(ret)
}

mergeToOneRegion = function(cR, eFreqs) {
  cR$x2[1] = max(cR$x2)
  cR$var[1] = sum(cR$var)
  cR$cov[1] = sum(cR$cov)
  cR$M[1] = sum(cR$M/cR$width^2)/sum(1/cR$width^2)
  cR$width[1] = 1/sqrt(sum(1/cR$width^2))
  cR = cR[1,]
  cf = correctedFrequency(cR, eFreqs)
  cR$f = cf[,'f']
  cR$ferr = cf[,'ferr']
  return(cR)
}

freqToDirectionProb = function(freq, freqX, f, clonality) {
  if ( f == 0.5 ) return(rep(0, length(freqX)))
  
  upProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    return(pBinom(freq$cov[j], freq$var[j], 0.5 + (0.5-f)*clonality))
  }
  downProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    return(pBinom(freq$cov[j], freq$var[j], 0.5 - (0.5-f)*clonality))
  }
  nullProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    return(pBinom(freq$cov[j], freq$var[j], 0.5))
  }
  u = 1e-5 + sapply(freqX, function(x) upProb(x))
  d = 1e-5 + sapply(freqX, function(x) downProb(x))
  n = 1e-5 + sapply(freqX, function(x) nullProb(x))
  return((1-n)*log(u/d))
}

#takes a dataframe of stories and groups them into subclone stories. returns a data frame of the subclone stories
#and a list of dataframes for the individual stories in each subclone.
storiesToCloneStories = function(stories, SNPs, storyList=as.list(rownames(stories)),
  minDistance=-qnorm(0.01), plotProgress=F, plot=F, cpus=1) {
  if ( length(storyList) < 2 )
    return(list('cloneStories'=stories, 'storyList'=storyList))

  #if a lot of mutations, merge them in batches, as the algorithm scales as O(N^2)
  #group less in this first pass (minDistance*0.5), and then group as specified last round.
  batchSize = 1000
  while ( length(storyList) > batchSize ) {
    catLog(length(storyList), ' stories. Merge first batch separately...', sep='')
    first300 = storiesToCloneStories(stories=stories, SNPs=SNPs, storyList=storyList[1:batchSize],
      minDistance=minDistance*0.5, plotProgress=plotProgress, plot=plot, cpu=cpus)
    storyList = c(first300$storyList, storyList[(batchSize+1):length(storyList)])
  }

  distance = do.call(rbind, mclapply(1:length(storyList), function(i) c(sapply(1:i, function(j) pairScore(stories, storyList[[i]], storyList[[j]])), rep(0, length(storyList)-i)), mc.cores=cpus))
  distance = distance + t(distance)
  distance = distance + (minDistance+1)*as.numeric(row(distance) == col(distance))
  
  while ( any(distance < minDistance) ) {
    merge = which(distance == min(distance), arr.ind=TRUE)[1,]
    if ( plotProgress ) plotStories(stories[c(storyList[[merge[1]]], storyList[[merge[2]]]),],
                                   col=c(rep('blue', length(storyList[[merge[1]]])), rep('red', length(storyList[[merge[2]]]))), SNPs, main=paste0('Distance = ', min(distance)))
    
    storyList[[merge[1]]] = c(storyList[[merge[1]]], storyList[[merge[2]]])
    distance[merge[1],] = distance[,merge[1]] = sapply(1:length(storyList), function(j)
                                      pairScore(stories, storyList[[merge[1]]], storyList[[j]]) + (minDistance+1)*as.numeric(j==merge[1]))
    distance = distance[-merge[2], -merge[2], drop=F]
    storyList = storyList[-merge[2]]
  }
  
  st = do.call(rbind, lapply(storyList, function(rows) {
    err = stories$errors[rows,,drop=F]
    if ( any(err<=0) ) err[err<=0] = rep(min(c(1,err[err>0]))/1e6, sum(err<=0))
    st = stories$stories[rows,,drop=F]
    w = t(1/t(err^2)/colsums(1/err^2))
    mean = colsums(st*w)
    ret = matrix(mean, nrow=1)
  }))
  err = do.call(rbind, lapply(storyList, function(rows) {
    err = stories$errors[rows,,drop=F]
    err = 1/sqrt(colsums(1/err^2))
    ret = matrix(err, nrow=1)
  }))
  rownames(err) = rownames(st) = 1:length(storyList)
  colnames(err) = colnames(st) = colnames(stories$stories)
  clusters = data.frame('call'=rep('clone', length(storyList)), 'x1'=rep(NA, length(storyList)), 'x2'=rep(NA, length(storyList)), row.names=1:length(storyList), stringsAsFactors=F)
  clusters$stories = st
  clusters$errors = err

  if ( plot ) plotStories(clusters, SNPs)
  names(storyList) = rownames(clusters)
  return(list('cloneStories'=clusters, 'storyList'=storyList))
}

#The metric on stories, used for clustering similar stories into subclones.
pairScore = function(stories, is, js) {
  if ( length(is) == 1 ) rms1 = noneg(0.5 - mean(stories$errors[is,]))
  else {
    err1 = stories$errors[is,]
    if ( any(err1<=0) ) err1[err1<=0] = rep(min(c(1,err1[err1>0]))/1e6, sum(err1<=0))
    st1 = stories$stories[is,]
    w1 = t(1/t(err1^2)/colsums(1/err1^2))
    mean1 = colSums(st1*w1)
    sigma1 = abs(t(t(st1)-mean1))/err1
    rms1 = max(sqrt(colmeans(sigma1^2)))
  }
  if ( length(js) == 1 ) rms2 = noneg(0.5 - mean(stories$errors[js,]))
  else {
    err2 = stories$errors[js,]
    if ( any(err2<=0) ) err2[err2<=0] = rep(min(c(1,err2[err2>0]))/1e6, sum(err2<=0))
    st2 = stories$stories[js,]
    w2 = t(1/t(err2^2)/colsums(1/err2^2))
    mean2 = colSums(st2*w2)
    sigma2 = abs(t(t(st2)-mean2))/err2
    rms2 = max(sqrt(mean(sigma2^2)))
  }
  unpairedRms = sqrt(rms1^2+rms2^2)

  err = stories$errors[c(is,js),]
  if ( any(err<=0) ) err[err<=0] = rep(min(c(1,err[err>0]))/1e6, sum(err<=0))
  if ( any(is.infinite(err)) ) err[is.infinite(err)] = rep(1e6, sum(is.infinite(err)))
  st = stories$stories[c(is,js),]
  w = t(1/t(err^2)/colsums(1/err^2))
  mean = colSums(st*w)
  sigma = abs(t(t(st)-mean))/err
  #rms = max(sqrt(colmeans(sigma^2)))

  p = min(apply(pnorm(-sigma)*2, 2, function(ps) min(p.adjust(ps, method='fdr'))))
  return(-qnorm(p/2))
  
  #return(rms + noneg(rms - unpairedRms) )
}


#helper function that decided which subclones are subclones of each other.
findCloneTree = function(cloneStories) {
  #add the purity as a story (this may or may not already be present)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(cloneStories$stories[,i]))
  #check which clones are subclones of which others, allowing an error bar from each clone.
  subcloneMx = as.matrix(apply(abs(cloneStories$stories) + cloneStories$errors*1.5, 1,
    function(story1) apply(abs(cloneStories$stories) - cloneStories$errors*1.5, 1,
                           function(story2) all(story1 > story2))))
  greaterSum = as.matrix(apply(abs(cloneStories$stories), 1,
    function(story1) apply(abs(cloneStories$stories), 1,
                           function(story2) sum(story1) > sum(story2))))
  subcloneMx = subcloneMx & greaterSum
  subcloneMx = enforceTransitive(subcloneMx)
  subcloneStories = cloneStories$stories
  rownames(subcloneMx) = colnames(subcloneMx) = rownames(subcloneStories) = rownames(cloneStories)
  
  cloneTree = list('all'=findChildren(subcloneMx, subcloneStories))

  return(cloneTree)
}
#helper function
findChildren = function(subcloneMx, subcloneStories) {
  #if no more subclones, return empty list
  if ( nrow(subcloneStories) == 0 ) return(list())

  #there will be at least one clone that is not a subclone (aka paranetless, the one with the largest clonality sum)
  is = which(rowSums(subcloneMx) == 0)
  cloneNames = rownames(subcloneStories)[is]
  #score the parentless subclones from the sum of clonalities.
  cloneScores = rowSums(abs(subcloneStories))[is]
  cloneScores = sort(cloneScores, decreasing=T)

  #go through the parentless clones in order of score, each one recurring with its subclones.
  ret = list()
  for ( clone in names(cloneScores) ) {
    subclones = subcloneMx[,clone]
    if ( !any(subclones) ) ret[[clone]] = list()
    else ret[[clone]] = findChildren(subcloneMx[subclones, subclones, drop=F], subcloneStories[subclones,,drop=F])
    subcloneMx = subcloneMx[!subclones, !subclones, drop=F]
    subcloneStories = subcloneStories[!subclones,,drop=F]
  }

  return(ret)
}

#helper function, switching A and B to show that the other allele is affected.
reverseCall = function(call) {
  call = gsub('A', 'a', call)
  call = gsub('B', 'A', call)
  call = gsub('a', 'B', call)
  return(call)
}

#enforces transitivity on a subclone matrix by taking the transitive closure
enforceTransitive = function(mx) {
  for ( col in 1:ncol(mx) ) {
    mx[,col] = mx[,col] | rowsums(mx[,mx[,col],drop=F]) > 0
  }
  return(mx)
}

#add the filtered stories to the clone stories if consistent.
mergeStories = function(clusteredStories, filteredStories) {
  distance = sapply(1:nrow(clusteredStories$cloneStories), function(cloneRow) {
    sapply(1:nrow(filteredStories), function(storyRow) {
      cloneStoryRMS(clusteredStories$cloneStories[cloneRow,], filteredStories[storyRow,])
    })
  })
  #if either as only one row, make it a 1-row or 1-column matrix.
  if ( class(distance) != 'matrix' )
    distance = matrix(distance, nrow = nrow(filteredStories))
  closestDistance = apply(distance, 1, min)
  toMerge = which(closestDistance <= 1)
  bestClone = sapply(toMerge, function(row) which(distance[row,] == closestDistance[row])[1])
  for ( i in 1:length(toMerge) )
    clusteredStories$storyList[[bestClone[i]]] = c(clusteredStories$storyList[[bestClone[i]]], rownames(filteredStories)[toMerge[i]])

  return(clusteredStories)
}

#distance between a story and a clone
cloneStoryRMS = function(clone, story, systematicVariance=0.02) {
  errors = sqrt(clone$errors^2+story$errors^2)+0.02
  sigmas = (clone$stories - story$stories)/errors
  rms = sqrt(mean(sigmas^2))
  return(rms)
}

combineStories = function(stories1, stories2) {
  ret = data.frame(row.names=c(rownames(stories1), rownames(stories2)), stringsAsFactors=F)
  ret$x1 = c(stories1$x1, stories2$x1)
  ret$x2 = c(stories1$x2, stories2$x2)
  ret$call = c(as.character(stories1$call), as.character(stories2$call))
  ret$stories = as.matrix(rbind(stories1$stories, stories2$stories))
  ret$errors = as.matrix(rbind(stories1$errors, stories2$errors))
  rownames(ret$stories) = rownames(ret$errors) = rownames(ret)
  return(ret)
}
