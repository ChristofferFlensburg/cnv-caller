
#summarises somatic SNVs and CNV calls into subclone evolution over samples
#and determines which subclones are subclones of which other subclones.
getStories = function(variants, normalVariants, cnvs, timeSeries, Rdirectory, plotDirectory, cpus=1, forceRedo=F) {
  setVariantLoss(normalVariants$variants)
  stories = list()
  if ( length(timeSeries) == 0 ) return(stories)
  saveFile = paste0(Rdirectory, '/stories.Rdata')
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading stories.\n')
    load(file=saveFile)
    return(stories)
  }
  for ( i in 1:length(timeSeries) ) {
    ts = timeSeries[[i]]
    name = names(timeSeries)[i]
    qs = variants$variants[ts]
    catLog('Tracking clonal evolution in ', name, '..', sep='')
    
    #select somatic SNPs
    somaticMx = do.call(cbind, lapply(qs, function(q) q$somaticP > 0.95))
    somatic = apply(somaticMx, 1, any)
    somaticQs = lapply(qs, function(q) q[somatic,])

    #set clonality of SNPs from frequency and local CNV, return stories
    catLog('SNVs..\n')
    snpStories = findSNPstories(somaticQs, cnvs[ts])
    catLog('Keeping', nrow(snpStories), 'SNV stories.\n')

    #combine CNV calls over sample into stories
    catLog('Tracking clonal evolution in CNVs..\n')
    cnvStories = getCNVstories(cnvs[ts])
    catLog('Keeping', nrow(cnvStories), 'CNV stories.\n')

    #cluster stories for SNPs, CNVs, both.
    catLog('merge events into clones..')
    allStories = data.frame(row.names=c(rownames(snpStories), rownames(cnvStories)))
    allStories$x1 = c(snpStories$x1, cnvStories$x1)
    allStories$x2 = c(snpStories$x2, cnvStories$x2)
    allStories$call = c(as.character(snpStories$call), as.character(cnvStories$call))
    allStories$stories = as.matrix(rbind(snpStories$stories, cnvStories$stories))
    allStories$errors = as.matrix(rbind(snpStories$errors, cnvStories$errors))
    rownames(allStories$stories) = rownames(allStories$errors) = rownames(allStories)
    clusteredStories = storiesToCloneStories(allStories, plotProgress=F, plot=F)
    
    #combine the clustered stories into clonal evolution (figure out which is subclone of which)
    catLog('decide subclone structure..')
    cloneTree = findCloneTree(clusteredStories$cloneStories)
    
    stories[[name]] = list('all'=allStories, 'clusters'=clusteredStories, 'cloneTree'=cloneTree)
    catLog('done!\n')
  }

  #return clustered and raw stories
  catLog('Saving stories.\n')
  save('stories', file=saveFile)
  return(stories)
}

findSNPstories = function(somaticQs, cnvs) {
  cov10 = rowMeans(do.call(cbind, lapply(somaticQs, function(q) q$cov))) >= 10
  somaticQs = lapply(somaticQs, function(q) q[cov10,])
  somaticQs = findSNPclonalities(somaticQs, cnvs)
  clonality = matrix(sapply(somaticQs, function(q) q$clonality), ncol=length(somaticQs))
  clonalityError = matrix(sapply(somaticQs, function(q) q$clonalityError), ncol=length(somaticQs))
  ret = data.frame(x1=somaticQs[[1]]$x, x2=somaticQs[[1]]$x, call=rownames(somaticQs[[1]]), row.names=rownames(somaticQs[[1]]))
  ret$stories = clonality
  ret$errors = clonalityError

  allSmall = rowSums(abs(ret$stories) < ret$errors*2 | is.na(ret$errors) |ret$stories < 0.3) == ncol(ret$stories)
  ret = ret[!allSmall,]
  uncertain = rowMeans(ret$errors) > 0.25
  ret = ret[!uncertain,]
  catLog('Filtered ', sum(allSmall) , ' small and ', sum(uncertain), ' uncertain stories.\n', sep='')

  return(ret)
}

findLocalCNV = function(qs, cnvs) {
  x = qs[[1]]$x
  for ( i in 1:length(cnvs) ) {
    call = sapply(x, function(X) c(cnvs[[i]]$clusters$call[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 'AB')[1])
    clonality = sapply(x, function(X) c(cnvs[[i]]$clusters$clonality[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 1)[1])
    clonalityError = sapply(x, function(X) c(cnvs[[i]]$clusters$clonalityError[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 0)[1])
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
    #calculate expected frequencies range if the SNP is on the A or B allele of the CNV, or in the AB background.
    cloneMin = pmax(0, q$CNVclonality-q$CNVclonalityError)
    cloneMax = pmin(1, q$CNVclonality+q$CNVclonalityError)
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
    closestFN = ifelse(f < fNmin, fNmin, ifelse(f > fNmax, fNmax, f))
    pA = pBinom(q$cov, q$var, closestFA)
    pB = pBinom(q$cov, q$var, closestFB)
    pN = pBinom(q$cov, q$var, closestFN)

    homeAllele = ifelse(q$CNV %in% c('AB', 'AB?', 'AB??', 'CL'), 'N', ifelse(pA > pmax(pB, pN), 'A', ifelse(pB > pN, 'B', 'N')))
    homeAllele[f==0] = 'N'
    clonality = pmin(1, ifelse(homeAllele == 'A', 2*f/(nA+f*(2-nA-nB)),
      ifelse(homeAllele == 'B', 2*f/(nB+f*(2-nA-nB)), 1 - (1 - f*2)/(f*(nA+nB-2)+1))))
    clonality[f==0] = 0

    #frequency error estimate: add up the poissonian width with the RIB as independent normal error sources.
    fErr = sqrt(frequencyError(q$var, q$cov, p0=0.15)^2 + q$RIB^2)
    fErr[q$cov==0] = 1
    #propagate to clonality
    clonalityHigh = ifelse(homeAllele == 'A', 2/(2+nA/(f+fErr)-nA-nB),
      ifelse(homeAllele == 'B', 2/(2+nB/(f+fErr)-nA-nB), 1 - (1 - (f+fErr)*2)/((f+fErr)*(nA+nB-2)+1)))
    #Not allowing the lower limit to go into negatives reduces the error estimate of low frequnecies
    #and is associated with a prior bias towards clonality 0 for low frequencies.
    clonalityLow = pmax(0, ifelse(homeAllele == 'A', 2/(2+nA/(f-fErr)-nA-nB),
      ifelse(homeAllele == 'B', 2/(2+nB/(f-fErr)-nA-nB), 1 - (1 - (f-fErr)*2)/((f-fErr)*(nA+nB-2)+1))))
    clonalityError = (clonalityHigh - clonalityLow)/2
    #clonalityError[f==0] = fErr[f==0]

    q$homeAllele = homeAllele
    q$clonality = clonality
    q$clonalityError = clonalityError
    return(q)
    })
  return(somaticQs)
}

frequencyError = function(var, cov, p0=0.15) {
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

  error = pmax(tsolutions-varOri, varOri-bsolutions)/covOri
  return(error)
}

getCNVstories = function(cnvs) {
  regions = splitRegions(cnvs)
  events = splitEvents(cnvs, regions)
  stories = cnvsToStories(cnvs, events)
  return(stories)
}

splitRegions = function(cnvs) {
  regions = data.frame('x1' = c(), 'x2' = c())
  cnvX1 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x1)))
  cnvX2 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x2)))
  x1 = x2 = -Inf
  while (x1 < max(cnvX1)) {
    x1 = min(cnvX1[cnvX1 >= x2])
    x2 = min(cnvX2[cnvX2 > x1])
    regions = rbind(regions, data.frame('x1'=x1, 'x2'=x2))
  }
  return(regions)
}

splitEvents = function(cnvs, regions) {
  events = data.frame('x1' = c(), 'x2' = c(), 'call' = c())
  cnvX1 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x1))
  cnvX2 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x2))
  cnvCalls = unlist(sapply(cnvs, function(cnv) cnv$cluster$call))
  for ( i in 1:nrow(regions) ) {
    calls = unique(cnvCalls[cnvX1 <= regions$x1[i] & cnvX2 >= regions$x2[i]])
    calls = calls[!(calls %in% c('AB', 'AB?', 'AB??'))]
    for ( call in calls ) events = rbind(events, data.frame('x1' = regions$x1[i], 'x2' = regions$x2[i], 'call' = call))
  }
  events = events[!grepl('\\?', events$call),]
  return(events)
}

cnvsToStories = function(cnvs, events) {
  stories = errors = data.frame()
  if ( nrow(events) > 0 ) {
    for ( i in 1:nrow(events) ) {
      call = as.character(events$call[i])
      clonalities = extractClonalities(cnvs, events[i,])
      stories = rbind(stories, clonalities$clonality)
      errors = rbind(errors, clonalities$clonalityError)
    }
  }
  else {
    events$stories=matrix(,nrow=0, ncol=length(cnvs))
    events$errors=matrix(,nrow=0,ncol=length(cnvs))
    return(events)
  }
  colnames(errors) = colnames(stories) = names(cnvs)
  ret = events
  rownames(ret) = make.names(paste0('chr', xToChr(ret$x1), '-', events$call), unique=T)
  ret$stories = as.matrix(stories)
  ret$errors = as.matrix(errors)
  ret$errors[is.na(ret$errors)] = 1

  allSmall = rowSums(ret$stories < 0.3) == ncol(ret$stories)
  ret = ret[!allSmall,]
  uncertain = rowMeans(ret$errors) > 0.4
  ret = ret[!uncertain,]
  notSignificant = rowSums(abs(ret$stories) < ret$errors*1.5 | is.na(ret$errors)) >= ncol(ret$stories)
  ret = ret[!notSignificant,]
  catLog('Filtered ', sum(allSmall) , ' all small, ', sum(uncertain), ' uncertain and ', sum(notSignificant), ' notsignificant stories.\n', sep='')
  return(ret)
}

extractClonalities = function(cnvs, event) {
  nSample = length(cnvs)
  subCR = do.call(rbind, lapply(cnvs, function(cs) mergeToOneRegion(cs$CR[cs$CR$x1 >= event$x1 & cs$CR$x2 <= event$x2,,drop=F])))
  subFreqs = lapply(cnvs, function(cs) cs$freqs[cs$freqs$x >= event$x1 & cs$freqs$x <= event$x2,])
  subFreqsMirror = lapply(cnvs, function(cs) {
    ret=cs$freqs[cs$freqs$x >= event$x1 & cs$freqs$x <= event$x2,]
    ret$var = mirrorDown(ret$var, ret$cov)
    return(ret)
    })
  fM = callTofM(as.character(event$call))
  ret = as.data.frame(do.call(rbind, lapply(1:nSample, function(sample)
    unlist(isCNV(cluster=subCR[sample,], freqs=subFreqsMirror[[sample]], M=log2(fM[2]), f=fM[1],
                 prior=callPrior(as.character(event$call)))))))
  direction = rep(0, nSample)

  if ( fM[1] != 0.5 ) {
    freqX = unique(unlist(lapply(subFreqs, function(freq) freq$x)))
    if ( length(freqX) > 0 ) {
      directionScores = do.call(cbind, lapply(1:nSample, function(i) freqToDirectionProb(subFreqs[[i]], freqX, fM[1], ret$clonality[i])))
      strongestSample = which(colSums(directionScores^2) == max(colSums(directionScores^2)))
      direction = sapply(1:nSample, function(i) scalarNorm(directionScores[,strongestSample], directionScores[,i]))
    }
  }

  ret$direction = direction

  ret$clonality[ret$sigma > 5] = 0
  ret$clonality[ret$sigma > 5] = 0
  ret$clonalityError[ret$sigma > 5] = sqrt(ret$clonalityError[ret$sigma > 5]^2 + 1/ret$sigma[ret$sigma > 5]^2)
  ret$direction[ret$sigma > 5] = 0
  ret$clonality[!is.na(ret$direction) & ret$direction < 0] = -ret$clonality[!is.na(ret$direction) & ret$direction < 0]

  return(ret)
}

mergeToOneRegion = function(cR) {
  cR$x2[1] = cR$x2[nrow(cR)]
  cR$var[1] = sum(cR$var)
  cR$cov[1] = sum(cR$cov)
  cR$M[1] = sum(cR$M/cR$width^2)/sum(1/cR$width^2)
  cR$width[1] = 1/sqrt(sum(1/cR$width^2))
  cR = cR[1,]
  cR$f = refUnbias(cR$var/cR$cov)
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
  u = sapply(freqX, function(x) upProb(x))
  d = sapply(freqX, function(x) downProb(x))
  n = sapply(freqX, function(x) nullProb(x))
  return((1-n)*log(u/d))
}

#takes a dataframe of stories and groups them into subclone stories. returns a data frame of the subclone stories
#and a list of dataframes for the individual stories in each subclone.
storiesToCloneStories = function(stories, minDistance=2, plotProgress=F, plot=F) {
  storyList = as.list(rownames(stories))
  if ( length(storyList) < 2 ) 
    return(list('cloneStories'=stories, 'storyList'=storyList))

  
  distance = do.call(rbind, lapply(1:length(storyList), function(is) sapply(1:length(storyList), function(js) pairScore(stories, is, js))))
  distance = distance + minDistance*as.numeric(row(distance) == col(distance))
  while ( any(distance < minDistance) ) {
    merge = which(distance == min(distance), arr.ind=TRUE)[1,]
    if ( plotProgress ) plotStories(stories[c(storyList[[merge[1]]], storyList[[merge[2]]]),],
                                   col=c(rep('blue', length(storyList[[merge[1]]])), rep('red', length(storyList[[merge[2]]]))))
    
    storyList[[merge[1]]] = c(storyList[[merge[1]]], storyList[[merge[2]]])
    distance[merge[1],] = distance[,merge[1]] = sapply(1:length(storyList), function(j)
                                      pairScore(stories, storyList[[merge[1]]], storyList[[j]]) + minDistance*as.numeric(j==merge[1]))
    distance = distance[-merge[2], -merge[2]]
    storyList = storyList[-merge[2]]
  }
  
  st = do.call(rbind, lapply(storyList, function(rows) {
    err = stories$errors[rows,,drop=F]
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
  clusters = data.frame('call'=rep('clone', length(storyList)), 'x1'=rep(NA, length(storyList)), 'x2'=rep(NA, length(storyList)), row.names=1:length(storyList))
  clusters$stories = st
  clusters$errors = err

  if ( plot ) plotStories(clusters)
  names(storyList) = rownames(clusters)
  return(list('cloneStories'=clusters, 'storyList'=storyList))
}

#The metric on stories, used for clustering similar stories into subclones.
pairScore = function(stories, is, js) {
  if ( length(is) == 1 ) wms1 = 0.5
  else {
    err1 = stories$errors[is,]
    st1 = stories$stories[is,]
    w1 = t(1/t(err1^2)/colsums(1/err1^2))
    mean1 = colSums(st1*w1)
    sigma1 = abs(t(t(st1)-mean1))/err1
    wms1 = sum(sigma1*w1)/sum(w1)
  }
  if ( length(js) == 1 ) wms2 = 0.5
  else {
    err2 = stories$errors[js,]
    st2 = stories$stories[js,]
    w2 = t(1/t(err2^2)/colsums(1/err2^2))
    mean2 = colSums(st2*w2)
    sigma2 = abs(t(t(st2)-mean2))/err2
    wms2 = sum(sigma2*w2)/sum(w2)
  }
  W1 = sum(1/stories$errors[is,]^2)
  W2 = sum(1/stories$errors[js,]^2)
  unpairedWms = (wms1*W1+wms2*W2)/(W1+W2)

  err = stories$errors[c(is,js),]
  st = stories$stories[c(is,js),]
  w = t(1/t(err^2)/colsums(1/err^2))
  mean = colSums(st*w)
  sigma = abs(t(t(st)-mean))/err
  wms = sum(sigma*w)/sum(w)

  return(wms + noneg(wms - unpairedWms) )
}

#defines a metric between two stories
storyDistance = function(story1, story2) {
  stories = rbind(story1, story2)
  err = stories$errors
  st = stories$stories
  w = t(1/t(err^2)/colsums(1/err^2))
  mean = colSums(st*w)
  sigma = abs(t(t(st)-mean))/err
  rms = sqrt(mean(sigma^2))

  return(rms)
}

#helper function that decided which subclones are subclones of each other.
findCloneTree = function(cloneStories) {
  #add the purity as a story (this may or may not already be present)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(cloneStories$stories[,i]))
  #check which clones are subclones of which others, allowing an error bar from each clone.
  subcloneMx = as.matrix(apply(abs(cloneStories$stories) + cloneStories$errors, 1,
    function(story1) apply(abs(cloneStories$stories) - cloneStories$errors, 1,
                           function(story2) all(story1 > story2))))
  greaterSum = as.matrix(apply(abs(cloneStories$stories), 1,
    function(story1) apply(abs(cloneStories$stories), 1,
                           function(story2) sum(story1) > sum(story2))))
  subcloneMx = subcloneMx & greaterSum
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
