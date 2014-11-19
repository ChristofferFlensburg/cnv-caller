
#ensures that all variants are present in both samples and normals.
#flags variants with suspicious behaviour in the normals.
#marks the somatic-looking variants in the samples
#ie non-db variants that are not present in the normals and not flagged as suspicious.
matchFlagVariants = function(variants, normalVariants, individuals, normals, Rdirectory, v=v, forceRedoMatchFlag=F) {
  saveFile = paste0(Rdirectory, '/allVariants.Rdata')
  if ( file.exists(saveFile) & !forceRedoMatchFlag ) {
    catLog('Loading final version of combined variants.\n')
    load(file=saveFile)
    return(allVariants)
  }
  variants = matchVariants(variants, normalVariants, v=v)
  normalVariants = matchVariants(normalVariants, variants, v=v)

  #Use the normals to flag variants that are noisy in the normals
  variants = flagFromNormals(variants, normalVariants, cpus=cpus, v=v)

  #mark somatic variants
  variants = markSomatics(variants, normalVariants, individuals, normals, cpus=cpus, v=v)
  
  allVariants = list('variants'=variants, 'normalVariants'=normalVariants)
  catLog('Saving final version of combined variants..')
  save('allVariants', file=saveFile)
  catLog('done.\n')
  return(allVariants)
}



#helper function that ensures that the first variants object have all the variants.
matchVariants = function(vs1, vs2, v='') {
  catLog('Matching variants..')      
  vs1Vars = rownames(vs1$variants[[1]])
  vs2Vars = rownames(vs2$variants[[1]])
  newV1 = setdiff(vs2Vars, vs1Vars)
  catLog(length(vs1Vars) , 'and', length(vs2Vars), 'combine to', length(vs1Vars)+length(newV1), 'variants..')
  is = which(vs2Vars %in% newV1)
  newVars = data.frame(
    x = as.numeric(gsub('[AGNTCagtcn+-].*$', '', newV1)),
    reference = vs2$variants[[1]]$reference[is],
    variant = vs2$variants[[1]]$variant[is],
    cov = rep(0, length(is)),
    ref = rep(0, length(is)),
    var = rep(0, length(is)),
    pbq = rep(1, length(is)),
    pmq = rep(1, length(is)),
    psr = rep(1, length(is)),
    RIB = rep(1, length(is)),
    flag = rep('', length(is)),
    db = vs2$variants[[1]]$db[is],
    somaticP = rep(1, length(is)),
    row.names = newV1)

  #make sure all the columns are in place
  for ( col in setdiff(names(newVars), names(vs1$variants[[1]])) ) {
    catLog('adding ', col, '..', sep='')
    vs1$variants = lapply(vs1$variants, function(q) {q[[col]] = NA; return(q)})
  }
  vs1$variants = lapply(vs1$variants, function(q) rbind(q, newVars))
  vs1$variants = lapply(vs1$variants, function(q) q[order(q$x, q$variant),])
  vs1$SNPs = shareSNPs(vs1$SNPs, vs2$SNPs)
  vs1$SNPs = vs1$SNPs[order(vs1$SNPs$x, vs1$SNPs$variant),]
  catLog('done.\n')
  return(vs1)
}

#helper function that flags variants that have suspicious behaviour in the pool of normals.
flagFromNormals = function(variants, normalVariants, cpus=1, v='') {  
  #check normals for recurring noise.
  setVariantLoss(normalVariants$variants, v=v)
  varN = do.call(cbind, lapply(normalVariants$variants, function(q) q$var))
  covN = do.call(cbind, lapply(normalVariants$variants, function(q) q$cov))
  db = normalVariants$variants[[1]]$db
  f = rowSums(varN)/rowSums(covN)
  f[rowSums(covN) == 0] = 0
  fs = varN/covN
  fs[covN == 0] = 0  
  psN = matrix(pBinom(as.integer(covN), as.integer(varN), rep(f, ncol(covN))), ncol=ncol(covN))
  pSameF = apply(psN, 1, fisherTest)[2,]
  non0 = f > 0.05 & rowMeans(varN) > 1
  isLow = fs < 0.1
  isHigh = fs > 0.9
  isHalf = matrix(pBinom(as.integer(covN), as.integer(varN), rep(refBias(0.5), nrow(covN)*ncol(covN))), ncol=ncol(covN)) > 0.01
  consistent = rowSums(isLow | isHigh | isHalf) == ncol(fs)   #are all samples 0, 0.5 or 1?
  normalNoise = (!db & non0 & pSameF > 0.01) | (db & non0 & pSameF > 0.01 & !consistent)
  catLog('Flagged', sum(normalNoise), 'out of', nrow(fs),
         'variants that are recurrently and consistently noisy in normals.\n')
  
  #check if the variants are consistent with the normal noise frequency, and flag Nnc or Nnn.
  variants$variants =
    lapply(variants$variants, function(q) {
        ps = pBinom(q$cov[normalNoise], q$var[normalNoise], f[normalNoise])
        flag = ifelse(ps > 0.01, 'Nnc', 'Nnn')        #normal noise (non)-consistent
        q$flag[normalNoise] = paste0(q$flag[normalNoise], flag)
        return(q)
    })

  #check for non-db SNPs that behave as polymoprhic db SNPs in the normals
  polymorphic = (!db &                                                   #db
                 rowSums(isLow | isHigh | isHalf) == ncol(fs) &          #all consistent
                 rowSums(isLow & !isHalf) > 0 &                          #one strictly ref 
                 rowSums(!isLow & !isHigh & isHalf) > 0 &                #one strictly het
                 pSameF < 0.01)                                          #not same frequency
  
  catLog('Flagged', sum(polymorphic), 'out of', sum(!db),
         ' non-db variants that are consistently polymorphic in normals.\n')

  variants$variants =
    lapply(variants$variants, function(q) {
      ps = pBinom(q$cov[polymorphic], q$var[polymorphic], f[polymorphic])
      q$flag[polymorphic] = paste0(q$flag[polymorphic], 'Pn')  #polymorphic Normal
      return(q)
    })
  
  return(variants)
}

#This helper function assign probabilities that variants are somatics.
# the probabilities that the variants are true somatic variants are added in a column 'somaticP'
markSomatics = function(variants, normalVariants, individuals, normals, cpus=cpus, v=v) {
  names = names(variants$variants)
  #pair up cancer normals
  correspondingNormal = findCorrespondingNormal(names, individuals, normals)
  CNs = which(!is.na(correspondingNormal))
  names(CNs) = names(variants$variants[CNs])

  #calculate and assign somatic p-values for the cancer samples
  #use the newVariant p-values
  for ( name in names(CNs) ) {
    catLog('Marking somatic mutations in ', name, ' using matching normal ', correspondingNormal[name], '..', sep='')
    q1 = variants$variants[[correspondingNormal[name]]]
    q2 = variants$variants[[name]]
    use = !q2$db & q2$flag == ''
    q1 = q1[use,]
    q2 = q2[use,]
    freq1 = q1$var/q1$cov
    freq1[is.na(freq1)] = -0.02
    freq2 = q2$var/q2$cov
    freq2[is.na(freq2)] = -0.02
  
    require(parallel)
    ps = unlist(mclapply(1:length(freq1), function(i) fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value, mc.cores=cpus))
  
    low = (freq1 < 0.2 | (pBinom(q1$cov, q1$var, 0.3) < 0.01 & freq1 < 0.3)) & !(freq1 == 0 & freq2 == 0)
    fdr = p.adjust(ps, method='fdr')
    fdr[low] = p.adjust(ps[low], method='fdr')
    normalOK = which(freq1 < freq2 & low)
    fdr[!normalOK] = 1
    pSampleOK = p.adjust(q2$pbq, method='fdr')*p.adjust(q2$pmq, method='fdr')*p.adjust(q2$psr, method='fdr')
    pZero = p.adjust(dbinom(q2$var, q2$cov, 0.001), method='fdr')   #the base quality cut on 30 correpsonds to 0.001 wrong base calls.
    
    variants$variants[[name]]$somaticP = 0
    variants$variants[[name]]$somaticP[use] = (1-fdr)*pSampleOK*(1-pZero)
    catLog('got roughly ', sum(variants$variants[[name]]$somaticP > 0.5), ' somatic variants.\n', sep='')
  }

  #require consistent and very low frequency between normals.
  varN = do.call(cbind, lapply(normalVariants$variants, function(q) q$var))
  covN = do.call(cbind, lapply(normalVariants$variants, function(q) q$cov))
  f = rowSums(varN)/rowSums(covN)
  fs = varN/covN
  fs[covN == 0] = 0  
  f[rowSums(covN) == 0] = 0
  psN = matrix(pBinom(as.integer(covN), as.integer(varN), rep(f, ncol(covN))), ncol=ncol(covN))
  pSameF = apply(psN, 1, fisherTest)[2,]
  
  #count the ratio of non-db SNPs that have the 'Pn' flag, ie polymorphic normal
  observedPolymorphic = length(grep('Pn',variants$variants[[1]]$flag))
  basePairs = 3e9
  nNormals = length(normalVariants$variants)
  polymorphicFrequency = 1-(1-observedPolymorphic/basePairs)^(1/nNormals)
  
  #calculate somatic p-values for the rest of the samples
  for ( name in names[!(names %in% CNs)] ) {
    catLog('Marking somatic mutations in ', name, ' using pool of normals..', sep='')
    q = variants$variants[[name]]
    #Require non-db SNP, and no flag
    use = !q$db & q$flag == ''
    q = q[use,]
    freq = q$var/q$cov
    freq[is.na(freq)] = -0.02
    normalFreq = f[use]
    
    
    #set p-value from both number of normals (for population-wide frequency uncertainty)
    #and difference in frequency (is the cancer sample really different)
    #should effectively be a cut on how certain it can be from number of normals, even with
    #perfect 0,0,0,0, 0.5 frequencies
    pPolymorphic = 1/(1+nrow(q)/(polymorphicFrequency*basePairs))
    pNormalFreq = pBinom(q$cov, q$var, normalFreq)
    normalOK = normalFreq < 0.01 & normalFreq < freq
    pSampleOK = p.adjust(q$pbq, method='fdr')*p.adjust(q$pmq, method='fdr')*p.adjust(q$psr, method='fdr')
    pZero = p.adjust(dbinom(q$var, q$cov, 0.01), method='fdr')   #the base quality cut on 30 correpsonds to 0.001 wrong base calls.

    somaticP = (1-pPolymorphic)*(1-pNormalFreq)*normalOK*pSampleOK*(1-pZero)
 
    variants$variants[[name]]$somaticP = 0
    variants$variants[[name]]$somaticP[use] = somaticP
    catLog('got roughly ', sum(variants$variants[[name]]$somaticP > 0.5), ' somatic variants.\n', sep='')
  }

  return(variants)
}

#helper function that matches samples with a normal from the same individual, if present.
findCorrespondingNormal = function(names, individuals, normals) {
  ret = sapply(1:length(names), function(name) {
    ind = individuals[name]
    isCancer = !normals[name]
    if ( !isCancer ) return(NA)
    hasNormal = any(normals & individuals == ind)
    if ( !hasNormal ) return(NA)
    return(which(normals & individuals == ind)[1]) #fix this to support replicate normals!
  })
  names(ret) = names
  return(ret)
}

#calculates the p-value from a two-tailed binomial.
pBinom = function(cov, var, f) {
  use = cov > 0 & var >= 0
  p = rep(1, length(cov))
  if ( length(f) == 1 ) f = rep(f, length(cov))
  if ( length(f) != length(cov) | length(var) != length(cov) )
    cat('Length of f', length(f), 'must match length of cov', length(cov), 'or be 1.\n')
  p[use] = pbinom(var[use], cov[use], f[use]) -  dbinom(var[use], cov[use], f[use])/2
  p[use] = 2*pmin(p[use], 1-p[use])
  return(p)
}

#helper function that matches variants between two SNPs objects
shareSNPs = function(SNPs1, SNPs2) {
  newX = setdiff(SNPs2$x, SNPs1$x1)
  newSNPs = SNPs2[as.character(newX),]
  newSNPs$reads = newSNPs$readsReference = newSNPs$readsVariant = newSNPs$frequency = 0
  newSNPs$referencePlus = newSNPs$referenceMinus = newSNPs$variantPlus = newSNPs$variantMinus = 0
  newSNPs$frequency = -0.02
  newSNPs$pValue = newSNPs$filterPValue = 1
  newSNPs$ref = newSNPs$het = newSNPs$hom = 0
  newSNPs$nc = max(SNPs1$ref+SNPs1$het+SNPs1$hom+SNPs1$nc)
  return(rbind(SNPs1, newSNPs))
}

#helper function combining p-values using fishers method.
fisherTest = function(p) {
  p = p[!is.na(p)]
  Xsq = -2*sum(log(p))
  pVal = 1-pchisq(Xsq, df = 2*length(p))
  return(c(Xsq = Xsq, pVal = pVal))
}
