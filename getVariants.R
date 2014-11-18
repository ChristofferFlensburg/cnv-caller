
require(GenomicRanges)

#Takes a set of vcf files and corresponding bam files and capture regions.
#The function filters variants outside of the capture regions, and a few other filters.
#The function the checks the variants in the bamfiles, and flags if the variant seems suspicious.
#Outputs a data frame for each sample with variant and reference counts, as well as some quality information.
getVariants = function(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset, Rdirectory, filterBoring=T, cpus, v, forceRedoSNPs=F, forceRedoVariants=F) {
  SNPsSaveFile = paste0(Rdirectory, '/SNPs.Rdata')
  if ( file.exists(SNPsSaveFile) & !forceRedoSNPs ) {
    catLog('Loading saved SNVs.\n')
    load(file=SNPsSaveFile)
  }
  else {
    catLog('Importing SNV positions..')
    SNPs = try(varScanSNPToMut(vcfFiles, genome=genome))
    if ( class(SNPs) == 'try-error' ) {
      catLog('Error in varScanSNPToMut.\n')
      stop('Error in varScanSNPToMut.')
    }
    if ( any(is.na(as.matrix(SNPs))) ) {
      catLog('NA in SNPs:\n')
      for ( row in which(sapply(1:nrow(SNPs), function(row) any(is.na(SNPs[row,])))) ) catLog(SNPs[row,], '\n')
      stop('NA in SNPs.')
    }
    catLog('done.\n')
    
    catLog('Tagging SNVs with padded capture region..')
    paddedCaptureRegions = captureRegions
    start(paddedCaptureRegions) = start(captureRegions) - 300
    end(paddedCaptureRegions) = end(captureRegions) + 300
    SNPs = inGene(SNPs, paddedCaptureRegions, genome=genome)
    catLog('done.\n')
    
    if ( filterBoring ) {
      boring = SNPs$het == 0 & SNPs$hom == 0
      catLog('Keeping ', sum(!boring), ' out of ', length(boring),
             ' (', round(sum(!boring)/length(!boring), 3)*100, '%) SNVs that are het or hom in at least one sample.\n', sep='')
      SNPs = SNPs[!boring,]
    }

    inCapture = SNP2GRanges(SNPs, genome=genome) %within% paddedCaptureRegions
    catLog('Keeping ', sum(inCapture), ' out of ', length(inCapture),
           ' (', round(sum(inCapture)/length(inCapture), 3)*100, '%) SNVs that are inside capture regions.\n', sep='')
    SNPs = SNPs[inCapture,]

    flag = flagStrandBias(SNPs)
    catLog('Keeping ', sum(!flag), ' out of ', length(flag),
           ' (', round(sum(!flag)/length(flag), 3)*100, '%) SNVs that have consistent strand ratios.\n', sep='')
    SNPs = SNPs[!flag,]
    
    catLog('Matching to dbSNPs.\n')
    SNPs = matchTodbSNPs(SNPs, genome=genome) #adding $db
    
    catLog('Saving SNVs..')
    save(SNPs, file=SNPsSaveFile)
    catLog('done.\n')
  }

  variantsSaveFile = paste0(Rdirectory, '/variants.Rdata')
  if ( file.exists(variantsSaveFile) & !forceRedoVariants ) {
    catLog('Loading saved variants..')
    load(file=variantsSaveFile)
    catLog('done.\n')
  }
  else {
    catLog('Calculating variants:\n')
    gc()
    if ( genome == 'mm10' ) SNPs$chr = paste0('chr', SNPs$chr)
    variants = lapply(bamFiles, function(file) {
      QCsnps(pileups=importQualityScores(SNPs, file, BQoffset, cpus=cpus, v=v)[[1]], SNPs=SNPs, cpus=cpus)})
    names(variants)=names
    variants = shareVariants(variants)
    if ( genome == 'mm10' ) SNPs$chr = gsub('chr', '', SNPs$chr)
    for ( i in 1:length(variants) )  variants[[i]]$db = SNPs[as.character(variants[[i]]$x),]$db
    catLog('Saving variants..')
    save(variants, file=variantsSaveFile)
    catLog('done.\n')
  }
  
  return(list(SNPs=SNPs, variants=variants))
}








#helper function that imports the variants from a vcf file.
varScanSNPToMut = function(files, genome='hg19') {
  #if more than one file, call each file separately and rbind the outputs.
  if ( length(files) > 1 ) {
    catLog('Found', length(files), 'files.', '\n')
    return(do.call(rbind, lapply(files, function(file) varScanSNPToMut(file, genome=genome))))
  }

  catLog('Reading file ', files, '...', sep='')
  raw = read.table(files, fill=T, skip=0, row.names=NULL, header=F, as.is=T)
  catLog('done. Processing data...')
  if ( nrow(raw) == 1 ) return(matrix(1, nrow=0, ncol=22))
  raw = raw[-1,]
  cons = strsplit(as.character(raw$V5), ':')
  strands = strsplit(as.character(raw$V6), ':')
  chrs = gsub('MT', 'M', gsub('chr', '', as.character(raw$V1)))
  ret = data.frame(
    chr = chrs,
    start = as.numeric(as.character(raw$V2)),
    end = as.numeric(as.character(raw$V2)),
    x = chrToX(chrs, as.numeric(as.character(raw$V2)), genome=genome),
    reference = raw$V3,
    variant = raw$V4,
    consensus = as.character(unlist(lapply(cons, function(v) v[1]))),
    reads = as.numeric(unlist(lapply(cons, function(v) v[2]))),
    readsReference = as.numeric(unlist(lapply(cons, function(v) v[3]))),
    readsVariant = as.numeric(unlist(lapply(cons, function(v) v[4]))),
    frequency = as.numeric(gsub('%', '',as.character(unlist(lapply(cons, function(v) v[5])))))/100,
    pValue = as.numeric(unlist(lapply(cons, function(v) v[6]))),
    filter = as.character(unlist(lapply(strands, function(v) v[1]))),
    referencePlus = as.numeric(unlist(lapply(strands, function(v) v[2]))),
    referenceMinus = as.numeric(unlist(lapply(strands, function(v) v[3]))),
    variantPlus = as.numeric(unlist(lapply(strands, function(v) v[4]))),
    variantMinus = as.numeric(unlist(lapply(strands, function(v) v[5]))),
    filterPValue = as.numeric(unlist(lapply(strands, function(v) v[6]))),
    ref = as.numeric(as.character(raw$V7)),
    het = as.numeric(as.character(raw$V8)),
    hom = as.numeric(as.character(raw$V9)),
    nc = as.numeric(as.character(raw$V10)), stringsAsFactors=F)

  rownames(ret) = ret$x
  catLog('done.\n')
  
  catLog('Returning data frame of dimension', dim(ret), '.\n')
  return(ret)
}

#helper function converting from chr+bp coordinates into a single coordinate x that runs over all chromosomes.
chrToX = function(chr, bp, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(humanChrLengths()), 'outside')
  return(prevChrL[gsub('chr', '', chr)] + bp)
}

#Helper function that takes a SNP data frame and returns an GRanges object for the SNPs.
SNP2GRanges = function(SNPs, genome=genome) {
  ir = IRanges(start = SNPs$start, end = SNPs$end, names = rownames(SNPs), width = rep(1, nrow(SNPs)))
  return(GRanges(seqnames = SNPs$chr, ranges = ir, seqlengths = chrLengths(genome)))
}

#helper function that flags variants with suspicious strand distribution.
flagStrandBias = function(SNPs) {
  p = sapply(1:nrow(SNPs), function(row) fisher.test(matrix(
    c(SNPs$referenceMinus[row], SNPs$variantMinus[row],
      SNPs$referencePlus[row], SNPs$variantPlus[row]), nrow=2))$p.value)
  p[is.na(p)] = 1
  fdr = p.adjust(p, method='fdr')
  
  flag = fdr < 0.05
  return(flag)
}
