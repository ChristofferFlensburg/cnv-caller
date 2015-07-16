

runVEP = function(variants, dir, cpus=1, forceRedoVEP=F) {
  catLog('VEP-ing..')
  for ( name in names(variants$variants) ) {
    infile = paste0(dir, '/', name, '.txt')
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    wd = getwd()
    setwd(dir)
    if ( !file.exists(VEPfile) | forceRedoVEP ) {
      catLog(name, ': ', sum(variants$variants[[name]]$somaticP > 0), ' variants.\n', sep='')
      if ( sum(variants$variants[[name]]$somaticP > 0) == 0 ) next
      call = paste0('vep -i ', basename(infile), ' -o ', basename(VEPfile), ' --everything --force_overwrite --fork ', cpus)
      cat(call, '\n')
      systemRet = system(call, intern=T)
      if ( !any(grepl('Finished', systemRet)) ) warning('VEP run didnt finish!')
    }
    setwd(wd)
  }
  catLog('done.\n')
  
  catLog('Importing VEP results:\n')
  for ( name in names(variants$variants) ) {
    catLog(name, ': ', sep='')
    if ( sum(variants$variants[[name]]$somaticP > 0) == 0 ) {
      catLog('no somatic variants.\n')
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
      next
    }
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    VEPdata = try(read.table(VEPfile, fill=T), silent=T)
    if ( class(VEPdata) == 'try-error' ) {
      catLog('failed to read VEP file.\n')
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
    }
    else {
      ID = as.numeric(as.character(VEPdata$V1))
      chr = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[1])
      pos = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[2])
      end = gsub('^.*-', '',pos)
      pos = gsub('-.*$', '',pos)
      var = as.character(VEPdata$V3)
      type = as.character(VEPdata$V7)
      sev = sapply(type, typeToSeverity)
      x = chrToX(chr, as.numeric(pos))

      lastCol = strsplit(as.character(VEPdata$V14), ';')
      polyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) ) return(as.numeric(gsub('\\)$', '', gsub('^.*\\(', '', strs[grep('PolyPhen=',strs)[1]]))))
        else return(0.05)
        })
      
      namevar = var
      namevar[nchar(namevar) > 1] = paste0('+', substring(var[nchar(var) > 1], 2))
      namevar[namevar == '-'] = paste0('-', pmax(1, as.numeric(end)[namevar == '-']-as.numeric(pos)[namevar == '-']))
      rowNames = paste0(x, namevar)
      qNames = rep('', length(unique(ID)))
      qNames[ID] = rowNames
      
      #find most severe effect
      mostSev = rep('unkown', length(unique(ID)))
      sevScore = rep(100, length(unique(ID)))
      for ( i in 1:nrow(VEPdata) ) {
        sevI = sev[i] - polyPhen[i]
        IDI = ID[i]
        if ( sevI < sevScore[IDI] ) {
          sevScore[IDI] = sevI - polyPhen[i]
          mostSev[IDI] = severityToType(sevI)
        }
      }
      
      #add effect and severity to q
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
      
      variants$variants[[name]][qNames,]$type = mostSev
      variants$variants[[name]][qNames,]$severity = sevScore
      catLog(length(mostSev), 'VEPed variants.\n')
    }
  }
  catLog('done.\n')
  return(variants)
}

typeToSeverity = function(type) {
  #first change some (older? newer?) notation into the ensemble SO terms
  type = gsub('non_coding_exon_variant', 'non_coding_transcript_exon_variant', type)
  type = gsub('nc_transcript_variant', 'non_coding_transcript_exon_variant', type)

  if ( grepl('transcript_ablation', type) ) return(1)
  if ( grepl('splice_acceptor_variant', type) ) return(2)
  if ( grepl('splice_donor_variant', type) ) return(3)
  if ( grepl('stop_gained', type) ) return(4)
  if ( grepl('frameshift_variant', type) ) return(5)
  if ( grepl('stop_lost', type) ) return(6)
  if ( grepl('initiator_codon_variant', type) ) return(7)
  if ( grepl('transcript_amplification', type) ) return(8)
  if ( grepl('inframe_insertion', type) ) return(9)
  if ( grepl('inframe_deletion', type) ) return(10)
  if ( grepl('missense_variant', type) ) return(11)
  if ( grepl('splice_region_variant', type) ) return(12)
  if ( grepl('incomplete_terminal_codon_variant', type) ) return(13)
  if ( grepl('stop_retained_variant', type) ) return(14)
  if ( grepl('synonymous_variant', type) ) return(15)
  if ( grepl('coding_sequence_variant', type) ) return(16)
  if ( grepl('mature_miRNA_variant', type) ) return(17)
  if ( grepl('5_prime_UTR_variant', type) ) return(18)
  if ( grepl('3_prime_UTR_variant', type) ) return(19)
  if ( grepl('non_coding_transcript_exon_variant', type) ) return(20)
  if ( grepl('intron_variant', type) ) return(21)
  if ( grepl('NMD_transcript_variant', type) ) return(22)
  if ( grepl('non_coding_transcript_variant', type) ) return(23)
  if ( grepl('upstream_gene_variant', type) ) return(24)
  if ( grepl('downstream_gene_variant', type) ) return(25)
  if ( grepl('TFBS_ablation', type) ) return(26)
  if ( grepl('TFBS_amplification', type) ) return(27)
  if ( grepl('TF_binding_site_variant', type) ) return(28)
  if ( grepl('regulatory_region_ablation', type) ) return(29)
  if ( grepl('regulatory_region_amplification', type) ) return(30)
  if ( grepl('regulatory_region_variant', type) ) return(31)
  if ( grepl('feature_elongation', type) ) return(32)
  if ( grepl('feature_truncation', type) ) return(33)
  if ( grepl('intergenic_variant', type) ) return(34)
  
  if ( grepl('unknown', type) ) return(100)
  cat('Dont know about the mutation type: ', type, '\n')
  return(100)
}

severityToType = function(severity) {
  if ( severity == 1 ) return('transcript_ablation')
  if ( severity == 2 ) return('splice_acceptor_variant')
  if ( severity == 3 ) return('splice_donor_variant')
  if ( severity == 4 ) return('stop_gained')
  if ( severity == 5 ) return('frameshift_variant')
  if ( severity == 6 ) return('stop_lost')
  if ( severity == 7 ) return('initiator_codon_variant')
  if ( severity == 8 ) return('transcript_amplification')
  if ( severity == 9 ) return('inframe_insertion')
  if ( severity == 10 ) return('inframe_deletion')
  if ( severity == 11 ) return('missense_variant')
  if ( severity == 12 ) return('splice_region_variant')
  if ( severity == 13 ) return('incomplete_terminal_codon_variant')
  if ( severity == 14 ) return('stop_retained_variant')
  if ( severity == 15 ) return('synonymous_variant')
  if ( severity == 16 ) return('coding_sequence_variant')
  if ( severity == 17 ) return('mature_miRNA_variant')
  if ( severity == 18 ) return('5_prime_UTR_variant')
  if ( severity == 19 ) return('3_prime_UTR_variant')
  if ( severity == 20 ) return('non_coding_transcript_exon_variant')
  if ( severity == 21 ) return('intron_variant')
  if ( severity == 22 ) return('NMD_transcript_variant')
  if ( severity == 23 ) return('non_coding_transcript_variant')
  if ( severity == 24 ) return('upstream_gene_variant')
  if ( severity == 25 ) return('downstream_gene_variant')
  if ( severity == 26 ) return('TFBS_ablation')
  if ( severity == 27 ) return('TFBS_amplification')
  if ( severity == 28 ) return('TF_binding_site_variant')
  if ( severity == 29 ) return('regulatory_region_ablation')
  if ( severity == 30 ) return('regulatory_region_amplification')
  if ( severity == 31 ) return('regulatory_region_variant')
  if ( severity == 32 ) return('feature_elongation')
  if ( severity == 33 ) return('feature_truncation')
  if ( severity == 34 ) return('intergenic_variant')
  if ( severity == 100 ) return('unknown')
  return('unknown')
}

mergeVEPs = function(files) {
  fileLines = lapply(files, function(file) {
        allFile = read.table(file, sep='\n')
        notComment = allFile[!grepl('^#', allFile)]
        return(notComment)
  })
  
  #remove the ID at first line?
  
  allLines = unique(unlist(fileLines))

  #print to file
}





postAnalyseVEP = function(Rdirectory, plotDirectory, parameters, forceRedo=F) {
  loadMethods()
  assign('.maxCov', parameters$maxCov, envir = .GlobalEnv)
  assign('.systematicVariance', parameters$systematicVariance, envir = .GlobalEnv)
  data = loadData(Rdirectory)
  if ( !('allVariants' %in% names(data)) ) {
    warning('Cant find a saved allVariants.\n')
    return()
  }

  backup = paste0(Rdirectory, '/allVariantsPreVEP.Rdata')
  if ( !file.exists(backup) ) {
    cat('Backing up variants (in case the vep run goes wrong) to', backup, '\n')
    allVariantsPreVEP = data$allVariants
    save('allVariantsPreVEP', file=backup)
  }

  variants = data$allVariants$variants
  sampleMetaData = try(importSampleMetaData(metaDataFile))
  normals = as.logical(gsub('YES', 'T', gsub('NO', 'F', sampleMetaData$NORMAL)))
  names = make.names(sampleMetaData$NAME, unique=T)
  names(normals) = names
  individuals = sampleMetaData$INDIVIDUAL
  timePoints = sampleMetaData$TIMEPOINT
  names(timePoints) = names(individuals) = names
  samplePairs = metaToSamplePairs(names, individuals, normals)
    
  source('runVEP.R')
  vcfDir = paste0(plotDirectory, '/somatics')
  variants = runVEP(variants, vcfDir, cpus=cpus, forceRedoVEP=forceRedo)
  allVariants = data$allVariants
  allVariants$variants = variants
  
  allVariantSaveFile = paste0(Rdirectory, '/allVariants.Rdata')
  save('allVariants', file=allVariantSaveFile)
  
  makeScatterPlots(variants, samplePairs, timePoints, plotDirectory, genome=genome, cpus=cpus, forceRedo=T)
  outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, forceRedo=T)
  outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, forceRedo=T)

}



getMoreVEPinfo = function(variants, plotDirectory) {
  dir = paste0(plotDirectory, '/somatics')
  catLog('Importing more VEP info:\n')
  for ( name in names(variants$variants) ) {
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    VEPdata = try(read.table(VEPfile, fill=T), silent=T)
    if ( class(VEPdata) == 'try-error' ) {
      catLog('failed to read VEP file.\n')
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
    }
    else {
      ID = as.numeric(as.character(VEPdata$V1))
      chr = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[1])
      pos = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[2])
      end = gsub('^.*-', '',pos)
      pos = gsub('-.*$', '',pos)
      var = as.character(VEPdata$V3)
      type = as.character(VEPdata$V7)
      sev = sapply(type, typeToSeverity)
      x = chrToX(chr, as.numeric(pos))
      AApos = VEPdata$V10
      AAbefore = gsub('\\/.*$', '', VEPdata$V11)
      AAafter = gsub('^.*\\/', '', VEPdata$V11)

      lastCol = strsplit(as.character(VEPdata$V14), ';')
      polyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) ) return(gsub('^PolyPhen=', '', strs[grep('PolyPhen=',strs)[1]]))
        else return('')
        })
      numericPolyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) )
          return(as.numeric(gsub('\\)$', '', gsub('^.*\\(', '', strs[grep('PolyPhen=',strs)[1]]))))
        else return(0.05)
        })

      symbol = sapply(lastCol, function(strs) {
        if ( any(grepl('SYMBOL=', strs)) ) return(gsub('^SYMBOL=', '', strs[grep('SYMBOL=',strs)[1]]))
        else return('')
        })
      exon = sapply(lastCol, function(strs) {
        if ( any(grepl('EXON=', strs)) ) return(gsub('^EXON=', '', strs[grep('EXON=',strs)[1]]))
        else return('')
      })
      sift = sapply(lastCol, function(strs) {
        if ( any(grepl('SIFT=', strs)) ) return(gsub('^SIFT=', '', strs[grep('SIFT=',strs)[1]]))
        else return('')
      })
      domain = sapply(lastCol, function(strs) {
        if ( any(grepl('DOMAINS=', strs)) ) return(gsub('\\,.*$', '', gsub('^DOMAINS=', '', strs[grep('DOMAINS=',strs)[1]])))
        else return('')
      })

      namevar = var
      namevar[nchar(namevar) > 1] = paste0('+', substring(var[nchar(var) > 1], 2))
      namevar[namevar == '-'] = paste0('-', pmax(1, as.numeric(end)[namevar == '-']-as.numeric(pos)[namevar == '-']))
      rowNames = paste0(x, namevar)
      qNames = rep('', length(unique(ID)))
      qNames[ID] = rowNames
      
      #find most severe effect
      mostSev = rep('unkown', length(unique(ID)))
      sevScore = rep(100, length(unique(ID)))
      polyPhenRet = rep('', length(unique(ID)))
      exonRet = rep('', length(unique(ID)))
      siftRet = rep('', length(unique(ID)))
      AAposRet = rep('', length(unique(ID)))
      AAbeforeRet = rep('', length(unique(ID)))
      AAafterRet = rep('', length(unique(ID)))
      domainRet = rep('', length(unique(ID)))
      for ( i in 1:nrow(VEPdata) ) {
        sevI = sev[i] - numericPolyPhen[i]
        IDI = ID[i]
        if ( sevI < sevScore[IDI] ) {
          polyPhenRet[IDI] = polyPhen[i]
          exonRet[IDI] = exon[i]
          siftRet[IDI] = sift[i]
          AAposRet[IDI] = AApos[i]
          AAbeforeRet[IDI] = AAbefore[i]
          AAafterRet[IDI] = AAafter[i]
          domainRet[IDI] = domain[i]
        }
      }
      
      #add effect and severity to q
      variants$variants[[name]]$polyPhen = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$sift = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$exon = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$AApos = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$AAbefore = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$AAafter = rep('', nrow(variants$variants[[name]]))
      variants$variants[[name]]$domain = rep('', nrow(variants$variants[[name]]))
      
      variants$variants[[name]][qNames,]$polyPhen = polyPhenRet
      variants$variants[[name]][qNames,]$sift = siftRet
      variants$variants[[name]][qNames,]$exon = exonRet
      variants$variants[[name]][qNames,]$AApos = AAposRet
      variants$variants[[name]][qNames,]$AAbefore = AAbeforeRet
      variants$variants[[name]][qNames,]$AAafter = AAafterRet
      variants$variants[[name]][qNames,]$domain = domainRet
      catLog(length(mostSev), 'VEPed variants.\n')
    }
  }
  catLog('done.\n')
  return(variants)
}
