

runVEP = function(variants, dir, forceRedoVEP=F) {

  catLog('VEP-ing..')
  for ( name in names(variants$variants) ) {
    infile = paste0(dir, '/', name, '.txt')
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    wd = getwd()
    setwd(dir)
    if ( !file.exists(VEPfile) | forceRedoVEP ) {
      catLog(name, '..', sep='')
      a=system(paste0('vep -i ', basename(infile), ' -o ', basename(VEPfile), ' --everything --force_overwrite'), intern=T)
      if ( !any(grepl('Finished', a)) ) warning('VEP run didnt finish!')
    }
    setwd(wd)
  }
  catLog('done.\n')
  
  catLog('Importing VEP results..')
  for ( name in names(variants$variants) ) {
    catLog(name, ': ', sep='')
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    VEPdata = try(read.table(VEPfile, fill=T), silent=T)
    if ( class(VEPdata) == 'try-error' ) {
      catLog('no variants..')
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
        sevI = sev[i]
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
      catLog(length(mostSev), 'VEPed variants..')
    }
  }
  catLog('done.\n')
  return(variants)
}


typeToSeverity = function(type) {
  if ( grepl('frameshift_variant', type) ) return(1)
  if ( grepl('stop_gained', type) ) return(2)
  if ( grepl('missense', type) ) return(3)
  if ( grepl('stop_lost', type) ) return(4)
  if ( grepl('inframe_deletion', type) ) return(5)
  if ( grepl('inframe_insertion', type) ) return(6)
  if ( grepl('splice_acceptor_variant', type) ) return(7)
  if ( grepl('splice_donor_variant', type) ) return(8)
  if ( grepl('splice_region_variant', type) ) return(9)
  if ( grepl('TF_binding_site_variant', type) ) return(10)
  if ( grepl('3_prime_UTR_variant', type) ) return(11)
  if ( grepl('non_coding_exon_variant', type) ) return(12)
  if ( grepl('5_prime_UTR_variant', type) ) return(13)
  if ( grepl('mature_miRNA_variant', type) ) return(14)
  if ( grepl('nc_transcript_variant', type) ) return(15)
  if ( grepl('regulatory_region_variant', type) ) return(16)
  if ( grepl('intron_variant', type) ) return(17)
  if ( grepl('stop_retained_variant', type) ) return(18)
  if ( grepl('synonymous_variant', type) ) return(19)
  if ( grepl('incomplete_terminal_codon_variant', type) ) return(20)
  if ( grepl('upstream_gene_variant', type) ) return(21)
  if ( grepl('downstream_gene_variant', type) ) return(22)
  if ( grepl('intergenic_variant', type) ) return(23)
  
  if ( grepl('unknown', type) ) return(100)
  cat('Dont know about the mutation type: ', type, '\n')
  return(100)
}

severityToType = function(severity) {
  if ( severity == 1 ) return('frameshift_variant')
  if ( severity == 2 ) return('stop_gained')
  if ( severity == 3 ) return('missense')
  if ( severity == 4 ) return('stop_lost')
  if ( severity == 5 ) return('inframe_deletion')
  if ( severity == 6 ) return('inframe_insertion')
  if ( severity == 7 ) return('splice_acceptor_variant')
  if ( severity == 8 ) return('splice_donor_variant')
  if ( severity == 9 ) return('splice_region_variant')
  if ( severity == 10 ) return('TF_binding_site_variant')
  if ( severity == 11 ) return('3_prime_UTR_variant')
  if ( severity == 12 ) return('non_coding_exon_variant')
  if ( severity == 13 ) return('5_prime_UTR_variant')
  if ( severity == 14 ) return('mature_miRNA_variant')
  if ( severity == 15 ) return('nc_transcript_variant')
  if ( severity == 16 ) return('regulatory_region_variant')
  if ( severity == 17 ) return('intron_variant')
  if ( severity == 18 ) return('stop_retained_variant')
  if ( severity == 19 ) return('synonymous_variant')
  if ( severity == 20 ) return('incomplete_terminal_codon_variant')
  if ( severity == 21 ) return('upstream_gene_variant')
  if ( severity == 22 ) return('downstream_gene_variant')
  if ( severity == 23 ) return('intergenic_variant')
  if ( severity == 100 ) return('unknown')
  return('unknown')
}
