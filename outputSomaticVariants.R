require(WriteXLS)


#prints the somatic variants to an excel sheet.
outputSomaticVariants = function(variants, genome, plotDirectory, cpus=cpus, v, forceRedo=forceRedoOutputSomatic) {
  outfile = paste0(plotDirectory, '/somaticVariants.xls')
  if ( (!file.exists(outfile) | forceRedo) ) {
    somatics = list()
    catLog('Printing somatic variants to ', outfile, '.\n', sep='')
    for ( sample in names(variants$variants) ) {
      catLog(sample, '..', sep='')
      q = variants$variants[[sample]]
      somaticP = q$somaticP
      toReturn = which(somaticP > 0)
      toReturn = toReturn[order(somaticP[toReturn], decreasing=T)]
      q = q[toReturn,]
      SNPs = variants$SNPs[variants$SNPs$x %in% q$x,]
      somatic = data.frame(
        chr=xToChr(q$x, genome),
        start=xToPos(q$x, genome),
        end=xToPos(q$x, genome),
        reference=q$reference,
        variant=q$variant,
        inGene=SNPs[as.character(q$x),]$inGene,
        f=q$var/q$cov,
        cov=q$cov,
        ref=q$ref,
        var=q$var,
        flag=q$flag,
        pbq=q$pbq,
        pmq=q$pmq,
        psr=q$psr,
        somaticP=q$somaticP,
        row.names=rownames(q))
      somatics[[sample]] = somatic
    }
    WriteXLS('somatics', outfile)
    catLog('done!\n')
  }
}
