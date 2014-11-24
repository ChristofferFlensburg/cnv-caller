
#imports metadata from a tab separated file
importSampleMetaData = function(sampleMetaDataFile) {
  catLog('Loading sample meta data from file...')
  metaData = read.table(sampleMetaDataFile, header=T, as.is=T, fill=T, sep='\t')
  if ( any(!(c('BAM', 'INDIVIDUAL', 'NAME', 'NORMAL') %in% colnames(metaData))) )
    stop('Could not find required columns BAM, INDIVIDUAL, NAME, NORMAL in sample meta data.\n
The meta data file should be a tab separated file with headings.\n')
  catLog('done.\n')
  return(metaData)
}

#helper function to extract sample pairs from meta data
metaToSamplePairs = function(names, individuals, normals) {
  catLog('Deciding which pairs to scatter plot..')
  pairs = list()
  for (individual in unique(individuals)) {
    rows = which(individual == individuals)
    if ( length(rows) < 2 ) next
    for ( row1 in rows ) {
      for ( row2 in rows[rows > row1] ) {
        if ( normals[row2] & !normals[row1] )
          pairs = c(pairs, list(c(names[row2], names[row1])))
        else
          pairs = c(pairs, list(c(names[row1], names[row2])))
      }
    }
  }
  catLog('done.\n')
  return(pairs)
}

#helper function to extract time series frmo emta data.
metaToTimeSeries = function(names, individuals, normals) {
  catLog('Deciding which time series to plot..')
  series = list()
  for (individual in unique(individuals)) {
    rows = which(individual == individuals)
    if ( length(rows) < 3 ) next
    series = c(series, list(names[rows]))
    names(series)[length(series)] = individual
  }
  catLog('done.\n')
  return(series)
}
