

#These are debugging functions, and other methods to handle errors and problems in general



dumpInput = function(Rdirectory, inputList) {
  catLog('Dumping input for debugging purposes..')
  debugDir = paste0(Rdirectory, '/debug')
  if ( !file.exists(debugDir) ) dir.create(debugDir)
  dumpDir = paste0(debugDir, '/inputDump-', Sys.time(), '-', abs(rnorm(1, 0, 0.1)))
  dir.create(dumpDir)
  for ( input in names(inputList) ) {
    catLog(input, '..', sep='')
    save(input, file=paste0(dumpDir, '-', input, '.Rdata'))
  }
  catLog('done.\n')
}
