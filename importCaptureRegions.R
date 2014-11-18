
#import a capture region with gc content and gene names.
importCaptureRegions = function(file='captureRegions.gc.bed', genome='hg19', gcColumn=5) {
  require(rtracklayer)
  cR = read.table(file, fill=T, quote='')
  chrL = chrLengths(genome)
  use = gsub('chr', '',cR$V1) %in% names(chrL)
  cR = cR[use,]
  if ( is.na(gcColumn) )
    gr = GRanges(ranges = IRanges(start=cR$V2, end = cR$V3),
      seqnames=gsub('chr', '',cR$V1), seqlengths=chrL, region = gsub(',', '', gsub("\"", '', as.character(cR$V4))))
  else
    gr = GRanges(ranges = IRanges(start=cR$V2, end = cR$V3),
      seqnames=gsub('chr', '',cR$V1), seqlengths=chrL, region = gsub(',', '', gsub("\"", '', as.character(cR$V4))),
      gc = cR[[gcColumn]])

  names(gr) = gr$region
  return(gr)
}

#helper functions that returns the lengths of the chromosomes of a genome.
humanAllChrLengths = function() {
lengths = c(249250621, 106433, 547496, 243199373, 198022430, 191154276, 590426, 189789, 191469, 180915260, 171115067, 4622290, 4795371,
4610396, 4683263, 4833398, 4611984, 4928567, 159138663, 182896, 146364022, 38914, 37175, 141213431, 90085, 169874, 187035, 36148,
135534747, 135006516, 40103, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 1680828, 37498, 81310, 174588, 41001,
78077248, 4262, 59128983, 92689, 159169, 63025520, 48129895, 27682, 51304566, 155270560, 59373566, 166566, 186858, 164239, 137718,
172545, 172294, 172149, 161147, 179198, 161802, 155397, 186861, 180455, 179693, 211173, 15008, 128374, 129120, 19913, 43691, 27386,
40652, 45941, 40531, 34474, 41934, 45867, 39939, 33824, 41933, 42152, 43523, 43341, 39929, 36651, 38154, 36422, 39786, 38502, 16571)

names(lengths) = c('1','1_gl000191_random','1_gl000192_random','2','3','4','4_ctg9_hap1','4_gl000193_random','4_gl000194_random',
'5','6','6_apd_hap1','6_cox_hap2','6_dbb_hap3','6_mann_hap4','6_mcf_hap5','6_qbl_hap6','6_ssto_hap7','7','7_gl000195_random',
'8','8_gl000196_random','8_gl000197_random','9','9_gl000198_random','9_gl000199_random','9_gl000200_random','9_gl000201_random',
'10','11','11_gl000202_random','12','13','14','15','16','17','17_ctg5_hap1','17_gl000203_random','17_gl000204_random',
'17_gl000205_random','17_gl000206_random','18','18_gl000207_random','19','19_gl000208_random',
'19_gl000209_random', '20', '21', '21_gl000210_random', '22', 'X', 'Y', 'Un_gl000211', 'Un_gl000212', 'Un_gl000213', 'Un_gl000214',
'Un_gl000215', 'Un_gl000216', 'Un_gl000217', 'Un_gl000218', 'Un_gl000219', 'Un_gl000220', 'Un_gl000221', 'Un_gl000222', 'Un_gl000223',
'Un_gl000224', 'Un_gl000225', 'Un_gl000226', 'Un_gl000227', 'Un_gl000228', 'Un_gl000229', 'Un_gl000230', 'Un_gl000231', 'Un_gl000232',
'Un_gl000233', 'Un_gl000234', 'Un_gl000235', 'Un_gl000236', 'Un_gl000237', 'Un_gl000238', 'Un_gl000239', 'Un_gl000240', 'Un_gl000241',
'Un_gl000242', 'Un_gl000243',  'Un_gl000244', 'Un_gl000245', 'Un_gl000246', 'Un_gl000247', 'Un_gl000248', 'Un_gl000249', 'M')

  return(lengths)
}
humanChrLengths = function() {
  humanAllChrLengths()[as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y', 'M'))]
}
chrLengths = function(genome='hg19') {
  if ( genome == 'hg19' ) return(humanChrLengths())
  else stop('chrLengths doesnt know about genome', genome, '\n')
}
