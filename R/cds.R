# Copyright 2019 by the authors.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#' @author Markus Gumbel


#' Check a fasta file for structure of coding sequence.
#'
#' @param filename Fasta file name.
ccmotif.checkCDS.fastafile = function(filename) {
  ff = read.fasta(filename, as.string = TRUE,
                  forceDNAtolower = FALSE)
  print(paste("Verify ", filename, " for CDS.", sep = ""))
  ccmotif.checkCDS.fasta(ff, filename, paste(filename, "_quality.txt"))
  print("Done.")
}

#' Check a fasta file entry for structure of coding sequence.
#'
#' @param ff Fasta file entry returned by `read.fasta`
#' @param source The source location of the fasta file.
#' @param logfilename File name for the logging data.
ccmotif.checkCDS.fasta = function(ff, source, logfilename = NULL) {
  cdsCounts = list(start = list(), stop = list())
  for (f in ff) {
    cdsCounts = ccmotif.checkCDS(f[[1]], attributes(f)$name, 
                                 source, cdsCounts, logfilename)
  }
  cat("Start codons:\n")
  print(unlist(cdsCounts$start))
  cat("Stop codons:\n")
  print(unlist(cdsCounts$stop))
}

ccmotif.checkCDS = function(seq, id, source, cdsCounts, filename = NULL) {
  startCounts = cdsCounts$start
  stopCounts = cdsCounts$stop
  firstCodon = substr(seq, 1, 3)
  n = nchar(seq)
  lastCodon = substr(seq, n - 2, n)
  
  k = if (is.null(startCounts[firstCodon][[1]])) 0 else startCounts[firstCodon][[1]]
  startCounts[firstCodon] = k + 1          
  l = if (is.null(stopCounts[lastCodon][[1]])) 0 else stopCounts[lastCodon][[1]]
  stopCounts[lastCodon] = l + 1
  
  if (!(firstCodon == "ATG" || firstCodon == "AUG")) {
    cdsInfo = paste("Unknown start codon: ", firstCodon, " in ",
                    source, ": ", id, "\n", sep = "")
    if (is.null(filename)) {
      print(cdsInfo)
    } else {
      cat(cdsInfo, file = filename, append = TRUE)
    }
  }
  
  if (!(lastCodon == "TAA" || lastCodon == "UAA"
      || lastCodon == "TAG" || lastCodon == "UAG"
      || lastCodon == "TGA" || lastCodon == "UGA")) {
    cdsInfo = paste("Unknown stop codon: ", lastCodon, " in ",
                    source, ": ", id, "\n", sep = "")
    if (is.null(filename)) {
      print(cdsInfo)
    } else {
      cat(cdsInfo, file = filename, append = TRUE)
    }
  }
  list(start = startCounts, stop = stopCounts)
}
