# Copyright 2018-19 by the authors.
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

#' Relative frequency of code usage.
#'
#' Calculates the relative frequency of the codon usage of a code X.
#'
#' @param ml Motif lengths.
#'
#' @return Vector with two elements:
#' 1. usage for codon in code, 2. usage for codons not in code X.
#' @export
ccmotif.codonfreq = function(ml) {
  c_in = sum(ml$incode)
  c_out = sum(ml$outcode)
  c_total = c_in + c_out
  c(`in` = c_in / c_total, out = c_out / c_total)
}

#' Analysis of motif lengths in a list of fasta sequences.
#'
#' The function scans over all DNA (RNA) sequences (\code{seqlist}) and
#' calculates the motif lengths for all codes (\code{codes}). Depending
#' on the parameters this might take a while.
#' This scan writes two comma separated files:
#' 1) table with relative frequency of codons usage. Filename is
#' \code{code_usage_<seqid>.csv}
#' 2) two files with the motif lengths of 1) the code and 2) codons not
#' in the code. Filename are
#' \code{motif_lengths_<seqid>.csv} and \code{motif_lengths_non_<seqid>.csv}
#'
#' @param ff Entries of a fasta file set.
#' @param codes List of codes to analyze.
#' @param frame Frame to analyze (0, 1 or 2). Default is 0.
#' @param seqid A short label (one word) for the sequences.
#' @param tmpDir Path to directory to save files. Default: "./".
#'
#' @export
ccmotif.scan.fasta = function(ff, codes, frame = 0,
                              label = "unkn", tmpDir = "./") {

  frameLabel = paste("_f", frame, sep = "")
  codenumbers = sapply(codes, function(c) c$id)

  codonfreqsin = c()
  codonfreqsout = c()
  #mlinpercode = list(c(1, 2), "a") # dummy
  #mloutpercode = list(c(1, 2), "a") # dummy
  mlinpercode = list() # dummy
  mloutpercode = list() # dummy
  
  i = 1
  for (code in codes) {
    ml = ccmotif.lengths.fasta(ff, code, frame)
    mlinpercode[[i]] = ml$incode
    mloutpercode[[i]] = ml$outcode

    codonfreq = ccmotif.codonfreq(ml)
    codonfreqsin = c(codonfreqsin, codonfreq["in"])
    codonfreqsout = c(codonfreqsout, codonfreq["out"])

    i = i + 1
  }

  # Write motif length distribution file:
  writeml = function(mlpercode, what) {
    # First, we need to have the motif lengths vector the same size:
    n = max(sapply(mlpercode, length)) # maximum length of motif lengths.
    i = 1
    for (ml in mlpercode) {
      length(mlpercode[[i]]) = n # Rezize the vectors in the list.
      i = i + 1
    }
    names(mlpercode) = codenumbers
    filename = paste(tmpDir, "motif_lengths", what, label,
                     frameLabel, ".csv.zip", sep = "")
    write.csv(mlpercode, gzfile(filename), row.names = FALSE)
    # write.csv(mlpercode, filename, row.names = FALSE)
  }
  writeml(mlinpercode, "_")
  writeml(mloutpercode, "_non_")

  # Write code usage file:
  rf = cbind(codonfreqsin, codonfreqsout)
  rownames(rf) = codenumbers
  colnames(rf) = c("in", "out")
  filename = paste(tmpDir, "code_usage_", label, frameLabel,
                   ".csv", sep = "")
  write.csv(rf, filename)
}

#' Motif lengths in sequences of a FASTA file.
#'
#' @param seqlist Fasta file set (as returned by seqinr fasta.read)
#' @param code
#' @param frame Frame to use (0, 1 or 2). Default: 0
#'
#' @return
#' @export
ccmotif.lengths.fasta = function(seqlist, code, frame = 0) {
  incode = c()
  outcode = c()
  for (s in seqlist) {
    sequ = s[[1]]
    # Remove leading elements to get correct frame:
    sequ = substr(sequ, (frame + 1), nchar(sequ))
    r = ccmotif.lengths(sequ, code = code)
    # r = ccmotif.rlengths(sequ, code = code) # TODO confirm
    incode = c(incode, r$incode)
    outcode = c(outcode, r$outcode)
  }
  list(incode = incode, outcode = outcode)
}

#' Map a sequence to classes "in code" (1) or "out of code" (0).
#'
#' @param seq Sequence
#' @param code Code as a list with entries id, and tuples,
#' e.g. list(id = "foo", codons = c("ACU", "GAG")).
#'
#' @return Vector of 0s and 1s.
#' @export
ccmotif.seq2cc = function(seq, code) {
  tuplelength = nchar(code$codons[1])
  # the indices where each substr will start
  starts = seq(1, nchar(seq), by = tuplelength)

  # chop it up
  s = sapply(starts, function(ii) substr(seq, ii, ii + tuplelength - 1))

  m = sapply(s, function(codon) {
    if (is.element(codon, code$codons)) 1 else 0
  })
  m
}

#' Motif lengths in a sequence.
#'
#'
#' @param seq Sequence
#' @param code Code as a list with entries id, and codons,
#' e.g. list(id = "foo", codons = c("ACU", "GAG")).
#'
#' @return Data structure with two elements: 1) incode is a vector of
#' motifs (run lengths) and 2) outcode is the vector of run lengths which are
#' not part of the code.
#' @export
ccmotif.lengths = function(seq, code) {
  cc = ccmotif.seq2cc(seq, code) # Sequence of 0,1.
  r = rle(cc) # Calulate run lengths
  # Remember/store only values but no additional meta data:
  x = as.vector(r$lengths[r$values == 1]) # in code
  n = as.vector(r$lengths[r$values == 0]) # not in code
  l = list(incode = x, outcode = n)
  class(l) = append("ccmotif.lengths", class(l))
  l
}

#' Motif lengths in list of a sequences.
#'
#'
#' @param seqlist List of sequences
#' @param code Code as a list with entries id, and codons,
#' e.g. list(id = "foo", codons = c("ACU", "GAG")).
#'
#' @return Data structure with two elements: 1) incode is a vector of
#' motifs (run lengths) and 2) outcode is the vector of run lengths which are
#' not part of the code.
#' @export
ccmotif.lengthslist = function(seqlist, code) {
  l = list(incode = c(), outcode = c())
  for (seq in seqlist) {
    h = ccmotif.lengths(seq, code)  
    l$incode = c(l$incode, h$incode)
    l$outcode = c(l$outcode, h$outcode)
  }
  class(l) = append("ccmotif.lengths", class(l))
  l
}

#' Read a motif lengths file into a data frame.
#'
#' @param s Species descriptor
#' @param results_path_prefix Prefix for building the the path to
#'  the results directory.
#' @param frame Which frame (0, 1 or 2). Default is 0.
#' @param rnd If true, motif lengths for random codes are read,
#' otherwise normal ones. Default is FALSE.
#' @param motif If true motif sequences are read, otherwise non-motif ones.
#'
#' @return Motif lengths as a data frame.
#' @export
ccmotif.readMotifLengths = function(s, results_path_prefix,
                                    frame = 0, rnd = FALSE, motif = TRUE,
                                    newZipVersion = FALSE) {
  frameLabel = paste("_f", frame, sep = "")
  rndLabel = if (rnd) "_rnd" else ""
  motifLabel = if (motif) "" else "non_"
  rangeIdx = paste("_", s$idxRange[1], "_", s$idxRange[2],
                   sep = "")
  mlFile = paste("motif_lengths_", motifLabel, s$fid, rndLabel, rangeIdx,
                 frameLabel, ".csv", sep = "")
  print(mlFile)
  # Zip file also possible. Takes 50% more time.
  # There seems to be a bug in unz() when reading freshly zip files.
  # The new version can read those files but it cannot read older ones.
  # Very strange. Let's see what was wrong..
  ml = if (newZipVersion) {
    read.csv(paste(results_path_prefix, mlFile, ".zip", sep = ""),
               header = TRUE, stringsAsFactors = FALSE)
  } else {
    read.csv(unz(paste(results_path_prefix, mlFile, ".zip", sep = ""), mlFile),
                header = TRUE, stringsAsFactors = FALSE)
  }
  # ml = read.csv(paste(results_path_prefix, mlFile, sep = ""),
  #              header = TRUE, stringsAsFactors = FALSE)
  ml
}

#' Create histogram classes for motifs (runs).
#'
#' @param ml List of motif lengths
#' @param p Probability (relative frequency) of code usage
#'
#' @return List of two elements: 1) sample histogram,
#' 2) geometric distribution histogram
#' @export
ccmotif.classes = function(ml, p, K = 8) {
  ml = na.omit(ml) # Samples
  n = length(ml) # Number of motif lengths
  q = 1 - p
  m = 1:(K - 1) # classes
  # Classes according to theoretical geometric distribution:
  n0 = sapply(m, function(k) n * q * (1 - q)^(k - 1))
  h0 = sapply(K:50, function(k) n * q * (1 - q)^(k - 1))
  n0 = c(n0, sum(h0)) # Remaing classes
  # Observed classes:
  N = sapply(m, function(k) length(which(ml %in% k)))
  N0 = sapply(K:50, function(k) length(which(ml %in% k)))
  N = c(N, sum(N0))
  list(sample = N, geom = n0)
}


#' Motif lengths in a sequence (Java).
#'
#'
#' @param seq Sequence
#' @param code Code as a list with entries id, and codons,
#' e.g. list(id = "foo", codons = c("ACU", "GAG")).
#'
#' @return Data structure with two elements: 1) incode is a vector of
#' motifs (run lengths) and 2) outcode is the vector of run lengths which are
#' not part of the code.
#' @export
ccmotif.jlengths = function(seq, code) {
  ccm = new(J("bio/gcat/ccmotif/CCMotif"))
  ja = ccm$motifLengthsR(seq, .jarray(code$codons), code$id)
  # ra = convertToJava(.jevalArray(ja))
  ra = .jevalArray(ja, simplify = TRUE)
  l1 = ra[1, ]
  l2 = ra[2, ]
  # Workaround for Java - R interface. Last value can be 0
  # indicating that this value does not exist:
  if (l1[length(l1)] == 0) {
    l1 = l1[-length(l1)]
  }
  if (l2[length(l2)] == 0) {
    l2 = l2[-length(l2)]
  }
  list(incode = l1, outcode = l2)
}


#' Next motif length.
#'
#'   + 1 2 3
#' 1 | 1 2 3
#' 2 | 1 3 2
#' 3 | 1 2 3
#' @param seq Sequence
#' @param code Code as a list with entries id, and codons,
#' e.g. list(id = "foo", codons = c("ACU", "GAG")).
#' @param max Maximal length. Default: 25.
#' @return
#' @export
ccmotif.nextLength = function(seq, code, max = 25, start = 1, fwd = 1) {
  cc = ccmotif.seq2cc(seq, code) # Sequence of 0,1.
  r = rle(cc) # Calulate run lengths

  M = matrix(0, nrow = max, ncol = max) # Result matrix of dim max x max.
  i0 = if (r$values[1] == 1) 1 else 2 # Start position of in-code.
  N = length(r$values) - i0 + 1
  if (N - fwd < i0) {
    return(M)
  }

  w = 0
  for (i in seq(i0, N - fwd, 2)) {
    k = r$lengths[i] # Motif length of current in-code
    l = r$lengths[i + fwd] # Motif length of next out-code
    if (k <= max && l <= max) {
      M[k, l] = M[k, l] + 1 # Increase this motif length
    } else {
      if (w == 0) {
        warning(paste("Max.length of ", max, " exceed: ", k, "; ", l))
      }
      w = w + 1
    }
  }
  M
}
