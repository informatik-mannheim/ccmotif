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

#' Workaround for GCATR
#'
#' @param k How many shifts (0, 1, 2)
#' @param codons
#'
#' @return
#' @export
shift_tuples = function(k, codons) {
  sapply(codons, function(codon) {
    n = nchar(codon)
    suffix = substr(codon, k + 1, n)
    prefix = substr(codon, 1, k)
    paste(suffix, prefix, sep = "")
  })
}


#' Create a code from a list of tuples.
#'
#' @param codons Vector of codons as Strings, e.g. c("AUC", "GCA").
#' @param id A brief description of the code.
#'
#' @return List with two elements: 1) id, 2) codons.
#' @export
codes.code = function(codons = c(), id = "unkn.") {
  code = list(id = id, codons = codons)
  class(code) = "codes.code"
  code
}

#' Pretty prints a code.
#'
#' @param c 
#'
#' @return
#' @export
print.codes.code = function(c) {
  s = paste(c$id, ": ", toString(unlist(c$codons)), sep = "")
  print.default(s)
}

#' Read a code from file.
#'
#' The structure of the file is: \cr
#' # header with short description \cr
#' 1, ACA, ACT, ... \cr
#' 2, ACT, ATT, ... \cr
#'
#' @param filename 
#' @param tupleSize The size of the tuples in the code. Default is 3.
#' @param skipCodeIDColumn If true the first column is left out.
#'
#' @return
#' @export
codes.readFile = function(filename, tupleSize = 3, 
                          skipCodeIDColumn = TRUE) {
  # TODO: stringAsFactors does not work (bug?).
  # read.csv reads factors, which are later manually converted
  # to characters.
  A = read.csv(filename,
               skip = 1, header = FALSE, # First line is a comment.
               stringsAsFactors = FALSE,
               comment.char = "#",
               sep = ",",
               strip.white = TRUE) # Codons as text.
  m = dim(A)[1] # Number of rows.
  n = dim(A)[2] # Number of columns.
  s = if (skipCodeIDColumn) 1 else 2 # Index if first column
  e = n # Index of last column
  codes = lapply(1:m, function(i) {
    cid = if (skipCodeIDColumn) i else A[i, 1]
    codes.code(id = cid, codons = as.character(A[i, s:e]))
  })
  codes
}


#' Write codes to a file.
#'
#' @param filename 
#' @param codes List of codes.
#'
#' @export
codes.writeFile = function(filename, codes, header = "# codes") {
  write(header, filename)
  res = lapply(codes, function(c) {
    s = paste(c$id, ", ", toString(unlist(c$codons)), sep = "")
    write(s, filename, append = TRUE, ncolumns = 1000)
  })
}

#' List of all 216 self complementary circular C3 codes.
#' @export
codes.c3 = {
  codes.readFile("data/C3-self-compl-circ-codes.txt",
                    tupleSize = 3, skipCodeIDColumn = FALSE)
}

#' Test whether a list of codons contains a stop codon.
#'
#' @param codonList List of codons (codons as strings).
#'
#' @return True if at least one codon is a stop codon, false otherwise.
#' @export
codes.containsStopCodon = function(codonList) {
  b = sapply(codonList, function(codon) # Map codons to true if stop codon.
    codon != "TAA" &
    codon != "TAG" &
    codon != "TGA"
  )
  !all(b) # If all entries are true, there is no stop codon.
}


#' List of all self complementary circular C3 codes which contains
#' a stop codon.
#'
#' @return
#' @export
codes.c3StopCodon = {
  s = sapply(codes.c3, function(code) codes.containsStopCodon(code$codons))
  codes.c3[s]
}

#' List of all self complementary circular C3 codes which contains
#' no stop codon.
#'
#' @return
#' @export
codes.c3NoStopCodon = {
  s = sapply(codes.c3, function(code) codes.containsStopCodon(code$codons))
  codes.c3[!s]
}

#' List of all self complementary circular C3 codes which contains
#' no stop codon in all three frames.
#'
#' @return
#' @export
codes.c3NoStopCodonAllFrames = {
  f0 = codes.c3NoStopCodon
  f1 = lapply(f0, function(code) {
    id = paste(code$id, "_1", sep = "")
    codes.code(shift_tuples(1, code$codons), id)
  })
  f2 = lapply(f0, function(code) {
    id = paste(code$id, "_2", sep = "")
    codes.code(shift_tuples(2, code$codons), id)
  })
  f = append(f0, f1)
  f = append(f, f2)
  f
}

#' Table for mapping of code numbers to equivalence classes.
#'
#' @return First column: C3 code number, second row: equivalence class.
#' @export
codes.c3.equivmatrix = {
  equiv = c(1, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9, 2, 3, 10, 11, 11, 10, 6, 1, 8, 12, 13, 9, 7, 9, 1, 7, 6, 12, 8, 10, 13, 9, 8, 7, 6, 10, 11, 3, 2, 11, 5, 4, 13, 12, 12, 13, 4, 5, 5, 4, 5, 13, 12, 4, 13, 12, 3, 2, 11, 7, 8, 9, 11, 6, 10, 9, 1, 7, 12, 6, 10, 13, 8, 6, 1, 10, 8, 9, 12, 7, 13, 3, 8, 7, 10, 2, 6, 9, 11, 11, 4, 5, 4, 5, 1, 1, 2, 3, 2, 3, 14, 14, 15, 16, 17, 16, 18, 15, 19, 17, 19, 18, 20, 20, 18, 15, 17, 16, 14, 14, 16, 18, 15, 17, 19, 20, 20, 19, 20, 19, 20, 19, 17, 18, 15, 16, 14, 16, 14, 17, 18, 15, 14, 15, 16, 14, 18, 16, 17, 17, 19, 15, 20, 19, 18, 20, 21, 21, 22, 23, 22, 23, 24, 21, 25, 25, 24, 26, 27, 24, 26, 27, 25, 24, 21, 25, 26, 27, 27, 26, 22, 23, 23, 22, 23, 22, 23, 27, 27, 22, 26, 26, 21, 25, 25, 24, 24, 26, 27, 24, 26, 27, 21, 24, 25, 25, 21, 22, 23, 21, 22, 23)
  cbind(1:216, equiv)
}

#' Equivalence class for a C3 code number.
#'
#' @param cid C3 code number.
#'
#' @return Its equivalence class.
#' @export
codes.c3.equiv = function(cid) codes.c3.equivmatrix[cid, 2]

#' All code numbers for given equivalance class
#'
#' @param eid Equivalance class.
#'
#' @return List of C3 codes numbers.
#' @export
codes.c3.inclass = function(eid) {
  codes.c3.equivmatrix[codes.c3.equivmatrix[, 2] == eid, 1]
}

#' Random code.
#'
#' @param size Number of codons in the code. Default is 20.
#' @param tuplelength Number of bases per tuple. Default is 3.
#' @param id Brief description of the code. Default is "unkn. rnd.".
#' @param uracil If false, all Us are replaced by Ts. Default is false.
#'
#' @return Vector or string based codons.
#' @export
codes.random = function(size = 20, tuplelength = 3,
                        id = "unkn. rnd", uracil = FALSE) {
  bases = c("A", "T", "C", "G")
  l = rep(list(bases), tuplelength) 
  codonsChars = expand.grid(l) # Create all combinations (as chars).
  # Make tuple as strings:
  allCodons = apply(codonsChars, 1, function(k) paste(k, collapse = ""))
  codons = sample(x = allCodons, size = size, replace = FALSE)
  if (uracil == FALSE) {
    codons = sapply(codons, function(codon) gsub("U", "T", codon))
    names(codons) = NULL
  }
  codes.code(id = id, codons = codons)
}

#' Calculates the code usage.
#'
#' @param cu Data frame for codon usages.
#' @param code A code.
#'
#' @return Percentage (0 to 1) of code usage.
#' @export
codes.usage = function(cu, code) {
  idx = match(code$codons, cu$codon) # Index position of codons of code.
  idx = idx[!is.na(idx)] # remove NA if not all codons are used.
  s = sum(cu$freq[idx]) / sum(cu$freq)
  s
}

#' Read a code usage file into a data frame.
#'
#' @param s Species descriptor
#' @param results_path_prefix Prefix for building the the path to
#' the results directory.
#' @param frame Which frame (0, 1 or 2). Default is 0.
#' @param rnd If true random codes are read. Otherwise C3 codes.
#'
#' @return Code usage as a data frame.
#' @export
codes.readCodeUsage = function(s, results_path_prefix,
                               frame = 0, rnd = FALSE) {
  frameLabel = paste("_f", frame, sep = "")
  rndLabel = if (rnd) "_rnd" else ""
  rangeIdx = paste("_", s$idxRange[1], "_", s$idxRange[2],
                   sep = "")
  cuFile = paste("code_usage_", s$fid, rndLabel, rangeIdx,
                 frameLabel, ".csv", sep = "")
  cu = read.csv(paste(results_path_prefix, cuFile, sep = ""),
                colClasses = c("character", "numeric", "numeric"))
  cu
}


codes.intersect = function(c1, c2) intersect(c1$codons, c2$codons)

#' Calculates the intersection of two C3 codes addressed by its id.
#'
#' @param i1 Id of first code (1 to 216)
#' @param i2 Id of second code (1 to 216)
#'
#' @return List of shared codons.
#' @export
codes.intersectByC3Id = function(i1, i2) {
  intersect(codes.c3[[i1]]$codons, codes.c3[[i2]]$codons)
}

#' Shows the differences of two C3 codes addressed by its id.
#'
#' @param i1 Id of first code (1 to 216)
#' @param i2 Id of second code (1 to 216)
#'
#' @return List of different codons.
#' @export
codes.complByC3Id = function(i1, i2) {
  u = union(codes.c3[[i1]]$codons, codes.c3[[i2]]$codons)
  i = intersect(codes.c3[[i1]]$codons, codes.c3[[i2]]$codons)
  list(i1 = setdiff(codes.c3[[i1]]$codons, i),
       i2 = setdiff(codes.c3[[i2]]$codons, i))
}

#' Plot histogram for code usage in all reading frames.
#'
#' @param cu0 Code usage for reading frame 0.
#' @param cu1 Code usage for reading frame 1.
#' @param cu2 Code usage for reading frame 2.
#' @param species A label for a species.
#' @param minCu Left interval border.
#' @param maxCu Right interval border.
#'
#' @return
#' @export
codes.usage.hist = function(cu0, cu1, cu2, species = "unkn.",
                            minCu = .125, maxCu = .5) {
    df = data.frame(code = cu0[, 1],
                    f0in = cu0[, 2],
                    f0out = cu0[, 3],
                    f1in = cu1[, 2],
                    f1out = cu1[, 3],
                    f2in = cu2[, 2],
                    f2out = cu2[, 3])
    breaks = seq(minCu, maxCu, by = .025)
    ggplot(df) +
      geom_histogram(aes(f0in, color = "frame 0"), 
                     breaks = breaks, fill = "blue", alpha = .5) +
      geom_histogram(aes(f1in, color = "frame 1"), linetype = "dashed",
                     breaks = breaks, fill = "red", alpha = .5) +
      geom_histogram(aes(f2in, color = "frame 2"),linetype = "dotted",
                     breaks = breaks, fill = "green", alpha = .5) +
      labs(x = "code usage", y = "count",
           title = paste("Code usage distribution in ", 
                         species, sep = "")) +
      scale_color_manual(values = c("frame 0" = "blue",
                                    "frame 1" = "red",
                                    "frame 2" = "green")) +
      theme_bw()
}