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

library(forcats)
library(dplyr)

#' Split string sequence into tupels (codons).
#'
#' If the last tuple has less than `size` bases it is skipped.
#'
#' @param seq Sequence as a string.
#' @param size Tuple size. Default: 3 for codons.
#'
#' @return Vector of tuples. Each tuple is a string.
#' @export
codon.split = function(seq, size = 3) {
  n = as.integer(nchar(seq) / 3) * 3 # number of tuples
  if (n < size) { # if there is no tuple at all...
    return(c()) # an empty vector is returned.
  }
  starts = seq(1, n, by = size)
  # chop it up:
  sapply(starts, function(ii) substr(seq, ii, ii + (size - 1)))
}

#' Split string sequences in a list into tupels (codons).
#'
#' The bases in every list are converted into tuples and
#' put together into a list of tuples (unflatten).
#'
#' @param seqlist List of sequences (e.g. from a fasta file).
#'
#' @return Vector of tuples. Each tuple is a string.
#' @export
codon.splitlist = function(seqlist, size = 3) {
  u = c()
  for (s in seqlist) {
    sequ = s[[1]]
    r = codon.split(sequ, size)
    u = c(u, r)
  }
  u
}

#' Frequencies of codon usage.
#'
#' Calculates the absolute frequencies of codons in a codon sequence.
#'
#' @param codSeq Sequence of tuples (codons)
#'
#' @return data frame with two columns:
#' 1. codon, 2. number of codons in sequence
#' @export
codon.usage = function(codSeq) {
  df = data.frame(cf = factor(codSeq))
  df = as.data.frame(table(df$cf))
  colnames(df) = c("codon", "freq")
  class(df) = append("codon.usage", class(df)) # order is important.
  df
}

#' Bar diagram of codon usage.
#'
#' @param cu Code usage
#' @param species The species. Default is "unknown".
#'
#' @export
plot.codon.usage = function(cu, species = "unknown") {
  # cu = cu[order(cu[, 1], decreasing = TRUE), ]
  #cu %>% 
  #  mutate(name = fct_reorder(codon, freq, .desc = F)) %>%
  p = ggplot(cu, aes(x = codon, y = freq)) +
    geom_bar(stat = "identity", fill = "darkgray") +
    labs(x = "codon", y = "codon count",
         title = paste("Codon frequencies for ", species, sep = "")) +
    # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_flip() + theme_bw()
  plot(p)
}
