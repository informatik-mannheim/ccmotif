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


#' Bar diagram with circular codes motif length distribution.
#'
#' @param ml
#' @param codeid
#'
#' @export
ccmotif.barplot = function(ml, codeid = "unkn. code") {
  x = ml
  df = data.frame(x = x)
  ggplot(df) +
    geom_bar(aes(x = x), color = "green", fill = "green") +
    labs(x = "motif length", y = "count",
         title = paste("Motif length distribution ", codeid, sep = "")) +
    theme_bw()
}

#' Bar diagram for comparison of two circular codes 
#' motif length distributions.
#' 
#' Usually, a distribution from biological data and one for
#' theoretical data is compared.
#'
#' @param ml Classes for motif lengths
#' @param tml Classes for theoretical distribution of motif lenghts
#' @param codeid A label for the code. Default: "unkn. code".
#'
#' @export
ccmotif.barplotDiff = function(ml, tml, codeid = "unkn. code") {
  m = min(length(ml), length(tml))
  x = 1:m
  df = data.frame(x = x, ml = ml[1:m], tml = tml[1:m])
  ggplot() +
    geom_bar(data = df, aes(x = x, y = ml), stat = "identity",
             color = "blue", fill = "blue", alpha = .5) +
    geom_bar(data = df, aes(x = x, y = tml), stat = "identity",
             color = "red", fill = "red", alpha = .5) +
    labs(x = "motif length", y = "count",
         title = paste("Motif length distribution ", codeid, sep = "")) +
    theme_bw() +
    theme(plot.title = element_text(size = 8))
}

#' Bar diagram with circular codes motif length distribution.
#'
#' @param ml
#' @param codeid
#'
#' @export
ccmotif.barplot2 = function(ml, codeid = "unkn. code") {
  x = ml$incode
  n = ml$outcode
  m = min(length(x), length(n))
  df = data.frame(x = x[1:m], n = n[1:m])
  ggplot(df) +
    geom_bar(aes(x = x), color = "green", fill = "green", alpha = .5) +
    geom_bar(aes(x = n), color = "red", fill = "red", alpha = .5) +
    labs(x = "motif length", y = "count",
         title = paste("Motif length distribution ", codeid, sep = "")) +
    theme_bw()
}

#' Plot histogram for motif lengths in all reading frames.
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
ccmotif.motif.hist = function(ml0, ml1, ml2, species = "unkn.",
                            minCu = 1, maxCu = 2) {
  df = data.frame(code = ml0[, 1],
                  ml0 = ml0[, 2],
                  ml1 = ml1[, 2],
                  ml2 = ml2[, 2])
  breaks = seq(minCu, maxCu, by = .05)
  ggplot(df) +
    geom_histogram(aes(ml0, color = "frame 0"), 
                   breaks = breaks, fill = "blue", alpha = .5) +
    geom_histogram(aes(ml1, color = "frame 1"), linetype = "dashed",
                   breaks = breaks, fill = "red", alpha = .5) +
    geom_histogram(aes(ml2, color = "frame 2"),linetype = "dotted",
                   breaks = breaks, fill = "green", alpha = .5) +
    labs(x = "code usage", y = "count",
         title = paste("Avg. motif length distribution in ", 
                       species, sep = "")) +
    scale_color_manual(values = c("frame 0" = "blue",
                                  "frame 1" = "red",
                                  "frame 2" = "green")) +
    theme_bw()
}
