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

#' Code usage in different positions and frames.
#'
#' @param seq Sequence to be analyzed.
#' @param codes List of three codes.
#'
#' @return List of size 9. The number indicates the offset 1 to 9. 
#' @export
ccmotif.allFrames = function(seq, codes) {
  W = 10 # minimal sequence length
  if (nchar(seq) < W) {
    #return(t(c(0, 0, 0)))
    stop("Sequence to short.")
  }
  n = nchar(seq)
  
  binSeq = function(i) {
    e = n - (W - i) # end index
    s = substr(seq, i, e)
    f = ((i - 1) %% 3) + 1 # frame indexed by 1, 2, 3
    ccmotif.seq2cc(s, codes[[f]])
  }
  
  b = lapply(1:9, function(i) binSeq(i))
  #b
  cbind(b[[1]], b[[2]], b[[3]], b[[4]], b[[5]], b[[6]],
        b[[7]], b[[8]], b[[9]])
}
