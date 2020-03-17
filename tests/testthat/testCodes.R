# Test for codes.

test_that("Test code X1", {
  c = codes.c3[[1]]
  expect_equal(c$id, "X1")
  expect_equal(length(c$codons), 20)
  expect_equal(c$codons[1], "AAC")
  expect_equal(c$codons[20], "TGA")
})

test_that("Test code X23", {
  c = codes.c3[[23]]
  expect_equal(c$id, "X23")
  expect_equal(length(c$codons), 20)
  expect_equal(c$codons[1], "AAC")
  expect_equal(c$codons[5], "ATT")
  expect_equal(c$codons[10], "GAC")
  expect_equal(c$codons[15], "GGT")
  expect_equal(c$codons[20], "TTC")
})

test_that("Test random code 20", {
  c = codes.random(20)
  expect_equal(c$id, "unkn. rnd")
  expect_equal(length(c$codons), 20)
  expect_equal(nchar(c$codons[1]), 3)
})

test_that("Test random code 10", {
  c = codes.random(10)
  expect_equal(c$id, "unkn. rnd")
  expect_equal(length(c$codons), 10)
})

test_that("Test random code of length 4", {
  c = codes.random(size = 12, tuplelength = 4)
  expect_equal(c$id, "unkn. rnd")
  expect_equal(length(c$codons), 12)
  expect_equal(nchar(c$codons[1]), 4)
})

test_that("Test contains no stop codon 1", {
  r = codes.containsStopCodon(c("AUG"))
  expect_equal(r, FALSE)
})

test_that("Test contains no stop codon 2", {
  r = codes.containsStopCodon(c("CCC", "GGG"))
  expect_equal(r, FALSE)
})

test_that("Test contains stop codon 1", {
  r = codes.containsStopCodon(c("TAA"))
  expect_equal(r, TRUE)
})

test_that("Test contains stop codon 2", {
  r = codes.containsStopCodon(c("CAG", "TAA", "AUC"))
  expect_equal(r, TRUE)
})

test_that("Test if X23 contains stop codon", {
  r = codes.containsStopCodon(codes.c3[[23]]$codons)
  expect_equal(r, FALSE)
})
