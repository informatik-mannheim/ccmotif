# Tests for codons.

test_that("Test codon 1", {
  seq = c("AAA")
  codons = codon.split(seq)
  expect_equal(length(codons), 1)
  expect_equal(codons[1], "AAA")
})

test_that("Test codon 2", {
  seq = c("AAATTT")
  codons = codon.split(seq)
  expect_equal(length(codons), 2)
  expect_equal(codons[1], "AAA")
  expect_equal(codons[2], "TTT")
})

test_that("Test codon incomplete 0.1", {
  seq = c("G")
  codons = codon.split(seq)
  expect_equal(length(codons), 0)
})

test_that("Test codon incomplete 2.1", {
  seq = c("AAATTTG")
  codons = codon.split(seq)
  expect_equal(length(codons), 2)
  expect_equal(codons[1], "AAA")
  expect_equal(codons[2], "TTT")
})

test_that("Test codon incomplete 2.2", {
  seq = c("AAATTTGG")
  codons = codon.split(seq)
  expect_equal(length(codons), 2)
  expect_equal(codons[1], "AAA")
  expect_equal(codons[2], "TTT")
})