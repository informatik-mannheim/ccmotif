# Tests for code usage.

test_that("Test code usage 1", {
  seq = c("AAA")
  cu = codon.usage(seq)
  expect_equal(length(cu[, 1]), 1)
  expect_equal(cu[1, 2], 1)
})

test_that("Test code usage 2", {
  seq = c("AAA", "AAA")
  cu = codon.usage(seq)
  expect_equal(length(cu[, 1]), 1)
  expect_equal(cu[1, 2], 2)
})

test_that("Test code usage 2 + 1", {
  seq = c("AAA", "AAA", "TTT")
  cu = codon.usage(seq)
  expect_equal(length(cu[, 1]), 2)
  expect_equal(cu[1, 2], 2)
  expect_equal(cu[2, 2], 1)
})

test_that("Test code usage 2 + 1 alternative 1", {
  seq = c("AAA", "TTT", "AAA")
  cu = codon.usage(seq)
  expect_equal(length(cu[, 1]), 2)
  expect_equal(cu[1, 2], 2)
  expect_equal(cu[2, 2], 1)
})

test_that("Test code usage 2 + 1 alternative 2", {
  seq = c("TTT", "AAA", "AAA")
  cu = codon.usage(seq)
  expect_equal(length(cu[, 1]), 2)
  expect_equal(cu[1, 2], 2)
  expect_equal(cu[2, 2], 1)
})

test_that("Test code usage with split", {
  seq = "TTTAAAAAA"
  cu = codon.usage(codon.split(seq, size = 3))
  expect_equal(length(cu[, 1]), 2)
  expect_equal(cu[1, 2], 2)
  expect_equal(cu[2, 2], 1)
})

test_that("Test code usage with motif lengths", {
  seq = "TTTAAAAAA"
  code = codes.code(c("TTT"))
  ml = ccmotif.lengths(seq, code)
  cu = ccmotif.codonfreq(ml)
  expect_equal(length(cu), 2)
  expect_equal(cu["in"], c("in" = 1 / 3))
  expect_equal(cu["out"], c("out" = 2 / 3))
})

test_that("Test code usage with motif lengths 2", {
  seq = "TTTAAAAAATTT"
  code = codes.code(c("TTT"))
  ml = ccmotif.lengths(seq, code)
  cu = ccmotif.codonfreq(ml)
  expect_equal(length(cu), 2)
  expect_equal(cu["in"], c("in" = 1 / 2))
  expect_equal(cu["out"], c("out" = 1 / 2))
})

test_that("Test code usage with motif lengths 3", {
  seq = "TTTAAAAAATTTGGG"
  code = codes.code(c("TTT"))
  ml = ccmotif.lengths(seq, code)
  cu = ccmotif.codonfreq(ml)
  expect_equal(length(cu), 2)
  expect_equal(cu["in"], c("in" = 2 / 5))
  expect_equal(cu["out"], c("out" = 3 / 5))
})

test_that("Test code usage with motif lengths 4", {
  seq = "TTTAAAAAATTTGGG"
  code = codes.code(c("TTT", "GGG"))
  ml = ccmotif.lengths(seq, code)
  cu = ccmotif.codonfreq(ml)
  expect_equal(length(cu), 2)
  expect_equal(cu["in"], c("in" = 3 / 5))
  expect_equal(cu["out"], c("out" = 2 / 5))
})