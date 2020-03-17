# Test

test_that("Test motif empty", {
  seq = "AAATTTCCC"
  code = codes.code(c("GGG"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 0)
  expect_equal(length(ml$outcode), 1)
  expect_equal(ml$outcode[1], 3)
})

test_that("Test motif 1", {
  seq = "AAATTTCCC"
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 1)
  expect_equal(length(ml$outcode), 1)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$outcode[1], 2)
})

test_that("Test motif 2", {
  seq = "AAATTTAAACCC"
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 2)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 1)
  expect_equal(ml$outcode[1], 1)
  expect_equal(ml$outcode[2], 1)
})

test_that("Test motif on odd sequence", {
  seq = "AAATTTAAACCCG" # not a multiple of 3
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 2)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 1)
  expect_equal(ml$outcode[1], 1)
  expect_equal(ml$outcode[2], 2) # or 1
})

test_that("Test motif on odd sequence 2", {
  seq = "AAATTTAAACCCGG" # not a multiple of 3
  code = codes.code(c("AAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 2)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 1)
  expect_equal(ml$outcode[1], 1)
  expect_equal(ml$outcode[2], 2) # or 1
})

# ---------------- Tests for tuple size 4 ----------------

test_that("Test C4 motif", {
  seq = "AAAATTTTAAAACCCC"
  code = codes.code(c("AAAA"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 2)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 1)
  expect_equal(ml$outcode[1], 1)
  expect_equal(ml$outcode[2], 1)
})

test_that("Test C4 motif 2", {
  seq = "AAAATTTTAAAACCCC"
  code = codes.code(c("AAAA", "CCCC"))
  ml = ccmotif.lengths(seq, code)
  expect_equal(length(ml), 2)
  expect_equal(length(ml$incode), 2)
  expect_equal(length(ml$outcode), 1)
  expect_equal(ml$incode[1], 1)
  expect_equal(ml$incode[2], 2)
  expect_equal(ml$outcode[1], 1)
})