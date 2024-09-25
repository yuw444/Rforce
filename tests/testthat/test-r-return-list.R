test_that("c return list works", {
  temp <- RC_ReturnList()
  expect_equal(length(temp), 3)
  expect_equal(temp[[1]], 0:2)
  expect_equal(temp[[2]], 0:3)
  expect_equal(temp[[3]], paste("string", 0:4))

  expect_true(is.numeric(temp[[1]]))
  expect_true(is.integer(temp[[2]]))
  expect_true(is.character(temp[[3]]))
})
