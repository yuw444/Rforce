test_that("c return matrix works", {
  x <- matrix(1:120, nrow = 30, ncol = 4)
  y1 <- RC_ColsPermute(x, c(1,3), 926)
  y2 <- RC_ColsPermute(x, c(1,3), 926)
  y3 <- RC_ColsPermute(x, c(1,3), 927)

  expect_equal(sum(x), sum(y1))
  expect_equal(sum(x), sum(y2))

  expect_equal(y1[,2], x[,2])
  expect_equal(y2[,2], x[,2])

  expect_true(any(y2 != y3))

})
