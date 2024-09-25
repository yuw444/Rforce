test_that("c return vector works", {

  expect_equal(RC_Cummax(c(1, 3, 2, 2, 5)), c(1, 3, 3, 3, 5))
  
})
