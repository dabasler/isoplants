context("plant18O_model class")
library(isoplants)

test_that("get_default_parameters returns data.frame", {
  ep <- get_default_parameters()
  expect_s3_class(ep, "data.frame")
})
