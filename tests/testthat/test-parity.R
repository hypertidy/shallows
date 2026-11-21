test_that("the package reproduces the parity file shared with the JS port", {
  skip_if_not_installed("jsonlite")
  f <- system.file("parity", "parity.json", package = "shallows")
  skip_if(f == "", "parity file not installed")
  ref <- jsonlite::fromJSON(f, simplifyVector = FALSE)$cases
  now <- shallows:::.parity_compute()
  expect_equal(length(now), length(ref))
  for (k in seq_along(now)) {
    expect_equal(now[[k]]$name, ref[[k]]$name)
    for (v in names(now[[k]]$expect)) {
      expect_equal(as.numeric(now[[k]]$expect[[v]]), as.numeric(unlist(ref[[k]]$expect[[v]])),
                   tolerance = 1e-12, info = paste(now[[k]]$name, v))
    }
  }
})
