# The helpers for forked processes, which prepare_for_ORFquant() uses and
# ORFquant imports. From ORFquant's tests (tests/testthat/test-helpers.R).

test_that("mc_lapply() is lapply() with each function in a forked process", {
  skip_on_os("windows")
  fs <- list(function() 1, function() "a", function() NULL)
  expected <- lapply(fs, function(f) f())
  expect_identical(mc_lapply(fs, 1), expected)
  expect_identical(mc_lapply(fs, 2), expected)

  # The warning of a forked process is given again, and its error stops
  expect_warning(res <- mc_lapply(list(function() { warning("careful"); 1 }, function() 2), 2),
                 "careful")
  expect_identical(res, list(1, 2))
  expect_error(suppressWarnings(mc_lapply(list(function() 1, function() stop("boom")), 2)), "boom")
  # A process that ends without a result, as when the system kills it for memory
  expect_error(suppressWarnings(mc_lapply(list(function() 1,
                                               function() system2("kill", c("-9", Sys.getpid()))), 2)),
               "ended without a result")
})

test_that("mc_job() evaluates in a forked process with more than 1 core", {
  skip_on_os("windows")
  expect_equal(mc_value(mc_job(1 + 1, 1)), 2)
  expect_equal(mc_value(mc_job(1 + 1, 2)), 2)
  expect_warning(res <- mc_value(mc_job({ warning("jobwarn"); 5 }, 2)), "jobwarn")
  expect_equal(res, 5)
  expect_error(mc_value(mc_job(stop("jobboom"), 2)), "jobboom")
  job <- mc_job({ Sys.sleep(2); 1 }, 2)
  expect_null(mc_value(job, wait = FALSE))   # not done yet
  expect_equal(mc_value(job), 1)
})
