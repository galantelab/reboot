############################################# TEST MODULE I - REGRESSION #############################################

test_that("rebootRegression runs on toy data", {
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")

  result <- rebootRegression(
    filein = example_file,
    bootstrap = 4,
    groupsize = 3,
    type = "transcript",
    force = TRUE,
    ncores = 1,
    seed = 123
  )

  expect_s3_class(result, "reboot_signature")
  expect_true(is.data.frame(result$signature))
  expect_true(nrow(result$signature) > 0)
  expect_true(all(c("feature", "coefficient") %in% colnames(result$signature)))

  expect_equal(result$params$ncores, 1)
  expect_equal(result$params$seed, 123)
})


test_that("rebootRegression is reproducible with the same seed", {
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")

  result1 <- rebootRegression(
    filein = example_file,
    bootstrap = 4,
    groupsize = 3,
    type = "transcript",
    force = TRUE,
    ncores = 1,
    seed = 123
  )

  result2 <- rebootRegression(
    filein = example_file,
    bootstrap = 4,
    groupsize = 3,
    type = "transcript",
    force = TRUE,
    ncores = 1,
    seed = 123
  )

  expect_identical(result1$signature, result2$signature)
})


test_that("rebootRegression gives identical results in serial and parallel modes", {
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")

  result_serial <- rebootRegression(
    filein = example_file,
    bootstrap = 4,
    groupsize = 3,
    type = "transcript",
    force = TRUE,
    ncores = 1,
    seed = 123
  )

  result_parallel <- rebootRegression(
    filein = example_file,
    bootstrap = 4,
    groupsize = 3,
    type = "transcript",
    force = TRUE,
    ncores = 2,
    seed = 123
  )

  expect_identical(result_serial$signature, result_parallel$signature)
})


test_that("rebootRegression validates ncores and seed", {
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")

  expect_error(
    rebootRegression(
      filein = example_file,
      bootstrap = 4,
      groupsize = 3,
      type = "transcript",
      force = TRUE,
      ncores = 0
    )
  )

  expect_error(
    rebootRegression(
      filein = example_file,
      bootstrap = 4,
      groupsize = 3,
      type = "transcript",
      force = TRUE,
      ncores = 1.5
    )
  )

  expect_error(
    rebootRegression(
      filein = example_file,
      bootstrap = 4,
      groupsize = 3,
      type = "transcript",
      force = TRUE,
      seed = 1.5
    )
  )
})

#######################################################################################################################

################################### TEST COMPLETE WORKFLOW - REGRESSION + SURVIVAL ###################################
test_that("rebootComplete runs on toy data", {
  # Temporary output directory
  tmp_dir <- tempdir()
  message("Using temporary directory: ", tmp_dir)
  outprefix <- file.path(tmp_dir, "rebootComplete")
  # outprefix <- "rebootComplete"

  # Loads toy datasets
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")
  clinical_file <- system.file("extdata", "toy_clinics.tsv", package = "Reboot")

  # Runs Reboot modules I + II - regression + survival
  result <- rebootComplete(
    filein = example_file,
    outprefix = outprefix,
    bootstrap = 10,
    groupsize = 10,
    percentagefilter = 0.3,
    variancefilter = 0.01,
    followup = NULL,
    type = "transcript",
    multivariate = TRUE,
    clinin = clinical_file,
    roc = TRUE,
    p.cutoff = 0.2,
    force = TRUE,
    plots = TRUE,
    table = TRUE,
    saveJSON = TRUE,
    saveRDS = TRUE,
    report = FALSE,
    log = TRUE,
    ncores = 1,
    seed = 123
  )

  # Main object
  expect_s3_class(result, "reboot_complete")

  # Internal modules
  expect_s3_class(result$regression, "reboot_signature")
  expect_s3_class(result$survival, "reboot_survival")

  # Regression outputs
  expect_true(is.data.frame(result$regression$signature))
  expect_true(nrow(result$regression$signature) > 0)

  # Survival outputs
  expect_true(is.list(result$survival$survival))
  expect_true(inherits(result$survival$survival$univariate$model, "coxph"))

  # Parameters
  expect_true(is.list(result$params))
  expect_equal(result$params$bootstrap, 10)
  expect_equal(result$params$ncores, 1)
  expect_equal(result$params$seed, 123)
  expect_equal(result$params$variancefilter, 0.01)
  expect_true(result$params$multivariate)
  expect_true(result$params$roc)

  # Metadata and reproducibility
  expect_true(!is.null(result$call))
  expect_true(is.list(result$metadata))
  expect_true(!is.null(result$metadata$package_version))
})
#######################################################################################################################
