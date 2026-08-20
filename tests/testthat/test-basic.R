############################################# TEST MODULE I - REGRESSION #############################################
# test_that("rebootRegression runs on toy data", {
#   # Temporary output directory
#   tmp_dir <- tempdir()
#   message("Using temporary directory: ", tmp_dir)
#   outprefix <- file.path(tmp_dir, "rebootRegression")
#   
#   # Loads toy dataset
#   example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")
#   
#   # Runs Reboot module I - regression
#   result <- rebootRegression(
#     filein = example_file,
#     outprefix = outprefix,
#     bootstrap = 10,
#     groupsize = 10,
#     percentagefilter = 0.3,
#     variancefilter = 0.01,
#     followup = NULL,
#     type = "transcript",
#     force = TRUE,
#     plots = TRUE,
#     table = TRUE,
#     saveJSON = TRUE,
#     saveRDS = FALSE,
#     report = TRUE,
#     log = TRUE
#   )
#   
#   # Main object
#   expect_s3_class(result, "reboot_signature")
#   
#   # Signature
#   expect_true(is.data.frame(result$signature))
#   expect_true(nrow(result$signature) > 0)
#   expect_true(all(c("feature", "coefficient") %in% colnames(result$signature)))
#   
#   # Metadata and reproducibility
#   expect_true(!is.null(result$call))
#   expect_true(is.list(result$metadata))
#   expect_true(!is.null(result$metadata$package_version))
# })
#######################################################################################################################

############################################## TEST MODULE II - SURVIVAL ##############################################
# test_that("rebootSurvival runs on toy data", {
#   # Temporary output directory
#   tmp_dir <- tempdir()
#   message("Using temporary directory: ", tmp_dir)
#   outprefix <- file.path(tmp_dir, "rebootSurvival")
#   
#   # Loads toy datasets
#   example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")
#   signature_file <- system.file("extdata", "toy_signature.txt", package = "Reboot")
#   clinical_file <- system.file("extdata", "toy_clinics.tsv", package = "Reboot")
#   
#   # Runs Reboot module II - survival
#   result <- rebootSurvival(
#     filein = example_file,
#     signature = signature_file,
#     outprefix = outprefix,
#     multivariate = TRUE,
#     clinin = clinical_file,
#     roc = TRUE,
#     variancefilter = 0.01,
#     followup = NULL,
#     p.cutoff = 0.2,
#     bootstrap = 10,
#     force = TRUE,
#     plots = TRUE,
#     table = TRUE,
#     saveJSON = TRUE,
#     saveRDS = FALSE,
#     report = TRUE,
#     log = TRUE
#   )
#   
#   # Main object
#   expect_s3_class(result, "reboot_survival")
#   
#   # Core components
#   expect_true(is.list(result$survival))
#   expect_true(is.data.frame(result$expression_data))
#   expect_true(is.data.frame(result$signature))
#   expect_true(is.data.frame(result$clinical_data))
#   
#   # Survival analyses
#   expect_true("univariate" %in% names(result$survival))
#   expect_true(inherits(result$survival$univariate$model, "coxph"))
#   
#   # ROC
#   expect_true(!is.null(result$roc_result))
#   expect_true(!is.null(result$ph_result))
#   
#   # Multivariate model
#   expect_true(!is.null(result$survival$multivariate))
#   
#   # Parameters
#   expect_true(result$params$multivariate)
#   expect_true(result$params$roc)
#   
#   # Metadata and reproducibility
#   expect_true(!is.null(result$call))
#   expect_true(is.list(result$metadata))
#   expect_true(!is.null(result$metadata$package_version))
# })
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
    report = TRUE,
    log = TRUE
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
  expect_equal(result$params$variancefilter, 0.01)
  expect_true(result$params$multivariate)
  expect_true(result$params$roc)

  # Metadata and reproducibility
  expect_true(!is.null(result$call))
  expect_true(is.list(result$metadata))
  expect_true(!is.null(result$metadata$package_version))
})
#######################################################################################################################

################################################# TEST SEED REPRODUCIBILITY ############################################
test_that("rebootRegression is reproducible for a fixed seed and varies otherwise", {
  example_file <- system.file("extdata", "toy_expression.tsv", package = "Reboot")

  # Two runs with the same explicit seed must yield an identical signature
  result_a <- rebootRegression(filein = example_file, bootstrap = 5, groupsize = 10,
                               type = "transcript", force = TRUE, seed = 123)
  result_b <- rebootRegression(filein = example_file, bootstrap = 5, groupsize = 10,
                               type = "transcript", force = TRUE, seed = 123)
  expect_identical(result_a$signature, result_b$signature)

  # The seed used is recorded in params for reproducibility bookkeeping
  expect_equal(result_a$params$seed, 123)

  # A different seed is not guaranteed (and, for this stochastic pipeline, not expected)
  # to produce the exact same signature
  result_c <- rebootRegression(filein = example_file, bootstrap = 5, groupsize = 10,
                               type = "transcript", force = TRUE, seed = 456)
  expect_false(isTRUE(identical(result_a$signature, result_c$signature)))

  # Setting a seed inside the pipeline must not leak into the caller's global RNG state
  old_seed <- .GlobalEnv$.Random.seed
  invisible(rebootRegression(filein = example_file, bootstrap = 5, groupsize = 10,
                             type = "transcript", force = TRUE, seed = 123))
  expect_identical(.GlobalEnv$.Random.seed, old_seed)
})
#######################################################################################################################
