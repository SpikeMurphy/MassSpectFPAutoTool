test_that("`processing_check_file` works as intended.", {
  # create correct temp file
  tmp_dir <- tempdir()
  tmp_file <- file.path(tmp_dir, "example.mzXML")
  file.create(tmp_file)
  
  expect_no_warning(processing_check_file(FILE = tmp_file))
  expect_no_error(processing_check_file(FILE = tmp_file))

  expect_error(processing_check_file(FILE = "~/dir/downloads/example.mzXML"))
  expect_error(processing_check_file(FILE = "example.mzXML"))
  expect_error(processing_check_file(FILE = "string"))
  expect_error(processing_check_file(FILE = ""))
  expect_error(processing_check_file(FILE = NA))
  expect_error(processing_check_file(FILE = 0))
  
  # create temp file with wrong extension
  tmp_dir <- tempdir()
  tmp_file <- file.path(tmp_dir, "example.xml")
  file.create(tmp_file)
  
  expect_error(processing_check_file(FILE = tmp_file))
})


test_that("`processing_check_cleavage` works as intended.", {
  # TODO
  expect_no_warning(processing_check_cleavage(CLEAVAGE = "Trypsin"))
  expect_no_error(processing_check_cleavage(CLEAVAGE = "Trypsin"))
  
  expect_warning(processing_check_cleavage(CLEAVAGE = ""))
  expect_warning(processing_check_cleavage(CLEAVAGE = c("Trypsin", "Trypsin")))
  
  expect_warning(processing_check_cleavage(CLEAVAGE = NA))
  expect_warning(processing_check_cleavage(CLEAVAGE = 5))
  expect_warning(processing_check_cleavage(CLEAVAGE = "string"))
})


test_that("`processing_check_contaminants` works as intended.", {
  # TODO
  expect_no_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = c()))
  expect_no_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE =c(1, 2, 3)))
  expect_no_error(processing_check_contaminants(KERATIN = FALSE, TAG = "RFP", EXCLUDE = c()))
  expect_no_error(processing_check_contaminants(KERATIN = FALSE, TAG = "RFP", EXCLUDE =c(1, 2, 3)))
  
  expect_warning(processing_check_contaminants(KERATIN = c(1, 2, 3), TAG = "GFP", EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = 5, TAG = "GFP", EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = NA, TAG = "GFP", EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = "KERATIN", TAG = "GFP", EXCLUDE = c()))
  
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "", EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = c("GFP", "RFP"), EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = NA, EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = TRUE, EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = FALSE, EXCLUDE = c()))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = 5, EXCLUDE = c()))
  
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = ""))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = "string"))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = NA))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = TRUE))
  expect_warning(processing_check_contaminants(KERATIN = TRUE, TAG = "GFP", EXCLUDE = FALSE))
})


test_that("`processing_check_prerequisites` works as intended.", {
  # TODO
  expect_no_warning(processing_check_prerequisites(SNR, PEAKS))
  expect_no_error(processing_check_prerequisites(SNR, PEAKS))
})


test_that("`processing_check_output` works as intended.", {
  # TODO
  expect_no_warning(processing_check_plots(PLOTS))
  expect_no_error(processing_check_plots(PLOTS))
})


test_that("`processing_check_processing` works as intended.", {
  # TODO
  expect_no_warning(processing_check_processing(TRANSFORM, SMOOTH, BASELINE, CALIBRATE))
  expect_no_error(processing_check_processing(TRANSFORM, SMOOTH, BASELINE, CALIBRATE))
})
