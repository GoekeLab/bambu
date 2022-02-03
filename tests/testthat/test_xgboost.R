context("XGBoost functionality")

# Comparison of previously saved JSON of model and JSON of current model is required 
# because an XGBoost model saved as JSON cannot be loaded properly back into an XGBoost object
test_that("fitXGBoostModel returns expected model for splice junction correction", {
  # Matrix representing 100,000 examples with 3 features
  data_train <- matrix(seq(1:300000), nrow=100000)
  # Test matrix representing 16,000 examples with 3 features
  data_test <- matrix(c(seq(1:28000), seq(280001:300000)), nrow=16000)
  labels_train <- c(rep(1,50000), rep(0,50000))
  
  xgb_model <- fitXGBoostModel(labels_train, data_train, show.cv=TRUE)
  # Compute the predictions
  xgb_predictions <- predict(xgb_model, data_test)
  
  # Dump the new model and compare it for equality with the previous JSON dump that is saved
  curr_model = xgb.dump(xgb_model, dump_format='json')

  # To save a new reference model and reference predictions. Can be removed.
  #xgb.dump(xgb_model, '~/Downloads/FYP/bambu_new/bambu/tests/testData/xgb_model_splice_junction_correction.txt', dump_format='json')
  #writeLines(as.character(xgb_predictions),'~/Downloads/FYP/bambu_new/bambu/tests/testData/xgb_predictions_splice_junction_correction.txt')


  # Load reference model and predictions
  # Have to load the reference where readLines will produce a list but need to join the elements back using \n to make it equivalent
  # to the original model
  xgb_model_ref <- paste(readLines(system.file('extdata','xgb_model_splice_junction_correction.txt', package='bambu')), collapse='\n')
  xgb_predictions_ref <- as.numeric(readLines(system.file('extdata','xgb_predictions_splice_junction_correction.txt', package='bambu')))
  expect_equal(curr_model, xgb_model_ref)
  expect_equal(xgb_predictions, xgb_predictions_ref)
}) 
