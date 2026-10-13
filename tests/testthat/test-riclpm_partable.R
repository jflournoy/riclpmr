pe_cols <- c('lhs', 'op', 'rhs', 'label', 'est', 'se', 'pvalue')

fit_both_ways <- function(constrain){
  d <- riclpm_test_data()
  model_text <- riclpm_text(riclpm_test_groups, constrain_over_waves = constrain)
  pt <- riclpm_partable(model_text)
  fit_text <- lavriclpm(riclpmModel = model_text, data = d)
  fit_pt <- lavriclpm(riclpmModel = pt, data = d)
  return(list(text = fit_text, pt = fit_pt))
}

expect_identical_fits <- function(fits){
  expect_true(lavaan::lavInspect(fits$text, 'converged'))
  expect_true(lavaan::lavInspect(fits$pt, 'converged'))
  pe_text <- lavaan::parameterEstimates(fits$text)
  pe_pt <- lavaan::parameterEstimates(fits$pt)
  expect_gt(nrow(pe_text), 0)
  expect_identical(nrow(pe_pt), nrow(pe_text))
  for(col in pe_cols){
    expect_identical(pe_pt[[col]], pe_text[[col]], info = col)
  }
  expect_identical(lavaan::fitMeasures(fits$pt, c('chisq', 'df', 'npar')),
                   lavaan::fitMeasures(fits$text, c('chisq', 'df', 'npar')))
}

test_that('fit from riclpm_partable() is identical to fit from text, constrained', {
  expect_identical_fits(fit_both_ways(constrain = TRUE))
})

test_that('fit from riclpm_partable() is identical to fit from text, unconstrained', {
  expect_identical_fits(fit_both_ways(constrain = FALSE))
})

test_that('the identity comparison can detect a different model', {
  fits_c <- fit_both_ways(constrain = TRUE)
  fits_u <- fit_both_ways(constrain = FALSE)
  pe_c <- lavaan::parameterEstimates(fits_c$pt)
  pe_u <- lavaan::parameterEstimates(fits_u$pt)
  expect_false(identical(pe_c$est, pe_u$est))
  expect_false(identical(lavaan::fitMeasures(fits_c$pt, 'df'),
                         lavaan::fitMeasures(fits_u$pt, 'df')))
})

test_that('riclpm_partable() returns a parameter table built with lavriclpm options', {
  pt_c <- riclpm_partable(riclpm_text(riclpm_test_groups, constrain_over_waves = TRUE))
  pt_u <- riclpm_partable(riclpm_text(riclpm_test_groups, constrain_over_waves = FALSE))

  expect_true(all(c('lhs', 'op', 'rhs', 'free', 'ustart', 'label') %in% names(pt_c)))

  ri_row <- which(pt_c$lhs == 'ri_x' & pt_c$op == '=~' & pt_c$rhs == 'x1')
  expect_length(ri_row, 1)
  expect_identical(as.numeric(pt_c$free[ri_row]), 0)
  expect_identical(as.numeric(pt_c$ustart[ri_row]), 1)

  expect_true('x_y' %in% pt_c$label)
  expect_true('y_x' %in% pt_c$label)
  expect_false('x_y' %in% pt_u$label)
  expect_false('y_x' %in% pt_u$label)

  manifest_resid <- which(pt_c$lhs == 'x1' & pt_c$op == '~~' & pt_c$rhs == 'x1')
  expect_length(manifest_resid, 1)
  expect_identical(as.numeric(pt_c$free[manifest_resid]), 0)
  expect_identical(as.numeric(pt_c$ustart[manifest_resid]), 0)

  d <- riclpm_test_data()
  fit <- do.call(lavaan::lavaan,
                 c(list(model = pt_c, data = d), riclpmr:::.riclpm_lavaan_opts))
  expect_true(lavaan::lavInspect(fit, 'converged'))
})

test_that('lavriclpm() and riclpm_partable() share one option list', {
  opts <- riclpmr:::.riclpm_lavaan_opts
  expect_identical(sort(names(opts)),
                   sort(c('int.ov.free', 'int.lv.free', 'auto.fix.first',
                          'auto.fix.single', 'auto.cov.lv.x', 'auto.cov.y',
                          'auto.var')))
  expect_true(all(vapply(opts, isFALSE, logical(1))))
})
