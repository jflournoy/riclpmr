#' a_random_intercept
#'
#' Create syntax for a single random intercept term based on a vector of
#' variable names and a variable name for the random intercept term. Used
#' internally by \code{\link{random_intercepts}}.
#'
#' @param vars Vector of variable names (column names in the target data.frame).
#' @param varname Used to name the latent random intercept.
#'
#' @return A list with the model text in \code{model_text}, random intercept
#'   name in \code{varname}, and \code{type} set to 'latent' (not used yet).
#'
#' @examples
a_random_intercept <- function(vars, varname){
  ri_name <- paste0('ri_', varname)
  right_side <- paste(paste0('1*', vars), collapse = ' + ')
  ri_text <- paste0(ri_name, ' =~ ', right_side)
  return(list(model_text = ri_text, varname = ri_name, type = 'latent'))
}

#' random_intercepts
#'
#' Creates random intercept model text for all sets of variables passed via
#' \code{var_groups}.
#'
#' @param var_groups A list of vectors of variable names. Each element of the
#'   list will be used to generate model text for a latent random intercept.
#'
#' @return A matrix with columns for each top-level element of \code{var_groups}
#'   and rows for each element of the list returned by
#'   \code{\link{a_random_intercept}}.
#'
#' @examples
random_intercepts <- function(var_groups){
  return(mapply(a_random_intercept,
                vars = var_groups,
                varname = names(var_groups)))
}

#' simple_vars
#'
#' Creates model text that covaries each element of \code{terms1} with each
#' element of \{terms2}.
#'
#' @param terms1 A vector of variable names.
#' @param terms2 A vector of variable names.
#' @param constrain Set all variances to equality. If not \code{FALSE}, provide
#'   the name for the equality constraint.
#'
#' @return Model text specifying covariances between terms.
#' @export
#'
simple_vars <- function(terms1, terms2, constrain = NULL){
  if(!is.null(constrain)){
    terms2 <- paste0(constrain, '*', terms2)
  }
  vars <- paste(paste0(terms1, ' ~~ ', terms2), collapse = '\n')
  return(vars)
}

#' residual_variances
#'
#' Creat model text specifying residual variances, possibly with constraints to
#' eqaulity across waves.
#'
#' @param var_groups A list of vectors of variable names. Each element of the
#'   list will be used to generate model text for residual variances, possibly
#'   with equality constraints across waves.
#'
#' @return
#' @export
#'
#' @examples
residual_variances <- function(var_groups, constrain = TRUE){
  residual_text <- mapply(function(vars, varname){
    if(constrain){
      first_res_vars <- simple_vars(vars[1], vars[1], constrain = NULL)
      const_res_vars <- simple_vars(vars[-1], vars[-1], constrain = paste0('e_', varname))
      res_vars <- paste(c(first_res_vars, const_res_vars), collapse = '\n')
    } else {
      res_vars <- simple_vars(vars, vars)
    }
    return(res_vars)
  }, var_groups, names(var_groups))

  return(paste(unlist(residual_text), collapse = '\n'))
}

#' varcovars
#'
#' Creates model text that specifies variances for each variable in
#' \code{terms}, as well as the full covariance matrix between those terms.
#'
#' @param terms Vector of variable names to create variance-covariance
#'   specification from.
#'
#' @return Model text specifying covariance between terms in \code{terms}.
#' @export
#'
varcovars <- function(terms){
  vars <- simple_vars(terms, terms)
  nterms <- length(terms)
  covars <- lapply(1:(nterms-1), function(i){
    right_side <- paste(terms[(i+1):nterms], collapse = ' + ')
    return(paste0(terms[i], ' ~~ ', right_side))
  })
  return(paste(c(vars, covars), collapse = '\n'))
}

#' manifest_ints
#'
#' Create model text that specifies intercepts for each term.
#'
#' @param terms Vector of variable names to create intercept specifications for.
#' @param constraint Whether to constrain the intercepts to a certain value. \code{'free'} is the default. Provide the specific value for the contraint.
#'
#' @return Model text for specifying intercepts
#' @export
#'
manifest_ints <- function(terms, constrain = 'free'){
  if(constrain != 'free'){
    mus <- paste0(constrain, '*1')
  } else {
    mus <- paste0(terms, '_mu*1')
  }
  ints <- paste0(terms, ' ~ ', mus)
  return(paste(ints, collapse = '\n'))
}

#' a_latent_resid_set
#'
#' Creates latent residual terms for each variable in \code{vars}, using
#' \code{varname} to name them.
#'
#' @param vars Vector of variable names (column names in the target data.frame).
#' @param varname Used to name the latent residuals.
#'
#' @return A list with the model text in \code{model_text}, latent residual
#'   names in \code{varname}, and \code{type} set to 'latent' (not used yet).
#' @export
#'
a_latent_resid_set <- function(vars, varname){
  nvars <- length(vars)
  right_side <- paste0('1*', vars)
  latent_names <- paste0('lat_', varname, 1:nvars)
  latent_resids <- paste0(latent_names, ' =~ ', right_side)
  return(list(model_text = latent_resids,
              varname = latent_names,
              type = 'latent'))

}

#' latent_resids
#'
#' Creates latent residuals model text for all sets of variables passed via
#' \code{var_groups}.
#'
#' @param var_groups A list of vectors of variable names. Each element of the
#'   list will be used to generate model text for a set of latent residuals
#'   based on the variables in that element.
#'
#' @return A matrix with columns for each top-level element of \code{var_groups}
#'   and rows for each element of the list returned by
#'   \code{\link{a_latent_resid_set}}.
#' @export
#'
#' @examples
#' some_vars <- list(x = c('x1', 'x2'), y = c('y1', 'y2'))
#' latent_resids(some_vars)
latent_resids <- function(var_groups){
  return(mapply(a_latent_resid_set,
                vars = var_groups,
                varname = names(var_groups)))
}

#' clpm
#'
#' Create model text specifying the cross-lag panel model portion of the
#' RI-CLPM, including all autoregressive and cross-lag effects. Set
#' \code{constrain} to \code{TRUE} to constrain these to be equal across all
#' waves.
#'
#' @param lat_resid_vars A list of vectors of variable names of the latent
#'   residual variables. Each element should be a vector of variable names that
#'   corresponds to the same construct measured repeatedly over waves.
#' @param constrain Constrain the path coefficiecients to be the same from
#'   wave-to-wave?
#'
#' @return Model text specifying the CLPM portion of the model.
#' @export
#'
clpm <- function(lat_resid_vars, constrain = TRUE){
  adf <- as.data.frame(lat_resid_vars)
  nsets <- dim(adf)[2]
  nvars <- dim(adf)[1]
  set_names <- names(adf)
  cross_lags <- lapply(1:nsets, function(i){
    left_side <- adf[-1, i]
    left_var <- set_names[i]
    if(constrain){
      right_side_adf <- data.frame(mapply( function(col, aname){
        constraint_name <- paste0(left_var, '_', aname)
        paste0(constraint_name, '*', col)
      }, adf[-nvars,], set_names))
    } else {
      right_side_adf <- adf
    }
    right_side <- apply(right_side_adf, 1, paste, collapse = ' + ')
    paste0(left_side, ' ~ ', right_side)
  })
  contemps <- contemp_covars(lat_resid_vars, constrain = constrain)
  return(paste(c(contemps, unlist(cross_lags)), collapse = '\n'))
}

#' contemp_covars
#'
#' Create model text to specify contemporaneous covariances between latent
#' residuals.
#'
#' @param lat_resid_vars A list of vectors of latent residual variable names.
#' @param constrain Constrain the covariances to be the same across waves?
#'
#' @return Model text specifying the contemporaneous covariances.
#' @export
#'
contemp_covars <- function(lat_resid_vars, constrain = TRUE){
  nsets <- length(lat_resid_vars)
  covars <- lapply(1:(nsets-1), function(i){
    left_side <- lat_resid_vars[[i]]
    if(constrain){
      lname <- names(lat_resid_vars)[[i]]
      rnames <- names(lat_resid_vars)[(i+1):nsets]
      constraint_names <- paste0('r_', lname, rnames)
      lat_resid_df <- data.frame(lat_resid_vars[(i+1):nsets], stringsAsFactors = FALSE)
      lat_resid_df_top_row <- lat_resid_df[1,]
      lat_resid_df_the_rest <- lat_resid_df[-1,]
      right_side_adf <- data.frame(mapply( function(col, cname){
        paste0(cname, '*', col)
      }, lat_resid_df_the_rest, constraint_names), stringsAsFactors = FALSE)
      right_side_adf <- rbind(lat_resid_df_top_row, right_side_adf)
    } else {
      right_side_adf <- data.frame(lat_resid_vars[(i+1):nsets])
    }
    right_side <- apply(right_side_adf, 1, paste, collapse = ' + ')
    return(paste0(left_side, ' ~~ ', right_side))
  })
  return(paste(unlist(covars), collapse = '\n'))
}

#' riclpm_text
#'
#' Given a list of sets of variables sampled repeatedly over equal time intervals, create a Random Intercept Cross Lagged Panel Model (RI-CLPM) specification for lavaan.
#'
#' @param var_groups A list of vectors of variable names.
#' @param constrain_over_waves Constrain regression coefficients, covariances, and residuals to be the same from wave to wave? Will not constrain variances and covariances of wave-1 latent residuals.
#' @param constrain_ints Constrain intercepts of manifest variables? Default is to free them. At the moment, passing other values is nonsense.
#'
#' @return Model text to be passed to a lavaan function.
#' @export
#'
riclpm_text <- function(var_groups, constrain_over_waves = TRUE, constrain_ints = 'free'){
  ri_terms <- random_intercepts(var_groups = var_groups)
  ri_text <- paste(ri_terms['model_text',], collapse = '\n')
  ri_varcov_text <- varcovars(ri_terms['varname',])
  man_ints_text <- manifest_ints(unlist(var_groups), constrain = constrain_ints)
  lat_resid_terms <- latent_resids(var_groups = var_groups)
  lat_resid_text <- paste0(unlist(lat_resid_terms['model_text',]), collapse = '\n')
  lat_resid_vars <- lat_resid_terms['varname',]
  clpm_text <- clpm(lat_resid_vars, constrain = constrain_over_waves)
  resid_vars_text <- residual_variances(lat_resid_vars, constrain = constrain_over_waves)

  lavmod <- paste(c(ri_text,
                    ri_varcov_text,
                    man_ints_text,
                    lat_resid_text,
                    clpm_text,
                    resid_vars_text),
                  collapse = '\n')
  return(lavmod)
}

#' lavriclpm
#'
#' This passes through a model and data to the lavaan::lavaan function with
#' proper constraints.
#'
#' @param riclpmModel The lavaan syntax for the RI-CLPM model
#' @param data The data
#' @param ... Other parameters passed to \code{\link[lavaan]{lavaan}}
#'
#' @return a fitted model
#' @export
#'
lavriclpm <- function(riclpmModel, data, ...){
  fit <- lavaan::lavaan(riclpmModel, data = data,
                        int.ov.free = F,
                        int.lv.free = F,
                        auto.fix.first = F,
                        auto.fix.single = F,
                        auto.cov.lv.x = F,
                        auto.cov.y = F,
                        auto.var = F, ...)
  return(fit)
}
