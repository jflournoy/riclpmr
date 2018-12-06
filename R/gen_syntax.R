#' make_riclpm_lavaan
#'
#' @param xname name of x variable
#' @param yname name of y variable
#' @param waveints integers specifying wave labels
#' @param constrain_all set all constraints via labels
#' @param alpha set label or value for alpha constraint
#' @param delta set label or value for delta constraint
#' @param beta set label or value for beta constraint
#' @param gamma set label or value for gamma constraint
#' @param uu set label or value for uu constraint
#' @param vv set label or value for vv constraint
#' @param uv set label or value for uv constraint
#' @param kk set label or value for kk constraint
#' @param oo set label or value for oo constraint
#' @param ko set label or value for ko constraint
#' @param pp set label or value for pp constraint
#' @param qq set label or value for qq constraint
#' @param pq set label or value for pq constraint
#'
#' @return A character string that is syntax for an RI-CLPM model as specified
#' @export
#'
make_riclpm_lavaan <- function(xname = 'x', yname = 'y', waveints, constrain_all = F,
                               alpha = '', delta = '', beta = '', gamma = '', uu = '', vv = '', uv = '',
                               kk = '', oo = '' , ko = '', pp = '', qq = '', pq = ''){
  waveids <- waveints
  if(constrain_all){
    b <- list(alpha = 'alpha',
              delta = 'delta',
              beta = 'beta',
              gamma = 'gamma',
              uu = 'uu',
              vv = 'vv',
              uv = 'uv',
              kk = 'kk',
              oo = 'oo',
              ko = 'ko',
              pp = 'pp',
              qq = 'qq',
              pq = 'pq')
  }
  b_ <- list(alpha = alpha,
             delta = delta,
             beta = beta,
             gamma = gamma,
             uu = uu,
             vv = vv,
             uv = uv,
             kk = kk,
             oo = oo,
             ko = ko,
             pp = pp,
             qq = qq,
             pq = pq)

  for(b_i in 1:length(b)) {
    if(b_[[b_i]] != ''){
      b[[b_i]] <- b_[[b_i]]
    }
    if(b[[b_i]] != ''){
      b[[b_i]] = paste0(b[[b_i]], '*')
    }
  }

  latent_means <- paste(paste0(c('kappa =~ ', 'omega =~ '),
                               lapply(lapply(c(xname, yname), paste0, waveids),
                                      function(vars){
                                        paste(paste0('1*', vars), collapse = ' + ')
                                      })),
                        collapse = '\n')

  intercepts <- paste(unlist(lapply(list(c(xname,'mu'), c(yname, 'pi')),
                                    function(pair){
                                      lapply(waveids, function(wave){
                                        paste0(pair[1], wave, ' ~ ', pair[2], wave, '*1')
                                      })
                                    })), collapse = '\n')

  latent_covar <- paste0('
                         kappa ~~ ',b$kk,'kappa #variance
                         omega ~~ ',b$oo,'omega #variance
                         kappa ~~ ',b$ko,'omega #covariance')

  latent_resids <- paste(unlist(lapply(list(c('p',xname), c('q', yname)),
                                       function(pair){
                                         lapply(waveids, function(wave){
                                           paste0(pair[1], wave, ' =~ ', '1*', pair[2], wave)
                                         })
                                       })), collapse = '\n')

  regressions <- paste(
    unlist(lapply(list(c('p','q', b$alpha, b$beta), c('q', 'p', b$delta, b$gamma)),
                  function(pair){
                    lapply(rev(waveids[-1]), function(wave){
                      paste0(pair[1], wave, ' ~ ', pair[3], pair[1], wave-1, ' + ', pair[4], pair[2], wave-1)
                    })
                  })), collapse = '\n')

  first_wave_varcovar <- paste0('
p', waveids[1], ' ~~ ', b$pp,'p', waveids[1], '
q', waveids[1], ' ~~ ', b$qq,'q', waveids[1], '
p', waveids[1], ' ~~ ', b$pq,'q', waveids[1])

  residvar <- paste(unlist(lapply(list(c('p', b$uu), c('q', b$vv)),
                                  function(avar){
                                    lapply(waveids[-1], function(wave){
                                      paste0(avar[1], wave, ' ~~ ', avar[2], avar[1], wave)
                                    })
                                  })), collapse = '\n')

  contemporaneous <- paste(lapply(waveids[-1],
                                  function(wave){
                                    paste0('p', wave, ' ~~ ', b$uv, 'q', wave)
                                  }),
                           collapse = '\n')


  generatingModel <- paste(latent_means,
                           intercepts,
                           latent_covar,
                           latent_resids,
                           regressions,
                           first_wave_varcovar,
                           residvar,
                           contemporaneous, sep = '\n')
  return(generatingModel)
}


#' lavriclpm
#'
#' This passes through a model and data to the lavaan::lavaan
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
