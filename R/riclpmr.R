#' riclpmr: A package for creating RI-CLPM syntax for lavaan
#'
#' This riclpmr package is a set of tools for creating lavaan-style model
#' specifications of Random Intercept Cross Lag Panel Models for arbitrarily
#' large numbers of constructs and measurement occations.
#'
#' It also contains a helpful function for passing this model specification to the appropriate \code{lavaan} or \code{blavaan} function.
#'
#' @docType package
#' @name riclpmr
#' @import lavaan
NULL

# lavaan options used by every function that hands RI-CLPM syntax to lavaan
# (lavriclpm, riclpm_partable, simriclpm). Kept in one place so a parameter
# table from riclpm_partable() is always built the same way lavriclpm() would
# build it from syntax.
.riclpm_lavaan_opts <- list(int.ov.free = FALSE,
                            int.lv.free = FALSE,
                            auto.fix.first = FALSE,
                            auto.fix.single = FALSE,
                            auto.cov.lv.x = FALSE,
                            auto.cov.y = FALSE,
                            auto.var = FALSE)
