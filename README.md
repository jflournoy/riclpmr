# riclpmr

This is an R package that should help you generate `lavaan` syntax for random intercept cross-lagged panel models (RI-CLPMs) for multivariate longitudinal designs. It is very much in alpha, so please check the output carefully and raise any issues. 

For documentation, such as it is, see [here](https://johnflournoy.science/riclpmr/). The main function is [`riclpm_text`](https://johnflournoy.science/riclpmr/reference/riclpm_text.html).

See [this tutorial](https://johnflournoy.science/riclpm-lavaan-demo.html) for more information about the models.

For a simulation that uses riclpmr to compare the RI-CLPM with the lagged multilevel models often used in its place, see [False positives in multilevel models for longitudinal data](https://johnflournoy.science/lagged_mlm/).

## Why the random intercept

The traditional cross-lagged panel model regresses each variable on every variable at the previous wave. When people differ in stable, trait-like ways, those stable differences load onto the autoregressive and cross-lagged paths, so a "lagged effect" of x on y can just be the fact that people who are generally high on x are also generally high on y (Hamaker et al., 2015). The RI-CLPM gives each construct a random intercept that absorbs those stable between-person differences, and puts the autoregressive and cross-lagged paths on the within-person deviations from them. The lagged paths then describe how a person's departure from their own usual level on one variable predicts their departure on another at the next wave.

Writing that syntax by hand gets long fast: k constructs over T waves need k random intercepts, k × T latent within-person deviations, and k² × (T − 1) lagged paths, plus the covariances among them. `riclpm_text()` generates it for any number of constructs and waves.

## What it generates

Given a named list with one vector of column names per construct, ordered by wave, `riclpm_text()` writes:

- a random intercept per construct, with every wave's loading fixed to 1, and the variances and covariances among the random intercepts;
- a freely estimated intercept for every observed variable;
- a latent within-person deviation for every observation, with its loading fixed to 1;
- autoregressive and cross-lagged paths among the deviations at adjacent waves, and covariances among the deviations within each wave.

By default (`constrain_over_waves = TRUE`) the lagged paths, the within-wave covariances, and the residual variances of the deviations are held equal from wave 2 on, while the wave-1 variances and covariances are left free. Those equality constraints assume the waves are equally spaced and the process is stationary. Set `constrain_over_waves = FALSE` to estimate them freely.

Observed-variable residual variances are not specified, so they are fixed to zero: each observation is fully decomposed into its random intercept and its within-person deviation, as in Hamaker et al. (2015). `lavriclpm()` fits the model with the `lavaan` options this needs (nothing is added automatically), and `simriclpm()` simulates data from the same syntax with `lavaan::simulateData()`. Any path you don't give a value to is simulated as zero, and in the default syntax that is every autoregressive and cross-lagged path, so edit the model text to set the values you want before simulating.

```r
remotes::install_github("jflournoy/riclpmr")

library(riclpmr)
var_groups <- list(
  x = c("x_t1", "x_t2", "x_t3"),
  y = c("y_t1", "y_t2", "y_t3"),
  z = c("z_t1", "z_t2", "z_t3"))
model_text <- riclpm_text(var_groups)
fake_data <- simriclpm(riclpmModel = model_text, sample.nobs = 500)
fit <- lavriclpm(riclpmModel = model_text, data = fake_data)
summary(fit)
```

`lavriclpm(..., blavaan = TRUE)` fits the same model with `blavaan` instead (install it separately).

## Status

Written in 2018–2020 and not developed since. In October 2026 the example above still ran and converged with lavaan 0.7.3. The package does not handle measurement models (multiple indicators per construct), unequal spacing between waves, or the later extensions of the RI-CLPM (Mulder & Hamaker, 2021); for those, write or adapt the syntax yourself.

## References

Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique of the cross-lagged panel model. *Psychological Methods, 20*(1), 102–116. https://doi.org/10.1037/a0038889

Mulder, J. D., & Hamaker, E. L. (2021). Three extensions of the random intercept cross-lagged panel model. *Structural Equation Modeling, 28*(4), 638–648. https://doi.org/10.1080/10705511.2020.1784738
