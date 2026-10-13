riclpm_test_groups <- list(x = paste0('x', 1:4), y = paste0('y', 1:4))

riclpm_test_population <- function(){
  waves <- 1:4
  lines <- c(
    paste0('ri_x =~ ', paste0('1*x', waves, collapse = ' + ')),
    paste0('ri_y =~ ', paste0('1*y', waves, collapse = ' + ')),
    'ri_x ~~ 1*ri_x',
    'ri_y ~~ 1*ri_y',
    'ri_x ~~ .3*ri_y',
    paste0('x', waves, ' ~ 0*1'),
    paste0('y', waves, ' ~ 0*1'),
    paste0('lat_x', waves, ' =~ 1*x', waves),
    paste0('lat_y', waves, ' =~ 1*y', waves),
    'lat_x1 ~~ 1*lat_x1',
    'lat_y1 ~~ 1*lat_y1',
    'lat_x1 ~~ .2*lat_y1',
    paste0('lat_x', waves[-1], ' ~ .4*lat_x', waves[-4], ' + .1*lat_y', waves[-4]),
    paste0('lat_y', waves[-1], ' ~ .2*lat_x', waves[-4], ' + .5*lat_y', waves[-4]),
    paste0('lat_x', waves[-1], ' ~~ .7*lat_x', waves[-1]),
    paste0('lat_y', waves[-1], ' ~~ .7*lat_y', waves[-1]),
    paste0('lat_x', waves[-1], ' ~~ .1*lat_y', waves[-1])
  )
  return(paste(lines, collapse = '\n'))
}

riclpm_test_data <- function(){
  return(simriclpm(riclpm_test_population(), sample.nobs = 300, seed = 4417))
}
