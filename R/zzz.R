.onLoad <- function(libname, pkgname) {
  if (length(names(stanmodels)) < 1) {
    return(invisible())
  }
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) Rcpp::loadModule(m, what = TRUE)
}
