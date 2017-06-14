.onLoad <- function(libname, pkgname) {
    stanmods <- names(stanmodels)
    if (length(stanmods != 0)) {
        modules <- paste0("stan_fit4", names(stanmodels), "_mod")
        for (m in modules) loadModule(m, what = TRUE)
    }
}
