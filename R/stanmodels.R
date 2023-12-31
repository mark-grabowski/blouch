# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("blouchBM_v1", "blouchOU1_v1", "blouchOUPredict_v1", "blouchOUReg_v1", "blouchOU_v1")

# load each stan module
Rcpp::loadModule("stan_fit4blouchBM_v1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4blouchOU1_v1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4blouchOUPredict_v1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4blouchOUReg_v1_mod", what = TRUE)
Rcpp::loadModule("stan_fit4blouchOU_v1_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})
