#' @title Package Initialization
#' @description Package startup functions
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  packageStartupMessage("")  
  packageStartupMessage("  scMetaboFlux - Single-Cell Energy Metabolism Inference Engine")
  packageStartupMessage("  Version: 1.0.0")
  packageStartupMessage("")
}

#' @title Package Attach
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("[scMetaboFlux] For citation information, run: citation(\"scMetaboFlux\")")
}

#' @title Package Unload
#' @keywords internal
.onUnload <- function(libpath) {
  library.dynam.unload("scMetaboFlux", libpath)
}
