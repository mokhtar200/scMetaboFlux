#' @title Package Attach
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("[scMetaboFlux] Single-Cell Energy Metabolism Inference Engine v1.0.0")
  packageStartupMessage("[scMetaboFlux] For citation, run: citation(\"scMetaboFlux\")")
}
