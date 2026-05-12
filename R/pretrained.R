# classydms: download and cache pretrained encoder files
# Uses the {piggyback} package to fetch large .pt files from GitHub
# Releases on demand, since they are too large to include in the package
# itself.

#' Storage path for cached pretrained encoders
#'
#' Returns the path inside the installed package where downloaded
#' encoders are cached. Created on first use.
#'
#' @keywords internal
pretrained_cache_dir <- function() {
  cache <- file.path(tools::R_user_dir("classydms", "data"), "pretrained")
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  cache
}

#' List available pre-trained encoder files
#'
#' Queries the GitHub Releases of `mljaniczek/classydms` and returns the
#' list of `.pt` and `.Rdata` files attached as release assets.
#'
#' @param repo GitHub repo in `"owner/name"` form (default "mljaniczek/classydms").
#' @param tag Release tag to query (default "pretrained-v0"); use NULL to
#'   list all releases.
#' @return Character vector of available file names.
#' @export
list_pretrained_encoders <- function(repo = "mljaniczek/classydms",
                                       tag = "pretrained-v0") {
  if (!requireNamespace("piggyback", quietly = TRUE)) {
    stop("Package 'piggyback' is required. Install with install.packages('piggyback').")
  }
  tryCatch({
    df <- piggyback::pb_list(repo = repo, tag = tag)
    df$file_name
  }, error = function(e) {
    message("Could not query releases for ", repo,
            " (tag=", tag, "): ", conditionMessage(e))
    character(0)
  })
}

#' Load a pre-trained encoder by name
#'
#' Downloads the requested encoder file from GitHub Releases (cached
#' locally) and returns the loaded torch module. Pass `with_autoencoder
#' = TRUE` to also fetch the full autoencoder (encoder + decoder), which
#' is needed for reconstruction QC but not for classification or saliency
#' mapping.
#'
#' @param name Base name of the encoder file (without `.pt` extension).
#' @param repo GitHub repo.
#' @param tag Release tag.
#' @param with_autoencoder If TRUE, also fetch the autoencoder file and
#'   return both as a list.
#' @param force_download Force re-download even if cached locally.
#' @return Either a torch encoder module (if `with_autoencoder = FALSE`)
#'   or a list with `encoder`, `autoencoder`, and `manifest` elements.
#' @examples
#' \dontrun{
#' enc <- load_pretrained_encoder("flu_v1")
#' bundle <- load_pretrained_encoder("flu_v1", with_autoencoder = TRUE)
#' bundle$encoder; bundle$autoencoder; bundle$manifest
#' }
#' @export
load_pretrained_encoder <- function(name,
                                      repo = "mljaniczek/classydms",
                                      tag = "pretrained-v0",
                                      with_autoencoder = FALSE,
                                      force_download = FALSE) {
  if (!requireNamespace("piggyback", quietly = TRUE)) {
    stop("Package 'piggyback' is required for downloading pretrained ",
         "encoders. Install with install.packages('piggyback').")
  }
  cache <- pretrained_cache_dir()
  encoder_file <- paste0(name, ".pt")
  encoder_path <- file.path(cache, encoder_file)
  if (force_download || !file.exists(encoder_path)) {
    message("Downloading ", encoder_file, " from ", repo, " (tag=", tag, ")...")
    piggyback::pb_download(file = encoder_file, repo = repo, tag = tag,
                            dest = cache, overwrite = TRUE)
  }
  encoder <- torch::torch_load(encoder_path)
  if (!with_autoencoder) return(encoder)

  ae_file <- paste0(name, "_autoencoder.pt")
  ae_path <- file.path(cache, ae_file)
  if (force_download || !file.exists(ae_path)) {
    message("Downloading ", ae_file, "...")
    tryCatch(
      piggyback::pb_download(file = ae_file, repo = repo, tag = tag,
                              dest = cache, overwrite = TRUE),
      error = function(e) {
        message("Autoencoder not available for ", name, "; returning NULL.")
      })
  }
  autoencoder <- if (file.exists(ae_path)) torch::torch_load(ae_path) else NULL

  manifest_file <- paste0(name, "_manifest.Rdata")
  manifest_path <- file.path(cache, manifest_file)
  if (force_download || !file.exists(manifest_path)) {
    tryCatch(
      piggyback::pb_download(file = manifest_file, repo = repo, tag = tag,
                              dest = cache, overwrite = TRUE),
      error = function(e) {
        message("Manifest not available for ", name, "; returning NULL.")
      })
  }
  manifest <- NULL
  if (file.exists(manifest_path)) {
    e <- new.env()
    load(manifest_path, envir = e)
    manifest <- e$training_manifest
  }
  list(encoder = encoder, autoencoder = autoencoder, manifest = manifest)
}
