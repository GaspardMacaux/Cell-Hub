# Keeps a download button disabled until the result it exports exists. available() returns that result: NULL, FALSE, an empty
# table or list, a req() / validate() failure or any error mean "nothing to download yet". The observer re-runs whenever the
# result changes, so the button is enabled as soon as the analysis has run and disabled again when the result is cleared.
# The "disabled" class is set too: on a download link (<a class="btn">) it is the class, not the attribute, that blocks the click.
bind_download_state <- function(id, available) {
  observe({
    ok <- tryCatch({
      x <- available()
      !is.null(x) && !isFALSE(x) && !(is.data.frame(x) && nrow(x) == 0) && !(is.list(x) && !is.data.frame(x) && length(x) == 0)
    }, error = function(e) FALSE)
    shinyjs::toggleState(id, condition = ok)
    shinyjs::toggleCssClass(id, "disabled", condition = !ok)
  })
}
