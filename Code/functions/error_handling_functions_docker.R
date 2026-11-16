# ============================================================================
# Error Handling Functions for Cell-Hub
# ============================================================================
# Shared error handling for all server modules (single, multiple, spatial,
# CellChat, trajectory).
#
# Why: in Shiny, an uncaught error inside observe() / observeEvent() closes the
# whole session (grey screen) and the user loses the analysis in progress.
# Shiny offers no global handler that keeps the session alive, so every
# observer body is wrapped in run_safely(), which turns any error into a
# visible notification instead.
#
# Rules:
# - No reactive code in this file (no reactiveVal, observe, output$, session$):
#   observers stay in the *_server.R files.
# - req() and validate() keep their usual silent behaviour.
# - Coverage is checked by tests/check_error_handling.R:
#     Rscript ../tests/check_error_handling.R         # audit
#     Rscript ../tests/check_error_handling.R --fix   # wraps unprotected observers
#
# Author: Gaspard
# Last updated: 2026-10-21
# ============================================================================

# Run an observer body so that an error never closes the session: req()/validate() stops stay silent,
# any other error closes the progress modal and is shown to the user (persistent for user actions)
run_safely <- function(context, expr, persistent = TRUE) {
  tryCatch(expr, error = function(e) {
    if (inherits(e, c("shiny.silent.error", "validation"))) return(invisible(NULL))
    msg <- conditionMessage(e)[1]
    message(sprintf("[%s] ERROR: %s", context, msg))
    removeModal()
    showNotification(sprintf("%s failed: %s", context, msg), type = "error", duration = if (persistent) NULL else 10)
    invisible(NULL)
  })
}