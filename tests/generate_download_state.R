# generate_download_state.R — lists every download button of the Cell-Hub servers and the result it needs, and writes the
# bind_download_state() lines that keep each button disabled until that result exists.
# The result is read from the code itself: the first req(...) of the download (content function of a downloadHandler, or the
# reactive given to createDownloadHandler); a createDownloadHandler fed with a reactive name uses that reactive.
# Only downloads registered directly in the server function are bound: an output defined inside an observer depends on that
# observer's local variables and is listed apart (left as it is).
# Usage: Rscript generate_download_state.R <app_dir> <out_dir>
# Output: <out_dir>/<server file>.download_state.R (block to paste before the closing brace of the server function) and
#         <out_dir>/download_inventory.csv (every download, its binding or the reason it has none)
args <- commandArgs(trailingOnly = TRUE)
app_dir <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
inventory <- list()
add <- function(file, id, kind, status, available = "", note = "") inventory[[length(inventory) + 1]] <<- data.frame(file = file, id = id, kind = kind, status = status,
                                                                                                                     available = available, note = note, stringsAsFactors = FALSE)
one_line <- function(e) paste(trimws(deparse(e, width.cutoff = 500L)), collapse = " ")

# First req(...) call found in an expression (depth-first, in code order); NULL when there is none
find_req <- function(e) {
  if (is.call(e)) {
    fn <- e[[1]]
    if (is.name(fn) && identical(as.character(fn), "req") || (is.call(fn) && identical(one_line(fn), "shiny::req"))) return(e)
    el <- as.list(e)
    for (i in seq_along(el)[-1]) {
      # Empty arguments (x[, 1]) are skipped without being bound to a variable
      if (is.null(el[[i]]) || is.pairlist(el[[i]]) || identical(el[[i]], quote(expr = ))) next
      r <- find_req(el[[i]]); if (!is.null(r)) return(r)
    }
  } else if (is.expression(e) || (is.list(e) && !is.pairlist(e))) {
    for (i in seq_along(e)) { if (identical(e[[i]], quote(expr = ))) next; r <- find_req(e[[i]]); if (!is.null(r)) return(r) }
  }
  NULL
}

# Output id of an assignment output$id <- ... / output[["id"]] <- ...; NULL otherwise
output_id <- function(lhs) {
  if (is.call(lhs) && identical(as.character(lhs[[1]]), "$") && identical(as.character(lhs[[2]]), "output")) return(as.character(lhs[[3]]))
  if (is.call(lhs) && identical(as.character(lhs[[1]]), "[[") && identical(as.character(lhs[[2]]), "output") && is.character(lhs[[3]])) return(lhs[[3]])
  NULL
}

# Names of the reactives defined in the server function (x <- reactiveVal / reactive / eventReactive / debounce / throttle)
reactive_names <- character(0)
reactive_bodies <- list()
collect_reactives <- function(top) {
  out <- character(0); bodies <- list()
  for (st in top) {
    if (is.call(st) && identical(as.character(st[[1]]), "<-") && length(st) == 3 && is.name(st[[2]]) && is.call(st[[3]])) {
      fn <- tail(as.character(st[[3]][[1]]), 1)
      if (fn %in% c("reactiveVal", "reactive", "eventReactive", "debounce", "throttle")) out <- c(out, as.character(st[[2]]))
      # Body of a computed reactive: its own first check is used instead of computing it (a session object or a plot is costly)
      if (fn == "reactive" && length(st[[3]]) >= 2) bodies[[as.character(st[[2]])]] <- st[[3]][[2]]
    }
  }
  reactive_bodies <<- bodies
  out
}

# First call x() of a server reactive found in an expression (depth-first, code order); NULL when there is none
first_reactive_call <- function(e) {
  if (is.call(e)) {
    if (is.name(e[[1]]) && as.character(e[[1]]) %in% reactive_names && length(e) == 1) return(e)
    el <- as.list(e)
    for (i in seq_along(el)[-1]) {
      if (is.null(el[[i]]) || is.pairlist(el[[i]]) || identical(el[[i]], quote(expr = ))) next
      r <- first_reactive_call(el[[i]]); if (!is.null(r)) return(r)
    }
  }
  NULL
}

# Calls never copied into an availability check (they write files, open dialogs or have other side effects)
side_effect_fns <- c("write.csv", "write.table", "ggsave", "saveRDS", "showModal", "removeModal", "showNotification", "png", "pdf", "dev.off",
                     "writeLines", "openxlsx::write.xlsx", "write.xlsx", "file.copy", "print", "message")
uses_side_effects <- function(e) {
  if (is.name(e)) return(identical(as.character(e), "file"))
  if (!is.call(e)) return(FALSE)
  if (one_line(e[[1]]) %in% side_effect_fns) return(TRUE)
  el <- as.list(e)
  any(vapply(seq_along(el), function(i) !is.null(el[[i]]) && !is.pairlist(el[[i]]) && !identical(el[[i]], quote(expr = )) && uses_side_effects(el[[i]]), logical(1)))
}

# Availability of a download body (content function body or reactive body):
# 1. statements up to the first top-level req(), when they have no side effect (local variables of the download are kept);
# 2. otherwise the first server reactive read by the body.
# validate(need(a, msg), need(b, msg)) rewritten as req(a, b)
as_req <- function(st) {
  if (!is.call(st)) return(NULL)
  fn <- tail(as.character(st[[1]]), 1)
  if (fn == "req") return(st)
  if (fn == "validate") {
    conds <- lapply(as.list(st)[-1], function(n) if (is.call(n) && tail(as.character(n[[1]]), 1) == "need") n[[2]] else NULL)
    conds <- Filter(Negate(is.null), conds)
    if (length(conds) > 0) return(as.call(c(as.name("req"), conds)))
  }
  NULL
}

body_availability <- function(b, resolve = TRUE) {
  stmts <- if (is.call(b) && identical(as.character(b[[1]]), "{")) as.list(b)[-1] else list(b)
  for (k in seq_along(stmts)) {
    st <- as_req(stmts[[k]])
    if (is.null(st)) next
    {
      prefix <- c(stmts[seq_len(k - 1)], list(st))
      if (!any(vapply(prefix, uses_side_effects, logical(1)))) {
        return(list(expr = if (k == 1) one_line(st) else paste0("{ ", paste(vapply(prefix, one_line, character(1)), collapse = "; "), " }"),
                    note = "statements up to the first req()"))
      }
      break
    }
  }
  r <- first_reactive_call(b)
  if (!is.null(r)) return(resolve_reactive(one_line(r), "first reactive read", resolve))
  list(expr = NULL, note = "no req() and no server reactive read")
}

# A check reduced to one computed reactive x() is replaced by that reactive's own first check (one level)
resolve_reactive <- function(expr, note, resolve = TRUE) {
  name <- sub("\\(\\)$", "", expr)
  if (resolve && grepl("^[A-Za-z0-9_.]+\\(\\)$", expr) && name %in% names(reactive_bodies)) {
    inner <- body_availability(reactive_bodies[[name]], resolve = FALSE)
    if (!is.null(inner$expr)) return(list(expr = inner$expr, note = paste0(note, ", checked through ", name, "'s own first check")))
  }
  list(expr = expr, note = note)
}

# Availability expression of one download call (downloadHandler or createDownloadHandler), or NULL with a note
availability <- function(call) {
  fn <- tail(as.character(call[[1]]), 1)
  args_call <- as.list(call)[-1]
  if (fn == "createDownloadHandler") {
    rd <- if ("reactive_data" %in% names(args_call)) args_call$reactive_data else args_call[[1]]
    if (is.name(rd)) return(resolve_reactive(paste0(as.character(rd), "()"), "reactive given to createDownloadHandler"))
    if (is.call(rd) && tail(as.character(rd[[1]]), 1) == "reactive") return(body_availability(rd[[2]]))
    return(list(expr = NULL, note = "reactive_data is neither a reactive name nor reactive({...})"))
  }
  content <- if ("content" %in% names(args_call)) args_call$content else if (length(args_call) >= 2) args_call[[2]] else NULL
  if (is.call(content) && identical(as.character(content[[1]]), "function")) return(body_availability(content[[3]]))
  list(expr = NULL, note = "content is not an inline function")
}

# Every output$x <- downloadHandler / createDownloadHandler anywhere in an expression (for the nested ones listed apart)
find_downloads <- function(e, acc = list()) {
  if (is.call(e)) {
    if (identical(as.character(e[[1]]), "<-") && length(e) == 3 && is.call(e[[3]])) {
      id <- output_id(e[[2]]); fn <- tail(as.character(e[[3]][[1]]), 1)
      if (!is.null(id) && fn %in% c("downloadHandler", "createDownloadHandler")) acc[[length(acc) + 1]] <- id
    }
    el <- as.list(e)
    for (i in seq_along(el)[-1]) {
      if (is.null(el[[i]]) || is.pairlist(el[[i]]) || identical(el[[i]], quote(expr = ))) next
      acc <- find_downloads(el[[i]], acc)
    }
  }
  acc
}

servers <- list.files(app_dir, pattern = "_server(_docker)?\\.R$", full.names = TRUE, recursive = FALSE)
servers <- servers[!grepl("^main_server", basename(servers))]
for (f in servers) {
  exprs <- parse(f, keep.source = FALSE)
  # Server function: <name> <- function(input, output, session) { ... }
  srv <- Filter(function(e) is.call(e) && identical(as.character(e[[1]]), "<-") && is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function"), as.list(exprs))
  if (length(srv) == 0) { add(basename(f), "", "", "no server function found"); next }
  body <- srv[[1]][[3]][[3]]
  top <- as.list(body)[-1]
  reactive_names <- collect_reactives(top)
  bound_ids <- character(0); lines <- character(0)
  for (st in top) {
    if (!(is.call(st) && identical(as.character(st[[1]]), "<-") && length(st) == 3 && is.call(st[[3]]))) next
    id <- output_id(st[[2]]); fn <- tail(as.character(st[[3]][[1]]), 1)
    if (is.null(id) || !fn %in% c("downloadHandler", "createDownloadHandler")) next
    bound_ids <- c(bound_ids, id)
    # Parameter exports and reports work at any time (they report what has been done so far)
    if (grepl("params|report", id)) { add(basename(f), id, fn, "always enabled", "", "parameter export / report"); next }
    a <- availability(st[[3]])
    if (is.null(a$expr)) { add(basename(f), id, fn, "NOT BOUND", "", a$note); next }
    lines <- c(lines, sprintf("  bind_download_state(\"%s\", function() %s)", id, a$expr))
    add(basename(f), id, fn, "bound", a$expr, a$note)
  }
  # Downloads registered inside observers: listed, not bound
  nested <- setdiff(unique(unlist(find_downloads(body))), bound_ids)
  for (id in nested) add(basename(f), id, "", "NOT BOUND", "", "defined inside an observer")
  block <- c("", "  ############################## Download buttons state ##############################",
             "  # Each download button stays disabled until the result it exports exists (generated from the first req() of each download)",
             lines)
  writeLines(block, file.path(out_dir, paste0(sub("\\.R$", "", basename(f)), ".download_state.R")))
}
inv <- do.call(rbind, inventory)
write.csv(inv, file.path(out_dir, "download_inventory.csv"), row.names = FALSE)
print(table(inv$file, inv$status))
cat("\nNot bound:\n"); print(inv[inv$status != "bound", c("file", "id", "note")], row.names = FALSE)
