sig_lines <- readLines('inst/stan/helper/functions/joinme_dynpred_partial.stanfunctions')
call_lines <- readLines('inst/stan/joinme_dynpred.stan')

sig_start <- grep('^real partial_draw\(', sig_lines)
sig_end <- sig_start + which(sig_lines[sig_start:length(sig_lines)] == ') {')[1] - 1
sig_block <- paste(sig_lines[sig_start:sig_end], collapse = '\n')

call_start <- grep('target \+= partial_draw\(', call_lines)
call_end <- call_start + which(trimws(call_lines[call_start:length(call_lines)]) == ');')[1] - 1
call_block <- paste(call_lines[call_start:call_end], collapse = '\n')

extract_sig_names <- function(txt) {
  lines <- strsplit(txt, '\n', fixed = FALSE)[[1]]
  lines <- trimws(lines)
  lines <- lines[!grepl('^real partial_draw\($|^\) \{$|^/\*|^\*/|^\*', lines)]
  lines <- gsub(',+$', '', lines)
  lines <- lines[nzchar(lines)]
  sub('.*?([A-Za-z_][A-Za-z0-9_]*)\s*(\[[^]]*\])?$', '\\1', lines)
}

extract_call_names <- function(txt) {
  lines <- strsplit(txt, '\n', fixed = FALSE)[[1]]
  lines <- trimws(lines)
  lines <- lines[!grepl('^target \+= partial_draw\($|^\);$|^/\*|^\*/|^\*', lines)]
  lines <- gsub(',+$', '', lines)
  lines <- lines[nzchar(lines)]
  lines
}

sig_names <- extract_sig_names(sig_block)
call_names <- extract_call_names(call_block)
cat('sig_n=', length(sig_names), '\n', sep='')
cat('call_n=', length(call_names), '\n', sep='')
for (i in seq_len(min(length(sig_names), length(call_names)))) {
  if (!identical(sig_names[[i]], call_names[[i]])) {
    cat('MISMATCH at ', i, ': sig=', sig_names[[i]], ' call=', call_names[[i]], '\n', sep='')
    quit(save='no', status=0)
  }
}
cat('NO_MISMATCH\n')
