#' Parse transformation composition list to Stan format
#' @param tf_composition List with keys: tf_cv_mean_comp, tf_cv_marker_comp, tf_cs_mean_comp, tf_cs_marker_comp
#' @return List with codes and n_codes for each composition
#' @keywords internal
#' @noRd
.parse_tf_composition_list <- function(tf_composition) {
  # Expand named composition list into per-term code arrays
  tf_names <- c("tf_cv_mean_comp", "tf_cv_marker_comp", "tf_cs_mean_comp", "tf_cs_marker_comp")
  result <- list()
  
  for (name in tf_names) {
    tf_str <- tf_composition[[name]] %||% "none"
    parsed <- .parse_tf_composition(tf_str)
    result[[name]] <- list(
      codes = as.integer(parsed$codes),
      n_codes = as.integer(parsed$n_codes),
      description = parsed$description
    )
  }
  
  result
}

#' Parse transformation composition string
#' @description
#' Convert strings like "2_3" or "log_sqrt" into array of transform IDs.
#' Transform codes: 0=none, 1=id, 2=log, 3=sqrt, 4=logit, 5=inv, 6=exp
#' 
#' @param tf_string Character transformation specification. Can be:
#'   - "none" or "" or 0: identity
#'   - Numeric codes: "2", "3", "2_3", "1_2_3"
#'   - Named: "log", "sqrt", "log_sqrt", "sqrt_log"
#' @param max_length Maximum array length (default 7)
#' @return List with: codes (integer array), names (character), description (character)
#' 
#' @keywords internal
#' @noRd
.parse_tf_composition <- function(tf_string, max_length = 7) {
  # Parse transform composition into numeric codes for Stan
  tf_string <- as.character(tf_string)
  tf_string <- tolower(trimws(tf_string))
  
  # Map named transforms
  tf_map <- list(
    "none" = c(),
    "id" = c(1L),
    "log" = c(2L),
    "sqrt" = c(3L),
    "logit" = c(4L),
    "inv" = c(5L),
    "exp" = c(6L),
    "log_sqrt" = c(2L, 3L),  # Apply sqrt first, then log: sqrt(x) then log(sqrt(x))
    "sqrt_log" = c(3L, 2L),  # Apply log first, then sqrt: log(x) then sqrt(log(x))
    "log_log" = c(2L, 2L),
    "sqrt_sqrt" = c(3L, 3L),
    "logit_log" = c(4L, 2L),
    "exp_log" = c(6L, 2L)
  )
  
  if (tf_string %in% names(tf_map)) {
    tf_codes <- tf_map[[tf_string]]
  } else if (tf_string == "0" || tf_string == "") {
    tf_codes <- c()
  } else if (grepl("^[0-6](_[0-6])*$", tf_string)) {
    # Parse numeric codes like "2_3"
    parts <- as.integer(strsplit(tf_string, "_")[[1]])
    if (any(parts < 0 | parts > 6)) {
      stop("Transform codes must be in 0-6, got: ", tf_string)
    }
    tf_codes <- parts
  } else {
    stop("Unknown transformation: ", tf_string, 
         ". Use codes (0-6) or names (none, log, sqrt, log_sqrt, etc)")
  }
  
  if (length(tf_codes) == 0) tf_codes <- 1L  # identity
  if (length(tf_codes) > max_length) {
    stop("Transformation composition too long (max ", max_length, " transforms)")
  }
  
  # Pad to max_length with 0 (will be skipped in Stan)
  tf_padded <- c(tf_codes, rep(0L, max_length - length(tf_codes)))
  
  # Create names
  tf_names <- c("none", "id", "log", "sqrt", "logit", "inv", "exp")
  tf_name_str <- paste(tf_names[tf_codes + 1], collapse = " then ")
  if (length(tf_codes) == 1 && tf_codes[1] == 1) tf_name_str <- "identity"
  
  list(
    codes = tf_padded,
    n_codes = length(tf_codes),
    description = tf_name_str
  )
}