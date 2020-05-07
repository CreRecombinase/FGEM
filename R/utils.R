verbose_message_factory <- function(verbose) {
    if (verbose) {
        function(...) {
            message(...)
        }
    } else {
        function(...) {
            invisible(NULL)
        }
    }
}
