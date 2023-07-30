a_old <- function(...) {
  A <- list(...)
  Anames <- names(A)
  Alength <- length(A)
  B <- A[[20]]
  invisible()
}
a_new <- function(...) {
  A <- list(...)
  Anames <- ...names()
  Alength <- ...length()
  B <- ..20
  invisible()
}

a_new2 <- function(...) {
  Anames <- ...names()
  Alength <- ...length()
  B <- ..20
  invisible()
}

L <- as.list(seq_len(100))
names(L) <- as.character(L)
profvis::profvis({
  bench::mark(
    old = do.call(a_old, L),
    new = do.call(a_new, L),
    new = do.call(a_new2, L)
  )
})
