# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sparse_plsa_3d <- function(aye, jay, kay, qty, num_a, num_b, num_c, nz, niter) {
    .Call('_HOPLSA_sparse_plsa_3d', PACKAGE = 'HOPLSA', aye, jay, kay, qty, num_a, num_b, num_c, nz, niter)
}

sparse_plsa_2d <- function(aye, jay, qty, num_a, num_b, nz, niter) {
    .Call('_HOPLSA_sparse_plsa_2d', PACKAGE = 'HOPLSA', aye, jay, qty, num_a, num_b, nz, niter)
}

sparse_plsa_4d <- function(aye, jay, kay, ell, qty, num_a, num_b, num_c, num_d, nz, niter) {
    .Call('_HOPLSA_sparse_plsa_4d', PACKAGE = 'HOPLSA', aye, jay, kay, ell, qty, num_a, num_b, num_c, num_d, nz, niter)
}

