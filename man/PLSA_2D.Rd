\name{PLSA_2D}
\alias{PLSA_2D}

\title{
2D PLSA
}
\description{
Computes a PLSA decomposition for a 2D tensor (matrix). The full tensor is never constructed, only the nonzeros are used for faster calculations. 
}
\usage{
PLSA_2D(flat_format, n_latent = 2, niter = 2, all_factors = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flat_format}{A data frame (or matrix) describing the tensor. The first 2 columns are indices for the tensor, the 3rd column is the count data.}
  \item{n_latent}{number of latent factors}
  \item{niter}{number of iterations}
  \item{all_factors}{Should the conditional latent probabilities also should be calculated.}
}
\details{
The code doesn't check for duplicates in flat_format and having duplicates rows has undefined behavors. The first 2 columns can be character values, the unique values will be used as the dimension names for the factors. Currently uses random initialization, so should use multiple re-tries. 
}
\value{
List of factors.
\item{Pz}{Vector of probabilities of latent groups z}
\item{Pa_z}{Matrix of probabilities of A (1st tensor dimension) conditional on latent group z}
\item{Pb_z}{Matrix of probabilities of B (2nd tensor dimension) conditional on latent group z}
\item{conv}{Vector of cost vs iteration}

If \code{all_factors = TRUE} also returns
\item{Pz_a}{matrix}
\item{Pz_b}{matrix}

}


\seealso{\code{\link{PLSA_4D}} , \code{\link{PLSA_3D}} }

\examples{
out = PLSA_2D(my_data_frame,2,2)
}

