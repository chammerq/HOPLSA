\name{PLSA_3D}
\alias{PLSA_3D}

\title{
3D PLSA
}
\description{
Computes a PLSA decomposition for a 3D tensor. The full tensor is never constructed, only the nonzeros are used for faster calculations. 
}
\usage{
PLSA_3D(flat_format, n_latent = 2, niter = 2, all_factors = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flat_format}{A data frame (or matrix) describing the tensor. The first 3 columns are indices for the tensor, the 4th column is the count data.}
  \item{n_latent}{number of latent factors}
  \item{niter}{number of iterations}
  \item{all_factors}{Should the conditional latent probabilities also should be calculated.}
}
\details{
The code doesn't check for duplicates in flat_format and having duplicates rows has undefined behavors. The first 3 columns can be character values, the unique values will be used as the dimension names for the factors. Currently uses random initialization, so should use multiple re-tries. 
}
\value{
List of factors.
\item{Pz}{Vector of probabilities of latent groups z}
\item{Pa_z}{Matrix of probabilities of A (1st tensor dimension) conditional on latent group z}
\item{Pb_z}{Matrix of probabilities of B (2nd tensor dimension) conditional on latent group z}
\item{Pc_z}{Matrix of probabilities of C (3rd tensor dimension) conditional on latent group z}
\item{conv}{Vector of cost vs iteration}

If \code{all_factors = TRUE} also returns
\item{Pz_a}{matrix}
\item{Pz_b}{matrix}
\item{Pz_c}{matrix}

}


\seealso{\code{\link{PLSA_2D}} , \code{\link{PLSA_4D}} }

\examples{
out = PLSA_3D(my_data_frame,2,2)
}

