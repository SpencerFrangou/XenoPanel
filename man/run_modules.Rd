\name{run_modules}
\alias{run_crossreactivity}
\alias{run_qc_de}
\alias{run_decomposition}
\alias{run_entropy_pathways}
\title{The four analysis modules}
\usage{
run_crossreactivity(cfg, project_dir, out_dir)

run_qc_de(cfg, m, samples, out_dir, mask = NULL)

run_decomposition(cfg, de, out_dir)

run_entropy_pathways(cfg, de, samples, out_dir, project_dir)
}
\arguments{
  \item{cfg}{The parsed YAML configuration (a list).}
  \item{project_dir}{Directory against which relative paths in the configuration resolve.}
  \item{out_dir}{Output directory.}
  \item{m}{Count matrix (features by samples) as returned by \code{read_counts} and
  \code{align_counts}.}
  \item{samples}{Sample sheet as returned by \code{read_samples}.}
  \item{mask}{Character vector of feature IDs to remove before differential expression, or
  \code{NULL}.}
  \item{de}{The list returned by \code{run_qc_de} (normalised matrices and fits).}
}
\value{
\code{run_crossreactivity} returns the exclusion mask (character vector); \code{run_qc_de}
returns a list with the normalised matrices, fits and contrast tables; the other two return
\code{NULL} invisibly. All four write CSV tables into \code{out_dir}.
}
\description{
The modules called in order by \code{\link{xenopanel_run}}; they can also be called one at a
time from R.
}
