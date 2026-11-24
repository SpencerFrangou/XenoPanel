\name{readers}
\alias{read_counts}
\alias{read_samples}
\alias{align_counts}
\title{Read and align the count matrix and sample sheet}
\usage{
read_counts(path, collapse = FALSE)

read_samples(path)

align_counts(m, samples)
}
\arguments{
  \item{path}{Path to a CSV or TSV file; see \code{docs/INPUT_FORMAT.md}.}
  \item{collapse}{Collapse duplicate gene symbols by summing their counts.}
  \item{m}{Count matrix as returned by \code{read_counts}.}
  \item{samples}{Sample sheet as returned by \code{read_samples}.}
}
\value{
\code{read_counts} returns an integer matrix (features by samples); \code{read_samples} a
data frame with at least the columns \code{sample} and \code{group}; \code{align_counts} the
matrix with columns ordered and named by the sample sheet.
}
\description{
Readers that fail loudly on malformed input.
}
