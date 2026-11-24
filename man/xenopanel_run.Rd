\name{xenopanel_run}
\alias{xenopanel_run}
\alias{xenopanel_version}
\title{Run the XenoPanel pipeline from a YAML configuration file}
\usage{
xenopanel_run(cfg_path)

xenopanel_version()
}
\arguments{
  \item{cfg_path}{Path to a YAML configuration file; see \code{config.example.yaml} in the
  repository and \code{docs/INPUT_FORMAT.md} for the layout.}
}
\value{
\code{xenopanel_run} returns the output directory invisibly after writing all tables and
\code{sessionInfo.txt} into it. \code{xenopanel_version} returns the version string.
}
\description{
Reads the configuration, runs the enabled modules in order (cross-reactivity, QC and
differential expression, decomposition, entropy and pathway activity) and writes every
output under \code{<project_dir>/<output_dir>}. Each module stops with a clear message if a
package it needs is not installed; see the package \code{Suggests} field.
}
\examples{
\dontrun{
xenopanel_run("config.example.yaml")
}
}
