% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEseq2_pairwise.R
\name{DEseq2_pairwise}
\alias{DEseq2_pairwise}
\title{Run post-hoc DESeq2 comparison at a given time point (Genus level)}
\usage{
DEseq2_pairwise(
  physeq,
  time_point,
  group = "treatment",
  time_var = "timeline",
  comparison = c("groupA", "groupB"),
  padj_cutoff = 0.05
)
}
\arguments{
\item{physeq}{A \code{phyloseq} object (raw, not yet aggregated).}

\item{time_point}{A string like "Week 0" or "Day 3" (must match sample_data()[\link{time_var}]).}

\item{group}{Name of grouping variable (default: "treatment").}

\item{time_var}{Name of time variable (default: "timeline").}

\item{comparison}{A character vector of length 2: c("groupA", "groupB").}

\item{padj_cutoff}{Adjusted p-value threshold (default: 0.05).}
}
\value{
A list with sig_table and plot.
}
\description{
Run post-hoc DESeq2 comparison at a given time point (Genus level)
}
