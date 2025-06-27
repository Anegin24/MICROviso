#' make_tree: Build and optionally plot a phylogenetic tree from ASV sequences in a phyloseq object
#'
#' @param ps A phyloseq object containing ASV sequences in the refseq slot (DNAStringSet)
#' @param msa_method Alignment method to use: "ClustalW", "Muscle", or "ClustalOmega"
#' @param optimize_ml Logical; whether to optimize the tree using a maximum likelihood model (default: TRUE)
#' @param model Substitution model used for ML tree optimization (default: "GTR")
#' @param optimize_gamma Logical; whether to optimize the gamma rate parameter (default: TRUE)
#' @param return_tree_only Logical; if TRUE, return only the tree (phylo object). If FALSE, return a phyloseq object with the tree attached.
#' @param plot_tree Logical; whether to plot the resulting tree using ggtree (default: TRUE)
#'
#' @return A `phylo` object (tree only) or an updated `phyloseq` object with the phylogenetic tree attached
#'
#' @examples
#' ps <- make_tree(ps, msa_method = "ClustalOmega", optimize_ml = TRUE)
#' tree <- make_tree(ps, optimize_ml = FALSE, return_tree_only = TRUE)
#' @export
make_tree <- function(ps,
                      msa_method = "ClustalW",
                      optimize_ml = TRUE,
                      model = "GTR",
                      optimize_gamma = TRUE,
                      return_tree_only = FALSE,
                      plot_tree = TRUE) {
  require(msa)
  require(Biostrings)
  require(phangorn)
  require(phyloseq)
  require(ape)
  require(ggtree)

  if (is.null(refseq(ps))) {
    stop("No reference sequences found in the phyloseq object. Please assign a DNAStringSet to refseq(ps).")
  }

  cat("Aligning sequences using method:", msa_method, "...\n")
  seqs <- refseq(ps)
  names(seqs) <- taxa_names(ps)
  alignment <- msa(seqs, method = msa_method)
  aligned <- as(alignment, "DNAStringSet")
  aligned_matrix <- as.matrix(aligned)

  cat("Converting alignment to phangorn format...\n")
  phangDat <- phyDat(aligned_matrix, type = "DNA")

  cat("Constructing NJ tree from distance matrix...\n")
  dm <- dist.ml(phangDat)
  treeNJ <- NJ(dm)

  if (optimize_ml) {
    cat("Optimizing tree using Maximum Likelihood with model:", model, "...\n")
    fit <- pml(treeNJ, data = phangDat)
    fit_opt <- optim.pml(fit,
                         model = model,
                         optInv = TRUE,
                         optGamma = optimize_gamma,
                         rearrangement = "stochastic",
                         control = pml.control(trace = 0))
    tree_final <- midpoint(fit_opt$tree)
  } else {
    cat("Skipping ML optimization; using NJ tree directly.\n")
    tree_final <- midpoint(treeNJ)
  }

  if (plot_tree) {
    cat("Plotting phylogenetic tree...\n")
    print(
      ggtree(tree_final) +
        geom_tiplab(size = 2) +
        theme_tree2() +
        ggtitle("Phylogenetic Tree")
    )
  }

  if (return_tree_only) {
    return(tree_final)
  } else {
    cat("Merging tree into the phyloseq object...\n")
    ps_out <- merge_phyloseq(ps, phy_tree(tree_final))
    return(ps_out)
  }
}
