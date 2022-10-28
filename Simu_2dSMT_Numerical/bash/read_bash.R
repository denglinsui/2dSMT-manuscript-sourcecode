#' Simulation for one case
if (!requireNamespace("argparse", quietly = TRUE)){
  install.packages("argparse")}

library("argparse")

parser <- ArgumentParser()

parser$add_argument("--mu_type", type = "character", default = "Medium", help = "Type of mu: Sparse, Medium or Dense.")
parser$add_argument("--Cov_type", type = "character", default = "Medium", help = "Type of Covariance Dependency of noise: Weak, Medium or Strong.")
parser$add_argument("--magnitude", type = "double", default = "1", help = "The magnitude of signal.")
parser$add_argument("--est_cov", type = "logical", default = TRUE, help = "Whether estimate covariance matrix in simulation.")

args <- parser$parse_args()

mu_type <- args$mu_type
Cov_type <- args$Cov_type
magnitude <- args$magnitude
est_cov <- args$est_cov
estcov <- est_cov
#estcov <- F