library(testthat)
library(Seurat)
# library(UCell_Seuratv4) <- will be available after install

test_that("AddModuleScore_UCell runs on simple Seurat v4 object", {
  mat <- matrix(rpois(200, lambda=2), nrow=20)
  rownames(mat) <- paste0("G", 1:20)
  colnames(mat) <- paste0("C", 1:10)
  obj <- Seurat::CreateSeuratObject(mat)
  gene.sets <- list(sig1 = c("G1","G2","G3"))
  # call the function directly after loading package in tests
  res <- AddModuleScore_UCell(obj, features = gene.sets)
  expect_true("sig1_UCell" %in% colnames(res[[]]))
  vals <- res[["sig1_UCell"]]
  expect_true(is.numeric(vals))
  expect_true(all(vals >= 0 & vals <= 1))
})
