# UCell_Seuratv4

This is a Seurat v4-only fork of UCell, adapted to use Seurat v4 slot API (GetAssayData(..., slot=...)).

Install from GitHub (after you push to your GitHub repo):
  devtools::install_github("shobhitagrawal1/UCell_Seurat_v4")

This fork pins Seurat to >=4, <5. Use this repository if you are running Seurat v4 and want the UCell scoring tools without Seurat v5 layer logic.
