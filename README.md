# UCell_Seuratv4

This is a Seurat v4-only fork of UCell, adapted to use Seurat v4 slot API (GetAssayData(..., slot=...)).

Install from GitHub (after you push to your GitHub repo):

devtools::install_github("shobhitagrawal1/UCell_Seurat_v4")

library(UCellSeuratv4)

usage:
seu_obj=AddModuleScore_UCell_v4(seu_obj, ... )

This fork pins Seurat to >=4, <5. Use this repository if you are running Seurat v4 and want the UCell scoring tools without Seurat v5 layer logic.

If you use this package, please cite the original UCell publication: 
UCell: robust and scalable single-cell gene signature scoring. Massimo Andreatta & Santiago J Carmona (2021) CSBJ https://doi.org/10.1016/j.csbj.2021.06.043

Revert to the original repo in case of doubt or errors:
https://github.com/carmonalab/UCell
