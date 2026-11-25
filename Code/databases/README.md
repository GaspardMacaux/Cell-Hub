# CellHub Database files

The Docker image `gaspardmacaux/cell-hub:2.0.0` already contains CellHub Database 2.1.2: nothing to do when you use the image.

To run the code outside Docker, download CellHub Database 2.1.2 from Zenodo ([doi:10.5281/zenodo.24197570](https://doi.org/10.5281/zenodo.24197570), all versions: take the latest) and copy the content of its `app_databases/` folder here, so that this folder holds one sub-folder per species:

```
Code/databases/human/CellHubDB_human_merged.rds
Code/databases/human/CellHubDB_human_CellChatDB_v2.rds
...
Code/databases/mouse/CellHubDB_mouse_merged.rds
...
```

The files are not versioned in this repository (size); their MD5 checksums are listed in the `MD5SUMS` file of the release.
