from importlib_resources import files

ensembl_file_37 = 'annotations/ensembl_genes_37.txt.gz'
ensembl_file_38 = 'annotations/ensembl_genes_38.txt.gz'

ensembl_path_37 = files('genmod').join_path(ensembl_file_37)
ensembl_path_38 = files('genmod').join_path(ensembl_file_38)