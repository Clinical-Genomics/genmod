import pkg_resources

ensembl_file_37 = 'annotations/ensembl_genes_37.txt.gz'
ensembl_file_38 = 'annotations/ensembl_genes_38.txt.gz'

ensembl_path_37 = pkg_resources.resource_filename('genmod', ensembl_file_37)
ensembl_path_38 = pkg_resources.resource_filename('genmod', ensembl_file_38)