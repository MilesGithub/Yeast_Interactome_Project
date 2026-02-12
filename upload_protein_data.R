
con_1 <- dbConnect(RMySQL::MySQL(),
                   dbname='YeRI',
                   host='localhost',
                   user='root',
                   password='')

yeast_gene_df<-read.delim('data/Saccharomyces_cerevisiae_genes.tsv')
uniprot_yeast_gene_df<-read.delim('data/uniprotkb_taxonomy_id_559292_AND_review_2026_02_10.tsv')


for (i in seq_len(nrow(yeast_gene_df))) {
  
  print(i)
  row <- yeast_gene_df[i, ]
  
  try({
    
    description <- uniprot_yeast_gene_df$description[
      uniprot_yeast_gene_df$ensembl_gene_id == row$ensembl_id
    ]
    if (length(description) == 0) description <- NA
    
    sequence <- uniprot_sub$Sequence[
      uniprot_sub$Entry == row$uniprot_id
    ]
    protein_name <- uniprot_sub$Protein.names[
      uniprot_sub$Entry == row$uniprot_id
    ]
    
    if (length(sequence) == 0) sequence <- NA
    if (length(protein_name) == 0) protein_name <- NA
    
    cols <- c()
    vals <- c()
    
    if (!is.na(row$gene_symbol) && row$gene_symbol != "") {
      cols <- c(cols, "`gene_name`")
      vals <- c(vals, shQuote(row$gene_symbol))
    }else{
      cols <- c(cols, "`gene_name`")
      vals <- c(vals, shQuote(row$ensembl_id))
    }
    
    if (!is.na(row$uniprot_id) && row$uniprot_id != "") {
      cols <- c(cols, "`uniprot_id`")
      vals <- c(vals, shQuote(row$uniprot_id))
    }
    
    if (!is.na(row$ensembl_id) && row$ensembl_id != "") {
      cols <- c(cols, "`ensembl_id`")
      vals <- c(vals, shQuote(row$ensembl_id))
    }
    
    if (!is.na(row$entrezgene_id)) {
      cols <- c(cols, "`entrez_id`")
      vals <- c(vals, shQuote(row$entrezgene_id))
    }
    
    if (!is.na(row$chromosome) && row$chromosome != "") {
      cols <- c(cols, "`chromosome`")
      vals <- c(vals, shQuote(row$chromosome))
    }
    
    if (!is.na(row$start)) {
      cols <- c(cols, "`start`")
      vals <- c(vals, row$start)
    }
    
    if (!is.na(row$end)) {
      cols <- c(cols, "`end`")
      vals <- c(vals, row$end)
    }
    
    if (!is.na(description) && description != "") {
      cols <- c(cols, "`description`")
      vals <- c(vals, shQuote(description[1]))
    }
    
    if (!is.na(sequence) && sequence != "") {
      cols <- c(cols, "`sequence`")
      vals <- c(vals, shQuote(sequence[1]))
    }
    
    if (!is.na(protein_name) && protein_name != "") {
      cols <- c(cols, "`protein_name`")
      vals <- c(vals, shQuote(protein_name[1]))
    }
    
    if (length(cols) == 0) next
    
    sql <- paste0(
      "INSERT INTO `protein` (",
      paste(cols, collapse = ", "),
      ") VALUES (",
      paste(vals, collapse = ", "),
      ")"
    )
    
    dbExecute(con, sql)
  })
}



