

con_1 <- dbConnect(RMySQL::MySQL(),
                   dbname='YeRI',
                   host='localhost',
                   user='root',
                   password='')


strain_ids_and_single_mutant_fitness<-read.delim('data/strain_ids_and_single_mutant_fitness.tsv')

SGA_ExE<-read.delim('data/SGA_ExE.txt')
SGA_ExN<-read.delim('data/SGA_ExN.txt')
SGA_NxN<-read.delim('data/SGA_NxN.txt')


sql<-paste(c("SELECT * FROM `protein`",";"),collapse='')
query_obj<-dbSendQuery(con_1, sql)
protein_df<-fetch(query_obj, n = -1)


for (i in seq_len(nrow(protein_df))) {
  
  print(i)
  
  try({
    
    protein_row <- protein_df[i, ]
    ensembl_id  <- protein_row$ensembl_id
    
    # Subset matching fitness rows
    fitness_rows <- strain_ids_and_single_mutant_fitness[
      strain_ids_and_single_mutant_fitness$Systematic.gene.name == ensembl_id,
    ]
    
    fitness_values <- as.numeric(fitness_rows$Single.mutant.fitness..30..)
    
    # Keep only valid numeric values
    fitness_values <- fitness_values[is.finite(fitness_values)]
    
    # If no valid values → default to 1
    if (length(fitness_values) == 0) {
      mean_fitness <- 1
    } else {
      mean_fitness <- mean(fitness_values)
    }
    
    # Parameterized query (safe)
    dbExecute(
      con_2,
      "UPDATE protein 
       SET single_mutant_fitness = ?
       WHERE ensembl_id = ?",
      params = list(mean_fitness, ensembl_id)
    )
    
  })
}

sql<-paste(c("SELECT * FROM `interaction`",";"),collapse='')
query_obj<-dbSendQuery(con_1, sql)
interaction_df<-fetch(query_obj, n = -1)

# Combine datasets
combined_sga <- rbind(SGA_ExE, SGA_ExN, SGA_NxN)

# Extract interactor IDs
combined_sga$query_interactor_id <- sub("_.*", "", combined_sga$query_strain_id)
combined_sga$array_interactor_id <- sub("_.*", "", combined_sga$array_strain_id)


for (i in seq_len(nrow(interaction_df))) {
  
  print(i)
  
  interaction_row   <- interaction_df[i, ]
  interaction_id    <- interaction_row$id
  interactor_A_id   <- interaction_row$interactor_A_id
  interactor_B_id   <- interaction_row$interactor_B_id
  
  # Match either orientation
  matching_rows <- combined_sga[
    (combined_sga$query_interactor_id == interactor_A_id &
       combined_sga$array_interactor_id == interactor_B_id) |
      (combined_sga$query_interactor_id == interactor_B_id &
         combined_sga$array_interactor_id == interactor_A_id),
  ]
  
  if (nrow(matching_rows) > 0) {
    
    # Select lowest p-value
    best_index <- which.min(matching_rows$p_value)
    best_match <- matching_rows[best_index, ]
    
    # Select required columns
    insert_values <- best_match[, c(
      "query_strain_id",
      "query_allele_name",
      "query_interactor_id",
      "array_strain_id",
      "array_allele_name",
      "array_interactor_id",
      "arraytype_temp",
      "genetic_interaction_score",
      "p_value",
      "query_single_mutant_fitness_smf",
      "array_smf",
      "double_mutant_fitness",
      "double_mutant_fitness_standard_deviation"
    )]
    
    # Clean NA / NaN properly
    insert_values[] <- lapply(insert_values, function(col) {
      if (is.numeric(col)) {
        col[!is.finite(col)] <- NA
      } else {
        col[col == "NaN"] <- NA
      }
      col
    })
    
    try({
      dbExecute(
        con_1,
        paste0(
          "INSERT INTO genetic_interaction (
            query_strain_id,
            query_allele_name,
            query_interactor_id,
            array_strain_id,
            array_allele_name,
            array_interactor_id,
            arraytype_temp,
            genetic_interaction_score,
            p_value,
            query_single_mutant_fitness_smf,
            array_smf,
            double_mutant_fitness,
            double_mutant_fitness_standard_deviation
          ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
        ),
        params = as.list(insert_values[1, ])
      )
    })
  }
}
