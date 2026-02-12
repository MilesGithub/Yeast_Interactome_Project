

con_1 <- dbConnect(RMySQL::MySQL(),
                   dbname='YeRI',
                   host='localhost',
                   user='root',
                   password='')

ValBin<-read.delim("C:/Users/Miles/Downloads/ValBin-25.tsv")

for(i in 1:nrow(ValBin)){
  
  print(i)
  row<-ValBin[i,]
  
  try({
    orf_a <- row$orf_name_a
    orf_b <- row$orf_name_b
    
    id_a <- protein_test[protein_test$ensembl_id == orf_a,]$id
    id_b <- protein_test[protein_test$ensembl_id == orf_b,]$id
    
    sql <- paste0(
      "INSERT INTO `interaction` (`interactor_A`, `interactor_B`, `interactor_A_id`, `interactor_B_id`) VALUES (",
      shQuote(id_a), ", ",
      shQuote(id_b), ", ",
      shQuote(orf_a), ", ",
      shQuote(orf_b), ")"
    )
    
    dbExecute(con_1, sql)
  })
}


sql<-paste(c("SELECT * FROM `interaction`",";"),collapse='')
query_obj<-dbSendQuery(con_1, sql)
interaction_df<-fetch(query_obj, n = -1)


# --------------------------------------------------
#
# Add join tables to datasets and categories
#
# --------------------------------------------------

for (i in seq_len(nrow(ValBin))) {
  
  print(i)
  
  val_row <- ValBin[i, ]
  
  try({
    
    # ORF identifiers
    orf_a <- val_row$orf_name_a
    orf_b <- val_row$orf_name_b
    
    # Protein IDs
    id_a <- protein_test[protein_test$ensembl_id == orf_a, ]$id
    id_b <- protein_test[protein_test$ensembl_id == orf_b, ]$id
    
    # Interaction lookup
    interaction_row <- interaction_df[
      interaction_df$interactor_A_id == orf_a &
        interaction_df$interactor_B_id == orf_b,
    ]
    
    if (nrow(interaction_row) > 0 && !is.na(interaction_row$id)) {
      
      interaction_id <- interaction_row$id
      
      # Flags from ValBin
      in_i3d_exp_24   <- val_row$in_I3D.exp.24
      in_lit_bm_24    <- val_row$in_Lit.BM.24
      in_y2h_union_25 <- val_row$in_Y2H.union.25
      in_afrf_core    <- val_row$in_AFRF.core
      
      # --------------------------------------------------
      # Literature BM 24
      # --------------------------------------------------
      if (in_lit_bm_24 == "True") {
        
        sql_insert_category_lit <- paste0(
          "INSERT INTO interaction_interaction_category ",
          "(interaction_category_id, interaction_id) VALUES (1, '",
          interaction_id, "')"
        )
        dbExecute(con_1, sql_insert_category_lit)
        
        sql_insert_dataset_lit <- paste0(
          "INSERT INTO interaction_dataset ",
          "(dataset_id, interaction_id) VALUES (1, '",
          interaction_id, "')"
        )
        dbExecute(con_1, sql_insert_dataset_lit)
      }
      
      # --------------------------------------------------
      # Y2H Union 25
      # --------------------------------------------------
      if (in_y2h_union_25 == "True") {
        
        sql_insert_category_y2h <- paste0(
          "INSERT INTO interaction_interaction_category ",
          "(interaction_category_id, interaction_id) VALUES (2, '",
          interaction_id, "')"
        )
        dbExecute(con_1, sql_insert_category_y2h)
        
        y2h_row <- Y2H_union_25[
          Y2H_union_25$orf_name_a == orf_a &
            Y2H_union_25$orf_name_b == orf_b,
        ]
        
        uetz_screen <- y2h_row$Uetz.screen
        ito_core    <- y2h_row$Ito.core
        ccsb_yi1    <- y2h_row$CCSB.YI1
        yeri        <- y2h_row$YeRI
        
        if (yeri == "TRUE") {
          dbExecute(con_1, paste0(
            "INSERT INTO interaction_dataset ",
            "(dataset_id, interaction_id) VALUES (2, '",
            interaction_id, "')"
          ))
        }
        
        if (ito_core == "TRUE") {
          dbExecute(con_1, paste0(
            "INSERT INTO interaction_dataset ",
            "(dataset_id, interaction_id) VALUES (3, '",
            interaction_id, "')"
          ))
        }
        
        if (uetz_screen == "TRUE") {
          dbExecute(con_1, paste0(
            "INSERT INTO interaction_dataset ",
            "(dataset_id, interaction_id) VALUES (4, '",
            interaction_id, "')"
          ))
        }
        
        if (ccsb_yi1 == "TRUE") {
          dbExecute(con_1, paste0(
            "INSERT INTO interaction_dataset ",
            "(dataset_id, interaction_id) VALUES (5, '",
            interaction_id, "')"
          ))
        }
      }
      
      # --------------------------------------------------
      # I3D Experimental 24
      # --------------------------------------------------
      if (in_i3d_exp_24 == "True") {
        
        dbExecute(con_1, paste0(
          "INSERT INTO interaction_interaction_category ",
          "(interaction_category_id, interaction_id) VALUES (3, '",
          interaction_id, "')"
        ))
        
        dbExecute(con_1, paste0(
          "INSERT INTO interaction_dataset ",
          "(dataset_id, interaction_id) VALUES (7, '",
          interaction_id, "')"
        ))
      }
      
      # --------------------------------------------------
      # AFRF Core
      # --------------------------------------------------
      if (in_afrf_core == "True") {
        
        dbExecute(con_1, paste0(
          "INSERT INTO interaction_interaction_category ",
          "(interaction_category_id, interaction_id) VALUES (4, '",
          interaction_id, "')"
        ))
        
        dbExecute(con_1, paste0(
          "INSERT INTO interaction_dataset ",
          "(dataset_id, interaction_id) VALUES (8, '",
          interaction_id, "')"
        ))
      }
      
    }
    
  })
}

# ------------------------------------------------------------------
#
# Add number of interactions in database and number of drug screens
#
# ------------------------------------------------------------------

sql<-paste(c("SELECT * FROM `interaction_genetic_interaction_df`",";"),collapse='')
query_obj<-dbSendQuery(con_1, sql)
interaction_genetic_interaction_df<-fetch(query_obj, n = -1)

for (i in seq_len(nrow(interaction_df))) {
  
  print(i)
  
  interaction_row <- interaction_df[i, ]
  interaction_id  <- interaction_row$id
  
  # Link table: interaction → genetic interaction
  interaction_genetic_links <- interaction_genetic_interaction_df[
    interaction_genetic_interaction_df$interaction_id == interaction_id,
  ]
  
  # Retrieve associated genetic interactions
  genetic_interactions <- genetic_interaction_df[
    genetic_interaction_df$id %in% interaction_genetic_links$genetic_interaction_id,
  ]
  
  # Order by smallest p-value
  genetic_interactions <- genetic_interactions[
    order(genetic_interactions$p_value),
  ]
  
  if (nrow(genetic_interactions) > 0) {
    
    best_genetic_interaction <- genetic_interactions[1, ]
    
    best_p_value  <- best_genetic_interaction$p_value
    best_gi_score <- best_genetic_interaction$genetic_interaction_score
    
    update_query <- paste0(
      "UPDATE interaction SET ",
      "genetic_interaction_score = ", best_gi_score, ", ",
      "genetic_interaction_p_value = ", best_p_value,
      " WHERE id = ", interaction_id
    )
    
    dbExecute(con_1, update_query)
  }
}

