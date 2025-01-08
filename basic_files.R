
# Creating basic files

VER <- "V4" # first run ""
ORTHOLOG_GROUP_FILE <- "Data/DATA/Orthogroups_20240823_clean.tsv" 


species_list <- c("Lodge", "Asp", "Nor","Scots","Birch","Cher")
combo <- data.frame(t(combn(species_list, 2)))

file_list_1 <- c()
for (i in 1:nrow(combo)){
  
  s1 <- combo[i, "X1"]
  s2 <- combo[i, "X2"]
  
  file_name <- paste0("Data/DATA/comparisonFiles/comparison-",s1,"_",s2,"-pearsonMR0.03no-table-",VER,".RData")
  file_list_1 <- append(file_list_1, file_name, after = length(file_list_1))
  
}

## COMBINING ALL COMPARISON FILES FROM COMPLEX (all_xpressed_genes)

expr_genes <- c()
for (x in file_list_1){
  if(file.exists(x)){
    
    # x <- file_list_1[1]
    load(x)
    
    key_word_s1 <- sapply(strsplit(x, "_"), "[",1) 
    key_word_s1 <- sapply(strsplit(key_word_s1, "-"), "[",2) 
    key_word_s2 <- sapply(strsplit(x, "_"), "[",2)
    key_word_s2 <- sapply(strsplit(key_word_s2, "-"), "[",1) 
    
    
    genes <- comparison_table %>%
      as.data.frame() %>% 
      mutate_at("Max.p.val", as.numeric) 
    
    expr_genes <- rbind(expr_genes, data.frame(
      OrthoGroup = genes$OrthoGroup,
      Species1 = rep(key_word_s1, nrow(genes)),
      GeneSpecies1 = genes$Species1,
      Species2 = rep(key_word_s2, nrow(genes)),
      GeneSpecies2 = genes$Species2,
      Species1pVal = genes$Species1.p.val, 
      Species2pVal = genes$Species2.p.val,
      Max.p.Val =  genes$Max.p.val))
    
    print(paste0(key_word_s1, "-", key_word_s2)) 
  }
}

expr_genes <- expr_genes %>% 
  arrange(OrthoGroup)

save(expr_genes, file = "Data/DATA/all_expressed_genes.RData")



## ORTHOGROUPS WITH BINARY INPUT FOR PRESENCE/ABSENCE OF EXPRESSED GENE PAIRS (ones_and_zeros_all)

ortholog_group_file <- read.delim(ORTHOLOG_GROUP_FILE, header = TRUE, sep = "\t")

ortho_general_filtering <- ortholog_group_file %>%
  mutate(Pinus_sylvestris_cp = Pinus_sylvestris) %>% 
  rename(
    Aspen = Populus_tremula,
    Birch = Betula_pendula,
    `Norway spruce` = Picea_abies,
    `Scots pine` = Pinus_sylvestris,
    Cherry = Prunus_avium,
    `Lodgepole pine` = Pinus_sylvestris_cp) %>% 
  select(OrthoGroup, all_of(species_list_full_names))

ones_and_zeros_all <- ortho_general_filtering 

for (x in file_list_1){
  if (file.exists(x)){ 
    load(x) 
    
    
    key_word_s1 <- sapply(strsplit(x, "_"), "[",1) 
    key_word_s1 <- sapply(strsplit(key_word_s1, "-"), "[",2) 
    key_word_s2 <- sapply(strsplit(x, "_"), "[",2)
    key_word_s2 <- sapply(strsplit(key_word_s2, "-"), "[",1)
    
    
    df <- comparison_table%>%
      select(Species1, Species2, OrthoGroup, Max.p.val) %>%
      mutate_at("Max.p.val", as.numeric) %>%
      group_by(OrthoGroup) %>%
      slice(1) 
    
    new_col_name <-as.character(paste0(key_word_s1, "-", key_word_s2))
    
    colnames(df)[1] = new_col_name
    
    df <- df %>%
      select(-c(Species2, Max.p.val))
    
    new_column <- ones_and_zeros_all %>%
      left_join(df, join_by(OrthoGroup == OrthoGroup)) %>%
      select(all_of(new_col_name))
    
    ones_and_zeros_all <- cbind(ones_and_zeros_all, new_column)
    
    print(new_col_name)
    
  }}

ones_and_zeros_all <- ones_and_zeros_all %>% 
  column_to_rownames(var = "OrthoGroup") %>% 
  select(-c(1:6))

ones_and_zeros_all[is.na(ones_and_zeros_all)] <- 0
ones_and_zeros_all[ones_and_zeros_all != 0] <-1

ones_and_zeros_all <- ones_and_zeros_all %>%
  mutate_if(is.character, as.integer)

col_order <- c(  "Asp-Birch",  "Asp-Cher", "Birch-Cher",  "Lodge-Asp", "Lodge-Birch", "Lodge-Cher",  "Asp-Nor","Nor-Birch", "Nor-Cher", "Scots-Birch", "Scots-Cher",  "Asp-Scots" ,  "Lodge-Scots", "Lodge-Nor", "Nor-Scots")

ones_and_zeros_all<- ones_and_zeros_all[, col_order]

save(ones_and_zeros_all, file = "Data/DATA/ones_and_zeros_all.RData")



