library(ggrepel)
library(msdap)
library(writexl)
library(openxlsx)
library(readxl)
library(dplyr)
library(VennDiagram)

######### 13. Loading Inputs from MSDAP2####
prot.input<-read.delim(paste0(FolderName,DateTimeStamp2,"/protein_abundance__input data as-is.tsv"))
dea<-readxl::read_excel(paste0(FolderName,DateTimeStamp2,"/differential_abundance_analysis.xlsx"), sheet = 2)
colnames(prot.input) <- c("Protein.ID", "FASTA", "Gene.ID",
                           "5_Ctrl","3_DR","4_DR","1_IP","2_IP")


prot.input <- prot.input[, c("Protein.ID", "FASTA", "Gene.ID", "1_IP", "2_IP", "3_DR", "4_DR", "5_Ctrl")]

dea <- dea[,c(1:3,6,10,14,16)] #this doesn't change
colnames(dea) <- c("protein_id","accessions",                                         
                    "fasta_headers","gene_symbols_or_id",                                
                    foldchange.colNam,
                    pvalue.colNam,
                    qvalue.colNam)

dea <- na.omit(dea)

assign(paste0("dea.", NameCond2, "vs", NameCond1), dea)

prot.input[[paste0("counter.vv.",NameCond2)]] <- rowSums(!is.na(prot.input[,c(col.start.cond2:col.end.cond2)]))
prot.input[[paste0("counter.vv.",NameCond1)]] <- rowSums(!is.na(prot.input[,c(col.start.cond1:col.end.cond1)]))

### PCA ###

### PCA 1
# prepare input matrix for PCA
PCA_mat <- prot.input[, 4:8]
rownames(PCA_mat) <- PCA_mat$Protein.ID
PCA_mat$Protein.ID <- NULL; PCA_mat$Gene.ID <- NULL; PCA_mat$FASTA <- NULL
PCA_mat <- na.omit(PCA_mat)

#perform pricipal component analyis, matrix transposed (input of prcomp: sample = rows)
pca <- prcomp(na.omit(as.matrix(t(PCA_mat))), scale = F) 


#calculate Variation for each experiment
#square of std.dev is variation
# calculate percentage of variation
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100,1) 

# generate a dataframe from the info to use for plotting
pca.var.data <- data.frame(Component = seq(1,length(pca.var.per),1),
                           Y = seq(1,length(pca.var.per),1),
                           X = pca.var.per)



#format for ggplot
pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1],
                       Y = pca$x[,2],
                       Z = pca$x[,3])


# assign pca groups test #assuming these are the groups you setup in the contrasts
pca.data$Group <- c(rep(NameCond1,cond1), rep(NameCond2,cond2), rep(NameControl, control))

#barplot(pca.var.per, main = "Scree plot", xlab = "Principal component", ylab = "Percent variation")
A <- ggplot(data = pca.var.data, aes(x = as.factor(Y), y = X)) + 
  geom_bar(stat = "identity")+
  ylab("Percentage of variation")+
  xlab("Principal component")+
  #ggtitle("Scree plot: Variance per PC")+
  theme_classic() + theme(legend.position = "none", axis.text=element_text(colour="black"))

cols <- c("#DCCB4E", "#E98905", "#3A9AB2")
B <- ggplot(data = pca.data, aes(x = X, y = Y, label = Sample)) +
  geom_point(aes(col = Group), size = 5) +
  geom_text(aes(label = Sample), hjust = -0.2, vjust = 0.5, size = 3) +
  #geom_point(shape = 1, size = 4, colour = "black") +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep =""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep =""))+
  #geom_label_repel(size = 2.8)+
  scale_color_manual(values = cols)+
  #ggtitle("PCA: log2 intensity")+
  theme_classic() + theme(legend.position = "none", axis.text=element_text(colour="black"))

(B|A)
PCAplot<-(B|A)
rm(A,B, PCA_mat, pca.data, pca, pca.var, pca.var.data)


######### 4. Protein Count and Intensities###################################
#prepare input and counter-vector
prot <- prot.input; prot_count <- c()
i = NULL
for(i in col.start.cond1:prot.col.end) { #change col numbers for sample number of runs
  temp <- subset.data.frame(prot, prot[,i] != "NA")
  temp <- nrow(temp)
  prot_count <- c(prot_count, temp)
  rm(temp)
}



### boxplot of intensties
par(mar = c(5, 4, 1, 2), cex.main = 0.9, mfrow = c(1,1), cex.axis = 0.9)
boxplot(prot[,col.start.cond1:prot.col.end],
        las = 2,
        col = c(rep("#DCCB4E" ,cond1), rep("#E98905",cond2), rep("#3A9AB2",control)), #change these to match number of replicates
        pch = 20,
        ylab = "log2 protein intensity")
Intensitiesplot<-recordPlot()

# barplot of IDs in GGPlot stylw
tmp <- data.frame(
  name=colnames(prot[,col.start.cond1:prot.col.end]) ,  #change columns to match sample numbers
  value=prot_count)

ggplot(tmp, aes(x=name, y=value)) + 
  geom_bar(stat = "identity", colour = "black", fill = c(rep("#DCCB4E" ,cond1), rep("#E98905",cond2), rep("#3A9AB2",control))) + 
  scale_x_discrete(limits=colnames(prot.input[,col.start.cond1:prot.col.end]))+                                                            #change numbers
  theme_classic()+ theme(legend.position = "none", axis.text=element_text(colour="black")) +
  ggtitle("")+
  ylab("protein groups") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Groupsplot<-recordPlot()

#----------------------------------------------------------------------------------------------------------
# List 1: Common proteins.
list1 <- dplyr::filter(prot.input, get(paste0("counter.vv.",NameCond1)) >= a & get(paste0("counter.vv.",NameCond2)) >= b)
write_xlsx(list1, path = file.path(FolderList, paste0("List1_", NameCond1, ".xlsx")))
# List 2: Enriched IP proteins.
list2.Fc1.5 <- list1[list1$Protein.ID %in% subset(dea, dea[[qvalue.colNam]] < pval & dea[[foldchange.colNam]] >= log2(foldlog2))$protein_id, ]
write_xlsx(list2.Fc1.5, path = file.path(FolderList, paste0("List2_", NameCond1, "FC1.5", ".xlsx")))
# List 2.1: Enriched DR proteins.
list2.1Fc1.5 <- list1[list1$Protein.ID %in% subset(dea, dea[[qvalue.colNam]] < pval & dea[[foldchange.colNam]] <= -log2(foldlog2))$protein_id, ]
write_xlsx(list2.1Fc1.5, path = file.path(FolderList, paste0("List2.1_", NameCond2, "FC1.5", ".xlsx")))
# List 3: Unique IP proteins
list3<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] == 0)
write_xlsx(list3, path = file.path(FolderList, paste0("List3_", NameCond1, ".xlsx")))
# List 3.1: Unique DR proteins
list3.1<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond2)]] >=a & prot.input[[paste0("counter.vv.",NameCond1)]] == 0)
write_xlsx(list3.1, path = file.path(FolderList, paste0("List3.1_", NameCond2, ".xlsx")))

#-----------------------------------------------------------------------------------------------------------


commons<-dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] >= b)
assign(paste0("commons.",NameCond1,"_", NameCond2),commons)
write_xlsx(commons, path = paste0(FolderList,"Commons_",NameCond2, "_", NameCond1, ".xlsx"))
rm(commons)

uniques<-dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] == 0)
assign(paste0("uniques.",NameCond1),uniques)
write_xlsx(uniques, path = paste0(FolderList,"List3", NameCond1, ".xlsx"))
rm(uniques)
uniques<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond2)]] >=b & prot.input[[paste0("counter.vv.",NameCond1)]] == 0)
assign(paste0("uniques.",NameCond2),uniques)
write_xlsx(uniques, path = paste0(FolderList,"List3_", NameCond1, ".xlsx"))
rm(uniques)

#-------------------------------------------------------------------------------------------------------------

# # Crear el diagrama de Venn
# venn.plot <- venn.diagram(
#   x = list(DataFrame1 = set1, DataFrame2 = set2),
#   filename = NULL, # para mostrar directamente en RStudio
#   fill = c("lightblue", "pink"),
#   alpha = 0.5,
#   cat.cex = 1.5,
#   cex = 2,
#   margin = 0.1
# )
# 
# # Dibujar diagrama en pantalla
# grid.newpage()
# grid.draw(venn.plot)

#-------------------------------------------------------------------------------------------------------------
# Count the rows for each category and condition
counts <- data.frame(
  Category = c("Commons", "Uniques", "Uniques"),
  Condition = c(paste0(NameCond1, "_", NameCond2), NameCond1, NameCond2),
  Count = c(
    nrow(get(paste0("commons.", NameCond1,"_",NameCond2))),
    nrow(get(paste0("uniques.", NameCond1))),
    nrow(get(paste0("uniques.", NameCond2)))
  ),
  Word = c("List1", "List3", "List3.1")
)


ggplot(counts, aes(x = Condition, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_text(aes(label = paste(Count, Word, sep = "\n")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,  # Move the labels above the bars
            size = 3.5) +  # Adjust label size as needed
  theme_classic() +
  theme(plot.title = element_text(margin = margin(b = 20))) +  # Add margin to title
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Extend y-axis limits
  labs(
    title = "Comparison of Exclusives, Commons, and Uniques",
    x = "Condition",
    y = "Count"
  ) +
  scale_fill_manual(values = c("Exclusives" = "#DCCB4E", "Commons" = "#E98905", "Uniques" = "#3A9AB2"))

Countsplot <- recordPlot()
rm(counts)


# Filter the most enriched proteins from both sides
enriched_proteins <- get(paste0("dea.", NameCond2, "vs", NameCond1)) %>%
  dplyr::filter(get(pvalue.colNam) <= pval & 
                  (get(foldchange.colNam) >= log2(foldlog2) | 
                     get(foldchange.colNam) <= -log2(foldlog2))) %>%
  dplyr::mutate(abs_foldchange = abs(get(foldchange.colNam))) %>%
  dplyr::arrange(desc(abs_foldchange)) %>%
  dplyr::group_by(sign = ifelse(get(foldchange.colNam) >= 0, "IP", "DR")) %>%
  dplyr::slice_head(n = 10) %>%  # Seleccionar los 10 principales
  dplyr::ungroup()


######### 14b. Generating Candidate Lists based on p-value and log############
# Initialize candidates column
dea_data <- get(paste0("dea.", NameCond2, "vs", NameCond1))
dea_data$candidates.p <- NA  # or c("")
x<-dea_data
for (n in 1:nrow(x)) {
  # Check for NA values before comparison
  if (!is.na(x[[foldchange.colNam]][n]) &&
      !is.na(x[[pvalue.colNam]][n]) &&
      abs(x[[foldchange.colNam]][n]) >= log2(foldlog2) && 
      x[[pvalue.colNam]][n] <= pval) {
    
    x$candidates.p[n] <- x$gene_symbols_or_id[n]
  } else {
    x$candidates.p[n] <- ""  # Set to empty string for non-matching rows
  }
}

dea_data<-x
assign(paste0("dea.", NameCond2, "vs", NameCond1, "."), dea_data)
rm(dea_data)
######### 14c. Generating Candidate Lists based on q-value and log############
# Initialize candidates column
dea_data <- get(paste0("dea.", NameCond2, "vs", NameCond1, "."))
dea_data$candidates.q <- NA  # or c("")
x<-dea_data
for (n in 1:nrow(x)) {
  # Check for NA values before comparison
  if (!is.na(x[[foldchange.colNam]][n]) &&
      !is.na(x[[qvalue.colNam]][n]) &&
      abs(x[[foldchange.colNam]][n]) >= log2(foldlog2) && 
      x[[qvalue.colNam]][n] <= qval) {
    
    x$candidates.q[n] <- x$gene_symbols_or_id[n]
  } else {
    x$candidates.q[n] <- ""  # Set to empty string for non-matching rows
  }
}

dea_data<-x
assign(paste0("dea.", NameCond2, "vs", NameCond1, "."), dea_data)
rm(dea_data)
################## 14.a PRINCIPAL VOLCANO ##########################
data = get(paste0("dea.", NameCond2, "vs", NameCond1, "."))
# Count the significant proteins on each side
count_positive <- nrow(data %>% dplyr::filter(get(pvalue.colNam) < pval & get(foldchange.colNam) >= log2(foldlog2)))
count_negative <- nrow(data %>% dplyr::filter(get(pvalue.colNam) < pval & get(foldchange.colNam) <= -log2(foldlog2)))

# Create the plot title with the counts in English and in the requested order
title_with_counts <- paste0(NameCond1, "vs", NameCond2,
                            " | Enriched proteins: IP (", count_positive, 
                            ") / DR (", count_negative, ")")

# Clean the data by removing rows with NA in the 'gene_symbols_or_id' column
data_clean <- data %>% filter(!is.na(gene_symbols_or_id))

# Create the plot without the 'NA' category in the legend
C <- ggplot(data = data_clean, 
            aes(x = get(foldchange.colNam), 
                y = -log10(get(pvalue.colNam)), 
                colour = ifelse(get(pvalue.colNam) <= pval & get(foldchange.colNam) >= log2(foldlog2), 
                                "IP", 
                                ifelse(get(pvalue.colNam) <= pval & get(foldchange.colNam) <= -log2(foldlog2), 
                                       "DR", 
                                       "Non-enriched proteins")))) + 
  geom_point(size = 0.9) +
  scale_colour_manual(name = NULL, # Remove legend title
                      values = c("IP" = '#003366', "DR" = '#FFA500', "Non-enriched proteins" = 'grey')) +
  theme_classic() + 
  theme(legend.position = "bottom", # Position the legend at the bottom
        legend.direction = "horizontal", # Display the legend horizontally
        axis.text = element_text(colour = "black")) + 
  ggtitle(title_with_counts) + 
  xlab("log2 fold change") + 
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(pval), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_vline(xintercept = log2(foldlog2), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_vline(xintercept = -log2(foldlog2), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_label_repel(data = enriched_proteins, 
                   aes(label = gene_symbols_or_id), 
                   size = 3, 
                   color = "black", 
                   fill = "lightgrey", 
                   box.padding = 0.5, 
                   point.padding = 0.5, 
                   max.overlaps = Inf) +
  guides(colour = guide_legend(nrow = 1)) # Make the legend horizontal below the title


DEAvolplot<-(C)
print(C)
rm(C)

#-------------------------------------------------------------------------------------------

# Crear los subconjuntos de proteínas enriquecidas
proteins_IP <- dea %>%
  filter(get(pvalue.colNam) < pval & get(foldchange.colNam) >= log2(foldlog2))

proteins_DR <- dea %>%
  filter(get(pvalue.colNam) < pval & get(foldchange.colNam) <= -log2(foldlog2))

# Verificar la estructura de los nuevos dataframes
print("Proteínas enriquecidas en IP:")
print(dim(proteins_IP))
print(head(proteins_IP))

print("Proteínas enriquecidas en DR:")
print(dim(proteins_DR))
print(head(proteins_DR))

# Crear el top 10 de proteínas más intensas en función de log2 fold change y -log10(p-value)
top10_IP <- proteins_IP %>%
  arrange(desc(get(foldchange.colNam)), desc(-log10(get(pvalue.colNam)))) %>%
  head(10)

top10_DR <- proteins_DR %>%
  arrange(desc(get(foldchange.colNam)), desc(-log10(get(pvalue.colNam)))) %>%
  head(10)

# Verificar los nuevos dataframes
print("Top 10 proteínas enriquecidas en IP:")
print(top10_IP)

print("Top 10 proteínas enriquecidas en DR:")
print(top10_DR)

top10_IP <- rename(top10_IP, "Protein.ID" = "protein_id","FASTA" = "fasta_headers",
                             "Gene.ID" = "gene_symbols_or_id")
top10_IP <- merge(top10_IP, list2.Fc1.5,by="Protein.ID")
top10_IP <- select(top10_IP, -"accessions", -paste0("counter.vv.",NameCond2), -paste0("counter.vv.",NameCond1),
                             -paste0(foldchange.colNam), -paste0(pvalue.colNam), -paste0(qvalue.colNam))

top10_DR <- rename(top10_DR, "Protein.ID" = "protein_id","FASTA" = "fasta_headers",
                   "Gene.ID" = "gene_symbols_or_id")
top10_DR <- merge(top10_DR, list2.1Fc1.5,by="Protein.ID")
top10_DR <- select(top10_DR, -"accessions", -paste0("counter.vv.",NameCond2), -paste0("counter.vv.",NameCond1),
                              -paste0(foldchange.colNam), -paste0(pvalue.colNam), -paste0(qvalue.colNam))

# proteins_plus_mHtt <- proteins_plus_mHtt[, !(names(proteins_plus_mHtt) %in% cols_to_remove)]
# top10_plus_mHtt <- top10_plus_mHtt[, !(names(top10_plus_mHtt) %in% cols_to_remove)]
# proteins_minus_mHtt <- proteins_minus_mHtt[, !(names(proteins_minus_mHtt) %in% cols_to_remove)]
# top10_minus_mHtt <- top10_minus_mHtt[, !(names(top10_minus_mHtt) %in% cols_to_remove)]

# Crear una lista con los dataframes y nombres de hoja
# lista_hojas <- list(
#   "proteins_plus_mHtt" = proteins_plus_mHtt,
#   "top10_plus_mHtt" = top10_plus_mHtt,
#   "proteins_minus_mHtt" = proteins_minus_mHtt,
#   "top10_minus_mHtt" = top10_minus_mHtt
# )

#------------------------------------------------------------------------------------------

calculate_foldchange <- function(df, 
                                 col_start_cond2 = col.start.cond2, col_end_cond2 = col.end.cond2, 
                                 col_start_cond1 = col.start.cond1, col_end_cond1 = col.end.cond1, 
                                 variable_name_cond2 = NameCond2, variable_name_cond1 = NameCond1, 
                                 foldchange_threshold = log2(1), ratio_threshold = 1) {
  
  # Calcular el promedio para las columnas cond2
  df[[paste0("mean_", variable_name_cond2)]] <- rowMeans(df[, col_start_cond2:col_end_cond2], na.rm = TRUE)
  
  # Calcular el promedio para las columnas cond1
  df[[paste0("mean_", variable_name_cond1)]] <- rowMeans(df[, col_start_cond1:col_end_cond1], na.rm = TRUE)
  
  # Calcular el fold change en log2
  df[[paste0("Log2_foldchange_", variable_name_cond1, "_", variable_name_cond2)]] <- 
    log2(df[[paste0("mean_", variable_name_cond1)]] / df[[paste0("mean_", variable_name_cond2)]])
  
  # Calcular el fold change sin log2
  df[[paste0("Ratio_", variable_name_cond1, "_", variable_name_cond1)]] <- 
    df[[paste0("mean_", variable_name_cond1)]] / df[[paste0("mean_", variable_name_cond2)]]
  
  # Determinar la significancia basada en el log2 fold change
  df[[paste0("significant_Log2_", variable_name_cond1)]] <- ifelse(
    df[[paste0("Log2_foldchange_", variable_name_cond2, "_", variable_name_cond1)]] > foldchange_threshold, 
    "IP", 
    ifelse(
      df[[paste0("Log2_foldchange_", variable_name_cond2, "_", variable_name_cond1)]] < foldchange_threshold, 
      "DR", 
      NA
    )
  )
  
  # Determinar la significancia basada en el ratio sin log2
  df[[paste0("significant_Ratio_", variable_name_cond2)]] <- ifelse(
    df[[paste0("Ratio_", variable_name_cond2, "_", variable_name_cond1)]] > ratio_threshold, 
    "IP", 
    ifelse(
      df[[paste0("Ratio_", variable_name_cond2, "_", variable_name_cond1)]] < (1 / ratio_threshold), 
      "DR", 
      NA
    )
  )
  
  # Contar los valores de IP y DR para ambas medidas
  count_IP_log2 <- sum(df[[paste0("significant_Log2_", variable_name_cond2)]] == "IP", na.rm = TRUE)
  count_DR_log2 <- sum(df[[paste0("significant_Log2_", variable_name_cond2)]] == "DR", na.rm = TRUE)
  
  count_IP_ratio <- sum(df[[paste0("significant_Ratio_", variable_name_cond2)]] == "DR", na.rm = TRUE)
  count_DR_ratio <- sum(df[[paste0("significant_Ratio_", variable_name_cond2)]] == "DR", na.rm = TRUE)
  
  # Guardar los resultados en un dataframe separado para el reporte
  summary_df <- data.frame(
    category = c("IP_Log2", "DR_Log2", "IP_Ratio", "DR_Ratio"),
    count = c(count_IP_log2, count_DR_log2, count_IP_ratio, count_DR_ratio)
  )
  
  # Devolver el dataframe con las nuevas columnas y el dataframe de resumen
  return(list(df = df, summary = summary_df))
}

# Aplicar la función a tu dataframe
result <- calculate_foldchange(list1)

# Obtener los datos principales (df) y el resumen de los conteos
df_ratio <- result$df
summary_df <- result$summary

# Filtrar y agregar las columnas adicionales como antes
protein_ids <- unique(top10_IP$Protein.ID)
result_filtered_IP <- subset(df_ratio, Protein.ID %in% protein_ids)


# Filtrar y agregar las columnas adicionales como antes
protein_ids <- unique(top10_DR$Protein.ID                                                )
result_filtered_DR <- subset(df_ratio, Protein.ID %in% protein_ids)
