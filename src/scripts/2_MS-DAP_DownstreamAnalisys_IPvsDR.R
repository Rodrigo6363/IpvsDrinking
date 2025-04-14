library(ggrepel)
library(ggplot2)
library(msdap)
library(writexl)
library(openxlsx)
library(readxl)
library(dplyr)
library(ggvenn)
library(plotly)
library(svglite)

######### 13. Loading Inputs from MSDAP2####
prot.input<-read.delim(paste0(FolderName,DateTimeStamp,"/protein_abundance__input data as-is.tsv"))
dea<-readxl::read_excel(paste0(FolderName,DateTimeStamp,"/differential_abundance_analysis.xlsx"), sheet = 2)
colnames(prot.input) <- c("Protein.ID", "FASTA", "Gene.ID",
                           "5_Ctrl","3_DR","4_DR","1_IP","2_IP")


prot.input <- prot.input[, c("Protein.ID", "FASTA", "Gene.ID", "1_IP", "2_IP", "3_DR", "4_DR", "5_Ctrl")]

dea <- dea[,c(1:3,6,10,14,16)] #this doesn't change
colnames(dea) <- c("Protein.ID","accessions",                                         
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
ggsave(paste0(Folderfig,"/PCA.png"), plot = PCAplot, width = 8, height = 6, dpi = 300)
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
IntensitiesPlot <- recordPlot()
dev.off()
# barplot of IDs in GGPlot stylw
tmp <- data.frame(
  name=colnames(prot[,col.start.cond1:prot.col.end]) ,  #change columns to match sample numbers
  value=prot_count)

barplot_group <- ggplot(tmp, aes(x = name, y = value)) + 
  geom_bar(stat = "identity", colour = "black", 
           fill = c(rep("#DCCB4E", cond1), rep("#E98905", cond2), rep("#3A9AB2", control))) +
  scale_x_discrete(limits = colnames(prot.input[, col.start.cond1:prot.col.end])) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(colour = "black")) +
  ylab("protein groups") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(filename = paste0(Folderfig, "/barplot_group.svg"),
       plot = barplot_group,
       width = 8, height = 6, units = "in")

#----------------------------------------------------------------------------------------------------------
# List 1: Common proteins.
list1 <- dplyr::filter(prot.input, get(paste0("counter.vv.",NameCond1)) == 2 & get(paste0("counter.vv.",NameCond2)) == 2)
write_xlsx(list1, path = file.path(FolderList, paste0("List1_", NameCond1,"_", NameCond2, ".xlsx")))
# List 2: Enriched IP proteins.
list2.Fc1.5 <- list1[list1$Protein.ID %in% subset(dea, dea[[qvalue.colNam]] < pval & dea[[foldchange.colNam]] >= log2(foldlog2))$Protein.ID, ]
list2.Fc1.5_ordered <- merge(list2.Fc1.5, dea, by="Protein.ID")
list2.Fc1.5_ordered <- list2.Fc1.5_ordered %>% select(-accessions, -fasta_headers, -gene_symbols_or_id)
list2.Fc1.5_ordered <- list2.Fc1.5_ordered %>%
  arrange(desc(-foldchange.log2.IP.DR), desc(pvalue.log2.DR.IP))
write_xlsx(list2.Fc1.5, path = file.path(FolderList, paste0("List2.1_", NameCond2, "_FC1.5", ".xlsx")))
write_xlsx(list2.Fc1.5_ordered, path = file.path(FolderList, paste0("List2.1_", NameCond2, "_FC1.5_ordered", ".xlsx")))
# List 2.1: Enriched DR proteins.
list2.1Fc1.5 <- list1[list1$Protein.ID %in% subset(dea, dea[[qvalue.colNam]] < pval & dea[[foldchange.colNam]] <= -log2(foldlog2))$Protein.ID, ]
list2.1Fc1.5_ordered <- merge(list2.1Fc1.5, dea, by="Protein.ID")
list2.1Fc1.5_ordered <- list2.1Fc1.5_ordered %>% select(-accessions, -fasta_headers, -gene_symbols_or_id)
list2.1Fc1.5_ordered <- list2.1Fc1.5_ordered %>%
  arrange(desc(-foldchange.log2.IP.DR), desc(pvalue.log2.DR.IP))
write_xlsx(list2.1Fc1.5, path = file.path(FolderList, paste0("List2.1_", NameCond2, "_FC1.5", ".xlsx")))
write_xlsx(list2.1Fc1.5_ordered, path = file.path(FolderList, paste0("List2.1_", NameCond2, "_FC1.5_ordered", ".xlsx")))
# List 3: Unique IP proteins
list3<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] == 0)
write_xlsx(list3, path = file.path(FolderList, paste0("List3_", NameCond1, ".xlsx")))
# List 3.1: Unique DR proteins
list3.1<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond2)]] >=a & prot.input[[paste0("counter.vv.",NameCond1)]] == 0)
write_xlsx(list3.1, path = file.path(FolderList, paste0("List3.1_", NameCond2, ".xlsx")))
# List 4: IP Proteome
list4 <- rbind(list2.Fc1.5, list3)
write_xlsx(list4, path = file.path(FolderList, paste0("List4_", NameCond1, ".xlsx")))
# List 4.1: DR Proteome
list4.1 <- rbind(list2.1Fc1.5, list3.1)
write_xlsx(list4.1, path = file.path(FolderList, paste0("List4.1_", NameCond2, ".xlsx")))

#-----------------------------------------------------------------------------------------------------------


commons<-dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] >= b)
assign(paste0("commons.",NameCond1,"_", NameCond2),commons)
#write_xlsx(commons, path = paste0(FolderList,"Commons_",NameCond2, "_", NameCond1, ".xlsx"))
rm(commons)

uniques<-dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond1)]] >=a & prot.input[[paste0("counter.vv.",NameCond2)]] == 0)
assign(paste0("uniques.",NameCond1),uniques)
#write_xlsx(uniques, path = paste0(FolderList,"List3", NameCond1, ".xlsx"))
rm(uniques)
uniques<- dplyr::filter(prot.input, prot.input[[paste0("counter.vv.",NameCond2)]] >=b & prot.input[[paste0("counter.vv.",NameCond1)]] == 0)
assign(paste0("uniques.",NameCond2),uniques)
#write_xlsx(uniques, path = paste0(FolderList,"List3_", NameCond1, ".xlsx"))
rm(uniques)

#-------------------------------------------------------------------------------------------------------------

# Define los sets de proteínas detectadas en cada condición
set_IP <- prot.input$Protein.ID[prot.input[[paste0("counter.vv.", NameCond1)]] == 2]
set_DR <- prot.input$Protein.ID[prot.input[[paste0("counter.vv.", NameCond2)]] == 2]

venn_data <- list(
  IP = set_IP,
  DR = set_DR
)

venn_plot <- ggvenn(venn_data,
                    fill_color = c("lightblue", "salmon"),
                    stroke_size = 0.5,
                    set_name_size = 8,
                    text_size = 5) +
  ggtitle("Venn diagram detected proteins") +
  
  # Añadir anotaciones debajo de cada área
  annotate("text", x = -1.2, y = -0.5, label = "list3", size = 7, fontface = "italic") +       # IP
  annotate("text", x =  1.2, y = -0.5, label = "list3.1", size = 7, fontface = "italic") +     # DR
  annotate("text", x =  0, y =  0.5, label = "list1", size = 7, fontface = "italic") +           # Intersección
  
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b=20)),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

ggsave(paste0(Folderfig,"/Venn_diagram.svg"), plot = venn_plot, width = 8, height = 6, dpi = 300)

# Mostrar el gráfico
print(venn_plot)

venn_data_com <- list(
  IP = list2.Fc1.5$Protein.ID,
  DR = list2.1Fc1.5$Protein.ID
)

venn_plot_com <- ggvenn(venn_data_com,
                    fill_color = c("lightblue", "salmon"),
                    stroke_size = 0.5,
                    set_name_size = 8,
                    text_size = 5) +
  ggtitle("Venn diagram detected proteins") +
  
  # Añadir anotaciones debajo de cada área
  annotate("text", x = -1.2, y = -0.5, label = "list2", size = 7, fontface = "italic") +       # IP
  annotate("text", x =  1.2, y = -0.5, label = "list2.1", size = 7, fontface = "italic") +     # DR
  #annotate("text", x =  0, y =  0.5, label = "list1", size = 7, fontface = "italic") +        # Intersección
  
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b=20)),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

print(venn_plot_com)
#-------------------------------------------------------------------------------------------------------------
# 1. Calcular los promedios de IP y DR
df_scatter <- prot.input %>%
  rowwise() %>%
  mutate(
    mean_IP = mean(c_across(contains("_IP")), na.rm = TRUE),
    mean_DR = mean(c_across(contains("_DR")), na.rm = TRUE)
  ) %>%
  ungroup()

# 2. Scatterplot de IP vs DR con ejes definidos, sin grid, puntos redondeados y mismo color para todos
scatter <-  ggplot(df_scatter, aes(x = mean_IP, y = mean_DR)) +
    geom_point(size = 3, shape = 21, fill = "#1f77b4", stroke = 0.3) +
    labs(
      title = "Scatterplot: IP vs DR",
      x = "Mean IP intensity",
      y = "Mean DR intensity"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black")
    ) +
    expand_limits(x = 0, y = 0)

ggsave(paste0(Folderfig,"/scatter.svg"), plot = scatter, width = 8, height = 6, dpi = 300)


print(scatter)
#-------------------------------------------------------------------------------------------------------------
calcular_media_cond <- function(df, sufijo_cond) {
  df %>%
    rowwise() %>%
    mutate(mean = mean(c_across(contains(sufijo_cond)), na.rm = TRUE)) %>%
    pull(mean) %>%
    mean(na.rm = TRUE)
}

mean_list2   <- calcular_media_cond(list2.Fc1.5, "_IP")
mean_list3   <- calcular_media_cond(list3, "_IP")
mean_list4   <- calcular_media_cond(list4, "_IP")
mean_list2.1 <- calcular_media_cond(list2.1Fc1.5, "_DR")
mean_list3.1 <- calcular_media_cond(list3.1, "_DR")
mean_list4.1 <- calcular_media_cond(list4.1, "_DR")


# Crear dataframe para el plot
bar_data <- data.frame(
  grupo = factor(c("list2", "list3", "list4", "list2.1", "list3.1", "list4.1"),
                 levels = c("list2", "list3", "list4", "list2.1", "list3.1", "list4.1")),  # orden fijo
  condicion = c("IP", "IP", "IP", "DR", "DR", "DR"),
  media = c(mean_list2, mean_list3, mean_list4, mean_list2.1, mean_list3.1, mean_list4.1)
)


# Gráfico de barras agrupado sin grid y con ejes bien marcados
int_list <- ggplot(bar_data, aes(x = grupo, y = media, fill = condicion)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("IP" = "#1f77b4", "DR" = "#ff7f0e")) +
    labs(
      title = "Mean intensity per group",
      x = "Group",
      y = "Mean intensity",
      fill = "Condition"
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

print(int_list)
ggsave(paste0(Folderfig,"/protein_intensity_lists.svg"), plot = int_list, width = 8, height = 6, dpi = 300)


#-------------------------------------------------------------------------------------------------------------
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
assign(paste0("dea.", NameCond2, "vs", NameCond1), dea_data)
rm(dea_data)
######### 14c. Generating Candidate Lists based on q-value and log############
# Initialize candidates column
dea_data <- get(paste0("dea.", NameCond2, "vs", NameCond1))
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
                            " | Enriched proteins: DR (", count_negative, 
                            ") / IP (", count_positive, ")")

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
                                       "Non-enriched proteins")),
                label = gene_symbols_or_id  # <- para tooltip interactivo
            )) + 
  geom_point(size = 0.9) +
  scale_colour_manual(name = NULL,
                      values = c("IP" = '#003366', "DR" = '#FFA500', "Non-enriched proteins" = 'grey')) +
  theme_classic() + 
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_text(colour = "black")) +
  ggtitle(title_with_counts) + 
  xlab("log2 fold change") + 
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(pval), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_vline(xintercept = log2(foldlog2), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_vline(xintercept = -log2(foldlog2), linetype = "dashed", color = "grey50", size = 0.6) +
  geom_label_repel(data = enriched_proteins,
                   aes(x = get(foldchange.colNam),
                       y = -log10(get(pvalue.colNam)),
                       label = gene_symbols_or_id),
                   size = 3,
                   color = "black",
                   fill = "lightgrey",
                   box.padding = 0.5,
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   inherit.aes = FALSE)

  guides(colour = guide_legend(nrow = 1))

DEAvolplot<-(C)
print(C)
rm(C)

ggsave(paste0(Folderfig,"/volcano.svg"), plot = DEAvolplot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------------

data_clean$log2fc <- data_clean[[foldchange.colNam]]
data_clean$log10pval <- -log10(data_clean[[pvalue.colNam]])

data_clean$grupo <- ifelse(
  data_clean[[pvalue.colNam]] <= pval & data_clean[[foldchange.colNam]] >= log2(foldlog2),
  "IP",
  ifelse(
    data_clean[[pvalue.colNam]] <= pval & data_clean[[foldchange.colNam]] <= -log2(foldlog2),
    "DR",
    "Non-enriched proteins"
  )
)

fig <- plot_ly(
  data = data_clean,
  x = ~log2fc,
  y = ~log10pval,
  type = 'scatter',
  mode = 'markers',
  color = ~grupo,
  colors = c("IP" = '#003366', "DR" = '#FFA500', "Non-enriched proteins" = 'grey'),
  text = ~paste("Protein ID:", Protein.ID,
                "<br>Gene:", gene_symbols_or_id,
                "<br>log2FC:", round(log2fc, 2),
                "<br>-log10(p):", round(log10pval, 2)),
  hoverinfo = 'text',
  marker = list(size = 5, opacity = 0.8)
)


fig <- fig %>%
  layout(
    title = list(text = title_with_counts),
    xaxis = list(title = "log2 fold change"),
    yaxis = list(title = "-log10(p-value)"),
    legend = list(orientation = "h", x = 0.3, y = -0.2),
    margin = list(b = 80)
  )


htmlwidgets::saveWidget(fig, "plotly_plot.html", selfcontained = TRUE)
browseURL("plotly_plot.html")
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
  arrange(desc(get(foldchange.colNam))) %>%
  head(10)

top10_DR <- proteins_DR %>%
  arrange(!!sym(foldchange.colNam)) %>%  # orden ascendente → más negativo primero
  head(10)

# Verificar los nuevos dataframes
print("Top 10 proteínas enriquecidas en IP:")
print(top10_IP)

print("Top 10 proteínas enriquecidas en DR:")
print(top10_DR)

top10_IP <- rename(top10_IP,"FASTA" = "fasta_headers",
                             "Gene.ID" = "gene_symbols_or_id")
top10_IP <- merge(top10_IP, list2.Fc1.5,by="Protein.ID")
top10_IP <- select(top10_IP, -"accessions", -paste0("counter.vv.",NameCond2), -paste0("counter.vv.",NameCond1),
                             -paste0(foldchange.colNam), -paste0(pvalue.colNam), -paste0(qvalue.colNam))

top10_DR <- rename(top10_DR,"FASTA" = "fasta_headers",
                   "Gene.ID" = "gene_symbols_or_id")
top10_DR <- merge(top10_DR, list2.1Fc1.5,by="Protein.ID")
top10_DR <- select(top10_DR, -"accessions", -paste0("counter.vv.",NameCond2), -paste0("counter.vv.",NameCond1),
                              -paste0(foldchange.colNam), -paste0(pvalue.colNam), -paste0(qvalue.colNam))


#Crear una lista con los dataframes y nombres de hoja
lista_hojas <- list(
  "proteins_IP" = proteins_IP,
  "top_IP" = top10_IP,
  "proteins_DR" = proteins_DR,
  "top10_DR" = top10_DR
)

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
  df[[paste0("Ratio_", variable_name_cond1, "_", variable_name_cond2)]] <- 
    df[[paste0("mean_", variable_name_cond1)]] / df[[paste0("mean_", variable_name_cond2)]]
  
  # Determinar la significancia basada en el log2 fold change
  df[[paste0("significant_Log2_", variable_name_cond1)]] <- ifelse(
    df[[paste0("Log2_foldchange_", variable_name_cond1, "_", variable_name_cond2)]] > foldchange_threshold, 
    "IP", 
    ifelse(
      df[[paste0("Log2_foldchange_", variable_name_cond1, "_", variable_name_cond2)]] < foldchange_threshold, 
      "DR", 
      NA
    )
  )
  
  # Determinar la significancia basada en el ratio sin log2
  df[[paste0("significant_Ratio_", variable_name_cond1)]] <- ifelse(
    df[[paste0("Ratio_", variable_name_cond1, "_", variable_name_cond2)]] > ratio_threshold, 
    "IP", 
    ifelse(
      df[[paste0("Ratio_", variable_name_cond1, "_", variable_name_cond2)]] < (1 / ratio_threshold), 
      "DR", 
      NA
    )
  )
  
  # Contar los valores de IP y DR para ambas medidas
  count_IP_log2 <- sum(df[[paste0("significant_Log2_", variable_name_cond1)]] == "IP", na.rm = TRUE)
  count_DR_log2 <- sum(df[[paste0("significant_Log2_", variable_name_cond1)]] == "DR", na.rm = TRUE)
  
  count_IP_ratio <- sum(df[[paste0("significant_Ratio_", variable_name_cond1)]] == "IP", na.rm = TRUE)
  count_DR_ratio <- sum(df[[paste0("significant_Ratio_", variable_name_cond1)]] == "DR", na.rm = TRUE)
  
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
protein_ids <- unique(top10_DR$Protein.ID)
result_filtered_DR <- subset(df_ratio, Protein.ID %in% protein_ids)


# Crear una lista con los dataframes
dataframes <- list(
  "df_ratio" = df_ratio, 
  "summary" = summary_df, 
  "volcano_IP" = result_filtered_IP, 
  "volcano_DR" = result_filtered_DR
)

# Filtrar dataframes vacíos (con 0 filas)
dataframes <- dataframes[sapply(dataframes, function(df) nrow(df) > 0)]

# Guardar en Excel solo si la lista no está vacía
if (length(dataframes) > 0) {
  write_xlsx(dataframes, path = paste0(FolderList,"df_ratio.xlsx"))
  print("Archivo Candidates.xlsx guardado correctamente.")
} else {
  print("No hay dataframes válidos para guardar.")
}

# Mensaje en consola para imprimir el número de +mHtt y -mHtt para ambas condiciones
message("Number of IP (Log2): ", summary_df$count[summary_df$category == "IP_Log2"])
message("Number of DR (Log2): ", summary_df$count[summary_df$category == "DR_Log2"])
message("Number of IP (Ratio): ", summary_df$count[summary_df$category == "IP_Ratio"])
message("Number of DR (Ratio): ", summary_df$count[summary_df$category == "DR_Ratio"])


######### 16. Printing all plots to PDF############
pdf(paste0(Folderfig,"/",NameCond1,"vs",NameCond2,"Downstream_Results.pdf"))
PCAplot
# IntensitiesPlot
GroupsPlot_recorded
venn_plot
scatter
DEAvolplot
dev.off()
