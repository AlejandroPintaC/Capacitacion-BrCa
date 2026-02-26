library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GO.db)
library(ggplot2)
library(AnnotationDbi) 
library(enrichplot)
url_raw <- "https://raw.githubusercontent.com/guillermodeandajauregui/WorkshopAdvancedBioinformatics2021/main/data/BasalBreastCancer.graphml"

g <- read_graph(url_raw, format = "graphml")


# Comunidades Louvain
comunidades <- cluster_louvain(g)
V(g)$comunidad <- membership(comunidades)


# Usar AnnotationDbi::select() para consultar la base de datos GO
go_mf <- AnnotationDbi::select(org.Hs.eg.db, 
                                keys = keys(org.Hs.eg.db), 
                                columns = c("GO", "ONTOLOGY"), 
                                keytype = "ENTREZID") %>%
  filter(ONTOLOGY == "MF") %>%           
  dplyr::select(term_id = GO, gene_id = ENTREZID) %>%  # Usar dplyr::select() para dataframes
  distinct()

# Verificar que funcionó
head(go_mf)

go_names <- AnnotationDbi::select(GO.db,
                                  keys = unique(go_mf$term_id),
                                  columns = "TERM",
                                  keytype = "GOID") %>%
  distinct()

term2name <- go_names %>%
  dplyr::rename(term_id = GOID, term_name = TERM)


# 5. Enriquecimiento por comunidad
resultados_lista <- list()

for (comm in sort(unique(V(g)$comunidad))) {
  
  genes_comunidad <- V(g)$name[V(g)$comunidad == comm]
  
  # Convertir a Entrez si son símbolos
  if (!all(grepl("^[0-9]+$", genes_comunidad))) {
    conv <- suppressWarnings(bitr(genes_comunidad, fromType = "SYMBOL", 
                                  toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
    if (!is.null(conv)) {
      genes_comunidad <- conv$ENTREZID
    } else {
      genes_comunidad <- genes_comunidad[grepl("^[0-9]+$", genes_comunidad)]
    }
  }
  
  universe <- V(g)$name
  if (!all(grepl("^[0-9]+$", universe))) {
    conv_univ <- suppressWarnings(bitr(universe, fromType = "SYMBOL", 
                                       toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
    if (!is.null(conv_univ)) {
      universe <- unique(conv_univ$ENTREZID)
    } else {
      universe <- universe[grepl("^[0-9]+$", universe)]
    }
  }
  
  enr <- enricher(gene = genes_comunidad,
                  TERM2GENE = go_mf,
                  TERM2NAME = term2name,
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  universe = universe,
                  minGSSize = 5,
                  maxGSSize = 500)
  
  if (!is.null(enr) && nrow(enr) > 0) {
    resultados_lista[[paste0("Comunidad_", comm)]] <- enr
    cat("Comunidad", comm, ":", nrow(enr), "términos\n")
  }
}

# 6. Exportar resultados
dir.create("resultados_enriquecimiento", showWarnings = FALSE)
for (nombre in names(resultados_lista)) {
  df_result <- as.data.frame(resultados_lista[[nombre]]) %>%
    arrange(p.adjust)
  write.csv(df_result, 
            file = file.path("resultados_enriquecimiento", paste0(nombre, "_MF.csv")), 
            row.names = FALSE)
}

#### Graficar ####
### Comunidad 1 ###
enr_com1 <- resultados_lista[["Comunidad_1"]]


# Dotplot básico (top 20 términos por FDR)
p <- dotplot(enr_com1, 
        showCategory = 20, 
        title = "Comunidad 1 - GO Molecular Function",
        font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # Fondo blanco [[6]]
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = file.path("/Users/alejandropintacastro/resultados_enriquecimiento", 
                            paste0(nombre, "_dotplot.png")),
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Comunidad 2

enr_com2 <- resultados_lista[["Comunidad_2"]]

p <- dotplot(enr_com2, 
             showCategory = 20, 
             title = "Comunidad 2 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_2_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

#Comunidad 3

enr_com3 <- resultados_lista[["Comunidad_3"]]

p <- dotplot(enr_com3, 
             showCategory = 20, 
             title = "Comunidad 3 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_3_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Comunidad 5

enr_com5 <- resultados_lista[["Comunidad_5"]]

p <- dotplot(enr_com5, 
             showCategory = 20, 
             title = "Comunidad 5 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_5_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Comunidad 7

enr_com7 <- resultados_lista[["Comunidad_7"]]

p <- dotplot(enr_com7, 
             showCategory = 20, 
             title = "Comunidad 7 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_7_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Comunidad 8

enr_com8 <- resultados_lista[["Comunidad_8"]]

p <- dotplot(enr_com8, 
             showCategory = 20, 
             title = "Comunidad 8 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_8_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

#Comunidad 9
enr_com9 <- resultados_lista[["Comunidad_9"]]

p <- dotplot(enr_com9, 
             showCategory = 20, 
             title = "Comunidad 9 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_9_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

#Comunidad 10
enr_com10 <- resultados_lista[["Comunidad_10"]]

p <- dotplot(enr_com10, 
             showCategory = 20, 
             title = "Comunidad 10 - GO Molecular Function",
             font.size = 10) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right"
  )

ggsave(filename = "/Users/alejandropintacastro/resultados_enriquecimiento/Comunidad_10_dotplot.png",
       plot = p,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")
