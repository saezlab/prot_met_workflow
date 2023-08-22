plot_reaction_network <- function (network_and_attributes, t_table, scores_df, column_index, 
          rbrewer_plalette_name = "RdBu", vis.height = 700, vis.degree = 2, node_font_size = 20) 
{
  edges <- network_and_attributes[[1]]
  edges <- edges[edges[, 1] %in% scores_df[, 1] | edges[, 
                                                        2] %in% scores_df[, 1], ]
  nodes <- network_and_attributes[[2]]
  nodes <- nodes[nodes[, 1] %in% edges[, 1] | nodes[, 1] %in% 
                   edges[, 2], ]
  names(nodes) <- c("id", "molecule_type")
  metab_stat <- t_table[, c(1, column_index + 1)]
  names(metab_stat) <- c("id", "metab_stat")
  enzyme_score <- scores_df[, c(1, column_index + 1)]
  names(enzyme_score) <- c("id", "enzyme_score")
  nodes <- merge(nodes, metab_stat, all.x = T)
  nodes <- merge(nodes, enzyme_score, all.x = T)
  nodes$metab_stat[is.na(nodes$metab_stat)] <- 0
  nodes$enzyme_score[is.na(nodes$enzyme_score)] <- 0
  nodes$stat <- nodes$metab_stat + nodes$enzyme_score
  if (rbrewer_plalette_name == "RdBu") {
    nodes_enzymes <- nodes[nodes$molecule_type %in% c("reaction_enzyme", 
                                                      "transporter"), ]
    nodes_metabolites <- nodes[nodes$molecule_type == "metabolite", 
    ]
    nodes_enzymes$color <- ocean:::make_discrete_palette(nodes_enzymes$stat * 
                                                   -1, rbrewer_plalette_name)
    nodes_metabolites$color <- ocean:::make_discrete_palette(nodes_metabolites$stat * 
                                                       -1, rbrewer_plalette_name)
    nodes <- as.data.frame(rbind(nodes_enzymes, nodes_metabolites))
  }
  else {
    nodes_enzymes <- nodes[nodes$molecule_type %in% c("reaction_enzyme", 
                                                      "transporter"), ]
    nodes_metabolites <- nodes[nodes$molecule_type == "metabolite", 
    ]
    nodes_enzymes$color <- make_discrete_palette(nodes_enzymes$stat, 
                                                 rbrewer_plalette_name)
    nodes_metabolites$color <- make_discrete_palette(nodes_metabolites$stat, 
                                                     rbrewer_plalette_name)
    nodes <- as.data.frame(rbind(nodes_enzymes, nodes_metabolites))
  }
  nodes$color <- ifelse(nodes$stat == 0, "#C1C1C1", nodes$color)
  nodes$label <- nodes$id
  nodes$shape <- ifelse(nodes$molecule_type == "reaction_enzyme", 
                        "square", "triangle")
  nodes$shape <- ifelse(grepl("transporter", nodes$id), "circle", 
                        nodes$shape)
  nodes$label <- gsub("^transporter.*", "", nodes$label)
  names(edges) <- c("from", "to")
  edges$arrows <- "to"
  mapping_vec <- mapping_table$metab
  names(mapping_vec) <- mapping_table$KEGG
  nodes$label <- sapply(nodes$label, function(x, mapping_vec) {
    if (grepl("cpd:", x)) {
      suffixe <- stringr::str_extract(x, "_.$")
      x <- gsub("_.$", "", x)
      x <- gsub("cpd:", "", x)
      if (x %in% names(mapping_vec)) {
        x <- mapping_vec[x]
      }
      if (!is.na(suffixe)) {
        x <- paste0(x, suffixe)
      }
      return(x)
    }
    else {
      return(x)
    }
  }, mapping_vec = mapping_vec, simplify = F, USE.NAMES = F)
  nodes$font.size = node_font_size
  visNetwork::visNetwork(nodes = nodes, edges = edges, width = "100%", 
                         height = vis.height) %>% visOptions(nodesIdSelection = TRUE, 
                                                             highlightNearest = list(enabled = T, degree = vis.degree, 
                                                                                     hover = T))
}
