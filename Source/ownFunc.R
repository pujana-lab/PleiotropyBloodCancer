# From dataframe with counts for state and subgroup variables get contingency table elements
dfContin <- function(df) {
  library(dplyr)
  dfc <- df %>% 
    left_join(df %>% group_by(state) %>% summarise(t=sum(n)),by=c("state")) %>% mutate(ny=t-n) %>% select(-t) %>% 
    left_join(df %>% group_by(subgroup) %>% summarise(t=sum(n)),by=c("subgroup")) %>% mutate(yn=t-n) %>% select(-t) %>%
    mutate(t=sum(df %>% pull(n))) %>% mutate(nn=t-ny-yn-n) %>% select(-t) %>%
    rename(yy=n)
  return(dfc)
}
# Diagram venn plot
vennDiagramPlot <- function(df,p) {
  library(ggVennDiagram)
  library(ggplot2)
  venn.data <- list(pl = unlist(df$pleio.gene),
                    ph = unlist(df$pheno.gene))
  plot <- ggVennDiagram(venn.data, label_alpha = 0,category.names = c("Pleiotropic",p),label="count") +
    scale_fill_continuous(low="white",high="white") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("OR = ",round(df$OR,2),", p-val = ",round(df$pval,3)))
  return(plot)
}