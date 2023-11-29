# From dataframe with counts for state and subgroup variables get contingency table elements
dfContin <- function(df) {
  box::use(dplyr[...])
  dfc <- df %>% 
    left_join(df %>% group_by(state) %>% summarise(t=sum(n)),by=c("state")) %>% mutate(ny=t-n) %>% select(-t) %>% 
    left_join(df %>% group_by(subgroup) %>% summarise(t=sum(n)),by=c("subgroup")) %>% mutate(yn=t-n) %>% select(-t) %>%
    mutate(t=sum(df %>% pull(n))) %>% mutate(nn=t-ny-yn-n) %>% select(-t) %>%
    rename(yy=n)
  return(dfc)
}