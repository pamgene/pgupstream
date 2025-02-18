# methods for c-scoring

#' Function for simple c-scoring on standard upstream DB
#' @import dplyr
#' @export
cscore0 = function(db, rank_table = data.frame( rank = c(0,1,6,12), score = c(0.9, 0.7, 0.5, 0)) ){
  db  = db %>%
    distinct(ID, Kinase_Name,PepProtein_PhosLink,  Database, .keep_all = TRUE) %>%
    mutate(cscore = rank_table$score[1])

  for(i in 2:nrow(rank_table)){
    db = db %>%
      mutate(cscore = ifelse(Kinase_Rank >= rank_table$rank[i], rank_table$score[i], cscore))
  }
  dbw = db %>%
    group_by(ID, Kinase_Name) %>%
    dplyr::summarise(w = 1 - prod(1 - cscore)) %>%
    ungroup()

  pepcount = dbw %>%
    filter(w > 0) %>%
    group_by(ID) %>%
    dplyr::summarise(ks_count = n(),
              ks_mean_weight = mean(w),
              ks_sum_weight = sum(w)) %>%
    ungroup()

  dbw = dbw %>%
    left_join(pepcount, by = "ID") %>%
    mutate(wn = w/ks_sum_weight)
}
