calculate_top1 <- function(ranked_list, true_list){
  
  order1 <- order(true_list)
  order2 <- order(ranked_list)
  mean(order2[1] %in% order1[1])
}

calculate_top3 <- function(ranked_list, true_list){
  
  order1 <- order(true_list)
  order2 <- order(ranked_list)
  mean(order2[1:3] %in% order1[1:3])
}