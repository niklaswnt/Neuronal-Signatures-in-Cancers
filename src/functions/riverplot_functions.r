scaleRiverplot = function(rp, min_rank = 3, max_rank = 4) {
  # function: scales the flows of a riverplot object obtained from the generateRiverplot function of the "riverplot" package

  # Input:
    # rp: the riverplot object
    # min_rank: the minimum rank of the riverplot object
    # max_rank: the maximum rank of the riverplot object

  # Output:
    # rp: the riverplot object with the rescaled flows

  for (sig in 1:min_rank) {
  outflow = rp$edges %>% filter(N1 == paste0(min_rank, "_S", sig))

  outflow_val = outflow[,"Value"]

  rp$edges[which(rp$edges$N1 == paste0(min_rank, "_S", sig)),"rescaled"] = rp$edges[which(rp$edges$N1 == paste0(min_rank, "_S", sig)),"Value"]

  }

  for (k in (min_rank + 1):max_rank) {
      for (sig in 1:k) {
          inflow = rp$edges %>% filter(N2 == paste0(k, "_S", sig)) 
          outflow = rp$edges %>% filter(N1 == paste0(k, "_S", sig))

          inflow_val = inflow[,"rescaled"]
          outflow_val = outflow[,"rescaled"]

          inflow
          outflow

          rp$edges[which(rp$edges$N1 == paste0(k, "_S", sig)),"rescaled"] = outflow_val / sum(outflow_val) * sum(inflow_val)
      }
  }

  colnames(rp$edges) = c("N1", "N2", "ignore", "ID", "Value")

  return(rp)
}