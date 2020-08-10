# ________________________________________________________________________
# Comparison of PLSC weights and PLSC-RP weights using similarity measures
# ________________________________________________________________________

# input:
# ------
# w_PLSC = weights of PLSC
# w_PLSC_RP = weights of PLSC-RP

# output:
# -------
# sim_corr = Pearson correlation
# sim_cos = Cosine measure
# sim_eJac = extended Jaccard similarity


library(proxy)

similarity <- function(w_PLSC, w_PLSC_RP){
  
  # check direction of weights: multiply with (-1) if reversed (necessary for correct eJaccard)
  if(cor(w_PLSC, w_PLSC_RP) < 0){
    w_PLSC_RP <- (-1) * w_PLSC_RP
  }
  
  sim_corr <- simil(w_PLSC, w_PLSC_RP, method="correlation", by_rows=F)
  sim_cos <- simil(w_PLSC, w_PLSC_RP, method="cosine", by_rows=F)
  sim_eJac <- simil(w_PLSC, w_PLSC_RP, method="eJaccard", by_rows=F)
  
  return(list(sim_corr=sim_corr, sim_cos=sim_cos, sim_eJac=sim_eJac))
}
