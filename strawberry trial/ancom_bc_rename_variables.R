#SEs
tab_se = res$se
colnames(tab_se) = col_name

#Test statistics
tab_w = res$W
colnames(tab_w) = col_name

#P-values
tab_p = res$p_val
colnames(tab_p) = col_name

#Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name

#Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name

#Bias-adjusted abundances
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 


