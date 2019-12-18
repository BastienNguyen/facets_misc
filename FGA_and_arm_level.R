###########################
#### nguyenb@mskcc.org ####
###########################
# Tue Apr 16 12:41:16 2019 ------------------------------
# Fri Jun 21 10:07:02 2019 ------------------------------ ADDED N_BREAKPOINTS and INTERGRAL PLOIDY
# Fri Jun 21 11:52:33 2019 ------------------------------ ADDED get_ARM_FGA
# Wed Dec 18 12:09:06 2019 ------------------------------ TO ADD calculate aneuploidy score (defined by the number of arm alteration) 
# Wed Dec 18 12:09:40 2019 ------------------------------ TO ADD try to call arm level event

load('~/Documents/DATA/Other/arm_position.Rdata')

# function to calculate the size of each integral copy
get_intergral_size = function(fit, method = 'em') {
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  if(method == 'cncf') {
    tt = table(cncf$tcn)
    }
  else {
    tt = table(cncf$tcn.em)
    }
  if(length(tt) == 1) {
    integral_size = sum(cncf$size, na.rm = T)
    names(integral_size) = names(tt)
  } else {
    tt2 = names(tt)
    integral_size=vector()
    for(y in 1:length(tt2)){
      if(method == 'cncf') {
        integral_size[y] = sum(cncf$size[which(cncf$tcn == tt2[y])])
      } else {
        integral_size[y] = sum(cncf$size[which(cncf$tcn.em == tt2[y])])
      }
    }
    names(integral_size) = tt2
  }
  names(integral_size) = as.numeric(names(integral_size))
  return(integral_size)
}

# main function to calculate FGA / GAIN / LOSS / LOH / WGD
get_FGA = function(facets_rdata, out = NULL, sampleID = NULL, method = 'em', include_loh = F) {
  load(facets_rdata)
  cncf = fit$cncf
  cncf$size = cncf$end - cncf$start
  integral_size = get_intergral_size(fit, method)
  cut_off = as.numeric(names(integral_size)) > 2
  is_WGD = sum(integral_size[cut_off], na.rm = T)/(sum(integral_size, na.rm = T)) > .5
  major_cn = as.numeric(names(which.max(integral_size)))
  gain = sum(integral_size[as.numeric(names(integral_size)) > major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  loss = sum(integral_size[as.numeric(names(integral_size)) < major_cn], na.rm = T)/sum(integral_size, na.rm = T)
  if(method == 'cncf') {
    LOH = sum(cncf$size[which(cncf$tcn == major_cn & cncf$lcn == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  } else {
    LOH = sum(cncf$size[which(cncf$tcn.em == major_cn & cncf$lcn.em == 0)], na.rm = T)/sum(integral_size, na.rm = T)
  }
  if(include_loh) {
    FGA = gain + loss + LOH
  } else {
    FGA = gain + loss
  }
  if(!is.null(out)) {
    sample_name = as.character(out$IGV[1,1])  
  } else {
    sample_name = sampleID
  }
  if(is.null(fit$emflags)) { fit$emflags = '' }
  ouput = data.frame('SAMPLE' = sample_name, 'PURITY' = fit$purity, 'PLOIDY' = fit$ploidy, 'FGA' = FGA, 'GAIN' = gain, 'LOSS'= loss, 'LOH' = LOH, 'N_BREAKPOINTS' = nrow(cncf), 'IS_WGD' = is_WGD, 'INTERGAL_PLOIDY' = major_cn, 'EM_FLAGS' = fit$emflags)
  return(ouput)
}

# function to compute FGA by chromosomal arm # need arm_position.Rdata 
get_ARM_FGA = function(facets_rdata, arm_position, method = 'em', include_loh = F, calls_threshold = 0.8) {
  load(facets_rdata)
  require(data.table)
  sample_name = as.character(out$IGV[1,1])
  integral_size = get_intergral_size(fit, method)
  major_cn = as.numeric(names(which.max(integral_size)))
  
  cncf = fit$cncf
  cncf_pos = data.table(chr=cncf$chrom, loc.start=cncf$start,
                        loc.end=cncf$end)
  setkey(cncf_pos, chr, loc.start, loc.end)
  arm_position = arm_position[ !arm_position$arm %in% c('13p', '14p', '15p', '21p', '22p','23p', '23q')] # remove acrocentric and sex chromosome
  
  setkey(arm_position, chr, start, end)
  fo_impact.idx <- foverlaps(arm_position, cncf_pos, nomatch=NA, which = T)
  fo_impact = foverlaps(arm_position, cncf_pos, nomatch=NA)
  fo_impact[, loc.start := ifelse(loc.start < start, start, loc.start)]
  fo_impact[, loc.end := ifelse(loc.end > end, end, loc.end)]
  fo_impact$size = fo_impact$loc.end - fo_impact$loc.start
  
  cncf = cncf[fo_impact.idx$yid,]
  cncf$arm = fo_impact$arm
  cncf$size = fo_impact$size
  
  n_breaks = cf_gain = gain = cf_loss = loss = vector()
  for(k in 1:length(arm_position$arm)) {
    tmp = cncf[which(cncf$arm == arm_position$arm[k]),]
    if( sum(is.na(tmp$size)) ==1 ) {
      n_breaks[k] = gain[k] = loss[k] = NA
    } else {
      n_breaks[k] = nrow(tmp)
      if(method == 'cncf') {
        gain[k] = sum(tmp$size[which(tmp$tcn > major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_gain_tmp =  if(length(tmp$cf[which(tmp$tcn > major_cn)]) == 0) { NA } else { tmp$cf[which(tmp$tcn > major_cn)]}
        cf_gain[k] =   if(length(cf_gain_tmp) == 1) { cf_gain_tmp} else { cf_gain_tmp[which.max(tmp$size[which(tmp$tcn > major_cn)])] }
        loss[k] = sum(tmp$size[which(tmp$tcn < major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_loss_tmp =  if(length(tmp$cf[which(tmp$tcn < major_cn)]) == 0) { NA } else { tmp$cf[which(tmp$tcn < major_cn)]}
        cf_loss[k] = if(length(cf_loss_tmp) == 1) { cf_loss_tmp } else { cf_loss_tmp[which.max(tmp$size[which(tmp$tcn < major_cn)])] } 
        
      } else {
        gain[k] = sum(tmp$size[which(tmp$tcn.em > major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_gain_tmp =  if(length(tmp$cf.em[which(tmp$tcn.em > major_cn)]) == 0) { NA } else { tmp$cf.em[which(tmp$tcn.em > major_cn)]}
        cf_gain[k] = if(length(cf_gain_tmp) == 1) { cf_gain_tmp} else { cf_gain_tmp[which.max(tmp$size[which(tmp$tcn.em > major_cn)])] }
        loss[k] = sum(tmp$size[which(tmp$tcn.em < major_cn)], na.rm = T)/sum(tmp$size, na.rm = T)
        cf_loss_tmp =  if(length(tmp$cf.em[which(tmp$tcn.em < major_cn)]) == 0) { NA } else { tmp$cf.em[which(tmp$tcn.em < major_cn)]}
        cf_loss[k] = if(length(cf_loss_tmp) == 1) { cf_loss_tmp } else { cf_loss_tmp[which.max(tmp$size[which(tmp$tcn.em < major_cn)])] } 
        
      }
    }
  }
  
  altered_gain = if(sum(gain > calls_threshold) == 0) {
    NA
  } else {
    paste0(arm_position$arm[which(gain > calls_threshold)], '_gain')
  }
  
  altered_loss = if(sum(loss > calls_threshold) == 0) {
    NA
  } else {
    paste0(arm_position$arm[which(loss > calls_threshold)], '_loss')
  }
  altered_arm = c(altered_gain, altered_loss)
  altered_arm = altered_arm[which(!is.na(altered_arm))]
  if(length(altered_arm) > 0) {
    output = data.frame('SAMPLE_ID' = sample_name, 'altered_arm' = altered_arm,
                        'altered_arm_cf' = c(as.numeric(cf_gain[which(gain > calls_threshold)]), as.numeric(cf_loss[which(loss > calls_threshold)])), 'purity' = as.numeric(fit$purity))
    output$altered_arm_ccf = output$altered_arm_cf/output$purity
    return(output)
  } 
}
