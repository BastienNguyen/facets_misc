###################################################################
#### This function is to rank genomic event based on clonality ####
###################################################################
# work in progress nguyenb@mskcc.org
# Wed Jan 15 18:08:13 2020 ------------------------------


DoBradleyTerry = function(mat) {
  if(sum(colnames(mat) %in% c('SAMPLE_ID', 'Alteration', 'is_clonal')) != 3) {
    stop('The matrix needs three columns; SAMPLE_ID, Alteration, is_clonal')
  }
  require('BradleyTerryScalable')
  output_list = list()
  uniq_sampleID = table(mat$SAMPLE_ID)
  uniq_sampleID = names(uniq_sampleID[uniq_sampleID > 1])
  
  sample_kept = vector()
  input_matrix = list()
  for(i in 1:length(uniq_sampleID)) {
    tmp = mat[mat$SAMPLE_ID %in% uniq_sampleID[i],]
    tmp = tmp[order(tmp$is_clonal, decreasing = T),]
    tmp = tmp[!duplicated(tmp$Alteration),]                             # keep only the clonal event in case there is two concomitent alteration
    if(nrow(tmp) > 1) {                                                 # discard sample that have only >1 sample alteration
      sample_kept[i] = T
      tmp = tmp[sample(nrow(tmp)),]
      clonal_alt = tmp$Alteration[tmp$is_clonal]
      input_mat = data.frame(t(combn(tmp$Alteration, m = 2)))
      input_mat$X1 = as.character(input_mat$X1)
      input_mat$X2 = as.character(input_mat$X2)
      
      player1 = (input_mat$X1 %in% clonal_alt)*1
      player2 = (input_mat$X2 %in% clonal_alt)*1
      
      input_mat$outcome = NA
      input_mat$outcome[which(player1 > player2)] = 'W1'
      input_mat$outcome[which(player2 > player1)] = 'W2'
      input_mat$outcome[which(player2 == player1)] = 'D'
      input_matrix[[i]] = input_mat
    } else {
      sample_kept[i] = F
    }
  }
  sample_kept = uniq_sampleID[sample_kept]
  mat = mat[mat$SAMPLE_ID %in% sample_kept,]
  mat_clonal = mat[mat$is_clonal,]
  mat_subclonal = mat[!mat$is_clonal,]
  clonal_freq = table(mat_clonal$SAMPLE_ID, mat_clonal$Alteration)
  clonal_freq = (clonal_freq > 0)*1
  clonal_freq = colSums(clonal_freq)
  
  subclonal_freq = table(mat_subclonal$SAMPLE_ID, mat_subclonal$Alteration)
  subclonal_freq = (subclonal_freq > 0)*1
  subclonal_freq = colSums(subclonal_freq) 
  
  input_matrix = do.call(rbind, input_matrix)
  colnames(input_matrix)[1:2] = c('player1', 'player2')

  input_matrix_4col <- codes_to_counts(input_matrix, c("W1", "W2", "D"))
  output_list[[1]] = sample_kept
  output_list[[2]] = input_matrix_4col
  input_mat_btdata <- btdata(input_matrix_4col) 
  print(summary(input_mat_btdata))
  
  MLE <- btfit(input_mat_btdata, 1)
  MAP <- btfit(input_mat_btdata, 1.1)
  MLE = summary(MLE, SE = T)
  MAP = summary(MAP, SE = T)
  
  output_list[[3]] = MLE
  output_list[[4]] = MAP
  
  incidence_barplot = rbind(clonal_freq[MLE$item_summary$item],   subclonal_freq[MLE$item_summary$item])
  row.names(incidence_barplot) = c('clonal', 'subclonal')
  incidence_barplot[is.na(incidence_barplot)] = 0
  
  output_list[[5]] = incidence_barplot
  
  names(output_list) = c('input_samples','input_matrix', 'MLE', 'MAP', 'incidence_barplot_MLE')
  return(output_list)
}

plot_BradleyTerry = function(item_summary, incidence_barplot = NULL, main = NA) {
  require(RColorBrewer)
  item_summary$estimate_CI_low = item_summary$estimate - item_summary$SE
  item_summary$estimate_CI_hi = item_summary$estimate + item_summary$SE
  
  if(!is.null(incidence_barplot)) {
    par(mar = c(5.1/10, 0.5/10, 4.1/10, .5/10), oma = c(.5/10,5/10,.5/10,3/10))
    par(mar = c(5.1, 0.5, 4.1, .5), oma = c(.5,7,.5,3))
    layout(mat = matrix(c(1,1,2,3), nrow = 2, byrow = T), widths = c(10,3), heights = c(2,12))
    par(mar = c(1,1,1,1), oma = par('oma'), xpd=F)
    plot.new()
    text(0.5,0.2,main,cex=1,font=2)
    #par(mar = c(5.1, 0.5, 4.1, .5), oma = c(.5,5,.5,3))
    
    plot(NA, NA, ylim = c( 0, length(item_summary$estimate)), xlim = extendrange(c(item_summary$estimate_CI_low,item_summary$estimate_CI_hi), f = .05), axes = F, ylab = '', xlab = '', yaxs = 'i')
    
    abline(h = 1:nrow(item_summary), lty = 3, col = adjustcolor(col = 'grey', alpha.f = .5))
    col_alt = rep(brewer.pal(n = 3, name = 'Set1')[3], length(item_summary$item))
    col_alt[grepl(item_summary$item, pattern = '_gain')] = brewer.pal(n = 3, name = 'Set1')[1]
    col_alt[grepl(item_summary$item, pattern = 'Amplification')] = '#9f1214'
    col_alt[grepl(item_summary$item, pattern = 'Deletion')] = '#25567d'
    col_alt[grepl(item_summary$item, pattern = '_loss')] = brewer.pal(n = 3, name = 'Set1')[2]
    col_alt[grepl(item_summary$item, pattern = '_fusion')] = '#8a00c4'
    
    points(rev(item_summary$estimate), 0.5:length(item_summary$estimate), pch = 19, col = rev(col_alt))
    axis(side = 2, at = (length(item_summary$estimate)-.5):.5, labels = item_summary$item, las = 2, cex.axis = .75)
    segments(x0 = rev(item_summary$estimate_CI_low), x1 = rev(item_summary$estimate_CI_hi), y0 = .5:length(item_summary$estimate), y1 = .5:length(item_summary$estimate), col = rev(col_alt))
    
    barplot(colSums(incidence_barplot[,rev(1:ncol(incidence_barplot))]), horiz = T, names.arg = rep(NA, ncol(incidence_barplot)),  border = NA, yaxs = 'i', col = adjustcolor(rev(col_alt), alpha.f = 0.5), axes = F)
    barplot(incidence_barplot[1,rev(1:ncol(incidence_barplot))], horiz = T, names.arg = rep(NA, ncol(incidence_barplot)),  border = NA, yaxs = 'i', col = rev(col_alt), add = T, axes = F)
    axis(side = 3, line = 0.25, cex.axis = .75)
  } else {
    plot(NA, NA, ylim = c( 0, length(item_summary$estimate)), xlim = extendrange(item_summary$estimate, f = .05), axes = F, ylab = '', xlab = '', yaxs = 'i')
    abline(h = 1:nrow(item_summary), lty = 3, col = adjustcolor(col = 'grey', alpha.f = .5))
    col_alt = rep(brewer.pal(n = 3, name = 'Set1')[3], length(item_summary$item))
    col_alt[grepl(item_summary$item, pattern = '_gain')] = brewer.pal(n = 3, name = 'Set1')[1]
    col_alt[grepl(item_summary$item, pattern = '_loss')] = brewer.pal(n = 3, name = 'Set1')[2]
    col_alt[grepl(item_summary$item, pattern = '_fusion')] = '#8a00c4'
    
    
    points(rev(item_summary$estimate), 0.5:length(item_summary$estimate), pch = 19, col = rev(col_alt))
    axis(side = 2, at = (length(item_summary$estimate)-.5):.5, labels = item_summary$item, las = 2, cex.axis = .75)
    segments(x0 = rev(item_summary$estimate_CI_low), x1 = rev(item_summary$estimate_CI_hi), y0 = .5:length(item_summary$estimate), y1 = .5:length(item_summary$estimate), col = rev(col_alt))
    
  }
}


  
