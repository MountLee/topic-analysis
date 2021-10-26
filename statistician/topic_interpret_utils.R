##### the function to search id of a given list of name elements #####
search_name <- function(name,name_list){
  potential <- c()
  for (i in 1:length(name_list)){
    x <- name_list[i]
    tf <- TRUE
    for (nm in name){
      tf <- tf & grepl(nm, x)
    }
    if (tf){
      potential <- append(potential,i)
    }
  }
  return(potential)
}

##### the function to get average topic weight of an author #####
get_pap_tw_author <- function(id,nm_vector = 0){
  pap_list <- unlist(author_pap_list[id])
  pap_list <- paper$wos[pap_list]
  Npaper <- length(pap_list)
  
  tf_label <- paper_time$wos %in% pap_list
  W_hat_interest <- as.matrix(W_hat.normed[,tf_label[valid_ix]])
  
  if (dim(W_hat_interest)[2] > 0){
    av_tw <- rowSums(W_hat_interest) / dim(W_hat_interest)[2]
    av_tw <- av_tw / sum(av_tw)
    av_tw <- av_tw - nm_vector
  }
  else{
    av_tw <- rep(0,dim(W_hat.normed)[1])
  }
  
  return(av_tw)
}

get_tw_author <- function(name_or_id, topics, use_name = TRUE, name_title = NULL, 
                          normal = 'divid', nm_vector = 0,save = TRUE,no_ylabel = FALSE,horiz = TRUE,topic_names = FALSE){
  if (use_name){
    x <- which(author_name == name_or_id)
    pap_list <- author_pap_list[[x]]
    if (is.null(name_title)){
      name_title <- name_or_id
    }
  }
  else{
    pap_list <- unlist(author_pap_list[name_or_id])
    pap_list <- paper$wos[pap_list]
    if (is.null(name_title)){
      name_title <- author_name[name_or_id]
    }    
  }
  
  Npaper <- length(pap_list)
  
  tf_label <- paper_time$wos %in% pap_list
  paper_interest <- paper_time[valid_ix,][tf_label[valid_ix],]
  W_hat_interest <- W_hat.normed[,tf_label[valid_ix]]
  
  av_tw <- rowSums(W_hat_interest) / dim(W_hat_interest)[2]
  av_tw <- av_tw / sum(av_tw)
  
  print(av_tw)
  
  ylabel <- "Average Weights"
  yrange <- c(0,0.12)
  
  prefix <- 'divid'
  if (normal == 'divid'){
    av_tw <- (av_tw - nm_vector) / nm_vector
    ylabel <- "Normalized Average Weights"
    yrange <- c(-1,1)
  }
  else{
    prefix <- 'minus'
    av_tw <- (av_tw - nm_vector)
    ylabel <- "Centralized Average Weights"
    yrange <- c(-0.05,0.1)
  }
  if(no_ylabel){
    ylabel <- ""
  }
  
  sort_topics <- sort(topics,decreasing = F,index.return = T)$ix
  
  if (!horiz){
    if (topic_names){
      par(mar = c(6.8, 4.1, 3.1, 0.1))
    }
    else{
      par(mar = c(6.8, 4.1, 3.1, 0.1))
    }
    plt <- barplot(av_tw[sort_topics],ylab = ylabel,ylim = yrange,main = name_title,las = 1,
                   xaxt = "n",cex.main = 2,cex.axis = 1.4)
    rot_angle <- 90
    
    if (topic_names){
      par(mar = c(6.8, 4.1, 3.1, 0.1))
      text(plt, par("usr")[3], labels = topics[sort_topics], srt = rot_angle, adj = c(1,0.5), xpd = TRUE, cex=1.3)  
    }
    else{
      par(mar = c(3, 4.1, 3.1, 0.1))
    }
  }
  else{
    # barplot(av_tw,names.arg = topics,xlab = ylabel,xlim = yrange,main = name_title,las = 1,
    # cex.axis = 1.4,horiz = horiz)
    if (topic_names){
      par(mar = c(2.1, 7.5, 3.1, 0.4))
    }
    else{
      par(mar = c(2.1, 7.5, 3.1, 0.4))
    }
    plt <- barplot(av_tw[sort_topics[length(sort_topics):1]],xlab = ylabel,xlim = yrange,main = name_title,las = 1,yaxt = "n",
                   cex.axis = 1.4,horiz = horiz,cex.main = 1.6)
    print(plt)
    
    y_pos <- plt
    y_pos <- y_pos[length(y_pos):1] 
    if (topic_names){
      text(par("usr")[1], y_pos, labels = topics[sort_topics], srt = 0, pos = 2, xpd = TRUE, cex=1.5)   
    }
  }
  
  fname_prefix <- gsub(" ", "",gsub(", ", "", gsub("\\.", "", name_title)))
  if (save){
    if (horiz){
      dev.copy2pdf(file=paste(fname_prefix,"_all_",prefix,".pdf",sep = ''),width = 4.77, height = 5.19)
    }
    else{
      dev.copy2pdf(file=paste(fname_prefix,"_all_",prefix,".pdf",sep = ''),width = 4.77, height = 5.19)
    }
    ## best size: 4.77 * 5.19
  }
  
}

ma <- function(x){
  m <- 1
  # the window width is 3 in the interior
  n <- length(x)
  res <- rep(0,n)
  for (i in 1:n){
    if (i <= m){
      window <- 1:(i + m)
      weight <- c(0.66,0.34)
    }
    else if (i > n - m){
      window <- (i - m):n
      weight <- c(0.34,0.66)
    }
    else {
      window <- (i - m):(i + m)
      weight <- c(0.25,0.5,0.25)
    }
    res[i] <- sum(x[window] * weight)
  }
  return(res)
}


get_lwd <- function(x,threshold){
  delta <- max(x) - min(x)
  if(delta > threshold){
    return(2)
  }
  else{
    return(0)
  }
}

get_pch <- function(x,threshold){
  delta <- max(x) - min(x)
  if(delta > threshold){
    return(1)
  }
  else{
    return(0)
  }
}

get_lty <- function(x,threshold){
  delta <- max(x) - min(x)
  if(delta > threshold){
    return(1)
  }
  else{
    return(2)
  }
}


draw_tw_path <- function(df_tw,name_title = 'Change of Topic Interest',ylabel = 'Topic Interest'){
  col.list <- c('red','green','blue','black','cyan2','goldenrod','blueviolet',
                'coral3','azure4','darkgoldenrod','purple','hotpink1','orange',
                'firebrick3','gold2')
  
  time_window <- df_tw$year
  ntopic <- dim(df_tw)[2] - 1
  
  par(mfrow = c(1,1))
  par(mar = c(3,3,2,7))
  par(xpd = F)
  
  tmp <- df_tw[,-which(names(df_tw) %in% c('year'))]
  
  vari_curve <- function(x){
    return(max(x) - min(x))
  }
  vari <- unlist(apply(tmp,MARGIN = 2,FUN = vari_curve))
  n_dif <- 5
  l.wd <- 2.5
  threshold <- sort(vari,decreasing = TRUE)[n_dif]
  dif.list <- which(vari >= threshold)
  lwd.list <- rep(1,11)
  lwd.list[dif.list] <- 2.5
  lty.list <- rep(2,11)
  lty.list[dif.list] <- 1
  # lwd.list <- apply(tmp,MARGIN = 2,FUN = get_lwd,threshold = threshold)
  # lty.list <- apply(tmp,MARGIN = 2,FUN = get_lty,threshold = threshold)
  pch.list <- apply(tmp,MARGIN = 2,FUN = get_lty,threshold = threshold)
  
  plot(time_window,(tmp[,1]),col = col.list[1],type = 'o',
       lwd = lwd.list[1],lty = lty.list[1],pch = 1,
       ylim = c(min(tmp),max(tmp)),xlab = '',ylab = ylabel,cex.axis = 1.4)
  title(main = name_title,cex.main = 2)
  for (i in 2:11){
    points(time_window,(tmp[,i]),col = col.list[i],pch = i,
           type = 'o',lwd = lwd.list[i],lty = lty.list[i])
  }
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
       lwd = par("lwd"), equilogs = TRUE)
  par(xpd = T)
  legend(legend = topics,col = col.list,x = c(2016,2021),y = c(min(tmp),max(tmp)),
         pch = 1:11,lwd = lwd.list,lty = lty.list,cex = 0.8)
}










get_tw_group_author_based <- function(group,id,topics, no_title = FALSE, name_title = '', yrange = c(),file_prefix = "",no_ylabel = T, normal = 'divid', nm_vector = 0,save = FALSE){
  if (no_title){
    name_title <- ""
  }
  else{
    name_title <- paste(name_title,"Community_",as.character(id),sep = '')
  }
  aut_id <- group[,1][group[,2] == id]
  
  W_aut <- t(author_av_tw[author_av_tw$author %in% aut_id,colnames(author_av_tw) %in% topics])
  av_tw <- rowSums(W_aut) / dim(W_aut)[2]
  #### topic weights of some authors may be 0, so re-normalize it
  av_tw <- av_tw / sum(av_tw)
  
  ylabel <- "Average Weights"
  
  prefix <- 'divid'
  if (normal == 'divid'){
    av_tw <- (av_tw - nm_vector) / nm_vector
    m <- max(abs(av_tw))
    ylabel <- "Normalized Average Weights"
    # yrange <- c(-m * 1.2,m * 1.2)
  }
  else{
    prefix <- 'minus'
    av_tw <- (av_tw - nm_vector)
    ylabel <- "Normalized Average Weights"
    # yrange <- c(-0.12,0.12)
  }
  
  if(no_ylabel){
    ylabel <- ""
  }
  
  sort_topics <- sort(topics,decreasing = F,index.return = T)$ix
  
  if (length(yrange) > 0){
    plt <- barplot(av_tw[sort_topics],ylab = ylabel,ylim = yrange,main = name_title,las = 1,
                   xaxt = "n",cex.main = 2,cex.axis = 1.4)
  }
  plt <- barplot(av_tw[sort_topics],ylab = ylabel,main = name_title,las = 1,
                 xaxt = "n",cex.main = 2,cex.axis = 1.4)
  
  rot_angle <- 90
  par(mar = c(6.8, 4.1, 1.1, 0.5))
  text(plt, par("usr")[3], labels = topics[sort_topics], srt = rot_angle, adj = c(1,0.5), xpd = TRUE, cex=1.3) 
  
  fname_prefix <- gsub(" ", "",gsub(", ", "", gsub("\\.", "", name_title)))
  if (save){
    dev.copy2pdf(file=paste(file_prefix,"_Community",id,"_",prefix,".pdf",sep = ''))
  }
  
}


get_tw_group_author_based_horiz <- function(group,id, name_title = '', yrange = c(),file_prefix = "",
                                            no_ylabel = T, y_factor = 1, normal = 'divid', nm_vector = 0,save = FALSE){
  if (name_title == 'community'){
    name_title <- paste(name_title,"Community_",as.character(id),sep = '')
  }
  aut_id <- group[,1][group[,2] == id]
  
  W_aut <- t(author_av_tw[author_av_tw$author %in% aut_id,colnames(author_av_tw) %in% topics])
  av_tw <- rowSums(W_aut) / dim(W_aut)[2]
  #### topic weights of some authors may be 0, so re-normalize it
  av_tw <- av_tw / sum(av_tw)
  
  
  ylabel <- "Average Weights"
  
  prefix <- 'divid'
  if (normal == 'divid'){
    av_tw <- (av_tw - nm_vector) / nm_vector
    m <- max(abs(av_tw))
    ylabel <- "Normalized Average Weights"
    # yrange <- c(-m * 1.2,m * 1.2)
  }
  else{
    prefix <- 'minus'
    av_tw <- (av_tw - nm_vector)
    ylabel <- "Normalized Average Weights"
    # yrange <- c(-0.12,0.12)
  }
  
  if(no_ylabel){
    ylabel <- ""
  }
  
  sort_topics <- sort(topics,decreasing = F,index.return = T)$ix
  av_tw <- av_tw * y_factor
  
  # if (length(yrange) > 0){
  #   plt <- barplot(av_tw[sort_topics],ylab = ylabel,ylim = yrange,main = name_title,las = 1,
  #                  xaxt = "n",cex.main = 2,cex.axis = 1.4)
  # }
  # plt <- barplot(av_tw[sort_topics],ylab = ylabel,main = name_title,las = 1,
  #                xaxt = "n",cex.main = 2,cex.axis = 1.4)
  # 
  # rot_angle <- 90
  # if (name_title == ""){
  #   par(mar = c(6.8, 3.1, 1.1, 0.5))
  # }
  # else{
  #   par(mar = c(6.8, 4.1, 4.1, 2.1))
  # }
  # text(plt, par("usr")[3], labels = topics[sort_topics], srt = rot_angle, adj = c(1,0.5), xpd = TRUE, cex=1.3) 
  
  plt <- barplot(av_tw[sort_topics[length(sort_topics):1]],xlab = ylabel,xlim = yrange,main = name_title,las = 1,yaxt = "n",
                 cex.axis = 1.4,horiz = TRUE,cex.main = 1.6)
  
  y_pos <- plt
  y_pos <- y_pos[length(y_pos):1] 
  text(par("usr")[1], y_pos, labels = topics[sort_topics], srt = 0, pos = 2, xpd = TRUE, cex=1.5) 
  par(mar = c(2.1, 7.5, 3.1, 0.4))
  
  fname_prefix <- gsub(" ", "",gsub(", ", "", gsub("\\.", "", name_title)))
  if (save){
    dev.copy2pdf(file=paste(file_prefix,"_Community",id,"_",prefix,".pdf",sep = ''))
  }
  
}


##### transform accent letters to english letters, need to be opened with encoding utf-8 to show normally
to_english <- function(s) {
  # 1 character substitutions
  old1 <- "ŚŠŞşśšŽŻŻŻżžźþÅȦàáăặắâãäąåḃČĆÇçćčĐèééėễėêëęȩİİìıíĭîïķľłŁðñņńØÖòóôõöőùůūúûůüÜýÿğĝǧḡeřţt"
  new1 <- "SSSsssZZZZzzzyAAaaaaaaaaaabCCCcccDeeeeeeeeeeIIiiiiiikllLdnnnOOoooooouuuuuuuUyyggggertt"
  s1 <- chartr(old1, new1, s)
  
  # 2 character substitutions
  old2 <- c("œ", "ß", "æ", "ø","a̧","g̃","t́")
  new2 <- c("oe", "ss", "ae", "oe","a","g","t")
  s2 <- s1
  for(i in seq_along(old2)) s2 <- gsub(old2[i], new2[i], s2, fixed = TRUE)
  
  return(s2)
}
