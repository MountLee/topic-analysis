library("ggplot2")
library("tictoc")
library('R.matlab')

setwd("E:/Dropbox/CMU/ADA/Data/data2019")

source(paste('mixSCORE.R',sep = ''))
source(paste('topic_interpret_utils.R',sep = ''))
load(paste("paper_topic.rdata",sep = ''))


##### compute the average topic weight of all papers #####
av_tw_all <- rowSums(W_hat.normed) / dim(W_hat.normed)[2]
av_tw_all <- av_tw_all / sum(av_tw_all)

############## static/global analysis ##############
name_title <- 'Average Weights of All Papers'
barplot(av_tw_all,names.arg = topics,ylab = "Average Weights",ylim = c(0,0.12),main = name_title,las = 2)  

##### some examples #####

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

par(mfrow = c(3,4))

author_name[search_name(c('Berger','J'),author_name)]
id <- search_name(c('Berger','J'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE,topic_names = TRUE)
# get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
#               normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name('Bickel, Peter',author_name)]
id <- search_name('Bickel, Peter',author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Carroll','R'),author_name)]
id <- search_name(c('Carroll','R'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Donoho'),author_name)]
id <- search_name(c('Donoho'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

id <- search_name(c('Fan, Jianqing'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id],
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE,topic_names = TRUE)
# get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
#               normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Hall, Peter'),author_name)]
id <- search_name(c('Hall, Peter'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Jordan, Michael'),author_name)]
id <- search_name(c('Jordan, Michael'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

id <- search_name(c('Lin, Xihong'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Liu', 'Jun'),author_name)]
id <- search_name(c('Liu', 'Jun'),author_name)[2]
##### check which is Jun Liu at Havard
# id <- search_name(c('Liu', 'Jun'),author_name)[2]
# author_pap_list[id[2]]
# author_pap_list[id[3]]
# b <- author_pap_list[[id[3]]]
# paper$title[b]
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id],
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE,topic_names = TRUE)
# get_tw_author(id, topics, use_name = FALSE, name_title = "Liu, Jun", 
#               normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Kathryn',"R"),author_name)]
id <- search_name(c('Kathryn',"R"),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Wasserman','L'),author_name)]
id <- search_name(c('Wasserman','L'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = author_name[id], 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)

author_name[search_name(c('Cun',"Zhang"),author_name)]
id <- search_name(c('Cun-hui'),author_name)
get_tw_author(id, topics, use_name = FALSE, name_title = "Zhang, Cun-Hui", 
              normal = 'minus', nm_vector = av_tw_all,save = FALSE,no_ylabel = TRUE)


dev.copy2pdf(file=paste("author-topic.pdf",sep = ''),width = 15,height = 9)



############ high-degree authors ###############
authorN <- length(author_name)
AllWOS <- as.character(paper$wos)
indexAuthorPaper = c()
for(i in 1:authorN){
  if(i%%1000==0) print(i)
  idx = (author_pap_list[[i]])
  thisIndex = cbind(rep(i,length(idx)),idx)
  indexAuthorPaper = rbind(indexAuthorPaper,thisIndex)
}
indexAuthorPaper = as.data.frame(indexAuthorPaper)
names(indexAuthorPaper) = c("idxAuthor","idxPaper")

Papercite <- read.table(paste("PaperCite.txt",sep = ''),header = TRUE,sep=',',stringsAsFactors = FALSE)
sumpap <- rowSums(Papercite[,-1])
citeNum <- sumpap
paper_id <- 1:dim(Papercite)[1]
paper_cite_num <- data.frame(idxPaper = paper_id,citeNum = citeNum)

author_paper_citeNum <- merge(x = indexAuthorPaper, y = paper_cite_num, by = "idxPaper")

author_citeNum <- aggregate(citeNum ~ idxAuthor, data = author_paper_citeNum, FUN = sum)

ranking <- order(author_citeNum$citeNum,decreasing = T)
a <- author_citeNum[ranking[1:10],]
author_name[a$idxAuthor]


author_select <- c()
id_select <- c()
show_n <- 80

#### citation
citat_deg <- author_citeNum$citeNum
id_deg_sort <- sort(citat_deg,decreasing = T,index.return = T)$ix

# show_n <- 45
author_select <- append(author_select,author_name[id_deg_sort[1:show_n]])
id_select <- append(id_select,id_deg_sort[1:show_n])


show_n <- 40
{cat("To be selected: \n")
cat("min citations: ", sort(citat_deg, decreasing = TRUE)[show_n]," \n")}


##### output selected authors to a table #####
citat_select <- citat_deg[id_select]

sort_alpha <- sort(author_select,decreasing = F,index.return = T)$ix
author_select <- author_select[sort_alpha]
id_select <- id_select[sort_alpha]
citat_select <- citat_select[sort_alpha]

author_df <- data.frame(name = author_select,
                        citation = citat_select)

write.table(author_df,file = "author_select.csv",quote = T,row.names = F,sep = ',')





##### output selected authors with topic weights to a table #####

av_tw_all <- rowSums(W_hat.normed) / dim(W_hat.normed)[2]
av_tw_all <- av_tw_all / sum(av_tw_all)

repre_author_df <- data.frame(name = c(),topic = c(),
                              weight = c(),topic2 = c(),
                              weight2 = c(),topic3 = c(),
                              weight3 = c(),topic4 = c(),
                              weight4 = c(),topic5 = c(),
                              weight5 = c())

for (ii in 1:length(author_select)){
  i <- id_select[ii]
  
  pap_list <- unlist(author_pap_list[i])
  pap_list <- paper$wos[pap_list]
  Npaper <- length(pap_list)
  
  tf_label <- paper_time$wos %in% pap_list
  paper_interest <- paper_time[valid_ix,][tf_label[valid_ix],]
  W_hat_interest <- W_hat.normed[,tf_label[valid_ix]]
  
  av_tw <- rowSums(W_hat_interest) / dim(W_hat_interest)[2]
  av_tw <- av_tw / sum(av_tw)
  
  W_tmp <- av_tw - av_tw_all
  weight_rank_id <- sort(W_tmp,decreasing = TRUE,index.return = TRUE)$ix
  weight_k <- W_tmp[weight_rank_id[1:5]]
  topic_k <- topics[weight_rank_id[1:5]]
  
  use_id <- which(weight_k >= weight_k[1] * 0.4)
  
  # df_this <- data.frame(name = author_select[ii],
  #                       topic = paste(topic_k[topic_k != ""],collapse = ', '))
  # repre_author_df <- rbind(repre_author_df,df_this)
  if(length(use_id) > 10){
    df_this <- data.frame(name = c(author_select[ii],""),
                          topic = c(paste(topic_k[use_id[1:3]],collapse = ', '),paste(topic_k[use_id[-c(1:3)]],collapse = ', ')))
    repre_author_df <- rbind(repre_author_df,df_this)
  }
  else{
    df_this <- data.frame(name = author_select[ii],
                          topic = paste(topic_k[use_id],collapse = ', '))
    repre_author_df <- rbind(repre_author_df,df_this)
  }
}
write.csv(repre_author_df, file = "author_select_topic_04.csv",
          row.names = FALSE, quote = TRUE)

