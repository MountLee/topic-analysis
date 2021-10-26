library("forcats")
library("tidytext")
library("topicmodels")
library("tm")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rARPACK")
library("tictoc")
library("corrplot")
library("grid") 
library("stringr")


############################################################
##### load data, set time window and range of journals #####
############################################################

setwd("E:/Dropbox/CMU/ADA/Data/data2019")
load("2019-12.RData")

paper$year <- as.numeric(paper$year)

min_year <- 1990
max_year <- max(paper$year)

time_window <- c(min_year:max_year)

paper_time <- paper[paper$year %in% time_window,]

stats <- journal$issn[!journal$shortName %in% c('AIHPP','AoP','PTRF')]
paper_time <- paper_time[paper_time$issn %in% stats,]


##### preprocess words in abstracts #####
common_words <- c("for", "well", "also", "can", "may", "the", "and", "with", "some", "via", "using", "vol",
                  "estimation","estimator","estimators","distributions","distribution","function","functions",
                  "data","paper","our","shown","point","idea","method","show","model","models",'modelling',"problem",
                  "used","use","new","propose","proposed","likelihood","sampling","processes","process",
                  "improve","sharp","reserved","strengthen")

other_words <- c("longrang","art","retriev","retail","emb","compani","northern",
                 "maxstabl","countabl","sudden","aftershock","arch","ultim","site","drive","disciplin",
                 "depart","royal","statistician","journal","canadian","review","preprint","springer","roy",
                 "annal","singh","assoc","ann","wiley","canada","societi","inst","soc","statist","lett","ser",
                 "sinica","length","let","pay","sas","masson","right","all","rth","kth","ith","edit","xnn",
                 "denot","exclud","larger","eta","liu","unrel","errorsinvari",
                 "scientifiqu","erron",'for','of','c',"appl",
                 "claim","come","comment","academ","accord","allow","answer",
                 "around","ask","assess","author","away","bad","becom","begin","belong","besid",
                 "build","built","call","come","conclud","conclus","confirm","consider",
                 "contact","copyright","creat","current","date","detail","discuss","done",
                 "dot","due","easi","eas","easier","enabl","end","enough","equal","etc","even",
                 "eventu","everi","exampl","except","extra","fall","fashion","five","four",
                 "futur","gave","get","give","given","good","goal","great","greater",
                 "greaterthanorequalto","greatest","guid","guidanc","hand","help","henc",
                 "hold","idea","iii","iid","impli","instead","interest","introduc","introduct",
                 "isanelementof","keep","key","know","knowledg","known","last","led",
                 "lee","left","lessthanorequalto","like","lie","lin","list","load","made","main","meet",
                 "mention","met","might","much","most","must","name","need","never",
                 "newli","next","nice","non","now","one","out","paper","per","perform",
                 "perhap","permit","play","plot","put","question",
                 "reason","regard","regardless","run","said","say","sci","see","seen",
                 "shall","set","show","sinc","six","stay","still","strict","suggest",
                 "take","taken","ten","therebi","thus","thought","toward","tool","togeth",
                 "twice","under","upon","usa","various","view","want","way","wish",
                 "work","xii","xxn","yet","york","zero","whole","whose","will",
                 "bickel","berger","fan","nth","xin",
                 "sampl","estim","evolv","tie","law","bar","math","person","prove","proven")

process_words <- function(mydoc){
  mydoc <- iconv(mydoc, "utf-8", "ASCII", sub="")
  mydoc <- Corpus(VectorSource(mydoc))
  # Convert the text to lower case
  mydoc <- tm_map(mydoc, content_transformer(tolower))
  # Remove numbers and punctuations
  mydoc <- tm_map(mydoc, removePunctuation)
  mydoc <- tm_map(mydoc, removeNumbers)
  # Remove english common stopwords
  mydoc <- tm_map(mydoc, removeWords, stopwords("english"))

  # Further modifications
  mydoc <- tm_map(mydoc, removeWords, common_words)

  mydoc <- tm_map(mydoc, content_transformer(gsub), pattern = "\\b(measurements|measurement|measuring)\\b", replacement = "measurement1")
  mydoc <- tm_map(mydoc, content_transformer(gsub), pattern = "\\b(testing error|test error)\\b", replacement = "testerror")
  mydoc <- tm_map(mydoc, content_transformer(gsub), pattern = "\\b(markov chain)\\b", replacement = "markovchain")
  mydoc <- tm_map(mydoc, content_transformer(gsub), pattern = "\\b(monte carlo)\\b", replacement = "montecarlo")

  # Text stemming
  mydoc <- tm_map(mydoc, PlainTextDocument)
  mydoc <- tm_map(mydoc, stemDocument)
  mydoc <- tm_map(mydoc, removeWords, other_words)

  mydoc <- strsplit(mydoc[[1]]$content, " ")

  mydoc <- unlist(mydoc, recursive = FALSE)
  mydoc <- paste(mydoc,collapse = " ")

  return(mydoc)
}


##### prepare the list of abstracts and titles #####

Npaper <- dim(paper_time)[1]
abstracts <- vector(mode = 'list',length = Npaper)
abs_len <- rep(0,Npaper)

tic('get_abstracts')
for (i in 1:Npaper){
  abstracts[[i]] <- process_words(paper_time$abstract[i])
  tmp <- strsplit(abstracts[[i]],' ')[[1]]
  abs_len[i] <- length(tmp[tmp != ''])
  if (i %% 10000 == 0){
    print(i)
  }
}
toc()

summary(abs_len)
sum(abs_len == 0)


###########################################
##### remove some words and documents #####
###########################################


topn <- 56500
valid_ix <- sort(abs_len,decreasing = T,index.return = T)$ix[1:topn]

abstracts_valid <- abstracts[valid_ix]

title_valid <- lapply(as.character(paper_time$title[valid_ix]),FUN = tolower)


#### make the corpus matrix
sum(abs_len[valid_ix] == 0)
max(valid_ix)

aa <- unlist(abstracts_valid)
tdm <- aa %>%
  VectorSource() %>%
  Corpus() %>%
  TermDocumentMatrix()

terms_stats <- tdm[[6]]$Terms
a <- table(tdm[[1]])

threshold_word <- 99
sum(a > threshold_word)
term_index_valid <- sort(unique(tdm[[1]]),decreasing = F)[a > threshold_word]
terms_valid <- tdm[[6]]$Terms[term_index_valid]

tdm_ws <- tdm[term_index_valid,]
terms_order <- order(terms_valid)
terms_valid <- terms_valid[terms_order]
tdm_ws <- tdm_ws[terms_order,]

m <- as.matrix(tdm_ws)
dim(m)

max(m)

###### check: it should be 0
colsum <- colSums(m)
sum(colsum == 0)
#######


###################
##### T-SCORE #####
###################


source(paste("T-SCORE_functions.R",sep = ''))


####### set random seed
rm(".Random.seed")
RNGversion(vstr = "3.6.0")
set.seed(seed = 10403,kind = "Mersenne-Twister")
.Random.seed[1:10]
  
tic("tscore")
k <- 11
result_Tscore <- norm_score(K = k, K0 = ceiling(1.5 * k), m = 10 * k, D = m, Mquantile=1, scatterplot=F, kmeans_start = 2000)
thet <- result_Tscore$theta
R <- result_Tscore$R
V <- result_Tscore$V
Pi <- result_Tscore$Pi
A_hat <- result_Tscore$A_hat
toc()

#########################################
##### anchor words by normalizing A #####
#########################################

r <- rowSums(A_hat)
sum(r == 0)
A.row_normed <- diag(1 / r) %*% A_hat


k <- 11
n0 <- 1
n <- 30
word_topic <- c()
ind_topic <- c()
for (i in 1:k){
  ind_topic[[i]] <- sort(A.row_normed[,i],decreasing = T,index.return = T)$ix[n0:n]
  word_topic[[i]] <- terms_valid[ind_topic[[i]]]
}

for (i in 1:k){
  print(paste(c('word_topic ',as.character(i)),sep = ' '))
  print(word_topic[[i]])
}

##### set the names for topics
topics <- c("Inference","Exp.Design",
            "Time Series","Mach.Learn.",
            "Bayes","Latent.Var.","Clinic.",
            "Math.Stats.","Regression","Hypo.Test","Bio./Med.")

sum(A.row_normed > 0.9)


anchor_df <- word_topic[[1]]
for (i in 2:k){
  anchor_df <- cbind(anchor_df,word_topic[[i]])
}
anchor_df <- data.frame(anchor_df)
colnames(anchor_df) <- topics
anchor_df <- anchor_df[,order(names(anchor_df))]
write.table(anchor_df, file = paste("anchor_words.csv"),
            append = FALSE,sep = ',',row.names = FALSE,quote = FALSE)



#####################################################################################
################ visualize A: anchorness (row-wise normalization) ###################
#####################################################################################


max_A <- apply(A_hat, MARGIN = 2, max)
min_A <- apply(A_hat, MARGIN = 2, min)

A_hat <- result_Tscore$A_hat
terms_tmp <- c()
topic_tmp <- c()
freq_tmp <- c()
word_topic <- c()
ind_topic <- c()
show_n <- 20
for (i in 1:k){
  ind_topic[[i]] <- sort(A.row_normed[,i],decreasing = T,index.return = T)$ix[1:show_n]
  word_topic[[i]] <- terms_valid[ind_topic[[i]]]
}

k <- 11
for (i in 1:k){
  topic_tmp <- c(topic_tmp, rep(topics[i],show_n))
  terms_tmp <- c(terms_tmp, word_topic[[i]])  
  # freq_tmp <- c(freq_tmp, A_hat[ind_topic[[i]],i])
  freq_tmp <- c(freq_tmp, A.row_normed[ind_topic[[i]],i])
}

tscore_topics <- tibble(term = terms_tmp,topic = topic_tmp, beta = freq_tmp)

tscore_terms <- tscore_topics %>%
  group_by(topic) %>%
  top_n(20, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)


tscore_terms %>%
  # mutate(term = reorder(term, beta)) %>%
  mutate(term = reorder(term,beta)) %>%
  ggplot(aes(term, beta, fill = factor(topic))) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),text = element_text(size=15))+
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free", nrow = 2) +
  coord_flip(ylim = c(0.35, 1))

dev.copy2pdf(file=paste("topic-anchor.pdf",sep = ''),width = 13,height = 9)
###################################################
################# Analysis by W ###################
###################################################
lambda = 0.03
W_hat2 <- solve(a = (t(A_hat) %*% A_hat + lambda * diag(rep(1,k))),b = t(A_hat) %*% m)
W_hat <- W_hat2
x <- c(W_hat)
sum(x<0)/length(x)


##### representative papers #####
W_hat.normed <- W_hat
for (i in 1:dim(W_hat)[2]){
  if (sum(W_hat[,i]) != 0){
    tmp <- W_hat[,i]
    tmp <- tmp * (tmp > 0)
    
    if (sum(tmp) != 0){
      W_hat.normed[,i] <- tmp / sum(tmp)
    }
    else{
      W_hat.normed[,i] <- 0
    }
  }
  else{
    W_hat.normed[,i] <- 0
  }
}


####### correlation between topics
order_topic <- order(topics)
A_hat_z <-  scale(A_hat[,order_topic], center = TRUE, scale = TRUE) / (dim(A_hat[,order_topic])[1] - 1)^0.5
colSums(A_hat_z)

G <- t(A_hat_z) %*% A_hat_z
sigma <- diag(diag(G)^(-0.5))
cor_A <- sigma %*% G %*% sigma

par(mfrow = c(1,2))
col1 <- colorRampPalette(c("blue", "white", "red")) 
corrplot(cor_A, method = "color", col = col1(100))


W_hat_z <- scale(t(W_hat.normed[order_topic,]), center = TRUE, scale = TRUE) / (dim(W_hat.normed)[1] - 1)^0.5
colSums(W_hat_z)

G <- t(W_hat_z) %*% W_hat_z
sigma <- diag(diag(G)^(-0.5))
cor_W <- sigma %*% G %*% sigma


col1 <- colorRampPalette(c("blue", "white", "red")) 
corrplot(cor_W, method = "color", col = col1(100))

###### representative papers
k <- 11
index.doc_topic <- vector(mode = 'list',length = k)
top_m <- 10
for (i in 1:k){
  index.doc_topic[[i]] <- sort(W_hat.normed[i,],decreasing = T,index.return = T)$ix[1:top_m]
}

for (i in 1:k){
  print(i)
  print(topics[i])
  for (j in 1:top_m){
    print(tolower(title_valid[index.doc_topic[[i]][j]]))
  }
}

##### output to tables
reorder_topic <- sort(topics,decreasing = FALSE,index.return = TRUE)$ix

index.doc_topic_custom <- index.doc_topic
length_title <- 400
sink("repre_papers_table.txt")
for (ii in 1:k){
  i <- reorder_topic[ii]
  index_i <- index.doc_topic_custom[[i]]
  show_n <- 3
  cat("\\hline \n")
  for (j in 1:show_n){
    id_j <- index_i[j]
    title_j <- tolower(title_valid[id_j])
    substr(title_j,1,1) <- toupper(substr(title_j,1,1))
    a <- nchar(title_j)
    if (a > length_title){
      title_j <- paste(substr(title_j, 1, length_title),"...")
    }
    weight_j <- round(W_hat.normed[i,id_j],digits = 2)
    if (j == 1){
      cat(paste("\\multirow{3}[2]{*}{",topics[i],"} & ", 
                title_j, " & ", weight_j, "\\\\",sep = ""))
      cat("\n")
    }
    else{
      cat(paste(" & ", title_j, " & ", weight_j, "\\\\",sep = ""))
      cat("\n")
    }
  }
}
cat("\\hline")
sink()

########################################
##### change of topic productivity #####
########################################
##### count papers in each class

W_hat <- W_hat.normed
dim(W_hat)
paper.class <- apply(W_hat, MARGIN = 2, which.max)
table(paper.class)

paper_class_num <- rep(0,k)
for (i in 1:k){
  paper_class_num[i] = sum(paper.class == i)
}
paper_class_num
topics
sum(paper_class_num)

#### time_window1 should be a subset of time_window
time_window1 <- 1990:2015
# time_window1 <- 1975:2015
nyear <- length(time_window1)
maxyear <- max(time_window1)
minyear <- min(time_window1)
ntopic <- length(topics)
paper_time$year <- as.numeric(paper_time$year)

paper_time_valid_year <- paper_time$year
paper_time_valid_year <- paper_time_valid_year[valid_ix]

npaper <- length(paper_time_valid_year)

W_hat <- W_hat.normed
topic_time <- matrix(rep(0, nyear * ntopic), ncol = nyear)
for (y in 1:nyear){
  id.year <- (1:npaper)[paper_time_valid_year == (y + minyear - 1)]
  tmp <- rowSums(W_hat[,id.year])
  tmp <- tmp / sum(tmp)
  topic_time[,y] <- t(tmp)
}

topic_time[,26]
topics

##### draw the plot for all topics
par(mfrow = c(1,1))
par(mar = c(3,3,2,7))
par(xpd = F)

col.list <- c("firebrick4","firebrick3","firebrick2","firebrick1",
              "coral1","royalblue4","royalblue3","steelblue3",
              "turquoise3","skyblue2","skyblue1")
pch.list <- c(1:ntopic)
rank_topic <- sort(topic_time[,dim(topic_time)[2]],decreasing = T,index.return = T)$ix 
lwd.list <- rep(2,ntopic)
lwd.list[c(3,6)] <- 3
plot(time_window1,topic_time[rank_topic[1],],col = col.list[1],
     type = 'o',lwd = lwd.list[1],pch = pch.list[1],
     ylim = c(min(topic_time),max(topic_time)),xlab = '',ylab = '')
for (i in 2:ntopic){
  points(time_window1,topic_time[rank_topic[i],],lwd = lwd.list[i],
         col = col.list[i],type = 'o',pch = pch.list[i])
}
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

abline(h = 1/11,col = 'black',lty = 2,lwd = 2) 
par(xpd = T)


legend(legend = c(topics[rank_topic],"Average"),lwd = c(lwd.list,2),col = c(col.list,"black"),
       x = c(2016,2021.5),y = c(min(topic_time),max(topic_time)),
       lty = c(rep(1,ntopic),2),pch = c(pch.list,NA),cex = 0.8)

dev.copy2pdf(file=paste("topic-dynamic.pdf",sep = ''),width = 8,height = 5.5)

#######################################################
###### topic level, patten of different journals ######
#######################################################

journal_name <- paper_time$issn
for (i in 1:length(journal_name)){
  journal_name[i] <- journal$shortName[journal$issn == paper_time$issn[i]]  
}
paper_time$journal <- journal_name


ma <- function(x, m = 1){
  # the window width is 2m + 1 in the interior
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

topic_multi_journal <- function(topic_id,y_adjust = 0,x_title = 0,y_title = 0){
  
  journal_to_draw <- c("AoS","Bka","JASA","JRSSB","Bcs","JMLR","Sini")
  n_j <- length(journal_to_draw)
  
  W_topic <- W_hat.normed[topic_id,]
  
  df_tw <- data.frame(tw = W_topic,journal = paper_time$journal[valid_ix],
                      year = paper_time$year[valid_ix])
  
  av_tw <- aggregate(. ~ year + journal, data = df_tw, FUN = mean)
  
  time_window1 <- as.numeric(as.character(unique(av_tw$year)))
  
  ylabel <- ""
  
  # nm_vector <- nm_vector[nm_vector$year %in% av_tw$year,-which(names(nm_vector) %in% c('year'))]
  # av_tw <- av_tw[,-which(names(av_tw) %in% c('year'))]
  
  
  par(mar = c(2,2,1,1))
  par(xpd = F)
  col.list <- c("firebrick3","steelblue3","goldenrod4",
                "darkorange","darkorchid3","forestgreen","navy")
  pch.list <- c(1:n_j)
  # rank_topic <- sort(av_tw[,dim(topic_time)[2]],decreasing = T,index.return = T)$ix 
  all_tw <- av_tw$tw[av_tw$journal %in% journal_to_draw]
  
  time_window2 <- as.numeric(as.character(unique(av_tw$year[av_tw$journal == journal_to_draw[1]])))
  # plot(time_window2,ma(av_tw$tw[av_tw$journal == journal_to_draw[1]]),col = col.list[1],
  #      type = 'o',lwd = 2,pch = pch.list[1],main = topics[topic_id],
  #      ylim = c(min(all_tw),max(all_tw)),xlab = '',ylab = '')
  
  y_top = max(all_tw) + y_adjust
  plot(time_window2,ma(av_tw$tw[av_tw$journal == journal_to_draw[1]]),col = col.list[1],
       type = 'o',lwd = 2,pch = pch.list[1],main = "",
       ylim = c(min(all_tw),y_top),xlab = '',ylab = '')
  for (i in 2:n_j){
    time_window2 <- as.numeric(as.character(unique(av_tw$year[av_tw$journal == journal_to_draw[i]])))
    points(time_window2,ma(av_tw$tw[av_tw$journal == journal_to_draw[i]]),lwd = 2,
           col = col.list[i],type = 'o',pch = pch.list[i])
  }
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
       lwd = par("lwd"), equilogs = TRUE)
  # abline(h = 1/11,col = 'grey',lty = 2,lwd = 2)
  abline(h = 1/11,col = 'black',lty = 2,lwd = 2)
  # legend("topright",
  #        legend = c(topics[topic_id]), 
  #        cex = 2, horiz = F)
  text(2010 + x_title, y_top + y_title,labels = topics[topic_id],cex=2,pos = 1)
  # par(xpd = T)
  # legend(legend = journal_to_draw,lwd = c(rep(2,n_j)),col = c(col.list),
  #        x = c(2016,2021),y = c(min(all_tw),max(all_tw)),
  #        lty = c(rep(1,n_j)),pch = c(pch.list),cex = 0.8)
  # 
  # dev.copy2pdf(file=paste("topic_dynamic_all_9015.pdf",sep = ''))
  # haha <- gsub("/", "_",gsub(", ", "", gsub("\\.", "", topics[topic_id])))
  # dev.copy2pdf(file=paste(haha,"_dynamic_some.pdf",sep = ''))
}


m <- matrix(c(1:12),nrow = 3,ncol = 4,byrow = TRUE)

layout(mat = m)

topic_order <- order(topics,decreasing = FALSE)
# for (i in 1:11){
#   ii <- topic_order[i]
#   topic_multi_journal(ii)
# }

topic_multi_journal(topic_order[1])
topic_multi_journal(topic_order[2],y_adjust = 0.01)
topic_multi_journal(topic_order[3],y_adjust = 0.005)
topic_multi_journal(topic_order[4],x_title = -2)
topic_multi_journal(topic_order[5],x_title = -1.5)
topic_multi_journal(topic_order[6],x_title = -1)
topic_multi_journal(topic_order[7],x_title = -2)
topic_multi_journal(topic_order[8],y_adjust = 0.02,x_title = -2)
topic_multi_journal(topic_order[9],x_title = -1)
topic_multi_journal(topic_order[10],x_title = -2)
topic_multi_journal(topic_order[11],x_title = -3)


col.list <- c("firebrick3","steelblue3","goldenrod4",
              "darkorange","darkorchid3","forestgreen","navy")
pch.list <- c(1:n_j)
journal_to_draw <- c("AoS","Bka","JASA","JRSSB","Bcs","JMLR","Sini")
n_j <- length(journal_to_draw)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottomright",
       legend = c(journal_to_draw,"Average"), 
       lwd = c(rep(2,n_j),2),col = c(col.list,"black"),
       lty = c(rep(1,n_j),2),pch = c(pch.list,NA),cex = 2, horiz = F)


dev.copy2pdf(file=paste("topic-journal.pdf",sep = ''),width = 10,height = 9)

############################################
###### citation matrix between topics ######
############################################

wos_time <- paper_time$wos[valid_ix]

W_hat = W_hat.normed

dim(W_hat)
paper.class <- apply(W_hat, MARGIN = 2, which.max)
table(paper.class)

cite_list <- read.csv(paste('citelist-wos.txt',sep = ''),header = TRUE,sep = ',',stringsAsFactors = FALSE)

id_1 <- cite_list$cited %in% paper_time$wos[valid_ix]
id_2 <- cite_list$citing %in% paper_time$wos[valid_ix]
cite_list <- cite_list[id_1 & id_2,]

k = dim(W_hat)[1]
cm_topic <- matrix(rep(0,k^2),ncol = k)

#### this would take long time

for (line in 1:dim(cite_list)[1]){
  if (cite_list[line,1] %in% wos_time & cite_list[line,2] %in% wos_time){
    citing <- paper.class[wos_time == cite_list[line,1]]
    cited <- paper.class[wos_time == cite_list[line,2]]
    cm_topic[citing,cited] <- cm_topic[citing,cited] + 1
  }
  if ( line %% 10000 ==0){
    print(as.character(line))
  }
}

sum(cm_topic)
cpm_topic <- cm_topic / sum(cm_topic)
sum(cpm_topic)




colnames(cpm_topic) <- topics
write.csv(cpm_topic[reorder_topic,reorder_topic],file = 'new_cite_mat_topic.csv',
          row.names = topics[reorder_topic])
colnames(cm_topic) <- topics
write.csv(cm_topic[reorder_topic,reorder_topic],file = 'new_cite_mat_topic_raw_count.csv',
          row.names = topics[reorder_topic])


##### check distribution of dimensions of W ####

W_hat = W_hat.normed
W_df = data.frame(t(W_hat))
names(W_df) = topics

sapply( names(W_df), function( y ) {
  quantile( x = unlist( W_df[,  y ] ), 
            c(.05, .10, .20, .30, .40, .50, .60, .70, .80, .90, .95),
            na.rm = TRUE )
})


#####################################################
###### weighted citation matrix between topics ######
#####################################################

wos_time <- paper_time$wos[valid_ix]

W_hat = W_hat.normed

dim(W_hat)

cite_list <- read.csv(paste('citelist-wos.txt',sep = ''),header = TRUE,sep = ',',stringsAsFactors = FALSE)

k = dim(W_hat)[1]
cm_topic <- matrix(rep(0,k^2),ncol = k)

#### this would take long time

wos_to_id <- 1:length(valid_ix)
names(wos_to_id) <- paper_time$wos[valid_ix]
paper_time$wos[valid_ix][1:6]

id_1 <- cite_list$cited %in% paper_time$wos[valid_ix]
id_2 <- cite_list$citing %in% paper_time$wos[valid_ix]
cite_list <- cite_list[id_1 & id_2,]

trun <- sapply( 1:11, function( y ) {
  quantile( x = unlist( W_hat[y, ] ), 
            c(.30),
            na.rm = TRUE )
})


W_hat.trun <- W_hat
W_hat.trun[W_hat <= trun %*% t(rep(1,dim(W_hat)[2]))] = 0

tmp = sapply( 1:11, function( y ) {
  quantile( x = unlist( W_hat.trun[y, ] ), 
            c(.29),
            na.rm = TRUE )
})
tmp


for (line in 1:dim(cite_list)[1]){
  if (cite_list[line,1] %in% wos_time & cite_list[line,2] %in% wos_time){
    citing <- W_hat.trun[, wos_to_id[cite_list[line,1]]]
    cited <- W_hat.trun[, wos_to_id[cite_list[line,2]]]
    cp <- citing %*% t(cited)
    cm_topic <- cm_topic + cp
  }
  if ( line %% 10000 ==0){
    print(as.character(line))
  }
}

sum(cm_topic)
cpm_topic <- cm_topic / sum(cm_topic)
sum(cpm_topic)




#### citing and cited barplot
citing_prop <- colSums(cpm_topic)
cited_prop <- rowSums(cpm_topic)

par(mfrow = c(1,1))

barplot(citing_prop,names.arg = topics,
        ylab = 'Citing Proportion',las = 2)

barplot(cited_prop,names.arg = topics,
        ylab = 'Cited Proportion',las = 2)

#####################################################################
############ output the cross-citation matrix of topics #############
#####################################################################

colnames(cpm_topic) <- topics
write.csv(cpm_topic[reorder_topic,reorder_topic],file = 'w_cite_mat_topic.csv',
          row.names = topics[reorder_topic])
colnames(cm_topic) <- topics
write.csv(cm_topic[reorder_topic,reorder_topic],file = 'w_cite_mat_topic_raw_count.csv',
          row.names = topics[reorder_topic])

adj <- cm_topic
adj <- adj/matrix(rowSums(adj), nrow=11,ncol=11)
colnames(adj) <- topics
write.csv(adj[reorder_topic,reorder_topic],file = 'w_P_cite_mat_topic.csv',
          row.names = topics[reorder_topic])


library(igraph)
adj <- cm_topic
vLabel <- topics

xadj <- adj/matrix(rowSums(adj), nrow=11,ncol=11)
xadj[xadj < 0.11] = 0
g <- graph.adjacency(xadj, mode="directed", 
                     diag=F,weighted=TRUE)


set.seed(0)
plot.igraph(g,layout=layout.fruchterman.reingold(g), 
            edge.width=E(g)$weight*40, 
            edge.arrow.size=1.5,
            edge.curved=TRUE,
            vertex.size=200*(colSums(cpm_topic)-diag(cpm_topic)), 
            vertex.label=vLabel,
            vertex.label.color="red",
            #    edge.label=round(E(g)$weight, 1),
            edge.color="blue"
)

dev.copy2pdf(file=paste("topic-citation.pdf",sep = ''))


P_raw = read.table(paste("P_cite_mat_topic.csv",sep = ''),header = TRUE,sep=',',stringsAsFactors = FALSE)
P_w = read.table(paste("w_P_cite_mat_topic.csv",sep = ''),header = TRUE,sep=',',stringsAsFactors = FALSE)

P_raw[P_raw < 0.09] = 0
P_w[P_w < 0.11] = 0




###################################
##### output as a .RData file #####
###################################

save("A_hat","journal","m","paper","paper_time","result_Tscore","Pi","k","paper_time_valid_year",
     "R","V","W_hat2","author_name","abs_len","abstracts","abstracts_valid","author_pap_list",
     "valid_ix","journal_name","other_words","common_words","paper_author","paper_ref",
     "time_window","time_window1","title_valid","terms_valid","topics","W_hat.normed",
     "cm_topic","cpm_topic","topic_time","paper_class_num","term_index_valid", 
     file = "paper_topic.RData")
