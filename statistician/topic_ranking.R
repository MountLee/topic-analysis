##### load data, set time window and range of journals #####
setwd("E:/Dropbox/CMU/ADA/Data/data2019")
load("paper_topic.RData")
load(paste("2019-12.rdata",sep = ''))


Citelist <- read.table(paste("citelist-wos.txt",sep = ''),header = TRUE,sep=',',stringsAsFactors = FALSE)

sum(!Citelist$cited %in% paper$wos)

# the two should be equal
sum(!paper$wos %in% Citelist$citing)
a <- unlist(lapply(paper_ref,length))
sum(a == 0)

id_1 <- Citelist$cited %in% paper_time$wos[valid_ix]
id_2 <- Citelist$citing %in% paper_time$wos[valid_ix]
Citelist_topic <- Citelist[id_1 & id_2,]


cedge_topic <- matrix(rep(0,2 * dim(Citelist_topic)[1]),ncol = 2)
wos_list <- paper_time$wos[valid_ix]
for (i in 1:dim(Citelist_topic)[1]){
  cedge_topic[i,1] <- which(wos_list == Citelist_topic$citing[i])
  cedge_topic[i,2] <- which(wos_list == Citelist_topic$cited[i])
}


X <- matrix(rep(0,11 * dim(Citelist_topic)[1]),ncol = 11)
for (i in 1:dim(X)[1]){
  id_citing <- cedge_topic[i,1]
  id_cited <- cedge_topic[i,2]
  X[i,] <- W_hat.normed[,id_cited] - W_hat.normed[,id_citing]
}
dim(X)

log_like <- function(mu){
  tmp <- (X %*% mu)[,1]
  return(sum(tmp - log(1 + exp(tmp))))
}


grad_hessian_l <- function(mu){
  tmp <- 1 / (1 + exp(X %*% mu)[,1])
  X.scale <- sweep(X, MARGIN=1, tmp, `*`)
  grad <- colSums(X.scale)
  
  tmp2 <- exp(X %*% mu)[,1] / (1 + exp(X %*% mu)[,1])
  X.scale2 <- sweep(X, MARGIN=1, tmp2, `*`)
  hessian <- -t(X.scale) %*% X.scale2
  
  return (list(grad = grad, hessian = hessian))
}



#### Newton's method
step_size <- 0.2
max_iter <- 1000
ths <- 1e-12
objective <- rep(0,max_iter)
mu <- rep(0,11)
objective[1] <- log_like(mu)
for (i in 2:max_iter){
  utils <- grad_hessian_l(mu)
  grad <- utils$grad[-1]
  hessian <- utils$hessian[-1,-1]
  mu[-1] <- mu[-1] - step_size * solve(hessian,grad)
  objective[i] <- log_like(mu)
  print(objective[i])
  if ( abs(objective[i] - objective[i - 1])<ths ){
    print("Converge.")
    break
  }
}
mu[sort(mu,decreasing = T,index.return = T)$ix]
topics[sort(mu,decreasing = T,index.return = T)$ix]

a <- mu[sort(mu,decreasing = T,index.return = T)$ix]
a <- a - a[6]
a

plt <- barplot(mu,ylab = "",main = "Export Scores of Topics",las = 1,
               xaxt = "n",cex.main = 2,cex.axis = 1.4)

rot_angle <- 90
par(mar = c(6.8, 4.1, 4.1, 2.1))
text(plt, par("usr")[3], labels = topics, srt = rot_angle, adj = c(1,0.5), xpd = TRUE, cex=1.3)



##### dispersion parameter #####
## note that for paper-to-paper citations, c_ij <= 1 and t_ij <= 2, 
## so we do not need to aggregate Citelist_topic

df_cedge = data.frame(citing = cedge_topic[,1], cited = cedge_topic[,2])
df_cedge = df_cedge[with(df_cedge, order(citing, cited)),]


cedge_topic_extend = rbind(cedge_topic, cedge_topic[,c(2,1)])
tmp = duplicated(cedge_topic_extend)
sum(tmp)

rev_cite = cedge_topic_extend[tmp,]

sum(rev_cite[,1] == rev_cite[,2])

node_set = unique(c(cedge_topic[,1],cedge_topic[,2]))
length(node_set)
## there are 6 self-cites and 212 / 2 = 106 pairs with reverse cites
m = dim(cedge_topic)[1] - 6 - 212 / 2
p = length(mu)

index = rep(TRUE,dim(cedge_topic)[1])
# index[which(cedge_topic[,1] == cedge_topic[,2])] = FALSE
dup = which(rev_cite[,1] >= rev_cite[,2])
for (i in 1:length(dup)){
  ix = which((cedge_topic[,1] == rev_cite[dup[i],1]) & (cedge_topic[,2] == rev_cite[dup[i],2]))
  index[ix] = FALSE
}
sum(index == FALSE)

tij_eq_2 = which(rev_cite[,1] < rev_cite[,2])
index_2 = rep(FALSE,sum(index))
cedge_topic_2 = cedge_topic[index,]
for (i in 1:length(tij_eq_2)){
  ix = which((cedge_topic_2[,1] == rev_cite[tij_eq_2[i],1]) & (cedge_topic_2[,2] == rev_cite[tij_eq_2[i],2]))
  index_2[ix] = TRUE
}


X_v = X[index,]
tmp <- (X_v %*% mu)[,1]
pi = exp(tmp) / (1 + exp(tmp))
tij = rep(1,dim(X_v)[1])
tij[index_2] = 2

sum(tij == 1)
sum(tij == 2)

cij = rep(1, length(tij))
phi = 1 / (m - p + 1) * sum((cij - tij * pi)^2 / tij * pi * (1 - pi))
length(tij)



mu[sort(mu,decreasing = T,index.return = T)$ix]
topics[sort(mu,decreasing = T,index.return = T)$ix]

a <- mu[sort(mu,decreasing = T,index.return = T)$ix]
a <- a - a[6]
a

rank = data.frame(topic = topics[sort(mu,decreasing = T,index.return = T)$ix], mu = a)

## check if there is self-cite
check = df_cedge$citing == df_cedge$cited
sum(check)
which(check)

df_cedge[check,]

check_origin = as.integer(row.names(df_cedge[check,]))
Citelist_topic[check_origin,]
wos_list[4148]

id = df_cedge$citing[check]
paper_time$wos[valid_ix][id]
paper_time[valid_ix,][id,]

i = id[3]
paper_time[valid_ix,][i,]

wos_ = "000337869500013"
ref_id = which(paper$wos == wos_)
paper[ref_id,]
ref_list = paper_ref[[ref_id]]
"000337869500013" %in% ref_list
paper_author[[ref_id]]
which(duplicated(paper$wos))



check_tmp = which(Citelist$citing == "000337869500013")
Citelist[check_tmp,]


tmp = duplicated(cedge_topic)
which(tmp)

## check if there is reverse-cite
cedge_topic_extend = rbind(cedge_topic, cedge_topic[,c(2,1)])


tmp = duplicated(cedge_topic_extend)
which(tmp)

rev_cite = cedge_topic_extend[tmp,]

sum(rev_cite[,1] == rev_cite[,2])

##
pair = rev_cite[1,]
pair_wos = c(wos_list[pair[1]],wos_list[pair[2]])

pair_ref_id = c(which(paper$wos == pair_wos[1]),which(paper$wos == pair_wos[2]))
paper[pair_ref_id,]
ref_list1 = paper_ref[[pair_ref_id[1]]]
ref_list2 = paper_ref[[pair_ref_id[2]]]

pair_wos[1] %in% ref_list2
pair_wos[2] %in% ref_list1

paper_author[[pair_ref_id[1]]]
paper_author[[pair_ref_id[2]]]
##
check = rep(0,dim(rev_cite)[1])
for (i in 1:dim(rev_cite)[1]){
  pair = rev_cite[i,]
  pair_wos = c(wos_list[pair[1]],wos_list[pair[2]])
  
  pair_ref_id = c(which(paper$wos == pair_wos[1]),which(paper$wos == pair_wos[2]))
  paper[pair_ref_id,]
  ref_list1 = paper_ref[[pair_ref_id[1]]]
  ref_list2 = paper_ref[[pair_ref_id[2]]]
  
  if ((pair_wos[1] %in% ref_list2) & (pair_wos[2] %in% ref_list1)){
    check[i] = 1
  }
}
sum(check)

##






### old stuff
check = df_cedge$citing > df_cedge$cited
sum(check)

id = df_cedge$citing[check]
id_rev = df_cedge$cited[check]

i = id[1]
i_rev = id_rev[1]
paper_time[valid_ix,][i,]

wos_ = "000309793400008"
ref_id = which(paper$wos == wos_)
paper[ref_id,]
ref_list = paper_ref[[ref_id]]

paper_author[[ref_id]]
which(duplicated(paper$wos))


i_rev = id_rev[1]
paper_time[valid_ix,][i_rev,]


"000293113300013" %in% ref_list
paper_author[[ref_id]]
which(duplicated(paper$wos))

wos_rev = "000293113300013"
ref_id_rev = which(paper$wos == wos_rev)
paper[ref_id_rev,]
ref_list_rev = paper_ref[[ref_id_rev]]

wos_ %in% ref_list_rev
