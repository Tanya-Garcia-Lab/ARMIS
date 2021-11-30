nsub<-6
n.time.points<-6
times<-c(0, 1, 2, 3, 4, 5)
center<-c(100, 102, 90, 120, 110, 130)
sd<-c(5, 7, 3, 4, 5, 4)
group<-1

resp<-vector("list", nsub)
list.data<-vector("list", nsub)

for(j in 1:n.time.points){
subid<-rep(j, n.time.points)
resp[[j]]<-rnorm(n.time.points, center[j], sd[j])
list.data[[j]]<-cbind(subid, group, times, resp[[j]])
}
pseudo_steer<-do.call("rbind.data.frame", list.data)
names(pseudo_steer)<-c("Steer", "Group", "Time", "citrulline")

library(usethis)
write.csv(pseudo_steer, file="data_raw/pseudo_steer.csv")

usethis::use_data(pseudo_steer, overwrite = TRUE)
save(pseudo_steer, file="data/pseudo_steer.rdata")
