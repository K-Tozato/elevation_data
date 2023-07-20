setwd("~/analysis/topo_data/")

##########################
## Make fname_ndata.txt ## ----------
##########################
fnames <- c("./xyz/53-5A.xyz", "./xyz/53-10B.xyz", 
            "./xyz/54-5A.xyz", "./xyz/54-10B.xyz")

nn <- c()
for(i in 1:length(fnames)){
  da <- read.table(fnames[i], header = F, sep=",")
  nn[i] <- nrow(da)
  if(i==1){dat <- da[,2:4]}
  if(i>=2){dat <- rbind(dat,da[,2:4])}
}

range(dat[,1])
range(dat[,2])

write.table(data.frame(c("",fnames),c(length(fnames),nn)), "./input/fname_ndata.txt", row.names = F, col.names = F)

##---------------------------------


#################
## Confirm Map ## ----------
#################
da <- as.matrix(read.table("./input/area.txt"))
nx <- (da[1,2]-da[1,1])/da[3,1] + 1
ny <- (da[2,2]-da[2,1])/da[3,1] + 1

coord <- as.matrix(read.table("./input/coordinate.txt"))
xx <- coord[1:nx,1]
yy <- coord[seq(1,nx*ny,nx),2]
zz <- matrix(coord[,3], nx,ny)

rgl::persp3d(xx,yy,zz, col="gray30")
rgl::aspect3d(max(xx)-min(xx),max(yy)-min(yy),max(zz)-min(zz)*3)
##-------------------------------