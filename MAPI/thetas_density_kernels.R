MA_thetas<-read.table("thetas.dat")
t_density<-density(MA_thetas$V1)
p_density<-density(MA_thetas$V2)
par(mfrow=c(2,2))
hist(MA_thetas$V1,breaks=25)
hist(MA_thetas$V2,breaks=25)
plot(t_density)
plot(p_density)
dev.print("thetas_density_kernels.pdf",device=pdf)
