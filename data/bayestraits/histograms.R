
table <- read.csv("outfile.txt", header = TRUE, sep = "\t")
p1 <- hist(table$q01,breaks=30)
p2 <- hist(table$q10,breaks=20)

# plot the first histogram: rates for scalariform -> simple (purplish color)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,15), ylim=c(0,200), xlab="Instantaneous transition rate" )

# add the second: rates for simple -> scalariform (salmon color)
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,15), ylim=c(0,200), add=T)  

# Two-sample Kolmogorov-Smirnov test
# data:  table$q01 and table$q10
# D^- = 0.923, p-value < 2.2e-16
# alternative hypothesis: the CDF of x lies below that of y
ks.test(table$q01,table$q10,alternative='l')

