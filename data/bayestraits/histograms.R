
table <- read.csv("outfile.txt", header = TRUE, sep = "\t")
set.seed(42)
p1 <- hist(table$q01,breaks=20)                     # centered at 4
p2 <- hist(table$q10,breaks=20)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,15), ylim=c(0,200))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,15), ylim=c(0,200), add=T)  # second
