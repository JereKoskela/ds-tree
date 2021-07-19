samples <- c("South", "Thistilfj")
likelihoods <- 1:2

hist_tags = c("a","b","c","d")
trace_tags = c("e","f","g","h")

# Length of burnin for histogram.
burnin <- 20000

ind <- 1
for (l in samples) {
  for (lik in likelihoods) {
    data <- read.table(paste0(l, "-gl", lik, "-abc.txt"))
    n <- dim(data)[1]

    pdf(paste0(l, "-gl", lik, "-c-trace.pdf"))
    par(mar=c(5,5,4,2))
    plot(1:n, data[1:n,1],type="l",xlab="Iteration",ylab="c",cex.lab=2.5,cex.axis=2,las=1)
    mtext(trace_tags[ind], side=3, adj=0, outer=FALSE, cex=3)
    dev.off()
  
    pdf(paste0(l, "-gl", lik, "-c-hist.pdf"))
    par(mar=c(5,5,4,2))
    hist(data[burnin:n,1],prob=TRUE,cex.lab=2.5,cex.axis=2,las=1,main="",xlab="c",ylab="")
    mtext(hist_tags[ind], side=3, adj=0, outer=FALSE, cex=3)
    dev.off()
    ind <- ind + 1
  }
}