samples <- c("South", "Thistilfj")
likelihoods <- 1:2

for (l in samples) {
  for (lik in likelihoods) {
    pdf(paste0(l, "-gl", lik, "-fit.pdf"))
    y_lim <- c(-9,3)
    par(mar=c(5,5,4,2))

    params <- 6:16
    if (l == "Thistilf") {
      params <- 5:14
    }
    for (r in params) {
      d <- read.table(paste0(l, "-expected-nsfs-", r, ".txt"))
        if (r == params[1]) {
          n <- length(d) + 1
          plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), d, type="l", ylim=y_lim, xlab="Logit mutant frequency", ylab="Logit normalised SFS", cex.lab=1.5, cex.axis=1.2, las=1) 
        } else {
          par(new=TRUE)
          plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), d, type="l", ylim=y_lim, xaxt="n", yaxt="n", xlab="", ylab="") 
        }
      }

      d <- read.table(paste0(l,"-gl", lik, "-nsfs.txt"))
      d <- log(d) - log(1 - d)
      par(new=TRUE)
      plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), d, ylim=y_lim, xaxt="n", yaxt="n", xlab="", ylab="", col="red")
      legend("topright", legend=c(paste("c =", paste(params, collapse = ", "), collapse = " "), "25kb fragments"), lty=c(1,NA), pch=c(NA,1), col=c("black", "red"), cex=1.2, lwd=c(2,2))
      dev.off()

      pdf(paste0(l, "-gl", lik, "-residuals.pdf"))
      y_lim = c(-1, 3)
      par(mar=c(5,5,4,2))

      param <- 8
      if (l == "Thistilf") {
        param <- 6
      }

      baseline <- read.table(paste0(l, "-expected-nsfs-", r, ".txt"))
      plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), baseline - baseline, type="l", ylim=y_lim, xlab="Logit mutant frequency", ylab="Residual", cex.lab=1.5, cex.axis=1.2, las=1) 

      d <- read.table(paste0(l,"-gl", lik, "-nsfs.txt"))
      d <- log(d) - log(1 - d)
      par(new=TRUE)
      plot(log(1:(n-1)/n) - log(1 - 1:(n-1)/n), d - baseline, ylim=y_lim, xaxt="n", yaxt="n", xlab="", ylab="", col="red")
      
      legend("topleft", legend="25kb fragments", lty=NA, pch=1, col="red", cex=1.2, lwd=2)
      dev.off()
    }
  }
}