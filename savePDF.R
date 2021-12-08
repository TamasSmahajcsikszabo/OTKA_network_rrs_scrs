savePDF <- function(plotCall, filename, ...) {
    # plotCall: barmely plot(...) 
    # filename: file neve, helye

    pdf(filename)
    plotCall
    dev.off()
}

x <- rnorm(100)
savePDF(plot(x), "rnorm.pdf")
