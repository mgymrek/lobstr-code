trf = read.csv("trf.coords.bed", sep="\t", header=T)

trf$length = trf$end - trf$start
trf = trf[trf$length <= 80,]

png("period_vs_length.png")
plot(trf$period, trf$length, pch = 20, col="lightgray")

periods = sort(unique(trf$period))
medians = c()
for (period in periods) {
    period_median = median(trf[trf$period == period,]$length)
    medians = c(medians, period_median)
}
points(periods, medians, col="red", pch = 20)
dev.off()