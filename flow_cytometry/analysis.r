library(ggplot2)

# preset parameters
# threshold of mKate positive cells
min_mKate <- 10^1.8
# range for mKate
mKate_bin <- 10^seq(1.8, 4.8, 0.2)
# sampling numbers for noise calculation
sample.num <- 100
# threshold for outlier removal
outlier.cut <- 0.05

# read files
path <- "./data/"
d.case <- read.table(paste0(path, "case.txt"), header = T)
d.ctrl <- read.table(paste0(path, "ctrl.txt"), header = T)

d.case$miRNA <- "miRNA"
d.ctrl$miRNA <- "ctrl"

d <- rbind(d.case, d.ctrl)
d <- d[d$EYFP > 0 & d$mKate > min_mKate, ]

# function to calculate noise
CalculateNoise <- function(d){
    bin.num <- length(mKate_bin) - 1
    mKate_mean <- c()
    mKate_sd <- c()
    EYFP_mean <- c()
    EYFP_sd <- c()
    EYFP_CV <- c()
    EYFP_CVsd <- c()
    cell_count <- c()
    
    # get noise
    for(i in 1:bin.num){
        this_data <- d[d$mKate > mKate_bin[i] & d$mKate < mKate_bin[i + 1], ]
        this_data <- this_data[order(this_data$EYFP), ]
        # filter out outliers
        cut_range_low <- outlier.cut * nrow(this_data)
        cut_range_high <- (1 - outlier.cut) * nrow(this_data)
        this_data <- this_data[cut_range_low:cut_range_high, ]
      
        mKate_mean_temp <- rep(NA, times = sample.num)
        EYFP_mean_temp <- rep(NA, times = sample.num)
        EYFP_CV_temp <- rep(NA, times = sample.num)
      
        for (j in 1:sample.num){
            bin_sample_temp <- this_data[sample(1:nrow(this_data), ceiling(nrow(this_data)/2)),]
            mKate_mean_temp[j] <- mean(bin_sample_temp$mKate)
            EYFP_mean_temp[j] <- mean(bin_sample_temp$EYFP)
            EYFP_CV_temp[j] <- sd(bin_sample_temp$EYFP) / mean(bin_sample_temp$EYFP)
        }
      
        mKate_mean[i] <- mean(mKate_mean_temp)
        mKate_sd[i] <- sd(mKate_mean_temp)
        EYFP_mean[i] <- mean(EYFP_mean_temp)
        EYFP_sd[i] <- sd(EYFP_mean_temp)
        EYFP_CV[i] <- mean(EYFP_CV_temp)
        EYFP_CVsd[i] <- sd(EYFP_CV_temp)
        cell_count[i] <- nrow(this_data)
    }
    
    this_data <- data.frame(mKate_mean = mKate_mean, mKate_sd = mKate_sd,
                            EYFP_mean = EYFP_mean, EYFP_sd = EYFP_sd, 
                            EYFP_CV = EYFP_CV, EYFP_CVsd = EYFP_CVsd,
                            cell_count = cell_count)

	return(this_data)
}

# start calculating noise
miRNA_all <- unique(d$miRNA)

d.result <- data.frame()
for(this.f in miRNA_all){
	print(this.f)
	this.d <- d[d$miRNA == this.f, ]
	cal.d <- CalculateNoise(this.d)
	cal.d$miRNA <- this.d$miRNA[1]
	d.result <- rbind(d.result, cal.d)
}

d.result <- d.result[d.result$cell_count >= 100, ]
write.table(d.result, "result/result.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# visualization
pdf("result/mean.pdf", width = 4, height = 3)
p <- ggplot(d.result, aes(x = log10(mKate_mean), y = log10(EYFP_mean),
  ymax = log10(EYFP_mean + EYFP_sd), ymin = log10(EYFP_mean - EYFP_sd)))
p <- p + geom_point(aes(color = miRNA)) + geom_line(aes(color = miRNA))
p <- p + theme_bw()
p <- p + ylim(1, 5)
print(p)
dev.off()

pdf("result/CV.pdf", width = 4, height = 3)
p <- ggplot(d.result, aes(x = log10(EYFP_mean), y = EYFP_CV,
  ymax = (EYFP_CV + EYFP_CVsd), ymin = (EYFP_CV - EYFP_CVsd)))
p <- p + geom_point(aes(color = miRNA)) + geom_line(aes(color = miRNA))
p <- p + ylim(0, 0.7)
p <- p + theme_bw()
print(p)
dev.off()

