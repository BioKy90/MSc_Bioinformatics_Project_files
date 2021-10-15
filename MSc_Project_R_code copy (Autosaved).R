library("ggplot2")
library("pylr")
library("scales")

# Create distribution plot of CDR3 sequence length
TRA <- c(13, 12, 16, 12, 14, 15, 15, 14, 15, 15, 14, 12, 12, 12, 12, 13, 11, 13, 16, 16, 15, 15, 15, 11, 11, 11, 13, 13, 16, 13, 13, 13, 15, 15, 11, 11, 11, 14, 13, 13, 15, 15, 15, 15, 13, 12, 11, 12, 15, 13, 16, 16, 15, 15, 15, 10, 15, 13, 13, 13, 13, 15, 15, 15, 13, 13, 13, 13, 5, 5, 13, 13, 5, 14)
TRB <- c(16, 15, 13, 13, 16, 16, 16, 14, 14, 14, 12, 14, 14, 16, 14, 14, 14, 14, 13, 13, 13, 11, 16, 16, 18, 15, 16, 15, 13, 13, 12, 12, 12, 13, 16, 13, 13, 14, 14, 16, 16, 16, 14, 13, 13, 14, 16, 14, 14, 11, 14, 13, 13, 12, 12, 12, 12, 15, 15, 16, 15, 13, 16, 15, 14, 15, 15, 15, 14, 15, 14, 15, 15, 15)
ggplot() + geom_histogram(aes(x=TRA, fill="TRA gene"), alpha=0.5,binwidth=1) + geom_histogram(aes(x=TRB, fill="TRB gene"), alpha=0.5,binwidth=1) + labs(title='Distribution of CDR3 sequence length in TRA and TRB genes', y='Frequency', x='CDR3 sequence length') +  theme(legend.title=element_blank())
CDR3_len_9mer <- c()

ks.test(TRA, TRB, alternative = c("two.sided"), exact = NULL, tol=1e-8, simulate.p.value=FALSE, B=2000)
# TRA 


library("ggplot2")
library(tidyverse)
CDR3_13mer <- c(8,8,9,9,9,9,10,10,10,13)
names(CDR3_13mer) <- rep("13mer", length(CDR3_13mer))

CDR3_14mer <- c(10,13,13)
names(CDR3_14mer) <- rep("14mer", length(CDR3_14mer))

CDR3_15mer <- c(9,9,9,9,9,9,9,10,10,11,11,11,11,13)
names(CDR3_15mer) <- rep("15mer", length(CDR3_15mer))

CDR3_16mer <- c(9,9,9)
names(CDR3_16mer) <- rep("16mer", length(CDR3_16mer))

all_CDR3 <- data.frame(
    CDR3_length = c(names(CDR3_13mer), names(CDR3_14mer), names(CDR3_15mer), names(CDR3_16mer)),
    
    Human_TRB_epitope_lengths = c(CDR3_13mer, CDR3_14mer, CDR3_15mer, CDR3_16mer),
    
    Epitope_length = factor(c(CDR3_13mer, CDR3_14mer, CDR3_15mer, CDR3_16mer),
                            levels = c(13:8)) 
)

ggplot(all_CDR3, aes(fill = Epitope_length, y = Human_TRB_epitope_lengths, x = CDR3_length)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw()

# TRB

CDR3_13mer <- c(9,9,9,9,9,9,9,10,10)
names(CDR3_13mer) <- rep("13mer", length(CDR3_13mer))

CDR3_14mer <- c(9,9,9,9,9,9,9,9,10,13,13,13)
names(CDR3_14mer) <- rep("14mer", length(CDR3_14mer))

CDR3_15mer <- c(9,9,9,10,10,10,10,10,13)
names(CDR3_15mer) <- rep("15mer", length(CDR3_15mer))

CDR3_16mer <- c(8,8,9,9,9,9,10,10)
names(CDR3_16mer) <- rep("16mer", length(CDR3_16mer))

all_CDR3 <- data.frame(
    CDR3_length = c(names(CDR3_13mer), names(CDR3_14mer), names(CDR3_15mer), names(CDR3_16mer)),
    
    Human_TRB_epitope_lengths = c(CDR3_13mer, CDR3_14mer, CDR3_15mer, CDR3_16mer),
    
    Epitope_length = factor(c(CDR3_13mer, CDR3_14mer, CDR3_15mer, CDR3_16mer),
                            levels = c(13:8)) 
)

ggplot(all_CDR3, aes(fill = Epitope_length, y = Human_TRB_epitope_lengths, x = CDR3_length)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw()


