library(corrplot)
library(readr)

d <- read.csv("Results/norm_table.csv",stringsAsFactors=F, na.strings="unknown")
data <- data.matrix(d)
data <- as.matrix(data)
a=prop.table(data,2)
a[is.nan(a)] <- 0
corrplot(a,col.lim = c(0,1),col = colorRampPalette(c("white", "lightskyblue3", "darkblue"))(20))

d <- read.csv("Results/canc_table.csv",stringsAsFactors=F, na.strings="unknown")
data <- data.matrix(d)
data <- as.matrix(data)
a=prop.table(data,2)
a[is.nan(a)] <- 0
corrplot(a,col.lim = c(0,1),col = colorRampPalette(c("white", "coral1", "indianred4"))(20))