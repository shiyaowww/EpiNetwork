library(openxlsx)
library(readxl) # NA on greatlakes
library(ggplot2)
library(reshape2) 

# define a function that extracts worksheets as data frames in a list
read_allsheets <- function(filename,tibble=FALSE){
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet=X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# extract matrices from sheets
matrices <- read_allsheets("/Users/shiyaow/Desktop/matrices.xlsx")
matrices

# vertical and hirizontal labels
assortative_mixing <- c(0,0.1,0.25,0.5,0.75,0.9,1)
susceptibility_dispaired <-c(1,1.25,1.5,1.75,2) 

# create melted dfs for plotting 
plot_dfs <- list()
for(i in 1:16){
  matrices[[i]] <- data.matrix(matrices[[i]])
  colnames(matrices[[i]]) <- susceptibility_dispaired
  rownames(matrices[[i]]) <- assortative_mixing
  plot_dfs[[i]] <- melt(matrices[[i]])
}

# create key 
key <- c("0.12lowses_prev_overall", "0.12lowses_prev_ses1", "0.12lowses_prev_ses0", "0.12lowses_PR",
         "0.12lowses_risk_overall", "0.12lowses_risk_ses1", "0.12lowses_risk_ses0", "0.12lowses_RR",
         "0.3lowses_prev_overall", "0.3lowses_prev_ses1", "0.3lowses_prev_ses0", "0.3lowses_PR",
         "0.3lowses_risk_overall", "0.3lowses_risk_ses1", "0.3lowses_risk_ses0", "0.3lowses_RR")

# create contour plots
contour <- list()
for(i in 1:16){
  contour[[i]] <-
    ggplot(plot_dfs[[i]], aes(x=Var1,y=Var2,z=value)) +
    geom_point(aes(color=value)) + 
    stat_contour(na.rm=T) +
    labs(x="assortative mixing scenarios",y="susceptibility disparity levels",fill="height", title=key[i])
}

# save 16 plots for 16 matrices
pdf("/Users/shiyaow/Desktop/16matrices.pdf")
contour
dev.off()
