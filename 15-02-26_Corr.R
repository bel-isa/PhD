#Additional computation following the random forest imputation

# Correlation matrix ------------------------------------------------------

library(Hmisc)

Metabs_p <- rcorr(Z, type="pearson") #can also be spearman
Metabs_s <- rcorr(Z, type="spearman")

# To file -----------------------------------------------------------------

df.Metabs_p.r = data.frame(Metabs_p$r)
df.Metabs_p.n = data.frame(Metabs_p$n)
df.Metabs_p.p = data.frame(Metabs_p$P)

df.Metabs_s.r = data.frame(Metabs_s$r)
df.Metabs_s.n = data.frame(Metabs_s$n)
df.Metabs_s.p = data.frame(Metabs_s$P)

write.csv(df.Metabs_p.r,'Full_Z_pearson.csv')
write.csv(df.Metabs_p.p,'Full_Z_pearson_sign.csv')

write.csv(df.Metabs_s.r,'Full_Z_spearman.csv')
write.csv(df.Metabs_s.p,'Full_Z_spearman_sign.csv')

# Make some nice plots ----------------------------------------------------

image(h)

image(k)

myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(10,10,2.5,6))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, las= 2, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7) 
  
  # Color Scale
  par(mar = c(5,3.5,3.5,5))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

myImagePlot(h, title=c("Pearson correlation coefficient"))
myImagePlot(k, title=c("Spearman rank coefficient"))