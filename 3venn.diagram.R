circle <- function(x, y, r, ...) {
  ang <- seq(0, 2*pi, length = 100)
  xx <- x + r * cos(ang)
  yy <- y + r * sin(ang)
  polygon(xx, yy, ...)
}

venndia <- function(A, B, C, getdata=FALSE, ...){
  cMissing <- missing(C)
  if(cMissing){ C <- c() }
  
  unionAB <- union(A, B)
  unionAC <- union(A, C)
  unionBC <- union(B, C)
  uniqueA <- setdiff(A, unionBC)
  uniqueB <- setdiff(B, unionAC)
  uniqueC <- setdiff(C, unionAB)
  intersAB <- setdiff(intersect(A, B), C)
  intersAC <- setdiff(intersect(A, C), B)
  intersBC <- setdiff(intersect(B, C), A)
  intersABC <- intersect(intersect(A, B), intersect(B, C))
  
  nA <- length(uniqueA)  
  nB <- length(uniqueB)
  nC <- length(uniqueC)
  
  nAB <- length(intersAB)
  nAC <- length(intersAC)
  nBC <- length(intersBC)
  
  nABC <- length(intersABC)	
  
  par(mar=c(2, 2, 0, 0))
  plot(-10, -10, ylim=c(0, 9), xlim=c(0, 9), axes=FALSE, ...)
  circle(x=3, y=6, r=3, col=rgb(1,0,0,.5), border=NA)
  circle(x=6, y=6, r=3, col=rgb(0,.5,.1,.5), border=NA)
  circle(x=4.5, y=3, r=3, col=rgb(0,0,1,.5), border=NA)
  
  text( x=c(1.2, 7.7, 4.5), y=c(7.8, 7.8, 0.8), c("A", "B", "C"), cex=3, col="gray90" )
  
  text(
    x=c(2, 7, 4.5, 4.5, 3, 6, 4.5), 
    y=c(7, 7, 2  , 7  , 4, 4, 5), 
    c(nA, nB, nC, nAB, nAC, nBC, nABC), 
    cex=2
  )
  
  if(getdata){
    list(A=uniqueA, B=uniqueB, C=uniqueC, 
         AB=intersAB , AC=intersAC , BC=intersBC , 
         ABC=intersABC
    )
  }
}

#create matrix venn from unions and intersects
vennMatrix = matrix("", ncol = length(venn), nrow = max(sapply(venn, length)))
colnames(vennMatrix) = names(venn)
for(n in names(venn)) vennMatrix[1:length(venn[[n]]), n] = venn[[n]]

#save vennMatrix
vennMatrix_fn = paste("vennList", dateTIME= format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt", sep="_")
write.table(vennMatrix, 
            file      = vennMatrix_fn,
            col.names = TRUE,
            quote     = FALSE, 
            sep       = '\t')
