#'
#' Do multi-criteria integration of 2 omics datasets. 
#' See paper: Couto Alves, Dysregulation of Complement System and CD4+ T Cell Activation Pathways Implicated in Allergic Response, plosone 2013.
#' 
#' In the paper, we integrated pathway enrichment analyses conducted on both *gene expression*
#' and *GWAS* data. We also provided an hypothesis test for the joint rejection of both null hypotheses.
#' The joint p-value is defined as:
#' $Pg = 1 - \prod_k{1-p_k}$ 
#' 
#' Arguments:
#' p - is a data.frame with 2 columns containing the p-values of the enrichment conducted on GWAS and gene expression
#' s - is a data.frame with 2 columns containing the odds-ration of the enrichment conducted on GWAS and gene expression
#' Notes:
#' * rownames should contain the pathway names.
#' * The labels of the x and y axis are the column names of p data.frame
#' 
#' 
plotParetoIntegration <- function(p,s,alpha=0.05,doPlot=T)
{

  #Truncate p-values < 1E-300 to avoid -inf
  p = apply(X = p,MARGIN = 2,FUN = function(x) ifelse(x<1E-300,1E-300,x))
  
  #Obtain complete records
  input = na.omit( merge(p,s,by="row.names") )

  #get log10 pvalues
  input = cbind(input,logx = -log10(input[,2]),logy = -log10(input[,3]))

  #Obtain p-value of the joint null hypothesis
  input$Pg = 1 - apply(1 - input[,c(2,3)],1,prod);

  #Obtain pareto front
  library(mco)
  pset = as.data.frame( paretoFilter(as.matrix(-input[,c("logx","logy")])) )
  
  #obtain the list of prioritized pathways
  pinput = input[ as.numeric(rownames(pset)),]
  
  if (doPlot == F){
    return(pinput)
  }
    

#---------------------
#Plot colour coded and size coded
#---------------------
#Get colour for each vertex based on p-value
m        = min(5,length(input$Pg))  #Set number of colors 
hc5      = rev(heat.colors(m))
Pg.max   = max(input$Pg)
imcolors = m - round( m * (input$Pg / Pg.max) ) + 1 
vcolors  = hc5[imcolors]

#Generate scale of colours and pvalues for the colorbar in the legend
Pleg = NULL #P value ticks for the colorbar
for (i in 1:m) {
  quantum = mean(input$Pg[imcolors == i])
  Pleg    = c(Pleg,quantum)
}
Cleg = rev(hc5) #Colour of the legends

radstat = rowMeans(input[,c(4,5)])
maxmag  = max(input$logx, input$logy)
maxrad  = max(radstat)
delta   = 0.05*maxmag
radius  = delta * radstat/maxrad

par(mar = c(5,5,3,6))   #space around c(bottom, left, top, right)

plot(input$logx,input$logy,
     xlim = c(min(input$logx) - delta,    delta*5 + max(input$logx)),
     ylim = c(min(input$logy) - delta,    delta*2 + max(input$logy)),
     col  = "white",
     xlab = colnames(p)[1], 
     ylab = colnames(p)[2])

symbols(input$logx,
        input$logy, 
        circles = radius, 
        inches  = F, 
        ann     = F, 
        bg      = vcolors, 
        fg      = vcolors, 
        add     = T )



#---------------------
#Plot Pareto Front
#---------------------
pinput = pinput[order(pinput$logx),]
lines(pinput$logx,pinput$logy,col="dimgrey",lwd = 4, lty = 3)
abline(a = log10(alpha),b = 0,lty = 2)

#Add pathway labels
for (i in 1:nrow(pinput)) {
  text(x = pinput$logx[i],
       y = pinput$logy[i], 
       labels = pinput$Row.names[i], adj = c(0,0),
       pos = NULL, offset = 0.5, vfont = NULL,
       cex = 0.8, col = NULL, font = NULL)
}

#---------------------
#Add colorbar
#---------------------
library(plotrix)
pcoord = par("usr")
color.legend(xl = pcoord[2] + (pcoord[2])/100,
             yb = pcoord[3],
             xr = pcoord[2] + (pcoord[2])/20,
             yt = pcoord[4]/3,
             legend   = sprintf("%1.3f",sort(Pleg,decreasing = T)),
             rect.col = rev(Cleg),
             gradient = "y",
             align    = "rb", 
             cex      = .8)

return(pinput)
}


#Run example
if(F){
  N = 100
  p = data.frame(p1 = runif(N),p2 = runif(N) )
  s = data.frame(s1 = exp(runif(N)),s2 = exp(runif(N)))
  rownames(p)=paste("Pathway",1:N,sep = "")
  rownames(s)=paste("Pathway",1:N,sep = "")
  colnames(p)=c(
    "GWAS enrichment\n-log10(p)",
    "-log10(p)\nGene expression enrichment"
  )
  colnames(s)=c(
    "GWAS enrichment\nOR",
    "Gene expression enrichment\nOR"
  )
  plotParetoIntegration(p,s,alpha=0.05)
}

