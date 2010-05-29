#################################################
#################################################
## Balance Statistics Output with Plots
## Mark Huberty
## 10 January 2010
## v0.4

## This code provides two functions to generate tables and plots of the balance statistics
## (p-values for the t-test for difference in means, and the ks statistic for difference in distributions)
## for the output of the MatchBalance function in the Matching() package for R. 

## This function combines code from the Balance Statistics functions written by Charlie Gibbons (http://cgibbons.berkeley.edu/Courses/PS236_F08/balanceTable.r)
## and the balance statistics plotting functions written by Rocio Titiunik (http://are.berkeley.edu/~rocio/R/graph.pval.public.R).
## The changes here are to make those two functions cooperate and enhance their error handling and flexibility.

## BETA.
## v0.2 corrected for column mismatch that made means incorrect.
## v0.3 added error handling and optional inclusion of post-matching means
## v0.4 updates the error handling and does some re-formatting to make the code easier to read

## TO DO: Provide options for inclusion of t-statistics and p-values (i.e., pval=TRUE, tstat=TRUE) for the matrix function

#################################################
#################################################

## Note: KS test statistics appear as 'NA' for dichotomous variables. Also,
##   'BM' is before matching and 'AM' is after matching.


# cov / covariates: the list of variable names corresponding to the BalanceMatrix given to GenMatch. This can either be the BalanceMatrix (data type "matrix") or a vector of covariate names (data type "character"). The latter is useful if BalanceMatrix contains categorical vectors with n categories, which R converts to n dummy variables for matching. In that case, the row count of the BalanceMatrix will not match the row count of the MatchBalance results, and thus cannot be used to generate variable names for MatchBalance output.

# bal.out:       the output of MatchBalance

# 'results':     a matrix whose rows are different variables; whose first two columns contain the means for treated and control;
#                and whose remaining columns have the pvalues to be plotted for every variable. In this function, "results" takes the
#                output from the initial balanceMat function.

# 'title':       title of the overall graph

# at1, at2,at3:  scalars which indicates where to locate the three differents groups (mean treatment, mean controls, graph area) in the                 figure area

# xlim1 :        the left limit of the x-axis; right limit is always set to 1

# textsize:      scalar indicating the size of text in the figure

# legend:        logical indicating whether the legend should be included

# legendx:       scalar indicating the x-coordinate of the legend's location

# legendy:       scalar indicating the y-coordinate of the legend's location

# parcex:        scalar setting cex parameter

# means:         logical indicating whether means should be included in the matrix or plot output

# title:         string for title to apply to the plot

# color:         logical for whether the plot should use color or b/w symbols


###############################################
## BEGIN CODE
###############################################


  balanceMat <- function(cov,
                         bal.out,
                         means=TRUE
                         ){

    # Check to make sure that the data is in the right format
    dtype <- class(cov)
    if(dtype!="character"&dtype!="matrix"&dtype!="data.frame")
      {stop("Covariate names must be supplied as a character vector, or matrix or data frame with column names.")}
         
    
    if(dtype=="matrix")
      {
        if(dim(cov)[2]!=length(bal.out$AfterMatching[]))
          {
            stop("Rowcount for the covariate matrix does not equal the covariate count from MatchBalance")
         }
      }
    if(dtype=="data.frame")
      {
        if(dim(cov)[2]!=length(bal.out$AfterMatching[]))
          {
            stop("Rowcount for the covariate data frame does not equal the covariate count from MatchBalance")
          }
      }
    if(dtype=="character"){
      if(length(cov)!=length(bal.out$AfterMatching[]))
        {
          stop("List of covariate names is not the same length as the list of covariates.")
        }
    }

    
  
                                      
    # Calculate the number of covariates
    n <- ifelse(dtype=="matrix"| dtype=="data.frame", dim(cov)[2], length(cov))

    # Determine how the covariate names are provided, and then grab them
    if(dtype=="matrix"|dtype=="data.frame")
      {
        rnames <- dimnames(cov)[[2]]
      }else{
        if(dtype=="character") rnames <- cov
      }
    
    # Construct the matrix of statistics from the MatchBalance data and attach it to the covariate names
    z <- t(sapply(1:n, function(x){
      c(rnames[x],
        round(bal.out$AfterMatching[[x]]$mean.Tr,3),
        round(bal.out$AfterMatching[[x]]$mean.Co,3),
        round(bal.out$BeforeMatching[[x]]$tt$p.value,2),
        round(bal.out$AfterMatching[[x]]$tt$p.value,2),
        round(bal.out$BeforeMatching[[x]]$tt$statistic,2),
        round(bal.out$AfterMatching[[x]]$tt$statistic,2),
        ifelse(is.null(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue) ==
               0,round(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue,2),
               NA),
        ifelse(is.null(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue) ==
               0, round(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue,2),
               NA)
        )
    }
                  )
           )
    
    z <- as.data.frame(z)
    ## print(z)
    
    z[,2:9] <- apply(z[,2:9], 2, function(x){as.numeric(x)})

    ## Determine if means should be included in the matrix
    ## and format appropriately.
    if(means==TRUE){
      mat <- z[,2:9]
      
      ## Apply the correct column names
      names(mat)<- c("Mean Tr.",
                     "Mean Con.",
                     "BM t p-value",
                     "AM t p-value",
                     "BM t stat",
                     "AM t stat",
                     "BM KS p-value",
                     "AM KS p-value")
      ## Apply the correct row names
      dimnames(mat)[[1]] <- z[,1]
      ## print(mat)

      mat
    }else{
      mat <- z[,4:9]
      
      ## Apply the correct column names
      names(mat)<- c("BM t p-value",
                     "AM t p-value",
                     "BM t stat",
                     "AM t stat",
                     "BM KS p-value",
                     "AM KS p-value")
      ## Apply the correct row names
      dimnames(mat)[[1]] <- z[,1]
      mat
    }

    
  } ## end balanceMat()
 

plot.pval <- function(covariates,
                      bal.out,
                      title=NULL,
                      means=TRUE,
                      color=TRUE,
                      legend=TRUE,
                      legendx=0.15,
                      legendy=2.2,
                      textsize=0.9,
                      parcex=0.8,
                      at1=-0.35,
                      at2=-0.15,
                      at3=-0.9,
                      xlim1=-0.85
                      ) {


  ## Take the function above and apply it to the data supplied in the command
  ## Note that here means is always FALSE to ensure correct formatting of the output data
  ## 'means' in the input to plot.pval only controls what's output to the plot
  results <- balanceMat(covariates, bal.out, means=means)

  ##print(results)
  ##return(results)

  ## set values of different parameters
  xlim = c(xlim1,1)
  pchset = c(21,24,22,23)
  
  ## Set the colors of the data points, if color is TRUE. If color is false, then the points are differentiated by whether they are colored in.
  if(color==TRUE)
    {
      pchcol = c("blue","red", "yellow", "darkgreen")
    }else{
      pchcol=c("black", "black", "black", "black")
    }
  if(color==TRUE)
    {
      pchbgcol = c("blue","red", "yellow", "darkgreen")
    }else{
      pchbgcol=c("white", "white", "black", "black")
    }


  # set margins and letter size
  par(cex=parcex, mai = c(0.5, 0.35, 1.1, 0.35))

  # set number of rows 
  ny = nrow(results)

  # create the empty figure
  if(!is.null(title))
    plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="", main=title)
  if(is.null(title))
    plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="")
  
  # add the 0, 0.05 and 0.1 vertical lines
  abline(v=c(0,0.05,0.1),lty=c(1,4,4), lwd=c(1,2,2))
  axis(side=1,at=c(0,0.05,0.1,1),tick=TRUE, las=2, cex.axis=0.7)

  # add labels on top of the three areas of the graph. Only add the means if desired
  # Can be useful to omit the means if the available plotting space is limited
  if(means==TRUE)
    axis(side=3,at=at1,labels="Mean\nTreated",tick=FALSE, padj=0.5,cex.axis=textsize)
  if(means==TRUE)
    axis(side=3,at=at2,labels="Mean\nControl",tick=FALSE, padj=0.5,cex.axis=textsize)

  axis(side=3,at=0.5,labels="P-values",tick=FALSE, padj=0.5,cex.axis=textsize)

  # Fill the figure with the information which is inside the 'results' matrix
  ## Add the p-values of the t-statistics as points
  if(means){
    t.cols <- c(3,4)
    ks.cols <- c(7,8)
  }else{
    t.cols <- c(1,2)
    ks.cols <- c(5,6)
  }
  
  #print(t.cols)
  #print(ks.cols)
  
  for(i in t.cols){
    points(results[,i],
           ny:1,
           pch = pchset[i-min(t.cols)+1],
           col = pchcol[i-min(t.cols)+1],
           bg = pchbgcol[i-min(t.cols)+1]
           )
  }
  ## Add the p-values of the ks statistics as points
  for(i in ks.cols){
    points(results[,i],
           ny:1,
           pch = pchset[i-min(ks.cols)+3],
           col = pchcol[i-min(ks.cols)+3],
           bg = pchbgcol[i-min(ks.cols)+3]
           )
  }

  #print("Points plotting successful")

  #print(dimnames(results))
  # Second, add each variable name and the means for treated and control
  for(i in 1:ny) {
    text(at3,ny-i+1,dimnames(results)[[1]][i],adj = 0,cex=textsize) # variable name
    if(means==TRUE)
      text(at1,ny-i+1,results[i,1], cex=textsize)        # treatment mean
    if(means==TRUE)
      text(at2,ny-i+1,results[i,2], cex=textsize)        # control mean
  }

  # Add dotted horizontal lines every two variables to make it prettier
  for(i in seq(2,by=2,length.out=floor((ny-1)/2)))
    abline(h = i+0.5, lty = 3)

  # Add legend
  if(legend==TRUE)
    legend(x=-1,y=0.5,
           c(colnames(results)[t.cols],
             colnames(results)[ks.cols]
             ),
           pch=pchset, pt.bg = pchbgcol,
           cex=0.8, ncol=2, xpd=NA
           )
  
}


#########################################
## END PLOTTING CODE
#########################################

## USE
## source(balStatsPlot.R)
## plot.pval(covariates, bal.out, title="")
