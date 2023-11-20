### from Wei Zhou
                                        # Some little R functions that might help

# Wrap long plot titles
wrap_strings <- function(vector_of_strings,width){
	sapply(vector_of_strings,FUN=function(x){paste(strwrap(x,width=width), collapse="\n")})
	}

# convert -log10(P) values to as.character(P)
log10toP <- function(log10P){
    log10P <- abs(as.numeric(log10P))
    if(is.na(log10P)) return(NA)
    if(log10P > 300){
        part1 <- log10P%/%100*100
        part2 <- log10P-part1
        P <- format(signif(10^-part2,3), scientific = T)
        P <- paste(as.numeric(gsub("e-.+","",P)),"E-",as.numeric(gsub(".+-","",P),sep="")+part1,sep="")
    } else {
        P <- signif(10^-log10P,3)
    }
    return(as.character(P))
}

qqplotdata <- function(logpvector,Lambda=NULL){
    denom<-qchisq(0.5, df=1)
    if(is.null(Lambda))  Lambda <- qchisq(10^-median(logpvector,na.rm=T), df=1, lower.tail = F) / denom
    #print(Lambda)
    
    o = sort(logpvector,decreasing=T)
    e = -log10(ppoints(length(o)))       
    
    qqdata <- data.frame(o,e)
    qqdata$o <- round(qqdata$o,3)
    qqdata$e <- round(qqdata$e,3)
    keepU <- which(!duplicated(qqdata))
    qqdata <- qqdata[keepU,]

    CHISQ <- qchisq(-abs(qqdata$o)*log(10), df=1, lower.tail = F,log.p=T)
    CHISQ_GC <- CHISQ/Lambda
    if(Lambda > 1) qqdata$o <- round(-pchisq(CHISQ_GC,df=1,lower.tail=F,log.p = T)/log(10),3)
    
    N <- length(logpvector) ## number of p-values
    ## create the confidence intervals
    qqdata$c95 <- NA
    qqdata$c05 <- NA

    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)

    for(k in 1:length(keepU)){
        j <- keepU[k]
        qqdata$c95[k] <- -log10(qbeta(1-0.05/2,j,N-j+1))
        qqdata$c05[k] <- -log10(qbeta(0.05/2,j,N-j+1))
    }
    
    return(list(qqdata,Lambda))
}

plotqq <- function(plotdata1,plotdata2,legendtext=NULL,maintitle="",break.top=16,top.size=0.25){	
	t_darkblue <- grDevices::rgb(t(grDevices::col2rgb("darkblue")), alpha=50, maxColorValue=255)
	t_darkgreen <- grDevices::rgb(t(grDevices::col2rgb("darkgreen")), alpha=50, maxColorValue=255)
	require(plotrix)
	y <- c(plotdata1$o,plotdata2$o)
	x <- c(plotdata1$e,plotdata2$e)
	#y <- c(plotdata2$o)	
	#x <- c(plotdata2$e)
	#cols <- c(rep("darkblue",dim(plotdata2)[1]))
	cols <- c(rep("darkgreen",dim(plotdata1)[1]),rep("darkblue",dim(plotdata2)[1]))
	xlim <- c(0,max(x,na.rm=T))
	ylim <- c(0,max(c(y,plotdata1$c05[1],plotdata2$c05[1]),na.rm=T))
	par(mar=c(5.1,5.1,4.1,1.1),las=1)

	if(ylim[2] > break.top+break.top*top.size){
		lab1 <- pretty(c(0,break.top),n=ceiling(12 * (1-top.size)))
		lab1 <- c(lab1[lab1 < break.top],break.top)
		lab2 <- pretty(c(break.top,max(y,na.rm=T)))
		top.range = break.top*(top.size)
		top.data = max(y)-break.top
		rescale = function(y) { break.top+(y-break.top)/(top.data/top.range)}
		rescaled.y <- ifelse(y>break.top,rescale(y),y)

		plot(0,0,
			ylim=c(min(y),break.top+break.top*top.size),xlim=xlim,axes=FALSE,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P[GC]),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
			main=maintitle,pch=19)

		polyY1 <- c(plotdata1$c95,plotdata1$c05[dim(plotdata1)[1]:1])
		polyY2 <- c(plotdata2$c95,plotdata2$c05[dim(plotdata2)[1]:1])	
    	polygon(c(plotdata1$e,plotdata1$e[dim(plotdata1)[1]:1]),
    		ifelse(polyY1>break.top,rescale(polyY1),polyY1),col=t_darkgreen,border = NA)
    	polygon(c(plotdata2$e,plotdata2$e[dim(plotdata2)[1]:1]),
    		ifelse(polyY2>break.top,rescale(polyY2),polyY2),col=t_darkblue,border = NA)

		points(x,rescaled.y,cex=1,col=cols,pch=19)

		axis(1,cex.axis=1.5,cex.lab=1.5)
		par(las=1)
		axis(side=2,at=c(lab1,rescale(lab2)),labels=c(lab1,lab2),cex.axis=1.5,cex.lab=1.5)

		box()
		axis.break(axis=2,breakpos=break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=break.top,style="zigzag",brw=0.02)
		lines(c(min(x),max(x)),c(break.top,break.top),col = "grey",lty = 6)
		lines(c(0,break.top),c(0,break.top),col="black",lty = 2)
		if(xlim[2] > break.top) lines(c(break.top,xlim[2]),c(break.top,rescale(xlim[2])),col="black",lty = 2)
	} else {
		plot(0,0,ylim=ylim,xlim=xlim,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P[GC]),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,col="transparent",
			main=maintitle,pch=19)

		polyY1 <- c(plotdata1$c95,plotdata1$c05[dim(plotdata1)[1]:1])
		polyY2 <- c(plotdata2$c95,plotdata2$c05[dim(plotdata2)[1]:1])
			
    	polygon(c(plotdata1$e,plotdata1$e[dim(plotdata1)[1]:1]),
    		polyY1,col=t_darkgreen,border = NA)
    	polygon(c(plotdata2$e,plotdata2$e[dim(plotdata2)[1]:1]),
    		polyY2,col=t_darkblue,border = NA)

		points(x,y,cex=1,col=cols,pch=19)
		lines(xlim,xlim,col="black",lty = 2)	
	}	
	if(!is.null(legendtext)) legend("topleft",legend=legendtext,col=c("darkgreen","darkblue"),pch=15,bty="n")
}

##merge overlap genomic regions
mergeOverlapRegion = function(data){
        data[,2][which(data[,2] <= 0)] = 1

        data_new=NULL
        a = data[1,]
        for (i in 2:nrow(data)){
                #cat("a0: ", a , "\n")
                #cat("data_new0: ", data_new , "\n")
                if (data[i,1] != a[1]){
                        data_new = rbind(data_new, a)
                        a = data[i,]
                }else{
                        if(data[i,2] <= a[3]){
                                a[3] = data[i,3]
                        }else{
                                data_new = rbind(data_new, a)
                                a = data[i,]
                        }
                }
                #cat("a: ", a , "\n")
                #cat("data_new: ", data_new , "\n")

        }
        data_new = rbind(data_new,a)
        return(data_new)
}
