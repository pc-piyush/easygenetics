#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


server <- function(input, output) {
  
  library(qtl) #if make into package, something we add in the documentation?
  
  #First Tab
observeEvent(input$go, {
  output$plot1 <- renderPlot({
    
    
    initial.n= input$initial.n
    growth= input$Growth
    f<- input$Rate
    num.generations<- input$Generations
    
    if (input$setseed == T) {
      set.seed(input$seed2)
    }
    
    for ( k in 1:length(initial.n)){
      for ( j in 1:length(f)){
        geno=rep(0,initial.n[k])
        sex=rbinom(initial.n[k],1,.5) # let's say sex=0 corresponds to male 
        allele.freq=rep(0,num.generations)
        allele.sum=rep(0,num.generations)
        allele.sum[1]=sum(geno)
        n=rep(0,num.generations)
        # SIMULATING GENOTYPES, CONT.
        for (i in (2:num.generations)){
          which.male=grep(T,sex==0)
          which.female=grep(T,sex==1) 
          num.matings=min(length(which.male),length(which.female)) 
          male.mate=sample(which.male,num.matings) 
          female.mate=sample(which.female,num.matings) 
          num.offspring=rpois(num.matings,growth) 
          n[i]=sum(num.offspring) 
          maternal.geno=rep(geno[female.mate],num.offspring) 
          paternal.geno=rep(geno[male.mate],num.offspring) 
          maternal.transmitted=rbinom(n[i],1,maternal.geno/2) 
          paternal.transmitted=rbinom(n[i],1,paternal.geno/2)
          # Mutation process
          na.m<-length(maternal.transmitted[maternal.transmitted==0]) 
          if (na.m >0 ) { 
            maternal.transmitted<-replace(maternal.transmitted, maternal.transmitted==0, rbinom(rep(1,na.m),1,f[j]))
          }           
          na.p<-length(paternal.transmitted[paternal.transmitted==0]) 
          if (na.p >0 ) { 
            paternal.transmitted<-replace(paternal.transmitted, paternal.transmitted==0, rbinom(rep(1,na.p),1,f[j]))
          }           
          geno=maternal.transmitted+paternal.transmitted 
          allele.sum[i]=sum(geno)
          sex=rbinom(n[i],1,.5) 
        }
        w=initial.n[k]+f[j]
        plot(n, main=paste0("Initial population Size:",  initial.n[k],"\n Mutation Rate: ", f[j]), xlab="Generation", ylab="Population", cex.axis=1.2, cex.lab=1.2)
        allele.freq=allele.sum/(2*n) 
        #plot(allele.freq,type="b", main=paste0("Initial population Size:",  initial.n[k],"\n Mutation Rate: ", f[j]), xlab="Generation", ylab="Allele Freq:", cex.axis=1.2, cex.lab=1.2)
      } 
    }
    
  })
    
  output$plot2 <- renderPlot({
    
    
    initial.n= input$initial.n
    growth= input$Growth
    f<- input$Rate
    num.generations<- input$Generations
    
    if (input$setseed == T) {
      set.seed(input$seed2)
    }
    
    for ( k in 1:length(initial.n)){
      for ( j in 1:length(f)){
        geno=rep(0,initial.n[k])
        sex=rbinom(initial.n[k],1,.5) # let's say sex=0 corresponds to male 
        allele.freq=rep(0,num.generations)
        allele.sum=rep(0,num.generations)
        allele.sum[1]=sum(geno)
        n=rep(0,num.generations)
        # SIMULATING GENOTYPES, CONT.
        for (i in (2:num.generations)){
          which.male=grep(T,sex==0)
          which.female=grep(T,sex==1) 
          num.matings=min(length(which.male),length(which.female)) 
          male.mate=sample(which.male,num.matings) 
          female.mate=sample(which.female,num.matings) 
          num.offspring=rpois(num.matings,growth) 
          n[i]=sum(num.offspring) 
          maternal.geno=rep(geno[female.mate],num.offspring) 
          paternal.geno=rep(geno[male.mate],num.offspring) 
          maternal.transmitted=rbinom(n[i],1,maternal.geno/2) 
          paternal.transmitted=rbinom(n[i],1,paternal.geno/2)
          # Mutation process
          na.m<-length(maternal.transmitted[maternal.transmitted==0]) 
          if (na.m >0 ) { 
            maternal.transmitted<-replace(maternal.transmitted, maternal.transmitted==0, rbinom(rep(1,na.m),1,f[j]))
          }           
          na.p<-length(paternal.transmitted[paternal.transmitted==0]) 
          if (na.p >0 ) { 
            paternal.transmitted<-replace(paternal.transmitted, paternal.transmitted==0, rbinom(rep(1,na.p),1,f[j]))
          }           
          geno=maternal.transmitted+paternal.transmitted 
          allele.sum[i]=sum(geno)
          sex=rbinom(n[i],1,.5) 
        }
        w=initial.n[k]+f[j]
        #plot(n, main=paste0("Initial population Size:",  initial.n[k],"\n Mutation Rate: ", f[j]), xlab="Generation", ylab="Population", cex.axis=1.2, cex.lab=1.2)
        allele.freq=allele.sum/(2*n) 
        plot(allele.freq,type="b", main=paste0("Initial population Size:",  initial.n[k],"\n Mutation Rate: ", f[j]), xlab="Generation", ylab="Allele Freq:", cex.axis=1.2, cex.lab=1.2)
      } 
    }
    
  })
})
  
  ######################################################
  #Second Tab Breeding value table
  
  output$estimates <- renderTable(bordered = T,{
    p <- input$p
    q = input$q
    #The midpoint m is calculated as:
    m = (input$A1A1 + input$A2A2) / 2; 
    
    # The additive effect is then:
    a <- ifelse(input$A1A1 > input$A2A2, input$A1A1 - m, input$A2A2 - m)
    # The dominacne effect is:
    d = input$A1A2 - m; 
    M = m + a*(p-q) + 2*d*p*q; 
    
    table2 <- rbind(
                    c("Midpoint (m)", m),
                    c("Population Mean (M)", M),
                    c("Additive Effect (a)", a),
                    c("Dominance Effect (d)",d))
    colnames(table2) <- c(" ", " ")
    table2
    
  })
  
  output$bvplot <- renderTable(spacing = "l" ,
                               bordered = T,
                               width = '100%',
                               {
    
    # Problem states that:
    p = input$p
    q = input$q
    
    #The midpoint m is calculated as:

    m = (input$A1A1 + input$A2A2) / 2; 
    
    # The additive effect is then:
    
    a <- ifelse(input$A1A1 > input$A2A2, input$A1A1 - m, input$A2A2 - m)
    # The dominacne effect is:
    d = input$A1A2 - m; 
    
    #The overall population mean is:
    M = m + a*(p-q) + 2*d*p*q; 
    
    
    # The average genetic effect is:
    alpha_1 = q*(a+d*(q-p)); # for allele a
    alpha_2 = -p*(a+d*(q-p)); # for allele A
    
    alpha_diff = alpha_1 - alpha_2; 
    
    # Breeding value is:
    bv_aa = 2*q*alpha_diff; 
    bv_Aa = (q-p)*alpha_diff; 
    bv_AA = -2*p*alpha_diff; 
    
    # Dominance deviation is:
    dd_aa = -2*q^2*d;
    dd_Aa = 2*p*q*d;
    dd_AA = -2*p^2*d;
    
    # Table:

    mytable<- rbind(c("Frequency", p^2, 2*p*q, q^2), 
                    c("Genotypic Value", a, d, -1*a), 
                    c("Breeding Value", bv_aa, bv_Aa, bv_AA), 
                    c("Dominance Deviation", dd_aa, dd_Aa, dd_AA))
    colnames(mytable) <- c("Genotype", "A1A1", "A1A2", "A2A2")
   
    
   # colnames(mytable) <- ifelse(input$aa > input$AA, c("Genotype", "aa", "Aa", "AA"),
  #                              c("Genotype", "AA", "Aa", "aa"))
    
    mytable
    
    
  })
  
  ##################################
  
  #Parent-Child tab:
  observeEvent(input$runcov, {
  output$PCCov <- renderText({
    PCCov <- function(n)
    {
      if (input$setseedcov == T) {
      set.seed(input$seedcov)
      }
      p1 <- rnorm(n,0,1) ##parental allelic 1 effect
      p2 <- rnorm(n,0,1) ##parental allelic 2 effect
      which.allele <- rbinom(n,1,0.5) # transmitted
      y1 <- p1 + p2  ##Parental additive genetic effect
      y2 <- p2 ## Child additive effect
      y2[which.allele == 1] <- p1[which.allele == 1]
      
      y1 <- y1 + rnorm(n,0,.5)  ##add environmental error: parental phenotype
      y2 <- y2 + rnorm(n,0,1) + rnorm(n,0,.5)  ##Child phenotype, another parent and environment
      
      cov1 <- cov(y1,y2)
      
      print(paste("The covariance between parents and offspring is", cov1))
      
      
    }
    
    PCCov(input$int)
  })
  
  })
  
  #####################################
  #Regression Tab:
  observeEvent(input$runhe, {
    output$HE <- renderPlot({
      
      HE <- function(plot,num.families,var.total,var.genetic,var.env,maf)
      {
        if (input$setseedhe == T) {
        set.seed(input$seedhe)
        }
        var.error = var.total - var.genetic - var.env
        father.effect1 = rnorm(num.families,0,sd=sqrt(var.genetic))
        father.effect2 = rnorm(num.families,0,sd=sqrt(var.genetic))
        mother.effect1 = rnorm(num.families,0,sd=sqrt(var.genetic))
        mother.effect2 = rnorm(num.families,0,sd=sqrt(var.genetic))
        enviro.effect = rnorm(num.families,0,sd=sqrt(var.env))
        error.sib1 = rnorm(num.families,0,sd=sqrt(var.error))
        error.sib2 = rnorm(num.families,0,sd=sqrt(var.error))
        
        
        from.father.sib1=rbinom(num.families,1,maf)+1
        from.mother.sib1=rbinom(num.families,1,maf)+1
        from.father.sib2=rbinom(num.families,1,maf)+1
        from.mother.sib2=rbinom(num.families,1,maf)+1
        father.effect.sib1=father.effect1*(from.father.sib1==1)+father.effect2*(from.father.sib1==2)
        mother.effect.sib1=mother.effect1*(from.mother.sib1==1)+mother.effect2*(from.mother.sib1==2)
        father.effect.sib2=father.effect1*(from.father.sib2==1)+father.effect2*(from.father.sib2==2)
        mother.effect.sib2=mother.effect1*(from.mother.sib2==1)+mother.effect2*(from.mother.sib2==2)
        pheno.sib1=father.effect.sib1+mother.effect.sib1+enviro.effect+error.sib1
        pheno.sib2=father.effect.sib2+mother.effect.sib2+enviro.effect+error.sib2
        IBD=(from.father.sib1==from.father.sib2)+(from.mother.sib1==from.mother.sib2)
        
        if (plot=="All")
        {  
          par(mfrow=c(2,3))
          hist(IBD)
          plot(IBD,(pheno.sib1-pheno.sib2)^2,pch=16,cex=.8)
          plot(pheno.sib1,pheno.sib2,pch=16,cex=.8)
          plot(pheno.sib1[IBD==0],pheno.sib2[IBD==0],pch=16,cex=.8)
          plot(pheno.sib1[IBD==1],pheno.sib2[IBD==1],pch=16,cex=.8)
          plot(pheno.sib1[IBD==2],pheno.sib2[IBD==2],pch=16,cex=.8)
        }
        
        
        else if (plot=="Histogram of IBD")
        {
          hist(IBD) 
        }
        else if (plot=="IBD vs Differences in Sibling Phenotypes")
        {
          plot(IBD,(pheno.sib1-pheno.sib2)^2,pch=16,cex=.8)
        }
        else if (plot=="Phenotype- Sibling 1 vs Sibling 2")
        {
          plot(pheno.sib1,pheno.sib2,pch=16,cex=.8)
        }
        
        
        else if (plot=="Phenotype IBD=0 Sibling 1 vs Sibling 2")
        {
          plot(pheno.sib1[IBD==0],pheno.sib2[IBD==0],pch=16,cex=.8)
        }
        else if (plot=="Phenotype IBD=1 Sibling 1 vs Sibling 2")
        {
          plot(pheno.sib1[IBD==1],pheno.sib2[IBD==1],pch=16,cex=.8)
          
        }
        else if (plot=="Phenotype IBD=2 Sibling 1 vs Sibling 2")
        {
          plot(pheno.sib1[IBD==2],pheno.sib2[IBD==2],pch=16,cex=.8)
        }
        
        #print(paste("The covariance between parents and offspring is", cov1))
      }
      
      
      HE(input$plothe, input$NumFamhe, input$TotVarhe, input$GenVarhe, input$EnvVarhe,input$mafhe)
    })
  })
    
    
    ######################################################
    #Multiple QTL Tab

    mQTL <- function(plot,len,nmark,space,type,pos1,eff1,pos2,eff2,ind)
    {
      if (input$setseedqtl == T) {
      set.seed(input$seedqtl)
      }
      tel <- FALSE
      #map1 <- sim.map(len, nmark, include.x=tel, eq.spacing=space) 
      if (tel=="FALSE"){
        if (space=="TRUE")
        {
          map1 <- sim.map(len, nmark, include.x=FALSE, eq.spacing=TRUE) 
        }
        else if (space=="FALSE")
        {
          map1 <- sim.map(len, nmark, include.x=FALSE, eq.spacing=FALSE) 
        }
      }
      
      
      met <- "hk"
      chr1 <- 1
      chr2 <- 1
      
      if (type=="f2")
      {
        fake1 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1,0)))
        fake1 <- calc.genoprob(fake1, step=1)
        tempone1 <- scanone(fake1, chr=chr1, pheno.col=1, method=met)
        qtl1 <- makeqtl(fake1, chr=chr1, pos=pos1, what='prob')
        getfit1 <- fitqtl(fake1, qtl=qtl1, get.ests=TRUE, method=met)
        
        fake2 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1,0), c(chr2,pos2,eff2,0)))
        fake2 <- calc.genoprob(fake2, step=1)
        tempone2 <- scanone(fake2, chr=chr1, pheno.col=1, method=met)
        temptwo <- scantwo(fake2, chr=chr1, pheno.col=1, method=met)
        #temptwo <- scantwo(fake2, chr=chr1, pheno.col=1, method=met)
        marker <- find.marker(fake2, chr=chr1, pos=pos1)
        qtl2 <- makeqtl(fake2, chr=c(chr1,chr2), pos=c(pos1,pos2), what='prob')
        getfit2 <- fitqtl(fake2, qtl=qtl2, get.ests=TRUE, method=met)
        
        fake3 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1,0), c(chr2,pos2,eff2,0)))
        fake3 <- calc.genoprob(fake3, step=1)
        tempone3 <- scanone(fake3, chr=chr1, pheno.col=1, method=met)
        mar <- find.marker(fake3, chr1, pos1)
        g <- pull.geno(fake3)[,mar]
        any(is.na(g))
        out2 <- scanone(fake3, addcovar=g)
      }
      
      
      
      if (type=="bc")
      {
        fake1 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1)))
        fake1 <- calc.genoprob(fake1, step=1)
        tempone1 <- scanone(fake1, chr=chr1, pheno.col=1, method=met)
        qtl1 <- makeqtl(fake1, chr=chr1, pos=pos1, what='prob')
        getfit1 <- fitqtl(fake1, qtl=qtl1, get.ests=TRUE, method=met)
        
        fake2 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1), c(chr2,pos2,eff2)))
        fake2 <- calc.genoprob(fake2, step=1)
        tempone2 <- scanone(fake2, chr=chr1, pheno.col=1, method=met)
        temptwo <- scantwo(fake2, chr=chr1, pheno.col=1, method=met)
        #temptwo <- scantwo(fake2, chr=chr1, pheno.col=1, method=met)
        marker <- find.marker(fake2, chr=chr1, pos=pos1)
        qtl2 <- makeqtl(fake2, chr=c(chr1,chr2), pos=c(pos1,pos2), what='prob')
        getfit2 <- fitqtl(fake2, qtl=qtl2, get.ests=TRUE, method=met)
        
        fake3 <- sim.cross(map1, type=type, n.ind=ind, model=rbind(c(chr1,pos1,eff1), c(chr2,pos2,eff2)))
        fake3 <- calc.genoprob(fake3, step=1)
        tempone3 <- scanone(fake3, chr=chr1, pheno.col=1, method=met)
        mar <- find.marker(fake3, chr1, pos1)
        g <- pull.geno(fake3)[,mar]
        any(is.na(g))
        out2 <- scanone(fake3, addcovar=g)
      }
      
      if (plot=="Genetic Map")
      {
        plot(map1)
      }
      
      else if (plot=="Single QTL")
      {
        # Single QTL, Large effect
        plot(tempone1)
        abline(v=c(pos1))
        summary(getfit1)
      }
      
      
      else if (plot=="Two QTLs")
      {
        plot(tempone2)
        abline(v=c(pos1,pos2))
        summary(getfit2)
      }
      
      
      else if (plot=="Heat Map-Two QTLs Comparison")
      {
        plot(temptwo)
      }
      
      else if (plot=="Two QTLs Comparison")
      {
        
        plotPXG(fake2, marker=marker)
        summary(temptwo)
      }
      
      
      
      else if (plot=="Composite Interval Mapping")
      {
        
        plot(tempone3, out2)
        summary(out2)
      }
      
    }
    
    observeEvent(input$runqtl, {
    output$mQTLP <- renderPlot({
      
      mQTL(input$plotqtl, input$len, input$nmark, input$space,input$type,
           input$pos1, input$eff1, input$pos2,input$eff2,
           input$ind)
    })
    
    output$mQTLT <- renderPrint({
      mQTL(input$plotqtl, input$len, input$nmark, input$space,input$type,
           input$pos1, input$eff1, input$pos2,input$eff2,
           input$ind)
   })
    })
    
    ########################################################
    #Interval_Mapping tab
    
    IntMap <- function(plot,len,nmark,space,type,qtl,met,perm,pos1,eff1,pos2,eff2,ind)
    {
      if (input$setseedim == T) {
      set.seed(input$seedim)
      }
      tel <- FALSE
      #map1 <- sim.map(len, nmark, include.x=tel, eq.spacing=space) 
      if (tel=="FALSE"){
        if (space=="TRUE")
        {
          map1 <- sim.map(len, nmark, include.x=FALSE, eq.spacing=TRUE) 
        }
        else if (space=="FALSE")
        {
          map1 <- sim.map(len, nmark, include.x=FALSE, eq.spacing=FALSE) 
        }
      }
      
      chr1 <- 1
      
      if (qtl==1)
      {
        if (type=="f2")
        {
          fakes <- sim.cross(map1, type="f2", n.ind=ind, model=rbind(c(chr1,pos1,eff1,0)))
          fake1 <- calc.genoprob(fakes, step=1)
        }
        
        else if (type=="bc")
        {
          fakes <- sim.cross(map1, type="bc", n.ind=ind, model=rbind(c(chr1,pos1,eff1)))
          fake1 <- calc.genoprob(fakes, step=1)
        }
      }
      
      
      if (qtl==2)
      {
        if (type=="f2")
        {
          fakes <- sim.cross(map1, type="f2", n.ind=ind, model=rbind(c(chr1,pos1,eff1,0),c(chr1,pos2,eff2,0)))
          fake1 <- calc.genoprob(fakes, step=1)
        }
        
        else if (type=="bc")
        {
          fakes <- sim.cross(map1, type="bc", n.ind=ind, model=rbind(c(chr1,pos1,eff1),c(chr1,pos2,eff2)))
          fake1 <- calc.genoprob(fakes, step=1)
        }
      }
      
      
      tempone1 <- scanone(fake1,method="em")#,n.perm = perm)
      tempone2 <- scanone(fake1,method="hk")#,n.perm = perm)
      tempone3 <- scanone(fake1,method="ehk")#,n.perm = perm)
      
      permone1 <- suppressMessages(scanone(fake1,method="em",n.perm = perm)) 
      permone2 <- suppressMessages(scanone(fake1,method="hk",n.perm = perm)) 
      permone3 <- suppressMessages(scanone(fake1,method="ehk",n.perm = perm)) 
      
      
      if (plot=="Properties of Cross")
      {
        print(summary(fakes))
      }
      
      t1 <- summary(tempone1)
      t2 <- summary(tempone2)
      t3 <- summary(tempone3)
      
      if (plot=="Plot")
      {
        if (met=="EM")
        {
          plot(tempone1,main="EM Interval Mapping", col="Red")
          abline( h = c(summary(permone1)[1,1]), col=c("Red"), lty=2)
          legend("topright",legend=c("EM"),col = c("Red"), lty = 1)
          
        }
        
        if (met=="HK")
        {
          plot(tempone2,main="HK Interval Mapping", col="Blue")
          abline( h = c(summary(permone2)[1,1]), col=c("Blue"), lty=2)
          legend("topright",legend=c("HK"),col = c("Blue"), lty = 1)
        }
        
        if (met=="EHK")
        {
          plot(tempone3,main="Extended HK Interval Mapping", col="Green")
          abline( h = c(summary(permone3)[1,1]), col=c("Green"), lty=2)
          legend("topright",legend=c("EHK"),col = c("Green"), lty = 1)
        }
        
        if (met=="All")
        {
          plot(tempone1,tempone2,tempone3,main="EM, HK and Extended HK Interval Mapping", col=c("Red","Blue","Green"))
          abline( h = c(summary(permone1)[1,1],summary(permone2)[1,1],summary(permone3)[1,1]), col=c("Red","Blue","Green"), lty=2)
          legend("topright",legend=c("EM","HK","EHK"),col = c("Red","Blue","Green"), lty = 1)
        }
        
      }
      
      if (plot=="LOD Thresholds")
      {
        if (met=="EM")
        {
          print("EM Interval Mapping Threshold")
          print(summary(permone1))
        }
        
        if (met=="HK")
        {
          print("HK Interval Mapping Threshold")
          print(summary(permone2))
        }
        
        if (met=="EHK")
        {
          print("EHK Interval Mapping Threshold")
          print(summary(permone3))
        }
        
        if (met=="All")
        {
          print("EM Interval Mapping Threshold")
          print(summary(permone1))
          print("HK Interval Mapping Threshold")
          print(summary(permone2))
          print("EHK Interval Mapping Threshold")
          print(summary(permone3))
        }
        
      }
      
      if (plot=="LOD Intervals")
      {
        if (met=="EM")
        {
          print("EM Interval Mapping LOD Score Interval")
          print(lodint(tempone1,chr = 1))
        }
        
        if (met=="HK")
        {
          print("HK Interval Mapping LOD Score Interval")
          print(lodint(tempone2,chr = 1))
        }
        
        if (met=="EHK")
        {
          print("EHK Interval Mapping LOD Score Interval")
          print(lodint(tempone3,chr = 1))
        }
        
        if (met=="All")
        {
          print("EM Interval Mapping LOD Score Interval")
          print(lodint(tempone1,chr = 1))
          print("HK Interval Mapping LOD Score Interval")
          print(lodint(tempone2,chr = 1))
          print("EHK Interval Mapping LOD Score Interval")
          print(lodint(tempone3,chr = 1))
        }
        
      }
      
    }
    
    observeEvent(input$runim, {
    output$IntMapP <- renderPlot({
      IntMap(input$plotim, input$lenim, input$nmarkim, input$spaceim,input$typeim,
             input$qtlim,input$metim,input$permim,input$pos1im, input$eff1im, input$pos2im,input$eff2im,
             input$indim)
    })
    
    output$IntMapT <- renderPrint({
      IntMap(input$plotim, input$lenim, input$nmarkim, input$spaceim,input$typeim,
             input$qtlim,input$metim,input$permim,input$pos1im, input$eff1im, input$pos2im,input$eff2im,
             input$indim)
    })
    })

  
  
  

  
  
  
   
}

#shinyApp(ui=ui, server = server)