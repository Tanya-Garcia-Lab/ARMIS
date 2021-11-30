#' Analysis of the repeated measured data in Animal Science
#'
#' @param data dataset you want to analyze whose response variables are repeatedly measured
#' @param num.aa number of amino acids, the default is 1
#' @param n number of animals
#' @param time.points number of level of time factor
#' @param subid subject id in your data set
#' @param group group name in your data set
#' @param resp.var repeated measured response variable in your data set (repeated measured)
#' @param amino.names name of the amino acid (response variables)
#' @param name.tt name of the time vector in plot 
#' @param name.steer name of the steer in plot
#' @param time time points at which the response variable are measured at
#' @param interv.length interval length between two grids to draw a y-axis
#' @param num.method number of ways to fit linear mixed effects model
#'                 1) with standard var-cov structure,
#'                 2) with heterogeneous variance and generalized AR(1),
#'                 3) with baseline covariates,
#'                 4) based on normalized data using the percentage-changes from baselines
#' @param corStruct  choice of correlation structures, e.g) "ar1", "gen.ar1", "compsymm", and "unstruct"
#' @param hetero TRUE for unequal variance for the within-subject errors.
#' @param file  data file name: .cvs file.  The data file format has to be csv file and long format with four columns: subject_id, group_id, Time and Outcomes for the repeated measured data
#'
#' @return one table of p-values to show whether time effect is statistically significant from three mixed effects model 
#'          and three figures: comparison of the original values of citrulline (amino acid) concentration over time and 
#'          the normalized values of citrulline using percentage based transformation. It shows that our strategies 
#'          capture the huge biological variations among all study subjects and account for the correlations 
#'          between measurements within individual subjects. 
#'
#' @importFrom nlme lme
#'
#' @examples
#'
#' library(ARMIS)
#' 
#' #generate dataset
#' data(pseudo_steer)
#' head(pseudo_steer)
#'
#' #specify the input parameters
#' data=data;
#' num.aa=1; n=6; time.points=6;
#' subid="Steer"; group="Group"; time="Time"; resp.var="citrulline"; amino.names="citrulline";
#' name.tt=c("time0", "time1", "time2", "time3", "time4", "time5");
#' name.steer=c("steer1","steer2","steer3","steer4","steer5","steer6");
#' interv.length=7; num.method=4; corStruct="gen.ar1"; hetero=TRUE;
#' file="cit_data.csv"
#' 
#' # analyze the data using anova.test() function
#' result<-anova.test(data=data, num.aa=1, n=6, time.points=6, subid="Steer",
#' group="Group", time="Time", resp.var="citrulline", amino.names="citrulline",
#' interv.length=7, num.method=4, corStruct="gen.ar1", hetero=TRUE,
#' file="rpaa_cit_data.csv")

anova.test<-function(data=data, num.aa=1, n=6, time.points=6, subid="Steer",
                     group="Group", time="Time", resp.var="citrulline", amino.names="citrulline",
                     name.tt=c("time0", "time1", "time2", "time3", "time4", "time5"),
                     name.steer=c("steer1","steer2","steer3","steer4","steer5","steer6"),
                     interv.length=7, num.method=4,
                     corStruct="gen.ar1", hetero=TRUE, file="rpaa_cit_data.csv"){

  #######################################################################################
  ## ###### repeated measured steer data ################################################
  ##
  #######################################################################################


  ## read the steer data
  steer.id<- data[,subid]
  group.name<-data[,group]
  tt<-data[,time]

  AA_list<-vector(mode='list', length=num.aa)
  for (i in 1:num.aa){

  AA_list[[i]]<-data[,resp.var[i]]
  }

  AA<-as.data.frame(do.call(cbind, AA_list))
  #names(AA)<-amino.names
  
  steer<-cbind(steer.id, group.name, tt, AA)
  names(steer)<-c("steer.id", "group", "tt", "true.AA")
  
  ## dummy variables for steers and time
  ## steers are six distinct steers
  ## time factor has 6 levles (t=0, 0.5, 1, 2, 4, 6),
  ##                  where 0 is the reference level

  tag<-factor(steer[,"steer.id"])
  time<-factor(steer[,"tt"])

  ## preparation steer data in the data.frame format for the analysis
  ## tag is the categorical variable with 6 distinct steers
  ## time is the catergorical variables with 6 levels

  data.steer<-as.data.frame(cbind(steer, tag, time))
  #head(data.steer)
  #names(data.steer)[3+num.aa]<-amino.names

  ########################################################
  # method 1 : Modeling variance-covariance structure   ##
  ########################################################

  ##{{

  ############################################
  ## initialization for the array of output ##
  ############################################



  lme.results<-array(0,dim=c(2, 1, num.aa), dimnames=list( c("intercept", "Time"),
                                                           c("P-value"), amino.names))

  lme.cov.results<-array(0,dim=c(2, 1, num.aa), dimnames=list( c("intercept", "Time"),
                                                               c("P-value"),amino.names ))

  ## model.select

  for(i in 1:num.aa){

    ########################################################
    # mixed effects model with Standard var-cov structure ##
    ########################################################

    model0<-lme(as.numeric(true.AA)~as.numeric(tt), data=data.steer, random=~1|tag)
    lme.anova0<-anova(model0)

    ## summary for anova results -with constant
    lme.anova.summary0<-lme.anova0[,"p-value"]
    lme.results[,,i]<-round(lme.anova.summary0, 3)



    ############################################################
    # Veriety of var-covariance structure can be implemented  ##
    ############################################################
    ###########################################################################################################
    ## To have unequal variances,                                                                           ###
    ## varIdent: variance function structure allows different variances for each level of a factor          ###
    ##           it can be used to fit heteroscedastic model                                                ###
    ##                                                                                                      ###
    ## To choose various variance-covariance structure,                                                     ###
    ## 1) ar1 : autoregressive with order 1 (AR(1))                                                         ###
    ## 2) gen.ar1:generalized AR(1),                                                                        ###
    ##            correlation decreases depending on distance between two observations from the same steer  ###
    ##            two observations are closer in time have higher correlation than those are further in time###
    ## 3) compsymm : compound symmetric covariance structure                                                ###
    ## 4) unstruct : unstructured covariance structure                                                      ###
    ## (not recommended if there is not enough number of data)                                              ###
    ###########################################################################################################



    if (corStruct=="ar1" && hetero==TRUE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corAR1(form = ~1|tag),  data=data.steer)

    }else if(corStruct=="gen.ar1" && hetero==TRUE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corExp(form = ~1|tag , nugget=TRUE), data=data.steer)

    }else if(corStruct=="compsymm"&& hetero==TRUE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corCompSymm(form = ~1|tag), data=data.steer)

    }else if(corStruct=="unstruct"&& hetero==TRUE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corSymm(form = ~1|tag ), data=data.steer)

    }else if(corStruct=="ar1"&& hetero==FALSE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag,
                  correlation=corAR1(form = ~1|tag), data=data.steer)

    }else if(corStruct=="gen.ar1"&& hetero==FALSE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag,
                  correlation=corExp(form = ~1|tag , nugget=TRUE), data=data.steer)

    }else if(corStruct=="compsymm"&& hetero==FALSE){

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag,
                  correlation=corCompSymm(form = ~1|tag), data=data.steer)

    }else{

      model1<-lme(as.numeric(true.AA)~as.numeric(tt), random=~1|tag,
                  correlation=corSymm(form = ~1|tag ), data=data.steer)
    }



    lme.anova1<-anova(model1)


    ## summary for anova results
    lme.anova.summary1<-lme.anova1[,"p-value"]
    lme.cov.results[,,i]<-round(lme.anova.summary1, 3)


    ##########################################################
    # model selection: standard cov vs. generalized AR(1)   ##
    ##########################################################
    model.select<-anova(model0, model1)

  }


  ###}}



  #############################################
  # method 2 : use baselines as  covariates  ##
  #############################################


  ## unique steer with tag number
  uniq.steer<-as.numeric(unique(data.steer[,"tag"]))


  ###{{

  ###############################################################
  # initalization list for dataset with baseline covariates    ##
  ###############################################################
  names(data.steer)[3+num.aa]<-amino.names
  trim.steer.base.cov<-vector("list", n)



  for(id in 1:length(uniq.steer)){

    subdata<-data.steer[data.steer[,"steer.id"]==uniq.steer[id],]
    base.amino<-matrix(100, nrow=(time.points-1), ncol=num.aa)

    ### get data structure for response variables at time t=0.5, 1, 2, 4, 6 and baseline covariates

    for (j in 1:num.aa){

      base.amino[,j]<-rep(subdata[ ,amino.names[j]][subdata[,"tt"]==0], (time.points-1))

    }

    trim.subdata<-subdata[-1,]
    colnames(base.amino)<-paste("base.", amino.names, sep="")
    trim.subdata.base.amino<-cbind(trim.subdata, base.amino)


    trim.steer.base.cov[[id]]<-trim.subdata.base.amino

  }


  ## create data.frame from list type of data
  trim.steer<-do.call("rbind.data.frame", trim.steer.base.cov)
  names(trim.steer)<-c("steer.id", "group.name", "tt", "AA", "tag", "time", "base.AA")



  ###########################################
  ## initalization for the array of output ##
  ###########################################



  lme.base.cov.results<-array(0,dim=c(3, 1, num.aa), dimnames=list( c("intercept", "baseline", "time"),
                                                                    c("p-value"), amino.names))


  ##AR(1) : ar1
  ##generalized AR(1): gen.ar1
  ##CompoundSymm: compsymm
  ##Unstructured: unsruct


  for(i in 1:num.aa){

    ## Fit model with baseline covariate with stpecific var-covariance structure


    if (corStruct=="ar1" && hetero==TRUE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corAR1(form = ~1|tag), data=trim.steer)
    }else if(corStruct=="gen.ar1" && hetero==TRUE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corExp(form = ~1|tag , nugget=TRUE), data=trim.steer)
    }else if(corStruct=="compsymm"&& hetero==TRUE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corCompSymm(form = ~1|tag), data=trim.steer)
    }else if(corStruct=="unstruct"&& hetero==TRUE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag, weights=varIdent(form=~1|tag),
                  correlation=corSymm(form = ~1|tag ), data=trim.steer)
    }else if(corStruct=="ar1"&& hetero==FALSE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag,
                  correlation=corAR1(form = ~1|tag), data=trim.steer)
    }else if(corStruct=="gen.ar1"&& hetero==FALSE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag,
                  correlation=corExp(form = ~1|tag , nugget=TRUE), data=trim.steer)
    }else if(corStruct=="compsymm"&& hetero==FALSE){
      model2<-lme(AA~1+base.AA+time, random=~1|tag,
                  correlation=corCompSymm(form = ~1|tag), data=trim.steer)
    }else{
      model2<-lme(AA~1+base.AA+time, random=~1|tag,
                  correlation=corSymm(form = ~1|tag ), data=trim.steer)
    }

    lme.anova2<-anova(model2)


    lme.anova.summary2<-lme.anova2[,"p-value"]
    lme.base.cov.results[,,i]<-round(lme.anova.summary2, 3)

  }

  ###}}





  ###################################################
  # method 3 : Percentage-changes from baselines   ##
  ###################################################

  ##{{

  #############################
  # Normalization for data   ##
  #############################

  #uniq.steer<-unique(steer$tag)
  no.base.steer.list<-vector("list", n)
  uniq.steer.list<-vector("list", n)

  for(id in 1:length(uniq.steer)){


    names(data.steer)[3+num.aa]<-amino.names
    subdata<-data.steer[data.steer[,"steer.id"]==uniq.steer[id],]

    ################################################
    # Initialize matrix for normalized amino acids
    ################################################

    normamino<-matrix(100, nrow=time.points, ncol=num.aa, dimnames=list(name.tt, amino.names))

    #########################################
    ## percentage-chages from baselines   ###
    #########################################

    for (j in 1:num.aa){


      ## Percentage-changes from baselines : divide all response variables by baselines

      baseamino<-subdata[ ,amino.names[j]][subdata[,"time"]==0]
      normamino[,j]<-{{(subdata[ ,amino.names[j]]/baseamino)}*100}-100
      subdata[,c(paste("norm.", amino.names, sep=""))]<-normamino[,j]

    }


    ## normalized baselines are all 0
    uniq.steer.list[[id]]<-subdata

    ## exclude the normalized baselines in the analysis
    no.base.subdata<-subdata[-1,]
    no.base.steer.list[[id]]<-no.base.subdata

  }


  #####################################################################
  ## noramlized data set by using percentage-changes from baselines  ##
  #####################################################################

  nobase.trans.steer<-do.call("rbind.data.frame", no.base.steer.list)

  ## exclude the original scaled response varibles
  trans.steer<-nobase.trans.steer[,-((3+1):(3+num.aa))]

  ## obtain specific name for normalized amino acids to be analysis
  norm.amino.names<-names(trans.steer)[6]

  names(trans.steer)<-c("steer.id", "group.name", "tt", "tag", "time", "norm.AA")
  #####################
  ## array of output ##
  #####################


  lme.percent.results<-array(0,dim=c(1, 1, num.aa), dimnames=list( c("Time"),
                                                                   c("P-value"), norm.amino.names))


  for(i in 1:num.aa){

    ## fit mixed effecs model to percentage-based from baselines
    model3<-lme(norm.AA~-1+time, data=trans.steer, random= ~1|tag )
    lme.anova3<-anova(model3)


    lme.anova.summary3<- lme.anova3[,"p-value"]
    lme.percent.results[,,i]<-round(lme.anova.summary3, 3)


  }


  ###}}



  #######################################################################################
  ## Combine all analysis results using linear mixed effects model                     ##
  ## 1) with standard covariance structure                                             ##
  ## 2) with heterogeneous variances and generalized AR(1) covariance structure        ##
  ## 3) with baselines covariates                                                      ##
  ## 4) with normalized data based on percentage-changes for the baselines             ##
  #######################################################################################


  ##{{ table 1

  combine.results<-array(0,dim=c(num.method, 1, 1), dimnames=list( c("lme.stand", "lme.cov", "lme.base.covr", "lme.precent")
                                                                   , c("p-value"),"Time"))

  ## obtain restuls  only for the time factor
  combine.results["lme.stand",,"Time"]<-lme.results[2,,amino.names]
  combine.results["lme.cov",,"Time"]<-lme.cov.results[2,,amino.names]
  combine.results["lme.base.covr",, "Time"]<-lme.base.cov.results[3,,amino.names]
  combine.results["lme.precent",,"Time"]<-lme.percent.results[1,,  norm.amino.names]

  zero.pval<-which(combine.results[,"p-value", "Time"]==0)

  combine.results[,"p-value", "Time"][c(zero.pval)]<-format.pval(combine.results[,"p-value", "Time"][c(zero.pval)], eps=.001, ditis=2)

  ### return results without quotation mark on p-value
  print(combine.results, quote=FALSE)

  ###}}






  ##################################################################################################
  ## This part is to test                                                                        ###
  ## if the rate of concentrations of citrulline at each time point is increasing.               ###
  ## Use independent t-test                                                                      ###
  ## First, produce summary statistics including averages and standar errors for each steer      ###
  ## Second, use independent t-test in order to test                                             ###
  ##         the average concentration of citrulline is equal to 0:                              ###
  ##                         H_0: mu_{t}=0 at each time point t                                  ###
  ##                         H_1: mu_{t}>0                                                       ###
  ##################################################################################################

  ###{{

  #######################################################################################
  ## Initialization to have data summary: average and standard errors for each steer   ##
  #######################################################################################
  trim.tt<-unique(tt)[-1]


  percent.avr.at.times<-array(0,dim=c(2, 5, num.aa), dimnames=list( c( "avr", "sterr"),
                                                                    trim.name.tt, norm.amino.names ))

  percent.steer.at.times<-array(0,dim=c(6, 5, num.aa), dimnames=list( c( name.steer), trim.name.tt, norm.amino.names))

  for(i in 1:num.aa){

    aa<-i+5


    second.pt<-which(trans.steer[,"time"]==trim.tt[1])
    steer.at.second.pt<-trans.steer[,aa][second.pt]
    avr.at.second.pt<-mean( steer.at.second.pt)
    str.at.second.pt<-sd(steer.at.second.pt)/sqrt(n)

    third.pt<-which(trans.steer[,"time"]==trim.tt[2])
    steer.at.third.pt<-trans.steer[,aa][third.pt]
    avr.at.third.pt<-mean(steer.at.third.pt)
    str.at.third.pt<-sd(steer.at.third.pt)/sqrt(n)

    fourth.pt<-which(trans.steer[,"time"]==trim.tt[3])
    steer.at.fourth.pt<-trans.steer[,aa][fourth.pt]
    avr.at.fourth.pt<-mean(steer.at.fourth.pt)
    str.at.fourth.pt<-sd(steer.at.fourth.pt)/sqrt(n)

    fifth.pt<-which(trans.steer[,"time"]==trim.tt[4])
    steer.at.fifth.pt<-trans.steer[,aa][fifth.pt]
    avr.at.fifth.pt<-mean(steer.at.fifth.pt)
    str.at.fifth.pt<-sd(steer.at.fourh)/sqrt(n)

    sixth.pt<-which(trans.steer[,"time"]==trim.tt[5])
    steer.at.sixth.pt<-trans.steer[,aa][sixth.pt]
    avr.at.sixth.pt<-mean(steer.at.sixth.pt)
    str.at.sixth.pt<-sd(steer.at.sixth.pt)/sqrt(n)

    ### average of outcomes at each time points ###
    avr<-cbind(avr.at.second.pt, avr.at.third.pt, avr.at.fourth.pt,  avr.at.fifth.pt, avr.at.sixth.pt)

    ### standard errors of outcomes at each time points ###
    se<-cbind(str.at.second.pt,  str.at.third.pt, str.at.fourth.pt,  str.at.fifth.pt, str.at.sixth.pt)


    stat.summary<-round(rbind(avr, se),4)

    percent.avr.at.times[,,i]<-stat.summary

    percent.steer.at.times[,,i]<-cbind(steer.at.second.pt,  steer.at.third.pt, steer.at.fourth.pt,  steer.at.fifth.pt,steer.at.sixth.pt)

  }


  ###########################################################################
  ### T-test if slope of each level of time factor is increasing or not    ##
  ###########################################################################

  summary.ttest<-array(0,dim=c(5, 4, num.aa), dimnames=list( trim.name.tt,
                                                             c( "outcome.avr", "T-value", "num.para","p.value"), "true.AA"))


  for (i in 1:num.aa){
    for (t in 1:(time.points-1)){
      test.time.eff<-t.test(percent.steer.at.times[,t,i], mu=0, alternative="two.sided", conf.level=0.95)

      summary.ttest[t,,i]<-as.numeric(cbind(test.time.eff$estimate, test.time.eff$statistic,
                                            test.time.eff$parameter, test.time.eff$p.value))

    }
  }

  summary.ttest<-round(summary.ttest,3)

  ###}}






  ########################################################################
  ## Produce plot (a) Individiual trajectories by original scaled data  ##
  ########################################################################


  ##{{

  for (aa in 1:num.aa){

    substeer<-data.steer

    ### sort tag numbers for distinct steers ###
    steer.tag<-unique(substeer[,"steer.id"])

    ### choose unique steer ###
    one.uniq.steer<-which(substeer[,"steer.id"]==steer.tag[1])
    one.uniq.steer.data<-substeer[one.uniq.steer,]

    ### obtain amino acids name ###
    ## amino.names: "CIT"

    min.amino<-min(substeer[,amino.names])
    max.amino<-max(substeer[,amino.names])
    ylim<-round(c(min.amino-5, max.amino+20), 0)


    ##############
    ### colors ###
    ##############

    cols<-c("black", "blue", "red", "seagreen", "purple", "orange")

    #########################################################################
    ## Create plot for the concentrations of amino acids for the first cow ##
    #########################################################################

    ###################################
    ### plot names for saving file  ###
    ###################################

    file0<-paste("", "true.AA", sep="")
    postscript(paste(file0,"_origin_scale.eps",sep=""))

    ##################################
    ## obtain names of amino acids ###
    ##################################

    yvalues<-one.uniq.steer.data[,amino.names]
    ylab.name<-paste("concentrations of ", "true.AA",sep="")


    ## plot (a)

    plot(one.uniq.steer.data[,"time"], yvalues,
         ylim=ylim, type="o", ylab=ylab.name, xlab="time",
         col=cols[1], cex.lab=1.5, lwd=2.5, axes=F)

    for (id in 2:n){

      uniq.steer<-which(substeer[,"steer.id"]==steer.tag[id])
      uniq.steer.data<-substeer[uniq.steer,]

      lines(uniq.steer.data[,"time"], uniq.steer.data[,amino.names], type="o", lty=id, col=cols[id], lwd=2.5)

    }

    ## x-axis and y-axis for plot (a)

    axis(1,at=c(one.uniq.steer.data[,"time"]),labels=c(one.uniq.steer.data[,"time"]), lwd=2.5)
    interval<-(ylim[2]-ylim[1])/interv.length
    ylabs<-round(seq(ylim[1],ylim[2], by=interval),0)
    axis(2,at=c(ylabs),labels=c(ylabs), lwd=2.5)

    ## title for plot (a)

    mtext('(a) Original Scaled Indvidual Trajectories', outer=F, side=3, line=1,cex=1.5)


    ## legend for plot (a)

    legend(0,max(ylabs),paste("Steer ", steer.tag, sep="")[1:3], col=cols[1:3],
           lty=c(1:3), bty="n", cex=1.3, lwd=2)

    legend(2,max(ylabs),paste("Steer ", steer.tag, sep="")[4:6], col=cols[4:6],
           lty=c(4:6), bty="n", cex=1.3, lwd=2)


    dev.off()

  }

  ###end plot (a) }}






  ########################################################################
  ## Produce plot (b) Individiual trajectories by normalized data       ##
  ## using percentage-changes from baselines                            ##
  ########################################################################

  ##{{

  full.trans.steer<-do.call("rbind.data.frame", uniq.steer.list)


  for (aa in 1:num.aa){

    j<-3+aa

    ## remmove original scaled Citrulline concentraions with dummy variables for tag and time
    trans.substeer<- full.trans.steer[,-((3+num.aa):(5+num.aa))]

    ### sort tag numbers for distinct steers ###
    steer.tag<-unique(trans.substeer[,"steer.id"])

    ### choose unique steer ###
    one.uniq.steer<-which( trans.substeer[,"steer.id"]==steer.tag[1])
    one.uniq.steer.data<- trans.substeer[one.uniq.steer,]

    ### obtain amino acids name ###
    amino.name<-as.character(names(one.uniq.steer.data)[4])

    min.amino<-min(trans.substeer[,amino.name])
    max.amino<-max(trans.substeer[,amino.name])
    ylim<-round(c(min.amino-5, max.amino+20), 0)


    ##############
    ### colors ###
    ##############
    # already defined colors previously
    #cols<-c("black", "blue", "red", "seagreen", "purple", "orange")

    #########################################################################
    ## Create plot for the concentrations of amino acids for the first cow ##
    #########################################################################

    ###################################
    ### plot names for saving file  ###
    ###################################


    file1<-paste("", "true.AA", sep="")
    postscript(paste(file1,"_percent_change.eps",sep=""))

    ##################################
    ## obtain names of amino acids ###
    ##################################

    norm.yvalues<-one.uniq.steer.data[,norm.amino.names]
    norm.ylab.name<-paste("concentrations of norm.", "true.AA",sep="")


    ## plot

    plot(one.uniq.steer.data[,"tt"], norm.yvalues,
         ylim=ylim, type="o", ylab=norm.ylab.name, xlab="time",
         col=cols[1], cex.lab=1.5, lwd=2.5, axes=F)

    for (id in 2:n){

      uniq.steer<-which(trans.substeer[,"steer.id"]==steer.tag[id])
      uniq.steer.data<-trans.substeer[uniq.steer,]

      lines(uniq.steer.data[,"tt"], uniq.steer.data[,norm.amino.names], type="o",
      lty=id, col=cols[id], lwd=2.5)

    }


    ## x-axis and y-axis for plot (b)

    axis(1,at=c(one.uniq.steer.data[,"tt"]),labels=c(one.uniq.steer.data[,"tt"]), lwd=2.5)
    interval<-(ylim[2]-ylim[1])/interv.length
    ylabs<-round(seq(ylim[1],ylim[2], by=interval),0)
    axis(2,at=c(ylabs),labels=c(ylabs), lwd=2.5)


    ## title for plot (b)
    mtext('(b) Individual Trajectories after Percentage-changes', outer=F, side=3, line=1,cex=1.5)



    ## legend for plot (b)
    legend(0,max(ylabs),paste("Steer ", steer.tag, sep="")[1:3], col=cols[1:3],
           lty=c(1:3), bty="n", cex=1.3, lwd=2)

    legend(2,max(ylabs),paste("Steer ", steer.tag, sep="")[4:6], col=cols[4:6],
           lty=c(4:6), bty="n", cex=1.3, lwd=2)



    dev.off()
  }


  ### end plot (b)}}






  ########################################################################
  ## Produce plot (c) average trajectory of concentrations of cit       ##
  ## for (1) original data                                              ##
  ##     (2) noramlized data using percentage-changes from baselines    ##
  ########################################################################




  #######################################################################################
  ## Initialization to have data summary: average and standard errors for each steer   ##
  #######################################################################################


  origial.scaled.summarystat<-array(0,dim=c(2, 6, num.aa), dimnames=list( c( "avr", "sterr"),
                                                                         name.tt, paste("orig.scaled.", "true.AA", sep="")))



  #######################################################################################################
  ##For original scaled data                                                                           ##
  ##Obtain average concentration of citrulline at each time point for all steers with standard errors  ##
  #######################################################################################################

  for (i in 1:num.aa){

    aa<-i+3
    uniq.tt<-unique(tt)

    init.ind<-which(data.steer[,"time"]==uniq.tt[1])
    avr.at.zeroh<-mean(data.steer[,aa][init.ind])
    sd.at.zeroh<-sd(data.steer[,aa][init.ind])/sqrt(n)

    second.pt<-which(data.steer[,"time"]==uniq.tt[2])
    avr.at.second.pt<-mean(data.steer[,aa][second.pt])
    sd.at.second.pt<-sd(data.steer[,aa][second.pt])/sqrt(n)

    
    third.pt<-which(data.steer[,"time"]==uniq.tt[3])
    avr.at.third.pt<-mean(data.steer[,aa][third.pt])
    sd.at.third.pt<-sd(data.steer[,aa][third.pt])/sqrt(n)
    
    
    fourth.pt<-which(data.steer[,"time"]==uniq.tt[4])
    avr.at.fourth.pt<-mean(data.steer[,aa][fourth.pt])
    sd.at.fourth.pt<-sd(data.steer[,aa][fourth.pt])/sqrt(n)
    
    
    fifth.pt<-which(data.steer[,"time"]==uniq.tt[5])
    avr.at.fifth.pt<-mean(data.steer[,aa][fifth.pt])
    sd.at.fifth.pt<-sd(data.steer[,aa][fifth.pt])/sqrt(n)
    

    
    sixth.pt<-which(data.steer[,"time"]==uniq.tt[6])
    avr.at.sixth.pt<-mean(data.steer[,aa][sixth.pt])
    sd.at.sixth.pt<-sd(data.steer[,aa][sixth.pt])/sqrt(n)
    

    avr<-cbind(avr.at.zeroh, avr.at.second.pt, avr.at.third.pt, avr.at.fourth.pt,  avr.at.fifth.pt, avr.at.sixth.pt
               )
    se<-cbind(sd.at.zeroh, sd.at.second.pt, sd.at.third.pt, sd.at.fourth.pt,  sd.at.fifth.pt, sd.at.sixth.pt)


    stat.summary<-round(rbind(avr, se),2)
    origial.scaled.summarystat[,,i]<-stat.summary


  }



  ###################################
  ## plot (c) for the original scaled
  ###################################

  ##{{

  file2<-paste("",  "true.AA", sep="")
  postscript(paste(file2,"_mean_concentrations.eps",sep=""))


  ### x and y values ###

  time.vec<-c(unique(data.steer[,"time"]))
  avr.raw.yvalues<-origial.scaled.summarystat[,,1]["avr", ]

  ## add 0 at time0 for the percentage-changes from baselines
  avr.norm.yvalues<-c(0, percent.avr.at.times[,,1]["avr", ])

  ### range of y values and name of y axis ###

  ylim<-c(-5, 160)
  ylab.name<-paste("mean concentrations of ", "true.AA",sep="")

  ### create color vector ####
  #cols<-c("black", "blue", "red", "seagreen", "purple", "orange")
  #colors are defined previously

  plot(time.vec, avr.raw.yvalues,
       ylim=ylim, type="o", ylab=ylab.name, xlab="time",
       col=cols[1], cex.lab=1.5, lwd=2.5, lty=1, axes=F)

  lines(time.vec,  avr.norm.yvalues, type="o",
        col=cols[2], cex.lab=1.5, lwd=2.5, lty=2)


  mtext('(c) Average Trajectory', outer=F, side=3, line=1,cex=1.5)

  legend(0,160, c("Original Scaled Data", "Percentage-chages from Baselines"), col=cols[1:2],
         lty=c(1, 2), bty="n", cex=1.3, lwd=2)


  axis(1,at=c(time.vec),labels=c(time.vec), lwd=2.5)
  interval<-(ylim[2]-ylim[1])/interv.length
  ylabs<-round(seq(ylim[1],ylim[2], by=interval),0)
  axis(2,at=c(ylabs),labels=c(ylabs), lwd=2.5)

  dev.off()

  ### end plot (c)}}


  ## returns list of results: table 1 and table 2 in the main paper
  ## returns plots (a), (b) and (c) in Figure 1 in your directory

  return(list(summary.ttest=summary.ttest))
}




##############################
##############################
# #Run function
# library(nlme)
# anova.test(file="rpaa_cit_data.csv")
# ###############################
