#' ysim.ppc.plot.code - Create plot for Prior and Poserior Predictive Checks
#'
#' @param dat  Data formatted for Blouch models by blouch.prep functions
#' @param post Posterior distribution of stanfit class
#' @param row.nums Rows to be sampled from distribution
#'
#' @return Plots in ggplot2 format
#' @export
#'
ysim.ppc.plot.code<-function(dat,post,row.nums){
  df<-data.frame(cbind(Y_obs=dat$Y_obs,t(data.frame(post$Y_sim)[row.nums*100,])))
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  #return(df)
  plots<-list()
  for (i in 1:max(row.nums)){
    name.column<-paste("X",i*100,sep="")
    prior.plot<-ggplot2::ggplot(data=df)+
      ggplot2::geom_density(ggplot2::aes(.data[[name.column]],fill="Simulated Y"),alpha=0.2)+
      ggplot2::geom_density(ggplot2::aes(.data$Y_obs,fill="Observed Y"),alpha=0.2)+
      ggplot2::theme_bw()+
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())+
      ggplot2::labs(title="",x="Y", y = "Density")+
      #ggsci::scale_fill_npg(name="",labels=c("Simulated Y","Observed Y"))
      ggplot2::scale_fill_manual(values=rev(mypal),name="",labels=c("Simulated Y","Observed Y"))

    prior.plot<-prior.plot + ggplot2::theme(legend.position = "none")
    plots[[i]]<-prior.plot
  }
  #modified_plots <- lapply(plots, function(plot) {
  #  plot +  ggsci::scale_fill_npg(name="",labels=c("Observed Y","Simulated Y"))
  #})
  return(plots)
}
#' hl.prior.plot.code- Create Plot for Half-life Prior Distribution
#'
#' @param hl.prior  Prior log mean and log sd for half-life
#'
#' @return Plots in ggplot2 format
#' @export
#'
hl.prior.plot.code<-function(hl.prior){
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  names(hl.sims)<-"prior.hl.sims"
  hl.prior.plot<-ggplot2::ggplot()+
    ggplot2::geom_density(ggplot2::aes(prior.hl.sims,fill="prior.hl.sims"),alpha=0.2,data=hl.sims)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Half-life", y = "Density")+
    #ggplot2::geom_vline(xintercept=c(hl),linetype=2)+
    ggsci::scale_fill_npg(name="",labels=c("Prior"))
  return(hl.prior.plot)
}
#' hl.prior.plot.code- Create Prior vs Posterior Plot for Half-life Prior Distribution for Simulated Dataset
#'
#' @param hl.prior  Prior log mean and log sd for half-life
#' @param post  Posterior distribution of stanfit class
#' @param hl True half-life value
#' @return Plots in ggplot2 format
#' @export
#'
hl.prior.post.plot.code<-function(hl.prior,post,hl){
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  mypal[2]<-palette()[1]

  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  names(hl.sims)<-"prior.hl.sims"
  hl.post<-data.frame(post$hl)
  names(hl.post)<-"post.hl.sims"
  df<-data.frame(cbind(hl.sims,hl.post))

  hl.prior.plot<-ggplot2::ggplot(data=df)+
    ggplot2::geom_density(ggplot2::aes(prior.hl.sims, fill=mypal[2]),alpha=0.2)+
    ggplot2::geom_density(ggplot2::aes(post.hl.sims, fill=mypal[1]),alpha=0.2)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Half-life", y = "Density")+
    ggplot2::geom_vline(xintercept=c(hl),linetype=2)+
    ggplot2::scale_fill_manual(values=mypal,name="",labels=c("Posterior","Prior"))
  return(hl.prior.plot)
}
#' hl.prior.post.emp.plot.code- Create Prior vs. Posterior Plot for Half-life Prior Distribution for Empirical Dataset
#'
#' @param hl.prior  Prior log mean and log sd for half-life
#' @param post  Posterior distribution of stanfit class
#' @return Plots in ggplot2 format
#' @export
#
hl.prior.post.emp.plot.code<-function(hl.prior,post){
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  mypal[2]<-palette()[1]
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  names(hl.sims)<-"prior.hl.sims"
  hl.post<-data.frame(post$hl)
  names(hl.post)<-"post.hl.sims"
  df<-data.frame(cbind(hl.sims,hl.post))

  hl.prior.plot<-ggplot2::ggplot(data=df)+
    ggplot2::geom_density(ggplot2::aes(prior.hl.sims, fill=mypal[2]),alpha=0.2)+
    ggplot2::geom_density(ggplot2::aes(post.hl.sims, fill=mypal[1]),alpha=0.2)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Half-life", y = "Density")+
    ggplot2::scale_fill_manual(values=mypal,name="",labels=c("Posterior","Prior"))
  return(hl.prior.plot)
}

#' vy.prior.plot.code- Create plot for Vy Prior Distribution
#'
#' @param vy.prior  Prior log mean and log sd for half-life
#' @return Plots in ggplot2 format
#' @export
#'
vy.prior.plot.code<-function(vy.prior){
  vy.sims<-rexp(n=1000,rate=vy.prior)
  vy.sims<-data.frame(vy.sims)
  names(vy.sims)<-"prior.vy.sims"
  vy.prior.plot<-ggplot2::ggplot()+
    ggplot2::geom_density(ggplot2::aes(prior.vy.sims,fill="prior.vy.sims"),alpha=0.2,data=vy.sims)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Vy", y = "Density")+
    #ggplot2::geom_vline(xintercept=log(vy),linetype=2)+
    ggsci::scale_fill_npg(name="",labels=c("Prior"))
  return(vy.prior.plot)
}
#' vy.prior.post.plot.code- Create Prior vs. Posterior Plot for Vy
#'
#' @param vy.prior  Prior scale parameter
#' @param post  Posterior distribution of stanfit class
#' @param vy True Vy
#' @return Plots in ggplot2 format
#' @export
#'
vy.prior.post.plot.code<-function(vy.prior,post,vy){
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  mypal[2]<-palette()[1]
  vy.sims<-rexp(n=1000,rate=vy.prior)
  vy.sims<-data.frame(vy.sims)
  names(vy.sims)<-"prior.vy.sims"
  vy.post<-data.frame(post$vy)
  names(vy.post)<-"post.vy.sims"
  df<-data.frame(cbind(vy.sims,vy.post))
  vy.prior.post.plot<-ggplot2::ggplot(data=df)+
    ggplot2::geom_density(ggplot2::aes(prior.vy.sims,fill="prior.vy.sims"),alpha=0.2,data=vy.sims)+
    ggplot2::geom_density(ggplot2::aes(post.vy.sims,fill="post.vy.sims"),alpha=0.2,data=vy.post)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Vy", y = "Density")+
    ggplot2::geom_vline(xintercept=c(vy),linetype=2)+
    ggplot2::scale_fill_manual(values=mypal,name="",labels=c("Posterior","Prior"))
  return(vy.prior.post.plot)
}
#' vy.prior.post.emp.plot.code- Create Prior vs. Posterior Plot for Vy for Empirical Dataset
#'
#' @param vy.prior  Prior log mean and log sd for half-life
#' @param post  Posterior distribution of stanfit class
#' @return Plots in ggplot2 format
#' @export
#'
vy.prior.post.emp.plot.code<-function(vy.prior,post){
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  mypal[2]<-palette()[1]
  vy.sims<-rexp(n=1000,rate=vy.prior)
  vy.sims<-data.frame(vy.sims)
  names(vy.sims)<-"prior.vy.sims"
  vy.post<-data.frame(post$vy)
  names(vy.post)<-"post.vy.sims"
  df<-data.frame(cbind(vy.sims,vy.post))
  vy.plot<-ggplot2::ggplot(data=df)+
    ggplot2::geom_density(ggplot2::aes(prior.vy.sims, fill=mypal[2]),alpha=0.2)+
    ggplot2::geom_density(ggplot2::aes(post.vy.sims, fill=mypal[1]),alpha=0.2)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Vy", y = "Density")+
    ggplot2::scale_fill_manual(values=mypal,name="",labels=c("Posterior","Prior"))
  return(vy.prior.post.plot)
}

#' sigma.prior.plot.code- Create Prior vs. Posterior Plot for Sigma
#'
#' @param sigma.prior  Prior for mean and sd parameters of sigma
#' @return Plots in ggplot2 format
#' @export
#'
sig.prior.plot.code<-function(sigma.prior){
  sigma.sims<-data.frame(abs(rnorm(n=1000,mean=sigma.prior[1],sd=sigma.prior[2])))
  names(sigma.sims)<-"prior.sigma.sims"
  sigma.prior.plot<-ggplot2::ggplot()+
    ggplot2::geom_density(ggplot2::aes(prior.sigma.sims,fill="prior.sigma.sims"),alpha=0.2,data=sigma.sims)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Sigma", y = "Density")+
    ggsci::scale_fill_npg(name="",labels=c("Prior"))
  return(sigma.prior.plot)
}
#' sigma.prior.post.plot.code- Create Prior vs. Posterior Plot for Sigma
#'
#' @param sigma.prior  Prior mean and standard deviation parmateres
#' @param post  Posterior distribution of stanfit class
#' @return Plots in ggplot2 format
#' @export
#'
sig.prior.post.plot.code<-function(sigma.prior,post){
  sigma.sims<-data.frame(abs(rnorm(n=1000,mean=sigma.prior[1],sd=sigma.prior[2])))
  names(sigma.sims)<-"prior.sigma.sims"
  sigma.post<-data.frame(post$sigma)
  names(sigma.post)<-c("post.sigma.optima","post.sigma.beta")

  sigma.plot<-ggplot2::ggplot()+
    ggplot2::geom_density(ggplot2::aes(prior.sigma.sims,fill="prior.sigma.sims"),alpha=0.2,data=sigma.sims)+
    ggplot2::geom_density(ggplot2::aes(post.sigma.optima,fill="post.sigma.optima"),alpha=0.2,data=sigma.post)+
    ggplot2::geom_density(ggplot2::aes(post.sigma.beta,fill="post.sigma.beta"),alpha=0.2,data=sigma.post)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Sigma", y = "Density")+
    ggsci::scale_fill_npg(name="",labels=c("Posterior Beta","Posterior Optima","Prior"))
  return(sigma.plot)
}
#' covariance.prior.direct.plot.code- Create Prior Covariance Plot for Direct Effect Model
#'
#' @param hl.prior  Prior mean and scale parameter for Half-life
#' @param vy.prior Prior scale parameter for Vy
#' @return Plots in ggplot2 format
#' @export
#'
covariance.prior.direct.plot.code<-function(hl.prior,vy.prior){
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  vy.sims<-rexp(n=1000,rate=vy.prior)
  a.sims<-log(2)/hl.sims
  sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims))
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  x<-seq(0,1,by=0.01)
  df<-data.frame(x)

  p<-ggplot2::ggplot(df,ggplot2::aes(x))
  for(i in 1:30){
    p<-p+ggplot2::stat_function(fun=function(x,i) {sigma2_y.sims[i,] /(2 * a.sims[i,]) * ((1 - exp(-2 * a.sims[i,] * (1-(x/2)))) * exp(-a.sims[i,] * x))},
                                args=list(i=i),alpha=0.2,lwd=2)}
  p<-p+ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Time Since MRCA", y = "Covariance")
  return(p)
}
#' covariance.prior.post.direct.plot.code- Create Prior vs. Posterior Covariance Plot for Direct Effect Model
#'
#' @param hl.prior  Prior mean and scale parameter for Half-life
#' @param vy.prior Prior scale parameter for Vy
#' @param post Posterior distribution in stanfit class
#' @return Plots in ggplot2 format
#' @export
#'
covariance.prior.post.direct.plot.code<-function(hl.prior,vy.prior,post){
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  vy.sims<-rexp(n=1000,rate=vy.prior)
  a.sims<-log(2)/hl.sims
  sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims))
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  x<-seq(0,1,by=0.01)
  df<-data.frame(x)
  p<-ggplot2::ggplot(df,ggplot2::aes(x))
  for(i in 1:30){
    p<-p+ggplot2::stat_function(fun=function(x,i) {sigma2_y.sims[i,] /(2 * a.sims[i,]) * ((1 - exp(-2 * a.sims[i,] * (1-(x/2)))) * exp(-a.sims[i,] * x))},
                                args=list(i=i),alpha=0.2,lwd=2)
    p<-p+ggplot2::stat_function(fun=function(x,i) {post$sigma2_y[i] /(2 * post$a[i]) * ((1 - exp(-2 * post$a[i] * (1-(x/2)))) * exp(-post$a[i] * x))},
                                args=list(i=i),alpha=0.2,lwd=2,color=mypal[1])
  }
  p<-p+ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Time Since MRCA", y = "Covariance")+
    ggplot2::scale_fill_manual(values=mypal,name="",labels=c("Prior","Posterior"))

  return(p)
}
#' covariance.prior.adapt.plot.code- Create Prior Covariance Plot for Adaptive Model
#'
#' @param hl.prior  Prior mean and scale parameter for Half-life
#' @param vy.prior Prior scale parameter for Vy
#' @param beta.prior Prior mean and scale parameter for beta/slope
#' @return Plots in ggplot2 format
#' @export
#'
covariance.prior.adapt.plot.code<-function(hl.prior,vy.prior,beta.prior){
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  vy.sims<-rexp(n=1000,rate=vy.prior)
  beta.sims<-data.frame(rnorm(n=1000,beta.prior[1],beta.prior[2]))
  a.sims<-log(2)/hl.sims
  sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims))
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  x<-seq(0,1,by=0.01)
  df<-data.frame(x)

  p<-ggplot2::ggplot(df,ggplot2::aes(x))
  for(i in 1:30){
    p<-p+ggplot2::stat_function(fun=function(x,i){calc_adaptive_cov_plot(a.sims[i,],sigma2_y.sims[i,],beta.sims[i,],x)},
                                args=list(i=i),alpha=0.2,lwd=2)
  }

  p<-p+ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Time Since MRCA", y = "Covariance")
  return(p)
}
#' covariance.prior.post.adapt.plot.code- Create Prior vs. Posterior Covariance Plot for Adaptive Model
#'
#' @param hl.prior  Prior mean and scale parameter for Half-life
#' @param vy.prior Prior scale parameter for Vy
#' @param beta.prior Prior mean and scale parameter for beta/slope
#' @param post Posterior distribution in stanfit class
#' @return Plots in ggplot2 format
#' @export
#'
covariance.prior.post.adapt.plot.code<-function(hl.prior,vy.prior,beta.prior,post){
  hl.sims<-data.frame(rlnorm(n=1000,meanlog=hl.prior[1],sdlog=hl.prior[2]))
  vy.sims<-rexp(n=1000,rate=vy.prior)
  beta.sims<-data.frame(rnorm(n=1000,beta.prior[1],beta.prior[2]))
  a.sims<-log(2)/hl.sims
  sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims))
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  x<-seq(0,1,by=0.01)
  df<-data.frame(x)

  p<-ggplot2::ggplot(df,ggplot2::aes(x))
  for(i in 1:30){
    p<-p+ggplot2::stat_function(fun=function(x,i){calc_adaptive_cov_plot(a.sims[i,],sigma2_y.sims[i,],beta.sims[i,],x)},
                                args=list(i=i),alpha=0.2,lwd=2)
    p<-p+ggplot2::stat_function(fun=function(x,i){calc_adaptive_cov_plot(post$a[i],post$sigma2_y[i],post$beta[i],x)},
                                args=list(i=i),alpha=0.2,color=mypal[1],lwd=2)


  }

  p<-p+ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(title="",x="Time Since MRCA", y = "Covariance")
  ggplot2::scale_color_manual(values=mypal,name="",labels=c("Posterior","Prior"))
  return(p)
}
############################################################################################################
#' direct.prior.plot.code- Create Prior Plot for Direct Effect Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
direct.prior.plot.code<-function(trdata,optima.prior,beta.prior){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  prior.slope.plot<-ggplot2::ggplot()+
    ggplot2:: geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=X_with_error))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  return(prior.slope.plot)
}
############################################################################################################
#' direct.prior.plot.code- Create Prior vs. Posterior Plot for Direct Effect Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
direct.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  optima.post<-data.frame(post$optima)
  names(optima.post)<-"post.optima"

  beta.post<-data.frame(post$beta)
  names(beta.post)<-"post.beta"

  mu.link<-function(x.seq){optima.post+x.seq*beta.post}
  x.seq <- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)
  mu <- sapply(x.seq , mu.link )
  mu.mean <-lapply( mu , mean )
  mu.mean<-data.frame(as.numeric(mu.mean))
  names(mu.mean)<-"mu.mean"

  mu.CI <- lapply( mu , rethinking::PI , prob=0.89 )
  mu.CI<-data.frame(t(data.frame(mu.CI)),x.seq)
  names(mu.CI)<-c("min.5.5","max.94.5","x.seq")

  df2<-data.frame(x.seq,mu.mean)

  prior.post.slope.plot<-ggplot2::ggplot()+
    ggplot2::geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=X_with_error))+
    ggplot2::geom_abline(intercept=optima,slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::geom_line(data=df2,ggplot2::aes(x=x.seq,y=mu.mean),linetype=1)+

    ggplot2::geom_ribbon(data=mu.CI,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  return(prior.post.slope.plot)
}
#' adapt.prior.plot.code- Create Prior Plot for Adaptive Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
adapt.prior.plot.code<-function(trdata,optima.prior,beta.prior){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  prior.slope.plot<-ggplot2::ggplot()+
    ggplot2:: geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=X_with_error))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("Adaptative trait")+
    ggsci::scale_color_npg()

  return(prior.slope.plot)
}
#' adapt.prior.post.plot.code- Create Prior vs. Posterior Plot for Adaptive Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
adapt.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)

  optima.post<-data.frame(post$optima)
  names(optima.post)<-"post.optima"

  beta.post<-data.frame(post$beta)
  names(beta.post)<-"post.beta"

  mu.link<-function(x.seq){optima.post+x.seq*beta.post}
  x.seq <- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)
  mu <- sapply(x.seq , mu.link )
  mu.mean <-lapply( mu , mean )
  mu.mean<-data.frame(as.numeric(mu.mean))
  names(mu.mean)<-"mu.mean"

  mu.CI <- lapply( mu , rethinking::PI , prob=0.89 )
  mu.CI<-data.frame(t(data.frame(mu.CI)),x.seq)
  names(mu.CI)<-c("min.5.5","max.94.5","x.seq")

  df2<-data.frame(x.seq,mu.mean)
  prior.post.slope.plot<-ggplot2::ggplot()+

    ggplot2::geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=X_with_error))+
    ggplot2::geom_abline(intercept=optima,slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::geom_line(data=df2,ggplot2::aes(x=x.seq,y=mu.mean),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("Adaptative trait")+
    ggsci::scale_color_npg()

  return(prior.post.slope.plot)
}
#' optima.prior.plot.code- Create Prior Plot for Multi-Optima Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @return Plots in ggplot2 format
#' @export
#'
optima.prior.plot.code<-function(trdata,optima.prior){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  optima.sims<-data.frame(optima.sims,Prior="Prior")
  df<-data.frame(Y=trdata$dat$Y_with_error,regimes=trdata$dat$regimes)
  names(df)<-c("Y","Regimes")
  optima.prior.plot<-ggplot2::ggplot(data=df, Mapping = aes(x = Y, y = Regimes))+
    ggplot2::geom_jitter(data=optima.sims,ggplot2::aes(y=optima.sims,x=Prior),alpha=0.25,width=0.15)+
    ggplot2::geom_jitter(data=df,ggplot2::aes(y=Y,x=Regimes),width=0.15)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    ggplot2::ylab("Optima") + ggplot2::xlab("")+
    ggsci::scale_color_npg()

  return(optima.prior.plot)
}
#' optima.prior.post.plot.code- Create Prior Plot for Multi-Optima Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @return Plots in ggplot2 format
#' @export
#'
optima.prior.post.plot.code<-function(trdata,optima.prior,post,optima){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  optima.sims<-data.frame(optima.sims,Prior="Prior")
  optima.post<-data.frame(post$optima)
  names(optima.post)<-c(paste("OU",1:dim(optima.post)[2],sep=""))

  mu.mean <-apply( optima.post, 2, mean )
  mu.mean<-data.frame(mu.mean)
  mu.mean<-data.frame(mu.mean,"Regimes"=1:dim(mu.mean)[1])
  mu.mean<-data.frame(mu.mean,optima)

  mu.CI <- apply( optima.post , 2, rethinking::PI , prob=0.89 )
  mu.CI<-data.frame(t(data.frame(mu.CI)),"Regimes"=1:dim(mu.CI)[2])
  names(mu.CI)<-c("min.5.5","max.94.5","Regimes")

  df<-data.frame(Y=trdata$dat$Y_with_error,regimes=trdata$dat$regimes)
  names(df)<-c("Y","Regimes")


  optima.prior.post.plot<-
    ggplot2::ggplot(data=df, Mapping = ggplot2aes(x = Y, y = Regimes))+

    ggplot2::geom_jitter(data=optima.sims,ggplot2::aes(y=optima.sims,x=Prior),alpha=0.25,width=0.15)+
    ggplot2::geom_segment(data=mu.mean,ggplot2::aes(x=Regimes-0.5,xend=Regimes+0.5,y=mu.mean,yend=mu.mean),linetype=1)+
    ggplot2::geom_segment(data=mu.mean,ggplot2::aes(x=Regimes-0.5,xend=Regimes+0.5,y=optima,yend=optima),linetype=2)+
    ggplot2::geom_linerange(data=mu.CI, mapping=ggplot2::aes(x=Regimes,ymin=min.5.5,ymax=max.94.5),linewidth=35,color="darkgrey",alpha=0.5)+
    ggplot2::geom_jitter(data=df,ggplot2::aes(y=Y,x=Regimes),width=0.15)+

    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    ggplot2::ylab("Optima") + ggplot2::xlab("")+
    ggsci::scale_color_npg()


  return(optima.prior.post.plot)
}
#' direct.adapt.prior.plot.code- Create Prior Plot for Direct Effect + Adaptive Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
direct.adapt.prior.plot.code<-function(trdata,optima.prior,beta.prior){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  prior.slope.plot.1<-ggplot2::ggplot()+
    ggplot2:: geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=Xd))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("Direct Effect Model - X")+
    ggsci::scale_color_npg()
  prior.slope.plot.2<-ggplot2::ggplot()+
    ggplot2:: geom_point(data=as.data.frame(trdata$dat),ggplot2::aes(y=Y_with_error,x=Xa))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("Adptation Model - X")+
    ggsci::scale_color_npg()
  return(list(prior.slope.plot.1,prior.slope.plot.2))
}
#' direct.adapt.prior.post.plot.code- Create Prior vs. Posterior Plot for Direct Effect + Adaptive Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
direct.adapt.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  mypal <- ggsci::pal_npg("nrc", alpha = 0.4)(2)
  optima.post<-data.frame(post$optima)
  names(optima.post)<-"post.optima"
  beta.post.1<-data.frame(post$beta[,1])
  names(beta.post.1)<-"post.beta.1"
  beta.post.2<-data.frame(post$beta[,2])
  names(beta.post.2)<-"post.beta.2"
  mu.link.1<-function(x.seq){optima.post+x.seq*beta.post.1}
  mu.link.2<-function(x.seq){optima.post+x.seq*beta.post.2}
  x.seq.d <- seq(from=min(trdata$dat$Xd), to=max(trdata$dat$Xd) , length.out=100)
  x.seq.a <- seq(from=min(trdata$dat$Xa), to=max(trdata$dat$Xa) , length.out=100)
  mu.1 <- sapply(x.seq.d , mu.link.1 )
  mu.2 <- sapply(x.seq.a , mu.link.2 )
  mu.mean.1 <-lapply( mu.1 , mean )
  mu.mean.2 <-lapply( mu.2 , mean )
  mu.mean.1<-data.frame(as.numeric(mu.mean.1))
  mu.mean.2<-data.frame(as.numeric(mu.mean.2))
  names(mu.mean.1)<-"mu.mean.1"
  names(mu.mean.2)<-"mu.mean.2"

  mu.CI.1 <- lapply( mu.1 , rethinking::PI , prob=0.89 )
  mu.CI.2 <- lapply( mu.2 , rethinking::PI , prob=0.89 )
  mu.CI.1<-data.frame(t(data.frame(mu.CI.1)),x.seq.d)
  mu.CI.2<-data.frame(t(data.frame(mu.CI.2)),x.seq.a)

  names(mu.CI.1)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.2)<-c("min.5.5","max.94.5","x.seq")

  df<-data.frame(Y=trdata$dat$Y_with_error,Xd=trdata$dat$Xd,Xa=trdata$dat$Xa)
  df2<-data.frame(x.seq.d,mu.mean.1)
  df3<-data.frame(x.seq.a,mu.mean.2)

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xd))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::geom_abline(intercept=optima,slope=beta[1],alpha=0.5,linetype=2)+
    ggplot2::geom_line(data=df2,ggplot2::aes(x=x.seq.d,y=mu.mean.1),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.1,ggplot2::aes(x=x.seq.d,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1

  slope.plot.2<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xa))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::geom_abline(intercept=optima,slope=beta[2],alpha=0.5,linetype=2)+
    ggplot2::geom_line(data=df3,ggplot2::aes(x=x.seq.a,y=mu.mean.2),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.2,ggplot2::aes(x=x.seq.a,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.2
  return(list(slope.plot.1,slope.plot.2))
}
#' reg.direct.prior.plot.code- Create Prior Plot for Multi-Optima Direct Effect Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.prior.plot.code<-function(trdata,optima.prior,beta.prior){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error)

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)

}
#' reg.adapt.prior.plot.code- Create Prior Plot for Multi-Optima Adaptive Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
reg.adapt.prior.plot.code<-function(trdata,optima.prior,beta.prior){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error)

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)


}
#' reg.direct.adapt.prior.plot.code- Create Prior Plot for Multi-Optima Direct Effect + Adaptive Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.adapt.prior.plot.code<-function(trdata,optima.prior,beta.prior){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  df<-data.frame(Y=trdata$dat$Y_with_error,Xd=trdata$dat$Xd,Xa=trdata$dat$Xa)

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xd))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1

  slope.plot.2<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xa))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.2

  return(list(slope.plot.1,slope.plot.2))

}
#' reg.direct.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima

  beta.post.1<-data.frame(post$beta[,1])
  names(beta.post.1)<-"post.beta.1"

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post.1}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post.1}
  x.seq <- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)

  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )
  mu.mean.11 <-lapply( mu.11 , mean )
  mu.mean.12 <-lapply( mu.12 , mean )

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- lapply( mu.11 , rethinking::PI , prob=0.89 )
  mu.CI.12 <- lapply( mu.12 , rethinking::PI , prob=0.89 )
  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)


  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)

  #slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_abline(intercept=optima[1],slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)

}
#' reg.direct.prior.post.emp.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect Model
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.prior.post.emp.plot.code<-function(trdata,optima.prior,beta.prior,post){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima

  beta.post.1<-data.frame(post$beta[,1])
  names(beta.post.1)<-"post.beta.1"

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post.1}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post.1}
  x.seq <- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)

  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )
  mu.mean.11 <-lapply( mu.11 , mean )
  mu.mean.12 <-lapply( mu.12 , mean )

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- lapply( mu.11 , rethinking::PI , prob=0.89 )
  mu.CI.12 <- lapply( mu.12 , rethinking::PI , prob=0.89 )
  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)


  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,Regimes=trdata$dat$BGS)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)

}
#' reg.adapt.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Adaptive Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.adapt.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima

  beta.post.1<-data.frame(post$beta[,1])
  names(beta.post.1)<-"post.beta.1"

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post.1}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post.1}
  x.seq <- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)

  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )
  mu.mean.11 <-lapply( mu.11 , mean )
  mu.mean.12 <-lapply( mu.12 , mean )

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- lapply( mu.11 , rethinking::PI , prob=0.89 )
  mu.CI.12 <- lapply( mu.12 , rethinking::PI , prob=0.89 )
  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)


  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)

  #slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_abline(intercept=optima[1],slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta,alpha=0.5,linetype=2)+
    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)

}
#' reg.direct.adapt.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect + Adaptive Model - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.adapt.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){
  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])
  optima.post<-post$optima

  beta.post.1<-data.frame(post$beta[,1])
  names(beta.post.1)<-"post.beta.1"

  beta.post.2<-data.frame(post$beta[,2])
  names(beta.post.2)<-"post.beta.2"

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post.1}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post.1}
  mu.link.21<-function(x.seq){optima.post[,1]+x.seq*beta.post.2}
  mu.link.22<-function(x.seq){optima.post[,2]+x.seq*beta.post.2}

  x.seq.Xd <- seq(from=min(trdata$dat$Xd), to=max(trdata$dat$Xd) , length.out=100)
  x.seq.Xa <- seq(from=min(trdata$dat$Xa), to=max(trdata$dat$Xa) , length.out=100)

  mu.11 <- sapply(x.seq.Xd , mu.link.11 )
  mu.12 <- sapply(x.seq.Xd , mu.link.12 )
  mu.21 <- sapply(x.seq.Xa , mu.link.21 )
  mu.22 <- sapply(x.seq.Xa , mu.link.22 )

  mu.mean.11 <-lapply( mu.11 , mean )
  mu.mean.12 <-lapply( mu.12 , mean )
  mu.mean.21 <-lapply( mu.21 , mean )
  mu.mean.22 <-lapply( mu.22 , mean )

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"
  mu.mean.21<-data.frame(as.numeric(mu.mean.21))
  mu.mean.22<-data.frame(as.numeric(mu.mean.22))
  names(mu.mean.21)<-"mu.mean.21"
  names(mu.mean.22)<-"mu.mean.22"

  mu.CI.11 <- lapply( mu.11 , rethinking::PI , prob=0.89 )
  mu.CI.12 <- lapply( mu.12 , rethinking::PI , prob=0.89 )
  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq.Xd)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq.Xd)

  mu.CI.21 <- lapply( mu.21 , rethinking::PI , prob=0.89 )
  mu.CI.22 <- lapply( mu.22 , rethinking::PI , prob=0.89 )
  mu.CI.21<-data.frame(t(data.frame(mu.CI.21)),x.seq.Xa)
  mu.CI.22<-data.frame(t(data.frame(mu.CI.22)),x.seq.Xa)

  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.21)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.22)<-c("min.5.5","max.94.5","x.seq")

  df<-data.frame(Y=trdata$dat$Y_with_error,Xd=trdata$dat$Xd,Xa=trdata$dat$Xa,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq.Xd,mu.mean.11)
  df12<-data.frame(x.seq.Xd,mu.mean.12)
  df21<-data.frame(x.seq.Xa,mu.mean.21)
  df22<-data.frame(x.seq.Xa,mu.mean.22)

  slope.plot.Xd<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xd,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_abline(intercept=optima[1],slope=beta[1],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[1],alpha=0.5,linetype=2)+

    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq.Xd,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq.Xd,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq.Xd,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq.Xd,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.Xd

  slope.plot.Xa<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xa,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+
    ggplot2::geom_abline(intercept=optima[1],slope=beta[2],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[2],alpha=0.5,linetype=2)+

    ggplot2::geom_line(data=df21,ggplot2::aes(x=x.seq.Xa,y=mu.mean.21),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.21,ggplot2::aes(x=x.seq.Xa,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df22,ggplot2::aes(x=x.seq.Xa,y=mu.mean.22),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.22,ggplot2::aes(x=x.seq.Xa,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.Xa

  return(list(slope.plot.Xd,slope.plot.Xa))

}
#' reg.direct.ve.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect Model - Varying Effects - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.ve.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  names(beta.post)<-c("post.beta.1","post.beta.2")

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}

  x.seq<- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)
  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )


  mu.mean.11<-colMeans(mu.11)
  mu.mean.12<-colMeans(mu.12)

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
  mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )

  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)

  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)


  #slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.1)+ #Prior

    ggplot2::geom_abline(intercept=optima[1],slope=beta[1],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[2],alpha=0.5,linetype=2)+

    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+

    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)
}
#' reg.direct.ve.prior.post.emp.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect Model - Varying Effects
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.ve.prior.post.emp.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){
  mypal <- ggsci::pal_aaas("default", alpha = 1)(3)

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  names(beta.post)<-c("post.beta.1","post.beta.2")

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}

  x.seq<- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)
  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )


  mu.mean.11<-colMeans(mu.11)
  mu.mean.12<-colMeans(mu.12)

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
  mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )

  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)

  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,BGS=trdata$dat$BGS)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)


  #slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=BGS))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.041)+ #Prior

    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+

    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("log Antler Volume (l)") + ggplot2::xlab("log Posterior Skull Length (cm)")+
    ggplot2::scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))
    #ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)
}
#' reg.adapt.ve.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Adaptive Model - Varying Effects - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.adapt.ve.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  names(beta.post)<-c("post.beta.1","post.beta.2")

  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}

  x.seq<- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)
  mu.11 <- sapply(x.seq , mu.link.11 )
  mu.12 <- sapply(x.seq , mu.link.12 )


  mu.mean.11<-colMeans(mu.11)
  mu.mean.12<-colMeans(mu.12)

  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
  mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )

  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)

  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq,mu.mean.11)
  df12<-data.frame(x.seq,mu.mean.12)


  #slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.1)+ #Prior

    ggplot2::geom_abline(intercept=optima[1],slope=beta[1],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[2],alpha=0.5,linetype=2)+

    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq,y=mu.mean.11),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+

    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  slope.plot.1
  return(slope.plot.1)
}
#' reg.direct.adapt.ve.prior.post.plot.code- Create Prior vs. Posterior Plot for Multi-Optima Direct Effect + Adaptive Model - Varying Effects - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param optima True optima
#' @param beta True beta
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.adapt.ve.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,optima,beta){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  names(beta.post)<-c("post.beta.1","post.beta.2","post.beta.3","post.beta.4")


  mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
  mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}

  mu.link.21<-function(x.seq){optima.post[,1]+x.seq*beta.post[,3]}
  mu.link.22<-function(x.seq){optima.post[,2]+x.seq*beta.post[,4]}

  x.seq.Xd <- seq(from=min(trdata$dat$Xd), to=max(trdata$dat$Xd) , length.out=100)
  x.seq.Xa <- seq(from=min(trdata$dat$Xa), to=max(trdata$dat$Xa) , length.out=100)

  mu.11 <- sapply(x.seq.Xd , mu.link.11 )
  mu.12 <- sapply(x.seq.Xd , mu.link.12 )
  mu.21 <- sapply(x.seq.Xa , mu.link.21 )
  mu.22 <- sapply(x.seq.Xa , mu.link.22 )


  mu.mean.11<-colMeans(mu.11)
  mu.mean.12<-colMeans(mu.12)
  mu.mean.21<-colMeans(mu.21)
  mu.mean.22<-colMeans(mu.22)


  mu.mean.11<-data.frame(as.numeric(mu.mean.11))
  mu.mean.12<-data.frame(as.numeric(mu.mean.12))
  names(mu.mean.11)<-"mu.mean.11"
  names(mu.mean.12)<-"mu.mean.12"

  mu.mean.21<-data.frame(as.numeric(mu.mean.21))
  mu.mean.22<-data.frame(as.numeric(mu.mean.22))
  names(mu.mean.21)<-"mu.mean.21"
  names(mu.mean.22)<-"mu.mean.22"



  mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
  mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )


  mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq.Xd)
  mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq.Xd)

  mu.CI.21 <- apply( mu.21 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
  mu.CI.22 <- apply( mu.22 , MARGIN=2, FUN=rethinking::PI , prob=0.89 )

  mu.CI.21<-data.frame(t(data.frame(mu.CI.21)),x.seq.Xa)
  mu.CI.22<-data.frame(t(data.frame(mu.CI.22)),x.seq.Xa)


  names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq.Xd")
  names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq.Xd")
  names(mu.CI.21)<-c("min.5.5","max.94.5","x.seq.Xa")
  names(mu.CI.22)<-c("min.5.5","max.94.5","x.seq.Xa")

  #df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
  df<-data.frame(Y=trdata$dat$Y_with_error,Xd=trdata$dat$Xd,Xa=trdata$dat$Xa,Regimes=trdata$dat$regimes)
  df11<-data.frame(x.seq.Xd,mu.mean.11)
  df12<-data.frame(x.seq.Xd,mu.mean.12)
  df21<-data.frame(x.seq.Xa,mu.mean.21)
  df22<-data.frame(x.seq.Xa,mu.mean.22)

  mypal <- ggsci::pal_npg("nrc", alpha = 0.7)(length(beta))

  slope.plot.1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xd,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_line(data=df11,ggplot2::aes(x=x.seq.Xd,y=mu.mean.11),linetype=1)+
    ggplot2::geom_abline(intercept=optima[1],slope=beta[1,1],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[2,1],alpha=0.5,linetype=2)+

    ggplot2::geom_ribbon(data=mu.CI.11,ggplot2::aes(x=x.seq.Xd,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df12,ggplot2::aes(x=x.seq.Xd,y=mu.mean.12),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.12,ggplot2::aes(x=x.seq.Xd,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()

  slope.plot.1

  slope.plot.2<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=Xa,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.15)+ #Prior
    ggplot2::geom_abline(intercept=optima[1],slope=beta[1,2],alpha=0.5,linetype=2)+
    ggplot2::geom_abline(intercept=optima[2],slope=beta[2,2],alpha=0.5,linetype=2)+

    ggplot2::geom_line(data=df21,ggplot2::aes(x=x.seq.Xa,y=mu.mean.21),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.21,ggplot2::aes(x=x.seq.Xa,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::geom_line(data=df22,ggplot2::aes(x=x.seq.Xa,y=mu.mean.22),linetype=1)+
    ggplot2::geom_ribbon(data=mu.CI.22,ggplot2::aes(x=x.seq.Xa,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+

    #ggtitle("Prior vs. Posterior for Intercept and Slope")+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()

  return(list(slope.plot.1,slope.plot.2))
}
#' reg.direct.mlm.ve.prior.post.plot.code- Create Prior vs. Posterior Plot for Multilevel Multi-Optima Direct Effect Model - Varying Effects - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param vary.effects Data frame with true optima and beta for each regime
#' @return Plots in ggplot2 format
#' @export
#'
reg.direct.mlm.ve.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,vary.effects){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  n_reg<-length(unique(trdata$dat$regimes))

  x.seq<- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)

  mu.link<-function(x.seq,i){optima.post[,i]+x.seq*beta.post[,i]}

  y.reg<-NULL
  y.reg.mu<-NULL
  saved.y.reg<-list()
  for(i in 1:n_reg){
    y.reg<-mapply(mu.link,x.seq,i)
    saved.y.reg[[i]]<-y.reg
    y.reg.mu<-rbind(y.reg.mu,colMeans(data.frame(y.reg)))
  }
  mu.dat<-cbind(x.seq,t(y.reg.mu))
  mu.dat<-data.frame(mu.dat)
  mu.CI.reg<-NULL
  mu.CIL.reg<-list()
  mu.CIU.reg<-list()
  for(i in 1:n_reg){
    names(mu.dat)[i+1]<-paste("OU",i,sep="")
    mu.CI.reg <- apply(saved.y.reg[[i]] , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
    mu.CIL.reg[[i]]<-mu.CI.reg[1,]
    mu.CIU.reg[[i]]<-mu.CI.reg[2,]
  }
  #return(mu.dat)

  df<-data.frame(Regimes=trdata$dat$regimes,Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error)
  sp1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.1) #Prior
  names.mu.dat<-names(mu.dat)
  j<-1
  for(i in names(mu.dat)[-1]){
      df.CI<-data.frame(x.seq=x.seq,Lower=mu.CIL.reg[[j]],Upper=mu.CIU.reg[[j]])
      sp1<-sp1+
        ggplot2::geom_line(data=mu.dat,ggplot2::aes(x=x.seq,y=.data[[i]]),alpha=0.5)+
        ggplot2::geom_ribbon(data=df.CI,ggplot2::aes(x=x.seq,ymin=Lower,ymax=Upper),linetype=2,alpha=0.25)
      j<-j+1
  }
  sp1<-sp1+
    ggplot2::geom_abline(data=vary.effects,ggplot2::aes(intercept=Intercept,slope=Slope),linetype=2,alpha=0.35)+
    ggplot2::theme_bw()+
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Direct Effect Model")+
    ggsci::scale_color_npg()
  return(sp1)
}
#' reg.adapt.mlm.ve.prior.post.plot.code- Create Prior vs. Posterior Plot for Multilevel Multi-Optima Adaptive Model - Varying Effects - Simulated Data
#' @param trdata Phylogeny and data in treeplyr format
#' @param optima.prior  Prior mean and scale parameter for optima
#' @param beta.prior Prior mean and scale parameter for half-life
#' @param post Posterior distribution in stanfit class
#' @param vary.effects Data frame with true optima and beta for each regime
#' @return Plots in ggplot2 format
#' @export
#'
reg.adapt.mlm.ve.prior.post.plot.code<-function(trdata,optima.prior,beta.prior,post,vary.effects){

  optima.sims<-rnorm(100,optima.prior[1],optima.prior[2])
  beta.sims<-rnorm(n=100,beta.prior[1],beta.prior[2])

  optima.post<-post$optima
  beta.post<-data.frame(post$beta)
  n_reg<-length(unique(trdata$dat$regimes))

  x.seq<- seq(from=min(trdata$dat$X_with_error), to=max(trdata$dat$X_with_error) , length.out=100)

  mu.link<-function(x.seq,i){optima.post[,i]+x.seq*beta.post[,i]}

  y.reg<-NULL
  y.reg.mu<-NULL
  saved.y.reg<-list()
  for(i in 1:n_reg){
    y.reg<-mapply(mu.link,x.seq,i)
    saved.y.reg[[i]]<-y.reg
    y.reg.mu<-rbind(y.reg.mu,colMeans(data.frame(y.reg)))
  }
  mu.dat<-cbind(x.seq,t(y.reg.mu))
  mu.dat<-data.frame(mu.dat)
  mu.CI.reg<-NULL
  mu.CIL.reg<-list()
  mu.CIU.reg<-list()
  for(i in 1:n_reg){
    names(mu.dat)[i+1]<-paste("OU",i,sep="")
    mu.CI.reg <- apply(saved.y.reg[[i]] , MARGIN=2, FUN=rethinking::PI , prob=0.89 )
    mu.CIL.reg[[i]]<-mu.CI.reg[1,]
    mu.CIU.reg[[i]]<-mu.CI.reg[2,]
  }
  #return(mu.dat)

  df<-data.frame(Regimes=trdata$dat$regimes,Y=trdata$dat$Y_with_error,X=trdata$dat$X_with_error)
  sp1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=df,ggplot2::aes(y=Y,x=X,color=Regimes))+
    ggplot2::geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.1) #Prior
  names.mu.dat<-names(mu.dat)
  j<-1
  for(i in names(mu.dat)[-1]){
    df.CI<-data.frame(x.seq=x.seq,Lower=mu.CIL.reg[[j]],Upper=mu.CIU.reg[[j]])
    sp1<-sp1+
      ggplot2::geom_line(data=mu.dat,ggplot2::aes(x=x.seq,y=.data[[i]]),alpha=0.5)+
      ggplot2::geom_ribbon(data=df.CI,ggplot2::aes(x=x.seq,ymin=Lower,ymax=Upper),linetype=2,alpha=0.25)
    j<-j+1
  }
  sp1<-sp1+
    ggplot2::geom_abline(data=vary.effects,ggplot2::aes(intercept=Intercept,slope=Slope),linetype=2,alpha=0.35)+
    ggplot2::theme_bw()+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::ylab("Y") + ggplot2::xlab("X - Adaptation Model")+
    ggsci::scale_color_npg()
  return(sp1)
}
