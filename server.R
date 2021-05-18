library(shiny)
library(shinyjs)
library(plyr)
library(car)
library(drc)
library(ggplot2)
library(tidyr)
library(gplots)
library(outliers)
library(scales)
library(rlist)
library(dplyr)

#### functions ######
mod <- function(df,norm) {
    conc_drug1=as.numeric(gsub(",",".",rownames(df)[2:nrow(df)]))
    conc_drug2=as.numeric(df[1,])
    cdrug2=rep(conc_drug2,8)
    cdrug1=c(rep(conc_drug1[1],24),rep(conc_drug1[2],24),rep(conc_drug1[3],24),
             rep(conc_drug1[4],24),rep(conc_drug1[5],24),rep(conc_drug1[6],24),
             rep(conc_drug1[7],24),rep(conc_drug1[8],24))
    fa=unlist(c(df[2,],df[3,],df[4,],df[5,],df[6,],df[7,],df[8,],df[9,]))

    dat=data.frame(drug1=cdrug1,drug2=cdrug2,fa=fa)
    dat=dat[!is.na(dat$fa),]
    dat<-dat[with(dat,order(dat[,1],dat[,2])),]

    ## outlier ##
    dat1=data.frame(drug1=numeric(),drug2=numeric(),fa=numeric())
    uni=unique(dat[,1:2])
    for (i in 1:nrow(uni)) {
        data=dat[dat$drug1==uni[i,1] & dat$drug2==uni[i,2],"fa"]
        if (length(data)>2) {
            if (grubbs.test(data,two.sided = T)$p.value<0.05) {
                data[-unlist(lapply(data, function(x) grep(x,grubbs.test(
                    data,two.sided = T)$alternative)))]
            }
        }
        new_data=data.frame(drug1=rep(uni[i,1],length(data)),
                            drug2=rep(uni[i,2],length(data)),
                            fa=data)
        dat1=rbind(dat1,new_data)
    }

    if (norm=="max") {maxSignal=max(dat1$fa)}
    if (norm=="zero") {maxSignal=mean(dat[dat1$drug1==0 & dat1$drug2==0,"fa"])}
    dat1$fa=1-((dat1$fa-min(dat1$fa))/(maxSignal-min(dat1$fa)))
    dat1$fa[dat1$fa<0]=0

    tmp<-dat1[,3]==0
    dat1[tmp,3]=0.0000001
    tmp<-dat1[,3]==1
    dat1[tmp,3]=0.9999999
    dat1=dat1[!(dat1$drug1==0 & dat1$drug2==0),]

    return(dat1)
}

mod_s = function(dosi,lettura,norm) {
    maxSignal=max(lettura)
    dat_drug=list()
    for (i in 1:ncol(dosi)) {
        conc_drug=as.numeric(dosi[,i])
        fa=as.numeric(lettura[,i])
        
        dat=data.frame(drug=conc_drug,fa=fa)
        dat<-dat[with(dat,order(dat[,1])),]
        
        ## outlier ##
        dat1=data.frame(drug=numeric(),fa=numeric())
        uni=unique(dat$drug)
        for (j in 1:length(uni)) {
            data=dat[dat$drug==uni[j],"fa"]
            if (length(data)>2) {
                if (grubbs.test(data,two.sided = T)$p.value<0.05) {
                    data[-unlist(lapply(data, function(x) grep(x,grubbs.test(
                        data,two.sided = T)$alternative)))]
                }
            }
            new_data=data.frame(drug=rep(uni[j],length(data)),fa=data)
            dat1=rbind(dat1,new_data)
        }
        
        if (norm=="zero") {maxSignal=mean(dat[dat1$drug==0,"fa"])}
        dat1$fa=(dat1$fa-min(dat1$fa))/(maxSignal-min(dat1$fa))
        dat1$fa[dat1$fa>1]=1
        
        tmp<-dat1[,2]==0
        dat1[tmp,2]=0.0000001
        tmp<-dat1[,2]==1
        dat1[tmp,2]=0.9999999
        dat1=dat1[!(dat1$drug==0),]
        
        dat_drug=list.append(dat_drug,dat1)
        
    }
    names(dat_drug)=colnames(dosi)
    return(dat_drug)
}

rsq = function (obs,pred,k) 1-((length(obs)-1)/(length(obs)-k-1))*
    (sum((obs-pred)^2)/sum((obs-(sum(obs)/length(obs)))^2))

s_drug= function(df) {

    colnames(df)=c("dose","viability")
    
    ## median effect (2 parameter) ##
    lm_med_eff <- lm(log(df$viability/(1-df$viability))~log(df$dose))
    dm_med_eff <- exp(-summary(lm_med_eff)$coef[1,1]/summary(lm_med_eff)$coef[2,1])
    m_med_eff=summary(lm_med_eff)$coef[2,1]
    x_med_eff=data.frame(fa=seq(0,1,length.out = 10000),
                         conc=dm_med_eff*(seq(0,1,length.out = 10000)/(1-seq(
                             0,1,length.out = 10000)))^(1/m_med_eff))
    fa_med_eff=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                          fa=1/(1+(dm_med_eff/seq(min(df$dose),max(df$dose),
                                                 length.out = 10000))^m_med_eff))
    r2_med_eff=rsq(df$viability,1/(1+(dm_med_eff/df$dose)^m_med_eff),2)
    ic50_med_eff=dm_med_eff*(0.5/(1-0.5))^(1/m_med_eff)
    
    ## log-logistic (4 parameter) ##
    fit=drm(df$viability ~ df$dose, data=df,fct = LL.4(names=c("m","min","max","dm")), 
            type = "continuous")
    m_log_log=fit$coefficients[1]
    min_log_log=fit$coefficients[2]
    max_log_log=fit$coefficients[3]
    dm_log_log=fit$coefficients[4]
    x_log_log=data.frame(fa=seq(0,1,length.out = 10000),conc=dm_log_log*(((
        max_log_log-min_log_log)/(seq(0,1,length.out = 10000)-min_log_log)
        )-1)^(1/m_log_log))
    fa_log_log=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                          fa=min_log_log+((max_log_log-min_log_log)/(1+(
                              seq(min(df$dose),max(df$dose),
                                  length.out = 10000)/dm_log_log)^(m_log_log))))
    
    r2_log_log=rsq(df$viability,min_log_log+((max_log_log-min_log_log)/(1+(
        df$dose/dm_log_log)^(m_log_log))),4)

    ic50_log_log=dm_log_log*(((max_log_log-min_log_log)/(0.5-min_log_log))-1
                             )^(1/m_log_log)
    
    ## log-logistic 0,1 (with 4 parameters coefficients) ##
    x_log_log_01=data.frame(fa=seq(0,1,length.out = 10000),
                            conc=dm_log_log*((1-seq(0,1,length.out = 10000))/
                                                 (seq(0,1,length.out = 10000)))^(
                                                     1/m_log_log))
    fa_log_log_01=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                             fa=1/(1+(seq(min(df$dose),max(df$dose),
                                          length.out = 10000)/dm_log_log)^m_log_log))
    
    r2_log_log_01=rsq(df$viability,1/(1+(df$dose/dm_log_log)^m_log_log),2)
    
    ic50_log_log_01=dm_log_log*((1-0.5)/(0.5))^(1/m_log_log)
    
    ## log-logistic 0 (3 parameter) ##
    fit_0=drm(df$viability ~ df$dose, data=df,fct = LL.3(names=c("m","max","dm")),
              type = "continuous")
    dm_log_log_0=summary(fit_0)$coeff[3,1]
    m_log_log_0=summary(fit_0)$coeff[1,1]
    max_log_log_0=summary(fit_0)$coeff[2,1]
    x_log_log_0=data.frame(fa=seq(0,1,length.out = 10000),
                           conc=dm_log_log_0*((max_log_log_0-seq(0,1,length.out=10000))/
                                                  (seq(0,1,length.out = 10000)))^(
                                                      1/m_log_log_0))
    fa_log_log_0=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                            fa=max_log_log_0/(1+(seq(min(df$dose),max(df$dose),
                                                     length.out = 10000)/dm_log_log_0
                                                 )^m_log_log_0))
    
    r2_log_log_0=rsq(df$viability,max_log_log_0/(1+(df$dose/dm_log_log_0)^m_log_log_0),3)

    ic50_log_log_0=dm_log_log_0*((max_log_log_0-0.5)/(0.5))^(1/m_log_log_0)
        
    ## log-logistic 1 (3 parameter) ##
    fit_1=drm(df$viability ~ df$dose, data=df,fct = LL.3u(names=c("m","min","dm")), 
              type = "continuous")
    dm_log_log_1=summary(fit_1)$coeff[3,1]
    m_log_log_1=summary(fit_1)$coeff[1,1]
    min_log_log_1=summary(fit_1)$coeff[2,1]
    x_log_log_1=data.frame(fa=seq(0,1,length.out = 10000),
                           conc=dm_log_log_1*(((1-min_log_log_1)/
                                                   (seq(0,1,length.out = 10000)-
                                                        min_log_log_1))-1)^(
                                                            1/m_log_log_1))
    fa_log_log_1=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                            fa=min_log_log_1+((1-min_log_log_1)/(1+(
                                seq(min(df$dose),max(df$dose),length.out = 10000)/
                                    dm_log_log_1)^(m_log_log_1))))
    
    r2_log_log_1=rsq(df$viability,min_log_log_1+((1-min_log_log_1)/(1+(
        df$dose/dm_log_log_1)^(m_log_log_1))),3)
 
   ic50_log_log_1=dm_log_log_1*(((1-min_log_log_1)/(0.5-min_log_log_1))-1)^(
       1/m_log_log_1)
    
    
    ## log-logistic 0,1 v2 (2 parameters) ##
    fit_01=drm(df$viability ~ df$dose, data=df,fct = LL.2(names=c("m","dm")),
               type = "continuous")
    dm_log_log_01=summary(fit_01)$coeff[2,1]
    m_log_log_01=summary(fit_01)$coeff[1,1]
    x_log_log_01_v2=data.frame(fa=seq(0,1,length.out = 10000),
                               conc=dm_log_log_01*((1-seq(0,1,length.out = 10000))/(
                                   seq(0,1,length.out = 10000)))^(1/m_log_log_01))
    fa_log_log_01_v2=data.frame(conc=seq(min(df$dose),max(df$dose),length.out = 10000),
                                fa=1/(1+(seq(min(df$dose),max(df$dose),length.out = 10000
                                             )/dm_log_log_01)^m_log_log_01))
    
    r2_log_log_01_v2=rsq(df$viability,1/(1+(df$dose/dm_log_log_01)^m_log_log_01),2)
    ic50_log_log_01_v2=dm_log_log_01*((1-0.5)/(0.5))^(1/m_log_log_01)
    
    return(list(df=df,
                med_eff=list("dm"=dm_med_eff,"m"=m_med_eff,"x"=x_med_eff,
                                   "fa"=fa_med_eff,"r2"=r2_med_eff,"ic50"=ic50_med_eff),
                log_log=list("m"=m_log_log,"min"=min_log_log,"max"=max_log_log,
                                   "dm"=dm_log_log,"x"=x_log_log,"fa"=fa_log_log,
                                   "r2"=r2_log_log,"ic50"=ic50_log_log),
                log_log_01=list("x"=x_log_log_01,"fa"=fa_log_log_01,
                                      "r2"=r2_log_log_01,"ic50"=ic50_log_log_01),
                log_log_0=list("dm"=dm_log_log_0,"m"=m_log_log_0,
                                     "max"=max_log_log_0,"x"=x_log_log_0,
                                     "fa"=fa_log_log_0,"r2"=r2_log_log_0,
                               "ic50"=ic50_log_log_0),
                log_log_1=list("dm"=dm_log_log_1,"m"=m_log_log_1,
                                     "min"=min_log_log_1,"x"=x_log_log_1,
                                     "fa"=fa_log_log_1,"r2"=r2_log_log_1,
                               "ic50"=ic50_log_log_1),
                log_log_01_v2=list("dm"=dm_log_log_01,"m"=m_log_log_01,
                                         "x"=x_log_log_01_v2,"fa"=fa_log_log_01_v2,
                                         "r2"=r2_log_log_01_v2,
                                   "ic50"=ic50_log_log_01_v2)))
}

plot_sdrug = function(mod_df,drug_name) {
    
    colors=c("blueviolet","coral","chartreuse2","lightpink",
             "darkturquoise","darkgoldenrod1")
    
    g_plot=ggplot(mod_df$df,aes(x=log(dose),y=viability)) +
                           geom_point(color='darkgray',size=2) +
                           labs(y="viability") + theme_light() +
                           scale_x_continuous(name="dose",
                                              breaks=log(sort(unique(mod_df$df$dose))),
                                              labels = sort(unique(mod_df$df$dose))) +
                           ggtitle(paste0("Drug ",drug_name)) + 
                           theme(plot.title=element_text(hjust=0.5,face="bold",size=24),
                                 axis.text=element_text(size=12),
                                 axis.title=element_text(size=14),
                                 legend.text=element_text(size=16),
                                 legend.title=element_blank(), 
                                 legend.position = "right",
                                 axis.text.x = element_text(angle = 35, hjust = 1),
                                 legend.box.background = element_rect(colour = "black"))+
                           geom_line(data=mod_df$log_log$fa,
                                     aes(x=log(conc),y=fa,colour="Log-logistic"),
                                     size=0.8) +
                           geom_line(data=mod_df$log_log_01$fa,
                                     aes(x=log(conc),y=fa,colour="Log-logistic[01]"),
                                     size=0.8) +
                           geom_line(data=mod_df$med_eff$fa,
                                     aes(x=log(conc),y=fa,colour="Median-effect"),
                                     size=0.8) +
                           geom_line(data=mod_df$log_log_01_v2$fa,
                                     aes(x=log(conc),y=fa,colour="Log-logistic[01]v2"),
                                     size=0.8) +
                           geom_line(data=mod_df$log_log_0$fa,
                                     aes(x=log(conc),y=fa,colour="Log-logistic[0]"),
                                     size=0.8) +
                           geom_line(data=mod_df$log_log_1$fa,
                                     aes(x=log(conc),y=fa,colour="Log-logistic[1]"),
                                     size=0.8) +
        geom_vline(xintercept=log(mod_df$log_log$ic50),
                   c,size=0.5,linetype="dashed", color="blueviolet") +
        geom_vline(xintercept=log(mod_df$log_log_0$ic50),
            c,size=0.5,linetype="dashed", color="coral") +
        geom_vline(xintercept=log(mod_df$log_log_01$ic50),
            c,size=0.5,linetype="dashed", color="chartreuse2") +
        geom_vline(xintercept=log(mod_df$log_log_01_v2$ic50),
                   c,size=0.5,linetype="dashed", color="lightpink") +
        geom_vline(xintercept=log(mod_df$log_log_1$ic50),
            c,size=0.5,linetype="dashed", color="darkturquoise") +
        geom_vline(xintercept=log(mod_df$med_eff$ic50),
            c,size=0.5,linetype="dashed", color="darkgoldenrod1") +
        labs(colour="Model") + scale_color_manual(values = colors)
    
    return(g_plot)
    }
    
plot_sdrug_d = function(mod_df,drug_name) {
    
    colors=c("blueviolet","coral","chartreuse2","lightpink",
             "darkturquoise","darkgoldenrod1")
    
    g_plot=ggplot(mod_df$df,aes(x=log(dose),y=viability)) +
        geom_point(color='darkgray',size=2) +
        labs(y="viability") + theme_light() +
        scale_x_continuous(name="dose",
                           breaks=log(sort(unique(mod_df$df$dose))),
                           labels = sort(unique(mod_df$df$dose))) +
        ggtitle(paste0("Drug ",drug_name)) + 
        theme(plot.title=element_text(hjust=0.5,face="bold",size=24),
              axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_blank(), 
              legend.position = "bottom",
              axis.text.x = element_text(angle = 35, hjust = 1),
              legend.box.background = element_rect(colour = "black"))+
        geom_line(data=mod_df$log_log$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0("Log-logistic\nR2=",
                                    round(mod_df$log_log$r2,2),
                                    "\nIC50=",
                                    signif(mod_df$log_log$ic50,3),
                                    "\n")),
                  size=0.8) +
        geom_line(data=mod_df$log_log_01$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0(
                          "Log-logistic[01]\nR2=",
                          round(mod_df$log_log_01$r2,2),"\nIC50=",
                          signif(mod_df$log_log_01$ic50,3),"\n")),
                  size=0.8) +
        geom_line(data=mod_df$med_eff$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0(
                          "Median-effect\nR2=",
                          round(mod_df$med_eff$r2,2),"\nIC50=",
                          signif(mod_df$med_eff$ic50,3))),
                  size=0.8) +
        geom_line(data=mod_df$log_log_01_v2$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0(
                          "Log-logistic[01] v2\nR2=",
                          round(mod_df$log_log_01_v2$r2,2),"\nIC50=",
                          signif(mod_df$log_log_01_v2$ic50,3))),
                  size=0.8) +
        geom_line(data=mod_df$log_log_0$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0(
                          "Log-logistic[0]\nR2=",
                          round(mod_df$log_log_0$r2,2),"\nIC50=",
                          signif(mod_df$log_log_0$ic50,3))),
                  size=0.8) +
        geom_line(data=mod_df$log_log_1$fa,
                  aes(x=log(conc),y=fa,
                      colour=paste0(
                          "Log-logistic[1]\nR2=",
                          round(mod_df$log_log_1$r2,2),"\nIC50=",
                          signif(mod_df$log_log_1$ic50,3),"\n")),
                  size=0.8) +
        geom_vline(xintercept=log(mod_df$log_log$ic50),
                   c,size=0.5,linetype="dashed", color="blueviolet") +
        geom_vline(xintercept=log(mod_df$log_log_0$ic50),
                   c,size=0.5,linetype="dashed", color="coral") +
        geom_vline(xintercept=log(mod_df$log_log_01$ic50),
                   c,size=0.5,linetype="dashed", color="chartreuse2") +
        geom_vline(xintercept=log(mod_df$log_log_01_v2$ic50),
                   c,size=0.5,linetype="dashed", color="lightpink") +
        geom_vline(xintercept=log(mod_df$log_log_1$ic50),
                   c,size=0.5,linetype="dashed", color="darkturquoise") +
        geom_vline(xintercept=log(mod_df$med_eff$ic50),
                   c,size=0.5,linetype="dashed", color="darkgoldenrod1") +
        labs(colour="Model") + scale_color_manual(values = colors)
    
    return(g_plot)
}

plot_cell=function(dosi1,cell_show,cell1,sep_ss1,dec_ss1,norm_ss,name1,cell2,sep_ss2,
                   dec_ss2,name2,cell3,sep_ss3,dec_ss3,name3,cell4,sep_ss4,dec_ss4,
                   name4,drug_list_s,rep_no,drug_ss) {
    dat_all=list()
    
    if ("cell_sh1" %in% cell_show) {
        infile2 <- cell1
        lettura=read.csv(infile2$datapath,sep=sep_ss1,dec=dec_ss1,
                         header = F,stringsAsFactors = F)
        lettura1=t(lettura)
        dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = norm_ss)
        dat_all=list.append(dat_all,dat_drug)
        names(dat_all)[length(dat_all)]=name1
    }
    
    if ("cell_sh2" %in% cell_show) {
        infile3 <- cell2
        lettura=read.csv(infile3$datapath,sep=sep_ss2,dec=dec_ss2,
                         header = F,stringsAsFactors = F)
        lettura1=t(lettura)
        dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = norm_ss)
        dat_all=list.append(dat_all,dat_drug)
        names(dat_all)[length(dat_all)]=name2
    }
    
    if ("cell_sh3" %in% cell_show) {
        infile4 <- cell3
        lettura=read.csv(infile4$datapath,sep=sep_ss3,dec=dec_ss3,
                         header = F,stringsAsFactors = F)
        lettura1=t(lettura)
        dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = norm_ss)
        dat_all=list.append(dat_all,dat_drug)
        names(dat_all)[length(dat_all)]=name3
    }
    
    if ("cell_sh4" %in% cell_show) {
        infile5 <- cell4
        lettura=read.csv(infile5$datapath,sep=sep_ss4,dec=dec_ss4,
                         header = F,stringsAsFactors = F)
        lettura1=t(lettura)
        dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = norm_ss)
        dat_all=list.append(dat_all,dat_drug)
        names(dat_all)[length(dat_all)]=name4
    }
    
    colors=c("blueviolet","coral","chartreuse2","lightpink")
    
    df_dat <- data.frame(drug=numeric(),fa=numeric(),cell_line=character(),
                         stringsAsFactors=FALSE)
    
    if("cell_sh1" %in% cell_show) df_dat=rbind(
        df_dat,cbind(dat_all[[name1]][[drug_list_s]],
                     cell_line=rep(name1,nrow(dat_all[[name1]][[drug_list_s]]))))
    if("cell_sh2" %in% cell_show) df_dat=rbind(
        df_dat,cbind(dat_all[[name2]][[drug_list_s]],
                     cell_line=rep(name2,nrow(dat_all[[name2]][[drug_list_s]]))))
    if("cell_sh3" %in% cell_show) df_dat=rbind(
        df_dat,cbind(dat_all[[name3]][[drug_list_s]],
                     cell_line=rep(name3,nrow(dat_all[[name3]][[drug_list_s]]))))
    if("cell_sh4" %in% cell_show) df_dat=rbind(
        df_dat,cbind(dat_all[[name4]][[drug_list_s]],
                     cell_line=rep(name4,nrow(dat_all[[name4]][[drug_list_s]]))))
    
    mod_df_all=list()
    if("cell_sh1" %in% cell_show) {
        df=dat_all[[name1]][[drug_list_s]]
        mod_df=s_drug(df)
        mod_df_all=list.append(mod_df_all,mod_df)
        names(mod_df_all)[length(mod_df_all)]=name1
    }
    
    if("cell_sh2" %in% cell_show) {
        df=dat_all[[name2]][[drug_list_s]]
        mod_df=s_drug(df)
        mod_df_all=list.append(mod_df_all,mod_df)
        names(mod_df_all)[length(mod_df_all)]=name2
    }
    
    if("cell_sh3" %in% cell_show) {
        df=dat_all[[name3]][[drug_list_s]]
        mod_df=s_drug(df)
        mod_df_all=list.append(mod_df_all,mod_df)
        names(mod_df_all)[length(mod_df_all)]=name3
    }
    
    if("cell_sh4" %in% cell_show) {
        df=dat_all[[name4]][[drug_list_s]]
        mod_df=s_drug(df)
        mod_df_all=list.append(mod_df_all,mod_df)
        names(mod_df_all)[length(mod_df_all)]=name4
    }
    
    if(rep_no=="rep") {cell_plot=ggplot(df_dat,aes(x=log(drug),y=fa))}
    if(rep_no=="mea") {
        df.summary2 <- df_dat %>%
            group_by(drug, cell_line) %>%
            summarise(
                sd = sd(fa),
                fa = mean(fa)
            )
        cell_plot=ggplot(df.summary2,aes(x=log(drug),y=fa)) +
            geom_errorbar(aes(ymin = fa-sd, ymax = fa+sd,colour=cell_line), 
                          data = df.summary2, width = 0.2)
    }
    
    cell_plot= cell_plot +
        geom_point(aes(colour=cell_line),size=2) +
        labs(y="viability") + theme_light() +
        scale_x_continuous(name="dose",
                           breaks=log(sort(unique(dat_all[[1]][[drug_list_s]]$drug))),
                           labels = sort(unique(dat_all[[1]][[drug_list_s]]$drug))) +
        ggtitle(paste0("Drug ",drug_list_s)) + 
        theme(plot.title=element_text(hjust=0.5,face="bold",size=24),
              axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_blank(), 
              legend.position = "right",
              axis.text.x = element_text(angle = 35, hjust = 1),
              legend.box.background = element_rect(colour = "black"))
    
    if ("cell_sh1" %in% cell_show) {
        cell_plot <- cell_plot + 
            geom_line(data=mod_df_all[[name1]][[drug_ss]]$fa, 
                      aes(x=log(conc),y=fa,colour=name1),size=0.8) +
            geom_vline(xintercept=log(mod_df_all[[name1]][[drug_ss]]$ic50),c,
                       size=0.5,linetype="dashed",color=colors[1])
    }
    
    if ("cell_sh2" %in% cell_show) {
        cell_plot <- cell_plot + 
            geom_line(data=mod_df_all[[name2]][[drug_ss]]$fa, 
                      aes(x=log(conc),y=fa,colour=name2),size=0.8) +
            geom_vline(xintercept=log(mod_df_all[[name2]][[drug_ss]]$ic50),c,
                       size=0.5,linetype="dashed",color=colors[2])
    }
    
    if ("cell_sh3" %in% cell_show) {
        cell_plot <- cell_plot + 
            geom_line(data=mod_df_all[[name3]][[drug_ss]]$fa, 
                      aes(x=log(conc),y=fa,colour=name3),size=0.8) +
            geom_vline(xintercept=log(mod_df_all[[name3]][[drug_ss]]$ic50),c,
                       size=0.5,linetype="dashed",color=colors[3])
    }
    
    if ("cell_sh4" %in% cell_show) {
        cell_plot <- cell_plot + 
            geom_line(data=mod_df_all[[name4]][[drug_ss]]$fa, 
                      aes(x=log(conc),y=fa,colour=name4),size=0.8) +
            geom_vline(xintercept=log(mod_df_all[[name4]][[drug_ss]]$ic50),c,
                       size=0.5,linetype="dashed",color=colors[4])
    }
    
    cols=c(colors[1:4])
    names(cols)=c(name1,name2,name3,name4)
    gg=c("cell_sh1","cell_sh2","cell_sh3","cell_sh4") %in% cell_show
    
    cell_plot <- cell_plot + 
        scale_colour_manual("Cell lines",
                            breaks = names(cols)[gg],values = cols[gg])
    return(cell_plot)
}

plot_r2=function(x,norm,drug_names) {
    dati=mod(x,norm)
    drug1=drug_names[1]
    drug2=drug_names[2]
    df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
    df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
    mod_df1=s_drug(df1)
    mod_df2=s_drug(df2)
    r2_df=data.frame(r2_drug1=c(mod_df1$log_log$r2,mod_df1$log_log_0$r2,
                                mod_df1$log_log_01$r2,mod_df1$log_log_01_v2$r2,
                                mod_df1$log_log_1$r2,mod_df1$med_eff$r2),
                     r2_drug2=c(mod_df2$log_log$r2,mod_df2$log_log_0$r2,
                                mod_df2$log_log_01$r2,mod_df2$log_log_01_v2$r2,
                                mod_df2$log_log_1$r2,mod_df2$med_eff$r2))
    
    r2_df$Model=c("Log-logistic","Log-logistic [0]","Log-logistic [01]",
                  "Log-logistic [01] v2","Log-logistic [1]","Median-effect")
    
    colors=c("blueviolet","coral","chartreuse2","lightpink",
             "darkturquoise","darkgoldenrod1")
    
    p=ggplot(r2_df,aes(x=r2_drug1, y=r2_drug2,colour=Model)) + geom_point(size=2) + 
        geom_vline(xintercept=r2_df$r2_drug1,c,size=0.4,linetype="dashed",
                   color=colors) +
        geom_hline(yintercept=r2_df$r2_drug2,c,size=0.4,linetype="dashed",
                   color=colors) +
        labs(x = paste0("R2 ",drug1),y = paste0("R2 ",drug2))  +
        scale_color_manual(values = colors) + theme_light() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_blank(), 
              legend.position = "bottom",
              legend.box.background = element_rect(colour = "black"))
    return(list(p,r2_df))
}

ci=function(df,mod1,mod2) {
    sum_data<-ddply(df[df$drug1!=0 & df$drug2!=0,],c(drug1="drug1",drug2="drug2"),
                      summarise,fa=mean(fa))
    sum_data<-sum_data[with(sum_data,order(sum_data[,1],sum_data[,2])),]
    
    ## median effect ##
    d_med_eff_1=mod1$med_eff$dm*(sum_data$fa/(1-sum_data$fa))^(1/mod1$med_eff$m)
    d_med_eff_2=mod2$med_eff$dm*(sum_data$fa/(1-sum_data$fa))^(1/mod2$med_eff$m)
    ci_med_eff=(sum_data$drug1/d_med_eff_1)+(sum_data$drug2/d_med_eff_2)
    dff_med_eff<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,ci=ci_med_eff)
    
    ## log-logistic ##
    d_log_log_1=mod1$log_log$dm*(((mod1$log_log$max-mod1$log_log$min)/(
        sum_data$fa-mod1$log_log$min))-1)^(1/mod1$log_log$m)
    d_log_log_2=mod2$log_log$dm*(((mod2$log_log$max-mod2$log_log$min)/(
        sum_data$fa-mod2$log_log$min))-1)^(1/mod2$log_log$m)
    ci_log_log=(sum_data$drug1/d_log_log_1)+(sum_data$drug2/d_log_log_2)
    dff_log_log<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,ci=ci_log_log)
    
    ## log-logistic 01 ##
    d_log_log_01_1=mod1$log_log$dm*((1-sum_data$fa)/(sum_data$fa))^(1/mod1$log_log$m)
    d_log_log_01_2=mod2$log_log$dm*((1-sum_data$fa)/(sum_data$fa))^(1/mod2$log_log$m)
    ci_log_log_01=(sum_data$drug1/d_log_log_01_1)+(sum_data$drug2/d_log_log_01_2)
    dff_log_log_01<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,
                               ci=ci_log_log_01)
    
    ## log-logistic 0 ##
    d_log_log_0_1=mod1$log_log_0$dm*((mod1$log_log_0$max-sum_data$fa)/(sum_data$fa))^(
        1/mod1$log_log_0$m)
    d_log_log_0_2=mod2$log_log_0$dm*((mod2$log_log_0$max-sum_data$fa)/(sum_data$fa))^(
        1/mod2$log_log_0$m)
    ci_log_log_0=(sum_data$drug1/d_log_log_0_1)+(sum_data$drug2/d_log_log_0_2)
    dff_log_log_0<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,ci=ci_log_log_0)
    
    ## log-logistic 1 ##
    d_log_log_1_1=mod1$log_log_1$dm*(((1-mod1$log_log_1$min)/(
        sum_data$fa-mod1$log_log_1$min))-1)^(1/mod1$log_log_1$m)
    d_log_log_1_2=mod2$log_log_1$dm*(((1-mod2$log_log_1$min)/(
        sum_data$fa-mod2$log_log_1$min))-1)^(1/mod2$log_log_1$m)
    ci_log_log_1=(sum_data$drug1/d_log_log_1_1)+(sum_data$drug2/d_log_log_1_2)
    dff_log_log_1<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,ci=ci_log_log_1)
    
    ## log-logistic 01 v2 ##
    d_log_log_01_v2_1=mod1$log_log_01_v2$dm*((1-sum_data$fa)/(sum_data$fa))^(
        1/mod1$log_log_01_v2$m)
    d_log_log_01_v2_2=mod2$log_log_01_v2$dm*((1-sum_data$fa)/(sum_data$fa))^(
        1/mod2$log_log_01_v2$m)
    ci_log_log_01_v2=(sum_data$drug1/d_log_log_01_v2_1)+(
        sum_data$drug2/d_log_log_01_v2_2)
    dff_log_log_01_v2<-data.frame(drug1=sum_data$drug1,drug2=sum_data$drug2,
                                  ci=ci_log_log_01_v2)
    
    return(list("Median-effect"=dff_med_eff,"Log-logistic"=dff_log_log,
                "Log-logistic [01]"=dff_log_log_01,"Log-logistic [0]"=dff_log_log_0,
                "Log-logistic [1]"=dff_log_log_1,
                "Log-logistic [01] v2"=dff_log_log_01_v2))
}

isobol_plot=function(df_ci,drug,model) {
    df_ci$syn<-rep(NA,length(df_ci$ci))
    df_ci$syn<-ifelse(df_ci$ci>=1,"Antagonism","Synergy")
    df_ci$syn[df_ci$ci<1.1 & df_ci$ci>0.9]<-"Additivity"
    df_ci$syn<-factor(df_ci$syn,levels=c("Synergy","Antagonism","Additivity"))
    isob=ggplot(df_ci,aes(x=factor(signif(drug1,2)), y=factor(signif(drug2,2)), 
                      col=factor(syn,levels=c("Synergy","Antagonism",
                                              "Additivity"))))+
        geom_point(aes(size=ifelse(!is.finite(ci),
                                   log(max(ci[!is.infinite(ci)],na.rm = T)),
                                   abs(log(ci)))),
                   shape=19,na.rm = T)+ guides(size = FALSE)+
        scale_color_manual("Drug-Drug interaction", 
                           breaks=c("Synergy","Antagonism","Additivity"),
                           values=c("Synergy"="dodgerblue",
                                    "Antagonism"="firebrick1",
                                    "Additivity"="black"))+
        scale_size_continuous("Interaction Strength", range=c(1,10))+
        xlab(paste0("\n",colnames(df_ci)[1]," concentration (",drug[1],")"))+
        ylab(paste0(colnames(df_ci)[2]," concentration (",drug[2],")\n"))+
        ggtitle("Isobologram",subtitle = paste0(model," plot\n"))+
        theme(
            panel.background = element_blank(),
            legend.text = element_text(size = 10),
            plot.title = element_text(face="bold",colour="black",hjust = 0.5),
            plot.subtitle = element_text(face="bold",colour="black",hjust = 0.5),
            legend.position = "bottom",
            panel.grid.major = element_line(linetype = "dotted",colour = "grey"),
            panel.grid.minor = element_line(linetype = "dotted",colour = "grey50")
        )
    
    return(isob)
}

heat_plot=function(m,drug,model) {
    m$ci=log(m$ci)
    
    p=spread(m[,1:3], drug1, ci, fill=0)
    rownames(p)=p$drug2
    p=p[,2:ncol(p)]
    
    p=p[order(as.numeric(rownames(p)),decreasing = T),]
    p=p[order(as.numeric(colnames(p))),]
    
    p=as.matrix(p)
    colnames(p)=signif(as.numeric(colnames(p)),2)
    rownames(p)=signif(as.numeric(rownames(p)),2)
    
    my_palette <- c(colorRampPalette(c("dodgerblue","dodgerblue","white"))(n = 50),
                    colorRampPalette(c("white", "gainsboro","white"))(n = 10),
                    colorRampPalette(c("white","firebrick1","firebrick1"))(n = 50))
    breaksList = c(seq(-max(abs(p[!is.infinite(p)]),na.rm = T),log(0.89),length.out= 50),
                   seq(log(0.9),log(1.1), length.out = 11),
                   seq(log(1.11),max(abs(p[!is.infinite(p)]),na.rm = T),length.out = 50))
    par(cex.main=1)
    heatmap.2(p, dendrogram = "none", scale="none", Rowv = F, Colv = F, 
              trace="none", density.info = "none",col = my_palette,notecol ="black",
              main=paste0("Log CI matrix\n\n ",model," plot"), symbreaks = T,
              symkey = T, cexCol = 1, cexRow = 1,cellnote = round(p,2),
              xlab = paste0("\n",colnames(m)[1]," concentration (",drug[1],")"),
              ylab = ylab(paste0(colnames(m)[2]," concentration (",drug[2],")\n")),
              margins=c(5,6.5),breaks = breaksList)
}

##### shiny #####
shinyServer(function(input, output, session) {
    
    v <- reactiveValues(data = NULL)
    v1 <- reactiveValues(data = NULL)
    v_s <- reactiveValues(data = NULL)
    observeEvent(input$go1, {v$data <- 100})
    observeEvent(input$go2, {v1$data <- 100})
    observeEvent(input$go_s, {v_s$data <- 100})
    
    output$fileUploaded <- reactive({
        return(is.null(input$file1))
    })
    outputOptions(output, 'fileUploaded',suspendWhenHidden=FALSE)

    output$fileUploaded <- reactive({
        return(is.null(input$file_s1))
    })
    outputOptions(output, 'fileUploaded',suspendWhenHidden=FALSE)

    output$fileUploaded <- reactive({
        return(is.null(input$file_s2))
    })
    outputOptions(output, 'fileUploaded',suspendWhenHidden=FALSE)
    
    output$fileUploaded <- reactive({
        return(is.null(input$dose1))
    })
    outputOptions(output, 'fileUploaded',suspendWhenHidden=FALSE)

    output$fileUploaded <- reactive({
        return(is.null(input$cell1))
    })
    outputOptions(output, 'fileUploaded',suspendWhenHidden=FALSE)
    
    observe({
        if (is.null(input$file_s1)) {return(NULL)}
        infile1 <- input$file_s1
        dosi=read.csv(infile1$datapath,sep=input$sep_s1,dec=input$dec_s1,header = F,
                      stringsAsFactors = F)
        x <- dosi$V1
        updateSelectizeInput(session, "drug_list",choices = x,selected = head(x, 1)
        )
    })
    
    observe({
        if (is.null(input$dose1)) {return(NULL)}
        infile1 <- input$dose1
        dosi=read.csv(infile1$datapath,sep=input$sep_ss0,dec=input$dec_ss0,header = F,
                      stringsAsFactors = F)
        x <- dosi$V1
        updateSelectizeInput(session, "drug_list_s",choices = x,selected = head(x, 1)
        )
    })

    ###### single drug #####
    output$drug_s1 <- renderPlot({
        if (is.null(input$file_s1)) {return(NULL)}
        if (is.null(input$file_s2)) {return(NULL)}
        withProgress(message = 'Loading Plot', value = 1, {
            infile1 <- input$file_s1
            infile2 <- input$file_s2
            dosi=read.csv(infile1$datapath,sep=input$sep_s1,dec=input$dec_s1,header = F,
                          stringsAsFactors = F)
            dosi1=t(dosi[,2:ncol(dosi)])
            colnames(dosi1)=dosi$V1
            lettura=read.csv(infile2$datapath,sep=input$sep_s2,dec=input$dec_s2,
                             header = F,stringsAsFactors = F)
            lettura1=t(lettura)
            dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = input$norm_s)
            colors=c("blueviolet","coral","chartreuse2","lightpink",
                     "darkturquoise","darkgoldenrod1")
            
            df=dat_drug[[input$drug_list]]
            mod_df=s_drug(df)
            drug_plot=plot_sdrug(mod_df,input$drug_list)
            plot(drug_plot)
        })
    })
    
    output$table_s1 <- renderTable({
        if (is.null(input$file_s1)) {return(NULL)}
        if (is.null(input$file_s2)) {return(NULL)}
        withProgress(message = 'Loading Table', value = 1, {
            infile1 <- input$file_s1
            infile2 <- input$file_s2
            dosi=read.csv(infile1$datapath,sep=input$sep_s1,dec=input$dec_s1,header = F,
                          stringsAsFactors = F)
            dosi1=t(dosi[,2:ncol(dosi)])
            colnames(dosi1)=dosi$V1
            lettura=read.csv(infile2$datapath,sep=input$sep_s2,dec=input$dec_s2,
                             header = F,stringsAsFactors = F)
            lettura1=t(lettura)
            dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = input$norm_s)
            df=dat_drug[[input$drug_list]]
            mod_df=s_drug(df)
            table=data.frame("Model"=c("Log-logistic","Log-logistic[0]",
                                       "Log-logistic[01]","Log-logistic[01]v2",
                                       "Log-logistic[1]","Median-effect"),
                             "IC50"=c(mod_df$log_log$ic50,mod_df$log_log_0$ic50,
                                      mod_df$log_log_01$ic50,mod_df$log_log_01_v2$ic50,
                                      mod_df$log_log_1$ic50,mod_df$med_eff$ic50),
                             "R2"=c(mod_df$log_log$r2,mod_df$log_log_0$r2,
                                    mod_df$log_log_01$r2,mod_df$log_log_01_v2$r2,
                                    mod_df$log_log_1$r2,mod_df$med_eff$r2))
            table
        })
    },rownames = F,colnames = T,digits = 6,sanitize.text.function=function(x){x})
    
    output$downloadAll_s <- downloadHandler(
        filename = function(){
            paste0(gsub(".csv$","",input$file_s2),".pdf")},
        content = function(file){
            if (is.null(input$file_s1)) {return(NULL)}
            if (is.null(input$file_s2)) {return(NULL)}
            withProgress(message = 'Preparing Download', min=0, max=1, {
                
                infile1 <- input$file_s1
                infile2 <- input$file_s2
                dosi=read.csv(infile1$datapath,sep=input$sep_s1,dec=input$dec_s1,header = F,
                              stringsAsFactors = F)
                dosi1=t(dosi[,2:ncol(dosi)])
                colnames(dosi1)=dosi$V1
                lettura=read.csv(infile2$datapath,sep=input$sep_s2,dec=input$dec_s2,
                                 header = F,stringsAsFactors = F)
                lettura1=t(lettura)
                dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = input$norm_s)
                colors=c("blueviolet","coral","chartreuse2","lightpink",
                         "darkturquoise","darkgoldenrod1")
                pdf(file)
                for (i in 1:length(dat_drug)) {
                    incProgress(1/length(dat_drug))
                    tryCatch({
                df=dat_drug[[i]]
                mod_df=s_drug(df)
                drug_plot=plot_sdrug_d(mod_df,names(dat_drug)[i])
                print(drug_plot)
                    }, error=function(e){})
                }
            dev.off()
            })
        })
    
    output$downloadAll_s1 <- downloadHandler(
        filename = function(){"Plot.pdf"},
        content = function(file){
            if (is.null(input$file_s1)) {return(NULL)}
            if (is.null(input$file_s2)) {return(NULL)}
            withProgress(message = 'Preparing Download', value=1, {
                infile1 <- input$file_s1
                infile2 <- input$file_s2
                dosi=read.csv(infile1$datapath,sep=input$sep_s1,dec=input$dec_s1,header = F,
                              stringsAsFactors = F)
                dosi1=t(dosi[,2:ncol(dosi)])
                colnames(dosi1)=dosi$V1
                lettura=read.csv(infile2$datapath,sep=input$sep_s2,dec=input$dec_s2,
                                 header = F,stringsAsFactors = F)
                lettura1=t(lettura)
                dat_drug=mod_s(dosi = dosi1,lettura = lettura1,norm = input$norm_s)
                colors=c("blueviolet","coral","chartreuse2","lightpink",
                         "darkturquoise","darkgoldenrod1")
                pdf(file)
                        df=dat_drug[[input$drug_list]]
                        mod_df=s_drug(df)
                        drug_plot=plot_sdrug_d(mod_df,input$drug_list)
                        print(drug_plot)
                dev.off()
            })
        })
    
    ###### drug comparison ####
    output$drug_cell <- renderPlot({
        if (is.null(input$dose1)) {return(NULL)}
        if (length(input$cell_show)==0) {return(NULL)}
        withProgress(message = 'Loading Plot', value = 1, {
            infile1 <- input$dose1
            dosi=read.csv(infile1$datapath,sep=input$sep_ss0,dec=input$dec_ss0,
                          header = F,stringsAsFactors = F)
            dosi1=t(dosi[,2:ncol(dosi)])
            colnames(dosi1)=dosi$V1
            if ("cell_sh1" %in% input$cell_show) if (is.null(input$cell1)) {return(NULL)}
            if ("cell_sh2" %in% input$cell_show) if (is.null(input$cell2)) {return(NULL)}
            if ("cell_sh3" %in% input$cell_show) if (is.null(input$cell3)) {return(NULL)}
            if ("cell_sh4" %in% input$cell_show) if (is.null(input$cell4)) {return(NULL)}
            cell_plot=plot_cell(dosi1,input$cell_show,input$cell1,input$sep_ss1,
                                input$dec_ss1,input$norm_ss,input$name1,input$cell2,
                                input$sep_ss2,input$dec_ss2,input$name2,input$cell3,
                                input$sep_ss3,input$dec_ss3,input$name3,input$cell4,
                                input$sep_ss4,input$dec_ss4,input$name4,
                                input$drug_list_s,input$rep_no,input$drug_ss)
            plot(cell_plot)
        })
        
        
    })
    
    output$do <- downloadHandler(
        filename = function(){
            paste0(input$drug_list_s,".pdf")},
        content = function(file){
            withProgress(message = 'Preparing Download', value=1, {
                infile1 <- input$dose1
                dosi=read.csv(infile1$datapath,sep=input$sep_ss0,dec=input$dec_ss0,
                              header = F,stringsAsFactors = F)
                dosi1=t(dosi[,2:ncol(dosi)])
                colnames(dosi1)=dosi$V1
                if ("cell_sh1" %in% input$cell_show) if (is.null(input$cell1)) {return(NULL)}
                if ("cell_sh2" %in% input$cell_show) if (is.null(input$cell2)) {return(NULL)}
                if ("cell_sh3" %in% input$cell_show) if (is.null(input$cell3)) {return(NULL)}
                if ("cell_sh4" %in% input$cell_show) if (is.null(input$cell4)) {return(NULL)}
                cell_plot=plot_cell(dosi1,input$cell_show,input$cell1,input$sep_ss1,
                                    input$dec_ss1,input$norm_ss,input$name1,input$cell2,
                                    input$sep_ss2,input$dec_ss2,input$name2,input$cell3,
                                    input$sep_ss3,input$dec_ss3,input$name3,input$cell4,
                                    input$sep_ss4,input$dec_ss4,input$name4,
                                    input$drug_list_s,input$rep_no,input$drug_ss)
                pdf(file,10,8)
                plot(cell_plot)
                dev.off()
            })
        })
    
    output$downloadAll_ss1 <- downloadHandler(
        filename = function(){"Drug_comp.pdf"},
        content = function(file){
            withProgress(message = 'Preparing Download', min=0, max=1, {
                infile1 <- input$dose1
                dosi=read.csv(infile1$datapath,sep=input$sep_ss0,dec=input$dec_ss0,
                              header = F,stringsAsFactors = F)
                dosi1=t(dosi[,2:ncol(dosi)])
                colnames(dosi1)=dosi$V1
                pdf(file,10,8)
                for (i in 1:length(colnames(dosi1))) {
                    incProgress(1/length(colnames(dosi1)))
                    tryCatch({
                if ("cell_sh1" %in% input$cell_show) if (is.null(input$cell1)) {return(NULL)}
                if ("cell_sh2" %in% input$cell_show) if (is.null(input$cell2)) {return(NULL)}
                if ("cell_sh3" %in% input$cell_show) if (is.null(input$cell3)) {return(NULL)}
                if ("cell_sh4" %in% input$cell_show) if (is.null(input$cell4)) {return(NULL)}
                cell_plot=plot_cell(dosi1,input$cell_show,input$cell1,input$sep_ss1,
                                    input$dec_ss1,input$norm_ss,input$name1,input$cell2,
                                    input$sep_ss2,input$dec_ss2,input$name2,input$cell3,
                                    input$sep_ss3,input$dec_ss3,input$name3,input$cell4,
                                    input$sep_ss4,input$dec_ss4,input$name4,
                                    colnames(dosi1)[i],input$rep_no,input$drug_ss)
                print(cell_plot)
                    }, error=function(e){})
                }
                dev.off()
            })
        })
    
    
    ##### drug combination #####
        output$table <- renderTable({
        if (is.null(input$file1)) {return(NULL)}
        withProgress(message = 'Loading Table', value = 1, {  
        infile <- input$file1
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        tab=mod(dati1,input$norm)
        mat=matrix(c(0,0,0,tab$fa),nrow = 8,ncol = 24,byrow = T)
        rownames(mat)=paste0("<strong>",gsub(",",".",rownames(dati1)[2:nrow(dati1)]),
                             "</strong>")
        colnames(mat)=dati1[1,]
        sum_dati1=sapply(unique(colnames(mat)), function(g) rowMeans(
            mat[,colnames(mat)==g,drop=FALSE]))
        })
    },rownames = T,colnames = T,digits = 6,sanitize.text.function=function(x){x})
    
    output$table_r2 <- renderTable({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        plot_r2(dati1,input$norm,c("drug1","drug2"))[[2]]
    })
    
    output$drug1 <- renderPlot({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        withProgress(message = 'Loading Plot', value = 1, {  
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        drug1=strsplit(rownames(dati1)[1],",")[[1]][1]
        dati=mod(dati1,input$norm)
        df=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
        mod_df=s_drug(df)
        p=plot_sdrug(mod_df,drug1)
        if (!is.null(v$data)) p=p+ggtitle(paste0("Drug ",input$name1))
        print(p)
        })
    })

    output$drug2 <- renderPlot({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        withProgress(message = 'Loading Plot', value = 1, {  
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        drug2=strsplit(rownames(dati1)[1],",")[[1]][2]
        dati=mod(dati1,input$norm)
        df=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
        mod_df=s_drug(df)
        p=plot_sdrug(mod_df,drug2)
        if (!is.null(v1$data)) p=p+ggtitle(paste0("Drug ",input$name2))
        print(p)
        })
    })

    output$r2_plot <- renderPlot({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        withProgress(message = 'Loading Plot', value = 1, {  
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            drug1=strsplit(rownames(dati1)[1],",")[[1]][1]
            drug2=strsplit(rownames(dati1)[1],",")[[1]][2]
            p=plot_r2(dati1,input$norm,c(drug1,drug2))[[1]]
            if (!is.null(v$data)) p=p+labs(x = paste0("R2 ",input$name1))
            if (!is.null(v1$data)) p=p+labs(y = paste0("R2 ",input$name2))
            print(p)
        })
    })
    
    output$downloadPlot <- downloadHandler(
        filename = function(){"R2_plot.png"},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            if (input$name1!="") {drug1=input$name1
            } else if (input$name1=="") {drug1=strsplit(rownames(dati1)[1],",")[[1]][1]}
            if (input$name2!="") {drug2=input$name2
            } else if (input$name2=="") {drug2=strsplit(rownames(dati1)[1],",")[[1]][2]}
            ggsave(file,plot_r2(dati1,input$norm,c(drug1,drug2))[[1]])
        }
    )
    
    output$downloadPlot1 <- downloadHandler(
        filename = function(){"drug1_plot.png"},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            if (input$name1!="") {drug1=input$name1
            } else if (input$name1=="") {drug1=strsplit(rownames(dati1)[1],",")[[1]][1]}
            dati=mod(dati1,input$norm)
            df=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
            mod_df=s_drug(df)
            ggsave(file,plot_sdrug(mod_df,drug1))
        }
    )
    output$downloadPlot2 <- downloadHandler(
        filename = function(){"drug2_plot.png"},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            if (input$name2!="") {drug2=input$name2
            } else if (input$name2=="") {drug2=strsplit(rownames(dati1)[1],",")[[1]][2]}
            dati=mod(dati1,input$norm)
            df=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
            mod_df=s_drug(df)
            ggsave(file,plot_sdrug(mod_df,drug2))
        }
    )
    
    output$plot1a <- renderPlot({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        withProgress(message = 'Loading Plot', value = 1, {  
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        dati=mod(dati1,input$norm)
        drug1=strsplit(rownames(dati1)[1],",")[[1]][1]
        drug2=strsplit(rownames(dati1)[1],",")[[1]][2]
        df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
        df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
        mod_df1=s_drug(df1)
        mod_df2=s_drug(df2)
        ci=ci(dati,mod_df1,mod_df2)
        p=isobol_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug)
        if (!is.null(v$data)) p=p+xlab(paste0("\ndrug1 concentration (",input$name1,")"))
        if (!is.null(v1$data)) p=p+ylab(paste0("drug2 concentration (",input$name2,")\n"))
        print(p)
        })
    })
    
    output$plot1b <- renderPlot({
        if (is.null(input$file1)) {return(NULL)}
        infile <- input$file1
        withProgress(message = 'Loading Plot', value = 1, {  
        dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                       dec = input$dec)
        dati=mod(dati1,input$norm)
        df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
        df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
        mod_df1=s_drug(df1)
        mod_df2=s_drug(df2)
        ci=ci(dati,mod_df1,mod_df2)
        if (is.null(v$data) & is.null(v1$data)) {
            drug1=strsplit(rownames(dati1)[1],",")[[1]][1]
            drug2=strsplit(rownames(dati1)[1],",")[[1]][2]
            print(heat_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug))
        }
        if (!is.null(v$data) & is.null(v1$data)) {
            drug1=input$name1
            drug2=strsplit(rownames(dati1)[1],",")[[1]][2]
            print(heat_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug))
        }
        if (is.null(v$data) & !is.null(v1$data)) {
            drug1=strsplit(rownames(dati1)[1],",")[[1]][1]
            drug2=input$name2
            print(heat_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug))
        }
        if (!is.null(v$data) & !is.null(v1$data)) {
            drug1=input$name1
            drug2=input$name2
            print(heat_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug))
        }
        })
    })
    
    output$downloadPlot1a <- downloadHandler(
        filename = function(){"ci_plot1.png"},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            dati=mod(dati1,input$norm)
            if (input$name1!="") {drug1=input$name1
            } else if (input$name1=="") {drug1=strsplit(rownames(dati1)[1],",")[[1]][1]}
            if (input$name2!="") {drug2=input$name2
            } else if (input$name2=="") {drug2=strsplit(rownames(dati1)[1],",")[[1]][2]}
            df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
            df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
            mod_df1=s_drug(df1)
            mod_df2=s_drug(df2)
            ci=ci(dati,mod_df1,mod_df2)
            ggsave(file,
                   isobol_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug))
        }
    )
    
    output$downloadPlot1b <- downloadHandler(
        filename = function(){"ci_plot2.png"},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            dati=mod(dati1,input$norm)
            if (input$name1!="") {drug1=input$name1
            } else if (input$name1=="") {drug1=strsplit(rownames(dati1)[1],",")[[1]][1]}
            if (input$name2!="") {drug2=input$name2
            } else if (input$name2=="") {drug2=strsplit(rownames(dati1)[1],",")[[1]][2]}
            df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
            df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
            mod_df1=s_drug(df1)
            mod_df2=s_drug(df2)
            ci=ci(dati,mod_df1,mod_df2)
            png(file,w=9, h=9, pointsize=15, res=300, units = "in")
            heat_plot(ci[[input$drug]],drug=c(drug1,drug2),model=input$drug)
            dev.off()
        }
    )
    
    output$downloadAll <- downloadHandler(
        filename = function(){
            paste0(gsub(".csv$","",input$file1),".pdf")},
        content = function(file){
            if (is.null(input$file1)) {return(NULL)}
            infile <- input$file1
            dati1=read.csv(infile$datapath,sep=input$sep,row.names = 1,header=F,
                           dec = input$dec)
            if (input$name1!="") {drug1=input$name1
            } else if (input$name1=="") {drug1=strsplit(rownames(dati1)[1],",")[[1]][1]}
            if (input$name2!="") {drug2=input$name2
            } else if (input$name2=="") {drug2=strsplit(rownames(dati1)[1],",")[[1]][2]}
            dati=mod(dati1,input$norm)
            df1=dati[dati$drug1!=0 & dati$drug2==0,c("drug1","fa")]
            df2=dati[dati$drug1==0 & dati$drug2!=0,c("drug2","fa")]
            mod_df1=s_drug(df1)
            mod_df2=s_drug(df2)
            ci=ci(dati,mod_df1,mod_df2)
            pdf(file)
            print(plot_sdrug(mod_df1,drug1))
            print(plot_sdrug(mod_df2,drug2))
            print(plot_r2(dati1,input$norm,c(drug1,drug2))[[1]])
            print(isobol_plot(ci[["Log-logistic"]],drug=c(drug1,drug2),
                              model="Log-logistic"))
            heat_plot(ci[["Log-logistic"]],drug=c(drug1,drug2),model="Log-logistic")
            print(isobol_plot(ci[["Log-logistic [0]"]],drug=c(drug1,drug2),
                              model="Log-logistic [0]"))
            heat_plot(ci[["Log-logistic [0]"]],drug=c(drug1,drug2),
                      model="Log-logistic [0]")
            print(isobol_plot(ci[["Log-logistic [01]"]],drug=c(drug1,drug2),
                              model="Log-logistic [01]"))
            heat_plot(ci[["Log-logistic [01]"]],drug=c(drug1,drug2),
                      model="Log-logistic [01]")
            print(isobol_plot(ci[["Log-logistic [01] v2"]],drug=c(drug1,drug2),
                              model="Log-logistic [01] v2"))
            heat_plot(ci[["Log-logistic [01] v2"]],drug=c(drug1,drug2),
                      model="Log-logistic [01] v2")
            print(isobol_plot(ci[["Log-logistic [1]"]],drug=c(drug1,drug2),
                              model="Log-logistic [1]"))
            heat_plot(ci[["Log-logistic [1]"]],drug=c(drug1,drug2),
                      model="Log-logistic [1]")
            print(isobol_plot(ci[["Median-effect"]],drug=c(drug1,drug2),
                              model="Median-effect"))
            heat_plot(ci[["Median-effect"]],drug=c(drug1,drug2),model="Median-effect")
            dev.off()
        })
    
    output$formula <- renderUI({
        withMathJax(
            if(input$drug=="Log-logistic") {
                helpText('$$D=D_m \\cdot \\left(\\frac{max-min}{f_a-min} -1
                         \\right)^\\frac1 m$$')},
            if(input$drug=="Log-logistic [0]") {
                helpText('$$D=D_m \\cdot \\left(\\frac{max-f_a}{f_a}
                         \\right)^\\frac1 m$$')},
            if(input$drug=="Log-logistic [01]") {
                helpText('$$D=D_m \\cdot \\left(\\frac{1-f_a}{f_a}
                         \\right)^\\frac1 m$$')},
            if(input$drug=="Log-logistic [01] v2") {
                helpText('$$D=D_m \\cdot \\left(\\frac{1-f_a}{f_a}
                         \\right)^\\frac1 m$$')},
            if(input$drug=="Log-logistic [1]") {
                helpText('$$D=D_m \\cdot \\left(\\frac{1-min}{f_a-min} -1
                         \\right)^\\frac1 m$$')},
            if(input$drug=="Median-effect") {
                helpText('$$D=D_m \\cdot \\left(\\frac{f_a}{1-f_a}
                         \\right)^\\frac1 m$$')}
        )
    })

    output$formula1 <- renderUI({
        withMathJax(
            if(input$drug=="Log-logistic") {
                helpText('$$f_a=min+ \\frac{max-min}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')},
            if(input$drug=="Log-logistic [0]") {
                helpText('$$f_a=\\frac{max}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')},
            if(input$drug=="Log-logistic [01]") {
                helpText('$$f_a=\\frac{1}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')},
            if(input$drug=="Log-logistic [01] v2") {
                helpText('$$f_a=\\frac{1}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')},
            if(input$drug=="Log-logistic [1]") {
                helpText('$$f_a=min+ \\frac{1-min}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')},
            if(input$drug=="Median-effect") {
                helpText('$$f_a=\\frac{1}{1+\\left(\\frac{D}{D_m}
                         \\right)^m}$$')}
        )
    })
    
})