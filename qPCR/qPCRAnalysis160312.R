
# Read in qpcr data, ``Cp'' values as spit out by Roche480 fit to max
# 2nd deriv, reformatted to make a nice dataframe
# q means glutamine pulse, w means just water
exp86gap1q  <- data.frame(read.csv("exp86rawdata/141105_darach_exp86gap1q.csv"),
		plate="141105")
exp86gap1w  <- data.frame(read.csv("exp86rawdata/141113_darach_exp86gap1w.csv"),
		plate="141113")
exp86dip5wq <- data.frame(read.csv("exp86rawdata/141129_darach_exp86dip5wq.csv"),
		plate="141129")
exp86 <- rbind(exp86gap1q, exp86gap1w, exp86dip5wq)

# Just to check sanity and data format
ggplot(exp86, aes(x=minutes, y=cp, col=primers:treat))+theme_bw()+
	geom_point()+theme(legend.position="bottom")

sum(is.na(exp86))

# Use 100% efficiency assumption, I suppose
efficez <- c(1,1,1,1)
names(efficez) <- c("1200","gap1","hta1","dip5")

# Calc concentration relative to cycles at threshold of max 2nd deriv
exp86$conc <- (1 + efficez[as.character(exp86$primers)] ) ^ -exp86$cp
ggplot(exp86, aes(x=minutes, y=conc, col=primers:treat))+theme_bw()+
	geom_point()+theme(legend.position="bottom")

# How many measurements, variance, is all good?
countz <- aggregate(conc~minutes:primers:treat:plate,
  data=exp86,length)
varz <- aggregate(conc~minutes:primers:treat:plate,
  data=exp86,var,na.action=na.pass)

# Remove that column where I loaded sample 5 into column 6 on 
# accident, except H6 because I realized it at the end of the 
# column and paniced and I didn't write what I put in that well, 
# but it's HTA1 so no need to worry about that right now
# Also screwed up H12 on 141113 and didn't have another well to put 
# it into, so drop that one
exp86 <- subset(exp86, !(plate == "141105" & grepl("H6", well)))
exp86 <- subset(exp86, !(plate == "141113" & grepl("E12", well))) 

# Ok good, let's pick the timepoints after the chase starts so we're 
# actually fitting models to the chase

exp86 <- subset(exp86, minutes > 5) 

# Calculate means and standard deviations for each measurement
meanz <- aggregate(conc~primers:minutes:treat, data=exp86, mean, na.action=na.pass)
sdz   <- aggregate(conc~primers:minutes:treat, data=exp86, sd  , na.action=na.pass)

# For each gene, subset those means and divide over the spikein means
rm(ratioz)
for (genez in c("gap1","dip5","hta1")) {
	for (treatz in unique(meanz$treat)) {
		# Take a subset dataframe to work on
		tmpmeanz <- meanz[meanz$treat==treatz&meanz$primers%in%c(genez,"1200"),]
		tmpsdz <- sdz[meanz$treat==treatz&meanz$primers%in%c(genez,"1200"),]
		# Make a new temporary data frame
		tmpratioz <- data.frame(
			gene    = genez, 
			treat   = treatz,
			minutes = subset(tmpmeanz, primers == genez)$minutes, 
			mean    = subset(tmpmeanz, primers == genez)$conc) 
		# Divide the means of gene by spikein
		tmpratioz$ratiomean <- tmpratioz$mean / 
			tmpmeanz[tmpmeanz$minutes%in%tmpratioz$minutes&
				tmpmeanz$primers=="1200","conc"]
		# Calculate the standard deviation divided by mean, of the ratios
		tmpratioz$ratiocv <- sqrt(
 			(	tmpsdz[tmpsdz$minutes%in%tmpratioz$minutes&
					tmpsdz$primers==genez,"conc"] /
				tmpmeanz[tmpmeanz$minutes%in%tmpratioz$minutes&
					tmpmeanz$primers==genez,"conc"]) ^2 +
 			(	tmpsdz[tmpsdz$minutes%in%tmpratioz$minutes&
					tmpsdz$primers=="1200","conc"] /
				tmpmeanz[tmpmeanz$minutes%in%tmpratioz$minutes&
					tmpmeanz$primers=="1200","conc"]) ^2 )
		# Divide the ratio mean by the first timepoint
		tmpratioz$ratiomean <- tmpratioz$ratiomean / 
				subset(tmpratioz,minutes == min(minutes))$ratiomean
		# Multiply the cv by mean to get the standard deviations
		tmpratioz$ratiosd <- tmpratioz$ratiocv * tmpratioz$ratiomean
		if (exists("ratioz")) {
			ratioz <- rbind(ratioz,tmpratioz)
		} else {
			ratioz <- tmpratioz
		}
	}
}
ratioz$gene <- factor(toupper(ratioz$gene),levels=c("GAP1","DIP5","HTA1"))

# Have I added the treatment (water or glutamine) yet?
ratioz$post <- ratioz$minutes > 13

# How about a linear model of the log mean ratio predicted by the minutes, and a term that takes into account the treatment and if it's been applied yet?
# Taking out measurements after 25 minutes in the glutamine condition because they're bottomed out and including those in a linear model is just going to over fit to noise around zero mRNA
summary(lm(data=subset(ratioz,gene=="GAP1"&!(minutes>25&treat=="q")),log(ratiomean)~minutes:treat:post+minutes)) 

# It likes a shift in rates. What's it look like, but also using stat_smooth to fit the same model as above?

# Just gap1
#pdf("150211exp86gap1.pdf",width=7,height=5)
ggplot(subset(ratioz,gene=="GAP1"), aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(cex=2)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	ylab("log ( mean ratio of GAP1 mRNA to spike-in ),\ndivided by first timepoint")+
	xlab("Minutes after label chase started")+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
		x=12.6,y=-2.5,angle=90,size=4)+
	scale_color_manual(name="Glutamine upshift at 13min?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat),
		data=subset(ratioz,!(minutes>25&treat=="q")),
		se=F,linetype="dotted")
dev.off()

# And the same for the DIP5 data. Only did 5 to 25 minutes, so no subsetting needed for the model

summary(lm(data=subset(ratioz,gene=="DIP5"),log(ratiomean)~minutes:treat:post+minutes)) 

#pdf("150211exp86dip5.pdf",width=7,height=5)
ggplot(subset(ratioz,gene=="DIP5"), aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(cex=2)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	ylab("log ( mean ratio of DIP5 mRNA to spike-in ),\ndivided by first timepoint")+
	xlab("Minutes after label chase started")+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
		x=12.6,y=-2.5,angle=90,size=4)+
	scale_color_manual(name="Glutamine upshift at 13min?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat),
		data=subset(ratioz,!(minutes>25&treat=="q")),
		se=F,linetype="dotted")
dev.off()

# Las dos

summary(lm(data=subset(ratioz,!(gene=="GAP1"&minutes>25&treat=="q")),log(ratiomean)~minutes:treat:post:gene+minutes:gene))

#pdf("150312exp86water.pdf",width=7,height=5)
ggplot(subset(ratioz), aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat,alpha=as.numeric(treat=="w")
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(cex=2)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene)+
	ylab("log ratio of mean mRNA to mean spike-in")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	scale_alpha_continuous(range=c(0,1))+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	#geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
	#	x=12,y=-2.5,angle=90,size=4)+
	#stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
	#	data=subset(ratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
	#	se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()

#pdf("150720exp86bothwater.pdf",width=8,height=6)
plotratioz <- subset(ratioz,gene!="HTA1"&minutes<25)
ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat,alpha=as.numeric(treat%in%c("w"))
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	#coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene)+
	ylab("log ratio of mean mRNA to mean spike-in\n")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	#geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
	#	x=12,y=-2.5,angle=90,size=4)+
	#stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
#		data=subset(plotratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
#		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()
#pdf("150719exp86bothextended.pdf",width=8,height=4)
plotratioz <- subset(ratioz,gene!="HTA1")
ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	#coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene,scales="free_x")+
	ylab("log ratio of mean mRNA to mean spike-in\n")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	#geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
	#	x=12,y=-2.5,angle=90,size=4)+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
		data=subset(plotratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()



summary(lm(data=plotratioz,log(ratiomean)~minutes:treat:post:gene+minutes:gene))





####

pdf("150819exp86_gap1.pdf",width=8,height=8)
plotratioz <- subset(ratioz,gene=="GAP1"&minutes<25)
ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat#,alpha=as.numeric(treat%in%c("w","q"))
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	ylab("log ratio of mean mRNA to mean spike-in\n")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	ggtitle("GAP1 mRNA")+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat),
		data=subset(plotratioz),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()

pdf("150819exp86_dip5.pdf",width=8,height=8)
plotratioz <- subset(ratioz,gene=="DIP5"&minutes<25)
ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat#,alpha=as.numeric(treat%in%c("w","q"))
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	ylab("log ratio of mean mRNA to mean spike-in\n")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	ggtitle("DIP5 mRNA")+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat),
		data=subset(plotratioz),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()






####

pdf("150903exp86both.pdf",width=5,height=3.5)
plotratioz <- subset(ratioz,gene!="HTA1"&minutes<25)
ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat#,alpha=as.numeric(treat%in%c("w","q"))
		))+theme_bw()+
	geom_vline(xintercept=13,col="black",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	#coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene)+
	ylab("log ( labeled transcript)")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	#geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
	#	x=12,y=-2.5,angle=90,size=4)+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
		data=subset(plotratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
dev.off()


plotratioz <- subset(ratioz,gene!="HTA1"&minutes<25)
g <- ggplot(plotratioz, aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat#,alpha=as.numeric(treat%in%c("w","q"))
		))+theme_dm(base_family="Helvetica")+
	geom_vline(xintercept=13,col="grey50",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	#coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene)+
	ylab("log ( labeled transcript )")+
	xlab("Minutes after label chase started")+
	guides(alpha=F)+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	#geom_text(label="Added water +/- glutamine",col="black",alpha=.2,
	#	x=12,y=-2.5,angle=90,size=4)+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
		data=subset(plotratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift at 13min (dashed line)?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
g
ggsave("150929exp86both.pdf",g,width=8,height=6)
ggsave("150929exp86bothhalf.pdf",g,width=4,height=4)


plotratioz <- subset(ratioz,gene!="HTA1")
g <- ggplot(plotratioz)+
	aes(x=minutes,y=log(ratiomean),
		ymax=log(ratiomean+ratiosd),ymin=log(ratiomean-ratiosd),
		col=treat)+
	theme_dmprez(base_family="Helvetica")+
	geom_vline(xintercept=13,col="grey50",linetype="dashed",alpha=.5)+
	geom_point(size=1.5)+geom_errorbar(width=.3)+
	theme(legend.position="bottom")+
	#coord_cartesian(ylim=c(-4.5,.5))+
	facet_grid(.~gene)+
	ylab("log ( labeled transcript )")+
	xlab("Minutes after label chase started")+
	#guides(alpha=F)+
	#ggtitle("Abundance of labeled GAP1 mRNA during chase,\nwith glutamine upshift")+
	geom_text(label="Added water +/- glutamine",col="white",alpha=.2,
		x=12,y=-2.5,angle=90,size=4)+
	stat_smooth(method="lm",formula="y~x",aes(group=factor(post):treat:gene),
		data=subset(plotratioz,!(gene=="GAP1"&minutes>25&treat=="q")),
		se=F,linetype="dotted")+
	scale_color_manual(name="Glutamine upshift?",
		labels=c("400uM","(blank)"),
		breaks=c("q","w"),
		values=c("#D55E00","#56B4E9"))
g
ggsave("151112exp86both.pdf",g,width=8,height=6)

####


plotratioz <- subset(ratioz,gene!="HTA1")
g <- ggplot(subset(plotratioz,gene=="GAP1"&treat=="w"))+
  aes(x=minutes,y=log(ratiomean),
    ymax=log(ratiomean+ratiosd),
    ymin=log(ratiomean-ratiosd),
    col=treat)+
  theme_wb()+
  geom_vline(xintercept=13,col="white",
    linetype="dashed",alpha=.5)+
  geom_point(size=1.5)+geom_errorbar(width=.3)+
  theme(legend.position="bottom")+
  coord_cartesian(ylim=c(-4.5,.5))+
  ylab("log ( labeled transcript )")+
  xlab("Minutes after label chase started")+
  stat_smooth(method="lm",formula="y~x",
    aes(group=factor(post):treat:gene),
    data=subset(plotratioz,gene=="GAP1"&
      treat=="w"),
    se=F,linetype="dotted")+
  scale_color_manual(name="Glutamine upshift?",
    labels=c("400uM","(blank)"),
    breaks=c("q","w"),
    values=c("#56B4E9","#D55E00"))
g
ggsave("/home/zed/lab/outz/figz/160211exp86water.png",
  g,width=8,height=6)

plotratioz <- subset(ratioz,gene!="HTA1")
g <- ggplot(subset(plotratioz,gene=="GAP1"))+
  aes(x=minutes,y=log(ratiomean),
    ymax=log(ratiomean+ratiosd),
    ymin=log(ratiomean-ratiosd),
    col=treat)+
  theme_wb()+
  geom_vline(xintercept=13,col="white",
    linetype="dashed",alpha=.5)+
  geom_point(size=1.5)+geom_errorbar(width=.3)+
  theme(legend.position="bottom")+
  coord_cartesian(ylim=c(-4.5,.5))+
  ylab("log ( labeled transcript )")+
  xlab("Minutes after label chase started")+
  stat_smooth(method="lm",formula="y~x",
    aes(group=factor(post):treat:gene),
    data=subset(plotratioz,gene=="GAP1"&
      !(treat=="q"&minutes>25)),
    se=F,linetype="dotted")+
  scale_color_manual(name="Glutamine upshift?",
    labels=c("400uM","(blank)"),
    breaks=c("q","w"),
    values=c("#D55E00","#56B4E9"))
g
ggsave("/home/zed/lab/outz/figz/160211exp86.png",
  g,width=8,height=6)

