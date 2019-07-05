options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/aso_survival')
}
library(survival)
library(sqldf)
library(drc)
library(reshape2)
library(beeswarm)
library(raster)
library(rgdal)

imgmode = 'png'
# imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

expand_range = function(x,by=0.5) {
  xmin = min(x,na.rm=T) - by
  xmax = max(x,na.rm=T) + by
  return (c(xmin, xmax))
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

plot_params = read.table('data/plot_params.tsv',sep='\t',header=T,comment.char='')

mouse_survival = read.table('data/mouse_survival.tsv',sep='\t',header=T)
mouse_survival$treatment = mouse_survival$cohort # alias for this variable



imgsave(paste('figures/script_generated/figure-1.',imgmode,sep=''),width=9*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,1,4,4,2,2,4,4,3,3,4,4),nrow=3,byrow=T)
layout(layout_matrix)

potency = read.table('data/potency.tsv',sep='\t',header=T)
potency$color = plot_params$color[match(potency$treatment, plot_params$treatment)]

psmry = sqldf("
              select   x, region, color, treatment, dose, avg(mrna) mean, stdev(mrna) sd, count(*) n
              from     potency
              group by 1, 2, 3, 4, 5
              order by 1, 2, 3, 4, 5
              ;")
psmry$l95 = psmry$mean - 1.96*psmry$sd/sqrt(psmry$n)
psmry$u95 = psmry$mean + 1.96*psmry$sd/sqrt(psmry$n)

psmry$color =  alpha(psmry$color, psmry$dose/max(psmry$dose,na.rm=T))
psmry$color[is.na(psmry$dose)] = plot_params$color[match(psmry$treatment[is.na(psmry$dose)], plot_params$treatment)]

txlabs = sqldf("
select   region, treatment, color, avg(x) xmid
from     potency
group by 1, 2, 3
order by 1, 2, 3
;")

reglabs = sqldf("
select   region, avg(x) xmid
from     potency
group by 1
order by 1
;")

par(mar=c(6,5,4,1))
xlims = c(min(potency$x)-1, max(potency$x) + 0.5)
plot(NA,NA,xlim=xlims,ylim=c(0,1.2),xaxs='i',yaxs='i',ann=FALSE,axes=FALSE)
points(psmry$x, psmry$mean, type='h', lwd=20, lend=1, col=psmry$color)
arrows(x0=psmry$x, y0=psmry$l95, y1=psmry$u95, angle=90, length=0.05, code=3, lwd=1, col='black')
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1, las=2)
abline(h=0,lwd=2)
abline(h=1, lwd=1, lty=3)
abline(v=min(xlims), lwd=2)
mtext(side=1, at=psmry$x, text=psmry$dose, line=0, col=psmry$color, cex=0.7, font=1)
par(xpd=T)
text(x=txlabs$xmid, y=rep(-.15,nrow(txlabs)), labels=txlabs$treatment, srt=45, adj=c(1,1), col=txlabs$color, font=1)
par(xpd=F)
mtext(side=3, at=reglabs$xmid, text=reglabs$region, col='black', cex=1, font=1, line=0)
mtext(side=2, line=3, text='relative mRNA')
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)


par(mar=c(4,5,4,1))

# blank panel for the Western blot to be added in PowerPoint
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

washout = read.table('data/washout.tsv',sep='\t',header=T)
washout$color = plot_params$color[match(washout$treatment, plot_params$treatment)]

wsmry = sqldf("
              select   weeks, region, color, treatment, avg(mrna) mean, stdev(mrna) sd, count(*) n
              from     washout
              group by 1, 2, 3, 4
              order by 1, 2, 3, 4
;")
wsmry$l95 = wsmry$mean - 1.96*wsmry$sd/sqrt(wsmry$n)
wsmry$u95 = wsmry$mean + 1.96*wsmry$sd/sqrt(wsmry$n)

par(mar=c(4,5,4,1))
xlims=c(0,20)
plot(NA, NA, xlim=xlims, ylim=c(0,1.2), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
subs1 = subset(wsmry, region=='cortex')
pp = subset(plot_params, !(treatment %in% c('uninfected','no treatment')))
legend(x=13,y=.7,legend=pp$treatment,lwd=1.5,pch=19,col=pp$color,text.col=pp$color,text.font=1,bty='n')
axis(side=1, at=(0:4)*4, lwd=0, lwd.ticks=1)
abline(h=0,lwd=2)
abline(h=1, lwd=1, lty=3)
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1, las=2)
abline(v=min(xlims), lwd=2)
mtext(side=1, line=2.5, text='weeks post-dose')
mtext(side=2, line=3, text='relative mRNA')
for (tx in subs1$treatment) {
  subs2 = subset(subs1, treatment==tx)
  points(subs2$weeks, subs2$mean, pch=19, type='b', lwd=1.5, col=subs2$color)
  arrows(x0=subs2$weeks, y0=subs2$l95, y1=subs2$u95, lwd=1, col=subs2$color, angle=90, code=3, length=0.05)
}
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

par(mar=c(4,5,4,1))
# empty plot for panel D to be inserted
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
mtext('D', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off()





imgsave(paste('figures/script_generated/figure-2.',imgmode,sep=''), width=10*resx, height=6*resx, res=resx)

layout_matrix = matrix(c(1,1,1,2,2,3,3,3,4,4,5,5,5,6,6,6),nrow=2,byrow=T)
layout(layout_matrix)

# read in NIH mouse data
nihp = subset(mouse_survival, expt == 'NIH prophylactic')

sf = survfit(Surv(nihp$onset_dpi, nihp$onset) ~ nihp$treatment)
sf$treatment = gsub('nihp\\$treatment=','',names(sf$strata))
sf$color = plot_params$color[match(sf$treatment, plot_params$treatment)]

par(mar=c(4,5,4,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-25,550), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=c(-1:11)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.01)
axis(side=1, at=c(0:5)*100, lwd=0, lwd.ticks=1, tck=-0.03)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(h=0)
abline(v=0)
par(xpd=TRUE)
text(x=c(-14,46,106),y=c(1.05,1.05,1.05),pos=3,labels=c('','300 µg doses',''))
points(x=c(-14,46,106),y=c(1.05,1.05,1.05),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='proportion asymptomatic')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='onset')

mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)

# onset delayed by 82 - 99% on median:
median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='active ASO 1']) / median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='saline'])
median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='active ASO 2']) / median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='saline'])

nihp$x = plot_params$x[match(nihp$treatment, plot_params$treatment)]
nihp$color = plot_params$color[match(nihp$treatment, plot_params$treatment)]
set.seed(1)
nihp$xjit = jitter(nihp$x,.5)

nihp_smry = sqldf("
select   x, color, avg(duration) mean, stdev(duration) sd, sum(case when duration is not null then 1 else 0 end) n
from     nihp
group by 1, 2
order by 1
")
nihp_smry$l95 = nihp_smry$mean - 1.96 * nihp_smry$sd / sqrt(nihp_smry$n)
nihp_smry$u95 = nihp_smry$mean + 1.96 * nihp_smry$sd / sqrt(nihp_smry$n)

par(mar=c(4,5,4,1))
plot(NA,NA,xlim=c(0.5,5.5),ylim=c(0,75),axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
axis(side=1, at=c(0,6), labels=NA, lwd=1, lwd.ticks=0)
par(xpd=T)
subs = subset(plot_params, treatment != 'uninfected')
text(y=rep(0,5), x=subs$x, adj=c(1,1), labels=subs$treatment, col=subs$color, font=1, srt=45)
par(xpd=F)
axis(side=2, at=(0:3)*25, las=2)
mtext(side=2, line=2.5, text='progression (days)')
# points(nihp$xjit, nihp$duration, pch=19, col=nihp$color)
# beeswarm looks better than jitter:
# beeswarm(duration ~ x, data=nihp, pwcol=color, pch=19, method='swarm', add=T)
# for this plot, beeswarm() works since all groups are represented, however, for consistency
# with Fig 3B we are plotting this way - see comment at 3B
for (treatment in unique(nihp$treatment[!is.na(nihp$duration)])) {
  rows = nihp$treatment == treatment & !is.na(nihp$duration)
  nihp$xbee[rows] = nihp$x[rows] + beeswarm(duration ~ x, data=nihp[rows,], do.plot=F)$x - 1
}
points(nihp$xbee, nihp$duration, pch=19, col=nihp$color)
# means & error bars per JCII requirements
bar_color = '#000000'
# mean bar
segments(x0=nihp_smry$x-.25, x1=nihp_smry$x+.25, y0=nihp_smry$mean, col=bar_color)
# l95 & u95 bars
segments(x0=nihp_smry$x-.125, x1=nihp_smry$x+.125, y0=nihp_smry$l95, col=bar_color)
segments(x0=nihp_smry$x-.125, x1=nihp_smry$x+.125, y0=nihp_smry$u95, col=bar_color)
# connect them all
segments(x0=nihp_smry$x, y0=nihp_smry$l95, y1=nihp_smry$u95, col=bar_color)

mtext(side=3, line=1.5, text='progression')

mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

sf = survfit(Surv(nihp$dpi, nihp$acm) ~ nihp$treatment)
sf$treatment = gsub('nihp\\$treatment=','',names(sf$strata))
sf$color = plot_params$color[match(sf$treatment, plot_params$treatment)]

par(mar=c(4,5,4,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-25,550), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=c(-1:11)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.01)
axis(side=1, at=c(0:5)*100, lwd=0, lwd.ticks=1, tck=-0.03)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)
par(xpd=TRUE)
text(x=c(-14,46,106),y=c(1.05,1.05,1.05),pos=3,labels=c('','300 µg doses',''))
points(x=c(-14,46,106),y=c(1.05,1.05,1.05),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='proportion alive (all cause)')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')

mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

# mortality delayed by 81 - 98% on median
median(nihp$dpi[nihp$acm==1 & nihp$treatment=='active ASO 1']) / median(nihp$dpi[nihp$acm==1 & nihp$treatment=='saline'])
median(nihp$dpi[nihp$acm==1 & nihp$treatment=='active ASO 2']) / median(nihp$dpi[nihp$acm==1 & nihp$treatment=='saline'])

# number of days that corresponds to: 116 and 140
median(nihp$dpi[nihp$acm==1 & nihp$treatment=='active ASO 1']) - median(nihp$dpi[nihp$acm==1 & nihp$treatment=='saline'])
median(nihp$dpi[nihp$acm==1 & nihp$treatment=='active ASO 2']) - median(nihp$dpi[nihp$acm==1 & nihp$treatment=='saline'])

# pseudopanel to hold master legend
par(mar=c(4,5,4,1))
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
legend(x=0, y=1, legend=plot_params$treatment, lwd=3, lty=plot_params$lty, col=plot_params$color, bty='n', text.col=plot_params$color, text.font=1, cex=1.4, title='cohorts', title.col='black')
rml_bounds = c(.28,.76)
axis(side=2, at=rml_bounds, labels=NA, line=1, tck=0.05)
mtext(side=2, at=mean(rml_bounds), line=1.5, text='prion-infected',cex=1)

# Broad data
wts = read.table('data/weights.tsv',sep='\t',header=T)
wts$color = plot_params$color[match(wts$treatment, plot_params$treatment)]

par(mar=c(4,5,4,1))
plot(NA, NA, xlim=c(120,315), ylim=c(14,35), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
for (treatment in plot_params$treatment) {
  subs = wts$treatment == treatment
  points(wts$dpi[subs], wts$mean[subs], col=wts$color[subs], pch=20, type='l', lty=plot_params$lty[plot_params$treatment==treatment], lwd=3)
  polygon(x=c(wts$dpi[subs],rev(wts$dpi[subs])), y=c(wts$u95[subs],rev(wts$l95[subs])), col=alpha(wts$color[subs],.35), border=NA)
}
axis(side=1, at=c(120,3:8*50), lwd=1, lwd.ticks=1)
axis(side=2, at=c(0:7)*5, lwd=1, lwd.ticks=1, las=2)
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=2, line=2.5, text='body weight (g)')
mtext(side=3, line=1.5, text='body weights')
par(xpd=TRUE)
par(xpd=FALSE)

mtext('D', side=3, cex=2, adj = -0.1, line = 0.5)

survival = subset(mouse_survival, expt=='Broad prophylactic')

# survival curves - all cause mortality
par(mar=c(4,5,4,1))
sf = survfit(Surv(survival$dpi, survival$acm) ~ survival$cohort)
sf$cohort = gsub('survival\\$cohort=','',names(sf$strata))
sf$color = plot_params$color[match(sf$cohort, plot_params$treatment)]
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-25,550), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=c(-1:11)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.01)
axis(side=1, at=c(0:5)*100, lwd=0, lwd.ticks=1, tck=-0.03)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
# abline(h=0)
abline(v=0)
par(xpd=TRUE)
text(x=mean(c(-14,76)),y=c(1.05),pos=3,labels='500 µg doses')
points(x=c(-14,76),y=c(1.05,1.05),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='proportion alive (all cause)')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')
par(xpd=TRUE)
#legend(x=max(dates$dpi)+10,y=1, legend=cohorts$alias, lwd=3, col=cohorts$color, bty='n', text.col=cohorts$color, text.font=1)
par(xpd=FALSE)

mtext('E', side=3, cex=2, adj = -0.1, line = 0.5)

# all-cause survival increased by how much
median(survival$dpi[survival$acm==1 & survival$cohort=='active ASO 1']) / median(survival$dpi[survival$acm==1 & survival$cohort=='saline'])
median(survival$dpi[survival$acm==1 & survival$cohort=='active ASO 2']) / median(survival$dpi[survival$acm==1 & survival$cohort=='saline'])

# number of days
median(survival$dpi[survival$acm==1 & survival$cohort=='active ASO 1']) - median(survival$dpi[survival$acm==1 & survival$cohort=='saline'])
median(survival$dpi[survival$acm==1 & survival$cohort=='active ASO 2']) - median(survival$dpi[survival$acm==1 & survival$cohort=='saline'])

dev.off() # end Figure 2


# how many standard deviations past the mean was the long-surviving control ASO mouse when it was culled?
mean_aso_c = mean(survival$dpi[survival$cause_of_death == 'prion disease' & survival$cohort=='control ASO' & survival$dpi < 300])
sd_aso_c = sd(survival$dpi[survival$cause_of_death == 'prion disease'  & survival$cohort=='control ASO' & survival$dpi < 300])
longest = max(survival$dpi[survival$cohort == 'control ASO'])
(longest - mean_aso_c) / sd_aso_c



# Figure 3
imgsave(paste('figures/script_generated/figure-3.',imgmode,sep=''), width=10*resx, height=3.5*resx, res=resx)

layout_matrix = matrix(c(1,1,1,2,2,3,3,3,3),nrow=1,byrow=T)
layout(layout_matrix)

# NIH 120 dpi data
nihp = subset(mouse_survival, expt == 'NIH 120 dpi')
nihp$treatment = nihp$cohort

sf = survfit(Surv(nihp$onset_dpi, nihp$onset) ~ nihp$treatment)
sf$treatment = gsub('nihp\\$treatment=','',names(sf$strata))
sf$color = plot_params$color[match(sf$treatment, plot_params$treatment)]

# in case some reviewer insists we add p values
survdiff(Surv(onset_dpi, onset) ~ treatment, data=subset(nihp, treatment %in% c('active ASO 1','saline'))) # p = 6e-5
survdiff(Surv(dpi, onset) ~ treatment, data=subset(nihp, acm==1 & treatment %in% c('active ASO 1','saline'))) # p = 1e-4

par(mar=c(5,5,4,1))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-14,300), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=c(-1:6)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.01)
axis(side=1, at=c(0:3)*100, lwd=0, lwd.ticks=1, tck=-0.03)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
# abline(h=0)
abline(v=0)
par(xpd=TRUE)
text(x=c(120),y=c(1.05),pos=3,labels='300 µg dose')
points(x=c(120),y=c(1.05),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='proportion asymptomatic')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='onset')

mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)

# percent by which onset delayed
median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='active ASO 1']) / median(nihp$onset_dpi[nihp$onset==1 & nihp$treatment=='saline'])

nihp$x = plot_params$x[match(nihp$treatment, plot_params$treatment)]
nihp$color = plot_params$color[match(nihp$treatment, plot_params$treatment)]
set.seed(1)
nihp$xjit = jitter(nihp$x,.5)

nihp_smry = sqldf("
select   x, color, avg(duration) mean, stdev(duration) sd, sum(case when duration is not null then 1 else 0 end) n
from     nihp
group by 1, 2
order by 1
")
nihp_smry$l95 = nihp_smry$mean - 1.96 * nihp_smry$sd / sqrt(nihp_smry$n)
nihp_smry$u95 = nihp_smry$mean + 1.96 * nihp_smry$sd / sqrt(nihp_smry$n)

par(mar=c(5,5,4,1))
xlims = c(1.5,4.5)
ylims = c(0,75)
plot(NA,NA,xlim=xlims,ylim=ylims,axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
axis(side=1, at=c(0,6), labels=NA, lwd=1, lwd.ticks=0)
#axis(side=1, at=plot_params$x, labels=plot_params$treatment, srt=45, lwd=0, lwd.ticks=0)
#mtext(side=1, line=0, at=plot_params$x, text=plot_params$treatment, srt=45, col=plot_params$color)
par(xpd=T)
subs = subset(plot_params, !(treatment %in% c('uninfected','no treatment','active ASO 2')))
text(y=rep(0,5), x=subs$x, adj=c(1,1), labels=subs$treatment, col=subs$color, font=1, srt=45)
par(xpd=F)
axis(side=2, at=(0:3)*25, las=2)
mtext(side=2, line=2.5, text='progression (days)')
mtext(side=3, line=1.5, text='progression')
# beeswarm(duration ~ x, data=nihp, pwcol=color, pch=19, method='swarm', at=2:4, add=T, axes=F, ann=F, glim=xlims, dlim=ylims)
# if you directly use beeswarm to plot, it assigns its own "x" values based on the number of groups
# with non-NA data, which doesn't match our own "x" values. you can override this with the "at" param
# but the docs make no guarantee that the ordering of the "at" will match the ordering of your own "x"
# so to be absolutely safe that the right data are plotted in the right column, need to do this:
for (treatment in unique(nihp$treatment[!is.na(nihp$duration)])) {
  rows = nihp$treatment == treatment & !is.na(nihp$duration)
  nihp$xbee[rows] = nihp$x[rows] + beeswarm(duration ~ x, data=nihp[rows,], do.plot=F)$x - 1
}
points(nihp$xbee, nihp$duration, pch=19, col=nihp$color)
# means & error bars per JCII requirements
bar_color = '#000000'
# mean bar
segments(x0=nihp_smry$x-.25, x1=nihp_smry$x+.25, y0=nihp_smry$mean, col=bar_color)
# l95 & u95 bars
segments(x0=nihp_smry$x-.125, x1=nihp_smry$x+.125, y0=nihp_smry$l95, col=bar_color)
segments(x0=nihp_smry$x-.125, x1=nihp_smry$x+.125, y0=nihp_smry$u95, col=bar_color)
# connect them all
segments(x0=nihp_smry$x, y0=nihp_smry$l95, y1=nihp_smry$u95, col=bar_color)

mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

sf = survfit(Surv(nihp$dpi, nihp$acm) ~ nihp$treatment)
sf$treatment = gsub('nihp\\$treatment=','',names(sf$strata))
sf$color = plot_params$color[match(sf$treatment, plot_params$treatment)]

par(mar=c(5,5,4,10))
plot(sf, ann=FALSE, axes=FALSE, xlim=c(-14,300), ylim=c(0,1.05), xaxs='i', yaxs='i', col=sf$color, lwd=c(3,3,3,3,3))
axis(side=1, at=c(-1:6)*50, labels=NA, lwd=1, lwd.ticks=1, tck=-0.01)
axis(side=1, at=c(0:3)*100, lwd=0, lwd.ticks=1, tck=-0.03)
axis(side=2, at=c(0:4)/4, labels=percent(0:4/4), lwd=0, lwd.ticks=1, las=2)
abline(v=0)
par(xpd=TRUE)
text(x=c(120),y=c(1.05),pos=3,labels='300 µg dose')
points(x=c(120),y=c(1.05),pch=25,col='black',bg='black')
par(xpd=FALSE)
mtext(side=2, line=3.5, text='proportion alive (all cause)')
mtext(side=1, line=2.5, text='days post-infection')
mtext(side=3, line=1.5, text='mortality')

subs = subset(plot_params, !(treatment %in% c('uninfected','no treatment')))
par(xpd=T)
legend(x=270, y=1, legend=subs$treatment, lwd=3, lty=subs$lty, col=subs$color, bty='n', text.col=subs$color, text.font=1, cex=1.4, title='cohorts', title.col='black')
par(xpd=F)

mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

median(nihp$dpi[nihp$acm==1 & nihp$treatment=='active ASO 1']) / median(nihp$dpi[nihp$acm==1 & nihp$treatment=='saline'])

dev.off() # end Figure 3





tolerab = read.table('data/tolerability.tsv',sep='\t',header=T,comment.char='',quote='')
tolerab$color = plot_params$color[match(tolerab$treatment,plot_params$treatment)]

screening = read.table('data/screening.tsv',sep='\t',header=T)
screening$color = plot_params$color[match(screening$treatment,plot_params$treatment)]


imgsave(paste('figures/script_generated/figure-s1.',imgmode,sep=''),width=9*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(rep(1,10),2,2,2,3,3,3,4,4,4,4,2,2,2,3,3,3,4,4,4,4),nrow=3,byrow=T)
layout(layout_matrix)

# gene diagram
# active aso 1 - CTATTTAATGTCAGTCT - GRCm38 mmu2 21901020 to 21901036 and 131938400 to 131938416
# active aso 2 - TTGCAATTCTATCCAAA - GRCm38 mmu2 119699923 to 119699939 and 131914306 to 131914322
# from UCSC table browser GRCm38: uc008mly.3	chr2	+	131909927	131938436	131936429	131937194	3	131909927,131912194,131936419,	131910003,131912292,131938436,	P04925	uc008mly.3

active_aso_1 = c(131938400, 131938416)
active_aso_2 = c(131914306, 131914322)

mouse_exon1_start = 131909927 # transcription start site
mouse_exon1_end = 131910003

mouse_exon2_start = 131912194
mouse_exon2_end = 131912292

mouse_exon3_start = 131936419
mouse_exon3_end = 131938436 # transcription end site

mouse_cds_start = 131936429
mouse_cds_end = 131937194

# make plots

ymin = -1.7
ymax = -0.5
yaso = -0.6
xmin = mouse_exon1_start - 5000
xmax = mouse_exon3_end + 5000
transcript = data.frame(starts = c(mouse_exon1_start),
                        ends = c(mouse_exon3_end))
exons = data.frame(starts = c(mouse_exon1_start, mouse_exon2_start, mouse_exon3_start),
                   ends = c(mouse_exon1_end, mouse_exon2_end, mouse_exon3_end))
cds = data.frame(starts = c(mouse_cds_start),
                 ends = c(mouse_cds_end))
# reference_length = mouse_exon3_end - mouse_exon1_start
# reference_length_kb = paste(formatC(reference_length/1000, format='f', digits=1), 'kb')

reference_y = -1.0

intron_lwd = 1
utr_lwd = 10
cds_lwd = 20

big_ticks = seq(131900000, 131950000, by=10000)
big_tick_labels = formatC(big_ticks, format='fg', big.mark=',')
small_ticks = seq(131900000, 131950000, by=1000)

par(mar=c(5,5,3,2))
plot(NA, NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
arrows(x0=xmin+600, x1=xmin+1600, y0=ymin + 0.15, y1=ymin + 0.15, length=0.025, code=3, angle=90)
text(x=xmin+1150, y=ymin+0.1, pos=3, labels='1 kb')
axis(side=1, at=small_ticks, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=big_ticks, labels=big_tick_labels, lwd=1, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='GRCm38 genomic coordinates')
# reference sequence
segments(y0=reference_y, x0=transcript$starts, x1=transcript$ends, lwd=intron_lwd, lend=1)
segments(y0=reference_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=reference_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
# label gene name
text(x=mouse_exon1_start, y=reference_y, pos=2, font=3, labels='Prnp')

# plot the ASOs
par(xpd=T)
arrows(x0=mean(active_aso_1), y0=yaso, y1=reference_y+0.2, code=2, angle=90, length=0.02, col=plot_params$color[plot_params$treatment=='active ASO 1'], lwd=2)
text(x=mean(active_aso_1), y=yaso, pos=3, labels='active ASO 1', col=plot_params$color[plot_params$treatment=='active ASO 1'], font=1)
arrows(x0=mean(active_aso_2), y0=yaso, y1=reference_y+0.2, code=2, angle=90, length=0.02, col=plot_params$color[plot_params$treatment=='active ASO 2'], lwd=2)
text(x=mean(active_aso_2), y=yaso, pos=3, labels='active ASO 2', col=plot_params$color[plot_params$treatment=='active ASO 2'], font=1)
par(xpd=F)
mtext('A', side=3, cex=2, adj = -0.025, line = 0.5)

ssmry = sqldf("
              select   x, region, color, treatment, avg(mrna) mean, stdev(mrna) sd, count(*) n
              from     screening
              group by 1, 2, 3, 4
              order by 1, 2, 3, 4
              ;")
ssmry$l95 = ssmry$mean - 1.96*ssmry$sd/sqrt(ssmry$n)
ssmry$u95 = ssmry$mean + 1.96*ssmry$sd/sqrt(ssmry$n)


txlabs = sqldf("
               select   s.region, p.treatment, s.color, avg(s.x) xmid
               from     screening s, plot_params p
               where    s.treatment = p.treatment
               group by 1, 2, 3
               order by 1, 2, 3
               ;")

reglabs = sqldf("
                select   region, avg(x) xmid
                from     screening
                group by 1
                order by 1
                ;")

par(mar=c(5,5,3,2))
xlims = c(min(screening$x)-1, max(screening$x) + 0.5)
plot(NA,NA,xlim=xlims,ylim=c(0,1.2),xaxs='i',yaxs='i',ann=FALSE,axes=FALSE)
points(ssmry$x, ssmry$mean, type='h', lwd=20, lend=1, col=ssmry$color)
arrows(x0=ssmry$x, y0=ssmry$l95, y1=ssmry$u95, angle=90, length=0.05, code=3, lwd=1, col='black')
axis(side=2, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1, las=2)
abline(h=0,lwd=2)
abline(h=1, lwd=1, lty=3)
abline(v=min(xlims), lwd=2)
mtext(side=1, at=ssmry$x, text=ssmry$dose, line=0, col=ssmry$color, cex=0.7, font=1)
par(xpd=T)
text(x=txlabs$xmid, y=rep(0,nrow(txlabs)), labels=txlabs$display, srt=45, adj=c(1,1), col=txlabs$color, font=1)
par(xpd=F)
mtext(side=3, at=reglabs$xmid, text=reglabs$region, col='black', cex=1, font=1, line=0)
mtext(side=2, line=3, text='relative mRNA')
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

par(mar=c(5,5,3,2))
plot(NA,NA,xlim=c(0,8),ylim=c(0,27.5),axes=FALSE,ann=FALSE,xaxs='i',yaxs='i')
baseline_colno = which(colnames(tolerab)=='bw0wk')
#bwmat = as.matrix(tolerab[,baseline_colno:(baseline_colno+8)])
#bwmat = bwmat / bwmat[,1]
for (i in 1:nrow(tolerab)) {
  points(0:8,tolerab[i,baseline_colno:(baseline_colno+8)],type='l',lwd=2,col=tolerab$color[i])
}
axis(side=1,at=0:8)
mtext(side=1,at=c(0,8),text=c('surgery','euthanasia'),line=2)
mtext(side=1,text='weeks',line=3,font=1)
axis(side=2,at=(0:6)*5,las=2)
mtext(side=2,line=2.5,text='body weight (g)')
pp = subset(plot_params, treatment %in% tolerab$treatment)
legend(x=3,y=15,legend=pp$treatment,lwd=2,col=pp$color,text.col=pp$color,text.font=1,bty='n')
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

markers = melt(tolerab[,c('treatment','cortex_iba1_pct','cord_iba1_pct','cord_cd68_pct')])
markers$value = markers$value / 100
markers$x = rep(c(1:3,5:7,9:11),each=4)
markers$color = plot_params$color[match(markers$treatment,plot_params$treatment)]

msmry = sqldf("
select   x, color, treatment, variable, avg(value) mean, stdev(value) sd, sum(case when value is not null then 1 else 0 end) n
from     markers m
group by 1, 2, 3, 4
order by 1, 2, 3, 4
;")

msmry$l95 = msmry$mean - 1.96 * msmry$sd / sqrt(msmry$n)
msmry$u95 = msmry$mean + 1.96 * msmry$sd / sqrt(msmry$n)

msmry$display = plot_params$treatment[match(msmry$treatment, plot_params$treatment)]

varlabs = sqldf("
select   variable, avg(x) xmid
from     msmry
group by 1
order by 1
;")
varlabs$display[varlabs$variable=='cord_cd68_pct'] = 'cord\nCD68'
varlabs$display[varlabs$variable=='cord_iba1_pct'] = 'cord\nIba1'
varlabs$display[varlabs$variable=='cortex_iba1_pct'] = 'cortex\nIba1'

xlims=c(0.5,11.6)
ylims=c(0,1.50)

par(mar=c(5,5,3,2))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
abline(h=0,lwd=2)
par(xpd=T)
text(x=msmry$x, y=rep(0,nrow(msmry)), labels=msmry$display, adj=c(1,1), srt=45, font=1, col=msmry$color)
par(xpd=F)
axis(side=2,at=(0:6)/4,labels=percent((0:6)/4),las=2)
mtext(side=2, line=3, text='relative mRNA')
mtext(side=3, line=0, at=varlabs$xmid, text=varlabs$display, font=1)
abline(v=c(4,8,11.5))
abline(h=1,lwd=1,lty=3)
points(msmry$x, msmry$mean, type='h', lwd=20, col=msmry$color, lend=1)
arrows(x0=msmry$x, y0=msmry$l95, y1=msmry$u95, code=3, angle=90, lwd=1, length=0.05)

mtext('D', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off()





#### Figure S5

imgsave(paste('figures/script_generated/figure-s5.',imgmode,sep=''),width=8*resx,height=4*resx,res=resx)
layout_matrix = matrix(c(1,2,2,3,3),nrow=1,byrow=T)
layout(layout_matrix)

par(mar=c(4,0,4,0))
plot(NA, NA, xlim=c(0,1), ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

posthisto = read.table('data/nih-post-histo-s5.tsv',sep='\t',header=T)
posthisto$y = -posthisto$x # try this plot horizontal, negative sign so left = top
posthisto$grp = paste(posthisto$plot, posthisto$y)
# for (grp in unique(posthisto$grp)) {
#   rows = posthisto$grp == grp
#   posthisto$xbee[rows] = posthisto$x[rows] + beeswarm(value ~ grp, data=posthisto[rows,], cex=5, method='hex', spacing=1, corralWidth=0.5, do.plot=F)$x - 1
# }
# posthisto$ybee = -posthisto$xbee 
# beeswarm is too frustrating, try jitter
set.seed(1)
posthisto$ybee = jitter(posthisto$y, 0.35)
posthisto$color = plot_params$color[match(posthisto$tx, plot_params$treatment)]

ylabs1 = sqldf("
               select   tx display, y
               from     posthisto
               group by 1, 2
               order by 1, 2
               ;")
ylabs2 = sqldf("
               select   stage, min(y)-0.25 miny, max(y)+0.25 maxy, avg(y) meany
               from     posthisto
               group by 1
               order by 1
               ")
ylims=expand_range(posthisto$y,by=0.5)

# left half
pdat = subset(posthisto, thing=='ASO')
xlims = c(0,15)
par(mar=c(4,0,4,2))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=2, line=0.5, at=ylabs1$y, text=ylabs1$display, las=2)
par(xpd=T)
axis(side=2, at=c(ylabs2$miny[1], ylabs2$maxy[1]), labels=NA, lwd=1, tck=0.07, line=10)
axis(side=2, at=c(ylabs2$miny[2], ylabs2$maxy[2]), labels=NA, lwd=1, tck=0.07, line=10)
mtext(side=2, at=ylabs2$meany, text=ylabs2$stage, line=10.5, font=1)
par(xpd=F)
axis(side=1, at=(0:4)*5, cex.axis=1.2)
mtext(side=1, line=2.5, text='pixel count')
points(pdat$value, pdat$ybee, col=pdat$color, pch=19, cex=1.5)
mtext(side=3, line=0, text='ASO', cex=1.2, font=1)
mtext('A', side=3, cex=2, adj = 0.0, line = 0.5)


# right half
pdat = subset(posthisto, thing=='PrP')
xlims = c(0,80)
par(mar=c(4,1,4,1))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
axis(side=1, at=(0:4)*20, cex.axis=1.2)
mtext(side=1, line=2.5, text='pixel count')
points(pdat$value, pdat$ybee, col=pdat$color, pch=19, cex=1.5)
mtext(side=3, line=0, text='PrP', cex=1.2, font=1)
mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off()







imgsave(paste('figures/script_generated/figure-s3.',imgmode,sep=''),width=8*resx,height=4*resx,res=resx)
layout_matrix = matrix(c(1,2,2,3,3),nrow=1,byrow=T)
layout(layout_matrix)

par(mar=c(4,0,4,0))
plot(NA, NA, xlim=c(0,1), ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

prophhisto = read.table('data/nih-proph-histo-s3.tsv',sep='\t',header=T)
prophhisto$y = -prophhisto$x # try this plot horizontal, negative sign so left = top
prophhisto$grp = paste(prophhisto$plot, prophhisto$y)
# for (grp in unique(prophhisto$grp)) {
#   rows = prophhisto$grp == grp
#   prophhisto$xbee[rows] = prophhisto$x[rows] + beeswarm(value ~ grp, data=prophhisto[rows,], cex=5, method='hex', spacing=1, corralWidth=0.5, do.plot=F)$x - 1
# }
# prophhisto$ybee = -prophhisto$xbee 
# beeswarm is too frustrating, try jitter
set.seed(1)
prophhisto$ybee = jitter(prophhisto$y, 0.35)
prophhisto$color = plot_params$color[match(prophhisto$tx, plot_params$treatment)]

ylabs1 = sqldf("
               select   tx display, y
               from     prophhisto
               group by 1, 2
               order by 1, 2
               ;")
ylabs2 = sqldf("
               select   stage, min(y)-0.25 miny, max(y)+0.25 maxy, avg(y) meany
               from     prophhisto
               group by 1
               order by 1
               ")
ylims=expand_range(prophhisto$y,by=0.5)

# left half
pdat = subset(prophhisto, thing=='ASO')
xlims = c(0,25)
par(mar=c(4,0,4,2))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=2, line=0.5, at=ylabs1$y, text=ylabs1$display, las=2)
par(xpd=T)
axis(side=2, at=c(ylabs2$miny[1], ylabs2$maxy[1]), labels=NA, lwd=1, tck=0.07, line=10)
axis(side=2, at=c(ylabs2$miny[2], ylabs2$maxy[2]), labels=NA, lwd=1, tck=0.07, line=10)
mtext(side=2, at=ylabs2$meany, text=ylabs2$stage, line=10.5, font=1)
par(xpd=F)
axis(side=1, at=(0:5)*5, cex.axis=1.2)
mtext(side=1, line=2.5, text='pixel count')
points(pdat$value, pdat$ybee, col=pdat$color, pch=19, cex=1.5)
mtext(side=3, line=0, text='ASO', cex=1.2, font=1)
mtext('A', side=3, cex=2, adj = 0.0, line = 0.5)

# right half
pdat = subset(prophhisto, thing=='PrP')
xlims = c(0,80)
par(mar=c(4,1,4,1))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
axis(side=1, at=(0:4)*20, cex.axis=1.2)
mtext(side=1, line=2.5, text='pixel count')
points(pdat$value, pdat$ybee, col=pdat$color, pch=19, cex=1.5)
mtext(side=3, line=0, text='PrP', cex=1.2, font=1)
mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off()



moparms = data.frame(mouse=1:3,pch=c(16,17,15)) # different pch for different mice

imgsave(paste('figures/script_generated/figure-s4.',imgmode,sep=''),width=8*resx,height=4*resx,res=resx)
layout_matrix = matrix(c(1,2,2,3,3),nrow=1,byrow=T) 
layout(layout_matrix)

prophwb = read.table('data/nih-proph-wb-s4.tsv',sep='\t',header=T)
prophwb$pch = moparms$pch[match(prophwb$mouse,moparms$mouse)]

par(mar=c(4,0,4,0))
plot(NA, NA, xlim=c(0,1), ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

prophwb$y = -prophwb$x # try this plot horizontal, negative sign so left = top
prophwb$grp = paste(prophwb$plot, prophwb$y)
# for (grp in unique(prophwb$grp)) {
#   rows = prophwb$grp == grp
#   prophwb$xbee[rows] = prophwb$x[rows] + beeswarm(value ~ grp, data=prophwb[rows,], cex=1, method='hex', corral='gutter', corralWidth=0.5, do.plot=F)$x - 1
# }
# prophwb$ybee = -prophwb$xbee
# beeswarm is too frustrating, try jitter
set.seed(1)
prophwb$ybee = jitter(prophwb$y, 0.7)
prophwb$color = plot_params$color[match(prophwb$tx, plot_params$treatment)]

ylabs1 = sqldf("
               select   tx display, y
               from     prophwb
               group by 1, 2
               order by 1, 2
               ;")
ylabs2 = sqldf("
               select   stage, min(y)-0.25 miny, max(y)+0.25 maxy, avg(y) meany
               from     prophwb
               where    stage != ''
               group by 1
               order by 1
               ")
ylims=expand_range(prophwb$y,by=0.5)

# left half
pdat = subset(prophwb, thing=='-PK')
xlims = c(0,3)
par(mar=c(4,0,4,2.5))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=2, line=0.5, at=ylabs1$y, text=ylabs1$display, las=2)
par(xpd=T)
axis(side=2, at=c(ylabs2$miny[1], ylabs2$maxy[1]), labels=NA, lwd=1, tck=0.07, line=10)
axis(side=2, at=c(ylabs2$miny[2], ylabs2$maxy[2]), labels=NA, lwd=1, tck=0.07, line=10)
mtext(side=2, at=ylabs2$meany, text=ylabs2$stage, line=10.5, font=1)
par(xpd=F)
axis(side=1, at=(0:15)/5, labels=NA, tck=-0.01)
axis(side=1, at=(0:3), labels=percent(0:3), cex.axis=1.2, tck=-0.03)
abline(v=1,lwd=0.5,lty=3)
mtext(side=1, line=2.5, text='intensity')
points(pdat$value, pdat$ybee, col=pdat$color, pch=pdat$pch, lwd=3, cex=1.5)
mtext(side=3, line=0, text='-PK', cex=1.2, font=1)
mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

# right half
pdat = subset(prophwb, thing=='+PK')
xlims = c(0,3)
par(mar=c(4,1,4,1.5))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
axis(side=1, at=(0:15)/5, labels=NA, tck=-0.01)
axis(side=1, at=(0:3), labels=percent(0:3), cex.axis=1.2, tck=-0.03)
abline(v=1,lwd=0.5,lty=3)
mtext(side=1, line=2.5, text='intensity')
points(pdat$value, pdat$ybee, col=pdat$color, pch=pdat$pch, lwd=3, cex=1.5)
mtext(side=3, line=0, text='+PK', cex=1.2, font=1)
mtext('C', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off()







imgsave(paste('figures/script_generated/figure-s6.',imgmode,sep=''),width=8*resx,height=4*resx,res=resx)
layout_matrix = matrix(c(1,2,2,3,3),nrow=1,byrow=T)
layout(layout_matrix)

par(mar=c(4,0,4,0))
plot(NA, NA, xlim=c(0,1), ylim=ylims, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')

postwb = read.table('data/nih-post-wb-s6.tsv',sep='\t',header=T)
postwb$pch = moparms$pch[match(postwb$mouse, moparms$mouse)]
postwb$y = -postwb$x # try this plot horizontal, negative sign so left = top
postwb$grp = paste(postwb$plot, postwb$y)
set.seed(1)
postwb$ybee = jitter(postwb$y, 0.7)
postwb$color = plot_params$color[match(postwb$tx, plot_params$treatment)]

ylabs1 = sqldf("
               select   tx display, y
               from     postwb
               group by 1, 2
               order by 1, 2
               ;")
ylabs2 = sqldf("
               select   stage, min(y)-0.25 miny, max(y)+0.25 maxy, avg(y) meany
               from     postwb
               where    stage != ''
               group by 1
               order by 1
               ")
ylims=expand_range(postwb$y,by=0.5)

# left half
pdat = subset(postwb, thing=='-PK')
xlims = c(0,2)
par(mar=c(4,0,4,2.5))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
mtext(side=2, line=0.5, at=ylabs1$y, text=ylabs1$display, las=2)
par(xpd=T)
axis(side=2, at=c(ylabs2$miny[1], ylabs2$maxy[1]), labels=NA, lwd=1, tck=0.07, line=10)
axis(side=2, at=c(ylabs2$miny[2], ylabs2$maxy[2]), labels=NA, lwd=1, tck=0.07, line=10)
mtext(side=2, at=ylabs2$meany, text=ylabs2$stage, line=10.5, font=1)
par(xpd=F)
axis(side=1, at=(0:15)/5, labels=NA, tck=-0.01)
axis(side=1, at=(0:3), labels=percent(0:3), cex.axis=1.2, tck=-0.03)
abline(v=1,lwd=0.5,lty=3)
mtext(side=1, line=2.5, text='intensity')
points(pdat$value, pdat$ybee, col=pdat$color, pch=pdat$pch, lwd=3, cex=1.5)
mtext(side=3, line=0, text='-PK', cex=1.2, font=1)
mtext('B', side=3, cex=2, adj = 0.0, line = 0.5)

# right half
pdat = subset(postwb, thing=='+PK')
xlims = c(0,2)
par(mar=c(4,1,4,1.5))
plot(NA,NA,xlim=xlims, ylim=ylims,xaxs='i',yaxs='i',ann=F,axes=F)
axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0, las=2)
axis(side=1, at=(0:15)/5, labels=NA, tck=-0.01)
axis(side=1, at=(0:3), labels=percent(0:3), cex.axis=1.2, tck=-0.03)
abline(v=1,lwd=0.5,lty=3)
mtext(side=1, line=2.5, text='intensity')
points(pdat$value, pdat$ybee, col=pdat$color, pch=pdat$pch, lwd=3, cex=1.5)
mtext(side=3, line=0, text='+PK', cex=1.2, font=1)
mtext('C', side=3, cex=2, adj = 0.0, line = 0.5)

dev.off()




# write out a Table S2
write.table(mouse_survival[,c('expt','cohort','dpi','cause_of_death','acm')],'figures/script_generated/table_s2.tsv',sep='\t',quote=F,row.names=F,col.names=T)



