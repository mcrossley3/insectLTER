
setwd()
source('./code/AR_reml.R')

th = 9 #require time series to include at least 4 years of non-missing data

# Note: here Cedar Creek Sweep1 and Sweep2, and North Temperate Lakes Crayfish2 are recurated to remove treatments not deemed useful for estimating insect abundance trends.
# All other datasets are curated as in matters.arising.code_m4_standard_subsite.R

#############################################################################################
# Cedar Creek Ecosystem Arthropods Sweep 1
# https://www.cedarcreek.umn.edu/research/data/dataset?arce153

data1 = read.table('./raw_data/e153_Arthropod sweepnet sampling.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='undet under'),]
data1$Number = data1$Specimens
data1$PEFB = paste(data1$Plot,data1$Exclosure,data1$Fertilized,data1$Burned,sep='_') #each plot has a specific treatment regime
plots = unique(data1$Plot)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,'Exclosure'=NA,'Fertilized'=NA,'Burned'=NA,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$Plot==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'CedarCreek' #LTER
			out[ox,2] = plots[i] #"Locale" (watershed*plant species*fire frequency)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (in this case, number of dates in which sampling is reported)
			out[ox,6] = mean(data4$Number,na.rm=T)
			out[ox,7] = unique(data4$Exclosure)
			out[ox,8] = unique(data4$Fertilized)
			out[ox,9] = unique(data4$Burned)
			ox = ox + 1
		}
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999,'Exclosure'=NA,'Fertilized'=NA,'Burned'=NA)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > th) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
			cmin = min(dat2$Abundance[which(dat2$Abundance>0)],na.rm=T)
			if (is.na(cmin)){
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 
			} else {
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 * cmin #replace zeroes with 0.5*minimum abundance value in this time series
			}
			dat2$log.Abundance = log(dat2$Abundance) #log-transform abundances
			if (length(unique(dat2$log.Abundance[which(!is.na(dat2$log.Abundance))]))==1){ #abundances are all = 1
				abundance.trend = 0
			} else {
				ys = sort(dat2$Year)
				ys.stretch = seq(ys[1],ys[length(ys)],1) #expand time series to include missing years
				X = match(ys.stretch,ys)
				X[which(!is.na(X))] = dat2$log.Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				if (nrow(dat2)==2){
					abundance.trend = (Z[2] - Z[1]) / 2 #simple slope for 2-yr time series
				} else {
					arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
					abundance.trend = arr.Z$coef[2]
				}
#					png(paste0('./plots/CedarCreek/sweep1/',s,'_',locales[l],'.png'))
#					plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'CedarCreek'
			trends[tx,2] = 'sweep1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			trends[tx,6] = unique(dat2$Exclosure)
			trends[tx,7] = unique(dat2$Fertilized)
			trends[tx,8] = unique(dat2$Burned)
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}

# none of these trends passed filtering criteria
write.table(trends,paste0('./summary_tables/m10_rmtreat/CedarCreek_Sweep1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m10_rmtreat/CedarCreek_Sweep1_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Cedar Creek Ecosystem Arthropods Sweep 2
# https://www.cedarcreek.umn.edu/research/data/dataset?aage120

data1 = read.table('./raw_data/e120_Main Plots All Arthropod Insect Sweepnet Sampling 1996-2006.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet_undet'),]
data1$Number = as.numeric(data1$Count)
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='unk unk' & data1$Species!='none none' & data1$Species!='na? na?' & data1$Species!='na na'),]
data1$PS = paste(data1$Plot,data1$NumSp,sep='_') #each plot has a specific plant diversity treatment
plots = unique(data1$Plot)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,'Plants'=NA,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$Plot==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'CedarCreek' #LTER
			out[ox,2] = plots[i] #"Locale" (plot)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (in this case, number of dates in which sampling is reported)
			out[ox,6] = mean(data4$Number,na.rm=T)
			out[ox,7] = unique(apply(data4[,15:32],1,function(x){paste0(x,collapse='')}))
			ox = ox + 1
		}
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999,'Plants'=NA)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > th) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
			cmin = min(dat2$Abundance[which(dat2$Abundance>0)],na.rm=T)
			if (is.na(cmin)){
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 
			} else {
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 * cmin #replace zeroes with 0.5*minimum abundance value in this time series
			}
			dat2$log.Abundance = log(dat2$Abundance) #log-transform abundances
			if (length(unique(dat2$log.Abundance[which(!is.na(dat2$log.Abundance))]))==1){ #abundances are all = 1
				abundance.trend = 0
			} else {
				ys = sort(dat2$Year)
				ys.stretch = seq(ys[1],ys[length(ys)],1) #expand time series to include missing years
				X = match(ys.stretch,ys)
				X[which(!is.na(X))] = dat2$log.Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				if (nrow(dat2)==2){
					abundance.trend = (Z[2] - Z[1]) / 2 #simple slope for 2-yr time series
				} else {
					arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
					abundance.trend = arr.Z$coef[2]
				}
#					png(paste0('./plots/CedarCreek/sweep2/',s,'_',locales[l],'.png'))
#					plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'CedarCreek'
			trends[tx,2] = 'sweep2'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			trends[tx,6] = unique(dat2$Plants)
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s],'Plants'=unique(dat2$Plants))),stringsAsFactors=F)
		}	
	}
}

plot(factor(trends$Plants),trends$Abundance.trend)

# Remove plots with any experimental seeding
out2 = out[which(out$Plants=='000000000000000000'),]
trends2 = trends[which(trends$Plants=='000000000000000000'),] #none of the time series for these plots met filtering criteria
write.table(trends2,paste0('./summary_tables/m10_rmtreat/CedarCreek_Sweep2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out2,paste0('./summary_tables/m10_rmtreat/CedarCreek_Sweep2_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Temperate Lakes - Crayfish 2
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.269.2

data1 = read.csv('./raw_data/sparkling_crayfish.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$YEAR4
species = c('Orconectes rusticus','Orconectes virilis')
spcol = match(c('Rusty CPUE','Virilis CPUE'),colnames(data1))
years = sort(unique(data1$Year))
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	for (j in 1:length(years)){
		data2 = data1[which(data1$Year==years[j]),]
		out[ox,1] = 'NorthTemperateLakes' #LTER
		out[ox,2] = 'SP' #"Locale"
		out[ox,3] = species[s]
		out[ox,4] = years[j] #Year
		out[ox,5] = sum(data2$TRAP_DAYS,na.rm=T) #"N.Obs" (number of trap days)
		out[ox,6] = mean(data2[,spcol[s]],na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > th) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
			cmin = min(dat2$Abundance[which(dat2$Abundance>0)],na.rm=T)
			if (is.na(cmin)){
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 
			} else {
				dat2$Abundance[which(dat2$Abundance==0)] = 0.5 * cmin #replace zeroes with 0.5*minimum abundance value in this time series
			}
			dat2$log.Abundance = log(dat2$Abundance) #log-transform abundances
			if (length(unique(dat2$log.Abundance[which(!is.na(dat2$log.Abundance))]))==1){ #abundances are all = 1
				abundance.trend = 0
			} else {
				ys = sort(dat2$Year)
				ys.stretch = seq(ys[1],ys[length(ys)],1) #expand time series to include missing years
				X = match(ys.stretch,ys)
				X[which(!is.na(X))] = dat2$log.Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				if (nrow(dat2)==2){
					abundance.trend = (Z[2] - Z[1]) / 2 #simple slope for 2-yr time series
				} else {
					arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
					abundance.trend = arr.Z$coef[2]
				}
#					png(paste0('./plots/NorthTemperateLakes/Crayfish2/',species[s],'_',locales[l],'.png'))
#					plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'NorthTemperateLakes'
			trends[tx,2] = 'Crayfish2'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}
# Remove O. rusticus data (because trapped individuals were removed while O. virilis individuals were released again)
out = out[which(out$Species.code=='Orconectes virilis'),]
trends = trends[which(trends$Species=='Orconectes virilis'),]
write.table(trends,paste0('./summary_tables/m10_rmtreat/NorthTemperateLakes_Crayfish2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m10_rmtreat/NorthTemperateLakes_Crayfish2_abundance.txt'),sep='\t',quote=F,row.names=F)


######################################################################################################
# Compile all trends into a single dataframe

files = list.files(paste0('./summary_tables/m10_rmtreat'),full.names=T)
files = files[grep('_trends.txt',files)]
trends = c()
for (i in 1:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)
	if (length(which(is.na(add$LTER)))==nrow(add)){ #skip datasets that had no time series retained
		print(noquote(files[i]))
	} else {
		trends = data.frame(rbind(trends,add),stringsAsFactors=F)
	}
}
write.table(trends,'./summary_tables/m10/time_trends_arthropods_m10_wDatasetNames.txt',sep='\t',row.names=F,quote=F)
# Relabel datasets for plotting
trends$Dataset[which(trends$Dataset=='sweep1' | trends$Dataset=='sweep2')] = 'Sweeps'
trends$Dataset[which(trends$Dataset=='AntsNantucket')] = 'Ants'
trends$Dataset[which(trends$Dataset=='lepidoptera1' | trends$Dataset=='lepidoptera2')] = 'Leptidoptera'
trends$Dataset[which(trends$Dataset=='BenthicMacroinvertebrates' | trends$Dataset=='PelagicMacroinvertebrates')] = 'Macroinvertebrates'
trends$Dataset[which(trends$Dataset=='Crayfish1' | trends$Dataset=='Crayfish2')] = 'Crayfish'
trends$Dataset[which(trends$Dataset=='pitfall1' | trends$Dataset=='pitfall2')] = 'Pitfalls'
trends$Dataset[which(trends$Dataset=='Invertebrates' & trends$Species=='Prokelisia')] = 'PlantHoppers'
trends$Dataset[which(trends$Dataset=='Invertebrates' & trends$Species=='Acrididae')] = 'Grasshoppers'
trends$Dataset[which(trends$Dataset=='BurrowingCrabs' | trends$Dataset=='FiddlerCrabs')] = 'Crabs'
write.table(trends,'./summary_tables/m10_rmtreat/time_trends_arthropods_m10.txt',sep='\t',row.names=F,quote=F)

# Compile abundances into a single file
files = list.files(paste0('./summary_tables/m10_rmtreat'),full.names=T)
files = files[grep('_abundance.txt',files)]
abundances = read.table(files[1],sep='\t',as.is=T,check.names=F,header=T)
abundances$Dataset = strsplit(strsplit(files[1],'/')[[1]][5],'_')[[1]][2]
for (i in 2:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)[,1:6]
	add$Dataset = strsplit(strsplit(files[i],'/')[[1]][5],'_')[[1]][2]
	abundances = data.frame(rbind(abundances,add),stringsAsFactors=F)
}
write.table(abundances,'./summary_tables/m10_rmtreat/PerSpecies_Abundance_LTER_m10.txt',sep='\t',row.names=F,quote=F)


######################################################################################################
# Mean & 95% CI

dat = read.table(paste0('./summary_tables/m10_rmtreat/time_trends_arthropods_m10.txt'),sep='\t',as.is=T,check.names=F,header=T)
dat$LTER.Dataset = paste(dat$LTER,dat$Dataset,sep=' ')
levels1 = c('Arctic StreamInsects','Coweeta StreamInvertebrates','NorthTemperateLakes Crayfish','NorthTemperateLakes Macroinvertebrates','BonanzaCreek BarkBeetles','BonanzaCreek Leafminers','CedarCreek Grasshoppers',
	'CedarCreek Sweeps','GeorgiaCoastal Crabs','GeorgiaCoastal Grasshoppers','HarvardForest Ticks','HubbardBrook Leptidoptera','KonzaPrairie Grasshoppers','Sevilleta Grasshoppers','Midwest Aphids','Phoenix Pitfalls')
dat$LTER.Dataset = factor(dat$LTER.Dataset,levels=levels1)
dat$Habitat = 'terrestrial'
dat$Habitat[which(dat$LTER=='Arctic' | dat$LTER=='Coweeta' | dat$LTER=='NorthTemperateLakes')] = 'aquatic'
dat$Habitat[which(dat$LTER=='Midwest' | dat$LTER=='Baltimore' | dat$LTER=='Phoenix')] = 'impacted'
dat$Habitat = factor(dat$Habitat,levels=c('aquatic','terrestrial','impacted'))

dat$Sub = paste(dat$LTER,dat$Dataset,dat$Locale,sep='_')
dat$Spe = paste(dat$LTER,dat$Dataset,dat$Species,sep='_')
Spe = unique(dat$Spe)
Sub = unique(dat$Sub)
m2 = apply(array(Spe),1,function(x){mean(dat$Abundance.trend[which(dat$Spe==x)],na.rm=T)}) #mean of species means
dat2 = dat[match(Spe,dat$Spe),]
dat2$Species.mean.trend = m2
# Split KonzaPrairie Grasshoppers by grazed and ungrazed
dat2 = dat2[which(dat2$LTER.Dataset!='KonzaPrairie Grasshoppers'),]
kg = dat[which(dat$LTER.Dataset=='KonzaPrairie Grasshoppers'),]
kg.sp = unique(kg$Species)
add1 = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999,'LTER.Dataset'=NA,'Habitat'=NA,'Sub'=NA,'Spe'=NA,'Species.mean.trend'=-999)
ax = 1
for (i in 1:length(kg.sp)){
	kg2 = kg[which(kg$Species==kg.sp[i]),]
	gpos = which(kg2$Locale=='N01A' | kg2$Locale=='N01B' | kg2$Locale=='N04A' | kg2$Locale=='N04D' | kg2$Locale=='N20A' | kg2$Locale=='N20B')
	upos = which(kg2$Locale=='001D' | kg2$Locale=='0SpB' | kg2$Locale=='002C' | kg2$Locale=='002D' | kg2$Locale=='004F' | kg2$Locale=='020B')
	if (length(gpos)>0){
		add1[ax,1] = 'KonzaPrairie'
		add1[ax,2] = 'Grasshoppers'
		add1[ax,3] = 'Grazed'
		add1[ax,4] = kg.sp[i]
		add1[ax,5] = NA
		add1[ax,6] = 'KonzaPrairie Grasshoppers'
		add1[ax,7] = 'terrestrial'
		add1[ax,8] = NA
		add1[ax,9] = paste('KonzaPrairie_Grasshoppers',kg.sp[i],sep='_')
		add1[ax,10] = mean(kg2$Abundance.trend[gpos])
		ax = ax + 1
	} else {
	}
	if (length(upos)>0){
		add1[ax,1] = 'KonzaPrairie'
		add1[ax,2] = 'Grasshoppers'
		add1[ax,3] = 'Ungrazed'
		add1[ax,4] = kg.sp[i]
		add1[ax,5] = NA
		add1[ax,6] = 'KonzaPrairie Grasshoppers'
		add1[ax,7] = 'terrestrial'
		add1[ax,8] = NA
		add1[ax,9] = paste('KonzaPrairie_Grasshoppers',kg.sp[i],sep='_')
		add1[ax,10] = mean(kg2$Abundance.trend[upos])
		ax = ax + 1
	} else {
	}

}
mean(add1$Species.mean.trend[which(add1$Locale=='Grazed')])
mean(add1$Species.mean.trend[which(add1$Locale=='Ungrazed')])
dat2 = data.frame(rbind(dat2,add1),stringsAsFactors=F)
L = unique(dat2$LTER)
means2 = apply(array(L),1,function(x){mean(dat2$Species.mean.trend[which(dat2$LTER==x)],na.rm=T)})

t2 = t.test(means2)

# Barplots
mn2 = t2[['estimate']]
uppers = t2[['conf.int']][2]
lowers = t2[['conf.int']][1]
png(paste0('./plots/m10_rmtreat/t-test_abundance_lters_barplot_Species.png'),width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mn2,ylim=c(-1.5,1.5),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
dev.off()


######################################################################################################
# Recreate violin plot (Figure 2) using updated trends

library(ggplot2)

# Species
dp2 = ggplot(dat2, aes(x=LTER.Dataset, y=Species.mean.trend, fill=Habitat)) + 
	geom_violin(trim=T,size=1) +
	scale_fill_manual(values=c('lightblue','white','darkorange')) +
	ylim(-6,6) +
	stat_summary(fun=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1) +
	labs(y="Species mean abundance trend (sd/yr)") +
	theme_classic() +
	theme(axis.ticks.y = element_line(size=1),
		axis.line.y = element_line(size=1),
		axis.text.y = element_text(size=16,colour='black'),
		axis.title.y = element_text(colour='black',size=20),
		axis.ticks.length.y = unit(.15, "cm")) +
	theme(axis.line.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.x = element_blank()) +
	theme(legend.position='none')
	ggsave(filename='./plots/m10_rmtreat/ALL_time_trends_boxplot_relaxed_noannotations_violin_Species.png',plot=dp2,dpi=600,unit='cm',height=8,width=24)
tcount = apply(array(levels(dat2$LTER.Dataset)),1,function(x){length(which(as.character(dat2$LTER.Dataset)==x))})
cbind(levels(dat2$LTER.Dataset),tcount)

data1 = read.table('./summary_tables/m10_rmtreat/PerSpecies_Abundance_LTER_m10.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$LTER.Dataset = paste(data1$LTER,data1$Dataset,sep='_')
ycount = apply(array(unique(data1$LTER.Dataset)),1,function(x){sort(unique(data1$Year[which(data1$LTER.Dataset==x)]))})
names(ycount) = unique(data1$LTER.Dataset)

