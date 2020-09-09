
setwd('C:/Users/mcros/Desktop/Postdoc UGA/LTER_insect_diversity')
source('./code/AR_reml.R')

#############################################################################################
# Konza Prairie grasshoppers
# https://doi.org/10.6073/pasta/7b2259dcb0e499447e0e11dfb562dc2f

data1 = read.csv('./raw_data/CGR022.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
# Fix species names
#write.csv(u.species,'KonzaPrairie_grasshoppers_duplicatenames.csv',quote=F,row.names=F)
# Merge duplicates that have different species name spelling
konza.duplicates = read.csv('./taxa_keys/KonzaPrairie_grasshoppers_duplicatenames.csv',as.is=T,check.names=F,header=F)
data1 = data1[which(data1$Species!='Unknown'),]
for (i in 1:nrow(data1)){
	name.swap = konza.duplicates[which(konza.duplicates[,1]==data1$Species[i]),2]
	if (length(name.swap)==0){
		name.swap = konza.duplicates[match(data1$Species[i],konza.duplicates[,2]),2]
	} else {
	}
	data1$Species[i] = name.swap

}
wsd = unique(data1$WATERSHED)
wsd.count = apply(array(wsd),1,function(x){length(which(data1$WATERSHED==x))})

# Match inconsistently spelled WATERSHED values
data1$WATERSHED[which(data1$WATERSHED=='001d')] = '001D'
data1$WATERSHED[which(data1$WATERSHED=='002c')] = '002C'
data1$WATERSHED[which(data1$WATERSHED=='002d')] = '002D'
data1$WATERSHED[which(data1$WATERSHED=='004b')] = '004B'
data1$WATERSHED[which(data1$WATERSHED=='004f')] = '004F'
data1$WATERSHED[which(data1$WATERSHED=='020b')] = '020B'
data1$WATERSHED[which(data1$WATERSHED=='0spb' | data1$WATERSHED=='0SPB')] = '0SpB'
data1$WATERSHED[which(data1$WATERSHED=='0sub' | data1$WATERSHED=='0SUB')] = '0SuB'
data1$WATERSHED[which(data1$WATERSHED=='n01a')] = 'N01A'
data1$WATERSHED[which(data1$WATERSHED=='n01b')] = 'N01B'
data1$WATERSHED[which(data1$WATERSHED=='n04a')] = 'N04A'
data1$WATERSHED[which(data1$WATERSHED=='n04d')] = 'N04D'
data1$WATERSHED[which(data1$WATERSHED=='n20a')] = 'N20A'
data1$WATERSHED[which(data1$WATERSHED=='n20b')] = 'N20B'
wsd = unique(data1$WATERSHED)
wsd.years = apply(array(wsd),1,function(x){unique(data1$Year[which(data1$WATERSHED==x)])}) #check years for undescribed WATERSHED values in metadata. 
ungrazed = c('002D','001D','0SuB','004F','020B','002C','0SpB','004B','n00b','010d','004g','000b','004d')
grazed = c('N20B','N01B','N04D','N04A','N01A','N20A')
wrdmy = paste(data1$WATERSHED,data1$REPSITE,data1$RECYEAR,data1$RECMONTH,data1$RECDAY,sep='_')
wrdmy.count = apply(array(wrdmy),1,function(x){length(which(wrdmy==x))})

watersheds = unique(data1$WATERSHED)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(watersheds)){
	data2 = data1[which(data1$WATERSHED==watersheds[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		species3 = sort(unique(data3$Species))
		days = paste(data3$RECMONTH,data3$RECDAY,sep='_')
		for (k in 1:length(species3)){
			out[ox,1] = 'KonzaPrairie' #LTER
			out[ox,2] = watersheds[i] #"Locale"
			out[ox,3] = species3[k] #Species.code (in this case, actual species name)
			out[ox,4] = years2[j] #Year
			out[ox,5] = length(unique(days)) #"N.Obs" (in this case, the number of days over which samples were taken)
			out[ox,6] = mean(data3$Number[which(data3$Species==species3[k])],na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Konza/grasshoppers/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'KonzaPrairie'
			trends[tx,2] = 'grasshoppers'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}


# Repeat, this time excluding any counts occuring in June
watersheds = unique(data1$WATERSHED)
out2 = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(watersheds)){
	data2 = data1[which(data1$WATERSHED==watersheds[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j] & data2$RECMONTH!=6),]
		species3 = sort(unique(data3$Species))
		days = paste(data3$RECMONTH,data3$RECDAY,sep='_')
		for (k in 1:length(species3)){
			out2[ox,1] = 'KonzaPrairie' #LTER
			out2[ox,2] = watersheds[i] #"Locale"
			out2[ox,3] = species3[k] #Species.code (in this case, actual species name)
			out2[ox,4] = years2[j] #Year
			out2[ox,5] = length(unique(days)) #"N.Obs" (in this case, the number of days over which samples were taken)
			out2[ox,6] = mean(data3$Number[which(data3$Species==species3[k])],na.rm=T)
			ox = ox + 1
		}
	}
}; str(out2)

out$LSY = paste(out$Locale,out$Species.code,out$Year,sep='_')
out2$LSY = paste(out2$Locale,out2$Species.code,out2$Year,sep='_')
merge1 = merge(out,out2,by='LSY',all.x=F)
plot(merge1$Abundance.x,merge1$Abundance.y,xlab='All months',ylab='June excluded')

# Estimate abundance trends
species = sort(unique(out2$Species.code))
trends2 = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost2 = c()
tx = 1
for (s in 1:length(species)){
	dat = out2[which(out2$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
			}
			trends2[tx,1] = 'KonzaPrairie'
			trends2[tx,2] = 'grasshoppers'
			trends2[tx,3] = locales[l]
			trends2[tx,4] = species[s]
			trends2[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost2 = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}

trends$LSY = paste(trends$Locale,trends$Species,sep='_')
trends2$LSY = paste(trends2$Locale,trends2$Species,sep='_')
merge2 = merge(trends,trends2,by='LSY',all.x=F)
plot(merge2$Abundance.trend.x,merge2$Abundance.trend.y,xlab='All months',ylab='June excluded')


# Repeat, this time excluding any counts occuring prior to 1996
data1. = data1[which(data1$Year>=1996),] #32% of data is removed at this step
watersheds = unique(data1.$WATERSHED)
out3 = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(watersheds)){
	data2 = data1.[which(data1.$WATERSHED==watersheds[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j] & data2$RECMONTH!=6),]
		species3 = sort(unique(data3$Species))
		days = paste(data3$RECMONTH,data3$RECDAY,sep='_')
		for (k in 1:length(species3)){
			out3[ox,1] = 'KonzaPrairie' #LTER
			out3[ox,2] = watersheds[i] #"Locale"
			out3[ox,3] = species3[k] #Species.code (in this case, actual species name)
			out3[ox,4] = years2[j] #Year
			out3[ox,5] = length(unique(days)) #"N.Obs" (in this case, the number of days over which samples were taken)
			out3[ox,6] = mean(data3$Number[which(data3$Species==species3[k])],na.rm=T)
			ox = ox + 1
		}
	}
}; str(out3)

# Estimate abundance trends
species = sort(unique(out3$Species.code))
trends3 = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost3 = c()
tx = 1
for (s in 1:length(species)){
	dat = out3[which(out3$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
			}
			trends3[tx,1] = 'KonzaPrairie'
			trends3[tx,2] = 'grasshoppers'
			trends3[tx,3] = locales[l]
			trends3[tx,4] = species[s]
			trends3[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost3 = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}
trends3$LSY = paste(trends3$Locale,trends3$Species,sep='_')
merge3 = merge(trends,trends3,by='LSY',all.x=F)
plot(merge3$Abundance.trend.x,merge3$Abundance.trend.y,xlab='All years',ylab='pre-1996 excluded')
abline(h=0); abline(v=0)


# Repeat, this time excluding any counts occuring in June and any counts prior to 1996
watersheds = unique(data1.$WATERSHED)
out4 = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(watersheds)){
	data2 = data1.[which(data1.$WATERSHED==watersheds[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j] & data2$RECMONTH!=6),]
		species3 = sort(unique(data3$Species))
		days = paste(data3$RECMONTH,data3$RECDAY,sep='_')
		for (k in 1:length(species3)){
			out4[ox,1] = 'KonzaPrairie' #LTER
			out4[ox,2] = watersheds[i] #"Locale"
			out4[ox,3] = species3[k] #Species.code (in this case, actual species name)
			out4[ox,4] = years2[j] #Year
			out4[ox,5] = length(unique(days)) #"N.Obs" (in this case, the number of days over which samples were taken)
			out4[ox,6] = mean(data3$Number[which(data3$Species==species3[k])],na.rm=T)
			ox = ox + 1
		}
	}
}; str(out4)

# Estimate abundance trends
species = sort(unique(out4$Species.code))
trends4 = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost2 = c()
tx = 1
for (s in 1:length(species)){
	dat = out4[which(out4$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
			}
			trends4[tx,1] = 'KonzaPrairie'
			trends4[tx,2] = 'grasshoppers'
			trends4[tx,3] = locales[l]
			trends4[tx,4] = species[s]
			trends4[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost2 = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}
trends4$LSY = paste(trends4$Locale,trends4$Species,sep='_')
merge4 = merge(trends3,trends4,by='LSY',all.x=F)
plot(merge4$Abundance.trend.x,merge4$Abundance.trend.y,xlab='All data',ylab='pre-1996 and June excluded')
abline(h=0); abline(v=0)

par(mfrow=c(2,2))
#a
hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*watershed',cex.lab=1.5,cex.axis=1.2,col='grey40',breaks=seq(-10,10,1),main='All years & months')
abline(v=0,col='pink',lwd=2)
#b
plot(merge3$Abundance.trend.x,merge3$Abundance.trend.y,xlab='All years & months',ylab='pre-1996 excluded',cex.lab=1.5,cex.axis=1.2,pch=16,col='grey40',cex=1.2)
abline(h=0); abline(v=0)
#c
hist(trends4$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*watershed',cex.lab=1.5,cex.axis=1.2,col='grey40',breaks=seq(-10,10,1),main='Excluding pre-1996 & June')
abline(v=0,col='pink',lwd=2)
dev.off()

trends4 = trends4[,1:5]
write.table(trends4,'./summary_tables/reWelti/Konza_Grasshoppers_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out4,'./summary_tables/reWelti/Konza_Grasshoppers_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Konza Prairie gall insects
# http://dx.doi.org/10.6073/pasta/b2ac9e918a66dbbb18c7a6b39dc1efab

data1 = read.csv('./raw_data/CGP011.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$RecYear
data1$Number = data1$GalledStems / data1$SampledStems #galls per stem
# Fix plant species name typoes
data1$Species[which(data1$Species=='ceanothu')] = 'ceanothus'
data1$Watershed.Plant.FireFreq = paste(data1$Watershed,data1$Species,data1$fireFrequency,sep='_')
wsd = unique(data1$Watershed.Plant.FireFreq)
wsd.count = apply(array(wsd),1,function(x){length(which(data1$Watershed.Plant.FireFreq==x))})

out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(wsd)){
	data2 = data1[which(data1$Watershed.Plant.FireFreq==wsd[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'KonzaPrairie' #LTER
		out[ox,2] = wsd[i] #"Locale" (watershed*plant species*fire frequency)
		out[ox,3] = 'gall insects'
		out[ox,4] = years2[j] #Year
		out[ox,5] = length(unique(data3$CensusReplicate)) #"N.Obs" (in this case, the number of sample replicates)
		out[ox,6] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}; str(out)

# Estimate abundance trends
locales = sort(unique(out$Locale))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (l in 1:length(locales)){
	dat2 = out[which(out$Locale==locales[l]),]
	if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
			X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(ys.stretch)
			t.scale = (t.scale-min(t.scale))/max(t.scale)
			arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
			abundance.trend = arr.Z$coef[2]
			png(paste0('./plots/Welti reanalysis/Konza/galls/',locales[l],'.png'))
			plot(t.scale,Z,main=locales[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
			abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
			dev.off()
		}
		trends[tx,1] = 'KonzaPrairie'
		trends[tx,2] = 'GallInsects'
		trends[tx,3] = locales[l]
		trends[tx,4] = NA
		trends[tx,5] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
		lost2 = data.frame(rbind(lost,data.frame('Locale'=locales[l])),stringsAsFactors=F)
	}
}

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. watershed*plant*fire frequency',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

trends$FireFreq = apply(array(trends$Locale),1,function(x){strsplit(x,'_')[[1]][3]})
trends$PlantSpecies = apply(array(trends$Locale),1,function(x){strsplit(x,'_')[[1]][2]})
plot(as.factor(trends$PlantSpecies),trends$Abundance.trend)
plot(as.factor(trends$FireFreq),trends$Abundance.trend)
par(oma=c(1,1,2,1),mar=c(5,5,2,2))
plot(factor(paste(trends$PlantSpecies,trends$FireFreq,sep='_'),
	levels=c('ceanothus_1','ceanothus_4','ceanothus_10','ceanothus_20',
	'solidago_1','solidago_4','solidago_10','solidago_20',
	'vernonia_1','vernonia_4','vernonia_10','vernonia_20')),
	trends$Abundance.trend,xaxt='n',xlab='Fire frequency',ylab='Abundance trend (sd/yr)',cex.lab=1.5); abline(v=c(4.5,8.5)); abline(h=0,col='pink',lwd=2)
mtext(text=rep(c(1,4,10,20),4),side=1,at=1:12,cex=1.2)
mtext(text=c('ceanothus','solidago','vernonia'),side=3,line=0,outer=F,at=c(2.5,6.5,10.5),cex=1.5)
dev.off()

par(mfrow=c(1,3))
hist(trends$Abundance.trend[which(trends$PlantSpecies=='ceanothus')],xlab='Abundance trend (sd/yr)',ylab='No. watershed*fire frequency',cex.lab=1.5,cex.axis=1.2,col='grey40',main='ceanothus'); abline(v=0,col='pink',lwd=2)
hist(trends$Abundance.trend[which(trends$PlantSpecies=='solidago')],xlab='Abundance trend (sd/yr)',ylab='No. watershed*fire frequency',cex.lab=1.5,cex.axis=1.2,col='grey40',main='solidago'); abline(v=0,col='pink',lwd=2)
hist(trends$Abundance.trend[which(trends$PlantSpecies=='vernonia')],xlab='Abundance trend (sd/yr)',ylab='No. watershed*fire frequency',cex.lab=1.5,cex.axis=1.2,col='grey40',main='vernonia'); abline(v=0,col='pink',lwd=2)
dev.off()

trends = trends[,1:5]
write.table(trends,'./summary_tables/reWelti/Konza_GallInsects_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/Konza_GallInsects_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################5
# Cedar Creek Ecosystem - Old Field Grasshopper Sampling
# https://www.cedarcreek.umn.edu/research/data/methods?e014

data1 = read.table('./raw_data/e014_Core Old Field Grasshopper Sampling.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1$Number = data1$Specimens
fields = sort(unique(data1$Field.num))
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(fields)){
	data2 = data1[which(data1$Field.num==fields[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'CedarCreek' #LTER
			out[ox,2] = fields[i] #"Locale" (watershed*plant species*fire frequency)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (in this case, number of months in which sampling is reported)
			out[ox,6] = mean(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/CedarCreek/grasshoppers/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'CedarCreek'
			trends[tx,2] = 'grasshoppers'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*fields',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/CedarCreek_Grasshoppers_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/CedarCreek_Grasshoppers_abundance.txt',sep='\t',quote=F,row.names=F)


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
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
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
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/CedarCreek/sweep1/',s,'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'CedarCreek'
			trends[tx,2] = 'sweep1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/CedarCreek_Sweep1_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/CedarCreek_Sweep1_abundance.txt',sep='\t',quote=F,row.names=F)


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
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
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
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/CedarCreek/sweep2/',s,'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'CedarCreek'
			trends[tx,2] = 'sweep2'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/CedarCreek_Sweep2_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/CedarCreek_Sweep2_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Central Arizona-Phoenix - Arthropod sweeps
# https://doi.org/10.6073/pasta/0669ee6a71b24abb1ae3827f4ee77f6d

data1 = read.csv('./raw_data/652_arthropods_e9c22403e0f7243f241ed64862e22e05.csv',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$arthropod_scientific_name
data1$Species = gsub('\\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
data1$Number = data1$number_of_arthropods
apply(array(unique(data1$Site)),1,function(x){unique(data1$sample_date[which(data1$Site==x)])})
plots = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$Site==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'Phoenix' #LTER
			out[ox,2] = plots[i] #"Locale" (plot)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = ceiling(nrow(data4) / length(unique(data4$arthropod_scientific_name))) #"N.Obs" (in this case, number of sweep samples with counts reported) Should =3, but zeroes are not reported.
			out[ox,6] = sum(data4$Number,na.rm=T) / 3 #metadata reports 3 sweep samples per plot
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Phoenix/sweep/',s,'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'Phoenix'
			trends[tx,2] = 'sweep'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/Phoenix_Sweep_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/Phoenix_Sweep_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Central Arizona-Phoenix - Arthropod Pitfall 1
# https://data.sustainability.asu.edu/cap-portal/mapbrowse?packageid=knb-lter-cap.41.16

data1 = read.csv('./raw_data/41_core_arthropods_1c843549f049757465771bf4c36e80ab.csv',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
dim(data1) #109,143 rows
data1 = data1[which(!is.na(data1$Species)),]; dim(data1) #105,945 rows
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]; dim(data1) #105,401 rows
plots = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$Site==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			n.traps = length(unique(data2$trap_name[which(data2$Year==years3[j])]))
			out[ox,1] = 'Phoenix' #LTER
			out[ox,2] = plots[i] #"Locale" (plot)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = n.traps #"N.Obs" (in this case, number of pitfall traps for which counts are reported at a site)
			out[ox,6] = sum(data4$Number,na.rm=T) / n.traps #metadata reports 3 sweep samples per plot
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Phoenix/pitfall1/',s,'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'Phoenix'
			trends[tx,2] = 'pitfall1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

# Repeat, removing any observations with an oddly high or low number of traps reported.
dim(out) #22,594 rows
out2 = out[which(out$N.Obs<=12 & out$N.Obs>=10),]; dim(out2) #20,539 rows

# Estimate abundance trends
species = sort(unique(out2$Species.code))
trends2 = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost2 = c()
tx = 1
for (s in 1:length(species)){
	dat = out2[which(out2$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
#				png(paste0('./plots/Welti reanalysis/Phoenix/pitfall1/',s,'_',locales[l],'.png'))
#				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#				dev.off()
			}
			trends2[tx,1] = 'Phoenix'
			trends2[tx,2] = 'pitfall1'
			trends2[tx,3] = locales[l]
			trends2[tx,4] = species[s]
			trends2[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost2 = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends2)

hist(trends2$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='Only counts from 10-12 traps')
abline(v=0,col='pink',lwd=2)

trends$SL = paste(trends$Locale,trends$Species,sep='_')
trends2$SL = paste(trends2$Locale,trends2$Species,sep='_')
merge1 = merge(trends,trends2,by='SL',all.y=F)
plot(merge1$Abundance.trend.x,merge1$Abundance.trend.y,xlab='All data',ylab='Only counts from 10-12 traps'); abline(v=0); abline(h=0)

trends2 = trends2[,1:5]
write.table(trends2,'./summary_tables/reWelti/Phoenix_Pitfall1_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out2,'./summary_tables/reWelti/Phoenix_Pitfall1_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Central Arizona-Phoenix - Arthropod Pitfalls 2
# https://sustainability.asu.edu/caplter/data/view/knb-lter-cap.643.2/

data1 = read.csv('./raw_data/643_mcdowell_pitfall_arthropods_76003ef56b85e55df4e9f2605c3e8492.csv',as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$display_name
data1$Species = gsub('\\(adult\\)','',data1$Species)
data1$Species = gsub('\\(winged\\)','',data1$Species)
data1$Species = gsub('\\(larvae\\)','',data1$Species)
data1$Species = gsub('\\(nymph\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(Winged Aphelinidae\\)','',data1$Species)
data1$Species = gsub(' \\(Aphelinidae\\)','',data1$Species)
data1$Species = gsub('_',' ',data1$Species)
data1$Species = gsub(' \\(Centipeds\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('Scale Insects \\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!=''),]
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
plots = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$Site==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			n.traps = length(unique(data2$trap_name[which(data2$Year==years3[j])]))
			out[ox,1] = 'Phoenix' #LTER
			out[ox,2] = plots[i] #"Locale" (plot)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = n.traps #"N.Obs" (in this case, number of pitfall traps for which counts are reported at a site)
			out[ox,6] = sum(data4$Number,na.rm=T) / n.traps #metadata reports 3 sweep samples per plot
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Phoenix/pitfall2/',s,'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'Phoenix'
			trends[tx,2] = 'pitfall2'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plots',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/Phoenix_Pitfall2_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/Phoenix_Pitfall2_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Harvard Forest ants 1
# https://doi.org/10.6073/pasta/7a6b956fb0960d7fe8bb048b1fe26956

data1 = read.csv('./raw_data/hf118-01-ants.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Year = data1$year
data1$Number = data1$abundance
#fix method labels
data1$trap.type[which(data1$trap.type=='Pitfall' | data1$trap.type=='pitfall' | data1$trap.type=='pit')]='pitfall'
data1$trap.type[which(data1$trap.type=='Litter' | data1$trap.type=='LItter' | data1$trap.type=='litter')]='litter'
data1$trap.type[which(data1$trap.type=='Hand sample' | data1$trap.type=='hand')]='hand'
data1$trap.type[which(data1$trap.type=='Bait' | data1$trap.type=='Biat' | data1$trap.type=='bait')]='bait'
data1$block.plot.treat.moose.type = paste(data1$block,data1$plot,data1$treatment,data1$moose.cage,data1$trap.type,sep='_') #identifies unique factors associated with each count
plots = unique(data1$block.plot.treat.moose.type)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(plots)){
	data2 = data1[which(data1$block.plot.treat.moose.type==plots[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'HarvardForest' #LTER
			out[ox,2] = plots[i] #"Locale" (unique combination of plot, treatment, trap type)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (in this case, number of traps & dates in which a count is reported. Absences are true zeroes)
			out[ox,6] = sum(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/HarvardForest/ants1/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'HarvardForest'
			trends[tx,2] = 'ants1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plot*method',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/HarvardForest_Ants1_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/HarvardForest_Ants1_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Harvard Forest ants - Nantucket
# https://doi.org/10.6073/pasta/3493424abf9fc36eac7b62b732e4ea55 

data1 = read.csv('./raw_data/hf147-08-nantucket-ants-2004-09.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Site = data1$site
data1$Year = data1$year
data1$Number = data1$qty
data1 = data1[which(!is.na(data1$Species)),]
sites = unique(data1$Site)
apply(array(sites),1,function(x){unique(data1$pitfalltrap[which(data1$Site==x)])}) #5 sites where different trapping method is indicated
apply(array(sites),1,function(x){unique(data1$management[which(data1$Site==x)])}) #3 sites report >1 (2) managements
data1$Method = data1$pitfalltrap
data1$Site.Method = paste(data1$Site,data1$Method,sep='_')
site.methods = unique(data1$Site.Method)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(site.methods)){
	data2 = data1[which(data1$Site.Method==site.methods[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'HarvardForest' #LTER
			out[ox,2] = site.methods[i] #"Locale" (site*method)
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (in this case, number of quantities reported)
			out[ox,6] = mean(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/HarvardForest/ants.Nantucket/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'HarvardForest'
			trends[tx,2] = 'antsNantucket'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*plot*method',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/HarvardForest_AntsNantucket_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/HarvardForest_AntsNantucket_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Harvard Forest - Ticks
# https://doi.10.6073/pasta/b29a97941c11ddf45540ea30066fde35 -link broken. Use: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF299

data1 = read.csv('./raw_data/hf299-01-survey.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$year
data1$Site = data1$location.name
data1 = data1[which(!is.na(data1$Site)),]
sites = unique(data1$Site)
species = c('Dermacentor variabilis','Ixodes scapularis')
scols = match(c('wood.found','deer.found'),colnames(data1))
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Site==sites[i]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		for (s in 1:length(species)){
			tick.counts = data3[,scols[s]]
			person.hours = sum(data3$hours,na.rm=T)
			out[ox,1] = 'HarvardForest' #LTER
			out[ox,2] = sites[i] #"Locale" ("location.name", as recorded by participants)
			out[ox,3] = species[s]
			out[ox,4] = years2[j] #Year
			out[ox,5] = person.hours #"N.Obs" (in this case, number of person hours contributing to tick count)
			out[ox,6] = sum(tick.counts,na.rm=T) / person.hours
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/HarvardForest/ticks/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'HarvardForest'
			trends[tx,2] = 'ticks'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/HarvardForest_Ticks_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/HarvardForest_Ticks_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Hubbard Brook - Lepidoptera 1
# https://doi.org/10.6073/pasta/5d2a8c67c5a3278032b2b14d66c09a7f

data1 = read.csv('./raw_data/lepidoptera_HB_1986-2018.csv',as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
data1 = data1[which(data1$Species!=0),]
data1 = data1[which(!is.na(data1$Species)),]
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'HubbardBrook' #LTER
		out[ox,2] = 'HBEF' #"Locale" "Hubbard Brook Experimental Forest"
		out[ox,3] = species[s]
		out[ox,4] = years2[j] #Year
		out[ox,5] = nrow(data3) #"N.Obs" (number of rows in dataframe, almost always = number of individuals) Effort reported as constant in metadata
		out[ox,6] = sum(data3$Number,na.rm=T)
		ox = ox + 1
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/HubbardBrook/lepidoptera1/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'HubbardBrook'
			trends[tx,2] = 'lepidoptera1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/HubbardBrook_Lep1_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/HubbardBrook_Lep1_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Hubbard Brook - Lepidoptera 2
# https://doi.org/10.6073/pasta/5d2a8c67c5a3278032b2b14d66c09a7f

data1 = read.csv('./raw_data/lepidoptera_WMNFsites.csv',as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
sites = unique(data1$Plot)
apply(array(sites),1,function(x){length(which(data1$Plot==x))})
data1$Site.Year = paste(data1$Plot,data1$Year,sep='_')
apply(array(unique(data1$Site.Year)),1,function(x){unique(data1$Date[which(data1$Site.Year==x)])})
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Plot==sites[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			sample.months = length(unique(data1$Date[which(data1$Plot==sites[i] & data1$Year==years3[j])]))
			out[ox,1] = 'HubbardBrook' #LTER
			out[ox,2] = sites[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = sample.months #"N.Obs" (number of months during which samples are reported in this plot*year)
			out[ox,6] = sum(data4$Number,na.rm=T) / sample.months
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/HubbardBrook/lepidoptera2/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'HubbardBrook'
			trends[tx,2] = 'lepidoptera2'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/HubbardBrook_Lep2_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/HubbardBrook_Lep2_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Northern Temperate Lakes - Pelagic Macroinvertebrates
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.14.30

data1 = read.csv('./raw_data/ntl14_v8.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Species = data1$taxon
data1$Number = ceiling(data1$avg_ind_per_m3)
sites = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Site==sites[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'NorthTemperateLakes' #LTER
			out[ox,2] = sites[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (number of rows in dataframe)
			out[ox,6] = sum(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/NorthTemperateLakes/PelagicMacroinvertebrates/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'NorthTemperateLakes'
			trends[tx,2] = 'PelagicMacroinvertebrates'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/NorthTemperateLakes_Pelagic_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/NorthTemperateLakes_Pelagic_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Northern Temperate Lakes - Benthic Macroinvertebrates
# https://doi.org/10.6073/pasta/cc25694cdde49853271df465a15007fb

data1 = read.csv('./raw_data/ntl11_1_v8.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Site[which(data1$Site=='sp')] = 'SP'
data1$Species = data1$taxon_code
data1$Number = data1$number_indiv
data1$Lake.Site = paste(data1$Site,data1$site,sep='_')
sites = unique(data1$Lake.Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Lake.Site==sites[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'NorthTemperateLakes' #LTER
			out[ox,2] = sites[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = nrow(data4) #"N.Obs" (number of rows in dataframe)
			out[ox,6] = sum(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/NorthTemperateLakes/BenthicMacroinvertebrates/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'NorthTemperateLakes'
			trends[tx,2] = 'BenthicMacroinvertebrates'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/NorthTemperateLakes_Benthic_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/NorthTemperateLakes_Benthic_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Temperate Lakes - Crayfish 1
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-ntl&identifier=3&revision=27

data1 = read.csv('./raw_data/ntl3_v7.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$spname
data1$Site = data1$lakeid
data1$Year = data1$year4
data1$Number = data1$total_caught
data1 = data1[which(data1$Species!='CRAYFISH'),] #remove totals
data1$Site.Method = paste(data1$Site,data1$gearid,sep='_')
site.methods = unique(data1$Site.Method)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(site.methods)){
	data2 = data1[which(data1$Site.Method==site.methods[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'NorthTemperateLakes' #LTER
			out[ox,2] = site.methods[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = data4$effort #"N.Obs" ("effort" refers to the number of traps)
			out[ox,6] = data4$Number / data4$effort
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/NorthTemperateLakes/Crayfish1/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'NorthTemperateLakes'
			trends[tx,2] = 'Crayfish1'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site*method',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/NorthTemperateLakes_Crayfish1_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/NorthTemperateLakes_Crayfish1_abundance.txt',sep='\t',quote=F,row.names=F)


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
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/NorthTemperateLakes/Crayfish2/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
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
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site*method',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/NorthTemperateLakes_Crayfish2_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/NorthTemperateLakes_Crayfish2_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Sevilleta - Grasshoppers
# https://doi.org/10.6073/pasta/c1d40e9d0ec610bb74d02741e9d22576

data1 = read.csv('./raw_data/sev106_hopperdynamics_20150826.txt',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$DATE),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Site = data1$SITE
data1$Species = data1$SPECIES
data1$Number = data1$CNT
sites = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Site==sites[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'Sevilleta' #LTER
			out[ox,2] = sites[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = length(unique(data4$DATE)) #"N.Obs" (number of dates with counts reported) Should = 2, but sometimes only 1 date is reported
			out[ox,6] = sum(data4$Number,na.rm=T)
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Sevilleta/grasshoppers/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'Sevilleta'
			trends[tx,2] = 'grasshoppers'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/Sevilleta_Grasshoppers_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/Sevilleta_Grasshoppers_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Sevilleta - Arthropod pitfalls
# https://doi.org/10.6073/pasta/9e7e6dc9c9d8f72e9e0bca07a1e76ccd

con = file('./raw_data/sev029_arthropop_02162009_0.txt','r')
lines1 = readLines(con)
close(con)
data1 = strsplit(lines1[1],',')[[1]]
for (i in 2:length(lines1)){
	add = gsub('-888',NA,gsub('na',NA,strsplit(lines1[i],',')[[1]]))
	data1 = data.frame(rbind(data1,add),stringsAsFactors=F)
}
colnames(data1) = data1[1,]; data1 = data1[-1,]; str(data1)
data1 = data1[which(data1$Year>1994),] #only use data after 1994 because the number of pitfall traps was reduced by half after 1994.
data1 = data1[which(data1$Year!=2002),] #remove single record from 2002, because it reports a count from only one trap (should be 15 traps)
data1$Species = paste(data1$Genus,data1$Species,sep='_')
data1$Number = as.numeric(data1$Count)
data1$Line.Trap = paste(data1$Line,data1$Trap,sep='_')
data1 = data1[which(data1$Species!='_' & data1$Species!='NA_NA'),]
sites = unique(data1$Site)
out = data.frame('LTER.site'=NA,'Locale'=NA,'Species.code'=NA,'Year'=-999,'N.Obs'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (i in 1:length(sites)){
	data2 = data1[which(data1$Site==sites[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			n.traps = length(unique(data1$Line.Trap[which(data1$Site==sites[i] & data1$Year==years3[j])]))
			out[ox,1] = 'Sevilleta' #LTER
			out[ox,2] = sites[i] #"Locale"
			out[ox,3] = species2[s]
			out[ox,4] = years3[j] #Year
			out[ox,5] = n.traps #"N.Obs" (number of traps)
			out[ox,6] = sum(data4$Number,na.rm=T) / n.traps
			ox = ox + 1
		}
	}
}; str(out)

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species==species[s]),]
	locales = sort(unique(dat$Locale))
	for (l in 1:length(locales)){
		dat2 = dat[which(dat$Locale==locales[l]),]
		if ( (length(which(!is.na(dat2$Abundance))) > 3) & (length(which(dat2$Abundance[!is.na(dat2$Abundance)] > 0)) > 0) ){ #greater than 3 non-missing abundance values, at least one non-zero
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
				X[which(!is.na(X))] = dat2$Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
				png(paste0('./plots/Welti reanalysis/Sevilleta/pitfall/',species[s],'_',locales[l],'.png'))
				plot(t.scale,Z,main=paste0(species[s],' ',locales[l]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
				dev.off()
			}
			trends[tx,1] = 'Sevilleta'
			trends[tx,2] = 'pitfall'
			trends[tx,3] = locales[l]
			trends[tx,4] = species[s]
			trends[tx,5] = abundance.trend
			tx = tx + 1
		} else { #time series did not meet filtering criteria
			lost = data.frame(rbind(lost,data.frame('Locale'=locales[l],'Species'=species[s])),stringsAsFactors=F)	
		}	
	}
}; str(trends)

hist(trends$Abundance.trend,xlab='Abundance trend (sd/yr)',ylab='No. species*site',cex.lab=1.5,cex.axis=1.2,col='grey40',main='')
abline(v=0,col='pink',lwd=2)

write.table(trends,'./summary_tables/reWelti/Sevilleta_Pitfall_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/reWelti/Sevilleta_Pitfall_abundance.txt',sep='\t',quote=F,row.names=F)


######################################################################################################
# Compile all trends into a single file

files = list.files('./summary_tables/reWelti',full.names=T)
files = files[grep('_trends.txt',files)]
trends = read.table(files[1],sep='\t',as.is=T,check.names=F,header=T)
for (i in 2:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)
	trends = data.frame(rbind(trends,add),stringsAsFactors=F)
}
# Relabel datasets
trends$Dataset[which(trends$Dataset=='grasshoppers')] = 'Grasshoppers'
trends$Dataset[which(trends$Dataset=='sweep1' | trends$Dataset=='sweep2')] = 'Sweeps'
trends$Dataset[which(trends$Dataset=='ants1' | trends$Dataset=='antsNantucket')] = 'Ants'
trends$Dataset[which(trends$Dataset=='ticks')] = 'Ticks'
trends$Dataset[which(trends$Dataset=='lepidoptera1' | trends$Dataset=='lepidoptera2')] = 'Leptidoptera'
trends$Dataset[which(trends$Dataset=='BenthicMacroinvertebrates' | trends$Dataset=='PelagicMacroinvertebrates')] = 'Macroinvertebrates'
trends$Dataset[which(trends$Dataset=='Crayfish1' | trends$Dataset=='Crayfish2')] = 'Crayfish'
trends$Dataset[which(trends$Dataset=='pitfall1' | trends$Dataset=='pitfall2' | trends$Dataset=='pitfall')] = 'Pitfalls'
trends$Dataset[which(trends$Dataset=='sweep')] = 'Sweep'
# Add original trends for sites not included in critique
add2 = read.csv('./summary_tables/reWelti/time_trends_arthropods_relaxed.csv',as.is=T,check.names=F,header=T)
keep = unlist(apply(array(c('Arctic','Baltimore','BonanzaCreek','Coweeta','GeorgiaCoastal','MidwestSTN')),1,function(x){which(add2$LTER==x)}))
add2 = add2[keep,]
add3 = data.frame('LTER'=add2$LTER,'Dataset'=NA,'Locale'=add2$Site,'Species'=add2$Species,'Abundance.trend'=add2$coef_slope,stringsAsFactors=F)
# Relabel datasets
add3$Dataset[which(add3$LTER=='Arctic')] = 'StreamInsects'
add3$Dataset[which(add3$LTER=='Baltimore')] = 'Mosquitoes'
add3$Dataset[which(add3$LTER=='BonanzaCreek' & add3$Locale=='bark beetle')] = 'BarkBeetles'
add3$Dataset[which(add3$LTER=='BonanzaCreek' & add3$Locale=='leafminer')] = 'Leafminers'
add3$Dataset[which(add3$Locale=='Burrowing Crab' | add3$Locale=='Fiddler Crab')] = 'Crabs'
add3$Dataset[which(add3$LTER=='GeorgiaCoastal' & add3$Locale=='grasshopper')] = 'Grasshoppers'
add3$Dataset[which(add3$LTER=='GeorgiaCoastal' & add3$Locale=='Salt Marsh' & add3$Species=='Prokelisia marginata')] = 'PlantHoppers'
add3$Dataset[which(add3$LTER=='GeorgiaCoastal' & add3$Locale=='Salt Marsh' & add3$Species=='Acrididae')] = 'Grasshoppers'
add3$Dataset[which(add3$LTER=='MidwestSTN')] = 'Aphids'
add3$Dataset[which(add3$LTER=='Coweeta')] = 'StreamInsects'
trends2 = data.frame(rbind(trends,add3),stringsAsFactors=F)
write.csv(trends2,'./time_trends_arthropods_reWelti.csv',row.names=F,quote=F)

# Compile abundances into a single file
files = list.files('./summary_tables/reWelti',full.names=T)
files = files[grep('_abundance.txt',files)]
abundances = read.table(files[1],sep='\t',as.is=T,check.names=F,header=T)
for (i in 2:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)
	abundances = data.frame(rbind(abundances,add),stringsAsFactors=F)
}
add2 = read.csv('PerSpecies_Abundance_LTER.csv',as.is=T,check.names=F,header=T)
keep = unlist(apply(array(c('Arctic','Baltimore','BonanzaCreek','Coweeta','GeorgiaCoastal','MidwestSTN')),1,function(x){which(add2$LTER==x)}))
add2 = add2[keep,]
abundances2 = data.frame(rbind(abundances,add2),stringsAsFactors=F)
abundances2 = abundances[,c(1:4,6)]
write.csv(abundances2,'./PerSpecies_Abundance_LTER_reWelti.csv',row.names=F,quote=F)


######################################################################################################
# Recreate violin plot (Figure 2) using updated trends

library(ggplot2)

trends = read.csv('./time_trends_arthropods_reWelti.csv',as.is=T,check.names=F,header=T)
trends$LL = paste(trends$LTER,trends$Dataset)
levels1 = c('Arctic StreamInsects','Coweeta StreamInsects','NorthTemperateLakes Crayfish','NorthTemperateLakes Macroinvertebrates','BonanzaCreek BarkBeetles','BonanzaCreek Leafminers','CedarCreek Grasshoppers',
	'CedarCreek Sweeps','GeorgiaCoastal Crabs','GeorgiaCoastal Grasshoppers','GeorgiaCoastal PlantHoppers','HarvardForest Ants','HarvardForest Ticks','HubbardBrook Leptidoptera','KonzaPrairie GallInsects',
	'KonzaPrairie Grasshoppers','Sevilleta Grasshoppers','Sevilleta Pitfalls','Baltimore Mosquitoes','MidwestSTN Aphids','Phoenix Pitfalls','Phoenix Sweep')
trends$LL = factor(trends$LL,levels=levels1)
trends$Habitat = 'terrestrial'
trends$Habitat[which(trends$LTER=='Arctic' | trends$LTER=='Coweeta' | trends$LTER=='NorthTemperateLakes')] = 'aquatic'
trends$Habitat[which(trends$LTER=='MidwestSTN' | trends$LTER=='Baltimore' | trends$LTER=='Phoenix')] = 'impacted'
trends$Habitat = factor(trends$Habitat,levels=c('aquatic','terrestrial','impacted'))

dp = ggplot(trends, aes(x=LL, y=Abundance.trend, fill=Habitat)) + 
	geom_violin(trim=T,size=1) +
	scale_fill_manual(values=c('lightblue','white','darkorange')) +
	ylim(-6,6) +
	stat_summary(fun.y=median, geom="point", shape=18, size=3, color="black") +
	geom_hline(yintercept=0,size=1) +
	labs(y="Abundance trend") +
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
	ggsave(filename='./plots/Welti Reanalysis/ALL_time_trends_boxplot_relaxed_noannotations_violin.png',plot=dp,dpi=600,unit='cm',height=8,width=24)

tcount = apply(array(levels(trends$LL)),1,function(x){length(which(as.character(trends$LL)==x))})
cbind(levels(trends$LL),tcount)


######################################################################################################
# Mean & 95% CI

trends = read.csv('./time_trends_arthropods_reWelti.csv',as.is=T,check.names=F,header=T)

# t-test: Is the mean of LTER mean trends different from zero?
uL = unique(trends$LTER)
means2 = apply(array(uL),1,function(x){mean(trends$Abundance.trend[which(trends$LTER==x)],na.rm=T)})
t2 = t.test(as.numeric(means2))

# Barplot
mean2 = t2[['estimate']]
uppers = t2[['conf.int']][2]
lowers = t2[['conf.int']][1]
png('./plots/Welti Reanalysis/t-test_abundance_lters_barplot_reWelti.png',width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mean2,ylim=c(-1.5,1.5),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
abline(h=0)
axis(2,lwd=2)
text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2)
arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/15)
dev.off()

