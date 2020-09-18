

setwd()
source('./code/AR_reml.R')

th = 3 #require time series to include at least 4 years of non-missing data

#############################################################################################
# Konza Prairie grasshoppers
# https://doi.org/10.6073/pasta/7b2259dcb0e499447e0e11dfb562dc2f

data1 = read.csv('./raw_data/CGR022.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$SPECIES
data1$Year = data1$RECYEAR
data1$Number = as.numeric(data1$TOTAL)
data1 = data1[which(data1$Species!='unknown' & data1$Species!='unknown ' & data1$Species!=''),]
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
data1$WATERSHED[which(data1$WATERSHED=='004d')] = '004D'
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
data1. = data1[which(data1$Year>=1996),]
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
years = sort(unique(data1.$Year))
for (j in 1:length(years)){
	data2 = data1.[which(data1.$Year==years[j] & data1.$RECMONTH!=6),]
	species2 = sort(unique(data2$Species))
	days = paste(data2$RECMONTH,data2$RECDAY,sep='_')
	for (k in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[k]),]
		out[ox,1] = 'KonzaPrairie' #LTER
		out[ox,2] = species2[k] #Species.code (in this case, actual species name)
		out[ox,3] = years[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	if ( (length(which(!is.na(dat$Abundance))) > th) & (length(which(dat$Abundance[!is.na(dat$Abundance)] > 0)) > 0) ){ #greater than th non-missing abundance values, at least one non-zero
		cmin = min(dat$Abundance[which(dat$Abundance>0)],na.rm=T)
		if (is.na(cmin)){
			dat$Abundance[which(dat$Abundance==0)] = 0.5 
		} else {
			dat$Abundance[which(dat$Abundance==0)] = 0.5 * cmin #replace zeroes with 0.5*minimum abundance value in this time series
		}
		dat$log.Abundance = log(dat$Abundance) #log-transform abundances
		if (length(unique(dat$log.Abundance[which(!is.na(dat$log.Abundance))]))==1){ #abundances are all = 1
			abundance.trend = 0
		} else {
			ys = sort(dat$Year)
			ys.stretch = seq(ys[1],ys[length(ys)],1) #expand time series to include missing years
			X = match(ys.stretch,ys)
			X[which(!is.na(X))] = dat$log.Abundance[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
			Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
			t.scale = 1:length(ys.stretch)
			t.scale = (t.scale-min(t.scale))/max(t.scale)
			if (nrow(dat)==2){
				abundance.trend = (Z[2] - Z[1]) / 2 #simple slope for 2-yr time series
			} else {
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				abundance.trend = arr.Z$coef[2]
			}
			png(paste0('./plots/m4_standard/KonzaGrasshoppers/',species[s],'.png'))
			plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
			abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
			dev.off()
		}
		trends[tx,1] = 'KonzaPrairie'
		trends[tx,2] = 'Grasshoppers'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Konza_Grasshoppers_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Konza_Grasshoppers_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Konza Prairie gall insects
# http://dx.doi.org/10.6073/pasta/b2ac9e918a66dbbb18c7a6b39dc1efab

data1 = read.csv('./raw_data/CGP011.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$RecYear
data1$Number = data1$GalledStems / data1$SampledStems #galls per stem
# Fix plant species name typoes
data1$Species[which(data1$Species=='ceanothu')] = 'ceanothus'
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'KonzaPrairie' #LTER
		out[ox,2] = unique(data3$Species)
		out[ox,3] = years2[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = unique(out$Species.code)
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#				png(paste0('./plots/Konza/galls/',locales[l],'.png'))
#				plot(t.scale,Z,main=locales[l],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#				dev.off()
		}
		trends[tx,1] = 'KonzaPrairie'
		trends[tx,2] = 'GallInsects'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}
}
write.table(trends,paste0('./summary_tables/m4_standard/Konza_GallInsects_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Konza_GallInsects_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################5
# Cedar Creek Ecosystem - Old Field Grasshopper Sampling
# https://www.cedarcreek.umn.edu/research/data/methods?e014

data1 = read.table('./raw_data/e014_Core Old Field Grasshopper Sampling.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1$Number = data1$Specimens
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
species = unique(data1$Species)
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'CedarCreek' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years2[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/CedarCreek/grasshoppers/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'CedarCreek'
		trends[tx,2] = 'Grasshoppers'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/CedarCreek_Grasshoppers_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/CedarCreek_Grasshoppers_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Cedar Creek Ecosystem Arthropods Sweep 1
# https://www.cedarcreek.umn.edu/research/data/dataset?arce153

data1 = read.table('./raw_data/e153_Arthropod sweepnet sampling.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='undet under'),]
data1$Number = data1$Specimens
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'CedarCreek' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years2[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/CedarCreek/sweep1/',s,'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'CedarCreek'
		trends[tx,2] = 'sweep1'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/CedarCreek_Sweep1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/CedarCreek_Sweep1_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Cedar Creek Ecosystem Arthropods Sweep 2
# https://www.cedarcreek.umn.edu/research/data/dataset?aage120

data1 = read.table('./raw_data/e120_Main Plots All Arthropod Insect Sweepnet Sampling 1996-2006.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Species = paste(data1$Genus,data1$Specific.epithet,sep=' ')
data1 = data1[which(data1$Species!='undet_undet'),]
data1$Number = as.numeric(data1$Count)
data1 = data1[which(data1$Species!='undet undet' & data1$Species!='unk unk' & data1$Species!='none none' & data1$Species!='na? na?' & data1$Species!='na na'),]
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'CedarCreek' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years2[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/CedarCreek/sweep2/',s,'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'CedarCreek'
		trends[tx,2] = 'sweep2'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/CedarCreek_Sweep2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/CedarCreek_Sweep2_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Central Arizona-Phoenix - Arthropod sweeps
# https://doi.org/10.6073/pasta/0669ee6a71b24abb1ae3827f4ee77f6d

data1 = read.csv('./raw_data/652_arthropods_e9c22403e0f7243f241ed64862e22e05.csv',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$sample_date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Site = data1$site_code
data1$Species = data1$arthropod_scientific_name
data1$Species = gsub('\\(immature\\)','',data1$Species)
data1$Species = gsub(' \\(immature\\)','',data1$Species)
data1$Species = gsub('\\(Coccoidea\\)','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
data1$Number = data1$number_of_arthropods
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years = sort(unique(data2$Year))
	for (j in 1:length(years)){
		data3 = data2[which(data2$Year==years[j]),]
		out[ox,1] = 'Phoenix' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = mean(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Phoenix/sweep/',s,'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Phoenix'
		trends[tx,2] = 'Sweep'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Phoenix_Sweep_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Phoenix_Sweep_abundance.txt'),sep='\t',quote=F,row.names=F)


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
data1$Species = gsub('\\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('\\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(data1$Number>=0),] #remove record = -2
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years = sort(unique(data2$Year))
	for (j in 1:length(years)){
		data3 = data2[which(data2$Year==years[j]),]
		d0 = data1[which(data1$Year==years[j]),]
		ds = unique(d0$Site)
		ds.traps = apply(array(ds),1,function(x){length(unique(d0$trap_name[which(d0$Site==x)]))})
		n.traps = sum(ds.traps)
		out[ox,1] = 'Phoenix' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n.traps #metadata reports 3 sweep samples per plot
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#				png(paste0('./plots/Phoenix/pitfall1/',s,'.png'))
#				plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#				abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#				dev.off()
		}
		trends[tx,1] = 'Phoenix'
		trends[tx,2] = 'pitfall1'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Phoenix_Pitfall1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Phoenix_Pitfall1_abundance.txt'),sep='\t',quote=F,row.names=F)


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
data1$Species = gsub('\\(Coccoidea\\)','',data1$Species)
data1$Species = gsub('\\(Coccoidea - Immature\\)','',data1$Species)
data1$Species = gsub(' \\(imm.\\)','',data1$Species)
data1$Species = gsub(' >10mm','',data1$Species)
data1$Species = gsub(' immature','',data1$Species)
data1$Species = trimws(data1$Species,which='right')
data1 = data1[which(data1$Species!=''),]
data1$Number = apply(data1[,13:17],1,function(x){sum(x,na.rm=T)})
data1 = data1[which(!is.na(data1$Species)),]
data1 = data1[which(data1$Species!='' & data1$Species!='Unknown'),]
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years = sort(unique(data2$Year))
	for (j in 1:length(years)){
		data3 = data2[which(data2$Year==years[j]),]
		d0 = data1[which(data1$Year==years[j]),]
		ds = unique(d0$Site)
		ds.traps = apply(array(ds),1,function(x){length(unique(d0$trap_name[which(d0$Site==x)]))})
		n.traps = sum(ds.traps)
		out[ox,1] = 'Phoenix' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n.traps #metadata reports 3 sweep samples per plot
		ox = ox + 1
	}
}


# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Phoenix/pitfall2/',s,'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Phoenix'
		trends[tx,2] = 'pitfall2'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Phoenix_Pitfall2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Phoenix_Pitfall2_abundance.txt'),sep='\t',quote=F,row.names=F)


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
traps = unique(data1$trap.type)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,'Trap.type'=NA,stringsAsFactors=F)
ox = 1
for (i in 1:length(traps)){
	data2 = data1[which(data1$trap.type==traps[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			bptm = paste(data4$block,data4$plot,data4$treatment,data4$moose.cage,sep='_')
			n.trts = length(unique(bptm))
			out[ox,1] = 'HarvardForest' #LTER
			out[ox,2] = species2[s]
			out[ox,3] = years3[j] #Year
			out[ox,4] = sum(data4$Number,na.rm=T) / n.trts
			out[ox,5] = traps[i]
			ox = ox + 1
		}
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999,'Trap.type'=NA)
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	traps = sort(unique(dat$Trap.type))
	for (l in 1:length(traps)){
		dat2 = dat[which(dat$Trap.type==traps[l]),]
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
#					png(paste0('./plots/HarvardForest/ants1/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'HarvardForest'
			trends[tx,2] = 'Ants'
			trends[tx,3] = species[s]
			trends[tx,4] = abundance.trend
			trends[tx,5] = traps[l]
			tx = tx + 1
		} else { #time series did not meet filtering criteria
		}	
	}
}
write.table(trends,paste0('./summary_tables/m4_standard/HarvardForest_Ants1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/HarvardForest_Ants1_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Harvard Forest ants - Nantucket
# https://doi.org/10.6073/pasta/3493424abf9fc36eac7b62b732e4ea55 

data1 = read.csv('./raw_data/hf147-08-nantucket-ants-2004-09.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$code
data1$Site = data1$site
data1$Year = data1$year
data1$Number = data1$qty
data1 = data1[which(!is.na(data1$Species)),]
data1$Method = data1$pitfalltrap
data1$Method[which(data1$Method!='Pitfall' & data1$Method!='Tuna' & data1$Method!='Hand' & data1$Method!='Litter')] = 'Pitfall'
traps = unique(data1$Method)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,'Trap.type'=NA,stringsAsFactors=F)
ox = 1
for (i in 1:length(traps)){
	data2 = data1[which(data1$Method==traps[i]),]
	species2 = unique(data2$Species)
	for (s in 1:length(species2)){
		data3 = data2[which(data2$Species==species2[s]),]
		years3 = sort(unique(data3$Year))
		for (j in 1:length(years3)){
			data4 = data3[which(data3$Year==years3[j]),]
			out[ox,1] = 'HarvardForest' #LTER
			out[ox,2] = species2[s]
			out[ox,3] = years3[j] #Year
			out[ox,4] = mean(data4$Number,na.rm=T)
			out[ox,5] = traps[i]
			ox = ox + 1
		}
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999,'Trap.type'=NA)
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	traps = sort(unique(dat$Trap.type))
	for (l in 1:length(traps)){
		dat2 = dat[which(dat$Trap.type==traps[l]),]
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
#					png(paste0('./plots/HarvardForest/ants.Nantucket/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'HarvardForest'
			trends[tx,2] = 'AntsNantucket'
			trends[tx,3] = species[s]
			trends[tx,4] = abundance.trend
			trends[tx,5] = traps[l]
			tx = tx + 1
		} else { #time series did not meet filtering criteria
		}	
	}
}
write.table(trends,paste0('./summary_tables/m4_standard/HarvardForest_AntsNantucket_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/HarvardForest_AntsNantucket_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Harvard Forest - Ticks
# https://doi.10.6073/pasta/b29a97941c11ddf45540ea30066fde35 -link broken. Use: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF299

data1 = read.csv('./raw_data/hf299-01-survey.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$year
data1$Site = gsub("'",'',data1$location.name)
data1 = data1[which(!is.na(data1$Site)),]
species = c('Dermacentor variabilis','Ixodes scapularis')
scols = match(c('wood.found','deer.found'),colnames(data1))
years = sort(unique(data1$Year))
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (j in 1:length(years)){
	data2 = data1[which(data1$Year==years[j]),]
	for (s in 1:length(species)){
		tick.counts = data2[,scols[s]]
		person.hours = sum(data2$hours,na.rm=T)
		out[ox,1] = 'HarvardForest' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = sum(tick.counts,na.rm=T) / person.hours
		ox = ox + 1
	}
}


# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/HarvardForest/ticks/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'HarvardForest'
		trends[tx,2] = 'Ticks'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/HarvardForest_Ticks_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/HarvardForest_Ticks_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Hubbard Brook - Lepidoptera 1
# https://doi.org/10.6073/pasta/5d2a8c67c5a3278032b2b14d66c09a7f

data1 = read.csv('./raw_data/lepidoptera_HB_1986-2018.csv',as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
data1 = data1[which(data1$Species!=0),]
data1 = data1[which(!is.na(data1$Species)),]
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years2 = sort(unique(data2$Year))
	for (j in 1:length(years2)){
		data3 = data2[which(data2$Year==years2[j]),]
		out[ox,1] = 'HubbardBrook' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years2[j] #Year
		out[ox,4] = sum(data3$Number,na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/HubbardBrook/lepidoptera1/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'HubbardBrook'
		trends[tx,2] = 'lepidoptera1'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/HubbardBrook_Lep1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/HubbardBrook_Lep1_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Hubbard Brook - Lepidoptera 2
# https://doi.org/10.6073/pasta/5d2a8c67c5a3278032b2b14d66c09a7f

data1 = read.csv('./raw_data/lepidoptera_WMNFsites.csv',as.is=T,check.names=F,header=T)
data1$Number = data1$NumberIndividuals
data1$Species = data1$Taxon
species = unique(data1$Species)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	data2 = data1[which(data1$Species==species[s]),]
	years = sort(unique(data2$Year))
	for (j in 1:length(years)){
		data3 = data2[which(data2$Year==years[j]),]
		d0 = data1[which(data1$Year==years[j]),]
		ps = paste(d0$Plot,d0$Date,sep='_')
		ps.count = apply(array(unique(ps)),1,function(x){length(which(ps==x))})
		sample.months = sum(ps.count,na.rm=T)
		out[ox,1] = 'HubbardBrook' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / sample.months
		ox = ox + 1
	}
}


# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/HubbardBrook/lepidoptera2/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'HubbardBrook'
		trends[tx,2] = 'lepidoptera2'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/HubbardBrook_Lep2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/HubbardBrook_Lep2_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Northern Temperate Lakes - Pelagic Macroinvertebrates
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.14.30

data1 = read.csv('./raw_data/ntl14_v8.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$year4
data1$Site = data1$lakeid
data1$Species = data1$taxon
data1$Number = ceiling(data1$avg_ind_per_m3)
years = sort(unique(data1$Year))
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	species = unique(data2$Species)
	lds = paste(data2$Site,data2$sta,data2$depth,sep='_')
	n = length(unique(lds))
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		out[ox,1] = 'NorthTemperateLakes' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n
		ox = ox + 1
	}
}


# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/NorthTemperateLakes/PelagicMacroinvertebrates/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'NorthTemperateLakes'
		trends[tx,2] = 'PelagicMacroinvertebrates'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Pelagic_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Pelagic_abundance.txt'),sep='\t',quote=F,row.names=F)


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
years = sort(unique(data1$Year))
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	species = unique(data2$Species)
	lds = paste(data2$Lake.Site,sep='_')
	n = length(unique(lds))
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		out[ox,1] = 'NorthTemperateLakes' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/NorthTemperateLakes/BenthicMacroinvertebrates/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'NorthTemperateLakes'
		trends[tx,2] = 'BenthicMacroinvertebrates'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Benthic_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Benthic_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Temperate Lakes - Crayfish 1
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-ntl&identifier=3&revision=27

data1 = read.csv('./raw_data/ntl3_v7.csv',as.is=T,check.names=F,header=T)
data1$Species = data1$spname
data1$Species[which(data1$Species=='PROPINQUUS')] = 'Orconectes propinquus'
data1$Species[which(data1$Species=='RUSTICUS')] = 'Orconectes rusticus'
data1$Species[which(data1$Species=='VIRILIS')] = 'Orconectes virilis'
data1$Site = data1$lakeid
data1$Year = data1$year4
data1$Number = data1$total_caught
data1 = data1[which(data1$Species!='CRAYFISH'),] #remove totals
traps = unique(data1$gearid)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,'Trap.type'=NA,stringsAsFactors=F)
ox = 1
for (i in 1:length(traps)){
	data2 = data1[which(data1$gearid==traps[i]),]
	years = unique(data2$Year)
	for (y in 1:length(years)){
		data3 = data2[which(data2$Year==years[y]),]
		species = unique(data3$Species)
		for (s in 1:length(species)){
			data4 = data3[which(data3$Species==species[s]),]
			out[ox,1] = 'NorthTemperateLakes' #LTER
			out[ox,2] = species[s]
			out[ox,3] = years[y] #Year
			out[ox,4] = mean(data4$Number / data4$effort,na.rm=T)
			out[ox,5] = traps[i]
			ox = ox + 1
		}
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999,'Trap.type'=NA)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat = out[which(out$Species.code==species[s]),]
	traps = sort(unique(dat$Trap.type))
	for (l in 1:length(traps)){
		dat2 = dat[which(dat$Trap.type==traps[l]),]
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
#					png(paste0('./plots/NorthTemperateLakes/Crayfish1/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
			}
			trends[tx,1] = 'NorthTemperateLakes'
			trends[tx,2] = 'Crayfish1'
			trends[tx,3] = species[s]
			trends[tx,4] = abundance.trend
			trends[tx,5] = traps[l]
			tx = tx + 1
		} else { #time series did not meet filtering criteria
		}	
	}
}
write.table(trends,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Crayfish1_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Crayfish1_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Temperate Lakes - Crayfish 2
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-ntl.269.2

data1 = read.csv('./raw_data/sparkling_crayfish.csv',as.is=T,check.names=F,header=T)
data1$Year = data1$YEAR4
species = c('Orconectes rusticus','Orconectes virilis')
spcol = match(c('Rusty CPUE','Virilis CPUE'),colnames(data1))
years = sort(unique(data1$Year))
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (s in 1:length(species)){
	for (j in 1:length(years)){
		data2 = data1[which(data1$Year==years[j]),]
		out[ox,1] = 'NorthTemperateLakes' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = mean(data2[,spcol[s]],na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/NorthTemperateLakes/Crayfish2/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'NorthTemperateLakes'
		trends[tx,2] = 'Crayfish2'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Crayfish2_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/NorthTemperateLakes_Crayfish2_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Sevilleta - Grasshoppers
# https://doi.org/10.6073/pasta/c1d40e9d0ec610bb74d02741e9d22576

data1 = read.csv('./raw_data/sev106_hopperdynamics_20150826.txt',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$DATE),1,function(x){strsplit(x,'/')[[1]][3]})
data1$Site = data1$SITE
data1$Species = data1$SPECIES
data1$Number = data1$CNT
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	species = unique(data2$Species)
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		ds = unique(paste(data3$DATE,data3$SITE,sep='_'))
		n = length(ds)
		out[ox,1] = 'Sevilleta' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data4$Number,na.rm=T) / n
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
lost = c()
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Sevilleta/grasshoppers/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Sevilleta'
		trends[tx,2] = 'Grasshoppers'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Sevilleta_Grasshoppers_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Sevilleta_Grasshoppers_abundance.txt'),sep='\t',quote=F,row.names=F)


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
colnames(data1) = data1[1,]; data1 = data1[-1,]
data1 = data1[which(data1$Year>1994),] #only use data after 1994 because the number of pitfall traps was reduced by half after 1994.
data1 = data1[which(data1$Year!=2002),] #remove single record from 2002, because it reports a count from only one trap (should be 15 traps)
data1$Species = paste(data1$Genus,data1$Species,sep='_')
data1$Number = as.numeric(data1$Count)
data1$Line.Trap = paste(data1$Line,data1$Trap,sep='_')
data1 = data1[which(data1$Species!='_' & data1$Species!='NA_NA'),]
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	slt = paste(data2$Site,data2$Line.Trap,sep='_')
	n.traps = length(unique(slt))
	species = unique(data2$Species)
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		out[ox,1] = 'Sevilleta' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n.traps
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Sevilleta/pitfall/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Sevilleta'
		trends[tx,2] = 'Pitfall'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,'./summary_tables/m4_standard/Sevilleta_Pitfall_trends.txt',sep='\t',quote=F,row.names=F)
write.table(out,'./summary_tables/m4_standard/Sevilleta_Pitfall_abundance.txt',sep='\t',quote=F,row.names=F)


#############################################################################################
# Arctic - stream insects
# https://arc-lter.ecosystems.mbl.edu/84-98hektot

data1 = read.csv('./raw_data/84-98hektot.csv',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){paste0(19,strsplit(x,'-')[[1]][3])})
data1 = data1[,which(colnames(data1)!='SNAILS')]
# Ensure consistent naming of stations
data1$Site = data1$`station (m)`
stations = unique(data1$Site)
data1$Site[which(data1$Site=='0.34K')] = 340
data1$Site[which(data1$Site=='3.0 K' | data1$Site=='3.0 K ' | data1$Site=='3.0K' | data1$Site=='3.0k')] = 3000
data1$Site[which(data1$Site=='2.05K' | data1$Site=='2.05 K ' | data1$Site=='2.05 K')] = 2050
data1$Site[which(data1$Site=='2.0k')] = 2000
data1$Site[which(data1$Site=='1.0K')] = 1000
data1$Site[which(data1$Site=='2.1K')] = 2100
data1$Site[which(data1$Site=='2.5k' | data1$Site=='2.5K')] = 2500
data1$Site[which(data1$Site=='1.6k')] = 1600
data1$Site[which(data1$Site=='2.6 K' | data1$Site=='2.6K')] = 2600
data1$Site[which(data1$Site=='1.4 K' | data1$Site=='1.4K')] = 1400
data1$Site[which(data1$Site=='1.2 K' | data1$Site=='1.2K')] = 1200
data1$Site[which(data1$Site=='.74 K' | data1$Site=='0.74K')] = 740
data1$Site[which(data1$Site=='.59K' | data1$Site=='0.59K')] = 590
data1$Site[which(data1$Site=='.14K' | data1$Site=='0.14K')] = 140
data1$Site[which(data1$Site=='.34K')] = 340
data1$Site[which(data1$Site=='.59K')] = 590
species = colnames(data1)[7:20]
years = unique(data1$Year)
spcol = 7:20
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	st = paste(data2$Site,data2$Trial,sep='_')
	n.trials = length(unique(st))
	for (s in 1:length(species)){
		out[ox,1] = 'Arctic' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data2[,spcol[s]],na.rm=T) / n.trials
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Arctic/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Arctic'
		trends[tx,2] = 'StreamInsects'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Arctic_StreamInsects_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Arctic_StreamInsects_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Baltimore - mosquitoes
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-bes.3500.100

data1 = read.csv('./raw_data/Biodiversity_-_Mosquito_-_ovitrap_mosquito_data.csv',as.is=T,check.names=F,header=T)
data1 = data1[which(!is.na(data1$site)),]
data1$Year = data1$year
data1$Site = data1$site
species = colnames(data1)[6:22]
years = unique(data1$Year)
spcol = 6:22
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	sw = paste(data2$Site,data2$week,sep='_')
	n.weeks = length(unique(sw))
	for (s in 1:length(species)){
		out[ox,1] = 'Baltimore' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data2[,spcol[s]],na.rm=T) / n.weeks
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Baltimore/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Baltimore'
		trends[tx,2] = 'Mosquitoes'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Baltimore_Mosquitoes_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Baltimore_Mosquitoes_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Bonanza Creek - Bark Beetles
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-bnz&identifier=35

data1 = read.csv('./raw_data/35_BNZ_Beetles_Werner_1975-2012.txt',as.is=T,check.names=F,header=T)
species = c('Spruce beetle','Ips beetle','Larch beetle')
spcol = 2:4
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
years = sort(unique(data1$Year))
for (s in 1:length(species)){
	for (j in 1:length(years)){
		data2 = data1[which(data1$Year==years[j]),]
		out[ox,1] = 'BonanzaCreek' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[j] #Year
		out[ox,4] = mean(data2[,spcol[s]],na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/BonanzaCreek/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'BonanzaCreek'
		trends[tx,2] = 'BarkBeetles'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/BonanzaCreek_BarkBeetles_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/BonanzaCreek_BarkBeetles_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Bonanza Creek aspen leaf miners
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-bnz&identifier=608

data1 = read.csv('./raw_data/608_ALM_eggSurvey_2004-2015.txt',as.is=T,check.names=F,header=T)
data1$Number = data1$Total_esm
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	st = paste(data2$Site,data2$Tree,sep='_')
	n.tree = length(unique(st))
	out[ox,1] = 'BonanzaCreek' #LTER
	out[ox,2] = 'leafminers'
	out[ox,3] = years[y] #Year
	out[ox,4] = sum(data2$Number,na.rm=T) / n.tree
	ox = ox + 1
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/BonanzaCreek/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'BonanzaCreek'
		trends[tx,2] = 'Leafminers'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/BonanzaCreek_Leafminers_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/BonanzaCreek_Leafminers_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Coweeta - aquatic invertebrates
# Link broken: http://coweeta.uga.edu/dbpublic/dataset_details.asp?accession=3023

data1 = read.csv('./raw_data/3023_invertebrate_1_6046.CSV',as.is=T,check.names=F,header=T)
data1 = data1[-c(1:2),c(2:4,grep('Abundance',colnames(data1)))]
data1 = data1[,-grep('Total',colnames(data1))]
data1$Year = data1$Starting_Year
species = gsub('_Abundance','',colnames(data1)[4:13])
spcol = 4:13
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	n = length(unique(data2$Site))
	for (s in 1:length(species)){
		out[ox,1] = 'Coweeta' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(as.numeric(data2[,spcol[s]])) / n
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Coweeta/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Coweeta'
		trends[tx,2] = 'StreamInvertebrates'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Coweeta_StreamInvertebrates_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Coweeta_StreamInvertebrates_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Georgia Coastal Ecosystems - Burrowing crabs
# http://gce-lter.marsci.uga.edu/public/app/dataset_details.asp?accession=INV-GCES-1609

data1 = read.table('./raw_data/INV-GCES-1609_1_0.TXT',sep='\t',as.is=T,check.names=F,header=T)
data1$Number = data1$Hole_Density
data1$Zone.Plot = paste(data1$Zone,data1$Plot,sep='_')
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	szp = paste(data1$Site,data1$Zone.Plot,sep='_')
	n = length(unique(szp))
	out[ox,1] = 'GeorgiaCoastal' #LTER
	out[ox,2] = 'BurrowingCrabs'
	out[ox,3] = years[y] #Year
	out[ox,4] = sum(data2$Number,na.rm=T) / n
	ox = ox + 1
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/GeorgiaCoastal/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'GeorgiaCoastal'
		trends[tx,2] = 'BurrowingCrabs'
		trends[tx,3] = 'BurrowingCrabs'
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/GeorgiaCoastal_BurrowingCrabs_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/GeorgiaCoastal_BurrowingCrabs_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Georgia Coastal Ecosystems - Fiddler crabs
# http://gce-lter.marsci.uga.edu/public/app/dataset_details.asp?accession=INV-GCED-1406

data1 = read.table('./raw_data/INV-GCED-1406_1_0.TXT',sep='\t',as.is=T,check.names=F,header=T)
data1$Number = (data1$Adult_Fiddler_crab_burrow_count / data1$Quadrat_Area_Adult_Fiddler_Crab_Burrow_Count) + (data1$Juvenile_Fiddler_crab_burrow_count / data1$Quadrat_Area_Juvenile_Fiddler_Crab_Burrow_Count)
data1$Site.Treatment = paste(data1$Site,data1$Experimental_Treatment,sep='_')
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	n = length(unique(data2$Site.Treatment))
	out[ox,1] = 'GeorgiaCoastal' #LTER
	out[ox,2] = 'FiddlerCrabs'
	out[ox,3] = years[y] #Year
	out[ox,4] = mean(data2$Number,na.rm=T) / n
	ox = ox + 1
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/GeorgiaCoastal/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'GeorgiaCoastal'
		trends[tx,2] = 'FiddlerCrabs'
		trends[tx,3] = 'FiddlerCrabs'
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/GeorgiaCoastal_FiddlerCrabs_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/GeorgiaCoastal_FiddlerCrabs_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Georgia Coastal Ecosystems - Grasshoppers
# https://gce-lter.marsci.uga.edu/public/app/data_catalog.asp

files = list.files('./raw_data/GeorgiaCoast_grasshoppers',full.names=T)
data1 = c()
for (f in 1:length(files)){ #merge yearly data into one dataset
	add.data = read.table(files[f],sep='\t',as.is=T,check.names=F,header=T)[,1:7]
	data1 = data.frame(rbind(data1,add.data))
}
data1$Species = data1$Species_code
data1$Number = data1$Count
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	st = paste(data2$Site,data2$Transect,sep='_')
	n = length(unique(st))
	species = unique(data2$Species)
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		out[ox,1] = 'GeorgiaCoastal' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/GeorgiaCoastal/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'GeorgiaCoastal'
		trends[tx,2] = 'Grasshoppers'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/GeorgiaCoastal_Grasshoppers_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/GeorgiaCoastal_Grasshoppers_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Georgia Coastal Ecosystems - Invertebrates
# http://gce-lter.marsci.uga.edu/public/app/dataset_details.asp?accession=PLT-GCEM-1610

data1 = read.table('./raw_data/PLT-GCEM-1610_Animals_3_0.TXT',sep='\t',as.is=T,check.names=F,header=T)
data1$Year = apply(array(data1$Date),1,function(x){strsplit(x,'-')[[1]][1]})
data1$Zone.Plot = paste(data1$Spartina_Zone,data1$Plot_Number,sep='_')
species = c('Acrididae','Prokelisia')
years = unique(data1$Year)
spcol = 6:7
sites = unique(data1$Zone.Plot)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	for (s in 1:length(species)){
		out[ox,1] = 'GeorgiaCoastal' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = mean(data2[,spcol[s]],na.rm=T)
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/GeorgiaCoastal/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'GeorgiaCoastal'
		trends[tx,2] = 'Invertebrates'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/GeorgiaCoastal_Invertebrates_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/GeorgiaCoastal_Invertebrates_abundance.txt'),sep='\t',quote=F,row.names=F)


#############################################################################################
# Midwest Aphid Suction Trap Network
# https://suctiontrapnetwork.org/

data1 = read.table('C:/Users/mcros/Desktop/Postdoc UGA/Aphid_STN/curated_data/STN_counts_curated_time-consistent.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$Number = apply(data1[,6:7],1,function(x){if(is.na(x[1]) & is.na(x[2])){NA}else{sum(x,na.rm=T)}})
years = unique(data1$Year)
out = data.frame('LTER.site'=NA,'Species.code'=NA,'Year'=-999,'Abundance'=-999,stringsAsFactors=F)
ox = 1
for (y in 1:length(years)){
	data2 = data1[which(data1$Year==years[y]),]
	sw = paste(data2$Site,data2$Date,sep='_')
	n.days = 7 * length(unique(sw))
	species = unique(data2$Species)
	for (s in 1:length(species)){
		data3 = data2[which(data2$Species==species[s]),]
		out[ox,1] = 'Midwest' #LTER
		out[ox,2] = species[s]
		out[ox,3] = years[y] #Year
		out[ox,4] = sum(data3$Number,na.rm=T) / n.days
		ox = ox + 1
	}
}

# Estimate abundance trends
species = sort(unique(out$Species.code))
trends = data.frame('LTER'=NA,'Dataset'=NA,'Species'=NA,'Abundance.trend'=-999)
tx = 1
for (s in 1:length(species)){
	dat2 = out[which(out$Species.code==species[s]),]
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
#					png(paste0('./plots/Midwest/',species[s],'.png'))
#					plot(t.scale,Z,main=species[s],xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
#					abline(a=arr.Z$coef[1],b=abundance.trend,lty=2,col='red',lwd=1.5)
#					dev.off()
		}
		trends[tx,1] = 'Midwest'
		trends[tx,2] = 'Aphids'
		trends[tx,3] = species[s]
		trends[tx,4] = abundance.trend
		tx = tx + 1
	} else { #time series did not meet filtering criteria
	}	
}
write.table(trends,paste0('./summary_tables/m4_standard/Midwest_Aphids_trends.txt'),sep='\t',quote=F,row.names=F)
write.table(out,paste0('./summary_tables/m4_standard/Midwest_Aphids_abundance.txt'),sep='\t',quote=F,row.names=F)


######################################################################################################
# Compile all trends into a single dataframe

files = list.files(paste0('./summary_tables/m4_standard'),full.names=T)
files = files[grep('_trends.txt',files)]
trends = c()
for (i in 1:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)
	if (ncol(add)>4){
		species1 = unique(add$Species)
		add2 = add[match(species1,add$Species),][,1:4]
		add2$Abundance.trend = apply(array(species1),1,function(x){mean(add$Abundance.trend[which(add$Species==x)],na.rm=T)})
	} else {
		add2 = add
	}
	trends = data.frame(rbind(trends,add2),stringsAsFactors=F)
}
write.table(trends,'./summary_tables/m4_standard/time_trends_arthropods_m4_wDatasetNames_standard.txt',sep='\t',row.names=F,quote=F)

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
trends$LTER.Dataset = paste(trends$LTER,trends$Dataset,sep=' ')
trends$Spe = paste(trends$LTER.Dataset,trends$Species,sep='_')
usp = unique(trends$Spe)
trends2 = trends[match(usp,trends$Spe),]
trends2$Abundance.trend = apply(array(usp),1,function(x){mean(trends$Abundance.trend[which(trends$Spe==x)])})
write.table(trends2,'./summary_tables/m4_standard/time_trends_arthropods_m4_standard.txt',sep='\t',row.names=F,quote=F)

# Compile abundances into a single file
files = list.files(paste0('./summary_tables/m4_standard'),full.names=T)
files = files[grep('_abundance.txt',files)]
abundances = read.table(files[1],sep='\t',as.is=T,check.names=F,header=T)
abundances$Method = NA
abundances$Dataset = strsplit(strsplit(files[1],'/')[[1]][5],'_')[[1]][2]
for (i in 2:length(files)){
	add = read.table(files[i],sep='\t',as.is=T,check.names=F,header=T)
	if (length(which(is.na(add$LTER)))>0){print(i)}else{}
	if (ncol(add)==4){
		add$Method = NA
	} else {
		colnames(add)[5] = 'Method'
	}
	add$Dataset = strsplit(strsplit(files[i],'/')[[1]][5],'_')[[1]][2]
	abundances = data.frame(rbind(abundances,add),stringsAsFactors=F)

}
write.table(abundances,'./summary_tables/m4_standard/PerSpecies_Abundance_LTER_m4_standard.txt',sep='\t',row.names=F,quote=F)


######################################################################################################
# Mean & 95% CI

dat = read.table(paste0('./summary_tables/m4_standard/time_trends_arthropods_m4_standard.txt'),sep='\t',as.is=T,check.names=F,header=T)
dat$LTER.Dataset = paste(dat$LTER,dat$Dataset,sep=' ')
levels1 = c('Arctic StreamInsects','Coweeta StreamInvertebrates','NorthTemperateLakes Crayfish','NorthTemperateLakes Macroinvertebrates','BonanzaCreek BarkBeetles','BonanzaCreek Leafminers','CedarCreek Grasshoppers',
	'CedarCreek Sweeps','GeorgiaCoastal Crabs','GeorgiaCoastal Grasshoppers','GeorgiaCoastal PlantHoppers','HarvardForest Ants','HarvardForest Ticks','HubbardBrook Leptidoptera','KonzaPrairie GallInsects','KonzaPrairie Grasshoppers',
	'Sevilleta Grasshoppers','Sevilleta Pitfall','Baltimore Mosquitoes','Midwest Aphids','Phoenix Pitfalls','Phoenix Sweep')
dat$LTER.Dataset = factor(dat$LTER.Dataset,levels=levels1)
dat$Habitat = 'terrestrial'
dat$Habitat[which(dat$LTER=='Arctic' | dat$LTER=='Coweeta' | dat$LTER=='NorthTemperateLakes')] = 'aquatic'
dat$Habitat[which(dat$LTER=='Midwest' | dat$LTER=='Baltimore' | dat$LTER=='Phoenix')] = 'impacted'
dat$Habitat = factor(dat$Habitat,levels=c('aquatic','terrestrial','impacted'))

dat$Spe = paste(dat$LTER,dat$Dataset,dat$Species,sep='_')
Spe = unique(dat$Spe)
L = unique(dat$LTER)
means1 = apply(array(L),1,function(x){mean(dat$Abundance.trend[which(dat$LTER==x)],na.rm=T)}) #mean of species*subplots

t1 = t.test(means1)
t1

mn1 = t1[['estimate']]
uppers = t1[['conf.int']][2]
lowers = t1[['conf.int']][1]
png(paste0('./plots/m4_standard/t-test_abundance_lters_barplot_Species.png'),width=200,height=500)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
bars = barplot(mn1,ylim=c(-1.5,1.5),names.arg='',ylab='Average change in abundance',main=NULL,lwd=2)
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
dp2 = ggplot(dat, aes(x=LTER.Dataset, y=Abundance.trend, fill=Habitat)) + 
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
	ggsave(filename='./plots/m4_standard/ALL_time_trends_boxplot_relaxed_noannotations_violin_Species.png',plot=dp2,dpi=600,unit='cm',height=8,width=24)
tcount = apply(array(levels(dat$LTER.Dataset)),1,function(x){length(which(as.character(dat$LTER.Dataset)==x))})
cbind(levels(dat$LTER.Dataset),tcount)

data1 = read.table('./summary_tables/m4_standard/PerSpecies_Abundance_LTER_m4.txt',sep='\t',as.is=T,check.names=F,header=T)
data1$LTER.Dataset = paste(data1$LTER,data1$Dataset,sep='_')
ycount = apply(array(unique(data1$LTER.Dataset)),1,function(x){sort(unique(data1$Year[which(data1$LTER.Dataset==x)]))})
names(ycount) = unique(data1$LTER.Dataset)


####################################################################################################
# Calculate diversity trends

setwd()
library(vegan)
source('./code/AR_reml.R')

# Vegan rarefy function that accepts non-integer values
rarefy1 = function (xx, sample1, se=FALSE, MARGIN=1, divfun){
    xx <- as.matrix(xx)
    minsample <- min(apply(xx, MARGIN, sum))
    rarefun <- function(y, sample2) {
        y <- y[y > 0]
        J <- sum(y)
        ldiv <- lchoose(J, sample2)
        p1 <- ifelse(J - y < sample2, 0, exp(lchoose(J - y, sample2) - ldiv)) #p1 is larger for rarer individuals
		out <- sum(1 - p1) #original output of rarefy function
		return(out)
	}
    S.rare <- apply(xx, MARGIN, rarefun, sample2 = sample1)
    attr(S.rare, "Subsample") <- sample1
    return(S.rare)
}

data1 = read.table('./summary_tables/m4_standard/PerSpecies_Abundance_LTER_m4_standard.txt',sep='\t',as.is=T,check.names=F,header=T)

# Remove dataset time series shorter than 4 years, and with fewer than 9 species
data1$LL = paste(data1$LTER,data1$Dataset,data1$Method,sep='_')
data1$LLY = paste(data1$LL,data1$Year,sep='_')
LL = unique(data1$LL) #39
LL.ycount = apply(array(LL),1,function(x){length(unique(data1$Year[which(data1$LL==x)]))})
LL2 = LL[which(LL.ycount>3)] #35
LL2.scount = apply(array(LL2),1,function(x){length(unique(data1$Species.code[which(data1$LL==x)]))})
LL3 = LL2[which(LL2.scount>=9)] #19
k = unlist(apply(array(LL3),1,function(x){which(data1$LL==x)}))
data2 = data1[k,]
data2$Abundance[which(data2$LTER.site=='Midwest')] = data2$Abundance[which(data2$LTER.site=='Midwest')] * 20 #transform aphids/day into aphids/20 days to improve rarefaction
data2$Abundance[which(data2$LTER.site=='Phoenix' & data2$Dataset=='Pitfall2')] = data2$Abundance[which(data2$LTER.site=='Phoenix' & data2$Dataset=='Pitfall2')] * 10 #transform arthropods/trap into arthropods/10 traps to improve rarefaction
LL4 = unique(data2$LL)
LL4 = LL4[-grep('Coweeta',LL4)] #exclude Coweeta datasets that report taxa by feeding type
LL4.counts = apply(array(LL4),1,function(x){ys=unique(data2$Year[which(data2$LL==x)]);apply(array(ys),1,function(x2){sum(data2$Abundance[which(data2$LL==x & data2$Year==x2)])})})
LL4s = apply(array(LL4),1,function(x){y=strsplit(x,'_')[[1]];return(paste(y[1],y[2],sep='_'))})
uLL4 = unique(LL4s)
out = data.frame('LTER.Dataset'=NA,'Minsamp'=-999,'Maxsamp'=-999,'P10'=-999,'P20'=-999,'length1'=-999)
for (i in 1:length(uLL4)){
	vals = c()
	ls1 = which(LL4s==uLL4[i])
	for (j in 1:length(ls1)){
		vals = c(vals,LL4.counts[[ls1[j]]])
	}
	out[i,1] = uLL4[i]
	out[i,2] = min(vals[which(vals>0)])
	out[i,3] = max(vals)
	out[i,4] = quantile(vals,probs=0.1)
	out[i,5] = quantile(vals,probs=0.2)
	out[i,6] = length(which(vals > out[i,4]))
}
out$minsamp = out$P10
out$minsamp[which(out$LTER.Dataset=='HarvardForest_Ants1')] = out$P20[which(out$LTER.Dataset=='HarvardForest_Ants1')] #use quantile=0.2 because quantile=0.1 gives a rarefaction sample = 0
out$minsamp[which(out$LTER.Dataset=='Baltimore_Mosquitoes')] = out$Minsamp[which(out$LTER.Dataset=='Baltimore_Mosquitoes')] #use minimum because quantile=0.1 removes 2 of the 5 years of data
LL4.minsamps = lapply(LL4s,function(x){ceiling(out$minsamp[which(out$LTER.Dataset==x)])})

# Diversity metrics per dataset per year
div = c()
for (i in 1:length(LL4)){
	print(noquote(paste0(i,' of ',length(LL4))))
	data3 = data2[which(data2$LL==LL4[i]),]
	years = sort(unique(data3$Year))
	species = unique(data3$Species.code)
	mat = matrix(0,nrow=length(years),ncol=length(species)) #community data matrix
	for (y in 1:length(years)){
		for (s in 1:length(species)){
			a = data3$Abundance[which(data3$Year==years[y] & data3$Species.code==species[s])]
			if (length(a)==0){
			} else {
				mat[y,s] = a
			}
		}
	}
	ky = which(rowSums(mat)>LL4.minsamps[[i]])
	if (length(ky)<4){ #skip datasets that do not have any years with abundance/density > 9
		print(noquote(LL4[i]))
	} else {
		mat2 = mat[ky,]
		rarefied.richness = rarefy1(mat2,LL4.minsamps[[i]])
		sdi = diversity(mat2,index='shannon')
		evenness = sdi / log(length(species)) #Pielou's evenness
		jaccards = as.matrix(1-betadiver(mat2,method="j")) #presence/absence-based, Jaccard (dis)similarity 
		jaccard = jaccards[row(jaccards)==col(jaccards)+1]
		div.add = data.frame('LTER'=rep(strsplit(LL4[i],'_')[[1]][1],length(ky)),'Dataset'=rep(LL4[i],length(ky)),'Year'=years[ky],'Rarefied.richness'=rarefied.richness,'Evenness'=evenness,'Jaccard'=c(NA,jaccard))
		div = data.frame(rbind(div,div.add),stringsAsFactors=F)
	}
}

# Diversity trends
datasets = unique(div$Dataset)
trends = data.frame('LTER'=NA,'Dataset'=NA,'Locale'=NA,'Rarefied.richness.trend'=-999,'Evenness.trend'=-999,'Jaccard.trend'=-999)
vars = c('Rarefied.richness','Evenness','Jaccard')
for (i in 1:length(datasets)){
	print(noquote(i))
	dat = div[which(div$Dataset==datasets[i]),]
	if (nrow(dat)<4){ #skip time series with < 4 data points
	} else {
		trends[i,1] = strsplit(datasets[i],'_')[[1]][1]
		trends[i,2] = datasets[i]
		trends[i,3] = strsplit(datasets[i],'_')[[1]][3]
		for (v in 1:length(vars)){
			vals = dat[,which(colnames(dat)==vars[v])]
			if (length(which(vals[!is.na(vals)]==0))>0){
				cmin = min(vals[which(!is.na(vals) & vals>0)])
				vals[which(vals==0)] = 0.5*cmin
			} else {
			}
			log.vals = log(vals[!is.na(vals)]) #remove NA from first year in beta diversity time series, and log-transform
			if (length(log.vals)<4){
				#skip instances where number of Jaccard values is < 4
			} else {
				ys = sort(dat$Year)
				ys.stretch = seq(ys[1],ys[length(ys)],1) #expand time series to include missing years
				X = match(ys.stretch,ys)
				X[which(!is.na(X))] = log.vals[X[which(!is.na(X))]] #fill-in abundance values in the expanded time series
				Z = (X – mean(X, na.rm=T))/sd(X, na.rm=T) # z-transform
				t.scale = 1:length(ys.stretch)
				t.scale = (t.scale-min(t.scale))/max(t.scale)
				arr.Z = AR_reml(Z ~ t.scale) #Z-transformed time trends
				trends[i,which(gsub('.trend','',colnames(trends))==vars[v])] = arr.Z$coef[2]
				png(paste0('./plots/m4_standard/diversity/',vars[v],'_',datasets[i],'.png'))
				plot(t.scale,Z,main=paste0(vars[v],' ',datasets[i]),xlab='Scaled time',ylab='Z-transformed abundance',lwd=2,cex=2,pch=16)
				abline(a=arr.Z$coef[1],b=arr.Z$coef[2],lty=2,col='red',lwd=1.5)
				dev.off()
			}
		}
	}
}
write.table(trends,'./summary_tables/m4_standard/arthropod_diversity_trends_m4_standard.txt',sep='\t',quote=F,row.names=F)


# Boxplots of trends
trends = read.table('./summary_tables/m4_standard/arthropod_diversity_trends_m4.txt',sep='\t',as.is=T,check.names=F,header=T)
trends$dataset = apply(array(trends$Dataset),1,function(x){y=strsplit(x,'_')[[1]];return(paste(y[1],y[2],sep='_'))})
trends$dataset[which(trends$dataset=='Phoenix_Pitfall1' | trends$dataset=='Phoenix_Pitfall2')] = 'Phoenix_Pitfalls'
trends$dataset[which(trends$dataset=='HarvardForest_Ants1' | trends$dataset=='HarvardForest_AntsNantucket')] = 'HarvardForest_Ants'
trends$dataset[which(trends$dataset=='CedarCreek_Sweep1' | trends$dataset=='CedarCreek_Sweep2')] = 'CedarCreek_Sweeps'
levels1 = c('Arctic_StreamInsects','NorthTemperateLakes_Benthic','CedarCreek_Grasshoppers','CedarCreek_Sweeps','HarvardForest_Ants','KonzaPrairie_Grasshoppers','Sevilleta_Grasshoppers','Sevilleta_Pitfall',
	'Baltimore_Mosquitoes','Midwest_Aphids','Phoenix_Pitfalls','Phoenix_Sweep')
trends$dataset = factor(trends$dataset,levels=levels1)
png('./plots/m4_standard/diversity_trends_boxplot_m4.png',res=300,width=480*2.5,height=480*4)
par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(1,5,1,1))
# Rarefied richness
plot(trends$dataset,trends$Rarefied.richness.trend,ylim=c(-5,5),xaxt='n',yaxt='n',axes=F,ylab='Rarefied richness trend (sd/yr)',col=c(rep('lightblue',2),rep('white',6),rep('darkorange',4)),cex.lab=2,cex.axis=1.5)
axis(2,lwd=2,at=seq(-5,5,2.5),labels=seq(-5,5,2.5),cex.axis=1.5)
abline(h=0,lty=1,col='black',lwd=2)
# Evenness
plot(trends$dataset,trends$Evenness.trend,ylim=c(-5,5),xaxt='n',yaxt='n',axes=F,ylab='Evenness trend (sd/yr)',col=c(rep('lightblue',2),rep('white',6),rep('darkorange',4)),cex.lab=2,cex.axis=1.5)
axis(2,lwd=2,at=seq(-5,5,2.5),labels=seq(-5,5,2.5),cex.axis=1.5)
abline(h=0,lty=1,col='black',lwd=2)
# Jaccard
plot(trends$dataset,trends$Jaccard.trend,ylim=c(-8,8),xaxt='n',yaxt='n',axes=F,ylab='Jaccard trend (sd/yr)',col=c(rep('lightblue',2),rep('white',6),rep('darkorange',4)),cex.lab=2,cex.axis=1.5)
axis(2,lwd=2,at=seq(-8,8,4),labels=seq(-8,8,4),cex.axis=1.5)
abline(h=0,lty=1,col='black',lwd=2)
dev.off()	

# Barplots
uL = unique(trends$dataset)
R.means = apply(array(uL),1,function(x){mean(trends$Rarefied.richness.trend[which(trends$dataset==x)],na.rm=T)})
R.t = t.test(as.numeric(R.means))
E.means = apply(array(uL),1,function(x){mean(trends$Evenness.trend[which(trends$dataset==x)],na.rm=T)})
E.t = t.test(as.numeric(E.means))
J.means = apply(array(uL),1,function(x){mean(trends$Jaccard.trend[which(trends$dataset==x)],na.rm=T)})
J.t = t.test(as.numeric(J.means))
png(paste0('./plots/m4_standard/t-test_diversity_richness.png'),width=200,height=450)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
mean1 = R.t[['estimate']]; uppers = R.t[['conf.int']][2]; lowers = R.t[['conf.int']][1]
bars = barplot(mean1,ylim=c(-4,4),names.arg='',ylab='Avg. richness trend',main=NULL,lwd=2); abline(h=0); axis(2,lwd=2) #text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2); arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()
png(paste0('./plots/m4_standard/t-test_diversity_evenness.png'),width=200,height=450)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
mean1 = E.t[['estimate']]; uppers = E.t[['conf.int']][2]; lowers = E.t[['conf.int']][1]
bars = barplot(mean1,ylim=c(-4,4),names.arg='',ylab='Avg. evenness trend',main=NULL,lwd=2); abline(h=0); axis(2,lwd=2) #text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2); arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()
png(paste0('./plots/m4_standard/t-test_diversity_Jaccard.png'),width=200,height=450)
par(cex.lab=2,cex.axis=1.5,oma=c(0,0,0,0),mar=c(15,8,2,2),lwd=2)
mean1 = J.t[['estimate']]; uppers = J.t[['conf.int']][2]; lowers = J.t[['conf.int']][1]
bars = barplot(mean1,ylim=c(-4,4),names.arg='',ylab='Avg. Jaccard trend',main=NULL,lwd=2); abline(h=0); axis(2,lwd=2) #text(x=bars,y=par("usr")[3]-0.05,srt=90,adj=1,labels='LTER',xpd=T,cex=2)
segments(bars,lowers,bars,uppers,lwd=2); arrows(bars,lowers,bars,uppers,lwd=2,angle=90,code=3,length=(par("usr")[2]-par("usr")[1])/4)
dev.off()

