library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library("ggsci")
source("toolbox_msk.R")
options(stringsAsFactors=FALSE)
  
changeNA2Zero<- function(dat) {	  
  for(i in c(1:dim(dat)[2])) {
	dat[which( is.na(dat[,i]) | is.na(as.character(dat[,i])) ), i]=0
  }
  return (dat)		
}
#----------------------Read in  & preprocess the data for plotting --------------------
kochdat=read.delim("input/Koch.Covid19.CH_min.txt", as.is=T, header=T)
clndat=unique(kochdat[,c(1:8)])

kochdat2 <- kochdat %>% mutate(
  ch_nondriver_az = case_when(
    ch_pd == 1 ~ 0, 
    ch_pd == 0 & CH_all == 1 ~ 1,
    CH_all == 0 ~ 0),
  ch_nondriver_3_az = case_when(
    ch_nondriver_az == 1 & TotVaf > 3  ~ 1, 
    T ~ 0
  ),
  ch_nondriver_4_az = case_when(
    ch_nondriver_az == 1 & TotVaf > 4  ~ 1, 
    T ~ 0
  ),
  ch_nondriver_5_az = case_when(
    ch_nondriver_az == 1 & TotVaf > 5  ~ 1, 
    T ~ 0
  ),
  ch_nondriver_10_az = case_when(
    ch_nondriver_az == 1 & TotVaf > 10  ~ 1, 
    T ~ 0
  ), 
  ch_nondriver_and_silent = case_when(
    ch_nondriver_az == 1 & CH_silent == 1 ~ 1, 
    T ~ 0
  )
)

D <- kochdat2 
D <- D %>% mutate(covid_cat=case_when(
  Severity==1 ~ "Severe COVID",
  Severity==0 ~ "Non-Severe COVID"
))

D <- D %>% mutate(covid_cat = factor(covid_cat, levels = c("Severe COVID","Non-Severe COVID")))

D <- D %>% filter(!is.na(covid_cat))


#----------------------Sup Fig 4 ----
pdf(paste("output/koCH_VariantFrequency-Severe.pdf",sep=""), width=10, height=16)
vafs=c(2)

for(i in c(1:length(vafs)) ) {

	minVaf=vafs[i]
	
	T=D %>% filter(TotVaf >= minVaf & ch_nondriver_az==1 & covid_cat=="Severe COVID")
	DF=data.frame(table(T$Gene))
	colnames(DF)=c("Gene","n")
	DF=DF[order(desc(DF$n), DF$Gene), ]
	
	T2=D %>% filter(TotVaf >= minVaf & ch_nondriver_and_silent==1 & covid_cat=="Severe COVID")
	DF2=data.frame(table(T2$Gene))
	colnames(DF2)=c("Gene","n")
	DF2=DF2[order(desc(DF2$n), DF2$Gene), ]	
	
	p1=ggplot(DF, aes(x=reorder(Gene, -n),y=n)) + geom_bar(stat="identity", position="dodge")+ xlab("Gene")+ ylab("n")+ scale_fill_nejm() + 
		theme_bw()+ theme(legend.position = "top")+ theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(),
						panel.grid.minor.y=element_blank(), panel.grid.minor.x=element_line(color="grey60", linetype="dashed")) +  ggtitle("Non-Driver gene frequency")
	p2=ggplot(DF2, aes(x=reorder(Gene, -n),y=n)) + geom_bar(stat="identity", position="dodge")+ xlab("Gene")+ ylab("n")+ scale_fill_nejm() +
		theme_bw()+ theme(legend.position = "top")+ theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(),
						panel.grid.minor.y=element_blank(), panel.grid.minor.x=element_line(color="grey60", linetype="dashed")) +  ggtitle("Silent gene frequency")


	grid.arrange(p1, p2, ncol=1, nrow=2)

}
dev.off()


#----------------------Sup Figure 6 ----

D_new <- kochdat2  
D_new <- D %>% mutate(covid_cat=case_when(
  Severity==1 & TotVaf>0 ~ "Severe COVID",
  Severity==0 ~ "Non-Severe COVID",
  Severity==0 ~ "Negative"
))

D_new <- D_new %>% mutate(covid_cat = factor(covid_cat, levels = c("Severe COVID","Non-Severe COVID","Negative")))

##Plot by age ##
minVaf=2
D_new <- D %>% filter(TotVaf >= minVaf)

changeAgelabels<-function(dat){
	agebreaks <- c(0,50,60,70,80,500)
	agelabels <- c("0-50","50-60","60-70","70-80","80+")
	
	dat$agegroups[dat$Age>=agebreaks[1] & dat$Age<agebreaks[2]]=agelabels[1]
	dat$agegroups[dat$Age>=agebreaks[2] & dat$Age<agebreaks[3]]=agelabels[2]
	dat$agegroups[dat$Age>=agebreaks[3] & dat$Age<agebreaks[4]]=agelabels[3]
	dat$agegroups[dat$Age>=agebreaks[4] & dat$Age<agebreaks[5]]=agelabels[4]
	dat$agegroups[dat$Age>=agebreaks[5]]=agelabels[5]
	
	return(dat)

}

D_new=changeAgelabels(D_new)

#plotAge(D_new)

	panel_theme = theme_bw() + theme(
		panel.border = element_blank(),
		legend.position = "top",
		panel.grid.minor = element_blank(),
		plot.subtitle = element_text(hjust = 0.5, size = 8),
		plot.title = element_text(face = 'bold', size = 12, hjust = 0, vjust = -11),
		panel.grid.major = element_blank(),
		strip.background = element_blank(),
		strip.text = element_text(size = 10),
		axis.text.y = element_text(size = 10),
		axis.text.x = element_text(size = 10),
		axis.title = element_text(size = 10),
		axis.line = element_line(),
		plot.margin = unit(c(0,0,0,0), 'pt')
	) 

	get_ch_grouped = function(D_new, CI = F) {

		CH_by_age_grouped = D_new %>% select(agegroups, CH) %>%
			mutate(CH = ifelse(is.na(CH), 0, CH)) %>%
			group_by(agegroups) %>%
			summarise(CH = sum(CH), total = n()) %>% 
			filter(!is.na(agegroups)) %>%
			mutate(freq = CH / total)
		
		if (CI) {
			CH_by_age_grouped = CH_by_age_grouped %>%
			cbind(
				apply(CH_by_age_grouped, 1, function(row) {
					CI = prop.test(row['CH'], row['total'], conf.level=0.95)$conf.int[1:2]
					return(c(lower = CI[1], upper = CI[2]))
				}) %>% t
			)
		}
		
		CH_by_age_grouped=data.frame(CH_by_age_grouped)
		
		return(CH_by_age_grouped)
	}

	font_size = 12

	age_curve_theme = 
	  theme(
		  legend.position = 'top',
		  legend.key.size = unit(5, 'mm'),
		  legend.title = element_blank(),
		  legend.direction = 'horizontal',
		  plot.title = element_text(hjust = -0.08),
		  axis.text.x = element_text(angle = 45, vjust = 0.5, size = font_size),
		  axis.text.y = element_text(size = font_size),
		  axis.title = element_text(size = font_size),
		  legend.text = element_text(size = font_size)
	  )

	## Age CH All
	aa1=get_ch_grouped(D_new %>% filter(covid_cat == "Severe COVID") %>% mutate(CH = CH_all)) %>% mutate(covid_cat = "Severe COVID")
	aa2=get_ch_grouped(D_new %>% filter(covid_cat == "Non-Severe COVID") %>% mutate(CH = CH_all)) %>% mutate(covid_cat = "Non-Severe COVID")

	agegroup_counts= rbind(aa1, aa2)

	agegroup_counts$covid_cat <- factor(agegroup_counts$covid_cat, levels = c("Non-Severe COVID", "Severe COVID"))

	agegroup_counts <- agegroup_counts %>% rowwise() %>% mutate(CI_lower = prop.test(CH, total, conf.level = .95)$conf.int[1], CI_upper = prop.test(CH, total, conf.level = .95)$conf.int[2]) %>% mutate(type = ifelse(grepl("Background", covid_cat), "simul", "actual"))

	pp1=ggplot(agegroup_counts, aes(x = agegroups, y = freq, group = covid_cat, color = covid_cat)) + 
	  geom_point(position = position_dodge(width = .4), aes(shape = type), show.legend = F) + 
	  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = .4), width = 0) + 
	  scale_color_manual(values = c("steelblue", "darkred", "lightgreen", "darkorange")) + 
	  ylab("Proportion of patients with CH mutations") + 
	  xlab("Age")
	  
	 pp1
	  
	ggsave("output/KoCH_age_vs_CH_with_covid_categories.pdf", width = 10, height = 4)


#--------------------- Sup Fig. 2 CH mutation number with severity------------
#--- CH by gene---#
D_wide <-D_new
D_long <-D_new #M

panel_theme = theme_bw() + theme(
    panel.border = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.title = element_text(face = 'bold', size = 12, hjust = 0, vjust = -11),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.line = element_line(),
    plot.margin = unit(c(0,0,0,0), 'pt')
) 

## Histogram by gene frequency
gene_list = D_long %>% count(Gene) %>% arrange(-n) %>% .$Gene %>% unique %>% .[1:10]

n_severe = D_wide %>% count(covid_cat) %>% filter(covid_cat == "Severe COVID") %>% pull(n)
n_non_severe = D_wide %>% count(covid_cat) %>% filter(covid_cat == "Non-Severe COVID") %>% pull(n)

# tally
D_new2 = D_long %>% filter(!is.na(covid_cat)) %>%
    reshape2::dcast(
        formula = Gene + covid_cat ~ .,
        value.var = 'Sample',
        fun.aggregate = function(MRNs) {length(unique(MRNs))}
    ) %>%
    dplyr::rename("n_patient" = ".") %>%
    mutate(
        prop_patient = case_when(
            covid_cat == "Severe COVID" ~ n_patient/n_severe,
            covid_cat == "Non-Severe COVID" ~ n_patient/n_non_severe
        )
    ) %>%
    filter(Gene %in% gene_list) %>%
    mutate(
        Gene = factor(Gene, gene_list)) %>%
    arrange(Gene)

D_new3 <- D_new2 %>% mutate(covid_cat = factor(covid_cat, levels = c("Severe COVID","Non-Severe COVID")))

  pp2=ggplot(
      D_new3,
      aes(x = Gene, y = prop_patient, fill = covid_cat)
  ) +
  geom_bar(stat = 'identity', position = "dodge", color = 'black', size = 0.25) +
  panel_theme +
  # theme(
  #     panel.grid.major = element_blank(), 
  #     panel.border = element_blank(),
  #     axis.line = element_line(colour = "black"),
  #     legend.title = element_blank(),
  #     legend.key.size = unit(5, 'mm'),
  #     legend.position = 'top',
  #     legend.direction = 'horizontal',
  #     axis.title = element_text(size = font_size),
  #     axis.text.x = element_text(angle = 45, hjust = 1, size = font_size),
  #     legend.text = element_text(size = font_size)
  # ) +
  ylab("Proportion with mutated Gene") +
  xlab('') +
  scale_fill_nejm()
  
  pp2

ggsave("output/KoCH_gene_vs_CH_with_covid_categories.pdf", width = 10, height = 4)




