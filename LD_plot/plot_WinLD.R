ld <- read.table("combine.ld",header = T)

ld$percentile <- ifelse(is.finite(ld$percentile),ld$percentile,0)
ggplot(data = ld,aes(w1,w2))+
  geom_tile(aes(fill=percentile))+
  scale_fill_gradientn(colours = c(brewer.pal(4,'Paired')[2],
                                   brewer.pal(6,'Set1')[6]),
                       breaks=c(0,0.5,1),limits=c(0,1))+
  theme(panel.grid = element_blank(),
        legend.position = '',
        panel.background = element_blank(),
        axis.line = element_line())+
  labs(x='Chr5 (Mbp)',y='Chr5 (Mbp)')
