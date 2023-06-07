# 群落组成柱状图 ----
#门水平
df <- as.data.frame(phylum.abun)
#属水平
df <- as.data.frame(genus.abun)
#
df$YC <-rowMeans(df[,str_sub(colnames(df),0,-3)=="A"])
df$YL <-rowMeans(df[,str_sub(colnames(df),0,-4)=="C"])
df$YH <-rowMeans(df[,str_sub(colnames(df),0,-4)=="D"])
df$OC <-rowMeans(df[,str_sub(colnames(df),0,-4)=="E"])
df$OL <-rowMeans(df[,str_sub(colnames(df),0,-4)=="F"])
df$OH <-rowMeans(df[,str_sub(colnames(df),0,-4)=="G"])
df_young <- df[,c(60,61,62)]
df_old <- df[,c(63,64,65)]
df_all <- df[,c(60,61,62,63,64,65)]
#做F/B，要单个样本丰度
#df_single=as.data.frame(phylum.abun)

##导出phylum水平的#分组#菌种丰度
#重新定义变量
data1=phylum.id
data2=df_all

#data2=df_single

#把行名变为第一列，列名为id
data2=data2 %>% rownames_to_column("id")
# 两表ID名取交集,取交集的列名必须一致
ids <- Reduce(intersect,list(data1$id,data2$id))
# 表1去除在表2中没有的ID对应的行
data1 <- data1[data1$id %in% ids,]
# 根据某一列的重复值，去除该重复值所在的整行（保留第一行）
data1 <- data1[!duplicated(data1$id),]
# 左合并,data1到data2中
data <- dplyr::left_join(data2,data1,by = "id")
data=data[,c(8,2,3,4,5,6,7)]
data3 <- rbind(colnames(data),data)
#导出phylum水平的菌种丰度
write.table(data3,file="phylum.菌门丰度.txt",sep="\t",quote=F,row.names = T,col.names = F)
##导出genus水平的#分组#菌种丰度
#重新定义变量
data1=genus.id
data2=df_all
#把行名变为第一列，列名为id
data2=data2 %>% rownames_to_column("id")
# 两表ID名取交集,取交集的列名必须一致
ids <- Reduce(intersect,list(data1$id,data2$id))
# 表1去除在表2中没有的ID对应的行
data1 <- data1[data1$id %in% ids,]
# 根据某一列的重复值，去除该重复值所在的整行（保留第一行）
data1 <- data1[!duplicated(data1$id),]
# 左合并,data1到data2中
data <- dplyr::left_join(data2,data1,by = "id")
data=data[,c(8,2,3,4,5,6,7)]
data3 <- rbind(colnames(data),data)
#导出genus水平的菌种丰度
write.table(data3,file="genus.菌属丰度.txt",sep="\t",quote=F,row.names = T,col.names = F)
draw.stack.bar2 <- function(abun.data,legend_name,id_level){
  
  # 根据输入的组别数据挑选对应的样本的数据
  abun.data.filter <- as.data.frame(abun.data)
  #abun.data.filter <- abun.data.filter[,colnames(abun.data.filter) %in% facet_group$sample]
  abun.data.filter$id <- rownames(abun.data.filter)
  abun.data.filter <- abun.data.filter[,c(ncol(abun.data.filter),1:(ncol(abun.data.filter)-1))]
  
  # 根据某种要求选取要画图的物种（如丰度前10个，其余物种都归类为others）
  bar.data<- reshape2::melt(abun.data.filter,id.vars = "id",variable.name="sample",value.name = "abundance")
  
  # 获取每一个样本中丰度前10的物种
  get.topn<-function(bar.data,n){
    sample <- unique(bar.data$sample)
    result <- data.frame()
    for (i in 1:length(sample)){
      simple.data <- bar.data[which(bar.data$sample==sample[i]),]
      simple.data <- simple.data %>% arrange(desc(abundance)) %>% lump_rows( id, abundance, n = n, other_level = "Others" ) 
      result <- rbind(result,simple.data)
    }
    return(result)
  }
  
  bar.data <- get.topn(bar.data,10)
  
  # 获取颜色
  getPalette <- colorRampPalette(c("#D9D9D9","#FB8072","#FFFFB3","#FDB462","#B3DE69","#8DD3C7","#80B1D3","#BEBADA","#BC80BD"))
  colors <- getPalette(length(unique(bar.data$id)))
  
  # 添加分面分组信息
  #bar.data <- merge(bar.data,facet_group,by="sample")
  bar.data <- merge(bar.data,id_level,by="id",all.x=T)
  bar.data$asv[is.na(bar.data$asv)] <- "Others"
  bar.data$asv[bar.data$asv=="unclassified"] <- "Others"
  
  bar.data$asv <- fct_relevel( fct_reorder( bar.data$asv, bar.data$abundance, .desc = F),"Others" )
  
  # 画图
  {
    stack_bar <- ggplot(bar.data,aes(sample,abundance,fill=asv))+
    #width=0.8是调节柱子宽度
    geom_bar(stat="identity",width=0.8) +
    #facet_wrap(group.x~group.y,scales = "free_x")+
    #ncol = 1,注释的列数
    guides(fill = guide_legend( ncol = 1))+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(name=legend_name,values = colors)+
    labs(x= "Group",y="Relative abundance(%)",title=legend_name)+
    # theme(legend.position="bottom",plot.title = element_text(hjust = 0.5,vjust=1,size=20),panel.background = element_rect(fill = 'white', colour = 'white'),axis.title = element_text(size=15))
    #legend.position="bottom"注释在下面，legend.position="right"注释在右边
    theme(text = element_text(size=30),legend.position="right",plot.title = element_text(hjust = 0.5,vjust=1,size=30),axis.text = element_text(size=30))
  
  return(stack_bar)
    }
  
  
  
}

phylum.stack.plot <- draw.stack.bar2(abun.data = df_all,id_level=phylum.id, legend_name = "phylum")
phylum.stack.plot <- draw.stack.bar2(abun.data = df_young,id_level=phylum.id, legend_name = "phylum")
phylum.stack.plot <- draw.stack.bar2(abun.data = df_old,id_level=phylum.id, legend_name = "phylum")
#species.stack.plot <- draw.stack.bar(abun.data = species.abun, facet_group = site.group1, id_level=species.id, legend_name = "Species")

ggsave(phylum.stack.plot, file="02.phylum.tissue_all_bottom.stackbar.mean.pdf", width=11, height=10)
ggsave(phylum.stack.plot, file="022.phylum.tissue_young.stackbar.mean.pdf", width=8, height=10)
ggsave(phylum.stack.plot, file="022.phylum.tissue_old2.stackbar.mean.pdf", width=8, height=10)
#ggsave(species.stack.plot, file="02.species.tissue.stackbar.pdf", width=25, height=20)

genus.stack.plot <- draw.stack.bar2(abun.data = df_all,id_level=genus.id, legend_name = "Genus")
genus.stack.plot <- draw.stack.bar2(abun.data = df_young,id_level=genus.id, legend_name = "Genus")
genus.stack.plot <- draw.stack.bar2(abun.data = df_old,id_level=genus.id, legend_name = "Genus")
#species.stack.plot <- draw.stack.bar(abun.data = species.abun, facet_group = site.group1, id_level=species.id, legend_name = "Species")

ggsave(genus.stack.plot, file="02.genus.tissue_all_bottom.stackbar.mean.pdf", width=11, height=10)
ggsave(genus.stack.plot, file="022.genus.tissue_young.stackbar.mean.pdf", width=8, height=10)
ggsave(genus.stack.plot, file="022.genus.tissue_old.stackbar.mean.pdf", width=9, height=10)
#ggsave(species.stack.plot, file="02.species.tissue.stackbar.pdf", width=25, height=20)

