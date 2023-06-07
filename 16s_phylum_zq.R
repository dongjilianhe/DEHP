# 1. load packages --------------------------------------------------------
library(readxl)
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(reshape2)
library(tidytidbits)
library(vegan)
library(ggsci)
library(ggpubr)
library(amplicon)
library(microbiome)
library(ggrepel)
library(ape)
library(phyloseq)
library(ggsignif)
library(patchwork)

# 2. 准备工作 -----------------------------------------------------------------
pro.dir <- "e:/五附院/学生/张强/new2/" #windows系统所用方式
setwd("e:/五附院/学生/张强/new2/phylum/")
# meta数据读入
metadata <- read_excel(paste0(pro.dir,"metadata_2.xlsx"))
# 读入phylum水平的asv数据
asv.phylum <- read.csv(paste0(pro.dir,"08_table_level_phylum.tsv"),header = T,sep="\t",quote="")
#shuzi-shuzhi
trans.to.numeric <- function(data){
  
  for(i in 2:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  return(data)
}
asv.phylum <- trans.to.numeric(data=asv.phylum)

# 将Taxonomy这一列分开来 ----
colnames(asv.phylum)[1] <- "asv"
asv.phylum <- separate(data=asv.phylum,col=asv,into=c('k','p'),sep=";")
asv.phylum.long <- gather(data=asv.phylum,key="sample",value="count",-c('k','p'))
asv.phylum.long$tax <- paste(asv.phylum.long$k,asv.phylum.long$p,sep=";")
asv.phylum.long <- asv.phylum.long[,c("sample",'k','p',"count","tax")]

#input asv.genus
#asv.genus2 <- asv.genus
#asv.genus2[1,]<- colnames(asv.genus2) 
#write.table(asv.genus2,file="asv.genus2.txt",sep="\t",quote=F,row.names = T,col.names = F)


{# 统计注释率 ----
  
  # 单个样本统计
  
  anno.rate <- data.frame()
  for(i in 1:59){
    sample.list <- unique(asv.phylum.long$sample) 
    asv.single.sample <- asv.phylum.long[asv.phylum.long$sample==sample.list[i],]
    asv.single.sample.nonzero <- asv.single.sample[asv.single.sample$count!=0,]
    anno.rate <- rbind(anno.rate,c(sample.list[i],sum(str_sub(asv.single.sample.nonzero$k,1,1)=="d")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$p,1,1)=="p")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$c,1,1)=="c")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$o,1,1)=="o")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$f,1,1)=="f")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$g,1,1)=="g")/nrow(asv.single.sample.nonzero)*100,sum(str_sub(asv.single.sample.nonzero$s,1,1)=="s")/nrow(asv.single.sample.nonzero)*100))
    
  }
  rm(sample.list,asv.single.sample.nonzero)
  colnames(anno.rate) <- c("sample","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
  anno.rate.long <- gather(anno.rate,Taxonomy,percent,-sample)
  anno.rate.long$percent <- as.numeric(anno.rate.long$percent)
  anno.rate.long$Taxonomy <- factor(anno.rate.long$Taxonomy,levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

  # 注释率箱线图
  colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(7)
  annot.rate.plot <- ggplot(anno.rate.long,aes(x=Taxonomy,y=percent,fill=Taxonomy))+
    geom_boxplot()+
    labs(x="Taxonomy",y="Ratio(%)")+
    scale_fill_manual(name="Taxonomy",values = colors)+
    theme(plot.title = element_text(hjust = 0.5,vjust=1,size=20),axis.title = element_text(size=15),axis.text  = element_text(size = 12))
  
  
  ggsave(annot.rate.plot, file="00.annot.rate.plot.pdf", width=8, height=8)
  
  
}  

{# phylum水平 ----
  # 数据整理
  phylum.asv <- asv.phylum.long[,c(1,3,4)]
  phylum.asv$p[str_sub(phylum.asv$p,1,1)!="p"] <- "p__unclassified"
  
  phylum.asv <- phylum.asv %>% group_by(sample,p) %>% summarise(count=sum(count))
  colnames(phylum.asv)[2] <- "Taxonomy"
  phylum.asv$Taxonomy <- gsub("p__","",phylum.asv$Taxonomy)
  
  
  # 为每一个phylum编号
  phylum.id <- data.frame(asv=unique(phylum.asv$Taxonomy))
  phylum.id$id <- paste0("phylum_",1:nrow(phylum.id))
  phylum.asv <- merge(phylum.asv,phylum.id,by.x="Taxonomy",by.y="asv")
  # rm(phylum.id)
  
  phylum.asv <- phylum.asv[,c(2,4,3)]
}

{# 将长数据形式的asv转换为宽数据 -----
  
  long.to.wide <- function(data){
    
    data <- spread(data,sample,count)
    data[is.na(data)] <- 0
    return(data)
    
  }
  
  phylum.table <- long.to.wide(data=phylum.asv)
 
  # 过滤掉测序量不正常的几个样本
  
  #remove.outlier <- function(data){
    # 	J3N、J33C、J37C、J3P、J31C
    #   J23C、J31N
    #outlier <- c("J3N","J33C","J37C","J3P","J31C","J23C","J31N","J5P")
    #data <- data[,-which(colnames(data) %in% outlier)]
  #}
  
  #genus.table2 <- remove.outlier(data=genus.table)
  #species.table2 <- remove.outlier(data=species.table)
  
  
}

{# 过滤低丰度 ----
  
  filter.low.count <- function(count){
    
    rownames(count) <- count[,1]
    count <- count[,-1]
    
    rm.list <- c()
    for (i in 1:nrow(count)){
      
      if(sum(count[i,]!=0)<=2 & max(count[i,])<10){
        
        rm.list <- c(rm.list,i)
        
      }
    }
    
    #count <- count[-rm.list,]
    return(count)
  }
  phylum.table.filt <- filter.low.count(count = phylum.table)
  #species.table.filt <- filter.low.count(count = species.table2)
  
}

{# 相对丰度计算 ----
  
  cal.rela.abun <- function(count,rowname.is.asv){
    
    if(rowname.is.asv){
      
      abun <- apply(count, 2, function(x){ x/sum(x) })
      
    }
    if(!rowname.is.asv){
      
      abun <- count
      abun[,2:ncol(abun)] <- apply(abun[,2:ncol(abun)], 2, function(x){ x/sum(x) })
      
    }
    return(abun)
    
  }
  
  phylum.abun <- cal.rela.abun(count=phylum.table.filt,rowname.is.asv=TRUE)
  #species.abun <- cal.rela.abun(count=species.table.filt,rowname.is.asv=TRUE)
  
}

# 3. 分组设置 -----------------------------------------------------------------
#
all.group <- data.frame(sample=colnames(phylum.table)[2:ncol(phylum.table)])
##tissue.group$group <- str_sub(tissue.group$sample,0,-3)
##tissue.group$group <- substring(tissue.group$sample,1,1)
#合并表格
all.group <- merge(all.group,metadata,by.x="sample",by.y="ID")
#将Group改为group
colnames(all.group) <- c("sample","group")  

young.group <- all.group[1:29,]
colnames(young.group) <- c("sample","group")
#young.group2 <- tissue.group[c(1:9,20:39),]#ACD
old.group <- all.group[30:49,]
colnames(old.group) <- c("sample","group")

##fenzu
##A.group <- data.frame(sample=tissue.group[tissue.group$group=="C","sample"])

df <- as.data.frame(phylum.abun)
df$C <-rowMeans(df[,str_sub(colnames(df),0,1)=="C"])
df$A <-rowMeans(df[,str_sub(colnames(df),0,1)=="A"])
df_CA <- df[,c(53,54)]
#做F/B，要单个样本丰度

##导出phylums水平的菌种丰度
#重新定义变量
data1=phylum.id
data2=df_CA
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
data=data[,c(4,2,3)]
data3 <- rbind(colnames(data),data)
#导出phylum水平的菌种丰度
write.table(data3,file="phylum.菌门丰度.txt",sep="\t",quote=F,row.names = T,col.names = F)

#导出phylums水平的df_CA
#write.table(df_CA,file="df_CA.txt",sep="\t",quote=F,row.names = T,col.names = F)
#导出phylums水平的phylum.id
#write.table(phylum.id,file="phylum.id.txt",sep="\t",quote=F,row.names = T,col.names = F)

#df1 <- as.data.frame(df)
#df1$C=df[,str_sub(colnames(df),0,-2)=="C"]
# 4. 后续分析 -----------------------------------------------------------------
{# 群落组成柱状图 ----
  
  draw.stack.bar <- function(abun.data,legend_name,id_level){
    
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
    
  
  
  phylum.stack.plot <- draw.stack.bar(abun.data = df_CA, id_level=phylum.id, legend_name = "phylum")
  ggsave(phylum.stack.plot, file="02.phylum.bottom.stackbar.mean.pdf", width=8, height=10)
  
  }
  
{# 数据抽平 ----
  
  data.t <- function(data,sample.list){
    
    if (sample.list=="all"){
      
      # rownames(data) <- data[,1]
      # data <- data[,-1]
      data <- as.data.frame(t(data))
      
    }else{
      
      # rownames(data) <- data[,1]
      # data <- data[,-1]
      data <- as.data.frame(t(data))
      
      data <- data[which(rownames(data) %in% sample.list$sample),]
    }
    
    
  }
  t.phylum.table.filt <- data.t(data=phylum.table.filt,sample.list = "all")
  #t.species.table.filt <- data.t(data=species.table.filt,sample.list = "all")
  
  # 抽平
  phylum.table.rarefy <- as.data.frame(rrarefy(t.phylum.table.filt,sample=min(rowSums(t.phylum.table.filt))))
  #species.table.rarefy <- as.data.frame(rrarefy(t.species.table.filt,sample=min(rowSums(t.species.table.filt))))
  
  # 检查抽平结果
  rowSums(phylum.table.rarefy)
  #rowSums(species.table.rarefy)
  
  
}

{# α多样性 ----
  
  {# α多样性计算 ----
    
    a.diversity <- function(count,sample.list,data.T){
      
      # 整理数据格式
      if(!data.T){
        
        rownames(count) <- count[,1]
        count <- count[,-1]
        count <- as.data.frame(t(count))
        count <- count[rownames(count) %in% sample.list$sample,]
        
      }
      if(data.T){
        
        count <- count[rownames(count) %in% sample.list$sample,]
        
      }
      observed_genus <- estimateR(count)[1,]
      shannon <- vegan::diversity(count, index = 'shannon', base = exp(1))
      Gini_simpson  <- vegan::diversity(count, index = 'simpson')
      chao1 <- estimateR(count)[2,]
      ace <- estimateR(count)[4,]
      goods_coverage <- 1 - rowSums(count == 1) / rowSums(count)
      diversity <- cbind(observed_genus,shannon,Gini_simpson,chao1,ace,goods_coverage)
      diversity <- merge(diversity,sample.list,by.x="row.names",by.y="sample")
      rownames(diversity) <- diversity[,1]
      diversity <- diversity[,-1]
      return(diversity)
    }
    
    #tissue2.group <- tissue.group
    #young2.group <- young.group
    #old2.group <- old.group
    
    #pdl2.group$sample <- gsub("C","N",pdl2.group$sample)
    
    
    phylum.a.diversity <- a.diversity(count = phylum.table.rarefy, sample.list = all.group, data.T = TRUE)
    phylum.a.diversity2 <- a.diversity(count = phylum.table.rarefy, sample.list = young.group, data.T = TRUE)
    phylum.a.diversity3 <- a.diversity(count = phylum.table.rarefy, sample.list = old.group, data.T = TRUE)
    #species.a.diversity <- a.diversity(count = species.table.rarefy, sample.list = pdl1.group, data.T = TRUE)
    
  }
  
  {# α多样性画图 ----
    
    form.comp <- function(data){
      
      pheno <-unique(data$group)
      result <- list()
      for (i in 1:(length(pheno)-1)){
        
        for(j in (i+1):length(pheno)){
          comp <- list(c(pheno[i],pheno[j]))
          result <- append(result,comp)
        }
        
      }
      return(result)
      
    }
    draw.a.diversity.plot <- function(data, index_name){
      
      mycomp <- form.comp(data=data)
      
      # colors <- brewer.pal(length(unique(data$group)), "Set3")
      # getPalette <- colorRampPalette(c("#ED000099","#42B54099","#00468B99"))
      # colors <- getPalette(length(unique(data$group)))
      colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(length(unique(data$group)))
      if(length(colors)==2){
        colors <- c(colors[2],colors[1])
        
      }
      
      colnames(data) <- c("Observed genus","Shannon","Simpson","Chao1","ACE","Goods coverage","group")
      index_col <- which(colnames(data)==index_name)
      
      symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
      
      plot <- ggplot(data,aes(x=group,y=data[,index_col],fill=group))+
        geom_boxplot()+
        scale_fill_manual(name="Group",values = colors)+
        labs(title=index_name, x="Group", y=index_name)+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))+
        # stat_compare_means(comparisons=mycomp,method = "wilcox")
        stat_compare_means(comparisons=mycomp,label = "p.signif",method = "wilcox",symnum.args=symnum.args)
      
      return(plot)
      
    }
    
    combine.a.plot <- function(data){
      
      obseved<- draw.a.diversity.plot(data=data,index_name = "Observed genus")
      shannon <- draw.a.diversity.plot(data=data,index_name = "Shannon")
      simpson <- draw.a.diversity.plot(data=data,index_name = "Simpson")
      chao1 <- draw.a.diversity.plot(data=data,index_name = "Chao1")
      ace <- draw.a.diversity.plot(data=data,index_name = "ACE")
      
      combined <- ggarrange(obseved,chao1,shannon,simpson,ace,ncol =3,nrow=2,common.legend = T,legend = "right")
      
      return(combined)
    }
    
    phylum.a.plot <- combine.a.plot(data=phylum.a.diversity)
    phylum.a.plot2 <- combine.a.plot(data=phylum.a.diversity2)
    phylum.a.plot3 <- combine.a.plot(data=phylum.a.diversity3)
    #species.a.plot <- combine.a.plot(data=species.a.diversity)
    
    ggsave(phylum.a.plot, file="03.phylum.all.a.diversity.pdf", width=9, height=11)
    ggsave(phylum.a.plot2, file="03.phylum.young.a.diversity.pdf", width=9, height=11)
    ggsave(phylum.a.plot3, file="03.phylum.old.a.diversity.pdf", width=9, height=11)
    #ggsave(species.a.plot, file="03.species.pdl1.a.diversity.pdf", width=9, height=11)
    
    
  }
  
  {# 稀释曲线 ----
    
    rare.plot <- function(data, group, type, rowname.is.null, row.is.asv){
      
      if(rowname.is.null){
        
        rownames(data) <- data[,1]
        data <- data[,-1]
        
      } 
      if(!row.is.asv){
        
        data <- as.data.frame(t(data))
        
      } 
      
      group.1 <- data.frame(Group=group$group)
      rownames(group.1) <- group$sample
      
      rare.plot <- alpha_rare_all(otu=data,map=group.1)
      rare.data <- rare.plot[[2]]
      
      colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(length(unique(group.1$Group)))
      if(length(colors)==2){
        colors <- c(colors[2],colors[1])
        
      }

      if(type=="sample"){
        
        p1 <- ggplot(data= rare.data,aes(x = i,y = index,group = ID,colour = Group)) +
          geom_smooth(span = 0.7, se = FALSE, method = "loess") +
          labs(x= "Number of sequences",y="Richness",title="Rarefaction curve")+
          theme_bw()+
          scale_color_manual(name="Group",values = colors)+
          theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
        return(p1)
        
      }
      if(type=="mean"){
        
        rare.data <- rare.data %>% group_by(Group,i) %>% summarise(mean(index), sd(index))
        colnames(rare.data) <- c(colnames(rare.data)[1:2],"mean","sd")
        rare.data <- na.omit(rare.data)
        
        p2 <- ggplot(data= rare.data,aes(x = i,y = mean,colour = Group)) +
          geom_smooth(span = 0.7,se = FALSE, method = "loess") +
          labs(x= "Number of sequences",y="Richness",title="Rarefaction curve(mean)")+
          theme_bw()+
          scale_color_manual(name="Group",values = colors)+
          theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
        return(p2)
        
      }
      if(type=="combined"){
        
        pheno <- unique(rare.data$Group)
        p1 <- ggplot(data= rare.data[rare.data$Group==pheno[1],],aes(x = i,y = index,group = ID,colour = Group)) +
          geom_smooth(span = 0.7, se = FALSE, method = "loess") +
          scale_color_manual(name="Group",values = "#ED0000B2")+
          labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[1],")"))+
          theme_bw()+
          theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
        
        p2 <- ggplot(data= rare.data[rare.data$Group==pheno[2],],aes(x = i,y = index,group = ID,colour = Group)) +
          geom_smooth(span = 0.7, se = FALSE, method = "loess") +
          scale_color_manual(name="Group",values = "#42B540B2")+
          labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[2],")"))+
          theme_bw()+
          theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))

        
        }
        
        rare.data <- rare.data %>% group_by(Group,i) %>% summarise(mean(index), sd(index))
        colnames(rare.data) <- c(colnames(rare.data)[1:2],"mean","sd")
        rare.data <- na.omit(rare.data)
        
         p3 <- ggplot(data= rare.data,aes(x = i,y = mean,colour = Group)) +
          geom_smooth(span = 0.7,se = FALSE, method = "loess") +
          scale_color_manual(name="Group",values = colors)+
          labs(x= "Number of sequences",y="Richness",title="Rarefaction curve(mean)")+
          theme_bw()+
          theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))

         p4 <- ggplot(data= rare.data[rare.data$Group==pheno[4],],aes(x = i,y = index,group = ID,colour = Group)) +
           geom_smooth(span = 0.7, se = FALSE, method = "loess") +
           scale_color_manual(name="Group",values = "#88419d")+
           labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[4],")"))+
           theme_bw()+
           theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
         
         p5 <- ggplot(data= rare.data[rare.data$Group==pheno[5],],aes(x = i,y = index,group = ID,colour = Group)) +
           geom_smooth(span = 0.7, se = FALSE, method = "loess") +
           scale_color_manual(name="Group",values = "#f768a1")+
           labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[5],")"))+
           theme_bw()+
           theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
         
         p6 <- ggplot(data= rare.data[rare.data$Group==pheno[6],],aes(x = i,y = index,group = ID,colour = Group)) +
           geom_smooth(span = 0.7, se = FALSE, method = "loess") +
           scale_color_manual(name="Group",values = "#bfd3e6")+
           labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[6],")"))+
           theme_bw()+
           theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
         
         #p7 <- ggplot(data= rare.data[rare.data$Group==pheno[7],],aes(x = i,y = index,group = ID,colour = Group)) +
         # geom_smooth(span = 0.7, se = FALSE, method = "loess") +
         # scale_color_manual(name="Group",values = "#4a1486")+
         # labs(x= "Number of sequences",y="Richness",title=paste0("Rarefaction curve(",pheno[7],")"))+
         # theme_bw()+
         # theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
         
    }
    
    rare.data <- rare.data %>% group_by(Group,i) %>% summarise(mean(index), sd(index))
    colnames(rare.data) <- c(colnames(rare.data)[1:2],"mean","sd")
    rare.data <- na.omit(rare.data)
    
    p7 <- ggplot(data= rare.data,aes(x = i,y = mean,colour = Group)) +
      geom_smooth(span = 0.7,se = FALSE, method = "loess") +
      scale_color_manual(name="Group",values = colors)+
      labs(x= "Number of sequences",y="Richness",title="Rarefaction curve(mean)")+
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5),text = element_text(size=15))
    
    #if(length(pheno)==7){
    
    #p9 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4,nrow=2)
    
    
    # return(p9)
    
    if(length(pheno)==6){
      
      p8 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,ncol=3,nrow=3)
      
      
      return(p8)
          
        
        
      
      
      
      
    }
    }
    
  phylum.rare.plot <- rare.plot(data=phylum.table.rarefy,group=all.group,type="combined",rowname.is.null=FALSE,row.is.asv = FALSE)
    #species.rare.plot <- rare.plot(data=species.table.rarefy,group=tissue.group,type="combined",rowname.is.null=FALSE,row.is.asv = FALSE)
    ggsave(phylum.rare.plot, file="04.phylum.rare.pdf", width=44, height=22)
    #ggsave(species.rare.plot, file="04.species.rare.pdf", width=22, height=11)
    
  }
  {# 物种累积曲线 ----
    
    spec.accum <- function(data,sample.list){
      
      pheno <- unique(sample.list$group)
      
      accumulate.plot <- specaccum(data[which(rownames(data) %in% sample.list$sample[which(sample.list$group==pheno[1])]),], method = 'random')
      accumulate.data <- data.frame(sample_number=accumulate.plot$sites,index=accumulate.plot$richness,sd=accumulate.plot$sd)
      p1 <- ggplot(data=accumulate.data,aes(x = sample_number,y = index,colour="#ED0000B2")) +
        geom_errorbar(data = accumulate.data,aes(ymin=index - sd, ymax=index + sd),alpha = 0.4, width=.2)+
        geom_line()+
        labs(x= "Number of samples",y="Richness",title=paste0("Accumulation curve(",pheno[1],")"))+
        theme_bw()+
        scale_color_manual(values = "#ED0000B2")+
        theme(legend.position = 'none',plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      accumulate.plot <- specaccum(data[which(rownames(data) %in% sample.list$sample[which(sample.list$group==pheno[2])]),], method = 'random')
      accumulate.data <- data.frame(sample_number=accumulate.plot$sites,index=accumulate.plot$richness,sd=accumulate.plot$sd)
      p2 <- ggplot(data=accumulate.data,aes(x = sample_number,y = index,colour="#42B540B2")) +
        geom_errorbar(data = accumulate.data,aes(ymin=index - sd, ymax=index + sd),alpha = 0.4, width=.2)+
        geom_line()+
        labs(x= "Number of samples",y="Richness",title=paste0("Accumulation curve(",pheno[2],")"))+
        theme_bw()+
        scale_color_manual(values = "#42B540B2")+
        theme(legend.position = 'none',plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      
      
      
      #accumulate.plot <- specaccum(data[which(rownames(data) %in% sample.list$sample[which(sample.list$group==pheno[7])]),], method = 'random')
      #accumulate.data <- data.frame(sample_number=accumulate.plot$sites,index=accumulate.plot$richness,sd=accumulate.plot$sd)
      #p7 <- ggplot(data=accumulate.data,aes(x = sample_number,y = index,colour="#4a1486")) +
      #  geom_errorbar(data = accumulate.data,aes(ymin=index - sd, ymax=index + sd),alpha = 0.4, width=.2)+
      #  geom_line()+
      #  labs(x= "Number of samples",y="Richness",title=paste0("Accumulation curve(",pheno[7],")"))+
      #  theme_bw()+
      #  scale_color_manual(values = "#4a1486")+
      #  theme(legend.position = 'none',plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      accumulate.plot <- specaccum(data[which(rownames(data) %in% sample.list$sample),], method = 'random')
      accumulate.data <- data.frame(sample_number=accumulate.plot$sites,index=accumulate.plot$richness,sd=accumulate.plot$sd)
      p3 <- ggplot(data=accumulate.data,aes(x = sample_number,y = index,colour="#ED0000FF")) +
        geom_errorbar(data = accumulate.data,aes(ymin=index - sd, ymax=index + sd),alpha = 0.4, width=.2)+
        geom_line()+
        labs(x= "Number of samples",y="Richness",title="Accumulation curve(all)")+
        theme_bw()+
        scale_color_manual(values = "black")+
        theme(legend.position = 'none',plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      plot <- ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = T,legend = "none")
      return(plot)
      
    }
    
    phylum.spec.plot <- spec.accum(data=phylum.table.rarefy,sample.list=all.group)
    #species.spec.plot <- spec.accum(data=species.table.rarefy,sample.list=tissue.group)
    
    ggsave(phylum.spec.plot, file="04.phylum.accumulation.pdf", width=44, height=22)
    #ggsave(species.spec.plot, file="04.species.accumulation.pdf", width=22, height=11)
    
    
  }
  

  
{# β多样性 ----

 {# 距离矩阵计算 ----
  
  # bray_curtis
  bray.phylum <- as.data.frame(as.matrix(vegdist(phylum.table.rarefy,method="bray")))
  #bray.species <- as.data.frame(as.matrix(vegdist(species.table.rarefy,method="bray")))
  
  # euclidean
  euclidean.phylum <- as.data.frame(as.matrix(vegdist(phylum.table.rarefy,method="euclidean")))
  #euclidean.species <- as.data.frame(as.matrix(vegdist(species.table.rarefy,method="euclidean")))
  
  # abund_jaccard
  jaccard.phylum <- as.data.frame(as.matrix(vegdist(phylum.table.rarefy,method="jaccard",binary = TRUE)))
  #jaccard.species <- as.data.frame(as.matrix(vegdist(species.table.rarefy,method="jaccard",binary = TRUE)))
  
  # unweighted_unifrac/weighted_unifrac
  unifrac.cal <- function(data){
    
    tree <- rtree(ncol(data), rooted = T, tip.label = colnames(data))
    unwei.unifrac <- as.data.frame(as.matrix(UniFrac(phyloseq(otu_table(data, taxa_are_rows = F), phy_tree(tree)), weighted=F)))
    wei.unifrac <- as.data.frame(as.matrix(UniFrac(phyloseq(otu_table(data, taxa_are_rows = F), phy_tree(tree)), weighted=T)))
    return(list(unwei=unwei.unifrac,wei=wei.unifrac))
    
  }
  unifrac.phylum <- unifrac.cal(data=phylum.table.rarefy)
  #unifrac.species <- unifrac.cal(data=species.table.rarefy)
  
 }

{# PCA ----
  
  pca.plot <- function( pca.data, group, labels){
    
    pca.result <- prcomp(pca.data,center = F,scale. = F)
    
    # 绘制主成分碎石图
    # fviz_screeplot(pca.result, addlabels = TRUE)
    
    # pca结果可视化
    pca_score <- as.data.frame(pca.result$x)
    pca_score$sample <- rownames(pca_score)
    
    pca_vis<- merge(pca_score,group,by="sample")
    rownames(pca_vis) <- pca_vis$sample
    pca_vis <- pca_vis[,-1]
    
    summ <- summary(pca.result)
    xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
    ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
    
    colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(length(unique(pca_vis$group)))
    if(length(colors)==3){
      colors <- c(colors[2],colors[3],colors[1])
    }
    
    if (labels==FALSE){
      
      p1 <- ggplot(data = pca_vis,aes(x=PC1,y=PC2,color=group))+
        stat_ellipse(level = 0.95)+
        geom_point()+
        labs(x=xlab,y=ylab,title="PCA")+
        # scale_fill_manual(groups = colors)+
        scale_colour_manual(values = colors)+
        guides(fill="none")+
        theme_classic()+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
      # geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")+
      # geom_vline(aes(xintercept=0),colour="#990000",linetype="dashed")
      
      
    }else{
      
      p1 <- ggplot(data = pca_vis,aes(x=PC1,y=PC2,color=group))+
        stat_ellipse(level = 0.95)+
        geom_point()+
        labs(x=xlab,y=ylab,title="PCA")+
        # scale_fill_manual(groups = colors)+
        scale_colour_manual(values = colors)+
        guides(fill="none")+
        theme_classic()+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))+
        # geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")+
        # geom_vline(aes(xintercept=0),colour="#990000",linetype="dashed")+
        geom_text_repel(aes(label = rownames(pca_vis)), size = 3, show.legend = FALSE, box.padding = unit(0.5, 'lines'))
      
      
    }
    
    return(p1)
    
  }
  
  pca1 <- pca.plot(pca.data=phylum.table.rarefy,group=all.group,labels=T)
  pca2 <- pca.plot(pca.data=phylum.table.rarefy,group=young.group,labels=T)
  pca3 <- pca.plot(pca.data=phylum.table.rarefy,group=old.group,labels=T)
  #pca2 <- pca.plot(pca.data=species.table.rarefy,group=pdl1.group,labels=T)
  
  ggsave(pca1, file="04.phylum.all.pca.pdf", width=11, height=8)
  ggsave(pca2, file="04.phylum.young.pca.pdf", width=11, height=8)
  ggsave(pca3, file="04.phylum.old.pca.pdf", width=11, height=8)
  #ggsave(pca2, file="04.species.pdl1.pca.pdf", width=11, height=8)
  
}

{# PCoA ----
  
  sign.comp <- function(group){
    
    pheno <-unique(group$group)
    result <- list()
    for (i in 1:(length(pheno)-1)){
      
      for(j in (i+1):length(pheno)){
        comp <- list(c(pheno[i],pheno[j]))
        result <- append(result,comp)
      }
      
    }
    return(result)
    
  }
  beta.pcoa <- function(dis_mat, group, groupID="group", ellipse=T, label=F, PCo=12) {
    
    # PCoA
    pcoa=cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
    points=as.data.frame(pcoa$points) # get coordinate string, format to dataframme
    eig=pcoa$eig
    points <- merge(points,group[,c("sample",groupID)],by.x="row.names",by.y="sample")
    rownames(points) <- points[,1]
    points <- points[,-1]
    colnames(points)=c("x", "y", "z","group")
    
    comp <- sign.comp(group=group)
    
    colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(length(unique(points$group)))
    if(length(colors)==3){
      colors <- c(colors[2],colors[3],colors[1])
      
    }
    
    # 按1、2轴绘图
    if (PCo == 12){
      p=ggplot(points, aes(x=x, y=y, color=group, shape=group))  +
        labs(x="",y="")+theme_bw()+guides(fill="none")
      # labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
      #      y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
    }
    # 按1、3轴绘图
    if (PCo == 13){
      p=ggplot(points, aes(x=x, y=z, color=group, shape=group)) +
        labs(x="",y="")+theme_bw()
      # labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
      #      y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""), color=groupID)
    }
    # 按2、3轴绘图
    if (PCo == 23){
      p=ggplot(points, aes(x=y, y=z, color=group, shape=group)) +
        labs(x="",y="")+theme_bw()
      # labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
      #      y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""), color=groupID)
    }
    
    p=p + geom_point(alpha=.7, size=2) + 
      labs(title="PCoA")+
      theme_bw()+
      geom_vline(xintercept = 0,lty="dashed")+
      geom_hline(yintercept = 0,lty="dashed")+
      theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))+
      scale_color_manual(name="Group",values = colors)+
      guides(shape=F)
    # stat_ellipse(data=points,geom = "polygon",aes(fill=group),alpha=0.3)
    
    # geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")+
    # geom_vline(aes(xintercept=0),colour="#990000",linetype="dashed")
    
    p.x <- ggplot(data=points,aes(x=group,y=x,fill=group))+
      geom_boxplot()+
      labs(y=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),x="")+
      scale_fill_manual(name="Group",values = colors)+
      theme_bw()+ coord_flip()+
      theme(text = element_text(size=15),axis.ticks.y = element_blank())+
      guides(fill="none")+scale_x_discrete(labels=NULL)+
      geom_signif(comparisons = comp,
                  step_increase = 0.1,
                  map_signif_level = T,
                  test = wilcox.test)
    
    p.y <- ggplot(data=points,aes(x=group,y=y,fill=group))+
      geom_boxplot()+
      labs(y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),x="")+
      scale_fill_manual(name="Group",values = colors)+
      theme_bw()+
      theme(text = element_text(size=15))+
      guides(fill="none")+scale_x_discrete(labels=NULL)+
      geom_signif(comparisons = comp,
                  step_increase = 0.1,
                  map_signif_level = T,
                  test = wilcox.test)
    p.z <- ggplot(data=points,aes(x=group,y=z,fill=group))+
      geom_boxplot()+
      labs(y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""),x="")+
      scale_fill_manual(name="Group",values = colors)+
      theme_bw()+
      theme(text = element_text(size=15),axis.ticks.x = element_blank())+
      guides(fill="none")+scale_x_discrete(labels=NULL)+
      geom_signif(comparisons = comp,
                  step_increase = 0.1,
                  map_signif_level = T,
                  test = wilcox.test)
    
    
    
    # 是否添加置信椭圆
    if (ellipse == T){
      p=p + stat_ellipse(level=0.95)
    }
    # 是否显示样本标签
    if (label == T){
      p=p + geom_text_repel(label=paste(rownames(points)), size=3)
    }
    
    if (PCo == 12){
      p.y <- p.y+theme(axis.ticks.x = element_blank())
      p.all<-(p.y+p+plot_layout(nrow = 1,widths = c(1,5)))-(plot_spacer()+p.x+plot_layout(nrow = 1,widths = c(1,5)))+plot_layout(ncol = 1,heights = c(5,1))
      
      
    }
    if (PCo == 13){
      p.all<-(p.z+p+plot_layout(nrow = 1,widths = c(1,5)))-(plot_spacer()+p.x+plot_layout(nrow = 1,widths = c(1,5)))+plot_layout(ncol = 1,heights = c(5,1))
      
      
    }
    if (PCo == 23){
      p.y <- p.y+coord_flip()+theme(axis.ticks.y = element_blank())
      p.all<-(p.z+p+plot_layout(nrow = 1,widths = c(1,5)))-(plot_spacer()+p.y+plot_layout(nrow = 1,widths = c(1,5)))+plot_layout(ncol = 1,heights = c(5,1))
      
      
    }
    
    return(p.all)
  }
  
  pcoa1 <- beta.pcoa(dis_mat = euclidean.phylum, group = all.group, groupID="group", ellipse=T, label=T, PCo=12)
  #pcoa2 <- beta.pcoa(dis_mat = euclidean.genus, group = young.group, groupID="group", ellipse=T, label=T, PCo=12)
  #pcoa3 <- beta.pcoa(dis_mat = euclidean.genus, group = old.group, groupID="group", ellipse=T, label=T, PCo=12)
  #pcoa2 <- beta.pcoa(dis_mat = euclidean.species, group = pdl1.group, groupID="group", ellipse=T, label=T, PCo=12)
  
  ggsave(pcoa1, file="04.phylum.all.euclidean.pcoa.pdf", width=11, height=8)
  #ggsave(pcoa2, file="04.genus.young.euclidean.pcoa.pdf", width=11, height=8)
  #ggsave(pcoa3, file="04.genus.old.euclidean.pcoa.pdf", width=11, height=8)
  #ggsave(pcoa2, file="plot/04.species.pdl1.euclidean.pcoa.pdf", width=11, height=8)
  
}
 
{# NMDS ----
  
  nmds.plot <- function(input, group, dist.name, dim=12, label=F){
    nmds <- metaMDS(input, k=3)
    nmds_point <- data.frame(nmds$points)
    
    nmds_point <- merge(nmds_point,group,by.x="row.names",by.y="sample")
    rownames(nmds_point) <- nmds_point[,1]
    nmds_point <- nmds_point[,-1]
    
    colors <- pal_lancet(palette = c("lanonc"), alpha = 0.7)(length(unique(nmds_point$group)))
    if(length(colors)==2){
      colors <- c(colors[2],colors[1])
    }
    
    
    if(dim==12){
      
      result <- ggplot(nmds_point, aes(x=MDS1, y=MDS2, color = group)) +
        geom_point() + 
        stat_ellipse(level = 0.95, show.legend = F, aes(fill = group)) +
        theme_classic() +  # 经典，只有xy轴
        labs(title = paste(dist.name,'NMDS  Stress =', round(nmds$stress, 4))) +
        scale_shape_manual(values=c(1, 2, 0)) +
        scale_color_manual(name="Group",values = colors)+
        guides(color = guide_legend(order=1), shape = guide_legend(order=2))+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      
    }
    if(dim==23){
      
      result <- ggplot(nmds_point, aes(x=MDS2, y=MDS3, color = group)) +
        geom_point() + 
        stat_ellipse(level = 0.95, show.legend = F, aes(fill = group)) +
        theme_classic() +  # 经典，只有xy轴
        labs(title = paste('Bray-Curtis  NMDS  Stress =', round(nmds$stress, 4))) +
        scale_shape_manual(values=c(1, 2, 0)) +
        scale_color_manual(name="Group",values = colors)+
        guides(color = guide_legend(order=1), shape = guide_legend(order=2))+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      
    }
    if(dim==13){
      
      result <- ggplot(nmds_point, aes(x=MDS1, y=MDS3, color = group)) +
        geom_point() + 
        stat_ellipse(level = 0.95, show.legend = F, aes(fill = group)) +
        theme_classic() +  # 经典，只有xy轴
        labs(title = paste('Bray-Curtis  NMDS  Stress =', round(nmds$stress, 4))) +
        scale_shape_manual(values=c(1, 2, 0)) +
        scale_color_manual(name="Group",values = colors)+
        guides(color = guide_legend(order=1), shape = guide_legend(order=2))+
        theme(plot.title=element_text(hjust=0.5),text = element_text(size=15))
      
      
    }
    if(label==T){
      
      result <- result+geom_text_repel(label=paste(rownames(nmds_point)), size=3)
      
    }
    
    
    return(result)
  }
  
  nmds1 <- nmds.plot(input=unifrac.phylum[[2]],dist.name="",group=all.group,dim=12,label=T)
  #nmds2 <- nmds.plot(input=unifrac.genus[[2]],dist.name="",group=young.group,dim=12,label=T)
  #nmds3 <- nmds.plot(input=unifrac.genus[[2]],dist.name="",group=old.group,dim=12,label=T)
  #nmds2 <- nmds.plot(input=unifrac.species[[2]],dist.name="",group=pdl1.group,dim=12,label=T)
  
  ggsave(nmds1, file="04.phylum.all.wei.unifrac.nmds.pdf", width=11, height=8)
  #ggsave(nmds2, file="04.genus.young.wei.unifrac.nmds.pdf", width=11, height=8)
  #ggsave(nmds3, file="04.genus.old.wei.unifrac.nmds.pdf", width=11, height=8)
  #ggsave(nmds2, file="04.species.pdl1.wei.jaccard.nmds.pdf", width=11, height=8)
  
  
}
}

{# 组间差异分析 ----
  
  {# anosim ----
    
    anosim.cal <- function(dist.mat, group){
      
      temp <- all.group
      rownames(temp) <- temp[,1]
      temp <- temp[rownames(temp) %in% colnames(dist.mat),]
      temp <- temp[,-1]
      result <- summary(anosim(dist.mat, temp, permutations=999))
      return(result)
      
    }
    
    anosim.cal(dist.mat = jaccard.phylum, group=all.group)
    #anosim.cal(dist.mat = jaccard.species, group=tissue.group)
    
    #plot
    
    #polt(anosim,col=c("#ED0000B2","#42B540B2","#00468BB2","#88419d","#f768a1","#bfd3e6","#4a1486","#ED0000FF"))
    
  }
  
  {# adonis ----
    
    adonis.cal <- function(dist.mat, group){
      
      result <- adonis(dist.mat~group,data=group,permutations = 999)
      return(result)
    }
    
    adonis.cal(dist.mat=unifrac.phylum[[1]],group=all.group)
    
  }
} 
 
# 5. 导出数据 -----------------------------------------------------------------
format.lefse <- function(data,group,level.id,row.is.asv,rowname.is.null,comp.group){
  
  data <- as.data.frame(data)
  
  if(row.is.asv){
    
    if(rowname.is.null){
      
      data <- merge(data,level.id,by="id")
      data <- data[,c(1,ncol(data),2:(ncol(data)-1))]
      # data$asv <- gsub("\\[","",data$asv)
      # data$asv <- gsub("\\]","",data$asv)
      
      rownames(data) <- data[,2]
      data <- data[,-c(1,2)]
      data <- data[,colnames(data) %in% group$sample]
      
      data <- as.data.frame(t(data))
      
      data$sample <- rownames(data)
      data <- merge(data,group,by="sample")
      data <- data[,c(1,ncol(data),2:(ncol(data)-1))]
      
      data <- as.data.frame(t(data))
      
    }else{
      
      data$id <- rownames(data)
      data <- merge(data,level.id,by="id")
      data <- data[,c(1,ncol(data),2:(ncol(data)-1))]
      # data$asv <- gsub("\\[","",data$asv)
      # data$asv <- gsub("\\]","",data$asv)
      
      rownames(data) <- data[,2]
      data <- data[,-c(1,2)]
      data <- data[,colnames(data) %in% group$sample]
      
      data <- as.data.frame(t(data))
      
      data$sample <- rownames(data)
      data <- merge(data,group,by="sample")
      data <- data[,c(1,ncol(data),2:(ncol(data)-1))]
      
      data <- as.data.frame(t(data))
      
      
    }
  }
  
  data <- data[,which(data["group",]==comp.group[1] | data["group",]==comp.group[2])]
  return(data)
  
}
phylum.output1 <- format.lefse(data=phylum.abun,group=all.group,comp.group=c("C","A"),row.is.asv=T,rowname.is.null=F,level.id=phylum.id)


#species.output <- format.lefse(data=species.abun,group=tissue.group,comp.group=c("N","P"),row.is.asv=T,rowname.is.null=F,level.id=species.id)

write.table(phylum.output1,file="phylum_C_A.txt",sep="\t",quote=F,row.names = T,col.names = F)


#write.table(species.output,file="./txt1/species_N_P.txt",sep="\t",quote=F,row.names = T,col.names = F)

#genus.output16 <- format.lefse(data=genus.abun,group=tissue.group,comp.group=c("YC","YL","YH","OC","OL","OH"),row.is.asv=T,rowname.is.null=F,level.id=genus.id)
#all.pdl1.species.output <- format.lefse(data=species.abun,group=pdl1.group,comp.group=c("positive","negative"),row.is.asv=T,rowname.is.null=F,level.id=species.id)

#write.table(genus.output16,file="all_genus_tissue.txt",sep="\t",quote=F,row.names = T,col.names = F)
#write.table(all.pdl1.species.output,file="./txt1/all_species_pdl1.txt",sep="\t",quote=F,row.names = T,col.names = F)



 