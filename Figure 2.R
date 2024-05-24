## F2A ----
load("Figure 2A1.Rdata")
#### my_pair.table ####  将一个列表中的多个数据框按照相同的列名和顺序进行对齐，并用0填充缺失的值。

#' @param feat_list 整合后的feat

my_pair.table=function(feat_list){
  library(tidyverse)
  
  
  feat_name=lapply(feat_list,function(data){colnames(data)})
  feat_ID=feat_name%>%unlist()%>%as.vector()%>%unique()
  b=lapply(feat_list, function(data){
    add=setdiff(feat_ID,colnames(data))
  })
  
  my_add=function(b,feat){
    if(length(b)!=0){
      b=as.character(b)
      x_add=data.frame(matrix(0,dim(feat)[1],length(b)))
      colnames(x_add)=b
      xtest=cbind(feat,x_add)
    }else{
      xtest=feat
    }
    return(xtest)
  }
  
  #add
  data_add=list()
  for (i in 1:length(feat_list)) {
    data_add[[i]]=my_add(b[[i]],feat_list[[i]])
  }
  names(data_add)=names(feat_list)
  
  feat_seq=sort(feat_ID)
  data_add=lapply(data_add, function(data){new_data=data[,c(feat_seq)]})
  
  return(data_add)
}

data_add <- my_pair.table(my_data$feat_list)
data_all <- rbind(data_add$PRJDB11203,data_add$PRJNA230363,data_add$PRJNA396840,
                  data_add$PRJNA678453,data_add$PRJNA717815,data_add$PRJNA932553)
rowSums(data_all)  #每个样本的所有菌种之和都是1
#### 所有菌在所有样本的相关性####  408个菌  223个样本
cor_data_all <- psych::corr.test(data_all,method = "spearman")#计算相关系数矩阵
#r值
cor_data_all_r <- cor_data_all$r
cor_data_all_r <- as.data.frame(cor_data_all_r)
ncr <- nrow(cor_data_all_r)
for (i in 1:ncr){
  for (j in 1:ncr){
    if (j>=i){
      cor_data_all_r[i,j] <- NA
    }
  }
}
cor_data_all_r <- as.data.frame(cor_data_all_r)
cor_data_all_r_long <- pivot_longer(rownames_to_column(cor_data_all_r), 
                           cols = -1,  # Assuming the first column is the identifier
                           names_to = "rowID",
                           values_to = "r") %>% 
  na.omit()
#p值
cor_data_all_p <- cor_data_all$p
cor_data_all_p <- as.data.frame(cor_data_all_p)
ncr <- nrow(cor_data_all_p)
for (i in 1:ncr){
  for (j in 1:ncr){
    if (j>=i){
      cor_data_all_p[i,j] <- NA
    }
  }
}
cor_data_all_p <- as.data.frame(cor_data_all_p)
cor_data_all_p_long <- pivot_longer(rownames_to_column(cor_data_all_p), 
                                    cols = -1,  # Assuming the first column is the identifier
                                    names_to = "rowID",
                                    values_to = "p_value") %>% 
  na.omit()

cor_data_all_pr_long_sig <- filter(cor_data_all_pr_long, p_value < 0.05 & abs(r)>0.5)


##绘制双环图数据
library(dplyr)
library(tidyr)
library(ggplot2)
#大环
diameter.1 <- 8
nodes.1 <- read.csv("Figure 2A2.csv") %>%  #读取数据
  filter(Part == "Disease") %>% #过滤数据，有两个分类，先保留大环的
  sample_n(n()) %>%#重新排序-乱序，如果有特定的顺序需要自行决定
  mutate(angle = seq(pi / 2, by = -2 * pi / nrow(.), length.out = nrow(.))) %>%  #计算等分角度
  mutate(x = (diameter.1/2) * cos(angle), #计算x坐标
         y = (diameter.1/2) * sin(angle) #计算y坐标    [大环的圆心在(0，0)]
  )
#小环
diameter.2 <- 3.5   # 小环直径
nodes.2 <- read.csv("Figure 2A2.csv") %>%
  filter(Part == "Healthy") %>%
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, by = -2 * pi / nrow(.), length.out = nrow(.))) %>%  #计算等分角度
  mutate(x = (diameter.2/2) * cos(angle) + 7, #计算x坐标
         y = (diameter.2/2) * sin(angle) + 2.25     #计算y坐标[大环的圆心在(7，2.25)]
  )

data.node.plot <- rbind(nodes.1, nodes.2)[-1] #合并节点数据

edges <- read.csv("Figure 2A3.csv", row.names = 1) #读取边数据

edges <- edges %>%
  mutate(From.x = data.node.plot$x[match(edges$source, data.node.plot$label)],
         From.y = data.node.plot$y[match(edges$source, data.node.plot$label)],
         To.x = data.node.plot$x[match(edges$target, data.node.plot$label)],
         To.y = data.node.plot$y[match(edges$target, data.node.plot$label)]) #根据节点坐标匹配连线数据

##得到Table S5
edges <- Table S5

##绘制双环
ggplot(data = data.node.plot) +
  geom_segment(data = edges, #设置连线数据
               aes_string(xend= "To.x", yend = "To.y", x = "From.x", y = "From.y", 
                          color = "factor(correlation)", linewidth = "factor(correlation)"),
               show.legend = F)+
  scale_linewidth_discrete(range = c(0.07, 0.25))+
  geom_text_repel(aes_string(x = "x", y = "y", label = "label"), # 用geom_label_repel函数绘制标签，指定位置，名称，填充颜色
                   show.legend = FALSE, size = 3, 
                  ) + # 调整标签的水平对齐方式) + # 不显示图例，指定大小，颜色，字体
  geom_point(aes_string(x = "x", y = "y", fill = "Part", size = "No..of.project"), # 映射size变量
             show.legend = F, shape = 21, color="black",stroke = 1)+ #设置节点数据
  coord_cartesian(xlim = c(-5, 15), ylim = c(-5, 5))+ # 调整坐标轴的范围
  scale_color_manual(values = c("Negative correlation" = "#21A675", 
                                "Positive correlation" = "#c93756")) + #为连线单独设置颜色
  scale_fill_manual(values = c("Disease" = "#67000D",
                               "Healthy" = "#08306B"))+  #为不同分类单独设置填充颜色
  scale_size(range = c(7, 16), guide = guide_legend()) + # 设置节点大小的范围和图例
  geom_label(data = data.frame(part = c("Disease", "Healthy"),
                               x = c(0,7),
                               y = c(0,2.25)),
             aes(label = part, x = x, y = y), size = 4.7)+ #生成圆心添加分类字体，添加背景色
  geom_label(data = data.frame(
    label = c(
      paste("Nodes:", nrow(data.node.plot)),
      paste("Edges:", nrow(edges)),
      paste("Negative edges:", sum(edges$correlation == "Negative correlation")) #统计节点、连线、正相关数据
    ),
    x = 5,
    y = c(-3,-3.5,-4)),
    aes(x = x, y = y, label = label),hjust = 0, size = 4.1)+ #添加背景色
  theme_void()  #空白绘图主题


## F2B ----
cor_pr_long_filter2 <- Table S5

from_count <- table(cor_pr_long_filter2$rowname)
from_count <- as.data.frame(from_count)
sum(from_count$Freq)
from_count <- column_to_rownames(from_count,"Var1")

to_count <- table(cor_pr_long_filter2$rowID)
to_count <- as.data.frame(to_count)
sum(to_count$Freq)
to_count <- column_to_rownames(to_count,"Var1")

feat_list <- list()
feat_list[["from_count"]] <- as.data.frame(t(from_count))
feat_list[["to_count"]] <- as.data.frame(t(to_count))
from_to_count <- my_pair.table(feat_list)
identical(colnames(from_to_count$from_count),colnames(from_to_count$to_count))
from_to_edge_count <- rbind(from_to_count$from_count,from_to_count$to_count)
from_to_edge_count <- as.data.frame(t(from_to_edge_count))
from_to_edge_count$sum <- rowSums(from_to_edge_count)
#添加group
from_to_edge_count$group <- "Health"
row_names <- rownames(from_to_edge_count)
from_to_edge_count$group[row_names %in% rownames(data_all_marker_P)] <- "Disease"
#添加degree
# 创建一个新的列degree并初始化为0
from_to_edge_count$degree <- 0
# 根据group的值计算degree
from_to_edge_count$degree[from_to_edge_count$group == "Disease"] <- round(from_to_edge_count$sum[from_to_edge_count$group == "Disease"] / 42,2)
from_to_edge_count$degree[from_to_edge_count$group == "Health"] <- round(from_to_edge_count$sum[from_to_edge_count$group == "Health"] / 12,2)


from_to_edge_count2 <- from_to_edge_count
from_to_edge_count2 <- rownames_to_column(from_to_edge_count2,"species")
from_to_edge_count2 <- arrange(from_to_edge_count2,group,sum)
# 选出健康组的数据
health_data <- filter(from_to_edge_count2, group == "Health")

# 给sum列乘以-1
health_data$sum <- health_data$sum * -1

# 把修改后的数据重新合并到原数据框中
from_to_edge_count2 <- rbind(filter(from_to_edge_count2, group == "Disease"), health_data)

ggplot(from_to_edge_count2, aes(x = reorder(species, -sum), y = sum)) +
  geom_bar(stat = "identity",aes(fill = group), width = 0.7,
           show.legend = T) +
  geom_text(aes(label = sum), hjust = 1.2, vjust = 0.5, size = 3) + # 添加这一行来标注sum
  #coord_flip() +
  scale_fill_manual(values = c("Disease" = "#C45C69","Health" = "#4D779B"))+
  theme_bw() + theme(panel.grid=element_blank(),#设置主题
                     panel.border = element_blank(),
                     axis.ticks.y = element_blank()) + #去掉外边边框
  #labs(x = "Scientific Name", y = "Project ID", fill = "LDA Score") +  #设置坐标轴标签
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  #设置坐标轴角度
  theme(legend.position = "right") #设置图例位置

## F2C ----
cor_pr_long_filter2 <- Table S5

from_count <- table(cor_pr_long_filter2$rowname)
from_count <- as.data.frame(from_count)
sum(from_count$Freq)
from_count <- column_to_rownames(from_count,"Var1")


to_count <- table(cor_pr_long_filter2$rowID)
to_count <- as.data.frame(to_count)
sum(to_count$Freq)
to_count <- column_to_rownames(to_count,"Var1")

feat_list <- list()
feat_list[["from_count"]] <- as.data.frame(t(from_count))
feat_list[["to_count"]] <- as.data.frame(t(to_count))
from_to_count <- my_pair.table(feat_list)
identical(colnames(from_to_count$from_count),colnames(from_to_count$to_count))
from_to_edge_count <- rbind(from_to_count$from_count,from_to_count$to_count)
from_to_edge_count <- as.data.frame(t(from_to_edge_count))
from_to_edge_count$sum <- rowSums(from_to_edge_count)
#添加group
from_to_edge_count$group <- "Health"
row_names <- rownames(from_to_edge_count)
from_to_edge_count$group[row_names %in% rownames(data_all_marker_P)] <- "Disease"
#添加degree
# 创建一个新的列degree并初始化为0
from_to_edge_count$degree <- 0
# 根据group的值计算degree
from_to_edge_count$degree[from_to_edge_count$group == "Disease"] <- round(from_to_edge_count$sum[from_to_edge_count$group == "Disease"] / 42,2)
from_to_edge_count$degree[from_to_edge_count$group == "Health"] <- round(from_to_edge_count$sum[from_to_edge_count$group == "Health"] / 12,2)
#箱线图
# 载入ggplot2包
library(ggplot2)
library(ggpubr)
# 用ggplot函数创建一个图层对象
ggplot(from_to_edge_count, aes(x = group, y = degree, fill = group))+ 
  geom_boxplot(width = 0.18,position = position_dodge(width = 1)) + 
  scale_fill_manual(values = c("#EE0000","#3B4992"),name = "Group")+
  #facet_wrap(~ subject, ncol = 3,scales = "free") +# 用facet_wrap函数按照diversity分面
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test",position = position_nudge(x = 0.5)) +
  theme_minimal() + theme(panel.grid=element_blank())+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")


## F2D ----
####连接强度大小  r的绝对值 箱线图
cor_pr_long_filter2 <- Table S5

onnectivity_species_P <- cor_pr_long_filter2[cor_pr_long_filter2$rowname %in% rownames(data_all_marker_P),]
connectivity_species_P <- connectivity_species_P[connectivity_species_P$rowID %in% rownames(data_all_marker_P),]
connectivity_species_P <- edges_filtered_P[,c(1,2,3)]
colnames(connectivity_species_P) <- c("species","group","r")
connectivity_species_P$group <- "disease"
connectivity_species_H <- edges_filtered_H[,c(1,2,3)]
colnames(connectivity_species_H) <- c("species","group","r")
connectivity_species_H$group <- "health"
connectivity_species <- rbind(connectivity_species_P,connectivity_species_H)

# 载入ggplot2包
library(ggplot2)
library(ggpubr)
# 用ggplot函数创建一个图层对象
ggplot(connectivity_species, aes(x = group, y = r, fill = group))+ 
  geom_boxplot(width = 0.18,position = position_dodge(width = 1)) + 
  scale_fill_manual(values = c("#EE0000","#3B4992"),name = "Group")+
  #facet_wrap(~ subject, ncol = 3,scales = "free") +# 用facet_wrap函数按照diversity分面
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test") +
  theme_bw() + theme(panel.grid=element_blank())+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
