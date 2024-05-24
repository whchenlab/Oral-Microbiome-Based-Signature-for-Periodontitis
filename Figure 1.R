## F1B ----
#加载包
library(ggplot2)
library(forcats)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)
load("Figure 01.Rdata")
feat_list <- my_data$feat_list
#读取数据
marker.adj=my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,is_plot=F,lda_cutoff=2,nproj_cutoff=1,level="species",change_name=F)
##Figure 1B.Rdata即为marker.adj  标志菌lefse结果 
load("Figure 1B.Rdata")
#处理数据
markers_data <- marker.adj$marker_data
markers_data$scientific_name <- gsub("s__","", markers_data$scientific_name)
class(markers_data$LDA)
markers_data$LDA <- as.numeric(markers_data$LDA)
#根据n_project(>=3)和LDA()>=2值进行过滤
markers_data_filter <- markers_data %>%
  filter(nrproj >= 3, abs(LDA) >= 2) %>% 
  filter(class == "adjust")
unique(markers_data_filter$scientific_name)  #57个
##绘制图处理数据
#计算菌种在项目的数目
dat <- markers_data_filter
tab <- table(dat$scientific_name, dat$project_id)

count <- rowSums(tab)
count <- data.frame(scientific_name = names(count), count = count, row.names = NULL)

#菌种名称按照LDA平均值进行排序
dat_mean <- dat %>%
  group_by(scientific_name) %>%
  summarise(mean_LDA = mean(LDA))
#class(dat_mean$scientific_name)
#dat_mean <- arrange(dat_mean,desc(mean_LDA))
#合并俩数据框
count_dat_mean <- merge(count,dat_mean)
class(count_dat_mean$count)
class(count_dat_mean$mean_LDA)
#先把count_dat_mean排序
count_dat_mean <- count_dat_mean[order(-count_dat_mean$mean_LDA),]
#提取健康和疾病,先排序
count_dat_mean_H <- count_dat_mean[count_dat_mean$mean_LDA < 0,]
count_dat_mean_H <- count_dat_mean_H[order(count_dat_mean_H$count,-count_dat_mean_H$mean_LDA),]
count_dat_mean_P <- count_dat_mean[count_dat_mean$mean_LDA > 0,]
count_dat_mean_P <- count_dat_mean_P[order(-count_dat_mean_P$count,-count_dat_mean_P$mean_LDA),]
#合并到一起
count_dat_mean_PH <- rbind(count_dat_mean_P,count_dat_mean_H)
class(count_dat_mean_PH)
class(count_dat_mean_PH$scientific_name)
#dat数据框的菌种按照合并一起的数据框菌种进行排序
dat$scientific_name <- factor(dat$scientific_name,
                              levels=rev(count_dat_mean_PH$scientific_name))
##绘制图完整代码
# 设置图形的宽度和高度为3英寸，分辨率为300dpi，不居中放置
pdf("F1B.pdf", width = 7, height = 12, pagecentre = FALSE)

定义颜色范围                                           
my_colors <- colorRampPalette(c("#08306B", "#4292C6", "#9ECAE1", "#C6DBEF", "#F6BDC0", "#F1959B", "#EA4C46", "#781426"))
#绘制热图

p1 <- ggplot(data = dat, aes(y = scientific_name, x = project_id, fill = LDA)) +
  geom_tile()  + #绘制单元格
  scale_fill_gradientn(colours = my_colors(20), breaks = c(-4, -3, -2, -1, 1, 2, 3, 4), 
                       labels = c("< -4", "-3", "-2", "-1", "1", "2", "3", "> 4")) + #设置颜色填充
  theme_bw() + theme(panel.grid=element_blank(),#设置主题
                     panel.border = element_blank(),
                     axis.ticks.y = element_blank()) + #去掉外边边框
  #labs(x = "Scientific Name", y = "Project ID", fill = "LDA Score") +  #设置坐标轴标签
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  #设置坐标轴角度
  theme(legend.position = "left") #设置图例位置
#coord_equal() #设置x轴和y轴比例为1:1
#柱状图
p2 <- ggplot(data = count_dat_mean_PH,aes(y=factor(scientific_name, levels=rev(unique(scientific_name))),x=count,fill=count)) +
  geom_bar(stat = 'identity',
           fill=ifelse(count_dat_mean_PH$mean_LDA>0,"#EA4C46","#4292C6"),
           width = 0.8) +#根据LDA平均值正负设置颜色
  coord_cartesian(xlim = c(0, 6))+ # 设置ylim参数为-2到2
  labs(y='',x='num_proj')+
  #coord_flip() + #横纵翻转的
  theme(panel.grid=element_blank(),
        panel.grid.major =element_blank(),  # 去除边框
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  #设置主题

# 把热图和柱状图组合在一起
plot_grid(p1, p2, align = "h",ncol = 2, rel_widths =c(0.9, 0.2),hjust = -1)

# 关闭pdf设备
dev.off()

## F1C ----
load("Figure 01.Rdata")
source("GMModels/CCM & SCM Fun.R")
source("GMModels/siamcat_models.R")
source("GMModels/siamcat_models_adj.R")
source("GMModels/my_lodo.R")
source("GMModels/my_lodo.adj.R")
method='randomForest' #method:c("lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest")
label='randomForest_adj_nc'
method_dir <- 'randomForest'
feat_list <- my_data$feat_list
meta_list <- my_data$meta_list

## CV_LODO建模 ----
minmax_normal_randomForest_CV_model.adj <- my_result_adj(feat_list,meta_list,
                              method=method,label=label,
                              models=NULL,
                              num.folds=5, num.resample=3,
                              feature.type = "normalized",
                              is_cross=T,do.con=T,
                              pca_plot = F,imbalance = T,max_index=3)#is_cross=T批次校正


##Figure 1C.Rdata即为minmax_normal_randomForest_CV_model.adj 建模结果
load("Figure 1C.Rdata")
#提取result_CV
result_CV <- minmax_normal_randomForest_CV_model.adj$result
result_CV <- result_CV[,-c(1,2)]
result_CV <-  as.data.frame(result_CV)
# 用美元符号访问数据框的r_average列,行平均值
result_CV$r_average <- vector(length = 6)

for (i in 1:6){
  result_CV$r_average[i] <- round((sum(result_CV[i,])-result_CV[i,i])/5,2)
}
#列平均值
result_CV["c_average",] <- vector(length = 7)
for (i in 1:6){
  result_CV["c_average",][i] <- round((sum(result_CV[,i])-result_CV[i,i])/5,2)
}
#行列对角线average处的平均值
sum = 0
for (i in 1:6){
  sum <- sum + result_CV[i,i]
}
sum
result_CV["c_average",][7] <- round(sum/6,2)

library(purrr)
minmax_normal_randomForest_LODO_model.adj <- my_lodo.adj(feat=feat_list,meta=meta_list,top=NULL,models=NULL,method=method,
                              num.folds=10, num.resample=3,
                              check.con.before=NULL,check.con.after=NULL,
                              add_group=T,is_combat=T,re_scale=T,is_cross=T,
                              pca_plot=F,imbalance=F,max_index=3,batch_group=T,do_adj=T,
                              do.fs=F,nest_top=NULL)#is_cross=T批次校正

#提取result_LODO
result_LODO <- minmax_normal_randomForest_LODO_model.adj$result
result_LODO <- result_LODO[-1,] %>% round(.,2)
result_LODO$r_average <- round(rowMeans(result_LODO, na.rm = TRUE),2)
#合并result_CV  result_LODO
merged_df <- rbind(result_CV,result_LODO)

df <- data.frame(merged_df)
df$subject <- factor(rownames(df), levels = c("PRJDB11203","PRJNA230363", "PRJNA396840","PRJNA678453", "PRJNA717815","PRJNA932553","c_average","lodo"))  
df_melt <- melt(df)
df_melt$subject <- factor(df_melt$subject, levels = c( "lodo","c_average","PRJNA932553", "PRJNA717815","PRJNA678453", "PRJNA396840","PRJNA230363","PRJDB11203"))
df_melt$variable <- factor(df_melt$variable, levels = c("PRJDB11203","PRJNA230363", "PRJNA396840","PRJNA678453", "PRJNA717815","PRJNA932553","r_average"))
##绘制热图
ggplot(df_melt, aes(x = variable, y = subject, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value),size = 5,colour = "black") +
  scale_fill_gradientn(colors = c("white", "red"), limits = c(0.2, 1.0)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.3),
        axis.text.y = element_text(hjust = 1)) +
  scale_x_discrete(position = "top") +# 把 x 轴的位置设置为上面
  xlab("testing cohorts") + # 添加横坐标标题
  ylab("training cohorts") + # 添加纵坐标标题
  theme(axis.title.y = element_text(margin = ggplot2::margin(0,20,0,0)),# 设置y轴标题和图的距离
        axis.title.x = element_text(margin = ggplot2::margin(50,0,0,0))) + # 设置x轴标题和图的距离
  geom_segment(
    x = 0 ,       # 分割线起点x坐标
    xend = as.numeric(max(as.numeric(df_melt$variable)))+0.5 ,    # 分割线终点x坐标
    y = as.numeric(min(as.numeric(df_melt$subject)))+0.5 ,       # 分割线y坐标
    yend = as.numeric(min(as.numeric(df_melt$subject)))+0.5 ,    # 分割线终点y坐标
    color = "white", linewidth = 5
  )+
  geom_segment(
    x = 0 ,       # 分割线起点x坐标
    xend = as.numeric(max(as.numeric(df_melt$variable)))+1.5 ,    # 分割线终点x坐标
    y = as.numeric(min(as.numeric(df_melt$subject)))+1.5 ,       # 分割线y坐标
    yend = as.numeric(min(as.numeric(df_melt$subject)))+1.5 ,    # 分割线终点y坐标
    color = "white", linewidth = 5
  )+  
  geom_segment(
    y = 0 ,       # 分割线起点y坐标
    yend = as.numeric(max(as.numeric(df_melt$variable)))+1.5 ,    # 分割线终点y坐标
    x = as.numeric(min(as.numeric(df_melt$subject)))+5.5 ,       # 分割线x坐标
    xend = as.numeric(min(as.numeric(df_melt$subject)))+5.5 ,    # 分割线终点x坐标
    color = "white", linewidth = 5
  )+ 
  theme(panel.background = element_blank())


## F1D ----
##读取数据
prediction_accuracy <- read.csv("Figure 1D.csv")

library(ggplot2)
prediction_accuracy_long <- reshape2::melt(prediction_accuracy)

  
# 绘制热图表格
ggplot(prediction_accuracy_long, aes(x = variable , y = X, fill = variable)) +
  geom_tile() +
  # 添加数值
  geom_text(aes(label = value), color = "white",size = 8) +
  labs(fill = "variable") +
  theme_bw() + theme(panel.grid=element_blank(),#设置主题
                     panel.border = element_blank(),
                     axis.ticks.y = element_blank()) + #去掉外边边框
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        # 设置上方坐标轴文本的角度和对齐方式
        axis.text.x.top = element_text(angle = 45, hjust = 1)) +  
  theme(legend.position = "right") +#设置图例位置
  # 使用数据中的颜色值
  scale_fill_manual(values = c("PRJDB6966" = "#83D475", "PRJNA552294" = "#F1959B")) +
  # 使用离散型x轴
  scale_x_discrete(position = "top")


## F1E ----
##读取数据
GL_accuracy <- read.csv("Figure 1E.csv")

##绘制小提琴图
ggplot(GL_accuracy,aes(x=Group,y=Accuracy,fill=Group)) +
 geom_violin(trim=FALSE,color="white") + #绘制小提琴图
 geom_boxplot(width=0.2,position=position_dodge(0.9))+
 stat_compare_means(aes(group=group),
 label = "p.signif",label.x = 1.5)+#添加P值
 labs(title = "Alpha diversity",x="",y="")+ #图片title x轴和y轴标题
 theme_classic()





