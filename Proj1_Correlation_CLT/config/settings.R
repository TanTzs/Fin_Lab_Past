# 载入相关的包
library(MASS)
library(Rfast)
library(data.table)

# 路径设置

## 数据路径
Data_A_path<-paste0(Root_path,'\\','data\\data_test1','\\','price_1')
Data_B_path<-paste0(Root_path,'\\','data\\data_test1','\\','price_2')
g_rho_path<-'D:\\tyw_part\\corr_CLT\\mean'
Result_path<-paste0(Root_path,'\\','assets','\\','result_1')

## 计算脚本路径
lib_path<-paste0(Root_path,'\\','lib')

# 参数设置
path_num=5 ## 模拟数据轨道数
l=10 ## 估计窗宽
n=288*5 ## 观测数
data_n = 1440 ## 生成数据每天的观测数
Delta=1/n 
kappa= seq(0,1,length.out = n)
Ln=5
N=0.8*n
N<-floor(N)
file_A_names<-paste0('LogPrice_I_Iter',1:path_num,'.txt')
file_B_names<-paste0('LogPrice_II_Iter',1:path_num,'.txt')

# 依赖函数
frequen_adj<-function(x,n){ # 调整观测值过多的数据为固定观测值
  b<-rep(NA,n) 
  batch_size = floor(length(x)/n)
  for(i in 1:n){
    b[i]<-sum(x[((i-1)*batch_size+1):(i*batch_size)])
  } 
  return(b)
}
