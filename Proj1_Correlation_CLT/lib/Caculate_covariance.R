# Caculate covariance function

# 保存数据
B_kappa_sum<-rep(0,nrow=length(kappa)) 
sigma_A_kappa_sum<-rep(0,length(kappa))
sigma_B_kappa_sum<-rep(0,length(kappa))
C_ij_kappa_sum<-array(0,dim = c(3,3,n,n))

g_rho_kappa_mat=matrix(NA,nrow=n,ncol=path_num) # 记录每条轨道
wenti = c()

# 逐个轨道计算
for(path_iternum in 1:path_num){
  
  print(paste0('processing: ',path_iternum,'/',path_num))
  
  iter_time_start = proc.time()
  # 读取数据
  data_A<-read.table(paste(Data_A_path,file_A_names[path_iternum],sep='\\'))
  data_A[,1]<-as.numeric(data_A[,1])
  data_A<-matrix(diff(log(data_A[,1]))[1:(data_n*N)],ncol=N) # 得到收益率数据
  data_A<-apply(data_A,2,frequen_adj,n=n)
  
  data_B<-read.table(paste(Data_B_path,file_B_names[path_iternum],sep='\\'))
  data_B[,1]<-as.numeric(data_B[,1])
  data_B<-matrix(diff(log(data_B[,1]))[1:(data_n*N)],ncol=N) # 得到收益率数据
  data_B<-apply(data_B,2,frequen_adj,n=n)
  
  # 曲线估计
  B_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N) 
  sigma_A_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N)
  sigma_B_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N)
  
  estcurve_time_start = proc.time()
  for(i in 1:N){
    if(i == 1){
      for(kappa_iternum in 1:length(kappa)){
        j_kappa<-floor(kappa[kappa_iternum]*n)
        j<-j_kappa-l+1
        
        ## local spot covariation and volatility
        if(j>=1){
          B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]*data_B[(j:j_kappa),i])
          sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]^2)
          sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_B[(j:j_kappa),i]^2)
        }else{
          B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_A[(1:l),i]*data_B[(1:l),i])
          )
          sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_A[(1:l),i]^2)
          )
          sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_B[(1:l),i]^2)
          )
        }
      }
    }else{
      # 对于i>=2的天数，可用前一天的数据进行估计
      for(kappa_iternum in 1:length(kappa)){
        
        j_kappa<-floor(kappa[kappa_iternum]*n)
        j<-j_kappa-l+1
        
        ## local spot covariation and volatility
        if(j>=1){
          B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]*data_B[(j:j_kappa),i])
          sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]^2)
          sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_B[(j:j_kappa),i]^2)
        }else{
          B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_A[(1:j_kappa),i]*data_B[(1:j_kappa),i])+ # 此为当天
              sum(data_A[((n+j):n),i-1]*data_B[((n+j):n),i-1]) # 此为前一天
          )
          sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_A[(1:j_kappa),i]^2)+ # 此为当天
              sum(data_A[((n+j):n),i-1]^2) # 此为前一天
          )
          sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
            sum(data_B[(1:j_kappa),i]^2)+ # 此为当天
              sum(data_B[((n+j):n),i-1]^2) # 此为前一天
          )
        }
      }
    }
  }
  estcurve_time_end = proc.time()
  # correlation diurnal pattern function estimate
  g_rho_kappa<-apply(B_i_kappa,1,mean)/(
    (apply(sigma_A_i_kappa,1,mean)^(1/2))*(apply(sigma_B_i_kappa,1,mean)^(1/2))
  )
  
  # Approximation of the Limiting Distribution
  
  ## Covariance function (matrix-valued function)
  A_m_kappa<-array(NA,dim=c(3,n,N))
  A_m_kappa[1,,]<-sigma_A_i_kappa-apply(sigma_A_i_kappa,1,mean) # A_1
  A_m_kappa[2,,]<-sigma_B_i_kappa-apply(sigma_B_i_kappa,1,mean) # A_2
  A_m_kappa[3,,]<-B_i_kappa-apply(B_i_kappa,1,mean) # A_3
  
  C_ij_kappa<-array(0,dim = c(3,3,n,n))
  
  estcov_time_start = proc.time()
  for(i in 1:3){
    for(j in i:3){
      component_1_mat<-(1/N)*A_m_kappa[i,,]%*%t(A_m_kappa[j,,])
      component_2_mat<-matrix(0,ncol=ncol(component_1_mat),
                              nrow=nrow(component_1_mat))
      for(h in 1:Ln){
        a2<-A_m_kappa[i,,1:(N-h)]%*%t(A_m_kappa[j,,(1+h):N])+
          A_m_kappa[i,,(1+h):N]%*%t(A_m_kappa[j,,1:(N-h)])
        component_2_mat<-component_2_mat+a2/(N-h)
      }
      C_ij_kappa[i,j,,]<-component_1_mat+component_2_mat
    }
  }

  for(i in 1:3){
    for(j in i:3){
      if(i==j){
        NULL
      }else if(j>i){
        C_ij_kappa[j,i,,]<-t(C_ij_kappa[i,j,,])
      }
    }
  }
  estcov_time_end = proc.time()
  
  # 保存变量
  if(identical(which(is.na(C_ij_kappa)), integer(0))){
    B_kappa_sum<-B_kappa_sum+apply(B_i_kappa,1,mean)
    sigma_A_kappa_sum<-sigma_A_kappa_sum+apply(sigma_A_i_kappa,1,mean)
    sigma_B_kappa_sum<-sigma_B_kappa_sum+apply(sigma_B_i_kappa,1,mean)
    C_ij_kappa_sum<-C_ij_kappa_sum+C_ij_kappa 
    
    g_rho_kappa_mat[,path_iternum]<-g_rho_kappa
  }else{
    B_kappa_sum<-B_kappa_sum
    sigma_A_kappa_sum<-sigma_A_kappa_sum
    sigma_B_kappa_sum<-sigma_B_kappa_sum
    
    g_rho_kappa_mat[,path_iternum]<-rep(NA,length(kappa))
    C_ij_kappa_sum<-C_ij_kappa_sum
    wenti=c(wenti,path_iternum) # 保留出问题的轨道
  }
  iter_time_end = proc.time()
  cat(
    ' 曲线估计所用时间:',(estcurve_time_end-estcurve_time_start)[1:3],'\n',
    '协方差估计所用时间:',(estcov_time_end-estcov_time_start)[1:3],'\n',
    '此轨道总耗时:',(iter_time_end-iter_time_start)[1:3],'\n'
      )
}

# 保留变量数据
data_list=list(
  'B_kappa_sum'=B_kappa_sum,
  'sigma_A_kappa_sum'=sigma_A_kappa_sum,
  'sigma_B_kappa_sum'=sigma_B_kappa_sum,
  'C_ij_kappa_sum'=C_ij_kappa_sum,
  'g_rho_kappa_mat'=g_rho_kappa_mat,
  'wenti' = wenti,
  'path_num' = path_num,
  'l' = l, 
  'n' = n ,
  'data_n' = data_n,
  'Delta' = Delta,
  'kappa' = kappa,
  'Ln' = Ln,
  'N' = N
)



