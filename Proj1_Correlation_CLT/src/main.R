# 导入配置
Root_path = 'E:\\Fin_Lab\\Proj1_Correlation_CLT'
config_path = paste0(Root_path,'\\','config')
setwd(config_path)
source('settings.R')

# 计算协方差函数
setwd(lib_path)
source('Caculate_covariance.R') # 计算根据模拟轨道得到的协方差函数
setwd(Result_path)
save(data_list,file = 'result_1')

# 计算极限分布
E_sigma_A = sigma_A_kappa_sum/(path_num-length(wenti))
E_sigma_B = sigma_B_kappa_sum/(path_num-length(wenti))
E_B_kappa = B_kappa_sum/(path_num-length(wenti))
C_ij_kappa = C_ij_kappa_sum/(path_num-length(wenti))

g_rho_kappa_hat = apply(g_rho_kappa_mat,1,mean)
yita_hat = sqrt(Delta*sum(g_rho_kappa_hat^2))
f_rho_kappa_hat = g_rho_kappa_hat/yita_hat

## 计算Cov_tilde
Cov_tilde =  
  (E_B_kappa%*%t(E_B_kappa)/(4*sqrt(E_sigma_A^3%*%t(E_sigma_A^3))*
                                         sqrt(E_sigma_B%*%t(E_sigma_B))))*C_ij_kappa[1,1,,]+
  (E_B_kappa%*%t(E_B_kappa)/(4*sqrt(E_sigma_A^3%*%t(E_sigma_A))*
                               sqrt(E_sigma_B%*%t(E_sigma_B^3))))*C_ij_kappa[1,2,,]-
  apply(C_ij_kappa[1,3,,],2,function(x) x*E_B_kappa)/(2*sqrt(E_sigma_A^3%*%t(E_sigma_A))*
                                                        sqrt(E_sigma_B%*%t(E_sigma_B)))+
  (E_B_kappa%*%t(E_B_kappa)/(4*sqrt(E_sigma_A%*%t(E_sigma_A^3))*
                               sqrt(E_sigma_B^3%*%t(E_sigma_B))))*C_ij_kappa[2,1,,]+
  (E_B_kappa%*%t(E_B_kappa)/(4*sqrt(E_sigma_A%*%t(E_sigma_A))*
                               sqrt(E_sigma_B^3%*%t(E_sigma_B^3))))*C_ij_kappa[2,2,,]-
  apply(C_ij_kappa[2,3,,],2,function(x) x*E_B_kappa)/(2*sqrt(E_sigma_A%*%t(E_sigma_A))*
                                                        sqrt(E_sigma_B^3%*%t(E_sigma_B)))-
  apply(C_ij_kappa[3,1,,],2,function(x) x*E_B_kappa)/(2*sqrt(E_sigma_A%*%t(E_sigma_A^3))*
                                                        sqrt(E_sigma_B%*%t(E_sigma_B)))-
  apply(C_ij_kappa[3,2,,],2,function(x) x*E_B_kappa)/(2*sqrt(E_sigma_A%*%t(E_sigma_A))*
                                                        sqrt(E_sigma_B%*%t(E_sigma_B^3)))+
  C_ij_kappa[3,3,,]/(sqrt(E_sigma_A%*%t(E_sigma_A))*sqrt(E_sigma_B%*%t(E_sigma_B)))
  
## 计算Cov_brave
Cov_brave_term1_abbr = rep(NA,n)
for(kappakappa in 1:n){
  Cov_brave_term1_abbr[kappakappa] = sum(Delta*g_rho_kappa_hat*Cov_tilde[kappakappa,])
}

Cov_brave_term2_abbr = sum(apply(Cov_tilde,1,function(x) sum(x*g_rho_kappa_hat*Delta))*
                             g_rho_kappa_hat*Delta)

Cov_brave = Cov_tilde - (f_rho_kappa_hat/yita_hat)%*%t(Cov_brave_term1_abbr)-
  t((f_rho_kappa_hat/yita_hat)%*%t(Cov_brave_term1_abbr))+
  (f_rho_kappa_hat%*%t(f_rho_kappa_hat)/yita_hat)*Cov_brave_term2_abbr

## 生成Z
Eigen_list = eigen(Cov_brave/(yita_hat^2))
