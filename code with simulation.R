library(psych)
library(matrixcalc)
library(maotai)
library(expm)
library(imager)

#open image ---------------------------------------------------------
#library(imager)
im <- load.image("flag_image_300x300.jpg")
str(im);
plot(im) 
grid(nx = 15, ny = 15, col = "black")
grid(nx = 15, ny = 15, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

#create vector (it will be by rows) of an image:
imgtest2<-as.data.frame(im) #convert to data.frame

imv<-imgtest2$value #there are 3 consecutive vectors of each cc in a 'value' column 

#as the picture has no color, they are 3 equal vectors and we need only one
len<-length(imv)/3;
imvByRows<-imv[1:len]


#amount of columns and rows
k<-sqrt(len);

#create matrix from vector
imagematrix<-matrix(imvByRows,nrow=k,byrow=TRUE);
#save(imagematrix, file = "imagematrix.Rdata") #-------------save


#range(imagematrix)



#divide on ROIs and create vector ---------------------------------------------------------
n<-k/20 
ROI_vector<-c()

for(j in 1:n){
  for(i in 1:n){
    ROI_vector[(i+(j-1)*n)]<-mean(imagematrix[(20*(i-1)+1):(20*i),(20*(j-1)+1):(20*j)])
  }
}

#save(ROI_vector, file = "ROI_vector.Rdata") #-----------save

#make changes ---------------------------------------------------------------------------------
changeM<-imagematrix

changeM[120:180, 120:180]<-changeM[120:180, 120:180]+0.005
#changeM[11:25, 11:25]<-changeM[11:25, 11:25]+0.05
#changeM[21:35, 21:35]<-changeM[21:35, 21:35]+0.005
#changeM[201:215, 201:215]<-changeM[201:215, 201:215]+0.005
#changeM<-changeM-0.0008

changed_vector<-c()
for(j in 1:n){
  for(i in 1:n){
    changed_vector[(i+(j-1)*n)]<-mean(changeM[(20*(i-1)+1):(20*i),(20*(j-1)+1):(20*j)])
  }
}
#save(changed_vector, file = "changed_vector+0_005(big).Rdata") 

# matrix W ---------------------------------------------------------
W <- matrix(0,n^2,n^2);

help_matrix <- matrix( c( rep(1:n,n), rep(1:n,rep(n,n))  ),n^2,2);
for( i in 1:n^2 ){
  for( j in i:n^2 ){
    distance <- sqrt( (help_matrix[i,][1]-help_matrix[j,][1])^2 + (help_matrix[i,][2]-help_matrix[j,][2])^2 );
    if (distance <= sqrt(18) ) { W[i,j] <- W[j,i] <- 1 };
  }
}
for (j in 1:n){
  W[j,j]<-0
}
for( k in 1:nrow(W) ) { W[k,] <- W[k,]/sum(W[k,]) };

#save(W, file = "new_W.Rdata") #-----------save

#--------------------------------------------------------------------------------------
A<-diag(0.5,n*n)
delt<-0.01
I<-diag(1,n*n)
I_minus_delt_W_inv<-solve(I-delt*W)
#library(expm)
Cov_of_residuals_begin<-sqrtm(I_minus_delt_W_inv%*%diag(0.005*0.005, n*n)%*%t(I_minus_delt_W_inv))
Cov_of_residuals_correct_ROI<-solve(Cov_of_residuals_begin)

#Gamma0---------------------------------------------------------------------------------

Phi<-I_minus_delt_W_inv%*%A
Phi_inv<-solve(Phi)
Minus_transp_Phi<-t(Phi)*(-1)
C_for_sylv<-Phi_inv%*%I_minus_delt_W_inv%*%diag(0.005*0.005, n*n )%*%t(I_minus_delt_W_inv)
#library(maotai)
Gamma0<-sylvester(Phi_inv,Minus_transp_Phi,C_for_sylv)


#function to calculate Gammas ---INDEXES ARE FROM 1---------------------------------
calculate_GAMMAs<-function(k, Gamma0, Phi){
  GAMMAs<-list()
  GAMMAs[[1]]<-Gamma0
  for(i in 1:k-1){
    GAMMAs[[i+1]]<-matrix.power(Phi,i)%*%Gamma0
  }
  return(GAMMAs)
}
#-----------------------------------------------------------------------

GAMMAs1_120<-calculate_GAMMAs(120, Gamma0, Phi)
#save(GAMMAs1_120, file = "GAMMAs1_120.Rdata")

#--------------------------------------Sigmas----------------------------------------------


#function to calculate Sigma_t_r------------------------------------
calculate_Sigma_t_r<-function(t,GAMMAs, lambda, n){
  I_minus_lambda<-(1-lambda)
  matrix_sum<-rep(0,(n*n))
  for(i in 0:(t-1)){
    for(j in 0:(t-1)){
      a<-diag((I_minus_lambda^i),(n*n))
      b<-diag((I_minus_lambda^j),(n*n))
      matrix_sum<-matrix_sum+(a%*%GAMMAs[[abs(i-j)+1]]%*%b)
    }
  }
  res<-diag(lambda,n*n)%*%matrix_sum%*%diag(lambda,n*n)
  return(res)
}
#---------------------------------------------------------------------

#library(matrixcalc)
lambda<-0.5
Sigmas<-list()
for (i in 1:120){
  Sigmas[[i]]<-calculate_Sigma_t_r(i,GAMMAs1_120, lambda, n)
}

#save(Sigmas, file = "Sigmas1_120_L05.Rdata")

#after some point all further Sigmas will be the same

#---------------------------------------------------------------------

#sigmas_l_r for CS_2
big_F<-I_minus_delt_W_inv%*%A
lambda<-0.2
Sigma_l_r_by_formula<-(lambda/(2-lambda))*(solve(I-(1-lambda)*big_F)%*%Gamma0+Gamma0%*%solve(I-(1-lambda)*t(big_F))-Gamma0)
#save(Sigma_l_r_by_formula, file = "Sigma_l_r_by_formula_lambda_0_2.Rdata")


lambda<-0.8
Sigma_l_r_by_formula<-(lambda/(2-lambda))*(solve(I-(1-lambda)*big_F)%*%Gamma0+Gamma0%*%solve(I-(1-lambda)*t(big_F))-Gamma0)
#save(Sigma_l_r_by_formula, file = "Sigma_l_r_by_formula_lambda_0_8.Rdata")


lambda<-0.5
Sigma_l_r_by_formula<-(lambda/(2-lambda))*(solve(I-(1-lambda)*big_F)%*%Gamma0+Gamma0%*%solve(I-(1-lambda)*t(big_F))-Gamma0)
#save(Sigma_l_r_by_formula, file = "Sigma_l_r_by_formula_lambda_0_5.Rdata")


lambda<-1
Sigma_l_r_by_formula<-(lambda/(2-lambda))*(solve(I-(1-lambda)*big_F)%*%Gamma0+Gamma0%*%solve(I-(1-lambda)*t(big_F))-Gamma0)
#save(Sigma_l_r_by_formula, file = "Sigma_l_r_by_formula_lambda_1.Rdata")
#----------------------------------------------------------------------
  

#--------------------------------------Simulation--------------------------------------------

Lambda<-diag(0.5, n*n) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#load("changed_vector+0.Rdata")
#load("changed_vector+0_02.Rdata")



#CS_1_t_r-----------------------------------------------------------------------------------
y<- list()
z<- list()

Z_t_minus_mu_minus_1<-rep(0,(n*n))

 
for (i in 1:20){ #some steps in the beginning should be excluded
  
  if (i==1){
    y[[1]]<-ROI_vector
    z[[1]]<-ROI_vector
  }else{
    eps<-rnorm(n*n , mean=0 , sd=0.005)
    z[[i]]<-y[[1]]+I_minus_delt_W_inv%*%A%*%(y[[i-1]]-y[[1]])
    y[[i]]<-z[[i]]+I_minus_delt_W_inv%*%eps
  }
  
  
} 

CS_1_t_r<-c()
for (t in 21:120){
  
  eps<-rnorm(n*n , mean=0 , sd=0.005)
  z[[t]]<-z[[1]]+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-z[[1]])
  y[[t]]<-changed_vector+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-changed_vector)+I_minus_delt_W_inv%*%eps
  
  
  Z_t_minus_mu<-(I-Lambda)%*%Z_t_minus_mu_minus_1+Lambda%*%(y[[t]]-y[[1]])
  Z_t_minus_mu_minus_1<-Z_t_minus_mu
  
  CS_1_t_r[t-20]<-(t(Z_t_minus_mu)%*%solve(Sigmas[[t-20]])%*%Z_t_minus_mu-n*n)/(sqrt(2*n*n))
}

plot(CS_1_t_r, type = "b", lwd=1.5)
abline(h= 2.337, col = "red") #control limit was estimated separately for ARL_0=100



#CS_3_t_r res-----------------------------------------------------------------------------------

lambda<-0.5
y<- list()
z<- list()

eta_t<-list()
CS_3t_res<-c()
Z_t_minus_1<-rep(0,(n*n)) 
O_t<-c() 


for (i in 1:20){ #some steps in the beginning should be excluded
  
  if (i==1){
    y[[1]]<-ROI_vector
    z[[1]]<-ROI_vector
  }else{
    eps<-rnorm(n*n , mean=0 , sd=0.005)
    z[[i]]<-y[[1]]+I_minus_delt_W_inv%*%A%*%(y[[i-1]]-y[[1]])
    y[[i]]<-z[[i]]+I_minus_delt_W_inv%*%eps
  }
  
  
} 

#start of calculating statistic
t<-20
for (t in 21:120){
  
  eps<-rnorm(n*n , mean=0 , sd=0.005)
  z[[t]]<-z[[1]]+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-z[[1]])
  y[[t]]<-changed_vector+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-changed_vector)+I_minus_delt_W_inv%*%eps
  
  temp_sum<-0 #continue calculating statistics
  for (j in 0:(t-20)){
    temp_sum<-temp_sum+((1-lambda)^(2*j))
  }
  recalculated_E_O_t<-(n*n)*((lambda)*(temp_sum)*(lambda))
  recalculated_Var_O_t<- 2*(n*n)*((lambda*(temp_sum)*lambda)^2) 
  
  eta_t[[t-20]]<-Cov_of_residuals_correct_ROI%*%(z[[t]]-y[[t]])
  Z_t<-(I-Lambda)%*%Z_t_minus_1+Lambda%*%eta_t[[t-20]] #<---------------------------eta
  Z_t_minus_1<-Z_t
  O_t[t-20]<-t(Z_t)%*%Z_t
  CS_3t_res[t-20]<-(O_t[t-20]-recalculated_E_O_t)/(sqrt(recalculated_Var_O_t))
  
  
  
} 
plot(CS_3t_res, type = "b")


#CS_2_t_r -----------------------------------------------------------------------------------

#load("Sigma_l_r_by_formula_lambda_0_5.Rdata")
Sigma_l_r_by_formula_inv<-solve(Sigma_l_r_by_formula)

y<- list()
z<- list()

#start of calculating statistic
Z_t_minus_mu_minus_1<-rep(0,(n*n))

CS_2_t_r<-c() 
for (i in 1:20){ #some steps in the beginning should be excluded
  
  if (i==1){
    y[[1]]<-ROI_vector
    z[[1]]<-ROI_vector
  }else{
    eps<-rnorm(n*n , mean=0 , sd=0.005)
    z[[i]]<-y[[1]]+I_minus_delt_W_inv%*%A%*%(y[[i-1]]-y[[1]])
    y[[i]]<-z[[i]]+I_minus_delt_W_inv%*%eps
  }
  
  
} 


for (t in 21:120){ 
  
  eps<-rnorm(n*n , mean=0 , sd=0.005)
  z[[t]]<-z[[1]]+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-z[[1]])
  y[[t]]<-changed_vector+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-changed_vector)+I_minus_delt_W_inv%*%eps
  
  
  Z_t_minus_mu<-(I-Lambda)%*%Z_t_minus_mu_minus_1+Lambda%*%(y[[t]]-y[[1]])
  Z_t_minus_mu_minus_1<-Z_t_minus_mu
  
  sigmas_product<-Sigma_l_r_by_formula_inv%*%Sigmas[[t-20]]
  
  CS_2_t_r[t-20]<-((t(Z_t_minus_mu)%*%Sigma_l_r_by_formula_inv%*%Z_t_minus_mu)-tr(sigmas_product))/(sqrt(2*tr((sigmas_product)^2)))
  
  
  
}

plot(CS_2_t_r, type = "b")



#CS_3_t_r -----------------------------------------------------------------------------------

y<- list()
z<- list()

#start of calculating statistic
Z_t_minus_mu_minus_1<-rep(0,(n*n))

CS_3_t_r<-c() 
for (i in 1:20){ #some steps in the beginning should be excluded
  
  if (i==1){
    y[[1]]<-ROI_vector
    z[[1]]<-ROI_vector
  }else{
    eps<-rnorm(n*n , mean=0 , sd=0.005)
    z[[i]]<-y[[1]]+I_minus_delt_W_inv%*%A%*%(y[[i-1]]-y[[1]])
    y[[i]]<-z[[i]]+I_minus_delt_W_inv%*%eps
  }
  
} 

for (t in 21:120){ 
  
  eps<-rnorm(n*n , mean=0 , sd=0.005)
  z[[t]]<-z[[1]]+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-z[[1]])
  y[[t]]<-changed_vector+I_minus_delt_W_inv%*%A%*%(y[[t-1]]-changed_vector)+I_minus_delt_W_inv%*%eps
  
  
  Z_t_minus_mu<-(I-Lambda)%*%Z_t_minus_mu_minus_1+Lambda%*%(y[[t]]-y[[1]])
  Z_t_minus_mu_minus_1<-Z_t_minus_mu
  
  CS_3_t_r[t-20]<-(t(Z_t_minus_mu)%*%Z_t_minus_mu-tr(Sigmas[[t-20]]))/(sqrt(2*tr(Sigmas[[t-20]]%*%Sigmas[[t-20]])))
  
  
}

plot(CS_3_t_r, type = "b")










