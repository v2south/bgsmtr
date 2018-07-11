
Wang_algorithm_FUN = function(X, Y, group, r1, r2){
  
  
  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE))) 
  
  eps=2.2204e-16 
  
  d = dim(X)[1]
  p = dim(Y)[1]
  group_set=unique(group) 
  group_num=length(group_set) 
  
  
  ob=rep(NA,group_num)
  ob2=rep(NA,d)
  
  W_Wang=matrix( 1 , d, p) 
  
  XX=X%*%t(X)
  Xy=X%*%t(Y)
  
  d1=rep(1,d)
  d2=rep(1,d)
  Wi=matrix(0,d,1)  
  obj = c()  
  Tol = 10e-4  
  iter = 0  
  W_change = 1
  
  
  while ((W_change > Tol)) {    
    
    iter = iter + 1
    W_Wang_0 = W_Wang     
    D1=diag(d1)   
    D2=diag(d2)
    W_Wang=solve(XX + r1*D1 + r2*D2)%*%(Xy)   
    
    for (k in 1:group_num){              
      idx=which(group==group_set[k])
      idx
      W.k=W_Wang[idx,]
      di=sqrt(sum(W.k*W.k)+eps)
      Wi[idx]=di  
      ob[k]=di  
    }
    
    Wi=c(Wi)    
    d1=0.5/(Wi)    
    Wi2=sqrt(rowSums(W_Wang*W_Wang)+eps)    
    d2=0.5/Wi2    
    ob2=sum(Wi2)
    
    W_change = norm((W_Wang-W_Wang_0), type = "F")/max( norm(W_Wang_0, type = "F"), 1)    
    obj[iter]=sum(diag((t(t(X)%*%W_Wang-t(Y))%*%(t(X)%*%W_Wang-t(Y)))))+r1*sum(ob)+r2*ob2  
    if(iter > 100){break}
  }
  
  list_to_return = list("W_Wang" = W_Wang, "Wang_obj_func" = obj, 'num_iterations' = iter)
  
  return(list_to_return)
  
}



Wang_CV_tuning_values_FUN = function(X, Y, group){
  
  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE))) 
  
  
  d = nrow(X)
  
  
  p = nrow(Y)
  
  gamma_grid=expand.grid( g_1 = c(10e-4, 10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2, 10e3),  
                          g_2 = c(10e-4, 10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2, 10e3))
  
  CV_tuning_results = data.frame(gamma_grid)
  
  F = 5  
  
  R = length(CV_tuning_results$g_1)  
  
  CV_tuning_results$RMSE_fold_1 = rep(NA, R)
  CV_tuning_results$RMSE_fold_2 = rep(NA, R)
  CV_tuning_results$RMSE_fold_3 = rep(NA, R)
  CV_tuning_results$RMSE_fold_4 = rep(NA, R)
  CV_tuning_results$RMSE_fold_5 = rep(NA, R)
  CV_tuning_results$RMSE_mean = rep(NA, R)
  
  #data = data.frame(cbind(t(Y), t(X)),tuning_id = rep(NA,nrow(t(Y))))
  data = as.data.frame(cbind(t(Y), t(X)))    
  data$tuning_id = rep(NA, nrow(data))
  
  data$tuning_id = sample(1:F, nrow(data), replace = TRUE)
  cvlist = 1:F
  
  
  for (r in 1:R){  
    
    for (f in 1:F){ 
      
      training_set = data[data$tuning_id != f,]
      test_set = data[data$tuning_id ==f,]
      
      X_train = t(training_set[, (p+1): (p+d) ])  
      X_test = t(test_set[, (p+1): (p+d)] )    
      
      Y_train = t(training_set[, 1:p ])  
      Y_test = t(test_set[, 1:p] ) 
      
      
      W_Wang_hat = Wang_algorithm_FUN( X = X_train, Y = Y_train, group = group,
                                       r1 = CV_tuning_results$g_1[r],   
                                       r2 = CV_tuning_results$g_2[r] )$W_Wang
      
      Y_predict = t(W_Wang_hat)%*%X_test
      
      RMSE = sqrt(mean((Y_test-Y_predict)^2))
      
      
      CV_tuning_results[r, f+2] = RMSE   
      
    }
    
  }
  
  CV_tuning_results$RMSE_mean = rowMeans(CV_tuning_results[, 3:(2+F)])
  
  gamma_values = CV_tuning_results[which(CV_tuning_results$RMSE_mean==min(CV_tuning_results$RMSE_mean)), c(1,2,8)]
  
  
  penalty_1 = as.numeric(gamma_values[1])  
  penalty_2 = as.numeric(gamma_values[2])  
  
  Wang_algo_results = Wang_algorithm_FUN( X, Y, group, penalty_1, penalty_2 )
  
  
  list_to_return = list('CV_tuning_results' = CV_tuning_results, 
                        'selected_gamma_values' = gamma_values,
                        'W_Wang_from_tuning_CV' = Wang_algo_results$W_Wang)
  
  return(list_to_return)
  
}


bgsmtr.waic = function( X, Y, group, lam_1_fixed, lam_2_fixed, WAIC_opt = TRUE, 
                        iter_num = 10000, burn_in = 5001){
  
  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE))) 
  
  
  a_sig_prior=3
  b_sig_prior=1
  
  mean_prior_sig=b_sig_prior/(a_sig_prior-1)
  var_prior_sig=b_sig_prior^2/((a_sig_prior-1)^2*(a_sig_prior-2))
  
  Gibbs_setup_return = list( 'iter_num' = iter_num, 'burn_in' = burn_in,
                             'a_sig_prior' = a_sig_prior, 'b_sig_prior' = b_sig_prior, 
                             'lam_1_fixed' = lam_1_fixed, 'lam_2_fixed' = lam_2_fixed)
  
  
  d = dim(X)[1]  
  n = dim(X)[2]  
  
  p = dim(Y)[1]  
  
  
  group_set=unique(group) 
  K = length(group_set)
  m=rep(NA, K) 
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }
  
  idx=list()
  for (k in 1:K){
    idx[[k]]=which(group==group_set[k])
  }
  
  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p)
  }
  
  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}    
  }
  
  
  Xnk=list()
  for (k in 1:K){
    Xnk[[k]]=X[-idx[[k]],]
  }
  
  
  Xk_Xk=list()
  for (k in 1:K){
    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }
  
  
  
  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE) 
  }
  
  
  Xk_Y=list()
  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }
  
  
  
  p_waic = rep(0,n)
  
  log_p_waic = rep(0,n)
  
  log_p2_waic = rep(0,n)
  
  waic_iter = 0
  
  
  
  W_est=array(NA, c(iter_num,d,p))
  tau=array(NA, c(iter_num, K))
  omega=array(NA, c(iter_num, d))
  sig=c(rep(NA, iter_num))
  
  tau_init_value = 1   
  omega_init_value = 1     
  sig_init_value = 1
  
  
  tau[1,1:K] = rep(tau_init_value ,K)
  
  omega[1,1:d] = rep( omega_init_value,d)
  
  sig[1] = sig_init_value
  
  stm <- proc.time()
  cat('Computing Initial Values \n')
  
  W_Wang = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = 0.5, r2 = 0.5)$W_Wang
  end_time<-proc.time() - stm  
  cat('time: ', end_time[3], 's \n')
  cat('Gibbs Sampler Initialized and Running \n')
  
  W_est[1,1:d,1:p]=W_Wang
  
  W_int=W_Wang
  
  
  
  Wk=list()
  for (k in 1:K){
    Wk[[k]]=W_int[idx[[k]],] 
  }
  
  
  Wnk=list()
  for (k in 1:K){
    Wnk[[k]]=W_int[-idx[[k]],]
  }
  
  
  
  mkp=list()
  Ak=list()
  mu_k=list()
  vec.Wk_t=list()    
  Wk_t=list()
  
  
  
  stm <- proc.time()
  for (iter in 1:(iter_num-1)){
    
    
    for (k in 1:K){   
      
      Wnk[[k]]=W_int[-idx[[k]],]
      
      
      mkp[[k]]<-(1/tau[iter,k]) +(1/omega[iter,idx_mkp[[k]]])
      
      Ak[[k]]=Xk_Xk[[k]]+.symDiagonal(n=m[k]*p,x=mkp[[k]]) 
      
      CH<-Cholesky((1/sig[iter])*Ak[[k]])
      
      mu_k[[k]]<-solve(Ak[[k]],Xk_Y[[k]]-Xk_Xnk[[k]]%*%as.vector(t(Wnk[[k]])))
      
      
      vec.Wk_t[[k]]=rmvn.sparse(n=1, mu=mu_k[[k]], CH=CH,prec=TRUE)
      
      
      Wk_t[[k]]=matrix((vec.Wk_t[[k]]), p, m[k])   
      Wk[[k]]=t(Wk_t[[k]])  
      
      W_int[idx[[k]],]<-Wk[[k]]
      
    }
    
    
    W_est[(iter+1),1:d,1:p]=W_int
    
    
    for (k in 1:K){
      tau[(iter+1),k]=(rinvgauss(1, sqrt((lam_1_fixed*sig[iter])/(t(as.vector(Wk[[k]]))%*%as.vector(Wk[[k]]))), lam_1_fixed ))^-1
    }
    
    for (i in 1:d){
      omega[(iter+1),i]=(rinvgauss(1, sqrt((lam_2_fixed*sig[iter])/((t(W_est[(iter+1),i,1:p])%*%W_est[(iter+1),i,1:p]))), lam_2_fixed))^(-1)
    }
    
    
    Wij2_vk=0
    for (k in 1:K){
      for (i in 1:m[k]){
        Wij2_vk=(t(W_est[(iter+1),idx[[k]][i],1:p])%*%(W_est[(iter+1), idx[[k]][i], 1:p]))*
          ((1/tau[(iter+1),k])+ (1/omega[(iter+1),idx[[k]][i]])) +Wij2_vk        
      }
    }
    
    a_sig=(p*n)/2 + (d*p)/2 + a_sig_prior 
    b_sig=(norm((Y-t(W_est[(iter+1),1:d,1:p])%*%X), 'F'))^2/2 + Wij2_vk/2+ b_sig_prior
    
    sig[(iter+1)]=rinvgamma(1, a_sig, scale = b_sig )
    
    
    if (WAIC_opt){
      if(iter >= burn_in){  
        waic_iter = waic_iter + 1
        
        lik.vec<-dmnorm(t(Y),crossprod(X,W_est[(iter+1),1:d,1:p]),varcov=(sig[(iter+1)]*diag(p)))
        p_waic<-p_waic + lik.vec
        log_p_waic<-log_p_waic+log(lik.vec)
        log_p2_waic<-log_p2_waic+(log(lik.vec)^2)
      }
    }        
    
  }   
  end_time<-proc.time() - stm  
  cat('time: ', end_time[3], 's \n')
  
  
  if (WAIC_opt){
    approx_lpd = sum( log ( (1/waic_iter)*p_waic ) )   
    
    approx_P_waic = sum(  (1/(waic_iter -1))*log_p2_waic  -
                            (1/(waic_iter*(waic_iter -1)))*(log_p_waic^2) )
    
    WAIC = -2*(approx_lpd - approx_P_waic)
  }
  
  
  mcmc_tau = mcmc(tau[burn_in:iter_num,1:K])
  mcmc_omega = mcmc(omega[burn_in:iter_num,1:d])
  mcmc_sig = mcmc(sig[burn_in:iter_num])
  mcmc_W_est = list()
  
  for(q in 1:d) {
    mcmc_W_est[[q]] = mcmc(as.matrix(W_est[burn_in:iter_num,q, 1:p]))  
  }
  
  
  mcmc_W_est_summaries=list()
  for (q in 1:d){
    mcmc_W_est_summaries[[q]]=summary(mcmc_W_est[[q]]) 
  }
  
  mcmc_tau_summary=summary(mcmc_tau)
  mcmc_omega_summary=summary(mcmc_omega)
  mcmc_sig_summary=summary(mcmc_sig)
  
  
  W_post_mean=matrix(NA, d, p)
  W_post_sd=matrix(NA, d, p)
  
  for (q in 1:d){
    W_post_mean[q,]=mcmc_W_est_summaries[[q]][[1]][,1]
    W_post_sd[q,]=mcmc_W_est_summaries[[q]][[1]][,2]
  }
  
  
  W_2.5_quantile=matrix(NA, d, p)
  W_97.5_quantile=matrix(NA, d, p)
  
  for (q in 1:d){
    W_2.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,1]
    W_97.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,5]
  }
  
  stm <- proc.time()
  cat('Computing Approximate Posterior Mode \n')
  r1.final<-2*mean(sqrt(sig[burn_in:iter_num]))*sqrt(lam_1_fixed)
  r2.final<-2*mean(sqrt(sig[burn_in:iter_num]))*sqrt(lam_2_fixed)
  W_post_mode = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = r1.final, r2 = r2.final)$W_Wang
  end_time<-proc.time() - stm  
  cat('time: ', end_time[3], 's \n')
  
  row.names(W_post_mode) = row.names(W_post_mean) = row.names(W_post_sd) = row.names(W_2.5_quantile) = row.names(W_97.5_quantile) = row.names(X)
  
  colnames(W_post_mode) = colnames(W_post_mean) = colnames(W_post_sd) = colnames(W_2.5_quantile) = colnames(W_97.5_quantile) = row.names(Y)
  
  
  
  Gibbs_W_summaries_return = list(  'W_post_mean' =  W_post_mean,
                                    'W_post_mode' = W_post_mode,
                                    'W_post_sd' = W_post_sd, 
                                    'W_2.5_quantile' =  W_2.5_quantile, 
                                    'W_97.5_quantile' =  W_97.5_quantile)
  
  
  if (WAIC_opt){ function_returns = list( 'WAIC' = WAIC,  'Gibbs_setup' =   Gibbs_setup_return,
                                          'Gibbs_W_summaries' = Gibbs_W_summaries_return)   } else { function_returns =
                                            list( 'Gibbs_setup' =   Gibbs_setup_return,
                                                  'Gibbs_W_summaries' = Gibbs_W_summaries_return) }
  
  
  
  return(function_returns)
  
}




bgsmtr.cv.mode = function( X, Y, group, iter_num = 10000, burn_in = 5001){
  
  
  WAIC_opt = FALSE
  
  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE))) 
  
  
  a_sig_prior=3
  b_sig_prior=1
  
  mean_prior_sig=b_sig_prior/(a_sig_prior-1)
  var_prior_sig=b_sig_prior^2/((a_sig_prior-1)^2*(a_sig_prior-2))
  
  
  Gibbs_setup_return = list( 'iter_num' = iter_num, 'burn_in' = burn_in,
                             'a_sig_prior' = a_sig_prior, 'b_sig_prior' = b_sig_prior, 
                             'lam_1_fixed' = NULL, 'lam_2_fixed' = NULL)
  
  
  
  d = dim(X)[1]  
  n = dim(X)[2]  
  
  p = dim(Y)[1]  
  
  
  group_set=unique(group) 
  K = length(group_set)
  m=rep(NA, K) 
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }
  
  
  idx=list()
  for (k in 1:K){
    idx[[k]]=which(group==group_set[k])
  }
  
  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p)
  }
  
  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}   
  }
  
  
  Xnk=list()
  for (k in 1:K){
    Xnk[[k]]=X[-idx[[k]],]
  }
  
  
  Xk_Xk=list()
  for (k in 1:K){
    
    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }
  
  
  
  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE) 
  }
  
  
  Xk_Y=list()
  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }
  
  
  
  p_waic = rep(0,n)
  
  log_p_waic = rep(0,n)
  
  log_p2_waic = rep(0,n)
  
  waic_iter = 0
  
  
  W_est=array(NA, c(iter_num,d,p))
  tau=array(NA, c(iter_num, K))
  omega=array(NA, c(iter_num, d))
  sig=c(rep(NA, iter_num))
  
  tau_init_value = 1   
  omega_init_value = 1     
  sig_init_value = 1
  
  
  tau[1,1:K] = rep(tau_init_value ,K)
  
  omega[1,1:d] = rep( omega_init_value,d)
  
  sig[1] = sig_init_value
  
  stm <- proc.time()
  cat('Computing Initial Values and Estimating Tuning Parameters Using Five-Fold Cross-Validation \n')
  
  CV_results = Wang_CV_tuning_values_FUN(X = X , Y = Y, group = group)
  W_Wang = CV_results$W_Wang_from_tuning_CV
  W_post_mode = W_Wang
  r1.hat = CV_results$selected_gamma_values[1]
  r2.hat = CV_results$selected_gamma_values[2]
  end_time<-proc.time() - stm  
  cat('time: ', end_time[3], 's \n')
  cat('Gibbs Sampler Initialized and Running \n')
  
  W_est[1,1:d,1:p]=W_Wang
  
  W_int=W_Wang
  
  
  lam_1_fixed = (r1.hat/(2*sqrt(sig[1])))^2
  lam_2_fixed = (r2.hat/(2*sqrt(sig[1])))^2
  
  
  
  Wk=list()
  for (k in 1:K){
    Wk[[k]]=W_int[idx[[k]],] 
  }
  
  
  Wnk=list()
  for (k in 1:K){
    Wnk[[k]]=W_int[-idx[[k]],]
  }
  
  
  
  mkp=list()
  Ak=list()
  mu_k=list()
  vec.Wk_t=list()    
  Wk_t=list()
  
  
  stm <- proc.time()
  for (iter in 1:(iter_num-1)){
    
    for (k in 1:K){   
      
      Wnk[[k]]=W_int[-idx[[k]],]
      
      mkp[[k]]<-(1/tau[iter,k]) +(1/omega[iter,idx_mkp[[k]]])
      
      Ak[[k]]=Xk_Xk[[k]]+.symDiagonal(n=m[k]*p,x=mkp[[k]]) 
      
      CH<-Cholesky((1/sig[iter])*Ak[[k]])
      
      
      mu_k[[k]]<-solve(Ak[[k]],Xk_Y[[k]]-Xk_Xnk[[k]]%*%as.vector(t(Wnk[[k]])))
      
      
      vec.Wk_t[[k]]=rmvn.sparse(n=1, mu=mu_k[[k]], CH=CH,prec=TRUE)
      
      
      Wk_t[[k]]=matrix((vec.Wk_t[[k]]), p, m[k])   
      Wk[[k]]=t(Wk_t[[k]])  
      
      W_int[idx[[k]],]<-Wk[[k]]
      
    }
    
    
    W_est[(iter+1),1:d,1:p]=W_int
    
    for (k in 1:K){
      tau[(iter+1),k]=(rinvgauss(1, sqrt((lam_1_fixed*sig[iter])/(t(as.vector(Wk[[k]]))%*%as.vector(Wk[[k]]))), lam_1_fixed ))^-1
    }
    
    for (i in 1:d){
      omega[(iter+1),i]=(rinvgauss(1, sqrt((lam_2_fixed*sig[iter])/((t(W_est[(iter+1),i,1:p])%*%W_est[(iter+1),i,1:p]))), lam_2_fixed))^(-1)
    }
    
    Wij2_vk=0
    for (k in 1:K){
      for (i in 1:m[k]){
        Wij2_vk=(t(W_est[(iter+1),idx[[k]][i],1:p])%*%(W_est[(iter+1), idx[[k]][i], 1:p]))*
          ((1/tau[(iter+1),k])+ (1/omega[(iter+1),idx[[k]][i]])) +Wij2_vk        
      }
    }
    
    a_sig=(p*n)/2 + (d*p)/2 + a_sig_prior 
    b_sig=(norm((Y-t(W_est[(iter+1),1:d,1:p])%*%X), 'F'))^2/2 + Wij2_vk/2+ b_sig_prior
    
    sig[(iter+1)]=rinvgamma(1, a_sig, scale = b_sig )
    
    
    lam_1_fixed = (r1.hat/(2*sqrt(sig[(iter+1)])))^2
    lam_2_fixed = (r2.hat/(2*sqrt(sig[(iter+1)])))^2
    
    
    
    if (WAIC_opt){
      if(iter >= burn_in){  # start calculations after burn_in
        waic_iter = waic_iter + 1
        
        lik.vec<-dmnorm(t(Y),crossprod(X,W_est[(iter+1),1:d,1:p]),varcov=(sig[(iter+1)]*diag(p)))
        p_waic<-p_waic + lik.vec
        log_p_waic<-log_p_waic+log(lik.vec)
        log_p2_waic<-log_p2_waic+(log(lik.vec)^2)
      }
    }   
    
    
  }   
  end_time<-proc.time() - stm  
  cat('time: ', end_time[3], 's \n')
  
  
  if (WAIC_opt){
    approx_lpd = sum( log ( (1/waic_iter)*p_waic ) )   
    
    approx_P_waic = sum(  (1/(waic_iter -1))*log_p2_waic  -
                            (1/(waic_iter*(waic_iter -1)))*(log_p_waic^2) )
    
    WAIC = -2*(approx_lpd - approx_P_waic)
  }
  
  
  
  mcmc_tau = mcmc(tau[burn_in:iter_num,1:K])
  mcmc_omega = mcmc(omega[burn_in:iter_num,1:d])
  mcmc_sig = mcmc(sig[burn_in:iter_num])
  mcmc_W_est = list()
  
  for(q in 1:d) {
    mcmc_W_est[[q]] = mcmc(as.matrix(W_est[burn_in:iter_num,q, 1:p]))  
  }
  
  
  mcmc_W_est_summaries=list()
  for (q in 1:d){
    mcmc_W_est_summaries[[q]]=summary(mcmc_W_est[[q]]) 
  }
  
  mcmc_tau_summary=summary(mcmc_tau)
  mcmc_omega_summary=summary(mcmc_omega)
  mcmc_sig_summary=summary(mcmc_sig)
  
  
  
  W_post_mean=matrix(NA, d, p)
  W_post_sd=matrix(NA, d, p)
  
  for (q in 1:d){
    W_post_mean[q,]=mcmc_W_est_summaries[[q]][[1]][,1]
    W_post_sd[q,]=mcmc_W_est_summaries[[q]][[1]][,2]
  }
  
  
  
  W_2.5_quantile=matrix(NA, d, p)
  W_97.5_quantile=matrix(NA, d, p)
  
  for (q in 1:d){
    W_2.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,1]
    W_97.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,5]
  }
  
  
  row.names(W_post_mode) = row.names(W_post_mean) = row.names(W_post_sd) = row.names(W_2.5_quantile) = row.names(W_97.5_quantile) = row.names(X)
  
  colnames(W_post_mode) = colnames(W_post_mean) = colnames(W_post_sd) = colnames(W_2.5_quantile) = colnames(W_97.5_quantile) = row.names(Y)
  
  
  Gibbs_W_summaries_return = list(  'W_post_mean' =  W_post_mean,
                                    'W_post_mode' = W_post_mode,
                                    'W_post_sd' = W_post_sd, 
                                    'W_2.5_quantile' =  W_2.5_quantile, 
                                    'W_97.5_quantile' =  W_97.5_quantile)
  
  
  if (WAIC_opt){ function_returns = list( 'WAIC' = WAIC,  'Gibbs_setup' =   Gibbs_setup_return,
                                          'Gibbs_W_summaries' = Gibbs_W_summaries_return)   } else { function_returns =
                                            list( 'Gibbs_setup' =   Gibbs_setup_return,
                                                  'Gibbs_W_summaries' = Gibbs_W_summaries_return) }
  
  
  
  return(function_returns)
  
}

#' Bayesian Group Sparse Multi-Task Regression for Imaging Genetics
#'
#' Runs the the Gibbs sampling algorithm to fit a Bayesian group sparse multi-task regression model. 
#' Tuning parameters can be chosen using either the MCMC samples and the WAIC (multiple runs) or using an approximation to 
#' the posterior mode and five-fold cross-validation (single run).
#'  
#'
#' @param X A d-by-n matrix; d is the number of SNPs and n is the number of subjects. Each row of X should correspond to a particular SNP
#' and each column should correspond to a particular subject. Each element of X should give the number of minor alleles for the corresponding
#' SNP and subject. The function will center each row of X to have mean zero prior to running the Gibbs sampling algorithm.
#' @param Y A c-by-n matrix; c is the number of phenotypes (brain imaging measures) and n is the number of subjects. Each row of 
#' Y should correspond to a particular phenotype and each column should correspond to a particular subject. Each element of Y should give 
#' the measured value for the corresponding phentoype and subject. The function will center and scale each row of Y to have mean zero and unit 
#' variance prior to running the Gibbs sampling algorithm.
#' @param group A vector of length d; d is the number of SNPs. Each element of this vector is a string representing a gene or group
#' label associated with each SNP. The SNPs represented by this vector should be ordered according to the rows of X.
#' @param tuning A string, either 'WAIC' or 'CV.mode'. If 'WAIC', the Gibbs sampler is run with fixed values of the tuning 
#' parameters specified by the arguments \emph{lam_1_fixed} and  \emph{lam_2_fixed} and the WAIC is computed based on the sampling output. This
#' can then be used to choose optimal values for \emph{lam_1_fixed} and \emph{lam_2_fixed} based on multiple runs with each run using different
#' values of \emph{lam_1_fixed} and \emph{lam_2_fixed}. This option is best suited for either comparing a small set of tuning parameter values or
#' for computation on a high performance computing cluster where different nodes can be used to run the function with different 
#' values of \emph{lam_1_fixed} and \emph{lam_2_fixed}. Posterior inference is then based on the run that produces the lowest value for the WAIC. 
#' The option 'CV.mode', which is the default, is best suited for computation using just a single processor. In this case the 
#' tuning parameters are chosen based on five-fold cross-validation over a grid of possible values with out-of-sample prediction based on an 
#' approximate posterior mode. The Gibbs sampler is then run using the chosen values of the tuning parameters. When tuning = 'CV.mode' the values 
#' for the arguments \emph{lam_1_fixed} and \emph{lam_2_fixed} are not required.
#' @param lam_1_fixed Only required if tuning = 'WAIC'. A positive number giving the value for the gene-specific tuning parameter. Larger values lead to a larger
#' degree of shrinkage to zero of estimated regression coefficients at the gene level (across all SNPs and phenotypes). 
#' @param lam_2_fixed Only required if tuning = 'WAIC'. A positive number giving the value for the SNP-specific tuning parameter. Larger values lead to a larger
#' degree of shrinkage to zero of estimated regression coefficients at the SNP level (across all phenotypes).
#' @param iter_num Positive integer representing the total number of iterations to run the Gibbs sampler. Defaults to 10,000.
#' @param burn_in Nonnegative integer representing the number of MCMC samples to discard as burn-in. Defaults to 5001.
#' 
#' @return A list with the elements
#' \item{WAIC}{If tuning = 'WAIC' this is the value of the WAIC computed from the MCMC output. If tuning = 'CV.mode' this component is excluded.}
#' \item{Gibbs_setup}{A list providing values for the input parameters of the function.}
#' \item{Gibbs_W_summaries}{A list with five components, each component being a d-by-c matrix giving some posterior summary of the regression parameter
#' matrix W, where the ij-th element of W represents the association between the i-th SNP and j-th phenotype. 
#' 
#' -Gibbs_W_summaries$W_post_mean is a d-by-c matrix giving the posterior mean of W. 
#' 
#' 
#' 
#' -Gibbs_W_summaries$W_post_mode is a d-by-c matrix giving the posterior mode of W. 
#' 
#' 
#' 
#' -Gibbs_W_summaries$W_post_sd is a d-by-c matrix giving the posterior standard deviation for each element of W.
#' 
#' 
#' 
#' -Gibbs_W_summaries$W_2.5_quantile is a d-by-c matrix giving the posterior 2.5 percent quantile for each element of W. 
#' 
#' 
#' 
#' -Gibbs_W_summaries$W_97.5_quantile is a d-by-c matrix giving the posterior 97.5 percent quantile for each element of W.'}
#' 
#' 
#' @author Farouk S. Nathoo, \email{nathoo@uvic.ca} 
#' @author Keelin Greenlaw  \email{keelingreenlaw@gmail.com}
#' @author Mary Lesperance  \email{mlespera@uvic.ca}
#' 
#' @examples
#' data(bgsmtr_example_data)
#' names(bgsmtr_example_data)
#' \dontshow{
#' ## Toy example with a small subset of the data for routine CRAN testing
#' ## reduce the number of SNPs to 5, subjects to 5, phenotypes to 5 
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data[1:5,1:5], Y = bgsmtr_example_data$BrainMeasures[1:5,1:5], 
#' group = bgsmtr_example_data$SNP_groups[1:5], tuning = 'WAIC', lam_1_fixed = 2, lam_2_fixed = 2,
#' iter_num = 5, burn_in = 1)
#' }
#' 
#' \dontrun{
#' ## test run the sampler for 100 iterations with fixed tunning parameters and compute WAIC
#' ## we recomend at least 5,000 iterations for actual use
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures, 
#' group = bgsmtr_example_data$SNP_groups, tuning = 'WAIC', lam_1_fixed = 2, lam_2_fixed = 2,
#' iter_num = 100, burn_in = 50)
#' ## posterior mean for regression parameter relating 100th SNP to 14th phenotype  
#' fit$Gibbs_W_summaries$W_post_mean[100,14]   
#' ## posterior mode for regression parameter relating 100th SNP to 14th phenotype
#' fit$Gibbs_W_summaries$W_post_mode[100,14]
#' ## posterior standard deviation for regression parameter relating 100th SNP to 14th phenotype
#' fit$Gibbs_W_summaries$W_post_sd[100,14]
#' ## 95% equal-tail credible interval for regression parameter relating 100th SNP to 14th phenotype
#' c(fit$Gibbs_W_summaries$W_2.5_quantile[100,14],fit$Gibbs_W_summaries$W_97.5_quantile[100,14])     
#'}
#'     
#'\dontrun{
#' ## run the sampler for 10,000 iterations with tuning parameters set using cross-validation
#' ## On a standard computer with a small numer of cores this is the recomended option
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures, 
#' group = bgsmtr_example_data$SNP_groups, tuning = 'CV.mode',iter_num = 10000, burn_in = 5000)     
#'}     
#'
#' @references Greenlaw, Keelin, Elena Szefer, Jinko Graham, Mary Lesperance, and Farouk S. Nathoo. "A Bayesian Group Sparse Multi-Task Regression Model for Imaging Genetics." arXiv preprint arXiv:1605.02234 (2016).   
#' @references Nathoo, Farouk S., Keelin Greenlaw, and Mary Lesperance. "Regularization Parameter Selection for a Bayesian Multi-Level Group Lasso Regression Model with Application to Imaging Genomics." arXiv preprint arXiv:1603.08163 (2016).    
#' 
#' @import Matrix mvtnorm
#' @importFrom sparseMVN rmvn.sparse
#' @importFrom statmod rinvgauss
#' @importFrom EDISON rinvgamma
#' @importFrom coda mcmc
#' @importFrom mnormt dmnorm  
#'          
#' @export
bgsmtr = function(X, Y, group, tuning = 'CV.mode', lam_1_fixed = NULL, lam_2_fixed = NULL, iter_num = 10000, burn_in = 5001)
{
  if (tuning == 'WAIC')
  {
    result = bgsmtr.waic( X=X, Y=Y, group=group, lam_1_fixed=lam_1_fixed, lam_2_fixed=lam_2_fixed, WAIC_opt = TRUE, 
                          iter_num = iter_num, burn_in = burn_in)
  }
  else 
  {
    result = bgsmtr.cv.mode( X=X, Y=Y, group=group, iter_num = iter_num, burn_in = burn_in)
  }
  return(result)
}

  
  




