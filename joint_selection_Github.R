
# install.packages("lars")

# Y : n X 1 data 
# X : n X p covariate matrix
DAGWishart <- function(Y, X, qn, alpha.plus, U, init.A=NULL, init.gam=NULL, Rj=NULL, a.hp, b.hp, tau2, sigma2, a0, b0, c0, niter, nburn, nadap=0){
   
   n = nrow(X)
   p = ncol(X)
   if(!exists("lars")) library(lars)
   
   # set Rj value
   if(is.null(Rj)) Rj = floor( n/(log(n, base=10)) )
   dj = NULL
   
   res = list()
   
   Sj = list()
   Sj.init = list()
   Sj.mat = list()
   sj.mat = list()
   Rsq.mat = list()
   
   logpost.mat = list()
   logpost.new = list()
   logpost.old = list()
   
   tilde.Xj = list()
   tilde.Xj.prod = list()
   Z.j = list()
   accept.ind = list(); accept.ind[[1]] = 0
   unique.index = list()
   logpost.uniq.mat = list()
   dj.uniq.mat = list()
   
   # set initial value for DAG
   for(j in 2:p){
     
     Sj.mat[[j]] = matrix(0, nrow=niter + nburn, ncol=j-1)
     sj.mat[[j]] = rep(0, niter + nburn)
     Rsq.mat[[j]] = rep(0, niter + nburn)
     logpost.mat[[j]] = rep(0, niter + nburn)
     
     # data (tilde.Xj) and design matrix (Z.j) for j-th dimension
     tilde.Xj[[j]] = as.matrix(X[, j])
     Z.j[[j]] = as.matrix(X[, 1:(j-1)])
     tilde.Xj.prod[[j]] = sum( tilde.Xj[[j]]^2 )
     
     # initialzation
     if(is.null(init.A)){
       # initial guess for Sj based on lasso
       ob.j = lars(x = Z.j[[j]], y = tilde.Xj[[j]], normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
       cv.j = cv.lars(x = Z.j[[j]], y = tilde.Xj[[j]], plot.it = FALSE, se = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
       las.est = coef(ob.j, s = cv.j$index[which.min(cv.j$cv)], mode = "fraction")
       if(sum(las.est != 0) > min(j-1, Rj)){
         las.est[ sample(which(las.est != 0), size = sum(las.est != 0) - min(j-1, Rj), replace = FALSE) ] = 0
       }
       Sj.init[[j]] = as.numeric(las.est != 0)
       cat("Initial value for S",j," is set. . . . . .\n")
     }else{
       Sj.init[[j]] = as.numeric(init.A[j, 1:(j-1)] != 0)
     }
     
     if(sum(Sj.init[[j]]) == 0) Sj.init[[j]][sample(1:(j-1), size = 1)] = 1
   }
   
   if(is.null(init.gam)){
     # set initial value for gamma
     ob.y = lars(x = X, y = Y, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
     cv.y = cv.lars(x = X, y = Y, plot.it = FALSE, se = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
     las.est = coef(ob.y, s = cv.y$index[which.min(cv.y$cv)], mode = "fraction")
     if(sum(las.est != 0) > min(j-1, Rj)){
       las.est[ sample(which(las.est != 0), size = sum(las.est != 0) - min(j-1, Rj), replace = FALSE) ] = 0
     }
     gamma.init = as.numeric(las.est != 0)
     cat("Initial value for gamma is set. . . . . .\n")
   }else{
     gamma.init = init.gam
   }
   
   gamma.mat = matrix(0, nrow=niter + nburn, ncol=p)
   gamma.mat[1,] = gamma.init
   
   
   
   for(j in 2:p){
      
      Sj[[j]] = Sj.init[[j]]
      logpost.old[[j]] = logpost.cao(X, j, Sj[[j]], qn, alpha.plus, U, gamma.init, b.hp)
      
      if(is.finite(logpost.old[[j]]$val) == FALSE){
         # initial guess for Sj based on lasso
         ob.j = lars(x = Z.j[[j]], y = tilde.Xj[[j]], normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         cv.j = cv.lars(x = Z.j[[j]], y = tilde.Xj[[j]], plot.it = FALSE, se = FALSE, normalize = FALSE, intercept = FALSE, use.Gram = FALSE)
         las.est = coef(ob.j, s = cv.j$index[which.min(cv.j$cv)], mode = "fraction")
         if(sum(las.est != 0) > min(j-1, Rj)){
            las.est[ sample(which(las.est != 0), size = sum(las.est != 0) - min(j-1, Rj), replace = FALSE) ] = 0
         }
         Sj.init[[j]] = as.numeric(las.est != 0)
         Sj[[j]] = Sj.init[[j]]
         logpost.old[[j]] = logpost.cao(X, j, Sj[[j]], qn, alpha.plus, U, gamma.init, b.hp)
      }
      
      Sj.mat[[j]][1, ] = Sj[[j]]
      sj.mat[[j]][1] = sum(Sj[[j]])
      Rsq.mat[[j]][1] = 1 - n*logpost.old[[j]]$dj / tilde.Xj.prod[[j]]
      logpost.mat[[j]][1] = logpost.old[[j]]$val
      logpost.uniq.mat[[j]] = logpost.old[[j]]$val
      dj.uniq.mat[[j]] = logpost.old[[j]]$dj
      unique.index[[j]] = 1
      logpost.new[[j]] = list()
      
      accept.ind[[j]] = 0 # number of acceptance in MH sampler
      
   }
   
   
   
   # MCMC sampling
   for(i in 2:(niter+nburn)){
      
      
      # MH part for DAG
      for(j in 2:p){
         
         if(i <= nadap){
            Sj.new = Sprop.adap(X, Sj[[j]], j, Rj) # new proposal of Sj (Sj.new)
            log.prop.kernel = lkernelprob.adap.cao(X, Sj[[j]], Sj.new, j) # log proportion of proposal kernel
         }else{
            Sj.new = Sprop(Sj[[j]], j, Rj) # new proposal of Sj (Sj.new)
            log.prop.kernel = lkernelprob.cao(Sj[[j]], Sj.new, j) # log proportion of proposal kernel
         }
         
         # calcultaion of logpost.new
         if(is.na(unique.index[[j]][2])){
            match.ind = prod(as.numeric(Sj.mat[[j]][1,] == Sj.new))
         }else{
            match.ind = compare.to.rows(as.matrix(Sj.mat[[j]][unique.index[[j]] ,]), Sj.new)
         }
         
         if(match.ind == 0){
            logpost.new[[j]] = logpost.cao(X, j, Sj.new, qn, alpha.plus, U, gamma.mat[i-1,], b.hp)
         }else{
            logpost.new[[j]]$val = logpost.uniq.mat[[j]][match.ind]
            logpost.new[[j]]$dj = dj.uniq.mat[[j]][match.ind]
         }
         
         
         # MH sampler for Sj index
         if(is.finite(logpost.new[[j]]$val) == FALSE) logpost.new[[j]]$val = -Inf
         if(runif(1) <= exp(logpost.new[[j]]$val - logpost.old[[j]]$val + log.prop.kernel)){
            Sj[[j]] = Sj.new
            logpost.old[[j]] = logpost.new[[j]]
            if(match.ind == 0){
               logpost.uniq.mat[[j]] = c(logpost.uniq.mat[[j]], logpost.old[[j]]$val)
               dj.uniq.mat[[j]] = c(dj.uniq.mat[[j]], logpost.old[[j]]$dj)
               unique.index[[j]] = c(unique.index[[j]], i)
            }
            accept.ind[[j]] = accept.ind[[j]] + 1
         }
         
         
         Sj.mat[[j]][i, ] = Sj[[j]]
         sj.mat[[j]][i] = sum(Sj[[j]])
         Rsq.mat[[j]][i] = 1 - n*logpost.old[[j]]$dj / tilde.Xj.prod[[j]]
         logpost.mat[[j]][i] = logpost.old[[j]]$val
      
      } # end of (j in 2:p) for loop 
      
      
      
      # MH part for gamma
      gamma.old = gamma.mat[i-1,]
      gamma.new = Gammaprop(gamma.old, Rj) # new proposal of gamma (gamma.new)
      
      new.square = 0
      old.square = 0
      adj.S.mat = matrix(0, p,p) # adjacent matrix for DAG
      for(jj in 2:p){
        adj.S.mat[jj, 1:(jj-1)] = Sj.mat[[jj]][i, ]
      }
      adj.S.mat = adj.S.mat + t(adj.S.mat)
      new.square = t(matrix(gamma.new)) %*% adj.S.mat %*% matrix(gamma.new)
      old.square = t(matrix(gamma.old)) %*% adj.S.mat %*% matrix(gamma.old)
      
      X.gam.new = X[, gamma.new*(1:p)]
      X.gam.old = X[, gamma.old*(1:p)]
      
      # case (a): sigma2 is known 
      # We set beta_gamma | gamma ~ N(0, tau2*I) for some tau2 > 0
      # logpost.gam.new = - a.hp*sum(gamma.new) + b.hp*new.square - t(Y)%*%solve( diag(n) + tau2*X.gam.new%*%t(X.gam.new) )%*%Y/(2*sigma2) 
      # logpost.gam.old = - a.hp*sum(gamma.old) + b.hp*old.square - t(Y)%*%solve( diag(n) + tau2*X.gam.old%*%t(X.gam.old) )%*%Y/(2*sigma2)
      
      # case (b): sigma2 is unknown
      # We set sigma2 ~ IG(a0, b0)
      # and beta_gamma | gamma, sigma2 ~ N(0, c0*sigma2*I) for some c0 > 0
      logpost.gam.new = - a.hp*sum(gamma.new) + b.hp*new.square - 0.5*(n+2*a0)*log( b0 + 0.5*t(Y)%*%solve( diag(n) + c0*X.gam.new%*%t(X.gam.new) )%*%Y ) - 0.5*as.numeric(determinant(diag(n) + c0*X.gam.new%*%t(X.gam.new), logarithm = T)$mod)
      logpost.gam.old = - a.hp*sum(gamma.old) + b.hp*old.square - 0.5*(n+2*a0)*log( b0 + 0.5*t(Y)%*%solve( diag(n) + c0*X.gam.old%*%t(X.gam.old) )%*%Y ) - 0.5*as.numeric(determinant(diag(n) + c0*X.gam.old%*%t(X.gam.old), logarithm = T)$mod)
      
      if(runif(1) <= exp(logpost.gam.new - logpost.gam.old)){
        gamma.mat[i,] = gamma.new
        gamma.old = gamma.new
        accept.ind[[1]] = accept.ind[[1]] + 1
      }else{
        gamma.mat[i,] = gamma.old
      }
      
      
      
      cat(i,"th iteration is completed. . . . . .\n")
      
   } # end of (i in 2:(niter+nburn)) for loop
   
   
   for(j in 2:p){
     Sj.mat[[j]] = Sj.mat[[j]][-(1:nburn),]
     sj.mat[[j]] = sj.mat[[j]][-(1:nburn)]
     Rsq.mat[[j]] = Rsq.mat[[j]][-(1:nburn)]
     logpost.mat[[j]] = logpost.mat[[j]][-(1:nburn)]
   } 
   
   res = list(Sj.mat = Sj.mat, sj.mat = sj.mat, gamma.mat = gamma.mat[-(1:nburn),],
              accept.ind = accept.ind,  Rsq.mat = Rsq.mat, logpost.mat = logpost.mat)
   
   return(res)
}


###################################################################
# Auxiliary functions
###################################################################

compare.to.rows <- function(SS, S){
  h <- function(v) sum(abs(v - S))
  o <- apply(SS, 1, h)
  if(all(o > 0)){
    return(0)
  }else{
    return(which(o == 0)[1])
  }
}

Sprop <- function(S, j, Rj, prop.type=2){
   s = sum(S)
   upper.ind = min(j-1, Rj)
   
   if(s == 0){ # if current S has no index
      S[sample(which(S == 0), 1)] = 1
   }else if(s == upper.ind){ # if current S has maximum index (upper.ind)
      S[sample(which(S > 0), 1)] = 0
   }else{
      
      if(prop.type==1){ # (1) sample any index with the same prob.
        ind = sample(x = 1:length(S), size = 1)
        S[ind] = abs(S[ind] - 1)
      }else{ # (2) sample 1's w.p. 0.5 or sample 0's w.p. 0.5
        
        if(runif(1) <= 0.5){ # introducing additional one 0 to current S
          if(s == 1){
            S[which(S == 1)] = 0
          }else{
            S[sample(which(S > 0), 1)] = 0
          }
          
        }else{ # introducing additional one 1 to current S
          if(j-1-s == 1){
            S[which(S == 0)] = 1
          }else{
            S[sample(which(S == 0), 1)] = 1
          }
        }
        
      }
   }
   return(S)
}
Sprop.adap <- function(X, S, j, Rj){
   n = nrow(X)
   s = sum(S)
   upper.ind = min(j-1, Rj)
   
   if(s == 0){ # if current S has no index
      S[sample(which(S == 0), 1)] = 1
   }else if(s == upper.ind){ # if current S has maximum index (upper.ind)
      S[sample(which(S > 0), 1)] = 0
   }else{
      
      if(runif(1) <= 0.5){ # introducing additional one 0 to current S
         if(s == 1){
            S[which(S == 1)] = 0
         }else{
            prob.of.new.ind = rep(0, s)
            for(k in 1:s){
               S.old = S
               S.old[which(S == 1)[k]] = 0
               prob.of.new.ind[k] = 1/dShat(X, j, S.old)
            }
            S[sample(which(S == 1), 1, prob = prob.of.new.ind)] = 0
         }
      }else{ # introducing additional one 1 to current S
         if(j-1-s == 1){
            S[which(S == 0)] = 1
         }else{
            prob.of.new.ind = rep(0, j-1-s)
            for(k in 1:(j-1-s)){
               S.old = S
               S.old[which(S == 0)[k]] = 1
               prob.of.new.ind[k] = 1/dShat(X, j, S.old)
            }
            S[sample(which(S == 0), 1, prob = prob.of.new.ind)] = 1
         }
      }
      
   }
   return(S)
}

Gammaprop <- function(S, Rj, prop.type=1){
   s = sum(S)
   upper.ind = min(p, Rj)
   
   if(s == 0){ # if current S has no index
      S[sample(which(S == 0), 1)] = 1
   }else if(s == upper.ind){ # if current S has maximum index (upper.ind)
      S[sample(which(S > 0), 1)] = 0
   }else{
      
     if(prop.type==1){ # (1) sample any index with the same prob.
       ind = sample(x = 1:length(S), size = 1)
       S[ind] = abs(S[ind] - 1)
     }else{ # (2) sample 1's w.p. 0.5 or sample 0's w.p. 0.5
       if(runif(1) <= 0.5){ # introducing additional one 0 to current S
         if(s == 1){
           S[which(S == 1)] = 0
         }else{
           S[sample(which(S > 0), 1)] = 0
         }
       }else{ # introducing additional one 1 to current S
         if(p-s == 1){
           S[which(S == 0)] = 1
         }else{
           S[sample(which(S == 0), 1)] = 1
         }
       }
     }
      
   }
   return(S)
}


lkernelprob.cao <- function(Sj, Sj.new, j){
   if(sum(Sj - Sj.new) == -1){
      return(log(j-1 - sum(Sj)) - log(sum(Sj.new)))
   }else{
      return(log(sum(Sj)) - log(j-1 - sum(Sj.new)))
   }
}
lkernelprob.adap.cao <- function(X, Sj, Sj.new, j){
   if(sum(Sj - Sj.new) == -1){
      prob.of.new.ind = rep(0, sum(Sj==0))
      for(k in 1:sum(Sj==0)){
         S.old = Sj
         S.old[which(Sj == 0)[k]] = 1
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.denom = prob.of.new.ind[which(Sj == 0) == which(Sj - Sj.new == -1)] / sum(prob.of.new.ind)
      prob.of.new.ind = rep(0, sum(Sj.new==1))
      for(k in 1:sum(Sj.new==1)){
         S.old = Sj.new
         S.old[which(Sj.new == 1)[k]] = 0
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.numer = prob.of.new.ind[which(Sj.new == 1) == which(Sj.new - Sj == 1)] / sum(prob.of.new.ind)
      
      return(log(q.numer) - log(q.denom))
   }else{
      prob.of.new.ind = rep(0, sum(Sj==1))
      for(k in 1:sum(Sj==1)){
         S.old = Sj
         S.old[which(Sj == 1)[k]] = 0
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.denom = prob.of.new.ind[which(Sj == 1) == which(Sj - Sj.new == 1)] / sum(prob.of.new.ind)
      prob.of.new.ind = rep(0, sum(Sj.new==0))
      for(k in 1:sum(Sj.new==0)){
         S.old = Sj.new
         S.old[which(Sj.new == 0)[k]] = 1
         prob.of.new.ind[k] = 1/dShat(X, j, S.old)
      }
      q.numer = prob.of.new.ind[which(Sj.new == 0) == which(Sj.new - Sj == -1)] / sum(prob.of.new.ind)
      
      return(log(q.numer) - log(q.denom))
   }
}

dShat <- function(X, j, S){
  n = nrow(X)
  s = sum(S)
  tilde.Xj = as.matrix(X[, j])
  
  if(s == 0){
    return( sum((tilde.Xj)^2)/n )
  }else{
    Zj = as.matrix(X[, 1:(j-1)])
    X.Sj = as.matrix(Zj[, S>0])
    VhatX = sum(tilde.Xj^2)/n   # scalar
    Covhat = t(X.Sj)%*%matrix(tilde.Xj)/n   # s X 1 matrix
    res = VhatX - t(Covhat)%*%solve( t(X.Sj)%*%X.Sj/n )%*%Covhat
    return( res )
  }
}


logpost.cao <- function(X, j, S, qn, alpha.plus, U, gamma, b.hp){
   n = nrow(X)
   p = ncol(X)
   s = sum(S)
   
   dS.hat = dShat(X, j, S)
   # logpij = s*log(qn) + (p-j-s)*log(1-qn)
   logpij = s*log(qn) + (j-1-s)*log(1-qn)
   
   S.plus.j = S
   S.plus.j[j] = 1
   Zj = as.matrix(X[, 1:j])
   Uj = U[1:j, 1:j]
   X.Sj = as.matrix(Zj[, S.plus.j>0])
   UnS.j = t(X.Sj) %*% X.Sj + Uj[S.plus.j>0, S.plus.j>0]
   
   logpost.val = logpij + (n+alpha.plus-3)/2 * as.numeric(determinant(as.matrix(UnS.j[-nrow(UnS.j), -ncol(UnS.j)]), logarithm = TRUE)$mod) + 
      (alpha.plus-2)/2 * as.numeric(determinant(Uj, logarithm = TRUE)$mod) - (n+alpha.plus-2)/2 * as.numeric(determinant(UnS.j, logarithm = TRUE)$mod) +
      b.hp * sum(S * gamma[1:(j-1)]) * gamma[j]
   
   
   return( list(val = logpost.val, dj = dS.hat) )
}







