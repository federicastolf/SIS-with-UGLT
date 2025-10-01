# function to simulate synthetic data
simulate_data = function(n, p, k, cov, sd_gamma, cp, sigma_sq, mseed){
  set.seed(mseed)
  q = dim(cov)[2]
  gamma = matrix(rnorm(q*p, sd = sd_gamma), q,p) # regression coefficient
  pilogit = cp*plogis(cov %*% gamma) 
  get_phi = simPhi(p, k, pilogit)
  Phi = get_phi$Phi
  theta_inv = rgamma(k,1,1)
  Lambda = matrix(0,p,k)
  for(j in 1:p){
    for(h in 1:k){
      Lambda[j,h]= Phi[j,h]/theta_inv[h]
    }
  }
  M = matrix(rnorm(n*k, 0, 1), ncol=k) # factor scores
  # data
  Y = M %*% t(Lambda) + sqrt(sigma_sq)*matrix(rnorm(n*p), nrow=n) 
  return(list("Y" = Y, "Lambda" = Lambda, "Phi"=Phi, "pivot"=get_phi$pivot,
              "delta"=get_phi$Delta))
}

# function for simulate Phi via gibbs-sampler
simPhi = function(p, k, pilogit, NiterGibbs = 200){
  # Initialization 
  phi = matrix(rbinom(p*k, 1, 0.2), p, k) 
  delta = phi
  l=rep(NA, k) # pivots
  l[1] = which(phi[,1]>0)[1]
  for(h in 2:k){
    w = which(phi[,h]>0) # counting 1's in the h column of Delta
    l[h] = w[!w %in% l[-h]][1] # take the first
    if(is.na(l[h])){
      delta[,h] = 0 # put all column to zeros
    }else{
      delta[1:l[h],h] = 0 # before pivots all zeros
      delta[l[h]:p, h] =  phi[l[h]:p,h] # from pivot and after equal to phi[,l]
    }
  }
  
  for(t in 2:NiterGibbs){
    for(h in 1:k){ 
      Lh = sort(l[-h])
      p_h = pilogit[,h]
      k = 0     
      # definisco prob di phi_h condizionatamente agli altri pivot
      for(pl in Lh){
        k = k+1
        if (pl-k>0){
          scale = 1-prod(1-pilogit[1:pl,h][-Lh[1:k]])
          p_h[pl] = min(p_h[pl]/scale, 1)
        }
        else{
          p_h[pl] = 1
        }
      }
      phi[,h] = rbinom(p,1,p_h) # generate phi
      w = which(phi[,h]>0)
      l[h] = w[!w %in% l[-h]][1]
      if(is.na(l[h])){
        delta[,h] = 0
      }else{
        delta[1:l[h],h] = 0
        delta[l[h]:p, h] =  phi[l[h]:p,h]
      }
    }
  }
  # return the last phi and delta
  return(list("Phi" = phi, "Delta"= delta, "pivot"=l))
}
