sim_mvOU<-function(phy,paths,z_0,OU_mat,F_mat,Sigma){ #Send phylogeny, function paths output, ancestral traits value, and optima values
  #Function to simulate multivariate OU traits, based on known parameters
  library(expm)
  library(ape)
  num_paths<-length(paths[[1]]) #Get number of paths to tips on phylogent
  N<-length(phy$tip.label)
  n_traits<-nrow(F_mat)
  I<-diag(n_traits)
  store_e_z_list<-list()

  mrca1 <- ape::mrca(phy) #Node numbers for MRCA of tips
  t_node <- ape::node.depth.edgelength(phy) #Time from root to node, starting with the tips
  t_root_MRCA <- matrix(t_node[mrca1], nrow=N, dimnames = list(phy$tip.label, phy$tip.label))# = Cij #Matrix with time from root to MRCA of pairs of tips - pulls out values of times that correspond with node numbers - integers
  t_tips <- t_node[1:N] #Times from root for tips
  t_MRCA_tips <- t_tips - t_root_MRCA #Ttimes from MRCA to tips

  P <- eigen(F_mat)$vectors
  lambdas <- eigen(F_mat)$values
  lambdas_matrix <- as.matrix(outer(lambdas,lambdas,"+"))
  inv_P <- solve(P)

  for(i in 1:num_paths){ #For each path back to root, calculate E[Z]
    #cat("Path:\n")
    #print(i)
    path<-paths[[1]][[i]] #Node numbers starting with root
    regimes<-paths[[2]][[i]] #Regimes for nodes starting with root
    branch_lengths<-paths[[3]][[i]] #Branch lengths between nodes starting with root
    e_z_des<-z_0 #First run of code start at ancestor, then use new descendent value as that of ancestor
    for(j in 1:length(branch_lengths)){ #For each segment of branch
      #Caclulation of E[Z]
      term1<-expm(-F_mat*branch_lengths[[j]])%*%as.matrix(e_z_des)
      term2<-I-expm(-F_mat*branch_lengths[[j]])

      term3<-t(as.matrix(OU_mat[regimes[[j]],]))
      #print(term3)
      e_z_des <- term1+term2%*%term3
      #cat("e_z :\n")
      #print(e_z_des)
    }
    store_e_z_list[[length(store_e_z_list)+1]]<-e_z_des
  }
  cov_M <- matrix(NA, nrow = N * n_traits, ncol= N * n_traits) #Number species * number traits

  for(i in 1:N){ #For each tip calculate covariance between tips
    for(j in 1:N){
      term1 <- expm(-F_mat*t_MRCA_tips[i,j]) #Effects of independent evolution of sp1
      term2 <- expm(t(-F_mat)*t_MRCA_tips[j,i]) #Effects of independent evolution of sp2

      #cat("t_root_MRCA[i,j]:\n")
      #print(t_root_MRCA[i,j])

      #cat("M_eigen:\n")
      if(t_root_MRCA[i,j] != 0){
        M_eigen <- ((1 / lambdas_matrix * (1 - exp(-lambdas_matrix*t_root_MRCA[i,j]))))
        }
      else(
        M_eigen <- matrix(0,nrow = n_traits, ncol = n_traits)
      )
      #print(M_eigen)
      #cat("M_transformed_Sigma:\n")
      M_transformed_Sigma <- inv_P %*% Sigma %*% t(Sigma) %*% t(inv_P) #This transforms the instantaneous rate matrix (R) into the eigen-space of F.
      #print(M_transformed_Sigma)

      M_integral <- P %*% (M_eigen * M_transformed_Sigma) %*% t(P)

      cov_traits_ij <- term1 %*% M_integral %*% term2 #Covariance between all traits for sp1 and sp2 - results in matrix
      #cat("Covariance matrix:\n")
      #print(paste(i,j))

      #print(cov_traits_ij)
      start_row <- (i - 1) * n_traits + 1
      end_row <- i * n_traits
      start_col <- (j - 1) * n_traits + 1
      end_col <- j * n_traits

      cov_M[start_row:end_row,start_col:end_col]<-cov_traits_ij #Stores matrix for sp 1 and sp2 comparison in larger matrix
    }

  }
  return(list(store_e_z_list,cov_M))
}
