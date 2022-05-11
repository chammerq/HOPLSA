
PLSA_3D = function(flat_format,n_latent=2,niter=2,all_factors=FALSE){
  
  
  # Get names of dimensions
  dimA = sort(unique(flat_format[,1]))
  dimB = sort(unique(flat_format[,2]))
  dimC = sort(unique(flat_format[,3]))
  
  if(any(is.na(dimA),is.na(dimB),is.na(dimC))){
    stop('Missing values not allowed in tensor indices')
  }
  
  # Get index (Consider sorting for better cache coherency)
  i = match(flat_format[,1],dimA)
  j = match(flat_format[,2],dimB)
  k = match(flat_format[,3],dimC)
  qty = flat_format[,4]
  
  # get size
  na = length(dimA)
  nb = length(dimB)
  nc = length(dimC)
  
  
  # check size of n_latent
  if(n_latent<1){
    stop("Must have a postive number of latent variables",call. =FALSE)
  }
  if(n_latent>min(c(na,nb,nc))){
    stop("Latent space needs to smaller than smallest tensor dimension (for now)",call.=FALSE)
  }
  
  # Run C++ code
  cout = .Call('_HOPLSA_sparse_plsa_3d', PACKAGE = 'HOPLSA', i,j,k, qty, na, nb, nc, n_latent, niter)
  
  # Rename
  rnames = paste("Group",1:n_latent)
  rownames(cout$Pa_z) = dimA 
  rownames(cout$Pb_z) = dimB 
  rownames(cout$Pc_z) = dimC
  colnames(cout$Pa_z) = rnames 
  colnames(cout$Pb_z) = rnames
  colnames(cout$Pc_z) = rnames
  names(cout$Pz) = rnames
  
  if(all_factors){
    dZ = diag(cout$Pz)
    cout$Pz_a = dZ%*%t(cout$Pa_z)
    cout$Pz_b = dZ%*%t(cout$Pb_z)
    cout$Pz_c = dZ%*%t(cout$Pc_z)
    cout$Pz_a = sweep(cout$Pz_a,2,colSums(cout$Pz_a),"/")
    cout$Pz_b = sweep(cout$Pz_b,2,colSums(cout$Pz_b),"/")
    cout$Pz_c = sweep(cout$Pz_c,2,colSums(cout$Pz_c),"/")
    colnames(cout$Pz_a) = dimA 
    colnames(cout$Pz_b) = dimB 
    colnames(cout$Pz_c) = dimC
    rownames(cout$Pz_a) = rnames 
    rownames(cout$Pz_b) = rnames
    rownames(cout$Pz_c) = rnames
  }

  return(cout)
}

PLSA_4D = function(flat_format,n_latent=2,niter=2,all_factors=FALSE){
  
  
  # Get names of dimensions
  dimA = sort(unique(flat_format[,1]))
  dimB = sort(unique(flat_format[,2]))
  dimC = sort(unique(flat_format[,3]))
  dimD = sort(unique(flat_format[,4]))
  
  if(any(is.na(dimA),is.na(dimB),is.na(dimC),is.na(dimD))){
    stop('Missing values not allowed in tensor indices')
  }
  
  # Get index
  i = match(flat_format[,1],dimA)
  j = match(flat_format[,2],dimB)
  k = match(flat_format[,3],dimC)
  l = match(flat_format[,4],dimC)
  qty = flat_format[,5]
  
  # get size
  na = length(dimA)
  nb = length(dimB)
  nc = length(dimC)
  nd = length(dimD)
  
  
  # check size of n_latent
  if(n_latent<1){
    stop("Must have a postive number of latent variables",call. =FALSE)
  }
  if(n_latent>min(c(na,nb,nc,nd))){
    stop("Latent space needs to smaller than smallest tensor dimension (for now)",call.=FALSE)
  }
  
  # Run C++ code
  cout = .Call('_HOPLSA_sparse_plsa_4d', PACKAGE = 'HOPLSA', i,j,k,l, qty, na, nb, nc,nd, n_latent, niter)
  
  # Rename
  rnames = paste("Group",1:n_latent)
  rownames(cout$Pa_z) = dimA 
  rownames(cout$Pb_z) = dimB 
  rownames(cout$Pc_z) = dimC
  rownames(cout$Pd_z) = dimD
  colnames(cout$Pa_z) = rnames 
  colnames(cout$Pb_z) = rnames
  colnames(cout$Pc_z) = rnames
  colnames(cout$Pd_z) = rnames
  names(cout$Pz) = rnames
  
  if(all_factors){
    dZ = diag(cout$Pz)
    cout$Pz_a = dZ%*%t(cout$Pa_z)
    cout$Pz_b = dZ%*%t(cout$Pb_z)
    cout$Pz_c = dZ%*%t(cout$Pc_z)
    cout$Pz_d = dZ%*%t(cout$Pd_z)
    cout$Pz_a = sweep(cout$Pz_a,2,colSums(cout$Pz_a),"/")
    cout$Pz_b = sweep(cout$Pz_b,2,colSums(cout$Pz_b),"/")
    cout$Pz_c = sweep(cout$Pz_c,2,colSums(cout$Pz_c),"/")
    cout$Pz_d = sweep(cout$Pz_d,2,colSums(cout$Pz_d),"/")

    
    colnames(cout$Pz_a) = dimA 
    colnames(cout$Pz_b) = dimB 
    colnames(cout$Pz_c) = dimC
    colnames(cout$Pz_d) = dimD
    rownames(cout$Pz_a) = rnames 
    rownames(cout$Pz_b) = rnames
    rownames(cout$Pz_c) = rnames
    rownames(cout$Pz_d) = rnames
  }
  
  return(cout)
}

PLSA_2D = function(flat_format,n_latent=2,niter=2,all_factors=FALSE){
  
  
  # Get names of dimensions
  dimA = sort(unique(flat_format[,1]))
  dimB = sort(unique(flat_format[,2]))
  
  if(any(is.na(dimA),is.na(dimB))){
    stop('Missing values not allowed in tensor indices')
  }
  
  # Get index
  i = match(flat_format[,1],dimA)
  j = match(flat_format[,2],dimB)
  qty = flat_format[,3]
  
  # get size
  na = length(dimA)
  nb = length(dimB)

  
  # check size of n_latent
  if(n_latent<1){
    stop("Must have a postive number of latent variables",call. =FALSE)
  }
  if(n_latent>min(c(na,nb))){
    stop("Latent space needs to smaller than smallest tensor dimension (for now)",call.=FALSE)
  }
  
  # Run C++ code
  cout = .Call('_HOPLSA_sparse_plsa_2d', PACKAGE = 'HOPLSA', i,j, qty, na, nb, n_latent, niter)
  
  # Rename
  rnames = paste("Group",1:n_latent)
  rownames(cout$Pa_z) = dimA 
  rownames(cout$Pb_z) = dimB 
  colnames(cout$Pa_z) = rnames 
  colnames(cout$Pb_z) = rnames
  names(cout$Pz) = rnames
  
  if(all_factors){
    dZ = diag(cout$Pz)
    cout$Pz_a = dZ%*%t(cout$Pa_z)
    cout$Pz_b = dZ%*%t(cout$Pb_z)
    cout$Pz_a = sweep(cout$Pz_a,2,colSums(cout$Pz_a),"/")
    cout$Pz_b = sweep(cout$Pz_b,2,colSums(cout$Pz_b),"/")
    colnames(cout$Pz_a) = dimA 
    colnames(cout$Pz_b) = dimB 
    rownames(cout$Pz_a) = rnames 
    rownames(cout$Pz_b) = rnames
  }
  
  return(cout)
}

