
#' Compute the density function
#'
#' \code{compute_den} renvoie la densité de proba correspondant à 'den_name'
#' sur l'interval 'interval'
#' calculée en nb_eval_pts points équidistants sur l'interval
#'
#'
#'
#' [m,s,p0] parameters are :
#' Dans le cas Gaussien : 2 scalaires, m pour la moyenne et s pour
#'   l'ecart-type.
#' Dans le cas Cauchy : 2 scalaires, m pour la position et s pour le facteur
#'   d'échelle.
#' Dans le cas Uniforme : 2 scalaires, m pour le min et s pour le max.
#' Dans le cas Exponentielle : 1 scalaires, s pour le facteur d'échelle.
#' Dans le cas melange Gaussien : 3 vecteurs, m pour les moyennes
#'   s pour les ecart-types et p0 pour les proportions.
#' Dans le cas melange uniforme, 3 vecteurs : m qui correspond aux minima,
#'   s qui correspond aux maxima et p0 qui correspond aux proportions des
#'   differentes uniformes.
#'
#'
#'
#'
#' @param den_name chaine de caractere definissant la loi d'intérêt.
#'   Peut être : 'Gauss', 'MixGauss', 'Unif', 'MixUnif', 'Cauchy', 'Exp'.
#' @param interval un vecteur de dimension 1 contenant 2 éléments, a et b,
#'   correspondant aux extrêmités de l'interval.
#' @param x0 on peut directement donner les points où calculer la densité via
#'   x0. Dans ce cas il ne faut pas donner d'interval.
#' @param nb_eval_pts nombre de bin où calculer la densité.
#' @param [m,s,p0] les paramètres de la loi d'intérêt.
#'   Voir les détails pour plus d'info.
#'
#'
#'
#'
#'
#' @return
#' une liste contenant 2 vecteurs de dimension 1 et de taille nb_eval_pts
#' le premier vecteur 'x' contient les nb_eval_pts points équidistants de
#' l'interval
#' le second vecteur 'y' contient les nb_eval_pts valeurs de la densité
#' calculée en x
#'
#'
#' @importFrom stats dcauchy dexp dnorm dunif
compute_den_1D <- function(den_name, interval = NULL, x0 = 0, nb_eval_pts = 100,
                           m = 0, s = 1, p0 = 1){
  
  
  
  if (is.null(interval)){
    x <- x0
  }else{
    
    # # bornes de l'interval
    # a <- interval[1];
    # b <- interval[2];
    #
    # # points où la densité est calculée (nb_eval_pts points équidistants de
    # # interval)
    # x <- 0:(nb_eval_pts-1);
    # x <- (b-a)*x/(nb_eval_pts-1) + a;
    
    # x <- seq(from = interval['min'], to = interval['max'], length.out = nb_eval_pts)
    
    x0 <- randtoolbox::sobol(nb_eval_pts)
    x <- x0 * (interval['max'] - interval['min']) + interval['min']
    
  }
  
  
  
  
  
  # les fonctions locales suivantes permettent de calculer la densité d'intérêt
  # en x
  
  myGaussDen <- function(x, m, s, p0){
    # densité  gaussiennes calculée aux points x
    return(stats::dnorm(x = x, mean = m, sd = s))
  }
  
  myMixGaussDen <- function(x, m, s, p0){
    # densité melange de gaussiennes calculée aux points x
    y = array(data = 0, dim = length(x))
    for (i in 1:length(p0)){
      y = y + p0[i]*stats::dnorm(x = x, mean = m[i], sd = s[i])
    }
    return(y)
  }
  
  myCauchyDen <- function(x, m, s, p0){
    # densité de Cauchy calculée aux points x
    return(stats::dcauchy(x = x, location = m, scale = s))
  }
  
  myUnifDen <- function(x, m, s, p0){
    # densité uniforme calculée aux points x
    return(stats::dunif(x = x, min = m, max = s))
  }
  
  myMixUnifDen <- function(x, m, s, p0){
    # densité melange d'e gaussiennes'uniformes calculée aux points x
    y <- array(data = 0, dim = length(x) )
    for (i in 1:length(p0)){
      y <- y + p0[i]*stats::dunif(x = x, min = m[i], max = s[i])
    }
    return(y)
  }
  
  myExpDen <- function(x, m, s, p0){
    # densité exponentielle calculée aux points x
    return(stats::dexp(x = x, rate = s))
  }
  
  myLaplaceDen <- function(x, m, s, p0){
    # densité de laplace calculée aux points x
    return(exp(-abs(x-m)/s)/(2*s))
  }
  
  myMixLaplaceDen <- function(x, m, s, p0){
    # densité melange de laplace calculée aux points x
    y <- array(data = 0, dim = length(x) )
    for (i in 1:length(p0)){
      y <- y + p0[i]*exp(-abs(x-m[i])/s[i])/(2*s[i])
    }
    return(y)
  }
  
  # appel de la fonction correspondant à la densité d'intérêt
  nomfunc <- paste('my', den_name, 'Den', sep = '')
  y <- do.call(nomfunc, args = list(x = x, m = m, s = s, p0 = p0))
  
  
  # resultats mis sous forme de matrice
  # cad densité et points où la densité est calculée
  # matfinale <- cbind(x, y)
  # names(matfinale) <- c('x', 'y')
  # return(matfinale)
  
  return(list('x'=x,'y'=y))
  
}



#' Generate samples from the desired density function
#'
#' \code{generate_ech} renvoie un échantillon de taille ech_size de la densité
#' de proba correspondant à 'den_name'
#'
#'
#'
#' [m,s,p0] parameters are :
#' Dans le cas Gaussien : 2 scalaires, m pour la moyenne et s pour
#'   l'ecart-type.
#' Dans le cas Cauchy : 2 scalaires, m pour la position et s pour le facteur
#'   d'échelle.
#' Dans le cas Uniforme : 2 scalaires, m pour le min et s pour le max.
#' Dans le cas Exponentielle : 1 scalaires, s pour le facteur d'échelle.
#' Dans le cas melange Gaussien : 3 vecteurs, m pour les moyennes
#'   s pour les ecart-types et p0 pour les proportions.
#' Dans le cas melange uniforme, 3 vecteurs : m qui correspond aux minima,
#'   s qui correspond aux maxima et p0 qui correspond aux proportions des
#'   differentes uniformes.
#'
#'
#'
#'
#' @param den_name chaine de caractere definissant la loi d'intérêt.
#'   Peut être : 'Gauss', 'MixGauss', 'Unif', 'MixUnif', 'Cauchy', 'Exp'.
#' @param ech_size taille de l'échantillon à générer.
#' @param [m,s,p0] les paramètres de la loi d'intérêt.
#'   Voir les détails pour plus d'info.
#'
#'
#'
#'
#' @importFrom stats rcauchy rexp rnorm runif rmultinom
#' @return
#' un vecteur de dimension 1 et de taille ech_size
#'
generate_ech_1D <- function(den_name, ech_size = 1000, m = 0, s = 1, p0 = 1){
  
  
  # les fonctions locales suivantes permettant de générer un échantillon de la
  # densité d'intérêt
  
  myGaussEch <- function(ech_size, m, s, p0){
    # échantillon d'une loi gaussiennes
    return(rnorm(n = ech_size, mean = m, sd = s))
  }
  
  myMixGaussEch <- function(ech_size, m, s, p0){
    # échantillon d'une loi melange de gaussiennes
    
    # génération d'un échantillon de taille ech_size selon une multinomiale
    ztilde = rmultinom(n = ech_size, size = 1, prob = p0)
    
    y <- array(data = 0, dim = ech_size)
    ind_deb <- 1
    ind_fin <- 0
    for (i in 1:length(p0)){
      
      lesi <- ztilde[i, (ztilde[i,]==1)]
      nbi <- length(lesi)
      if (nbi==0){
        print(paste('nbi Gauss = ',nbi))
      }
      
      if (nbi>0){
        ind_fin <- ind_fin + nbi
        y[ind_deb:ind_fin] <- rnorm(n = nbi, mean = m[i], sd = s[i]);
        ind_deb <- ind_deb + nbi
        
      }
      
    }
    return(y)
  }
  
  myCauchyEch <- function(ech_size, m, s, p0){
    # densité de Cauchy
    return(rcauchy(n = ech_size, location = m, scale = s))
  }
  
  myUnifEch <- function(ech_size, m, s, p0){
    # densité uniforme
    return(runif(n = ech_size, min = m, max = s))
  }
  
  myMixUnifEch <- function(ech_size, m, s, p0){
    # densité melange d'uniformes
    
    # génération d'un échantillon de taille ech_size selon une multinomiale
    ztilde = rmultinom(n = ech_size, size = 1, prob = p0)
    
    y <- array(data = NA, dim = ech_size)
    
    ind_deb <- 1
    ind_fin <- 0
    for (i in 1:length(p0)){
      
      lesi <- ztilde[i, (ztilde[i,]==1)]
      
      nbi <- length(lesi)
      
      
      if (nbi>0){
        ind_fin <- ind_fin + nbi
        y[ind_deb:ind_fin] <- runif(n = nbi, min = m[i], max = s[i])
        ind_deb <- ind_deb + nbi
      }
      
    }
    return(y)
  }
  
  myExpEch <- function(ech_size, m, s, p0){
    # densité exponentielle
    return(rexp(n = ech_size, rate = s))
  }
  
  
  myLaplaceEch <- function(ech_size, m, s, p0){
    # densité uniforme
    ui <- runif(n = ech_size, min = -0.5, max = 0.5)
    return( m  - s*sign(ui)*log(1-2*abs(ui)))
  }
  
  
  
  
  
  # appel de la fonction correspondant à la densité d'intérêt
  nomfunc <- paste('my', den_name, 'Ech', sep = '')
  y <- do.call(nomfunc, args = list(ech_size = ech_size, m = m, s = s, p0 = p0))
  
  
  
  # densité et points où la densité est calculée
  return(y)
  
}





#' renvoie un échantillon de taille n de la densité de proba correspondant à 'nom' 
#'
#' @param
#' den_name : chaine de caractere definissant la loi d'intérêt. Peut être :
#'            'Gauss', 'MixGauss', 'Unif', 'MixUnif', 'Cauchy', 'Exp'
#' ech_size : taille de l'échantillon à générer
#' m, s, p0 : les paramètres de la loi d'intérêt
#'            Dans le cas Gaussien : 2 scalaires, m pour la moyenne et s pour l'ecart-type
#'            Dans le cas Cauchy : 2 scalaires, m pour la position et s pour le facteur 
#'              d'échelle
#'            Dans le cas Uniforme : 2 scalaires, m pour le min et s pour le max
#'            Dans le cas Exponentielle : 1 scalaire, s pour le facteur d'échelle
#'            Dans le cas melange Gaussien : 3 vecteurs, m pour les moyennes
#'              s pour les ecart-types et p0 pour les proportions
#'            Dans le cas melange uniforme, 3 vecteurs : p0 qui correspond aux 
#'              proportions des differentes uniformes, m qui correspond aux minima 
#'              et s qui correspond aux maxima
#'  
#' @output
#' un vecteur de dimension 1 et de taille n
generate_ech_mD <- function(den_name, ech_size = 1000, m = 0, s = 1, p0 = 1){
  
  if (length(dim(m)) <= 1){
    d <- length(m)
  }else{
    d <- length(m[1,])
  }
  
  
  # les fonctions locales suivantes permettant de générer un échantillon de la densité 
  # d'intérêt
  
  myGaussEch <- function(n, m, s, p0){
    # échantillon d'une loi gaussiennes multidim
    return(LaplacesDemon::rmvn(n=n, mu=m, Sigma=s))
  }
  
  myCauchyEch <- function(n, m, s, p0){
    # échantillon d'une loi de Cauchy multidim
    return(LaplacesDemon::rmvc(n=n, mu=m, S=s))
  }
  
  # myFunc1Ech <- function(n, m, s, p0){
  #   yi <- array(data=NA, dim=c(n, 2))
  #   # X suit une Exp(lambda=1) 
  #   xi <- rexp(n = n, rate = 1)
  #   yi[, 1] <- xi
  #   # et Y sachant X=x suit une uniforme sur [0, x]
  #   #yi <- sapply(xi, function(x) runif(n=n, min=0, max=x))
  #   yi[, 2] <- xi*runif(n=n, min=0, max=1)
  #   # print(xi)
  #   # print(yi)
  #   #return(array(c(xi,yi)))
  #   return(yi)
  #   
  # }
  
  myUnifCercleEch <- function(n, m, s, p0){
    d <- length(m)
    x_max <- m + s
    x_min <- m - s
    
    yi <- array(data=NA, dim=c(n, d))
    x_tmp <- array(data=NA, dim=d)
    nbpts <- 0
    
    while (nbpts < n){
      for (no_d in 1:d){
        x_tmp[no_d] <- runif(n=1, min=x_min[no_d], max=x_max[no_d])
      }
      
      
      if (sum((x_tmp - m)^2) <= s^2){
        nbpts <- nbpts + 1
        yi[nbpts, ] <- x_tmp
      }
      
    }
    return(yi)
  }
  
  
  myMixGaussEch <- function(n, m, s, p0){
    # échantillon d'une loi melange de gaussiennes
    
    # génération d'un échantillon de taille n selon une multinomiale
    ztilde = rmultinom(n = n, size = 1, prob = p0)
    
    yi <- array(data = 0, dim = c(n, d))
    
    
    ind_deb <- 1
    ind_fin <- 0
    for (i in 1:length(p0)){
      
      lesi <- ztilde[i, (ztilde[i,]==1)]
      nbi <- length(lesi)
      
      
      if (nbi>0){
        ind_fin <- ind_fin + nbi
        yi[ind_deb:ind_fin,] <- MASS::mvrnorm(n = nbi, mu = m[i,], Sigma = s[,,i])
        #yi[ind_deb:ind_fin,] <- rmvn(n=nbi, mu=m[i,], Sigma=s[,,i])
        ind_deb <- ind_deb + nbi
        
      }
      
    }
    return(yi)
  }
  
  
  
  myMixCauchyEch <- function(n, m, s, p0){
    # échantillon d'une loi melange de Cauchy
    
    # génération d'un échantillon de taille n selon une multinomiale
    ztilde = rmultinom(n = n, size = 1, prob = p0)
    
    yi <- array(data = 0, dim = c(n, d))
    
    ind_deb <- 1
    ind_fin <- 0
    for (i in 1:length(p0)){
      
      lesi <- ztilde[i, (ztilde[i,]==1)]
      nbi <- length(lesi)
      
      if (nbi>0){
        ind_fin <- ind_fin + nbi
        yi[ind_deb:ind_fin,] <- LaplacesDemon::rmvc(n = nbi, mu = m[i,], S = s[,,i])
        ind_deb <- ind_deb + nbi
      }
      
    }
    return(yi)
  }
  
  
  
  
  
  
  # mySpiralEch <- function(n, m, s, p0){
  #   t <- runif(n = n, min = 0, max = 14)
  #   w <- rmvn(n=n, mu=c(0,0,0), Sigma=diag(c(1/4, 1/4, 1/4)))
  #   x1 <- (13-0.5*t)*cos(t)
  #   x2 <- -(13-0.5*t)*sin(t)
  #   x3 <- t
  #   x <- array(data=NA, dim=c(n,3))
  #   x[, 1] <- x1
  #   x[, 2] <- x2
  #   x[, 3] <- x3
  #   
  #   return(x+w)
  # }
  
  
  
  
  # appel de la fonction correspondant à la densité d'intérêt
  nomfunc <- paste('my', den_name, 'Ech', sep = '')
  y <- do.call(nomfunc, args = list(n = ech_size, m = m, s = s, p0 = p0))
  
  
  
  # densité et points où la densité est calculée
  return(y)
  
}









#' Generate N n-samples of all test laws in dimension d
#'
#' \code{gen_ech} Generate N n-samples of all test laws in dimension d
#'
#' @param n sample size
#' @param N Number of samples
#' @param d the dimension of the samples
#'
#' @return genere
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @importFrom LaplacesDemon is.positive.definite as.positive.definite
#' @export
gen_N_ech <- function(n, d, N){
  # cette fonction initialise les paramètres de toutes les lois test et génère
  # les échantillons associés
  
  ech_file_name <- sprintf('ech_d%iN%in%i', d, N, n)
  if (!file.exists(ech_file_name)){
    if (d==1){
      nom_loi <- c('Gauss', 'Cauchy', 'Unif', 'Exp', 'Mix Gauss', 'skewed',
                   'strong skewed', 'kurtotic', 'outlier', 'bimodal',
                   'separate bimodal', 'skewed bimodal', 'trimodal', 'Bart',
                   'double Bart', 'asymetric Bart', 'asymetric double Bart',
                   'smooth comb', 'discrete comb', 'Mix Unif')
      nom_loi_abr <- c('G', 'C','U', 'E', 'MG', 'Sk', 'Sk+', 'K', 'O', 'Bi',
                       'SBi', 'SkBi', 'T', 'B', 'DB', 'AB', 'ADB', 'SC', 'DC',
                       'MU')
      
      # paramètres des lois pour les tests
      param_loi <- list()
      
      
      param <- list()
      param$nom_func <- 'Gauss'
      param$m <- 0
      param$s <- 1
      # param$supp <- c(-3*param$s + param$m, 3*param$s + param$m)
      param_loi[['Gauss']] <- param
      
      param <- list()
      param$nom_func <- 'Cauchy'
      param$m <- 0
      param$s <- 1
      # param$supp <- c(-100, 100)
      param_loi[['Cauchy']] <- param
      
      param <- list()
      param$nom_func <- 'Unif'
      param$m <- 0
      param$s <- 1
      # param$supp <- c(-1, 2)
      param_loi[['Unif']] <- param
      
      param <- list()
      param$nom_func <- 'Exp'
      param$s <- 1
      # param$supp <- c(-0.5, 3)
      param_loi[['Exp']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(0, 3)
      param$s <- c(1, 1/3)
      param$p0 <- c(0.5, 0.5)
      # param$supp <- c(-3, 5)
      param_loi[['Mix Gauss']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(0, 0.5, 13/12)
      param$s <- c(1, 2/3, 5/9)
      param$p0 <- c(1/5, 1/5, 3/5)
      # param$supp <- c(-3, 3)
      param_loi[['skewed']] <- param
      
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- 0:7
      param$m <- 3*((2/3)^j - 1)
      param$s <- (2/3)^j
      param$p0 <- rep(x = 1/8, times = 8)
      # param$supp <- c(-3, 3)
      param_loi[['strong skewed']] <- param
      rm(list='j')
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(0, 0)
      param$s <- c(1, 1/10)
      param$p0 <- c(2/3, 1/3)
      # param$supp <- c(-3, 3)
      param_loi[['kurtotic']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(0, 0)
      param$s <- c(1, 1/10)
      param$p0 <- c(1/10, 9/10)
      # param$supp <- c(-3, 3)
      param_loi[['outlier']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(-1, 1)
      param$s <- c(2/3, 2/3)
      param$p0 <- c(1/2, 1/2)
      # param$supp <- c(-3, 3)
      param_loi[['bimodal']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(-3/2, 3/2)
      param$s <- c(1/2, 1/2)
      param$p0 <- c(1/2, 1/2)
      # param$supp <- c(-3, 3)
      param_loi[['separate bimodal']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(0, 3/2)
      param$s <- c(1, 1/3)
      param$p0 <- c(3/4, 1/4)
      # param$supp <- c(-3, 3)
      param_loi[['skewed bimodal']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- c(-6/5, 6/5, 0)
      param$s <- c(3/5, 3/5, 1/4)
      param$p0 <- c(9/20, 9/20, 1/10)
      # param$supp <- c(-3, 3)
      param_loi[['trimodal']] <- param
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- 0:4
      param$m <- c(0, j/2 - 1)
      param$s <- c(1, 0.1, 0.1, 0.1, 0.1, 0.1)
      param$p0 <- c(0.5, 0.1, 0.1, 0.1, 0.1, 0.1)
      # param$supp <- c(-3, 3)
      param_loi[['Bart']] <- param
      rm(list='j')
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- 0:6
      param$m <- c(-1, 1, (j-3)/2)
      param$s <- c(2/3, 2/3, rep(1/100, times=7))
      param$p0 <- c(49/100, 49/100, rep(1/350, times=7))
      # param$supp <- c(-3, 3)
      param_loi[['double Bart']] <- param
      rm(list='j')
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- -2:2
      param$m <- c(0, j + 1/2)
      param$s <- c(1, 2^(-j)/10)
      param$p0 <- c(1/2, 2^(1-j)/31)
      # param$supp <- c(-3, 3)
      param_loi[['asymetric Bart']] <- param
      rm(list='j')
      
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- 1:3
      param$m <- c(-1, 1, -j/2, j/2)
      param$s <- c(2/3, 2/3, rep(1/100, times=3), rep(x = 7/100, times = 3))
      param$p0 <- c(46/100, 46/100, rep(1/300,times=3), rep(7/300,times=3))
      # param$supp <- c(-3, 3)
      param_loi[['asymetric double Bart']] <- param
      rm(list='j')
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j <- 0:5
      param$m <- (65-96*(0.5^j))/21
      param$s <- (32/63)/(2^j)
      param$p0 <- 2^(5-j)/63
      # param$supp <- c(-3, 3)
      param_loi[['smooth comb']] <- param
      rm(list='j')
      
      
      param <- list()
      param$nom_func <- 'MixGauss'
      j1 <- 0:2
      j2 <- 8:10
      param$m <- c((12*j1-15)/7, 2*j2/7)
      param$s <- c(rep(2/7, times=3), rep(1/21, times=3))
      param$p0 <- c(rep(x = 2/7, times=3), rep(x = 1/21, times=3))
      # param$supp <- c(-3, 3)
      param_loi[['discrete comb']] <- param
      rm(list=c('j1', 'j2'))
      
      param <- list()
      param$nom_func <- 'MixUnif'
      param$p0 <- c(0.04, 0.145, 0.085, 0.05, 0.14, 0.2, 0.14, 0.2);
      param$m  <- c(0, 0.15, 0.20, 3/8, 4/8, 0.60, 0.80, 7/8)
      param$s <- c(0.15, 0.20, 3/8, 4/8, 0.60, 0.80, 7/8, 1)
      # param$supp <- c(0, 1)
      param_loi[['Mix Unif']] <- param
      
      
      
      # param <- list()
      # param$nom_func <- 'MixLaplace'
      # param$p0 <- c(0.04, 0.145, 0.085, 0.05, 0.14, 0.2, 0.14, 0.2);
      # param$m  <- c(0, 0.15, 0.20, 3/8, 4/8, 0.60, 0.80, 7/8)
      # param$s <- c(0.15, 0.20, 3/8, 4/8, 0.60, 0.80, 7/8, 1)
      # param$supp <- c(-1, 1)
      
      
      
      rm(param)
      
      
      
      
      print('fin d\'initialisation des paramètres')
      
      
      
      # Pour chacune des densités test on génère N échantillons composés de
      # n observations de la loi correspondante
      
      print('Génération des échantillons...')
      
      nomdim_echs <- list(NULL, 1:N, nom_loi )
      echs <- array(data = NA, dim = c(n, N, length(nom_loi)))
      dimnames(echs) <- nomdim_echs
      
      for (no_loi in 1:length(nom_loi)){
        print(nom_loi[no_loi])
        type <- param_loi[[no_loi]]$nom_func
        print(type)
        m <- param_loi[[no_loi]]$m
        s <- param_loi[[no_loi]]$s
        p0 <- param_loi[[no_loi]]$p0
        
        echs[, , no_loi] <- replicate(N, generate_ech_1D(den_name = type, ech_size = n,
                                                         m = m, s = s, p0 = p0))
        
      }
      
      
      
      
      
    }# en dimension > 1 : NON TESTEE
    else{
      nom_loi <- c('Gauss UnCorr', 'Gauss Corr','Cauchy', 'Mix Cauchy 1',
                   'Uniform', 'Strong Skewed', 'Skewed 2', 'Dumbbell',
                   'Kurtotic', 'Bimodal', 'Separate Bimodal',
                   'Asymmetric Bimodal', 'Trimodal', 'Fountain',
                   'Double Fountain', 'Asymmetric Fountain')
      nom_loi_abr <- c('UG', 'CG','C', 'MC', 'U', 'Sk+', 'Sk', 'D', 'K', 'Bi',
                       'SBi', 'ABi', 'T', 'F', 'DF', 'AF')
      
      
      # paramètres des lois pour les tests
      param_loi <- list()
      
      
      
      print('loi Gauss non correlee')
      # loi gaussienne multidim non correlee
      param <- list()
      param$nom_func <- 'Gauss'
      param$m <- rep(0, times=d)
      
      if (d%%2==0){
        param$s <- diag(c(rep(1/4, times=d/2 ), rep(1, times=d/2)))
        
      }else{
        param$s <- diag(c(rep(1/4, times=(d-1)/2 ), rep(1, times=(d+1)/2)))
      }
      
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*param$s[no_d, no_d] + param$m[no_d]
        param$supp[no_d, 2] <-  3*param$s[no_d, no_d] + param$m[no_d]
      }
      param_loi[['Gauss UnCorr']] <- param
      
      
      
      print('loi Gauss correlee')
      # loi gaussienne multidim correlee
      param <- list()
      param$nom_func <- 'Gauss'
      param$m <- rep(0, times=d)
      
      param$s <- array(data=9/10, dim =c(d, d))
      diag(param$s) <- 1
      
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*param$s[no_d, no_d] + param$m[no_d]
        param$supp[no_d, 2] <-  3*param$s[no_d, no_d] + param$m[no_d]
      }
      param_loi[['Gauss Corr']] <- param
      
      # print('loi Gauss correlee')
      # # loi gaussienne multidim correlee
      # param <- list()
      # param$nom_func <- 'Gauss'
      # param$m <- rep(0, times=d)
      #
      # param$s <- array(data=NA, dim =c(d, d))
      # param$s[1,1] <- 6.22
      # param$s[2,2] <- 2.93
      # param$s[3,3] <- 4.58
      #
      # param$s[1,2] <- -1.56
      # param$s[1,3] <- -2.05
      # param$s[2,3] <- 0.20
      # param$s[2,1] <- -1.56
      # param$s[3,1] <- -2.05
      # param$s[3,2] <- 0.20
      # if (!is.positive.definite(param$s)){
      #   param$s <- as.positive.definite(param$s)
      # }
      # param$supp <- matrix(data=NA, ncol=2, nrow=d)
      # for (no_d in 1:d){
      #   param$supp[no_d, 1] <- -3*param$s[no_d, no_d] + param$m[no_d]
      #   param$supp[no_d, 2] <-  3*param$s[no_d, no_d] + param$m[no_d]
      # }
      # param_loi[['Gauss Corr']] <- param
      
      
      
      
      
      
      # print('loi Spiral')
      # # il faut trouver la formule exacte pour la densité
      # param <- list()
      # param$nom_func <- 'Spiral'
      # param_loi[['Spiral']] <- param
      
      
      print('loi Cauchy')
      # loi Cauchy multidim
      param <- list()
      param$nom_func <- 'Cauchy'
      param$m <- rep(0, times=d)
      param$s <- array(data=9/10, dim =c(d, d))
      diag(param$s) <- 1
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*param$s[no_d, no_d] + param$m[no_d]
        param$supp[no_d, 2] <-  3*param$s[no_d, no_d] + param$m[no_d]
      }
      param_loi[['Cauchy']] <- param
      
      # print('loi Cauchy')
      # # loi Cauchy multidim
      # param <- list()
      # param$nom_func <- 'Cauchy'
      # param$m <- rep(0, times=d)
      # # param$s <- array(data=9/10, dim =c(d, d))
      # # diag(param$s) <- 1
      # param$s[1,1] <- 6.65
      # param$s[2,2] <- 6.33
      # param$s[3,3] <- 6.79
      #
      # param$s[1,2] <- 0.10
      # param$s[1,3] <- 0.25
      # param$s[2,3] <- 0.09
      # param$s[2,1] <- 0.10
      # param$s[3,1] <- 0.25
      # param$s[3,2] <- 0.09
      # param$supp <- matrix(data=NA, ncol=2, nrow=d)
      # for (no_d in 1:d){
      #   param$supp[no_d, 1] <- -3*param$s[no_d, no_d] + param$m[no_d]
      #   param$supp[no_d, 2] <-  3*param$s[no_d, no_d] + param$m[no_d]
      # }
      # param_loi[['Cauchy']] <- param
      
      
      
      # melange de Cauchy
      print('melange de Cauchy 1')
      param <- list()
      param$nom_func <- 'MixCauchy'
      param$m <- array(data=NA, dim=c(3,d))
      # param$m[1, ] <- rep(0, times=d)
      # param$m[2, ] <- rep(1/2, times=d)
      # param$m[3, ] <- rep(13/12, times=d)
      param$m[1, ] <- rep(0, times=d)
      param$m[2, ] <- rep(5, times=d)
      param$m[3, ] <- rep(10, times=d)
      param$s <- array(data=NA, dim =c(d, d, 3))
      param$s[, , 1] <- diag(rep(1, times=d))
      param$s[, , 2] <- diag(rep(4/9, times=d))
      param$s[, , 3] <- diag(rep(25/81, times=d))
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param$p0 <- c(1/5, 1/5, 3/5)
      param_loi[['Mix Cauchy 1']] <- param
      
      
      # Uniforme sur un disque
      # equation du cercle de centre (a,b) et de rayon r : (x-a)^2 + (y-b)^2 +(z-c)^2 = r^2
      param <- list()
      param$nom_func <- 'UnifCercle'
      # le centre du cercle
      param$m <- rep(2, times=d)
      # le rayon
      param$s <- 1
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- param$m[no_d] - param$s
        param$supp[no_d, 2] <- param$m[no_d] + param$s
      }
      #param$supp <- c(param$m - param$s, param$m + param$s)
      param_loi[['Uniform']] <- param
      # pb pour la loi uniforme lors du calcul de la vraie densité
      
      
      
      print('Strong skewed')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(8,d))
      for (no_d in 1:d){
        param$m[, no_d] <- (1 - (4/5)^(0:7))*(-1)^(no_d+1)*3
      }
      param$s <- array(data=NA, dim =c(d, d, 8))
      for (l in 1:8){
        param$s[, , l] <- array(data=-9*(4/5)^(2*(l-1))/10, dim=c(d,d))
      }
      for (no_d in 1:d){
        param$s[no_d, no_d, ] <- (4/5)^(2*(0:7))
      }
      for (l in 1:8){
        if (!LaplacesDemon::is.positive.definite(param$s[, , l])){
          param$s[, , l] <- LaplacesDemon::as.positive.definite(param$s[, , l])
        }
      }
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param$p0 <- rep(1/8, times=8)
      param_loi[['Strong Skewed']] <- param
      
      
      
      
      print('Skewed 2')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(3,d))
      param$m[1, ] <- rep(0, times=d)
      param$m[2, ] <- rep(5, times=d)
      param$m[3, ] <- rep(10, times=d)
      param$s <- array(data=NA, dim =c(d, d, 3))
      param$s[, , 1] <- diag(rep(1, times=d))
      param$s[, , 2] <- diag(rep(4/9, times=d))
      param$s[, , 3] <- diag(rep(25/81, times=d))
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param$p0 <- c(1/5, 1/5, 3/5)
      param_loi[['Skewed 2']] <- param
      
      
      
      
      
      print('Dumbbell')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(3, d))
      param$m[1,] <- (-1)^(1:d)*3/2
      param$m[2,] <- (-1)^(0:(d-1))*3/2
      param$m[3,] <- rep(0, times=d)
      param$s <- array(data=NA, dim =c(d, d, 3))
      param$s[,, 1] <- diag(rep(9/16, times=d))
      param$s[,, 2] <- diag(rep(9/16, times=d))
      param$s[,, 3] <- array(data = -9*18/(16*25), dim=c(d,d))
      diag(param$s[,, 3]) <- 9*4/(16*5)
      for (l in 1:3){
        if (!LaplacesDemon::is.positive.definite(param$s[, , l])){
          param$s[, , l] <- LaplacesDemon::as.positive.definite(param$s[, , l])
        }
      }
      param$p0 <- c(4/11, 4/11, 3/11)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Dumbbell']] <- param
      
      
      
      
      print('Kurtotic')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2, d))
      param$m[1,] <- rep(0, times=d)
      param$m[2,] <- rep(0, times=d)
      param$s <- array(data=NA, dim =c(d, d, 2))
      param$s[,, 1] <- array(data = 1, dim=c(d,d))
      diag(param$s[,,1]) <- c(1, 4, rep(4, times=d-2))
      param$s[,, 2] <- array(data = -1/3, dim=c(d,d))
      diag(param$s[,, 2]) <- 4/9
      for (l in 1:2){
        if (!LaplacesDemon::is.positive.definite(param$s[, , l])){
          param$s[, , l] <- LaplacesDemon::as.positive.definite(param$s[, , l])
        }
      }
      param$p0 <- c(2/3, 1/3)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Kurtotic']] <- param
      
      
      
      
      print('bimodal')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2,d))
      param$m[1,] <- c(-1, rep(0, times=d-1))
      param$m[2,] <- c(1, rep(0, times=d-1))
      param$s <- array(data=2/9, dim =c(d, d, 2))
      diag(param$s[,, 1]) <- 4/9
      diag(param$s[,, 2]) <- 4/9
      # param$s[, , 1] <- array(data=c(4/9, 2/9, 2/9, 4/9), dim=c(2,2))
      # param$s[, , 2] <- array(data=c(4/9, 2/9, 2/9, 4/9), dim=c(2,2))
      param$p0 <- rep(1/2, times=2)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Bimodal']] <- param
      
      
      
      print('Separate bimodal')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2,d))
      param$m[1,] <- c(-1, rep(1, times=d-1))
      param$m[2,] <- rep(0, times=d)
      param$s <- array(data=NA, dim =c(d, d, 2))
      param$s[, , 1] <- array(data=1/3, dim=c(d,d))
      diag(param$s[, , 1]) <- 4/9 #array(data=c(4/9, 1/3, 1/3, 4/9), dim=c(2,2))
      param$s[, , 2] <- diag(rep(4/9, times=d)) #array(data=c(4/9, 0, 0, 4/9), dim=c(2,2))
      for (l in 1:2){
        if (!LaplacesDemon::is.positive.definite(param$s[, , l])){
          param$s[, , l] <- LaplacesDemon::as.positive.definite(param$s[, , l])
        }
      }
      param$p0 <- rep(1/2, times=2)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Separate Bimodal']] <- param
      
      
      
      
      print('Asymmetric Bimodal')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2,d))
      param$m[1,] <- (-1)^(0:(d-1))
      param$m[2,] <- (-1)^(1:d)
      param$s <- array(data=NA, dim =c(d, d, 2))
      param$s[, , 1] <- array(data=14/45, dim=c(d,d))
      diag(param$s[,,1]) <- 4/9
      param$s[, , 2] <- diag(rep(4/9, times=d))
      param$p0 <- rep(1/2, times=2)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Asymmetric Bimodal']] <- param
      
      
      
      print('trimodal')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(3,d))
      param$m[1,] <- c(-1, rep(0, times=d-1))
      param$m[2,] <- c(1, rep(2/sqrt(3), times=d-1))
      param$m[3,] <- c(1, rep(-2/sqrt(3), times=d-1))
      param$s <- array(data=NA, dim =c(d, d, 3))
      param$s[, , 1] <- array(data=63/250, dim=c(d,d))#array(data=c((3/5)^2, 3*3*7/(5*5*10), 3*3*7/(5*5*10), (7/10)^2), dim=c(2,2))
      diag(param$s[,,1]) <- c(9/25, rep(49/100, times=d-1))
      param$s[, , 2] <- diag(c(9/25, rep(49/100, times=d-1)))#array(data=c((3/5)^2, 0, 0, (7/10)^2), dim=c(2,2))
      param$s[, , 3] <- diag(c(9/25, rep(49/100, times=d-1)))#array(data=c((3/5)^2, 0, 0, (7/10)^2), dim=c(2,2))
      param$p0 <- c(3/7, 3/7, 1/7)
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Trimodal']] <- param
      
      
      print('Fountain')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2^d + 2, d))
      param$m[1,] <- rep(0, times=d)
      param$m[2,] <- rep(0, times=d)
      for (no_d in 1:d){
        param$m[3:(2^d + 2), no_d] <- rep(x = c(-1,1), each=2^(no_d-1), length.out=2^d)
      }
      
      param$s <- array(data=NA, dim =c(d, d, 2^d + 2))
      param$s[, , 1] <- diag(nrow=d)
      for (no_g in 2:(2^d + 2)){
        param$s[, , no_g] <- diag(nrow=d)/16
      }
      
      param$p0 <- c(0.5, rep(0.5/(2^d +1), times=2^d+1))
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Fountain']] <- param
      
      
      
      
      print('Double Fountain')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(9, d))
      param$m[1,] <- c(-3/2, rep(0, times=d-1))
      param$m[2,] <- c(3/2, rep(0, times=d-1))
      param$m[6,] <- rep(0, times=d)
      for (i in -1:1){
        param$m[i+4,] <- c(i-3/2, rep(i, times=d-1))
        param$m[i+8,] <- c(i+3/2, rep(i, times=d-1))
      }
      
      param$s <- array(data=NA, dim =c(d, d, 9))
      for (i in 1:2){
        param$s[,,i] <- array(data=4/15, dim=c(d,d))
        diag(param$s[,,i]) <- rep(4/9, times=d)
      }
      
      for (i in 3:5){
        param$s[,,i] <- array(data=1/375, dim=c(d,d))
        diag(param$s[,,i]) <- rep(1/225, times=d)
        param$s[,,i+4] <- array(data=1/375, dim=c(d,d))
        diag(param$s[,,i+4]) <- rep(1/225, times=d)
      }
      
      param$s[,,6] <- array(data=3/45, dim=c(d,d))
      diag(param$s[,,6]) <- rep(1/9, times=d)
      
      param$p0 <- c(12/25, 12/25, rep(1/350, times=3), 8/350, rep(1/350, times=3))
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Double Fountain']] <- param
      
      
      
      
      print('Asymmetric Fountain')
      param <- list()
      param$nom_func <- 'MixGauss'
      param$m <- array(data=NA, dim=c(2^d + 2, d))
      param$m[1,] <- rep(0,times=d)
      param$m[2,] <- rep(0, times=d)
      for (no_d in 1:d){
        param$m[3:(2^d + 2), no_d] <- rep(x = c(-1,1), each=2^(no_d-1), length.out=2^d)
      }
      
      param$s <- array(data=NA, dim =c(d, d, 2^d+2))
      param$s[,,1] <- diag(nrow = d)
      mat1 <- array(data=-9/10, dim=c(d,d))
      diag(mat1) <- rep(1, times=d)
      param$s[,,2] <- mat1/16
      
      param$s[,, 3] <- mat1/4
      param$s[,, 4] <- diag(nrow = d)/8
      i0 <- 5
      i1 <- 6
      for (no_d in 1:(2^(d-1) -1)){
        param$s[,, i0] <- param$s[,, i0 - 2]/2
        param$s[,, i1] <- param$s[,, i1 - 2]/2
        i0 <- i0 + 2
        i1 <- i1 + 2
      }
      for (l in 1:(2^d+2)){
        if (!LaplacesDemon::is.positive.definite(param$s[, , l])){
          param$s[, , l] <- LaplacesDemon::as.positive.definite(param$s[, , l])
        }
      }
      
      param$p0 <- c(1/2, 3/40, 1/5, rep((1-0.5-3/40-1/5)/(2^d - 1), times=2^d - 1))
      param$supp <- matrix(data=NA, ncol=2, nrow=d)
      for (no_d in 1:d){
        param$supp[no_d, 1] <- -3*max(param$s[no_d, no_d, ]) + min(param$m[, no_d])
        param$supp[no_d, 2] <-  3*max(param$s[no_d, no_d, ]) + max(param$m[, no_d])
      }
      param_loi[['Asymmetric Fountain']] <- param
      rm(list='mat1')
      
      rm(list=c('i', 'i0', 'i1', 'l', 'no_d', 'no_g', 'param'))
      
      
      print('fin d\'initialisation des paramètres')
      
      
      
      # Pour chacune des densités test on génère N échantillons composés de
      # n observations de la loi correspondante
      
      print('Génération des échantillons...')
      nomdim <- list(NULL, 1:d,  1:N, nom_loi )

      echs <- array(data = NA, dim = c(n, d, N, length(nom_loi)), dimnames = nomdim)


      for (loi in 1:length(nom_loi)){
        print(nom_loi[loi])

        echs[, , , loi] <- replicate(N, generate_ech_mD(den_name = param_loi[[loi]]$nom_func, ech_size = n,
                                                   m = param_loi[[loi]]$m, s = param_loi[[loi]]$s,
                                                   p0 = param_loi[[loi]]$p0))

      }
      
      
      
      
      
      
      
      
      
      
      print('génération des échantillons d observations : ok')
      rm(list=c('loi', 'nomdim'))
      
      
    }
    
    print('génération des échantillons d observations : ok')
    # rm(list=c('no_loi', 'nomdim_echs'))
    # print(object.size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))), units="Mb")
    
    
    
    # A CALCULER DANS PRE_PCO
    # # intervalle d'estimation de la densite
    #
    # # cut = facteur multiplicatif de hmin pour déterminer l'interval sur
    # # lequel estimer la densité
    # # Plus précisément, l'interval (a,b) est défini par a=min(x_i)-cut*hmin et
    # # b=max(x_i)+cut*hmin
    # # Dans les fonctions de R, par défaut ce n'est pas hmin mais h.
    # # Le problème pour moi de prendre h c'est qu'ensuite les estimations ne
    # # sont pas toutes évaluées aux mêmes points et le calcul de la différence
    # # à la vraie densité oblige à recalculer la vraie
    # # densité autant de fois que d'estimations possibles
    # cut <- 3
    # h_cut <- 0.1
    #
    # nomdim_d <- outer(X='d', Y=1:d, FUN = function(x,y) paste(x,y))
    # nomdim_a = list(nomdim_d, c('min', 'max'), NULL, nom_loi)
    # a <- array(data=NA, dim=c(d, 2, N, length(nom_loi)), dimnames = nomdim_a)
    # a[,1,,] <- apply(X=densityech_all, MARGIN=2:4, FUN=min) - cut*h_cut
    # a[,2,,] <- apply(X=densityech_all, MARGIN=2:4, FUN=max) + cut*h_cut
    #
    # print('interval d\'estimation : ok')
    # # print(ls(envir = sys.frame(sys.nframe())))
    # # print(object.size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))), units="Mb")
    # rm(list=c('nomdim0', 'nomdim_a'))
    #
    
    
    
    print('sauvegarde des data...')
    
    
    #save.image(file = paste(rep_data, file_name, sep=''))
    #save.image(file = ech_file_name) # save.image enregistre le workspace de l'utilisateur pas de la fonction
    save(nom_loi, nom_loi_abr, ech_file_name, param_loi, echs, file=ech_file_name)
  }
  
  
}









compute_den_mD <- function(den_name, interval = NULL, x0 = 0, n = 100, m = 0, s = 1, p0 = 1){
  # renvoie la densité de proba correspondant à 'nom' sur l'interval 'interval'
  # calculée en n points équidistants sur l'interval.
  #
  # @param
  # nom      : chaine de caractere definissant la loi d'intérêt. Peut être :
  #            'Gauss', 'MixGauss', 'Unif', 'MixUnif', 'Cauchy', 'Exp'
  # interval : un vecteur de dimension 1 contenant 2 éléments, a et b, 
  #            correspondant aux extrêmités de l'interval
  # n        : nombre de bin où calculer la densité
  # m, s, p0 : les paramètres de la loi d'intérêt
  #            Dans le cas Gaussien : 2 scalaires, m pour la moyenne et s pour l'ecart-type
  #            Dans le cas Cauchy : 2 scalaires, m pour la position et s pour le facteur 
  #              d'échelle
  #            Dans le cas Uniforme : 2 scalaires, m pour le min et s pour le max
  #            Dans le cas Exponentielle : 1 scalaires, s pour le facteur d'échelle
  #            Dans le cas melange Gaussien : 3 vecteurs, m pour les moyennes
  #              s pour les ecart-types et p0 pour les proportions
  #            Dans le cas melange uniforme, 3 vecteurs : p0 qui correspond aux 
  #              proportions des differentes uniformes, m qui correspond aux minima 
  #              et s qui correspond aux maxima
  #  
  # @output
  # une liste contenant 2 vecteurs de dimension 1 et de taille n
  # le premier vecteur 'x' contient les n points équidistants de l'interval
  # le second vecteur 'y' contient les n valeurs de la densité calculée en x
  
  
  if (is.null(interval)){
    x <- x0
  }else{
    if (length(dim(m)) <= 1){
      #print(m)
      d <- length(m)
      #print(d)
    }else{
      d <- length(m[1,])
      #print(m)
      #print(d)
    }
    
    
    # chaque axe est decoupé en n points équidistants dans l'interval
    x <- array(data=NA, dim=c(n, d))
    for (no_d in 1:d){
      a <- interval['min', no_d]
      b <- interval['max', no_d]
      xtmp <- seq(from=a, to=b, length.out=n)
      x[, no_d] <- xtmp
      
    }
  }
  
  
  
  
  
  # les fonctions locales suivantes permettant de calculer la densité d'intérêt en x
  
  myGaussDen <- function(x, m, s, p0){
    
    # points où la densité est calculée (toutes les combinaisons possibles de coordonnées)
    xtmp <- array(data=NA, dim=c(n^d, d))
    
    for (no_d in 1:d){
      xtmp[, no_d] <- rep(x[, no_d], each=n^(no_d-1), len=n^d) 
    }
    
    U <- chol(s)
    # densité  gaussiennes calculée aux points x
    #y <- dmvn(x=xtmp, mu=m, Sigma=s)
    y <- LaplacesDemon::dmvnc(x=xtmp, mu=m, U=U)
    
    # reshape de la solution
    y <- array(data=y, dim=rep(n, times=d))
    
    return(y)
  }
  
  
  
  myCauchyDen <- function(x, m, s, p0){
    # points où la densité est calculée (toutes les combinaisons possibles de coordonnées)
    xtmp <- array(data=NA, dim=c(n^d, d))
    
    for (no_d in 1:d){
      xtmp[, no_d] <- rep(x[, no_d], each=n^(no_d-1), len=n^d) 
    }
    
    y <- LaplacesDemon::dmvc(x=xtmp, mu=m, S=s)
    
    # reshape de la solution
    y <- array(data=y, dim=rep(n, times=d))
    
    return(y)
  }
  
  
  myUnifCercleDen <- function(x, m, s, p0){
    d <- length(m)
    # points où la densité est calculée (toutes les combinaisons possibles de coordonnées)
    xtmp <- array(data=NA, dim=c(n^d, d))
    for (no_d in 1:d){
      xtmp[, no_d] <- rep(x[, no_d], each=n^(no_d-1), len=n^d)
    }
    x_max <- m + s
    x_min <- m - s
    #contrainte <- ((xtmp[, 1] - m[1])^2 + (xtmp[, 2] - m[2])^2 <= s^2)
    centre <- t(apply(xtmp, 1, function(x,m) x-m, m=m))
    contrainte <- rowSums((centre)^2) <= s^2 #((xtmp[, 1] - m[1])^2 + (xtmp[, 2] - m[2])^2 <= s^2)
    # fxy <- 1/((x_max[1]-x_min[1])*(x_max[2]-x_min[2]))*contrainte
    
    # fxy <- contrainte/prod(x_max-x_min)
    
    # calcul du volume d'une sphere en dimension d 
    S0 <- 2
    V0 <- 1
    Si <- S0
    Vi <- V0
    for (no_d in 1:d){
      Sip1 <- 2*pi*Vi
      Vip1 <- Si/no_d
      Si <- Sip1
      Vi <- Vip1
    }
    volume <- Vi*s^d
    fxy <- contrainte/volume
    return(array(data=fxy, dim=rep(n, times=d)))
  }
  
  myMixGaussDen <- function(x, m, s, p0){
    # densité melange de gaussiennes calculée aux points x
    d <- length(m[1,])
    
    
    # points où la densité est calculée (toutes les combinaisons possibles de coordonnées)
    # xtmp <- array(data=NA, dim=c(n^d, d))
    # for (no_d in 1:d){
    #   xtmp[, no_d] <- rep(x[, no_d], each=n^(no_d-1), len=n^d)
    # }
    # # 8 octets/double : 8*d*n^d
    # en dim 3, n=100 : 24Mo
    # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) # en dim 3, n=100, 24Mo
    # rm(xtmp)
    
    # quand la dimension augmente, xtmp devient beaucoup trop gros
    # au lieu de travailler directement sur x on utilise des tableaux d'indices de manière 
    # diviser la taille de xtmp par 2
    x_ind <- array(data=1:(n*d), dim=dim(x))
    xtmp_ind <- array(data=NA, dim=c(n^d, d))
    for (no_d in 1:d){
      xtmp_ind[, no_d] <- rep(x_ind[, no_d], each=n^(no_d-1), len=n^d)
    }
    # 4 octets/int : 4*d*n^d
    # en dim 3, n=100 : 12Mo
    #print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) 
    rm(x_ind)
    
    y <- array(data=0, dim=n^d)
    #print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))# en dim 3, n=100, 20Mo
    
    
    # on fait nb_k boucles sur n^d/nb_k points pour limiter la taille du vecteur passé en paramètre
    # de dmvn
    nb_k <- 10
    k <- n^d%/%nb_k
    for (no_k in 1:(nb_k-1)){
      xtmp <- x[as.vector(xtmp_ind[((no_k-1)*k + 1):(no_k*k),])]
      dim(xtmp) <- c(k, d)
      # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) 
      
      # on boucles sur les gaussiennes du melange
      for (no_G in 1:length(p0)){
        
        y[((no_k-1)*k + 1):(no_k*k)] <- y[((no_k-1)*k + 1):(no_k*k)] + p0[no_G]*LaplacesDemon::dmvn(x=xtmp, mu=m[no_G,], S=s[,,no_G])
        
      }
    }
    
    # no_k <- no_k + 1
    xtmp <- x[as.vector(xtmp_ind[((no_k*k) + 1):(n^d),])]
    dim(xtmp) <- c(n^d - no_k*k , d)
    #    print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
    for (no_G in 1:length(p0)){
      
      y[((no_k*k) + 1):(n^d)] <- y[((no_k*k) + 1):(n^d)] + p0[no_G]*LaplacesDemon::dmvn(x=xtmp, mu=m[no_G,], S=s[,,no_G])
      
    }
    rm(xtmp)
    rm(xtmp_ind)
    gc(verbose=TRUE)
    #   print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
    # on reshape y
    dim(y) <- rep(n, times=d)
    
    return(y)
  }
  
  
  
  myMixCauchyDen <- function(x, m, s, p0){
    # densité melange de gaussiennes calculée aux points x
    d <- length(m[1,])
    #print(d)
    # points où la densité est calculée (toutes les combinaisons possibles de coordonnées)
    # xtmp <- array(data=NA, dim=c(n^d, d))
    # for (no_d in 1:d){
    #   xtmp[, no_d] <- rep(x[, no_d], each=n^(no_d-1), len=n^d)
    # }
    
    
    # # 8 octets/double : 8*d*n^d
    # en dim 3, n=100 : 24Mo
    # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) # en dim 3, n=100, 24Mo
    # rm(xtmp)
    
    # quand la dimension augmente, xtmp devient beaucoup trop gros
    # au lieu de travailler directement sur x on utilise des tableaux d'indices de manière 
    # diviser la taille de xtmp par 2
    x_ind <- array(data=1:(n*d), dim=dim(x))
    xtmp_ind <- array(data=NA, dim=c(n^d, d))
    for (no_d in 1:d){
      xtmp_ind[, no_d] <- rep(x_ind[, no_d], each=n^(no_d-1), len=n^d)
    }
    # 4 octets/int : 4*d*n^d
    # en dim 3, n=100 : 12Mo
    # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) 
    rm(x_ind)
    
    y <- array(data=0, dim=n^d)
    # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))# en dim 3, n=100, 20Mo
    
    
    # on fait nb_k boucles sur n^d/nb_k points pour limiter la taille du vecteur passé en paramètre
    # de dmvn
    nb_k <- 10
    k <- n^d%/%nb_k
    for (no_k in 1:(nb_k-1)){
      xtmp <- x[as.vector(xtmp_ind[((no_k-1)*k + 1):(no_k*k),])]
      dim(xtmp) <- c(k, d)
      #  print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe()))))) 
      
      # on boucles sur les gaussiennes du melange
      for (no_G in 1:length(p0)){
        
        y[((no_k-1)*k + 1):(no_k*k)] <- y[((no_k-1)*k + 1):(no_k*k)] + p0[no_G]*LaplacesDemon::dmvc(x=xtmp, mu=m[no_G,], S=s[,,no_G])
        
      }
    }
    
    # no_k <- no_k + 1
    xtmp <- x[as.vector(xtmp_ind[((no_k*k) + 1):(n^d),])]
    dim(xtmp) <- c(n^d - no_k*k , d)
    #print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
    for (no_G in 1:length(p0)){
      
      y[((no_k*k) + 1):(n^d)] <- y[((no_k*k) + 1):(n^d)] + p0[no_G]*LaplacesDemon::dmvc(x=xtmp, mu=m[no_G,], S=s[,,no_G])
      
    }
    rm(xtmp)
    rm(xtmp_ind)
    # print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
    # on reshape y
    dim(y) <- rep(n, times=d)
    
    # y = array(data = 0, dim=rep(n, times=d))
    # for (no_G in 1:length(p0)){
    #   
    #   dtmp <- dmvc(x=xtmp, mu=m[no_G,], S=s[,,no_G])
    #   #print(dim(dtmp))
    #   # reshape de la solution
    #   y_tmp <- array(data=dtmp, dim=rep(n, times=d))
    #   y <- y + p0[no_G]*y_tmp
    # }
    return(y)
  }
  
  
  #   
  #   myCauchyDen <- function(x, m, s, p0){
  #     # densité de Cauchy calculée aux points x
  #     return(dcauchy(x = x, location = m, scale = s))
  #   }
  #   
  #   myUnifDen <- function(x, m, s, p0){
  #     # densité uniforme calculée aux points x
  #     return(dunif(x = x, min = m, max = s))
  #   }
  #   
  #   myMixUnifDen <- function(x, m, s, p0){
  #     # densité melange d'e gaussiennes'uniformes calculée aux points x
  #     y <- array(data = 0, dim = length(x) )
  #     for (i in 1:length(p0)){
  #       y <- y + p0[i]*dunif(x = x, min = m[i], max = s[i])
  #     }
  #     return(y)
  #   }
  #   
  #   myExpDen <- function(x, m, s, p0){
  #     # densité exponentielle calculée aux points x
  #     return(dexp(x = x, rate = s))
  #   }
  #   
  #   myLaplaceDen <- function(x, m, s, p0){
  #     # densité de laplace calculée aux points x
  #     return(exp(-abs(x-m)/s)/(2*s))
  #   }
  #   
  #   myMixLaplaceDen <- function(x, m, s, p0){
  #     # densité melange d'e gaussiennes'uniformes calculée aux points x
  #     y <- array(data = 0, dim = length(x) )
  #     for (i in 1:length(p0)){
  #       y <- y + p0[i]*exp(-abs(x-m[i])/s[i])/(2*s[i])
  #     }
  #     return(y)
  #   }
  
  # appel de la fonction correspondant à la densité d'intérêt
  nomfunc <- paste('my', den_name, 'Den', sep = '')
  y <- do.call(nomfunc, args = list(x = x, m = m, s = s, p0 = p0))
  
  
  # resultats (cad densité et points où la densité est calculée) mis sous forme de liste
  obj <- list()
  for (no_d in 1:d){
    nom_var <- paste('x_', no_d, sep = '')
    obj[[nom_var]] <- x[, no_d]
  }
  obj$z <- y
  
  return(obj)
  
  
}







#' Computes the density of all test laws in dimension d
#'
#' \code{gen_ech} Generate N n-samples of all test laws in dimension d
#'
#' @param n sample size
#' @param N Number of samples
#' @param d the dimension of the samples
#' @param nQMC le nombre de points où évaluer la densité
#' @param nh cardinal de H
#' @param nc cardinal de C
#' @param cmin min(C)
#' @param cmax max(C)
#'
#' @return calcul la densite en nQMC points issus d'une suite de Sobol sur
#' l'intervale [a, b] = [-5*max(H) + min(x_i); 5*max(H) + max(x_i)]
#' Les parametre nc cmin et cmax ne sont utiles que pour aller lire le fichier
#' pre_PCO
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @export
compute_N_den <- function(n, d, N, nQMC, maxH){
  # cette fonction calule la densite des lois test en nQMC points sur
  # l'intervalle [a, b]
  
  # on a besoin des paramètres qui ont servis à generer les échantillons
  # (en particulier n) pour pouvoir aller chercher le fichier contenant les
  # échantillons et ainsi en déduire l'intervalle de calcul [a,b]
  # a et b dépendent de H il faut donc aussi pre_PCO
  
  ech_file_name <- sprintf('ech_d%iN%in%i', d, N, n)
  load(file = ech_file_name)
  
  
  
  
  print('calcul de la vraie densite en cours...')
  
  
  
  if (d==1){
    
    a <- array(data=NA, dim = c(2, N, length(nom_loi)))
    dimnames(a) <- list(c('min', 'max'), 1:N, nom_loi)
    
    for (no_loi in 1:length(nom_loi)){
      print(nom_loi[no_loi])
      param <- param_loi[[no_loi]]
      for (no_ech in 1:N){
        f_file_name <- sprintf('f_d%iN%in%inloc%i_%s_%i', d, N, n, nQMC, nom_loi_abr[no_loi], no_ech)
        if (!file.exists(f_file_name)){
          x_i <- echs[, no_ech, no_loi]
          a['min', no_ech, no_loi] <- -5*maxH + min(x_i)
          a['max', no_ech, no_loi] <- 5*maxH + max(x_i)
          interval <- a[, no_ech, no_loi]
          f <- list()
          type <- param$nom_func
          m <- param$m
          s <- param$s
          p0 <- param$p0
          
          # compute_den('Gauss', interval = NULL, x0 = x, nb_eval_pts = nQMC,
          # m = 0, s = 1, p0 = 1)
          # compute_den va générer une suite de sobol dans interval et calculer
          # la densite en ces points
          f[[nom_loi[no_loi]]] <- compute_den_1D(den_name = type, m = m, s = s, p0 = p0,
                                                 interval = interval,
                                                 nb_eval_pts = nQMC)
          
          
          
          # print('sauvegarde des data...')
          save(list=c('f', 'f_file_name', 'interval'), file = f_file_name)
          # #print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
          # rm(f)
          # #gc(verbose=TRUE)
          # print(gc(verbose=TRUE))
        }# fin du if !file.exists
      }# fin de boucle sur no_ech
    }# fin de boucle sur no_loi
    
    
    
    
  }else{ # on est dans le cas d>1 : 
    
    
    
    a <- array(data=NA, dim = c(2, d, N, length(nom_loi)))
    dimnames(a) <- list(c('min', 'max'), 1:d, 1:N, nom_loi)
    
    a['min', , , ] <- apply(X=echs, MARGIN=2:4, FUN=min) - 0.3 # cut*h_cut = 0.3
    a['max', , , ] <- apply(X=echs, MARGIN=2:4, FUN=max) + 0.3
    
    for (no_loi in 1:length(nom_loi)){
      print(nom_loi[no_loi])
      param <- param_loi[[no_loi]]
      for (no_ech in 1:N){
        f_file_name <- sprintf('../f_d%iN%in%inloc%i_%s_%i', d, N, n, nQMC, nom_loi_abr[no_loi], no_ech)
        if (!file.exists(f_file_name)){
          x_i <- echs[, , no_ech, no_loi]
          
          interval <- a[, , no_ech, no_loi]
          f <- list()
          type <- param$nom_func
          m <- param$m
          s <- param$s
          p0 <- param$p0
          
          # compute_den('Gauss', interval = NULL, x0 = x, nb_eval_pts = nQMC,
          # m = 0, s = 1, p0 = 1)
          # compute_den va générer une suite de sobol dans interval et calculer
          # la densite en ces points
          f[[nom_loi[no_loi]]] <- compute_den_mD(den_name = type, m = m, s = s, p0 = p0,
                                                 interval = interval,
                                                 n = nQMC)
          
          
          
          # print('sauvegarde des data...')
          save(list=c('f', 'f_file_name', 'interval'), file = f_file_name)
          # #print(object_size(x=lapply(ls(), function(x) get(x, envir = sys.frame(sys.nframe())))))
          # rm(f)
          # #gc(verbose=TRUE)
          # print(gc(verbose=TRUE))
        }# fin du if !file.exists
      }# fin de boucle sur no_ech
    }# fin de boucle sur no_loi
    
    
    
    
  }
  
  print('calcul de la vraie densite : ok')
  
  
}



#' Generates the file that contains all h values
#'
#' \code{generate_H} Generates the file that contains all h values
#'
#' @param d dimension
#' @param nh number of h values
#' @param hmin_term the hmin value in 1D, and the diagonal term of hmin for multivariate case
#'
#' @return Generates the file that contains all h values. In multivariate case, the output is a list of vectors
#' that are later transformated (in the PCO_L2_loi_ech function) to obtain full h matrices.
#' Please note that hmin is not a parameter for epanechnikov and biweight kernel
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @export
generate_H(d, nh, hmin_term){
  if (d==1){
    Htmp <- sobol(n = nh - 2, dim = 1)
    # hmin <- eval(hmin_s)
    Htmp <- Htmp * (1 - hmin_term) + hmin_term
    H <- c(hmin_term, Htmp, 1)
    H <- sort(H)
    # print(sprintf("hmin_s = %s", hmin_s))
    # print(sprintf("hmin = %f", hmin))
    
    save(list=c('H'), file = "H")
    
  }else{
    
      H <- list()
      
      
      # hmin_term correspond à 1 terme de la diagonale qui sera répété d fois
      
      # suite de sobol ramenée à l'interval [1/(||K||_inf/n^(1/d)); 1]^d
      # avec en plus une matrice diagonale de determiant égal à ||K||_inf/n
      
      normesupK <- 1 / sqrt(2*pi) # gaussian kernel
      
      Htmp <- randtoolbox::sobol(n = nh - 1, dim = d) # nh-1 vecteurs de dimension d
      # cst_diag <- sqrt(2 * pi) * n ^ (1 / d)
      # hmin_diag_term <- eval(hmin_s)
      # hmin_diag_term <- hmin_term
      Htmp <- Htmp * (1 - hmin_term) + hmin_term
      Htmp <- rbind(rep(hmin_term, times = d), Htmp)
      # H_vect <- list()
      for (i in 1:nh){
        # H[ , , i] <- diag(Htmp[i, ])
        # H[[i]] <- diag(Htmp[i, ]) # potentiellement peut etre remplace par une liste de vecteur au lieu d'une liste de matrices
        # H_vect[[i]] <- Htmp[i, ]
        H[[i]] <- Htmp[i, ]
      }
      
      # produit des éléments de la diagonale de h
      hi_prod <- apply(Htmp, 1, prod)
      hi_prod_indexsorted <- sort(hi_prod, index.return=TRUE)$ix
      
      hi_prod_sorted <- sort(hi_prod)
      
      # hmin <- H[ , , hi_prod_indexsorted[1]]
      hmin <- H[[hi_prod_indexsorted[1]]]
      
      # H est trié en fonction du produit des éléments diagonaux de chaque matrice (des matrices diagonales)
      H <- H[hi_prod_indexsorted]
      # H_vect <- H_vect[hi_prod_indexsorted]
      Htmp_sorted <- Htmp[hi_prod_indexsorted]
      
      # print(paste('min(H) = ', min(hi_prod)))
      # print(paste('||K||_inf/n = 1/hmin_diag_term^d = ', 1/hmin_diag_term^d))
      # #print(paste('det(hmin) = ', det(hmin)))
      # print(paste('1/n = ', 1/n))
      
      
      save(list=c('H'), file = "H")
      
    
    
  }
}









#' Calcul le critere PCO_test_hmin
#'
#' \code{gen_ech} Calcul du critere PCO_test_hmin
#'
#' @param n sample size
#' @param d the dimension of the samples
#' @param nQMC
#' @param nh cardinal de H
#' @param nc cardinal de C
#' @param cmin min(C)
#' @param cmax max(C)
#'
#' @return ne retourne rien mais génere un fichier contenant C, H et pen
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @export
PCO_L2_loi_ech <- function(n, d, N, H, no_loi, no_ech, K_name = "gaussian", diagonal=TRUE){
  
  
  
  ech_file_name <- sprintf('ech_d%iN%in%i', d, N, n)
  load(file = ech_file_name)
  
  if (K_name == "gaussian"){
    K_letter <- "g"
  }
  if (K_name == "biweight"){
    K_letter <- "b"
  }
  if (K_name == "epanechnikov"){
    K_letter <- "e"
  }
  # pre_PCO_file_name <- sprintf('pre_PCO_L2_d%in%inh%inc%icm%icM%i_K%s', d, n, nh,
  #                              K_letter)
  # if (file.exists(pre_PCO_file_name)){
  #   load(file = pre_PCO_file_name)
  # }else{
  #   pre_PCO_L2(n, d, nh, K_name)
  #   load(file = pre_PCO_file_name)
  # }
  
  # print(H)
  
  
  PCO_file_name <- sprintf('PCO_L2_K%sd%iN%in%i_%s_%i', K_letter, d, N, n,
                           nom_loi_abr[no_loi], no_ech)
  if (!file.exists(PCO_file_name)){
    density_hat_PCO <- list()
    density_hat_PCO[[nom_loi[no_loi]]] <- list()
    
    # la penalite a ete calculee dans pre_PCO
    
    if (d==1){
      # il faut donc maintenant calculer la perte
      x_i <- echs[, no_ech, no_loi]
      # print(x_i)
      if (K_name=="gaussian"){
        # u <- outer(x_i, x_i, function(x,y) x-y)
        
        # perte <- loss_L2_GK_1d_exact_e3(x_i, H)
        crit <- crit_GK_exact_1D(x_i, H)
        
        hmin <- min(H)
        # ptm0 <- proc.time()
        # u <- outer(x_i, x_i, function(x,y) x-y)
        # crit0 <- crit_L2_GK_1d_exact_e(x_i, u, H, n)
        # t0 <- proc.time() - ptm0
        #
        # ptm0 <- proc.time()
        # crit1 <- crit_L2_GK_1d_exact_e3(x_i, H, n)
        # t1 <- proc.time() - ptm0
        
        
        no_h_opts <- which.min(crit)
        
        density_hat_PCO[[nom_loi[no_loi]]]$bw <- H[no_h_opts]
        density_hat_PCO[[nom_loi[no_loi]]]$no_bw <- no_h_opts
        
      }else{
        if (K_name=="epanechnikov"){
          
          
          # u <- outer(x_i, x_i, function(x,y) x-y)
          # perte <- loss_L2_EK_1d_exact(x_i, u, H, n)
          #
          #
          # ptm0 <- proc.time()
          # u <- outer(x_i, x_i, function(x,y) x-y)
          # perte0 <- loss_L2_EK_1d_exact(x_i, u, H, n)
          # t0 <- proc.time() - ptm0
          # #
          # ptm0 <- proc.time()
          # perte1 <- loss_L2_EK_1d_exact_e1(x_i, H, n)
          # t1 <- proc.time() - ptm0
          #
          # ptm0 <- proc.time()
          # perte2 <- loss_L2_EK_1d_exact_e2(x_i, H, n)
          # t2 <- proc.time() - ptm0
          
          # ptm0 <- proc.time()
          # perte3 <- loss_L2_EK_1d_exact_e3(x_i, H, n)
          crit <- crit_EK_exact_1D(x_i, H)
          # t3 <- proc.time() - ptm0
          
          
          # ptm0 <- proc.time()
          # u <- outer(x_i, x_i, function(x,y) x-y)
          # crit0 <- crit_L2_EK_1d_exact(x_i, u, H, n)
          # t0 <- proc.time() - ptm0
          #
          # ptm0 <- proc.time()
          # crit1 <- crit_L2_EK_1d_exact_e3(x_i, H, n)
          # t1 <- proc.time() - ptm0
          
          
        }else{
          if (K_name=="biweight"){
            # u <- outer(x_i, x_i, function(x,y) x-y)
            # perte <- loss_L2_BK_1d_exact(x_i, u, H, n)
            #
            #
            # ptm0 <- proc.time()
            # u <- outer(x_i, x_i, function(x,y) x-y)
            # perte0 <- loss_L2_BK_1d_exact(x_i, u, H, n)
            # t0 <- proc.time() - ptm0
            #
            # ptm0 <- proc.time()
            # perte1 <- loss_L2_BK_1d_exact_e1(x_i, H, n)
            # t1 <- proc.time() - ptm0
            #
            # ptm0 <- proc.time()
            # perte2 <- loss_L2_BK_1d_exact_e2(x_i, H, n)
            crit <- crit_BK_exact_1D(x_i, H)
            # t2 <- proc.time() - ptm0
            
            
            # ptm0 <- proc.time()
            # u <- outer(x_i, x_i, function(x,y) x-y)
            # crit0 <- crit_L2_BK_1d_exact(x_i, u, H, n)
            # t0 <- proc.time() - ptm0
            #
            # ptm0 <- proc.time()
            # crit1 <- crit_L2_BK_1d_exact_e3(x_i, H, n)
            # t1 <- proc.time() - ptm0
            
            
          }
        }
        
      }
    }else{
      x_i <- echs[, , no_ech, no_loi]
      if (diagonal){
        # H <- list()
        # 
        # 
        # 
        # # suite de sobol ramenée à l'interval [1/(||K||_inf/n^(1/d)); 1]^d
        # # avec en plus une matrice diagonale de determiant égal à ||K||_inf/n
        # 
        # Htmp <- randtoolbox::sobol(n = nh - 1, dim = d) # nh vecteurs de dimension d
        # ndethmin <- sqrt(2 * pi) * n ^ (1 / d)
        # Htmp <- Htmp * (1 - 1 / (ndethmin)) + 1 / (ndethmin)
        # Htmp <- rbind(rep(1/(ndethmin), times = d), Htmp)
        # H_vect <- list()
        # for (i in 1:nh){
        #   # H[ , , i] <- diag(Htmp[i, ])
        #   H[[i]] <- diag(Htmp[i, ]) # potentiellement peut etre remplace par une liste de vecteur au lieu d'une liste de matrices
        #   H_vect[[i]] <- Htmp[i, ]
        # }
        # 
        # # rm(i)
        # 
        # # produit des éléments de la diagonale de h
        # hi_prod <- apply(Htmp, 1, prod)
        # hi_prod_indexsorted <- sort(hi_prod, index.return=TRUE)$ix
        # 
        # hi_prod_sorted <- sort(hi_prod)
        # 
        # # hmin <- H[ , , hi_prod_indexsorted[1]]
        # hmin <- H[[hi_prod_indexsorted[1]]]
        # 
        # # H est trié en fonction du produit des éléments diagonaux de ses éléments (des matrices diagonales)
        # H <- H[hi_prod_indexsorted]
        # H_vect <- H_vect[hi_prod_indexsorted]
        # Htmp_sorted <- Htmp[hi_prod_indexsorted]
        # 
        # print(paste('min(H) = ', min(hi_prod)))
        # print(paste('||K||_inf/n = ', 1/ndethmin^d))
        # print(paste('det(hmin) = ', det(hmin)))
        # print(paste('1/n = ', 1/n))
        # 
        
        # crit_GK_exact_mD_diag(Eigen::MatrixXd x_i, List H, Eigen::VectorXd hmin)
        hmin <- H[[1]]
        crit <- crit_GK_exact_mD_diag(x_i, H, hmin)
        
        
        no_h_opts <- which.min(crit)
        
        density_hat_PCO[[nom_loi[no_loi]]]$bw <- H[[no_h_opts]]
        density_hat_PCO[[nom_loi[no_loi]]]$no_bw <- no_h_opts
        
        
      }else{
        
        print("diagonal = FALSE")
        # H est une liste contenant des vecteurs pareil à ceux du cas diagonal
        # Pour avoir des matrices pleines on applique une transformation
        # H_new correspond à P %*% H %*% Pinv où P est la matrice de passage 
        # obtenue lors de la diagonalisation de la matrice de covariance des xi
        H_new <- list()
        # 
        # hi_tmp <- randtoolbox::sobol(n = nh - 1, dim = d)
        # ndethmin <- sqrt(2 * pi) * n ^ (1 / d)
        # hi_tmp <- hi_tmp * (1 - 1 / (ndethmin)) + 1 / (ndethmin)
        # 
        S <- array(data=NA, dim=c(d, d))
        S <- cov(x_i)
        # Diagonalisation de S
        S.eig <- eigen(S)
        P <- S.eig$vectors
        P_inv <- solve(S.eig$vectors)

        nh <- length(H)
        # H_new[[1]] <- P %*% daig(H[[1]]) %*% P_inv
        for (i in 1:nh){
          H_new[[i]] <- P %*% diag(H[[i]]) %*% P_inv
        }
        # H est trie en fonction du produit des elements diagonaux 
        # (càd le determinant car H ne contient que des matrices doagonales) 
        # la transformation P %*% D %*% P_inv ne modifie pas le determinant. 
        # Donc il n'est pas nécessaire de retrier H_new en fonction du determinant
        
        # # 
        # # # calcul du determinant de chaque matrice
        # # # hi_det est un tableau de dimension h_size
        # hi_det <- unlist(lapply(X = H_new, det), use.names=FALSE)
        # # # tri de H en fonction du determinant
        # # # hi_det_indexsorted est un tableau de dimension h_size
        # hi_det_indexsorted <- sort(hi_det, index.return=TRUE)$ix
        # #
        # # hmin <- H_new[[hi_det_indexsorted[1]]]
        # #
        # H_new <- H_new[hi_det_indexsorted]
        hmin <- H_new[[1]]
        #
        # # Calcul des inverses de H
        # # H_inv <- list()
        # # H_inv <- lapply(X = H, solve)
        #
        #
        #
        # print(paste('min(H) = ', min(hi_det)))
        # # print(paste('||K||_inf/n = ', 1/ndethmin^d))
        # print(paste('det(hmin) = ', det(hmin)))
        # print(paste('1/n = ', 1/n))

        # crit_GK_exact_mD_full(Eigen::MatrixXd x_i, Eigen::MatrixXd S, List H, Eigen::VectorXd hmin_diag)
        # hmin <- H[[1]]
        crit <- crit_GK_exact_mD_full(x_i, H_new, hmin)
        
        
        no_h_opts <- which.min(crit)
        
        density_hat_PCO[[nom_loi[no_loi]]]$bw <- H_new[[no_h_opts]]
        density_hat_PCO[[nom_loi[no_loi]]]$no_bw <- no_h_opts
        
      }
      
    }
    
    
    
    
    
    # crit <- pen + perte
    
    
    # 
    density_hat_PCO[[nom_loi[no_loi]]]$hmin <- hmin
    
    density_hat_PCO[[nom_loi[no_loi]]]$crit <- crit
    
    # density_hat_PCO[[nom_loi[no_loi]]]$perte <- perte
    
    
    save(list=c('density_hat_PCO', 'PCO_file_name'), file = PCO_file_name)
    
  }
  
}


















#' Calcul le critere PCO_test_hmin
#'
#' \code{gen_ech} Calcul du critere PCO_test_hmin
#'
#' @param n sample size
#' @param d the dimension of the samples
#' @param nQMC
#' @param nh cardinal de H
#' @param nc cardinal de C
#' @param cmin min(C)
#' @param cmax max(C)
#'
#' @return ne retourne rien mais génere un fichier contenant C, H et pen
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @export
compute_risks_Lp <- function(p, n, d, N, nQMC, nh, K_name = "gaussian"){
# compute_risks_Lp <- function(p, n, d, N, nQMC, nh, nc, cmin, cmax, K_name = "gaussian"){
  
  
  print('calcul du risque')
  nc <- 1
  
  ech_file_name <- sprintf('ech_d%iN%in%i', d, N, n)
  load(ech_file_name)
  
  if (K_name == "gaussian"){
    K_letter <- "g"
  }
  if (K_name == "biweight"){
    K_letter <- "b"
  }
  if (K_name == "epanechnikov"){
    K_letter <- "e"
  }
  # pre_PCO_file_name <- sprintf('pre_PCO_L%i_d%in%inh%inc%icm%icM%i_K%s', p, d, n, nh,
  #                              nc, cmin, cmax, K_letter)
  # load(pre_PCO_file_name)
  
  
  # risk_file_name <- sprintf('risksL%i_n%id%iN%inh%inc%icmin%icmax%i', p, n, d, N, nh,
                            # nc, cmin, cmax)
  risk_file_name <- sprintf('risksL%i_n%id%iN%inh%i', p, n, d, N, nh)
  if (file.exists(risk_file_name)){
    load(risk_file_name)
  }else{
    risks <- array(data=NA, dim=c(nc, N, length(nom_loi)))
    dimnames(risks) <- list(1:nc, 1:N, nom_loi)
    dLp_fhmin <- array(data=NA, dim=c(nc, N, length(nom_loi)))
    dimnames(dLp_fhmin) <- list(1:nc, 1:N, nom_loi)
    no_c_opt <- array(data=NA, dim = c(N, length(nom_loi)))
    dimnames(no_c_opt) <- list(1:N, nom_loi)
    
  }
  
  
  
  
  for (no_loi in 1:length(nom_loi)){
    print(nom_loi[no_loi])
    param <- param_loi[[no_loi]]
    for (no_ech in 1:N){
      
      # on charge la vraie densite
      f_file_name <- sprintf('../f_d%iN%in%inloc%i_%s_%i', d, N, n, nQMC, nom_loi_abr[no_loi], no_ech)
      if (file.exists(f_file_name)){
        load(f_file_name)
        # print('f_file loaded')
        # on charge les h optimaux de PCO_test_hmin
        if (is.na(risks[, no_ech, no_loi])){
          # print("is.na(risk) = TRUE")
          PCO_file_name <- sprintf('PCO_L%i_K%sd%iN%in%i_%s_%i', p, K_letter, d, N, n, nom_loi_abr[no_loi], no_ech)
          if (file.exists(PCO_file_name)){
            # print("PCO_file exists")
            load(PCO_file_name)
            # print("PCO_file loaded")
            hopt_c <- density_hat_PCO[[nom_loi[no_loi]]]$bw
            no_hopt_c <- density_hat_PCO[[nom_loi[no_loi]]]$no_bw
            
            
            if (d==1){
              x_i <- echs[, no_ech, no_loi]
              a <- interval['min']
              b <- interval['max']
              # risk <- risque_Lp(p, f[[1]]$y, x_i, C, f[[1]]$x, hopt_c, a, b) # risque avec noyau gaussien
              risk <- risque_Lp(p, f[[1]]$y, x_i, f[[1]]$x, hopt_c, a, b) # risque avec noyau gaussien
              # no_c_opt[no_ech, no_loi] <- which.min(risk)
              risks[, no_ech, no_loi] <- risk
              # dLp_fhmin[, no_ech, no_loi] <- density_hat_PCO[[nom_loi[no_loi]]]$perte[no_hopt_c]
              
            }else{
              # print("d != 1")
              
              # x_test = array(data = NA, dim=c(2,2)) #dim = (nloc=3, d=2)
              # oi_test = array(data = 1, dim=c(4,2)) #dim = (n=5, d=2)
              # 
              # 
              # d <- 4
              # m <- rep(0, times=d)
              # sig_test <- array(data=9/10, dim =c(d, d))
              # diag(sig_test) <- 1
              # 
              # n <- 1000
              # nloc <- 40
              
              # o_i <- generate_ech_mD('Gauss', ech_size = n, m = m, s=sig_test)
              x_i <- echs[, , no_ech, no_loi]
              # print("x_i ok")
              a <- interval['min', ]
              b <- interval['max', ]
              # print("a et b ok")
              # a_min <- apply(X=x_i, MARGIN=2:4, FUN=min) - 0.3 # cut*h_cut = 0.3
              # a_max <- apply(X=x_i, MARGIN=2:4, FUN=max) + 0.3 # cut*h_cut = 0.3
              # 
              # # a['min', , no_ech, no_loi] <- apply(X=echs, MARGIN=2:4, FUN=min) - 0.3 # cut*h_cut = 0.3
              # # a['max', , no_ech, no_loi] <- apply(X=echs, MARGIN=2:4, FUN=max) + 0.3
              # interval <- array(data=NA, dim=c(2,d), dimnames = list(c('min', 'max'), NULL))
              # interval['min', ] <- a_min
              # interval['max', ] <- a_max
              # f_test <- compute_den_mD(den_name = 'Gauss', interval = interval, n = nloc, m = m, s = sig_test)
              
              vec_f <- f[[1]]$z
              # print("vec_f ok")
              # dim(vec_f) <- nloc^d#nloc^d
              dim(vec_f) <- nQMC^d#nloc^d
              # print("dim vec_f ok")
              
              # risque_Lp_mD(int p, Eigen::VectorXd f, Eigen::MatrixXd o_i, 
              #              Eigen::MatrixXd x, double h, double a, double b)
              
              # H_test <- diag(c(0.1,0.1,0.1))
              # H_test <- sig_test
              # # H_test <- diag(c(0.1,0.1))
              # H2_test <- H_test %*% H_test
              
              # nb_bench <- 100
              # ptm0 <- proc.time()
              # for (i in 1:nb_bench){
              #   r_1 <- risque_Lp_mD_1(p = 2, f = vec_f, o_i = o_i, x = x_test, 
              #                 h = H_test, 
              #                 a = a_min, b = a_max)
              # }
              # t_r1_mD <- (proc.time() - ptm0)/nb_bench
              
              # ptm0 <- proc.time()
              # for (i in 1:nb_bench){
              h_tmp <- unlist(hopt_c)
              if (is.null(dim(h_tmp))){
                # h_tmp est un vecteur
                h_pco <- diag(h_tmp)
              }else{
                h_pco <- h_tmp
              }
              r_2 <- risque_Lp_mD_3(p = 2, f = vec_f, o_i = x_i, nQMC = nQMC,
                                    h = h_pco, a = a, b = b)
              
              
              risks[, no_ech, no_loi] <- r_2
              # }
              # t_r2_mD <- (proc.time() - ptm0)/nb_bench
              
              # somme = 0
              # for (i in 1:n){
              #   # a_min_test <- matrix(a_min, ncol = length(a_min))
              #   # tmp <- backsolve(chol(H2_test), t(a_min_test) - o_i[i, ], transpose = TRUE)
              #   # print(tmp)
              #   # rss <- colSums(tmp^2)
              #   # logretval <- -sum(log(diag(chol(H2_test)))) - 0.5 * d * log(2 * 
              #   #                                                     pi) - 0.5 * rss
              #   # print(exp(logretval))
              #   print(mvtnorm::dmvnorm(x = a_min, mean = o_i[i, ], sigma = H2_test))
              # }
              # # mvtnorm::dmvnorm(x = a_min, mean = o_i[1,], sigma = H_test)
              # test_kde <- ks::kde(x=o_i, H=H2_test, xmin = a_min, xmax =a_max, binned=FALSE, gridsize = c(2,2,2))
              
            }#fin du else du if d==1
            
          }#fin du if PCO_file exists
        }# fin du if f_file exists
        
      }# fin du if is.na(risk)
      
    }#fin de boucle sur ech
  }#fin de boucle sur loi
  
  
  save(risks, file = risk_file_name)
  
  
}



#' Calcul le critere PCO_test_hmin
#'
#' \code{gen_ech} Calcul du critere PCO_test_hmin
#'
#' @param n sample size
#' @param d the dimension of the samples
#' @param nQMC
#' @param nh cardinal de H
#' @param nc cardinal de C
#' @param cmin min(C)
#' @param cmax max(C)
#'
#' @return ne retourne rien mais génere un fichier contenant C, H et pen
#'
#' @useDynLib PCOtesthmin
#' @importFrom Rcpp sourceCpp
#' @export
compute_risks_Lp_loi_ech <- function(p, n, d, nQMC, nh, no_ech, no_loi){
  print('calcul du risque')
  nc <- 1
  
  ech_file_name <- sprintf('ech_d%iN%in%i', d, N, n)
  load(ech_file_name)
  
  # if (K_name == "gaussian"){
  K_letter <- "g"
  # }
  # if (K_name == "biweight"){
  #   K_letter <- "b"
  # }
  # if (K_name == "epanechnikov"){
  #   K_letter <- "e"
  # }
  # pre_PCO_file_name <- sprintf('pre_PCO_L%i_d%in%inh%inc%icm%icM%i_K%s', p, d, n, nh,
  #                              nc, cmin, cmax, K_letter)
  # load(pre_PCO_file_name)
  
  
  # risk_file_name <- sprintf('risksL%i_n%id%iN%inh%inc%icmin%icmax%i', p, n, d, N, nh,
  # nc, cmin, cmax)
  # risk_file_name <- sprintf('risksL%i_n%id%iN%inh%i', p, n, d, N, nh)
  # if (file.exists(risk_file_name)){
  #   load(risk_file_name)
  # }else{
  #   risks <- array(data=NA, dim=c(nc, N, length(nom_loi)))
  #   dimnames(risks) <- list(1:nc, 1:N, nom_loi)
  #   dLp_fhmin <- array(data=NA, dim=c(nc, N, length(nom_loi)))
  #   dimnames(dLp_fhmin) <- list(1:nc, 1:N, nom_loi)
  #   no_c_opt <- array(data=NA, dim = c(N, length(nom_loi)))
  #   dimnames(no_c_opt) <- list(1:N, nom_loi)
  #   
  # }
  
  
  
  
  # for (no_loi in 1:length(nom_loi)){
  print(nom_loi[no_loi])
  param <- param_loi[[no_loi]]
  # for (no_ech in 1:N){
  
  # on charge la vraie densite
  f_file_name <- sprintf('../f_d%iN%in%inloc%i_%s_%i', d, N, n, nQMC, nom_loi_abr[no_loi], no_ech)
  if (file.exists(f_file_name)){
    
    # print('f_file loaded')
    # on charge les h optimaux de PCO_test_hmin
    # if (is.na(risks[, no_ech, no_loi])){
    # print("is.na(risk) = TRUE")
    PCO_file_name <- sprintf('PCO_L%i_K%sd%iN%in%i_%s_%i', p, K_letter, d, N, n, nom_loi_abr[no_loi], no_ech)
    if (file.exists(PCO_file_name)){
      # print("PCO_file exists")
      
      load(PCO_file_name)
      if (is.null(density_hat_PCO[[nom_loi[no_loi]]]$risk)){
        load(f_file_name)
        # print("PCO_file loaded")
        hopt_c <- density_hat_PCO[[nom_loi[no_loi]]]$bw
        no_hopt_c <- density_hat_PCO[[nom_loi[no_loi]]]$no_bw
        
        
        if (d==1){
          x_i <- echs[, no_ech, no_loi]
          a <- interval['min']
          b <- interval['max']
          # risk <- risque_Lp(p, f[[1]]$y, x_i, C, f[[1]]$x, hopt_c, a, b) # risque avec noyau gaussien
          risk <- risque_Lp(p, f[[1]]$y, x_i, f[[1]]$x, hopt_c, a, b) # risque avec noyau gaussien
          # no_c_opt[no_ech, no_loi] <- which.min(risk)
          risks[, no_ech, no_loi] <- risk
          # dLp_fhmin[, no_ech, no_loi] <- density_hat_PCO[[nom_loi[no_loi]]]$perte[no_hopt_c]
          
        }else{
          # print("d != 1")
          
          # x_test = array(data = NA, dim=c(2,2)) #dim = (nloc=3, d=2)
          # oi_test = array(data = 1, dim=c(4,2)) #dim = (n=5, d=2)
          # 
          # 
          # d <- 4
          # m <- rep(0, times=d)
          # sig_test <- array(data=9/10, dim =c(d, d))
          # diag(sig_test) <- 1
          # 
          # n_test <- 1000
          # nloc_test <- 40
          # 
          # o_i <- generate_ech_mD('Gauss', ech_size = n_test, m = m, s=sig_test)
          # 
          
          x_i <- echs[, , no_ech, no_loi]
          # print("x_i ok")
          a <- interval['min', ]
          b <- interval['max', ]
          # print("a et b ok")
          # a_min <- apply(X=x_i, MARGIN=2:4, FUN=min) - 0.3 # cut*h_cut = 0.3
          # a_max <- apply(X=x_i, MARGIN=2:4, FUN=max) + 0.3 # cut*h_cut = 0.3
          # 
          # # a['min', , no_ech, no_loi] <- apply(X=echs, MARGIN=2:4, FUN=min) - 0.3 # cut*h_cut = 0.3
          # # a['max', , no_ech, no_loi] <- apply(X=echs, MARGIN=2:4, FUN=max) + 0.3
          # interval <- array(data=NA, dim=c(2,d), dimnames = list(c('min', 'max'), NULL))
          # interval['min', ] <- a_min
          # interval['max', ] <- a_max
          # f_test <- compute_den_mD(den_name = 'Gauss', interval = interval, n = nloc_test, m = m, s = sig_test)
          # vec_f_test <- f_test$z
          # dim(vec_f_test) <- nloc_test^d#nloc^d
          
          
          # ptm <- proc.time()
          # r_2 <- risque_Lp_mD_2(p = 2, f = vec_f_test, o_i = o_i, nQMC = nloc_test,
          #                       h = h_pco, a = a, b = b)
          # t2 <- proc.time() - ptm
          # 
          # ptm <- proc.time()
          # r_3 <- risque_Lp_mD_3(p = 2, f = vec_f_test, o_i = o_i, nQMC = nloc_test,
          #                       h = h_pco, a = a, b = b)
          # t3 <- proc.time() - ptm
          
          
          
          
          
          vec_f <- f[[1]]$z
          # print("vec_f ok")
          # dim(vec_f) <- nloc^d#nloc^d
          dim(vec_f) <- nQMC^d#nloc^d
          # print("dim vec_f ok")
          
          # risque_Lp_mD(int p, Eigen::VectorXd f, Eigen::MatrixXd o_i, 
          #              Eigen::MatrixXd x, double h, double a, double b)
          
          # H_test <- diag(c(0.1,0.1,0.1))
          # H_test <- sig_test
          # # H_test <- diag(c(0.1,0.1))
          # H2_test <- H_test %*% H_test
          
          # nb_bench <- 100
          # ptm0 <- proc.time()
          # for (i in 1:nb_bench){
          #   r_1 <- risque_Lp_mD_1(p = 2, f = vec_f, o_i = o_i, x = x_test, 
          #                 h = H_test, 
          #                 a = a_min, b = a_max)
          # }
          # t_r1_mD <- (proc.time() - ptm0)/nb_bench
          
          # ptm0 <- proc.time()
          # for (i in 1:nb_bench){
          # h_tmp <- unlist(hopt_c)
          if (is.null(dim(hopt_c))){
            # h_tmp est un vecteur
            h_pco <- diag(hopt_c)
            print("diagonal = TRUE")
            
          }else{
            h_pco <- hopt_c
            print("diagonal = FALSE")
          }
          print(h_pco)
          # h_pco <- hopt_c
          r_3 <- risque_Lp_mD_3(p = 2, f = vec_f, o_i = x_i, nQMC = nQMC,
                                h = h_pco, a = a, b = b)
          
          
          density_hat_PCO[[nom_loi[no_loi]]]$risk <- r_3
          # risks[, no_ech, no_loi] <- r_2
          # }
          # t_r2_mD <- (proc.time() - ptm0)/nb_bench
          
          # somme = 0
          # for (i in 1:n){
          #   # a_min_test <- matrix(a_min, ncol = length(a_min))
          #   # tmp <- backsolve(chol(H2_test), t(a_min_test) - o_i[i, ], transpose = TRUE)
          #   # print(tmp)
          #   # rss <- colSums(tmp^2)
          #   # logretval <- -sum(log(diag(chol(H2_test)))) - 0.5 * d * log(2 * 
          #   #                                                     pi) - 0.5 * rss
          #   # print(exp(logretval))
          #   print(mvtnorm::dmvnorm(x = a_min, mean = o_i[i, ], sigma = H2_test))
          # }
          # # mvtnorm::dmvnorm(x = a_min, mean = o_i[1,], sigma = H_test)
          # test_kde <- ks::kde(x=o_i, H=H2_test, xmin = a_min, xmax =a_max, binned=FALSE, gridsize = c(2,2,2))
          
          
          save(list=c('density_hat_PCO', 'PCO_file_name'), file = PCO_file_name)
          
        }#fin du else du if d==1
      }# fin du if is.null(density_hat_PCO$risk)
    }#fin du if PCO_file exists
    
  }# fin du if f_file exists
  
  # }# fin du if is.na(risk)
  
  # }#fin de boucle sur ech
  # }#fin de boucle sur loi
  
  
  # save(risks, file = risk_file_name)
  
  
  
  
  
}





