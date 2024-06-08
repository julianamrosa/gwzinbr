#' @title Golden Section Search
#'
#' @description Runs a Golden Section Search (GSS) algorithm for determining the optimum bandwidth for the geographically weighted zero inflated negative binomial regression and other spatial regression models.
#'
#' @param data name of the dataset.
#' @param formula regression model formula as in \code{lm}.
#' @param xvarinf name of the covariates for the zero inflated part of the model, default value is \code{NULL}.
#' @param weight name of the variable containing the sample weights, default value is \code{NULL}.
#' @param lat name of the variable containing the latitudes in the dataset.
#' @param long name of the variable containing the longitudes in the dataset.
#' @param globalmin logical value indicating whether to find a global minimum in the optimization process, default value is \code{TRUE}.
#' @param method indicates the method to be used for the bandwidth calculation (\code{adaptive_bsq} or \code{fixed_g}).
#' @param model indicates the model to be used for the regression (\code{zinb}, \code{zip}, \code{negbin}, \code{poisson}), default value is\code{"zinb"}.
#' @param bandwidth indicates the criterion to be used for the bandwidth calculation (\code{cv}, \code{aic}), default value is \code{"cv"}.
#' @param offset name of the variable containing the offset values, if null then is set to a vector of zeros, default value is \code{NULL}.
#' @param force logical value indicating whether to force the indicated model even if it is not the best fit for the data, default value is \code{FALSE}.
#' @param maxg integer indicating the maximum number of iterations for the zero inflated part of the model, default value is \code{100}.
#' @param distancekm logical value indicating whether to calculate the distances in km, default value is \code{FALSE}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{h_values} - Initial values tested for the bandwidth.
#' \item \code{iterations} - All bandwidth values tested and respective cv/aic results for each Golden Section Search executed.
#' \item \code{gss_results} - Optimum bandwidth found and respective cv/aic.
#' \item \code{min_bandwidth} - Optimum bandwidth.
#' }
#'
#' @examples
#' ## Data
#'
#'
#' data(southkorea_covid19)
#'
#'
#' ## GSS algorithm
#'
#' gss <- Golden(data = southkorea_covid19[1:122, ],
#' formula = n_covid1~diff_sd,
#' xvarinf = NULL, weight = NULL, lat = "y", long = "x",
#' offset = NULL, model = "poisson", method = "fixed_g",
#' bandwidth = "cv", globalmin = FALSE, distancekm = FALSE,
#' force=FALSE)
#'
#' ## Bandwidth
#' gss$min_bandwidth
#'
#' ## Iterations
#' gss$iterations
#'
#' @importFrom sp spDistsN1
#'
#' @importFrom stats model.extract model.matrix dist
#'
#'
#' @export

Golden <- function(data, formula, xvarinf=NULL, weight=NULL,
                   lat, long, globalmin=TRUE,
                   method, model="zinb", bandwidth="cv", offset=NULL,
                   force=FALSE, maxg=100, distancekm=FALSE){ #flag -> nulls
  output <- list()
  E <- 10
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  XVAR <- attr(mt, "term.labels")
  y <- model.extract(mf, "response")
  N <- length(y)
  x <- model.matrix(mt, mf)
  if (is.null(xvarinf)){
    G <- matrix(1, N, 1)
    lambdag <- matrix(0, ncol(G), 1)
  }
  else{
    G <- as.matrix(data[, xvarinf])
    G <- cbind(rep(1, N), G)
  }
  wt <- rep(1, N)
  if (!is.null(weight)){
    wt <- unlist(data[, weight])
  }
  Offset <- rep(0, N)
  if (!is.null(offset)){
    Offset <- unlist(data[, offset])
  }
  nvar <- ncol(x)
  yhat <- rep(0, N)
  yhat2 <- rep(0, N)
  pihat <- rep(0, N)
  alphai <- rep(0, N)
  S <- rep(0, N)
  Si <- rep(0, N)
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  pos0 <- which(y==0)
  pos02 <- which(y==0)
  pos1 <- which(y>0)

  #### global estimates ####
  uj <- (y+mean(y))/2
  nj <- log(uj)
  parg <- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  cont3 <- 0
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (model == "zip" | model == "poisson"){
      parg <- 1/E^(-6)
      alphag <- 1/parg
    }
    if (model == "zinb" | model == "negbin"){
      if (cont>1){
        parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }
      while (abs(dpar)>0.0001 & cont1<200){
        if (parg<0){
          parg <- 0.00001
        }
        parg <- ifelse(parg<E^-10, E^-10, parg)
        gf <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hess <- ifelse(hess==0, E^-23, hess)
        par0 <- parg
        parg <- par0-as.vector(solve(hess, tol=E^-60))*gf
        if (parg>E^5){
          dpar <- 0.0001
          cont3 <- cont3+1
          if (cont3==1){
            parg <- 2
          }
          else if (cont3==2) {
            parg <- E^5
          }
          else if (cont3==3){
            parg <- 0.0001
          }
        }
        else{
          dpar <- parg-par0
          cont1 <- cont1+1
        }
        if (parg>E^6){
          parg <- E^6
          dpar <- 0
        }
      }
      alphag <- 1/parg
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<100){
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))<E^-60){
        bg <- rep(0,ncol(x))
      }
      else{
        bg <- solve(t(x)%*%(Ai*x), tol=E^-60)%*%t(x)%*%(Ai*zj)
      }
      nj <- as.vector(x%*%bg+Offset)
      nj <- ifelse(nj>700,700,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      uj <- ifelse(uj>100000,100000,uj)
      tt <- y/uj
      tt <- ifelse(tt==0,E^-10,tt)
      devg <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
      if (cont2>100){
        ddev <- 0.0000001
      }
      else{
        ddev <- devg-olddev
      }
      cont2 <- cont2+1
    }
    cont <-cont+1
    ddpar <- parg-parold
  }
  if (!is.null(xvarinf)){
    lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
    if (lambda0<=0){
      lambdag <- rep(0, ncol(G))
    }
    else{
      lambda0 <- log(lambda0/(1-lambda0))
      lambdag <- c(lambda0, rep(0, ncol(G)-1)) #flag
    }
    pargg <- parg
    ujg <- uj
    if (length(pos0)==0 | !any(lambdag)){
      if (length(pos0)==0){
        pos0 <- pos1
        if (model=="zinb" | model=="zip"){
          model <- "negbin"
        }
      }
      if (!force){
        model <- "negbin"
      }
    }
  }
  njl <- G%*%lambdag
  if (model!="zip" & model!="zinb"){
    zkg <- 0
  }
  else{
    zkg <- 1/(1+exp(-G%*%lambdag)*(parg/(parg+uj))^parg)
    zkg <- ifelse(y>0,0,zkg)
  }
  dllike <- 1
  llikeg <- 0
  j <- 0
  contador <- 0
  while (abs(dllike)>0.00001 & j<600){
    contador <- contador + 1
    ddpar <- 1
    cont <- 1
    contador2 <- 0
    while (abs(ddpar)>0.000001 & cont<100){
      contador2 <- contador2+1
      dpar <- 1
      parold <- parg
      aux1 <- 1
      aux2 <- 1
      aux3 <- 1
      cont3 <- 0
      int <- 1
      if (model=="zip" | model=="poisson"){
        alphag <- E^-6
        parg <- 1/alphag
      }
      else{
        if (j>0){
          parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
        }
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            parg <- 0.00001
          }
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hess <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <- as.vector(par0-solve(hess, tol=E^-60)%*%gf)
          if (aux2>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg <- 2
            }
            else if (cont3==2){
              parg <- E^5
            }
            else if (cont3==3){
              parg <- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          if (parg>E^6){
            parg <- E^6
            dpar <- 0
          }
          aux2 <- aux2+1
        }
        alphag <- 1/parg
      }
      devg <- 0
      ddev <- 1
      nj <- x%*%bg+Offset
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- as.vector((1-zkg)*((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2)))))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        uj <- ifelse(uj<E^-150, E^-150, uj)
        zj <- (nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2))))*(1+alphag*uj)))-Offset
        if (det(t(x)%*%(Ai*x))<E^-60){
          bg <- rep(0, nvar)
        }
        else{
          bg <- solve(t(x)%*%(Ai*x), tol=E^-60)%*%t(x)%*%(Ai*zj)
        }
        nj <- x%*%bg+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        olddev <- devg
        gamma1 <- (uj/(uj+parg))^y*(parg/(uj+parg))^parg #(gamma(par+y)/(gamma(y+1)#gamma(par)))#
        gamma1 <- ifelse(gamma1<=0, E^-10, gamma1)
        devg <- sum((1-zkg)*(log(gamma1)))
        ddev <- devg-olddev
        aux1 <- aux1+1
      }
      ddpar <- parg-parold
      cont <- cont+1
    }
    if (model == "zip" |model == 'zinb'){
      devg <- 0
      ddev <- 1
      njl <- G%*%lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < (-maxg),-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      contador3 <- 0
      while (abs(ddev)>0.000001 & aux3<100){
        contador3 <- contador3 + 1
        Ai <- as.vector(pig*(1-pig))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if (det(t(G*Ai)%*%G)<E^-60){
          lambdag <- matrix(0, ncol(G), 1)
        }
        else{
          lambdag <- solve(t(G*Ai)%*%G, tol=E^-60)%*%t(G*Ai)%*%zj
        }
        njl <- G%*%lambdag
        njl <- ifelse(njl > maxg, maxg, njl)
        njl <- ifelse(njl < (-maxg),-maxg, njl)
        pig <- exp(njl)/(1+exp(njl))
        olddev <- devg
        devg <- sum(zkg*njl-log(1+exp(njl)))
        ddev <- devg-olddev
        aux3 <- aux3+1
      }
    }
    zkg <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
    zkg <- ifelse(y>0, 0, zkg)
    if (model != 'zip' & model != 'zinb'){
      zkg <- 0
    }
    oldllike <- llikeg
    llikeg <- sum(zkg*(njl)-log(1+exp(njl))+(1-zkg)*(log(gamma1)))
    dllike <- llikeg-oldllike
    j <- j+1
  }
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  COORD <- matrix(c(long, lat), ncol=2)
  sequ <- 1:N
  cv <- function(h){
    for (i in 1:N){
      seqi <- rep(i, N)
      dx <- sp::spDistsN1(COORD,COORD[i,])
      distan <- cbind(seqi, sequ, dx)
      if (distancekm){
        distan[,3] <- distan[,3]*111
      }
      u <- nrow(distan)
      w <- rep(0, u)
      for (jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        if (method=="fixed_bsq"){
          w[jj] <- (1-(distan[jj,3]/h)^2)^2
        }
        if (bandwidth=="cv"){
          w[i] <- 0
        }
      }
      if (method=="fixed_bsq"){
        position <- which(distan[,3]<=h)
        w[position] <- 0
      }
      if (method=="adaptive_bsq"){
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, N, 2)
        hn <- distan[h,3]
        for (jj in 1:N){
          if (distan[jj,4]<=h){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
        }
        if (bandwidth=="cv"){
          w[which(w[,2]==i)] <- 0
        }
        w <- w[order(w[, 2]), ]
        w <- w[,1]
      }
      b <- bg
      nj <- x%*%b+Offset
      uj <- exp(nj)
      par <- parg
      lambda <- lambdag
      njl <- G%*%lambda
      njl <- ifelse(njl>maxg, maxg, njl)
      njl <- ifelse(njl<(-maxg), -maxg, njl)
      if (model!="zip" & model!="zinb"){
        zk <- 0
      }
      else{
        lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
        if (lambda0>0){
          lambda0 <- log(lambda0/(1-lambda0))
          lambda <- c(lambda0, rep(0, ncol(G)-1)) #flag
          njl <- G%*%lambda
        }
        zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
        zk <- ifelse(y>0, 0, zk)
      }
      dllike <- 1
      llike <- 0
      j <- 1
      contador4 <- 0
      while (abs(dllike)>0.00001 & j<=600){
        contador4 <- contador4+1
        ddpar <- 1
        #while (abs(ddpar)>0.000001)
        dpar <- 1
        parold <- par
        aux1 <- 1
        aux2 <- 1
        int <- 1
        if (model=="zip" | model=="poisson"){
          alpha <- E^-6
          par <- 1/alpha
        }
        else{
          if (par<=E^-5){
            if (i>1){
              par <- 1/alphai[i-1,2]
            }
          }
          if (par>=E^6){
            par <- E^6
            dpar <- 0
            alpha <- 1/par
            b <- bg
            uj <- exp(x%*%b+Offset)
            lambda <- lambdag
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
            zk <- ifelse(y>0, 0, zk)
            if (any(lambda)==0){
              zk <- 0
            }
          }
          while (abs(dpar)>0.000001 & aux2<200){
            par <- ifelse(par<E^-10, E^-10, par)
            gf <- sum(w*wt*(1-zk)*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
            hess <- sum(w*wt*(1-zk)*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
            par0 <- par
            par <- as.vector(par0-solve(hess, tol=E^-60)%*%gf)
            dpar <- par-par0
            if (par>=E^6){
              par <- E^6
              dpar <- 0
              alpha <- 1/par
              b <- bg
              uj <- exp(x%*%b+Offset)
              lambda <- lambdag
              njl <- G%*%lambda
              njl <- ifelse(njl>maxg, maxg, njl)
              njl <- ifelse(njl<(-maxg),-maxg,njl)
              zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
              zk <- ifelse(y>0,0,zk)
              if (any(lambda)==0){
                zk <- 0
              }
            }
            aux2 <- aux2+1
          }
          if (par<=E^-5){
            par <- E^6
            b <- bg
            uj <- exp(x%*%b+Offset)
            lambda <- lambdag
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
            zk <- ifelse(y>0, 0, zk)
            if (any(lambda)==0){
              zk <- 0
            }
          }
          alpha <- 1/par
        }
        dev <- 0
        ddev <- 1
        nj <- x%*%b+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        contador5 <- 0
        while (abs(ddev)>0.000001 & aux1<100){
          contador5 <- contador5+1
          uj <- ifelse(uj>E^100,E^100,uj)
          Ai <- as.vector((1-zk)*((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2)))))
          Ai <- ifelse(Ai<=0,E^-5,Ai)
          uj <- ifelse(uj<E^-150,E^-150,uj)
          denz <- (((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))*(1+alpha*uj))
          denz <- ifelse(denz==0,E^-5,denz)
          zj <- (nj+(y-uj)/denz)-Offset
          if (det(t(x)%*%(w*Ai*x*wt))<E^-60){
            b <- rep(0, nvar)
          }
          else{
            b <- solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x)%*%(w*Ai*wt*zj)
          }
          nj <- x%*%b+Offset
          nj <- ifelse(nj>700, 700, nj)
          nj <- ifelse(nj<(-700), -700, nj)
          uj <- exp(nj)
          olddev <- dev
          uj <- ifelse(uj>E^10, E^10, uj)
          uj <- ifelse(uj==0, E^-10, uj)
          temp <- (uj/(uj+par))
          temp <- ifelse(temp<(E^-307), 0, temp)
          if (par==E^6){
            gamma1 <- temp^y*exp(-uj)
            gamma1 <- ifelse(temp==0 & y==0, NA, gamma1)
          }
          else{
            gamma1 <- temp^y*(par/(uj+par))^par
            gamma1 <- ifelse(temp==0 & y==0, NA, gamma1)
          }
          gamma1 <- ifelse(gamma1<=0 | is.na(gamma1), E^-10, gamma1)
          dev <- sum((1-zk)*(log(gamma1)))
          ddev <- dev-olddev
          aux1 <- aux1+1
        }
        ddpar <- par-parold
        if (model=="zip" | model=="zinb"){
          if (j==1){
            alphatemp <- alpha
            lambdatemp <- lambda[1]
          }
          else{
            alphatemp <- c(alphatemp, alpha)
            lambdatemp <- c(lambdatemp, lambda[1])
          }
          alphatemp <- round(alphatemp, 7)
          ambdatemp <- round(lambdatemp, 7)
          if (model=="zinb"){
            condition <- (j>300 & length(alphatemp)>length(unique(alphatemp)) & length(lambdatemp)>length(unique(lambdatemp)))
          }
          else if (model=="zip"){
            condition <- (j>300 & length(lambdatemp)>length(unique(lambdatemp)))
          }
          if (condition){
            lambda <- rep(0, ncol(G))
            njl <- G%*%lambda
            zk <- rep(0, N)
          }
          else{
            aux3 <- 1
            dev <- 0
            ddev <- 1
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            pi <- exp(njl)/(1+exp(njl))
            contador6 <- 0
            while (abs(ddev)>0.000001 & aux3<100){
              contador6 <- contador6+1
              Aii <- as.vector(pi*(1-pi))
              Aii <- ifelse(Aii<=0, E^-5, Aii)
              zj <- njl+(zk-pi)/Aii
              if (det(t(G*Aii*w*wt)%*%G)<E^-60){
                lambda <- matrix(0, ncol(G), 1)
              }
              else{
                lambda <- solve(t(G*Aii*w*wt)%*%G, tol=E^-60)%*%t(G*Aii*w*wt)%*%zj
              }
              njl <- G%*%lambda
              njl <- ifelse(njl>maxg, maxg, njl)
              njl <- ifelse(njl<(-maxg), -maxg, njl)
              pi <- exp(njl)/(1+exp(njl))
              olddev <- dev
              dev <- sum(zk*njl-log(1+exp(njl)))
              ddev <- dev-olddev
              aux3 <- aux3+1
            }
          }
        }
        njl <- G%*%lambda
        njl <- ifelse(njl>maxg, maxg, njl)
        njl <- ifelse(njl<(-maxg), -maxg, njl)
        zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
        zk <- ifelse(y>0, 0, zk)
        if (any(lambda)==0){
          zk <- rep(0, N)
        }
        if (model!="zip" & model!="zinb"){
          zk <- 0
        }
        oldllike <- llike
        llike <- sum(zk*(njl)-log(1+exp(njl))+(1-zk)*(log(gamma1)))
        dllike <- llike-oldllike
        j <- j+1
      }
      yhat[i] <- uj[i]
      pihat[i] <- njl[i]
      alphai_ <- get("alphai")
      alphai_[i] <- alpha
      assign("alphai", alphai_, envir=parent.frame())
      alphai[i] <-  alpha
      if (det(t(x)%*%(w*Ai*x*wt))<E^-60){
        S[i] <- 0
      }
      else{
        S[i] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x*w*Ai*wt))[i]
      }
      if (model=="zip" | model=="zinb"){
        yhat[i] <- (uj*(1-exp(njl)/(1+exp(njl))))[i]
        yhat2[i] <- uj[i]
        if (det(t(G)%*%(w*Aii*G*wt))<E^-60){
          Si[i] <- 0
        }
        else{
          Si[i] <- (G[i,]%*%solve(t(G)%*%(w*Aii*G*wt), tol=E^-60)%*%t(G*w*Aii*wt))[i]
          if (any(lambda)==0){
            Si[i] <- 0
          }
        }
      }
    }
    CV <- t((y-yhat)*wt)%*%(y-yhat)
    par_ <- 1/alphai
    if (model=="zinb" | model == "zip"){
      npar <- sum(S)+sum(Si)
      if (bandwidth=="aic"){
        if (any(lambda)==0){
          ll <- sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+(par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0]))+
            sum(-log(0+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                  y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1])))

          llnull1 <- sum(-log(1+zk[pos0])+log(zk[pos0]+(par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0]))+
            sum(-log(1+zk[pos1])+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                  y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1])))

          llnull2 <- sum(-log(1+0)+log(0+(par_/(par_+mean(y)))^par_))+
            sum(-log(1+0)+lgamma(par_+y)-lgamma(y+1)-lgamma(par_)+y*log(mean(y)/(par_+mean(y)))+par_*log(par_/(par_+mean(y))))
        }
        else{
          ll <- sum(-log(1+exp(pihat[pos0]))+log(exp(pihat[pos0])+(par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0]))+
            sum(-log(1+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                  y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1])))

          llnull1 <- sum(-log(1+zk[pos0])+log(zk[pos0]+(par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0]))+
            sum(-log(1+zk[pos1])+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                  y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1])))
        }
        dev <- 2*(llnull1-ll)
        npar <- sum(S)+sum(Si)
        AIC <- 2*npar-2*ll
        AICc <- AIC+2*(npar*(npar+1)/(N-npar-1))
        if (model=="zinb"){
          AIC <- 2*(npar+npar/(ncol(x)+ncol(G)))-2*ll
          AICc <- AIC+2*((npar+npar/(ncol(x)+ncol(G)))*((npar+npar/(ncol(x)+ncol(G)))+1)/(N-(npar+npar/(ncol(x)+ncol(G)))-1))
        }
      }
    }
    else if (model=="poisson" | model=="negbin"){
      npar <- sum(S)
      if (bandwidth=="aic"){
        if (length(pos02)==0){
          pos0 <- pos1
          pos0x <- 1
          pos0xl <- 1
        }
        else {
          pos0x <- (par_[pos0]/(par_[pos0]+yhat[pos0]))^par_[pos0]
          pos0xl <- (par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0]
          pos0x <- ifelse(pos0x==0, E^-10, pos0x)
        }
        ll <- sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+pos0x))+
          sum(-log(0+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                y[pos1]*log(yhat[pos1]/(par_[pos1]+yhat[pos1]))+
                par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat[pos1])))

        llnull1 <- sum(-log(1+zk)+log(zk+pos0xl))+ sum(-log(1+zk)+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                                                         y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1])))
        dev <- 2*(llnull1-ll)
        npar <- sum(S)
        AIC <- 2*npar-2*ll
        AICc <- AIC+2*(npar*(npar+1)/(N-npar-1))
        if (model == "negbin"){
          AIC <- 2*(npar+npar/ncol(x))-2*ll
          AICc <-AIC+2*(npar+npar/ncol(x))*(npar+npar/ncol(x)+1)/(N-(npar+npar/ncol(x))-1)
        }
      }
    }
    if (bandwidth == "aic"){
      CV <- AICc
    }
    res <- cbind(CV, npar)
    return (res)
  }

  #### defining golden section search parameters ####
  if (method=="fixed_g" | method=="fixed_bsq"){
    ax <- 0
    bx <- as.integer(max(dist(COORD))+1)
    if (distancekm){
      bx <- bx*111
    }
  }
  if (method=="adaptive_bsq"){
    ax <- 5
    bx <- N
  }
  r <- 0.61803399
  tol <- 0.1
  if (!globalmin){
    lower <- ax
    upper <- bx
    xmin <- matrix(0, 1, 2)
  }
  else{
    lower <- cbind(ax, (1-r)*bx, r*bx)
    upper <- cbind((1-r)*bx, r*bx, bx)
    xmin <- matrix(0, 3, 2)
  }
  for (GMY in 1:3){
    ax1 <- lower[GMY]
    bx1 <- upper[GMY]
    h0 <- ax1
    h3 <- bx1
    h1 <- bx1-r*(bx1-ax1)
    h2 <- ax1+r*(bx1-ax1)
    if (GMY==1){
      h_values <- data.frame('h0'=h0, 'h1'=h1, 'h2'=h2, 'h3'=h3)
    }
    else{
      h_values <- rbind(h_values, c(h0, h1, h2, h3))
    }
    ################################
    res1 <- cv(h1)
    CV1 <- res1[1]
    res2 <- cv(h2)
    CV2 <- res2[1]
    if (GMY==1){
      band <- data.frame('GSS_count'=GMY, 'h1'=h1, 'cv1'=CV1, 'h2'=h2, 'cv2'=CV2)
    }
    else{
      band <- rbind(band, c(GMY, h1, CV1, h2, CV2))
    }
    int <- 1
    while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200){
      if (CV2<CV1){
        h0 <- h1
        h1 <- h3-r*(h3-h0)
        h2 <- h0+r*(h3-h0)
        CV1 <- CV2
        res2 <- cv(h2)
        CV2 <- res2[1]
      }
      else{
        h3 <- h2
        h1 <- h3-r*(h3-h0)
        h2 <- h0+r*(h3-h0)
        CV2 <- CV1
        res1 <- cv(h1)
        CV1 <- res1[1]
      }
      band <- rbind(band, c(GMY, h1, CV1, h2, CV2))
      int <- int+1
    }
    if (CV1<CV2){
      golden <- CV1
      xmin[GMY,1] <- golden
      xmin[GMY,2] <- h1
      npar <- res1[2]
      if (method=="adaptive_bsq"){
        xmin[GMY,2] <- floor(h1)
      }
    }
    else{
      golden <- CV2
      xmin[GMY,1] <- golden
      xmin[GMY,2] <- h2
      npar <- res2[2]
      if (method=="adaptive_bsq"){
        xmin[GMY,2] <- floor(h2)
      }
    }
    if (!globalmin){
      break
    }
  }
  min_bandwidth <- as.data.frame(xmin)
  names(min_bandwidth) <- c(bandwidth, 'bandwidth')
  output <- append(output, list(h_values))
  names(output)[length(output)] <- "h_values"
  output <- append(output, list(band))
  names(output)[length(output)] <- "iterations"
  output <- append(output, list(min_bandwidth))
  names(output)[length(output)] <- "gss_results"
  if (globalmin){
    message('Global Minimum (Da Silva and Mendes, 2018)')
  }
  hh <- min_bandwidth[which(unlist(min_bandwidth[, bandwidth])==min(unlist(min_bandwidth[, bandwidth]))), 'bandwidth']
  output <- append(output, list(hh))
  names(output)[length(output)] <- "min_bandwidth"
  message('Bandwidth: ', hh)
  invisible(output)
}
