#' @title Geographically Weighted Zero Inflated Negative Binomial Regression
#'
#' @description Fits a geographically weighted regression model using zero inflated probability distributions. Has the zero inflated negative binomial distribution (zinb) as default, but also accepts the zero inflated Poisson (zip), negative binomial (negbin) and Poisson distributions. Can also fit the global versions of each regression model.
#'
#' @param data name of the dataset.
#' @param formula regression model formula as in \code{lm}.
#' @param xvarinf name of the covariates for the zero inflated part of the model, default value is \code{NULL}.
#' @param weight name of the variable containing the sample weights, default value is \code{NULL}.
#' @param lat name of the variable containing the latitudes in the dataset.
#' @param long name of the variable containing the longitudes in the dataset.
#' @param grid name of the dataset containing the coordinates for the model locations, default value is \code{NULL}.
#' @param method indicates the method to be used for the bandwidth calculation (\code{adaptive_bsq} or \code{fixed_g}).
#' @param model indicates the model to be used for the regression (\code{zinb}, \code{zip}, \code{negbin}, \code{poisson}), default value is\code{"zinb"}.
#' @param offset name of the variable containing the offset values, if null then is set to a vector of zeros, default value is \code{NULL}.
#' @param distancekm logical value indicating whether to calculate the distances in km, default value is \code{FALSE}.
#' @param force logical value indicating whether to force the indicated model even if it is not the best fit for the data, default value is \code{FALSE}.
#' @param int_inf logical value indicating whether to include an intercept in the zero inflated part of the model, default value is \code{TRUE}.
#' @param maxg integer indicating the maximum number of iterations for the zero inflated part of the model, default value is \code{100}.
#' @param h integer indicating the bandwidth value (obtained from \code{golden()}), default value is \code{NULL}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{bandwidth} - Bandwidth value.
#' \item \code{measures} - Goodness of fit statistics and other measures.
#' \item \code{qntls_gwr_param_estimates} - Quantiles of GWR parameter estimates.
#' \item \code{descript_stats_gwr_param_estimates} - Descriptive statistics of GWR parameter estimates.
#' \item \code{t_test_gwr_param_estimates} - Results for the parameters significance t tests.
#' \item \code{qntls_gwr_se} - Quantiles of GWR standard errors.
#' \item \code{descript_stats_gwr_se} - Descriptive statistics of GWR standard errors.
#' \item \code{qntls_gwr_zero_infl_param_estimates} - Quantiles of GWR zero inflated parameter estimates.
#' \item \code{descript_stats_gwr_zero_infl_param_estimates} - Descriptive statistics of GWR zero inflated parameter estimates.
#' \item \code{t_test_gwr_zero_infl_param_estimates} - Results for the zero inflated parameters significance t tests.
#' \item \code{qntls_gwr_zero_infl_se} - Quantiles of GWR zero inflated standard errors.
#' \item \code{descript_stats_gwr_zero_infl_se} - Descriptive statistics of GWR zero inflated standard errors.
#' \item \code{non_stationary_test} - Results for the Non-Stationary Test for GWR parameter estimates.
#' \item \code{non_stationary_test_zero_infl} - Results for the Non-Stationary Test for GWR zero inflated parameter estimates.
#' \item \code{global_param_estimates} - Parameter estimates for the global model.
#' \item \code{analysis_max_like_zero_infl_param_estimated} - Analysis of Maximum Likelihood Zero Inflation Parameter Estimates.
#' \item \code{analysis_max_like_gof_measures} - Goodness of fit measures for the Analysis of Maximum Likelihood Zero Inflation Parameter Estimates.
#' \item \code{variance_covariance_matrix} - Variance-covariance matrix.
#' \item \code{residuals} - Model residuals.
#' \item \code{param_estimates_grid} - GWR parameter estimates using grid dataset.
#' \item \code{alpha_estimates} - Estimates for the alpha parameter (for zinb and negbin).
#' \item \code{gwr_param_estimates} - GWR parameter estimates.
#' }
#'
#' @examples
#' ## Data
#'
#'
#' data(southkorea_covid19)
#'
#'
#' ## Model
#'
#' mod <- gwzinbr(data = southkorea_covid19,
#' formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
#' diff_sd+Crowding+Migration+Health_behavior,
#' lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
#' model = "negbin", distancekm = TRUE, h=230, force=TRUE)
#'
#' ## Bandwidth
#' mod$bandwidth
#'
#' ## Goodness of fit measures
#' mod$measures
#'
#' @importFrom sp spDistsN1
#'
#' @importFrom stats model.extract model.matrix pnorm pt qt quantile pf
#'
#'
#' @export

gwzinbr <- function(data, formula, xvarinf=NULL, weight=NULL,
                    lat, long, grid=NULL, method, model = "zinb",
                    offset=NULL, distancekm=FALSE, force=FALSE, int_inf=TRUE,
                    maxg=100, h=NULL){
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
    if (int_inf){
      G <- cbind(rep(1, N), G)
    }
  }
  yhat <- rep(0, N)
  yhat2 <- rep(0, N)
  pihat <- rep(0, N)
  nvar <- ncol(x)
  wt <- rep(1, N)
  if (!is.null(weight)){
    wt <- unlist(data[, weight])
  }
  Offset <- rep(0, N)
  if (!is.null(offset)){
    Offset <- unlist(data[, offset])
  }
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  Iy2 <- Iy
  pos0 <-  which(y==0)
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
        hessg <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hessg <- ifelse(hessg==0, E^-23, hessg)
        par0 <- parg
        parg <- par0-as.vector(solve(hessg))*gf
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
    while (abs(ddev)>0.000001 & cont2<200){
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0) {
        bg <- rep(0,ncol(x))
      }
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
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
    cont <- cont+1
    ddpar <- parg-parold
  }
  if (!is.null(xvarinf)){
    lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
    if (lambda0 <= 0) {
      lambdag <- rep(0, ncol(G))
      message("NOTE: Expected number of zeros (", round(sum((parg/(uj+parg))^parg), 2),
              ") >= number of zeros (", length(pos0), "). No Need for Zero Model.")
    }
    else{
      message("NOTE: Expected number of zeros (", round(sum((parg/(uj+parg))^parg), 2),
              ") < number of zeros (", length(pos0), "). Zero Model Used.")
      lambda0 <- log(lambda0/(1-lambda0))
      lambdag <- rbind(lambda0, rep(0, ncol(G)-1))
    }
    pargg <- parg
    ujg <- uj
    if (length(pos0)==0 | any(lambdag)==0){
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
  while (abs(dllike)>0.00001 & j<600){
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
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
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            parg <- 0.00001
          }
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hessg <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hessg <- ifelse(hessg==0, E^-23, hessg)
          par0 <- parg
          parg <- as.vector(par0-solve(hessg)%*%gf)
          if ( aux2 > 50 |parg > E^5) {
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1) {
              parg <- 2
            }
            else if (cont3==2) {
              parg <- E^5
            }
            else if (cont3==3) {
              parg <- 0.0001
            }
          }
          else {
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
      nj <- ifelse(nj>700, 700, nj)
      nj <- ifelse(nj<(-700), -700, nj)
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- as.vector((1-zkg)*((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2)))))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        uj <- ifelse(uj<E^-150, E^-150, uj)
        zj <- (nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2))))*(1+alphag*uj)))-Offset
        if (det(t(x)%*%(Ai*x))==0){
          bg <- rep(0, nvar)
        }
        else{
          bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
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
    if (model == "zip" | model == "zinb"){
      devg <- 0
      ddev <- 1
      njl <- G%*%lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < (-maxg),-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      while (abs(ddev)>0.000001 & aux3<100){
        Ai <- as.vector(pig*(1-pig))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if(det(t(G)%*%(Ai*G))==0){
          lambdag <- matrix(0, ncol(G), 1)
        }
        else {
          lambdag <- solve(t(G)%*%(Ai*G))%*%t(G)%*%(Ai*zj)
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
  g1x <- parg/(parg+uj)
  g2x <- uj/(parg+uj)
  hgx <- exp(njl)+g1x^parg
  daa <- zkg*((g1x^parg*(log(g1x)+g2x))^2*(1-1/hgx)/hgx+g1x^parg*(g2x^2/parg)/hgx)+(1-zkg)*(trigamma(parg+y)-trigamma(parg)-2/(uj+parg)+1/parg+(y+parg)/(uj+parg)^2)
  dab <- zkg*(g1x^(2*parg+1)*uj*(log(g1x)+g2x)/hgx^2-g1x^parg*(-g2x^2+parg*g2x*(log(g1x)+g2x))/hgx)+(1-zkg)*(g2x*(y-uj)/(uj+parg))
  dal <- -zkg*(exp(njl)*g1x^parg*(log(g1x)+g2x)/hgx^2)
  daa <- daa*parg^4
  dab <- dab*parg^2
  if (any(lambdag)==0){
    Iy <- matrix(0, length(y), 1)
  }
  dll <- as.vector(Iy*(exp(njl)*g1x^parg/hgx^2)-exp(njl)/(1+exp(njl))^2)
  dbb <- as.vector(Iy*(-(parg*g1x^parg*g2x/hgx)^2+parg^2*g1x^parg*g2x^2*(1-1/uj)/hgx)-(1-Iy)*(parg*g2x*(1+(y-uj)/(parg+uj))))
  dlb <- as.vector(Iy*(parg*exp(njl)*g1x^parg*g2x/hgx^2))
  I1 <- matrix(1, length(y), 1)
  II <- rbind(cbind(-(t(I1 * daa)) %*% I1, -(t(I1 * dab)) %*% x, -(t(I1 * dal)) %*% G),
              cbind(-(t(x) %*% (dab * I1)), -t(x*dbb) %*% x, -t(x * dlb) %*% G),
              cbind(-(t(G) %*% (dal * I1)), -t(G) %*% (x * dlb), -t(G * dll) %*% G))
  if (all(lambdag)>0 & alphag==E^-6){
    II <- II[2:nrow(II), 2:nrow(II)]
  }
  else if(any(lambdag)==0 & alphag>E^-6){
    II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
  }
  else if(any(lambdag)== 0 & alphag==E^-6){
    II <- II[2:(ncol(x)+1), 2:(ncol(x)+1)]
  }
  varabetalambdag <- diag(solve(II))
  stdabetalambdag <- sqrt(abs(varabetalambdag))
  varcovg <- solve(II)
  ##################
  output <- append(output, list(h))
  names(output)[length(output)] <- "bandwidth"
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  COORD <- matrix(c(lat, long), ncol=2)
  if (is.null(grid)){
    POINTS <- matrix(c(lat, long), ncol=2)
  }
  else{
    long2 <- unlist(grid[, long])
    lat2 <- unlist(grid[, lat])
    POINTS <- matrix(c(lat2, long2), ncol=2)
  }
  mm <- nrow(COORD)
  bi <- matrix(0, ncol(x)*mm, 4)
  li <- matrix(0, ncol(G)*mm, 4)
  alphai <- matrix(0, mm, 3)
  BB <- matrix(0, ncol(x)*N, N)
  BBl <- matrix(0, ncol(G)*N, N)
  sumwi <- matrix(0, mm, 1)
  varbi <- matrix(0, ncol(x)*mm, 1)
  varli <- matrix(0, ncol(G)* mm, 1)
  S <- matrix(0, mm, 1)
  Si <- matrix(0, mm, 1)
  S2 <- matrix(0, mm, 1)
  biT <- matrix(0, mm, ncol(x)+1)
  ym <- y-mean(y)

  #### calculating distance ####
  sequ <- 1:N
  for (i in 1:mm){
    seqi <- rep(i, N)
    dx <- sp::spDistsN1(COORD, COORD[i,])
    distan <- cbind(seqi, sequ, dx)
    if (distancekm){
      distan[,3] <- distan[,3]*111
    }
    u <- nrow(distan)
    w <- rep(0, u)
    for(jj in 1:u){
      if(method=="fixed_g"){
        w[jj] <- exp(-0.5*(distan[jj, 3]/h)^2)
      }
      else if(method=="fixed_bsq"){
        w[jj] <- (1 -(distan[jj, 3]/h)^2)^2
      }
    }
    if(method=="adaptive_bsq"){
      distan <- distan[order(distan[,3]),]
      distan <- cbind(distan, 1:nrow(distan))
      w <- matrix(0,N,2)
      hn <- distan[h,3]
      for (jj in 1:N){
        if(distan[jj, 4] <= h){
          w[jj, 1] <- (1 -(distan[jj, 3]/hn)^2)^2
        }
        else{
          w[jj, 1] <- 0
        }
        w[jj, 2] <- distan[jj, 2]
      }
      w <- w[order(w[, 2]), 1]
    }
    #### model selection ####
    Iy <- Iy2
    b <- bg
    b2 <- b
    nj <- x%*%b+Offset
    uj <- exp(nj)
    par <- parg
    par2 <- par
    lambda <- lambdag
    njl <- ifelse(njl > maxg, maxg, njl)
    njl <- ifelse(njl < (-maxg), -maxg, njl)
    if(model != "zip" & model != "zinb" ){
      zk <- 0
    }
    else{
      lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
      if (lambda0 > 0){
        lambda0 <- log(lambda0/(1-lambda0))
        lambda <- matrix(c(lambda0, rep(0, ncol(G)-1)), ncol=1)
        njl <- G%*%lambda
      }
      zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
      zk <- ifelse(y>0, 0, zk)
    }
    dllike <- 1
    llike <- 0
    j <- 1
    while (abs(dllike) > 0.00001 & j <= 600){
      ddpar <- 1
      dpar <- 1
      parold <- par
      aux1 <- 1
      aux2 <- 1
      int <- 1
      if (model == "zip" || model == "poisson") {
        alpha <- E^-6
        par <- 1/alpha
      }
      if(model == "zinb" || model == "negbin"){
        if(par >= E^6){
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
          zk <- ifelse(y>0, 0, zk)
          if (any(lambda) == 0){
            zk <- 0
          }
        }
        while (abs(dpar)>0.000001 & aux2<200){
          par <- ifelse(par < E^-10, E^-10, par)
          gf <- sum(w*wt*(1-zk)*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj) - (par+y)/(par+uj)))
          hess <- sum(w*wt*(1-zk)*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
          hess <- ifelse(hess == 0, E^-23, hess)
          par0 <- par
          par <- as.vector(par0 - solve(hess)%*%gf)
          dpar <- par - par0
          if (par >= E^6){
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
            zk <- ifelse(y>0, 0, zk)
            if (any(lambda) == 0){
              zk <- 0
            }
          }
          aux2 <- aux2 + 1
        }
        if (par <= E^-5){
          par <- E^6
          b <- bg
          uj <- exp(x%*%b+Offset)
          lambda <- lambdag
          njl <- G %*% lambda
          njl <- ifelse(njl>maxg, maxg, njl)
          njl <- ifelse(njl<(-maxg),-maxg,njl)
          zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
          zk <- ifelse(y>0, 0, zk)
          if (any(lambda) == 0){
            zk <- 0
          }
        }
        alpha <- 1/par
      }
      dev <- 0
      ddev <- 1
      nj <- x%*%b+Offset
      nj <- ifelse(nj>700,700,nj)
      nj <- ifelse(nj<(-700), -700, nj)
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- as.vector((1-zk)*((uj/(1+alpha*uj) + (y-uj) * (alpha*uj/(1+2*alpha*uj+alpha^2*uj^2)))))
        Ai <- ifelse(Ai <=0, E^-5, Ai)
        uj <- ifelse(uj < E^-150, E^-150, uj)
        denz <- (((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))*(1+alpha*uj))
        denz <- ifelse(denz == 0, E^-5, denz)
        zj <- (nj+(y-uj)/denz)-Offset
        if(det(t(x) %*% (w*Ai*x*wt)) == 0){
          b <- matrix(0, nvar, 1)
        }
        else {
          b <- solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x)%*%(w*Ai*wt*zj)
        }
        nj <- x%*%b+Offset
        nj <- ifelse(nj>700,700,nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        olddev <- dev
        uj <- ifelse(uj > E^10, E^10, uj)
        uj <- ifelse(uj == 0, E^10, uj)
        if (par == E^6){
          gamma1 <- (uj/(uj+par))^y*exp(-uj)
        }
        else{
          gamma1 <- (uj/(uj+par))^y*(par/(uj+par))^par
        }
        gamma1 <- ifelse(gamma1 <=0, E^-10, gamma1)
        dev <- sum((1-zk)*log(gamma1))
        ddev <- dev - olddev
        aux1 <- aux1 + 1
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
        lambdatemp <- round(lambdatemp, 7)
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
          while (abs(ddev)>0.000001 & aux3<100){
            Aii <- as.vector(pi*(1-pi))
            Aii <- ifelse(Aii<=0, E^-5, Aii)
            zj <- njl+(zk-pi)/Aii
            if (det(t(G*Aii*w*wt)%*%G)==0){
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
    #### computing variance of beta and lambda ####
    if (det(t(x)%*%(w*Ai*x*wt))==0){
      C <- matrix(0, ncol(x),nrow(x))
    }
    else{
      C <- t(t(solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x))*(w*Ai*wt))
    }
    Ci <- matrix(0, ncol(G), 1)
    if (model == "zip" || model == "zinb"){
      if (det(t(G)%*%(w*Aii*G*wt))==0){
        Ci <- matrix(0, ncol(G), nrow(G))
      }
      else{
        Ci <- t(t(solve(t(G)%*%(w*Aii*G*wt), tol=E^-60)%*%t(G))*(w*Aii*wt))
      }
      if(any(lambda)==0){
        Ci <- matrix(0, ncol(G), N)
      }
    }
    g1x <- par/(par+uj)
    g2x <- uj/(par+uj)
    hgx <- exp(njl)+g1x^par
    hgx <- ifelse(hgx > E^10, E^10, hgx)
    daa <- as.vector(w*wt*(zk*((g1x^par*(log(g1x)+g2x))^2*(1-1/hgx)/hgx+g1x^par*(g2x^2/par)/hgx)+
                             (1-zk)*(trigamma(par+y)-trigamma(par)-2/(uj+par)+1/par+(y+par)/(uj+par)^2)))
    dab <- as.vector(w*wt*(zk*(g1x^(2*par + 1)*uj*(log(g1x) + g2x) / hgx^2 - g1x^par * (-g2x^2+par*g2x*(log(g1x) + g2x)) / hgx) +
                             (1-zk) * (g2x*(y-uj) / (uj+par))))
    dal <- as.vector(-w*wt*zk*(exp(njl)*g1x^par*(log(g1x)+g2x) / hgx^2))
    daa <- daa*par^4
    dab <- dab*par^2
    if (any(lambda)==0) {
      Iy <- matrix(0, length(y), 1)
    }
    exphx <- 1 + exp(njl)
    exphx <- ifelse(exphx>E^90, E^90, exphx)
    dll <- as.vector(w*wt*(Iy*(exp(njl)*g1x^par/hgx^2)-exp(njl)/(exphx)^2))
    dbb <- as.vector(sqrt(w)*wt*(Iy*(-(par*g1x^par*g2x/hgx)^2+par^2*g1x^par*g2x^2*(1 - 1/uj)/hgx) -
                                   (1 - Iy)*(par*g2x*(1 + (y-uj)/(par+uj)))))
    dlb <- as.vector(w*wt*Iy*(par*exp(njl)*g1x^par*g2x/hgx^2))
    dll <- ifelse(is.na(dll), E^100, dll)
    daa <- ifelse(is.na(daa), E^100, daa)
    dab <- ifelse(is.na(dab), E^100, dab)
    dal <- ifelse(is.na(dal), E^100, dal)
    dbb <- ifelse(is.na(dbb), E^100, dbb)
    dlb <- ifelse(is.na(dlb), E^100, dlb)
    I1 <- matrix(1, length(y), 1)
    if(any(b)==0 & any(lambda)==0){
      II <- matrix(0, ncol(x)+ncol(G)+1, ncol(x)+ncol(G)+1)
    }
    else if(any(lambda)==0){
      dbb <- as.vector(w*wt*(Iy*(-(par*g1x^par*g2x/hgx)^2 + par^2*g1x^par*g2x^2*(1 - 1/uj)/hgx)-(1 - Iy)*(par*g2x*(1 + (y-uj) / (par+uj)))))
      if(det(t(x*dbb*dbb/Ai) %*% x)==0){
        II <- rbind(
          cbind(-(t(I1*daa))%*%I1, -(t(I1*dab))%*%x, -(t(I1*dal))%*%G),
          cbind(-t(x)%*%(dab*I1), -(t(x*dbb)%*%x), -t(x*dlb)%*%G),
          cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -(t(G*dll)%*%G))
        )
      }
      else{
        II <- rbind(
          cbind(-t(I1*daa)%*%I1, -t(I1*dab)%*%x, -t(I1*dal)%*%G),
          cbind(-t(x)%*%(dab*I1), -(t(x*dbb)%*%x)%*%solve(t(x*dbb*dbb/Ai)%*%x)%*%t(x*dbb)%*%x, -(t(x*dlb)%*%G)),
          cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -t(G*dll)%*%G)
        )
      }
    }
    else{
      II <- rbind(
        cbind(-t(I1*daa)%*%I1, -t(I1*dab)%*%x, -t(I1*dal)%*%G),
        cbind(-t(x)%*%(dab*I1), -(t(x*dbb)%*%x), -t(x*dlb)%*%G),
        cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -t(G*dll)%*%G)
      )
    }
    if (all(lambda) > 0 & alpha == E^-6){
      II <- II[2:nrow(II), 2:nrow(II)]
    }
    else if (any(lambda)==0 & alpha > E^-6){
      II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
    }
    else if (any(lambda) == 0 & alpha == E^-6){
      II <- II[2:(ncol(x)+1), 2:(ncol(x)+1)]
    }
    if (det(II) == 0) {
      if (all(lambda) > 0 & alpha == E^-6) {
        II <- II[1:ncol(x), 1:ncol(x)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        }
        else {
          varabetalambda <- rbind(diag(solve(II, tol=E^-60)), matrix(0, ncol(G), 1))
        }
      }
      else if (any(lambda) == 0 & alpha == E^-6) {
        II <- II[1:ncol(x), 1:ncol(x)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        }
        else {
          varabetalambda <- rbind(diag(solve(II,tol=E^-60)), matrix(0, ncol(G), 1))
        }
      }
      else {
        II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        }
        else {
          varabetalambda <- rbind(diag(solve(II, tol=E^-60)), matrix(0, ncol(G), 1))
        }
      }
    }
    else{
      varabetalambda <- diag(solve(II, tol=E^-60))
    }
    if (all(lambda) > 0 & alpha > E^-6){
      varb <- varabetalambda[2:(ncol(x)+1)]
      varl <- varabetalambda[(ncol(x)+2):length(varabetalambda)]
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(abs(varabetalambda[1]))
    }
    else if (all(lambda) > 0 & alpha == E^-6) {
      varb <- varabetalambda[1:ncol(x)]
      varl <- varabetalambda[(ncol(x)+1):length(varabetalambda)]
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(1/abs(-(t(I1*daa))%*%I1))
    }
    else if (any(lambda) == 0 & alpha > E^-6) {
      varb <- varabetalambda[2:(ncol(x)+1)]
      varl <- matrix(0, ncol(G), 1)
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(abs(varabetalambda[1]))
    }
    else if (any(lambda) == 0 & alpha == E^-6) {
      varb <- varabetalambda[1:ncol(x)]
      varl <- matrix(0, ncol(G), 1)
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(1/abs(-(t(I1*daa))%*%I1))
    }
    #############
    m1 <- (i-1)*ncol(x)+1
    m2 <- m1+(ncol(x)-1)
    bi[m1:m2, 1] <- i
    bi[m1:m2, 2] <- b
    bi[m1:m2, 3] <- COORD[i, 1]
    bi[m1:m2, 4] <- COORD[i, 2]
    varbi[m1:m2, 1] <- varb

    if (model == "zip" || model == "zinb") {
      m1 <- (i-1)*ncol(G)+1
      m2 <- m1 + (ncol(G)-1)
      li[m1:m2, 1] <- i
      li[m1:m2, 2] <- lambda
      li[m1:m2, 3] <- COORD[i, 1]
      li[m1:m2, 4] <- COORD[i, 2]
      varli[m1:m2, 1] <- varl
    }
    if (is.null(grid)) {
      r <- x[i, ]%*%C
      S[i] <- r[i]
      S2[i] <- r%*%t(r)
      yhat[i] <- uj[i]
      pihat[i] <- njl[i]
    }
    if (model == "zip" | model == "zinb") {
      ri <- G[i, ]%*%Ci
      Si[i] <- ri[i]
      yhat2[i] <- uj[i]
      yhat[i] <- (uj*(1-exp(njl)/(1+exp(njl))))[i]
    }
    #### creating non-stationarity matrix ####
    if (method != "adaptive_bsq"){
      CCC <- cbind(x, w, wt)
      m1 <- (i-1)*ncol(x)+1
      m2 <- m1 + (ncol(x)-1)
      if(det(t(CCC[, 1:ncol(x)])%*%(CCC[, ncol(CCC)-1]*CCC[, 1:ncol(x)]*CCC[, ncol(CCC)])) == 0){
        BB[m1:m2, ] <- matrix(0, ncol(x), nrow(x))
      }
      else{
        BB[m1:m2, ] <- t(t(solve(t(CCC[, 1:ncol(x)])%*%(CCC[, ncol(CCC)-1]*CCC[, 1:ncol(x)]*CCC[, ncol(CCC)]))%*%t(CCC[, 1:ncol(x)]))*(CCC[, ncol(CCC)-1]*CCC[, ncol(CCC)]))
      }
      if (model == "zip" || model == "zinb"){
        CCCl <- cbind(G, w, wt)
        m1 <- (i-1)*ncol(G)+1
        m2 <- m1+(ncol(G)-1)

        if (det(t(CCCl[, 1:ncol(G)])%*%(CCCl[, ncol(CCCl)-1]*CCCl[, 1:ncol(G)]*CCCl[, ncol(CCCl)])) == 0) {
          BBl[m1:m2, ] <- matrix(0, ncol(G), nrow(G))
        }
        else{
          BBl[m1:m2, ] <- t(t(solve(t(CCCl[, 1:ncol(G)])%*%(CCCl[, ncol(CCCl)-1] * CCCl[, 1:ncol(G)]*CCCl[, ncol(CCCl)]))%*%t(CCCl[, 1:ncol(G)]))*(CCCl[, ncol(CCCl)-1]*CCCl[, ncol(CCCl)]))
        }
      }
    }
    ############
    w_ <- w
    w_ <- w_[order(w_)]
    sumwi[i] <- sum(w_[1:length(w_)])
    if (i==1){
      W_f <- cbind(w, 1:length(w))
    }
    else{
      W_f <- rbind(W_f, cbind(w, 1:length(w)))
      #View(W_f)
    }
  }
  if (is.null(grid)){
    v1 <- sum(S)+sum(Si)
    v11 <- sum(S)+sum(Si)
    v2 <- sum(S2)
    nparmodel <- N-v11
    if (v11 < v2){
      v1 <- v11
    }
    res <- y-yhat
    rsqr1 <- t(res*wt)%*%res
    ym <- t(y*wt)%*%y
    rsqr2 <- ym-sum((y*wt)^2)/sum(wt)
    rsqr <- 1-rsqr1/rsqr2
    rsqradj <- 1-((N-1)/(N-v1))*(1-rsqr)
    sigma2 <- N*rsqr1/((N-v1)*sum(wt))
    root_mse <- sqrt(sigma2)
    measures <- c(sigma2, root_mse, v1, nparmodel, v2)
    influence <- as.vector(S)
    resstd <- res/(sqrt(sigma2)*sqrt(abs(1-influence)))
    CooksD <- resstd^2*influence/(v1*(1-influence))
    df <- N-(nvar+ncol(G))
    stdbi <- sqrt(abs(varbi))
    tstat <- bi[, 2]/stdbi
    probt <- 2*(1-pt(abs(tstat), df))
    malpha <- 0.05*(nvar/v1)
    if (model == "zip" || model == "zinb"){
      stdli <- sqrt(abs(varli))
      tstati <- li[, 2]

      for (j in 1:nrow(stdli)) {
        if (stdli[j] == 0){
          tstati[j] <- 0
        } else {
          tstati[j] <- li[j, 2]/stdli[j]
        }
      }
      tstati <- ifelse(is.na(tstati), 0, tstati)
      probti <- 2*(1-pt(abs(tstati), df))
      malpha <- 0.05*((nvar+ncol(G))/v1)
    }
    t_critical <- abs(qt(malpha/2, df))
    par_ <- 1/alphai[, 2]

    if (model == "zinb" || model == "zip") {
      if (any(lambda) == 0) {
        ll <- sum(-log(0+exp(pihat[pos0])) + log(0*exp(pihat[pos0]) + (par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0])) +
          sum(-log(0+exp(pihat[pos1])) + lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1] * log(par_[pos1]/(par_[pos1]+yhat2[pos1])))

        llnull1 <- sum(-log(1+zk[pos0]) + log(zk[pos0]+(par_[pos0]/(par_[pos0] + y[pos0]))^par_[pos0])) +
          sum(-log(1+zk[pos1])+lgamma(par_[pos1] + y[pos1])-lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(y[pos1]/(par_[pos1] + y[pos1])) + par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1])))

        llnull2 <- sum(-log(1+0) + log(0+(par_/(par_ + mean(y)))^par_)) +
          sum(-log(1+0) + lgamma(par_+y) - lgamma(y+1) - lgamma(par_) +
                y*log(mean(y)/(par_ + mean(y))) + par_*log(par_/(par_+mean(y))))
      } else {
        ll <- sum(-log(1+exp(pihat[pos0])) + log(exp(pihat[pos0]) + (par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0])) +
          sum(-log(1+exp(pihat[pos1])) + lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1])) + par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1])))

        llnull1 <- sum(-log(1+zk[pos0]) + log(zk[pos0] + (par_[pos0]/(par_[pos0] + y[pos0]))^par_[pos0])) +
          sum(-log(1 + zk[pos1]) + lgamma(par_[pos1] + y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1])) + par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1])))

        llnull2 <- sum(-log(1+0) + log(0+(mean(par_)/(mean(par_)+mean(y)))^mean(par_))) +
          sum(-log(1+0) + lgamma(mean(par_)+y) - lgamma(y+1) - lgamma(mean(par_)) +
                y*log(mean(y)/(mean(par_)+mean(y))) + mean(par_)*log(mean(par_)/(mean(par_) + mean(y))))
      }
      dev <- 2*(llnull1-ll)
      pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
      AIC <- 2*v1-2*ll
      AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
      adjpctll <- 1-(llnull1-ll+v1+0.5)/(llnull1-llnull2)

      if(model == "zinb"){
        AIC <- 2*(v1+v1/(ncol(x)+ncol(G)))-2*ll
        AICc <- AIC + 2*((v1+v1/(ncol(x)+ncol(G)))*((v1+v1/(ncol(x)+ncol(G)))+1)/ (N-(v1+v1/(ncol(x)+ncol(G))) - 1))
        adjpctll <- 1 - (llnull1-ll+(v1+v1/(ncol(x) + ncol(G))) + 0.5) / (llnull1 - llnull2)
      }
    }
    if(model == "poisson" || model == "negbin"){
      if (length(pos02)==0){
        pos0 <- pos1
        pos0x <- rep(1, length(pos1))
        pos0xl <- rep(1, length(pos1))
      } else {
        pos0x <- (par_[pos0]/(par_[pos0] + yhat[pos0]))^par_[pos0]
        pos0xl <- (par_[pos0]/(par_[pos0] + y[pos0]))^par_[pos0]
        pos0x <- ifelse(pos0x==0, E^-10, pos0x)
      }
      ll <- sum(-log(0+exp(pihat[pos0])) + log(0*exp(pihat[pos0]) + pos0x)) + sum(-log(0+exp(pihat[pos1])) +
                                                                                    lgamma(par_[pos1] + y[pos1]) - lgamma(y[pos1] + 1) - lgamma(par_[pos1]) + y[pos1]*log(yhat[pos1]/(par_[pos1] + yhat[pos1])) +
                                                                                    par_[pos1]*log(par_[pos1]/(par_[pos1] + yhat[pos1])))

      llnull1 <- sum(-log(1+zk) + log(zk+pos0xl)) + sum(-log(1+zk) +
                                                          lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1])+y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1])) +
                                                          par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1])))

      pos1xx <- 0+(par_[pos0]/(par_[pos0] + mean(y)))*par_[pos0]
      pos1xx <- ifelse(pos0x <= 0, E^-100, pos0x)

      llnull2 <- sum(-log(1+0) + log(pos1xx)) +
        sum(-log(1+0) + lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
              y[pos1]*log(mean(y)/(par_[pos1]+ mean(y))) + par_[pos1]*log(par_[pos1]/(par_[pos1]+ mean(y))))
      dev <- 2*(llnull1-ll)
      AIC <- -2*ll+2*v1
      AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
      pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
      adjpctll <- 1-(llnull1-ll+v1+0.5)/(llnull1-llnull2)
      if(model == "negbin"){
        AIC <- 2*(v1+v1/ncol(x))-2*ll
        AICc <- AIC+2*(v1+v1/ncol(x))*(v1+v1/ncol(x)+1)/(N-(v1+v1/ncol(x))-1)
        adjpctll <- 1-(llnull1-ll+v1+v1/ncol(x)+0.5)/(llnull1-llnull2)
      }
    }
    measures <- c(measures, dev, ll, pctll, adjpctll, AIC, AICc)
    names(measures) <- c("sigma2","root_mse", "n_params", "n_local_params",
                         "tot_variance", "deviance", "full_ll", "pct_ll", "adj_pct_ll", "AIC", "AICc")
    output <- append(output, list(measures))
    names(output)[length(output)] <- "measures"
    beta_ <- matrix(bi[, 2], nrow = N, byrow=T)
    beta2_ <- beta_
    if(model == "negbin" || model == "zinb"){
      alpha_= matrix(alphai[, 1:2], N)
      beta2_= cbind(beta_, alpha_[, 2])
    }
    qntl <- apply(beta2_, 2, quantile, c(0.25, 0.5, 0.75))
    IQR <- qntl[3, ]-qntl[1, ]
    qntl <- rbind(qntl, IQR)
    descriptb <- rbind(apply(beta2_, 2, "mean"), apply(beta2_, 2, "min"), apply(beta2_, 2, "max"))
    rownames(descriptb) <- c('Mean', 'Min', 'Max')
    if (model=='negbin' | model=="zinb"){
      colnames(qntl) <- c('Intercept', XVAR, 'alpha')
      colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
    }
    else{
      colnames(qntl) <- c('Intercept', XVAR)
      colnames(descriptb) <- c('Intercept', XVAR)
    }
    output <- append(output, list(qntl))
    names(output)[length(output)] <- "qntls_gwr_param_estimates"
    output <- append(output, list(descriptb))
    names(output)[length(output)] <- "descript_stats_gwr_param_estimates"
    stdbeta_ <- matrix(stdbi, N, byrow=TRUE)
    stdbeta2_ <- stdbeta_
    if (model=="negbin" | model=="zinb"){
      stdalpha_ <- matrix(alphai[, 3], N, byrow=TRUE)
      stdbeta2_ <- cbind(stdbeta_, stdalpha_)
    }
    qntls <- apply(stdbeta2_, 2, quantile, c(0.25, 0.5, 0.75))
    IQR <- qntls[3, ]-qntls[1, ]
    qntls <- rbind(qntls, IQR)
    descripts <- rbind(apply(stdbeta2_, 2, "mean"), apply(stdbeta2_, 2, "min"), apply(stdbeta2_, 2, "max"))
    rownames(descripts) <- c('Mean', 'Min', 'Max')
    if (model=='negbin' | model=="zinb"){
      colnames(qntls) <- c('Intercept', XVAR, 'alpha')
      colnames(descripts) <- c('Intercept', XVAR, 'alpha')
    }
    else{
      colnames(qntls) <- c('Intercept', XVAR)
      colnames(descripts) <- c('Intercept', XVAR)
    }
    t_test_gwr_param_estimates <- c(malpha, t_critical, df)
    names(t_test_gwr_param_estimates) <- c("p_value", "t_critical", "df")
    output <- append(output, list(t_test_gwr_param_estimates))
    names(output)[length(output)] <- "t_test_gwr_param_estimates"
    output <- append(output, list(qntls))
    names(output)[length(output)] <- "qntls_gwr_se"
    output <- append(output, list(descripts))
    names(output)[length(output)] <- "descript_stats_gwr_se"
    if (model=="zip" | model=="zinb"){
      lambda_ <- matrix(li[, 2], N, byrow=TRUE)
      lambda2_ <- lambda_
      qntl <- apply(lambda2_, 2, quantile, c(0.25, 0.5, 0.75))
      IQR <- qntl[3, ]-qntl[1, ]
      qntl <- rbind(qntl, IQR)
      descriptl <- rbind(apply(lambda2_, 2, "mean"), apply(lambda2_, 2, "min"), apply(lambda2_, 2, "max"))
      rownames(descriptl) <- c('Mean', 'Min', 'Max')
      if (int_inf){
        colnames(qntl) <- c('Intercept', xvarinf)
        colnames(descriptl) <- c('Intercept', xvarinf)
      }
      else{
        colnames(qntl) <- c(xvarinf)
        colnames(descriptl) <- c(xvarinf)
      }
      output <- append(output, list(qntl))
      names(output)[length(output)] <- "qntls_gwr_zero_infl_param_estimates"
      output <- append(output, list(descriptl))
      names(output)[length(output)] <- "descript_stats_gwr_zero_infl_param_estimates"
      stdlambda_ <- matrix(stdli, N, byrow=TRUE)
      stdlambda2_ <- stdlambda_
      qntls <- apply(stdlambda2_, 2, quantile, c(0.25, 0.5, 0.75))
      IQR <- qntls[3, ]-qntls[1, ]
      qntls <- rbind(qntls, IQR)
      descriptls <- rbind(apply(stdlambda2_, 2, "mean"), apply(stdlambda2_, 2, "min"), apply(stdlambda2_, 2, "max"))
      rownames(descriptls) <- c('Mean', 'Min', 'Max')
      if (int_inf){
        colnames(qntls) <- c('Intercept', xvarinf)
        colnames(descriptls) <- c('Intercept', xvarinf)
      }
      else{
        colnames(qntls) <- c(xvarinf)
        colnames(descriptls) <- c(xvarinf)
      }
      t_test_gwr_zero_infl_param_estimates <- c(malpha, t_critical, df)
      names(t_test_gwr_zero_infl_param_estimates) <- c("p_value", "t_critical", "df")
      output <- append(output, list(t_test_gwr_zero_infl_param_estimates))
      names(output)[length(output)] <- "t_test_gwr_zero_infl_param_estimates"
      #####
      output <- append(output, list(qntls))
      names(output)[length(output)] <- "qntls_gwr_zero_infl_se"
      output <- append(output, list(descriptls))
      names(output)[length(output)] <- "descript_stats_gwr_zero_infl_se"
    }
  }
  #### Non-Stationarity Test ####
  if (is.null(grid)){
    if(method!="adaptive_bsq"){
      BBk <- matrix(0, N, N)
      Vk <- rep(0, ncol(x))
      df1k <- rep(0, ncol(x))
      df2k <- rep(0, ncol(x))
      for (k in 1:ncol(x)){
        ek <- rep(0, ncol(x))
        ek[k] <- 1
        for (i in 1:N){
          m1 <- (i-1)*ncol(x)+1
          m2 <- m1+(ncol(x)-1)
          BBk[i, ] <- t(ek)%*%BB[m1:m2, ]
        }
        Vk[k] <- (t(y)*(1/N))%*%t(BBk)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBk%*%y
        df1k[k] <- sum(diag((1/N)*t(BBk)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBk))
        df2k[k] <- sum(diag(((1/N)*t(BBk)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBk)%*%((1/N)*t(BBk)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBk)))
      }
      Vk <- ifelse(abs(Vk)<=E^-8, 0, Vk)
      Fk <- (Vk/df1k)/sigma2
      ndf <- df1k^2/df2k
      ddf <- N-v1
      ddf <- rep(ddf, ncol(x))
      probf <- 1-pf(Fk, ndf, ddf)
      non_stat_test <- cbind(Vk, Fk, ndf, ddf, probf)
      rownames(non_stat_test) <- c('Intercept', XVAR)
      colnames(non_stat_test) <- c("V", "F", "ndf", "ddf", "Pr > F")
      output <- append(output, list(non_stat_test))
      names(output)[length(output)] <- "non_stationarity_test"
      if (model=="zip" | model=="zinb"){
        BBkl <- matrix(0, N, N)
        Vkl <- rep(0, ncol(G))
        df1kl <- rep(0, ncol(G))
        df2kl <- rep(0, ncol(G))
        for (k in 1:ncol(G)){
          ekl <- matrix(0, ncol(G), 1)
          ekl[k] <- 1
          for (i in 1:N){
            m1 <- (i-1)*ncol(G)+1
            m2 <- m1+(ncol(G)-1)
            BBkl[i, ] <- t(ekl)%*%BBl[m1:m2, ]
          }
          Vkl[k] <- (t(y)*(1/N))%*%t(BBkl)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBkl%*%y
          df1kl[k] <- sum(diag((1/N)*t(BBkl)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBkl))
          df2kl[k] <- sum(diag(((1/N)*t(BBkl)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBkl)%*%((1/N)*t(BBkl)%*%(diag(N)-(1/N)*matrix(1, N, N))%*%BBkl)))
        }
        Vkl <- ifelse(abs(Vkl)<=E^-8, 0, Vkl)
        Fkl <- (Vkl/df1kl)/sigma2
        ndfl <- df1kl^2/df2kl
        ddfl <- N-v1
        ddfl <- rep(ddfl, ncol(G))
        probfl <- 1-pf(Fkl, ndfl, ddfl)
        non_stat_test_zi <- cbind(Vkl, Fkl, ndfl, ddfl, probfl)
        if (int_inf){
          rownames(non_stat_test_zi) <- c('Intercept', xvarinf)
        }
        else{
          rownames(non_stat_test_zi) <- xvarinf
        }
        colnames(non_stat_test_zi) <- c("V", "F", "ndfl", "ddfl", "Pr > F")
        output <- append(output, list(non_stat_test_zi))
        names(output)[length(output)] <- "non_stationarity_test_zero_infl"
      }
    }
  }
  #### global estimates ####
  if (is.null(weight)){
    dfg <-N-nvar
  }
  if (model=="zinb"){
    b2 <- rbind(bg, alphag)
    if (alphag==E^-6){
      stdg <- c(stdabetalambdag[1:ncol(x)], (sqrt(1/abs(hessg))/(parg^2)))
    }
    else{
      stdg <- c(stdabetalambdag[2:(ncol(x)+1)], stdabetalambdag[1])
    }
    tg <- b2/stdg
    dfg <- length(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    lambdag <- lambdag
    if (alphag==E^-6){
      stdlambdag <- stdabetalambdag[(ncol(x)+1):length(stdabetalambdag)]
    }
    else{
      stdlambdag <- stdabetalambdag[(ncol(x)+2):length(stdabetalambdag)]
    }
    tlambdag <- lambdag/stdlambdag
    dflg <- length(y)-ncol(G)
    probtlambdag <- 2*(1-pt(abs(tlambdag), dflg))
    p <- ncol(x)+1+ncol(G)
  }
  if (model=="zip"){
    b2 <- bg
    stdg <- stdabetalambdag[1:ncol(x)]
    tg <- b2/stdg
    dfg <- length(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    lambdag <- lambdag
    stdlambdag <- stdabetalambdag[(ncol(x)+1):length(stdabetalambdag)]
    tlambdag <- lambdag/stdlambdag
    dflg <- length(y)-ncol(G)
    probtlambdag <- 2*(1-pt(abs(tlambdag), dflg))
    p <- ncol(x)+ncol(G)
  }
  if (model=="negbin"){
    b2 <- rbind(bg, alphag)
    if(alphag==E^-6){
      stdg <- c(stdabetalambdag[1:length(stdabetalambdag)], (sqrt(1/abs(hessg))/(parg^2)))
    }
    else{
      stdg <- c(stdabetalambdag[2:length(stdabetalambdag)], stdabetalambdag[1])
    }
    tg <- b2/stdg
    dfg <- length(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    p <- ncol(x)+1
  }
  if (model=="poisson"){
    b2 <- bg
    stdg <- stdabetalambdag
    tg <- b2/stdg
    dfg <- length(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    p <- ncol(x)
  }
  global_ests <- cbind(b2, stdg, tg, probtg)
  if (model=="negbin" | model=="zinb"){
    rownames(global_ests) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    rownames(global_ests) <- c('Intercept', XVAR)
  }
  colnames(global_ests) <- c("Par. Est.", "Std Error", "t Value", "Pr > |t|")
  output <- append(output, list(global_ests))
  names(output)[length(output)] <- "global_param_estimates"
  message("NOTE: The denominator degrees of freedom for the global model t tests is ", dfg, ".")
  if (model=="zip" | model=="zinb"){
    if (!is.null(xvarinf)){
      varnamezip1 <- xvarinf
      if (int_inf){
        varnamezip1 <- c("Intercept", xvarinf)
      }
    }
    else{
      if(int_inf){
        varnamezip1 <- "Intercept"
      }
    }
    analysis_max_like_zero_infl_param_estimates <- cbind(lambdag, stdlambdag, tlambdag, probtlambdag)
    rownames(analysis_max_like_zero_infl_param_estimates) <- varnamezip1
    colnames(analysis_max_like_zero_infl_param_estimates) <- c("Estimate", "Std Error", "t Value", "Pr > |t|")
    output <- append(output, list(analysis_max_like_zero_infl_param_estimates))
    names(output)[length(output)] <- "analysis_max_like_zero_infl_param_estimates"
    message("NOTE: The denominator degrees of freedom for the analysis of maximum likelihood zero inflation t tests is ", dflg, ".")
    ll <- sum(-log(1+exp(G[pos0, ]%*%lambdag))+
                log(exp(G[pos0, ]%*%lambdag)+
                      (parg/(parg+exp(x[pos0, ]%*%bg+Offset[pos0])))^parg))+
      sum(-log(1+exp(G[pos1, ]%*%lambdag))+
            lgamma(parg+y[pos1])-lgamma(y[pos1]+1)-
            lgamma(parg)+y[pos1]*log(exp(x[pos1, ]%*%bg+Offset[pos1])/(parg+exp(x[pos1, ]%*%bg+Offset[pos1])))+
            parg*log(parg/(parg+exp(x[pos1, ]%*%bg+Offset[pos1]))))
    AIC <- 2*p-2*ll
    AICc <- AIC+2*(p*(p+1)/(N-p-1))
    llnull1 <- sum(-log(1+zkg[pos0])+log(zkg[pos0]+
                                           (parg/(parg+y[pos0]))^parg))+
      sum(-log(1+zkg[pos1])+lgamma(parg+y[pos1])-
            lgamma(y[pos1]+1)-lgamma(parg)+
            y[pos1]*log(y[pos1]/(parg+y[pos1]))+
            parg*log(parg/(parg+y[pos1])))
    llnull2 <- sum(-log(1+0)+log(0+(parg/(parg+mean(y)))^parg))+
      sum(-log(1+0)+lgamma(parg+y)-lgamma(y+1)-lgamma(parg)+
            y*log(mean(y)/(parg+mean(y)))+parg*log(parg/(parg+mean(y))))
    devg <- 2*(llnull1-ll)
    pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
    adjpctll <- 1-(llnull1-ll+p+0.5)/(llnull1-llnull2)
    analysis_max_like_gof_measures <- c(devg, ll, pctll, adjpctll, AIC, AICc)
    names(analysis_max_like_gof_measures) <- c("deviance", "full_ll", "pct_ll", "adj_pct_ll", "AIC", "AICc")
    output <- append(output, list(analysis_max_like_gof_measures))
    names(output)[length(output)] <- "analysis_max_like_gof_measures"
  }
  if (model=="poisson" | model=="negbin"){
    yhatg <- exp(x%*%bg+Offset)
    if (length(pos02)==0){
      pos0 <- pos1
      pos0x <- 1
      pos0xl <- 1
    }
    else{
      pos0x <- (parg/(parg+exp(x[pos0, ]%*%bg+Offset[pos0])))^parg
      pos0xl <- (parg/(parg+y[pos0]))^parg
    }
    ll <- sum(-log(0+exp(G[pos0, ]%*%lambdag))+
                log(0*exp(G[pos0, ]%*%lambdag)+pos0x))+
      sum(-log(0+exp(G[pos1, ]%*%lambdag))+
            lgamma(parg+y[pos1])-lgamma(y[pos1]+1)-
            lgamma(parg)+y[pos1]*log(exp(x[pos1, ]%*%bg+
                                           Offset[pos1])/(parg+exp(x[pos1, ]%*%bg+
                                                                     Offset[pos1])))+
            parg*log(parg/(parg+exp(x[pos1, ]%*%bg+Offset[pos1]))))
    llnull1 <- sum(-log(1+zkg)+log(zkg+pos0xl))+
      sum(-log(1+zkg)+lgamma(parg+y[pos1])-
            lgamma(y[pos1]+1)-lgamma(parg)+
            y[pos1]*log(y[pos1]/(parg+y[pos1]))+
            parg*log(parg/(parg+y[pos1])))
    devg <- 2*(llnull1-ll)
    AIC <- -2*ll+2*nvar
    AICc <- -2*ll+2*nvar*(N/(N-nvar-1))
    tt2 <- y/mean(y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(y*log(tt2)-(y-mean(y)))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvar))*(1-pctdevg)
    if (model=="negbin"){
      AIC <- -2*ll+2*(nvar+1)
      AICc <- -2*ll+2*(nvar+1)*(N/(N-(nvar+1)-1))
      tt2 <- y/mean(y)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devnullg <- 2*sum(y*log(tt2)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*mean(y))))
      pctdevg <- 1-devg/devnullg
      adjpctdevg <- 1-((N-1)/(N-(nvar+1)))*(1-pctdevg)
    }
    analysis_max_like_gof_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(analysis_max_like_gof_measures) <- c("deviance", "full_ll", "pct_devg", "adj_pct_devg", "AIC", "AICc")
    output <- append(output, list(analysis_max_like_gof_measures))
    names(output)[length(output)] <- "analysis_max_like_gof_measures"
  }
  if (model=="zinb"){
    rownames(varcovg) <- c("alpha", "Intercept", XVAR, paste0(varnamezip1, "_xvarinf"))
    colnames(varcovg) <- c("alpha", "Intercept", XVAR, paste0(varnamezip1, "_xvarinf"))
  }
  else if (model=="negbin"){
    rownames(varcovg) <- c("alpha", "Intercept", XVAR)
    colnames(varcovg) <- c("alpha", "Intercept", XVAR)
  }
  else if (model=="zip"){
    rownames(varcovg) <- c("Intercept", XVAR, paste0(varnamezip1, "_xvarinf"))
    colnames(varcovg) <- c("Intercept", XVAR, paste0(varnamezip1, "_xvarinf"))
  }
  else{
    rownames(varcovg) <- c("Intercept", XVAR)
    colnames(varcovg) <- c("Intercept", XVAR)
  }
  output <- append(output, list(varcovg))
  names(output)[length(output)] <- "variance_covariance_matrix"
  ##################
  if (is.null(grid)){
    pihat <- exp(pihat)/(1+exp(pihat))
    res_ <- cbind(wt, y, yhat, res, resstd, influence, CooksD, sumwi, pihat)
    colnames(res_) <- c("wt", "y", "yhat", "res", "resstd", "influence", "cooksD", "sumwi", "pihat")
    output <- append(output, list(res_))
    names(output)[length(output)] <- "residuals"
    #View(res_)
    beta_out <- bi
    colnames(beta_out) <- c("id", "B", "x", "y")
    bistdt <- cbind(bi, stdbi, tstat, probt)
    parameters_ <- bistdt
    colnames(parameters_) <- c("id", "B", "x", "y", "stdbi", "tstat", "probt")
    tstat_ <- matrix(0, N, ncol(x))
    for (j in 1:nrow(stdbeta_)){
      for (k in 1:ncol(stdbeta_)){
        if (stdbeta_[j, k]==0){
          tstat_[j, k] <- 0
        }
        else{
          tstat_[j, k] <- (beta_[j, k])/stdbeta_[j, k]
        }
      }
    }
    tstat_ <- ifelse(is.na(tstat_), 0, tstat_)
    probt_ <- 2*(1-pt(abs(tstat_), df))
    sig_ <- matrix("not significant at 90%", N, ncol(x))
    for (i in 1:N){
      for (j in 1:ncol(x)){
        if (probt_[i, j]<0.01*(nvar/v1)){
          sig_[i, j] <- "significant at 99%"
        }
        else if (probt_[i, j]<0.05*(nvar/v1)){
          sig_[i, j] <- "significant at 95%"
        }
        else if (probt_[i, j]<0.1*(nvar/v1)){
          sig_[i, j] <- "significant at 90%"
        }
        else{
          sig_[i, j] <- "not significant at 90%"
        }
        if (is.na(probt_[i, j])){
          sig_[i, j] <- "not significant at 90%"
        }
      }
    }
    bistdt_ <- cbind(COORD, beta_, stdbeta_, tstat_, probt_)
    colname1 <- c("Intercept", XVAR)
    label_ <- c(rep("std_", ncol(x)), rep("tstat_", ncol(x)), rep("probt_", ncol(x)))
    colname <- c(c("x", "y"), colname1, paste0(label_, rep(colname1, 3)))
    if (model=="zip" | model=="zinb"){
      tstatl_ <- lambda_
      for (j in 1:nrow(stdlambda_)){
        for(k in 1:ncol(stdlambda_)){
          if (stdlambda_[j, k]==0){
            tstatl_[j, k] <- 0
          }
          else{
            tstatl_[j, k] <- lambda_[j, k]/stdlambda_[j, k]
          }
        }
      }
      probtl_ <- 2*(1-pt(abs(tstatl_), df))
      sigl_ <- matrix("not significant at 90%", N, ncol(G))
      for (i in 1:N){
        for (j in 1:ncol(G)){
          if (probtl_[i, j]<0.01*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 99%"
          }
          else if (probtl_[i, j]<0.05*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 95%"
          }
          else if (probtl_[i, j]<0.1*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 90%"
          }
          else{
            sigl_[i, j] <- "not significant at 90%"
          }
          if (is.na(probtl_[i, j])){
            sigl_[i, j] <- "not significant at 90%"
          }
        }
      }
      bistdt_ <- cbind(bistdt_, lambda_, stdlambda_, tstatl_, probtl_)
      colname2 <- xvarinf
      if (int_inf){
        colname2 <- c("Intercept", xvarinf)
      }
      label3_ <- rep("Inf_", ncol(G))
      colname3 <- paste0(label3_, colname2)
      label2_ <- c(rep("Inf_std_", ncol(G)), rep("Inf_tstat_", ncol(G)), rep("Inf_probt_", ncol(G)))
      colname <- c(c("x", "y"), colname1, paste0(label_, rep(colname1, 3)), colname3, paste0(label2_, rep(colname2, 3)))
      labell_ <- rep("sig_Inf_", ncol(G))
      colnamel <- paste0(labell_, colname2)
      sig_inf_parameters2 <- sigl_
      colnames(sig_inf_parameters2) <- colnamel
    }
    parameters2_ <- bistdt_
    colnames(parameters2_) <- colname
    label_ <- rep("sig_", ncol(x))
    colname_ <- paste0(label_, colname1)
    sig_parameters2 <- sig_
    colnames(sig_parameters2) <- colname_
    if (model=="negbin" | model=="zinb"){
      atstat <- alphai[, 2]
      for(j in 1:nrow(alphai)){
        if (alphai[j, 3]==0){
          atstat[j] <- 0
        }
        else{
          atstat[j] <- alphai[j, 2]/alphai[j, 3]
        }
      }
      atstat <- ifelse(is.na(atstat), 0, atstat)
      aprobtstat <- 2*(1-pnorm(abs(atstat)))
      siga_ <- matrix("not significant at 90%", N, 1)
      for (i in 1:N){
        if (aprobtstat[i]<0.01*(nvar/v1)){
          siga_[i] <- "significant at 99%"
        }
        else if (aprobtstat[i]<0.05*(nvar/v1)){
          siga_[i] <- "significant at 95%"
        }
        else if (aprobtstat[i]<0.1*(nvar/v1)){
          siga_[i] <- "significant at 90%"
        }
        else{
          siga_[i] <- "not significant at 90%"
        }
      }
      alphai <- cbind(alphai, atstat, aprobtstat)
      alpha_ <- alphai
      colnames(alpha_) <- c("id", "alpha", "std", "tstat", "probt")
      sig_alpha_ <- siga_
      colnames(sig_alpha_) <- "sig_alpha"
    }
  }
  else{
    beta_out <- bi
    colnames(beta_out) <- c("id", "B", "x", "y")
    stdbi <- sqrt(abs(varbi))
    tstat <- bi[, 2]/stdbi
    for (i in 1:nrow(tstat)){
      if (is.na(tstat[i])){
        tstat[i] <- 0
      }
    }
    bistdt <- cbind(bi, stdbi, tstat)
    parameters_ <- bistdt
    colnames(parameters_) <- c("id", "B", "x", "y", "stdbi", "tstat")
    if (model=="negbin" | model=="zinb"){
      atstat <- alphai[, 2]
      for(j in 1:nrow(alphai)){
        if (alphai[j, 3]==0){
          atstat[j] <- 0
        }
        else
          atstat[j] <- alphai[j, 2]/alphai[j, 3]
      }
      atstat <- ifelse(is.na(atstat), 0, atstat)
      aprobtstat <- 2*(1-pnorm(abs(atstat)))
      alphai <- cbind(POINTS, alphai, atstat, aprobtstat)
      alpha_ <- alphai
      colnames(alpha_) <- c("x", "y", "id", "alpha", "std", "tstat", "probt")
    }
    beta_ <- matrix(bi[, 2], mm, byrow=TRUE)
    stdbeta_ <- matrix(stdbi, mm, byrow=TRUE)
    tstat_ <- beta_/stdbeta_
    bistdt_ <- cbind(POINTS, beta_, stdbeta_, tstat_)
    colname1 <- c("Intercept", XVAR)
    label_ <- c(rep("std_", nvar), rep("tstat_", nvar))
    colname <- c(c("x", "y"), colname1, paste0(label_, rep(colname1, 2)))
    parameters_grid_ <- bistdt_
    colnames(parameters_grid_) <- colname
    output <- append(output, list(parameters_grid_))
    names(output)[length(output)] <- "param_estimates_grid"
    #View(parameters_grid_)
  }
  if (is.null(grid)){
    parameters2 <- cbind(parameters2_, sig_parameters2)
    if (model=="zip" | model=="zinb"){
      parameters2 <- cbind(parameters2, sig_inf_parameters2)
    }
    if (model=="negbin" | model=="zinb"){
      alpha_ <- cbind(alpha_, sig_alpha_)
      output <- append(output, list(alpha_))
      names(output)[length(output)] <- "alpha_estimates"
      #View(alpha_)
    }
    output <- append(output, list(parameters2))
    names(output)[length(output)] <- "parameter_estimates"
    #View(parameters2)
  }
  invisible(output)
}
