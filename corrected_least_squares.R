library("ncvreg")
library("logistf")

corrected_least_squares<-function(X, Y, etaHat){
  
  # A FUNCTION TO COMPUTE THE CORRECTED LEAST SQUARES ESTIMATOR OF Lewis & Battey (2024)
  # https://academic.oup.com/biomet/advance-article/doi/10.1093/biomet/asad065/7338235?login=true
  
  # Inputs:
  # X is the n x p design matrix (rows are observations)
  # etaHat is an nx1 matrix giving the estimate of X%*%betaStar, eg: from a LASSO regression
  # Y is the nx1 response vector (taking values in {-1,1})
  
  # Output: a data frame. The first column gives the estimates and the second gives the standard 
  # errors based on an approximate Gaussian distribution. Note that no guarantees have 
  # been provided for the standard errors to date.
  
  XTX=t(X)%*%X
  invXTX=solve(XTX)
  Px=X%*%invXTX%*%t(X)
  
  #obtain least-squares estimate
  betaHat0=invXTX%*%t(X)%*%Y
  
  #estimate c
  cHat=t(etaHat)%*%tanh(etaHat/2)/norm(etaHat, type="2")^2
  
  #estimate delta
  Phat=Px-etaHat%*%t(etaHat)/norm(etaHat, type="2")^2
  deltaHat=invXTX%*%t(X)%*%Phat%*%tanh(etaHat/2)
  
  #compute correction to least squares
  betaTildeStar_OLS=(betaHat0-deltaHat)/as.vector(cHat)
  
  #estimate of asymptotic variance (no guarantees)
  GammaMat=diag(as.vector(1-tanh(etaHat/2)^2))
  se_betaTildeStar=sqrt(diag(invXTX%*%t(X)%*%GammaMat%*%X%*%invXTX))/as.vector(cHat)
  
  
  return(data.frame(CLS_estimates=betaTildeStar_OLS, CLS_standard_errors=se_betaTildeStar))
  
}

#An example
set.seed(123)
n=700 
p=0.1*n
s=3 
sigStrength=0.9
rho=0.5
sigma=rho*matrix(1,p,p)+(1-rho)*diag(1,p,p)
sigma[1,2]=0.95
sigma[2,1]=0.95
betaStar=sigStrength*rbind(matrix(1, s,1), matrix(0, p-s,1))
betaStar[2]=-betaStar[2]

#generate design matrix
Xn1=mvtnorm::rmvnorm(n,matrix(0,p,1), sigma)
X=Xn1

#generate response vector
logisticProb=exp(X%*%betaStar)/(1+exp(X%*%betaStar))
Y01=rbinom(n,1,logisticProb)
Y=2*Y01-1

#estimate X*betaStar using a SCAD penalised regression
scadRegress=cv.ncvreg(Xn1, Y01, family="binomial", penalty="SCAD")
etaHat=as.vector(coef(scadRegress))[1]+X%*%as.vector(coef(scadRegress))[-1]

CLS=corrected_least_squares(X, Y, etaHat)

#Apply Firth's estimator to this dataset
firthRegress=logistf::logistf(Y01~Xn1-1)
betaHatStar_Firth=t(t(coef(firthRegress)))

results=data.frame(betaStar=betaStar, SCAD=coef(scadRegress)[-1],Firth=betaHatStar_Firth,  CLS)
results
