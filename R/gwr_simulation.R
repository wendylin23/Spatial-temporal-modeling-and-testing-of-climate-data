N = 10
p = 1/N
Nsim = 100
T.sim = 50
T.seq = c(1:T.sim - mean(1:T.sim))/10
# generate some points
set.seed(1234)
x.coord = 1:N
y.coord = 1:N
N.loc = N*N
points = expand.grid(x.coord,y.coord)

# distance matrix between points
Dd = as.matrix(dist(points))

# cor ~ dist regression
## temperature
pmean <- function(x,y) (x+y)/2
t=2
tmp.cor = as.vector(cor.list[[t]])
vec.dist = as.vector(distance)
dist.seq = seq(0, max(vec.dist), by=0.1)
tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist,
                     log.cor = log(tmp.cor),dist2 = vec.dist^2)
fit_lm1 = lm(I(tmp.cor - 1) ~ 0 + vec.dist, data = tmp.dat) 
beta_cor = coef(fit_lm1)[1]
sigma_cor = summary(fit_lm1)$sigma
rho_sim = array(0, dim = c(Nsim, dim(Dd)[1], dim(Dd)[2]))
y_sim = array(0, dim = c(Nsim, N.loc, T.sim))
set.seed(2022)
beta.sim = rnorm(N.loc, mean(beta_est[,t]), sd(beta_est[,t]))
b0.sim = rnorm(N.loc, mean(b0_est[,t]), sd(b0_est[,t]))
for(nsim in 1:Nsim)
{
  rho_t = 1 + beta_cor * Dd + matrix(rnorm(N.loc,0, sigma_cor),dim(Dd)[1], dim(Dd)[2])
  diag(rho_t) = 1
  rho_t[] = pmean(rho_t, matrix(rho_t, nrow(rho_t), byrow=TRUE))
  rho_1 = nearPD(rho_t,corr = TRUE, keepDiag = TRUE, maxit = 1000)
  rho_sim[nsim,,] = as.matrix(rho_1$mat)
  err_sim = mvrnorm(n = T.sim, rep(0,N.loc), rho_sim[nsim,,])
  for(tsim in 1:T.sim)
  {
    y_sim[nsim,,tsim] = as.vector(b0.sim + beta.sim * T.seq[tsim] + err_sim[tsim,])
  }
}

source(paste0(code_path,"gwr.R"))
d_bis = seq(1,5,by = 0.5)
w2 = gwr.lm(Dd,d_bis[4])
w0 = diag(1, N.loc)
fixed_beta = matrix(c(1,1,0), nrow = 3, ncol = 1)
slope_lm = slope_gwr = matrix(0, nrow = Nsim, ncol = N.loc)
b0_lm = b0_gwr = matrix(0, nrow = Nsim, ncol = N.loc)
slope_se_lm = slope_se_gwr = matrix(0, nrow = Nsim, ncol = N.loc)

for(nsim in 1:Nsim)
{
  y.w = array(0, dim = c(T.sim,1,N.loc))
  x.w = array(0, dim = c(T.sim,3,N.loc))
  for(i in 1:N.loc)
  {
    lm.fit = lm(y_sim[nsim, i, ] ~ T.seq)
    slope_lm[nsim, i] = summary(lm.fit)$coefficient[2,1]
    b0_lm[nsim, i] = summary(lm.fit)$coefficient[1,1]
    slope_se_lm[nsim, i] = summary(lm.fit)$coefficient[2,2]
    y.w[,1,i] = y_sim[nsim, i, ]
    x.w[,,i] = cbind(rep(1,T.sim),T.seq, T.seq)
  }
  gwr.fit = varw_fixed(y.w,p=0,xt = x.w,w = w2,fixed = fixed_beta)
  slope_gwr[nsim,] = gwr.fit$coef[2,1,]
  b0_gwr[nsim,] = gwr.fit$coef[1,1,]
  slope_se_gwr[nsim,] = gwr.fit$se.coef[2,1,]
}

## MSE of beta
mse.lm = apply(slope_lm, 2, function(x) mean((x - beta.sim)^2))
mse.gwr = apply(slope_gwr, 2, function(x) mean((x - beta.sim)^2))

## MSE of predicted values
pred.lm = pred.gwr = matrix(0, nrow = N.loc, ncol = T.sim)
y.mse.lm = y.mse.gwr = rep(0, Nsim)
for(nsim in 1:Nsim)
{
  for(i in 1:N.loc)
  {
    pred.lm[i,] = b0_lm[nsim, i] + slope_lm[nsim, i] * T.seq
    pred.gwr[i,] = b0_gwr[nsim, i] + slope_gwr[nsim, i] * T.seq
  }
  y.mse.lm[nsim] = mean((pred.lm - y_sim[nsim,,])^2)
  y.mse.gwr[nsim] = mean((pred.gwr - y_sim[nsim,,])^2)
}

## Zscore comparison
zscore.lm = apply(slope_lm/slope_se_lm, 2, mean)
zscore.gwr = apply(slope_gwr/slope_se_gwr, 2, mean)

sim_dat = list(y_sim = y_sim, b0.sim = b0.sim, beta.sim = beta.sim,
               slope_lm = slope_lm, slope_se_lm = slope_se_lm,
               slope_gwr = slope_gwr, slope_se_gwr = slope_se_gwr,
               b0_lm = b0_lm, b0_gwr = b0_gwr,
               mse.lm = mse.lm, mse.gwr = mse.gwr,
               y.mse.lm = y.mse.lm, y.mse.gwr = y.mse.gwr,
               zscore.lm = zscore.lm, zscore.gwr = zscore.gwr)

save(sim_dat, file = paste0("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/simulation/sim_dat_N",N,"_Nsim",Nsim,"_T",T.sim,".rdata"))

## precipitation

gwr_gamma = function(b0,b1, logshape) {
  loglik = rep(0, N.loc)
  for(j in 1:N.loc)
  {
    linear_predictor = x[[j]][,1]*b0 + x[[j]][,2]*b1
    sum_j = sum((exp(logshape)-1) *log(y[j,]) - log(gamma(exp(logshape))) - 
                  exp(logshape)*(linear_predictor - logshape) - exp(logshape) * y[j,]/exp(linear_predictor))
    loglik[j] = w_i[j] * sum_j
  }
  # rate = shape / mean:
  # sum of negative log likelihoods:
  -sum(loglik)
}

gr_gamma = function(b0,b1, logshape) {
  loglik = rep(0, N.loc)
  for(j in 1:N.loc)
  {
    linear_predictor = x[[j]][,1]*b0 + x[[j]][,2]*b1
    sum_j = sum((exp(logshape)-1) *log(y[j,]) - log(gamma(exp(logshape))) - 
                  exp(logshape)*(linear_predictor - logshape) - exp(logshape) * y[j,]/exp(linear_predictor))
    loglik[j] = w_i[j] * sum_j
  }
  # rate = shape / mean:
  # sum of negative log likelihoods:
  -sum(loglik)
}

pmean <- function(x,y) (x+y)/2
t=2
tmp.cor = as.vector(cor.glm[[t]])
vec.dist = as.vector(distance)
dist.seq = seq(0, max(vec.dist), by=0.1)
tmp.dat = data.frame(tmp.cor = tmp.cor, vec.dist = vec.dist,
                     log.cor = log(tmp.cor),dist2 = vec.dist^2)
fit_lm1 = lm(I(tmp.cor - 1) ~ 0 + vec.dist, data = tmp.dat) 
beta_cor = coef(fit_lm1)[1]
sigma_cor = summary(fit_lm1)$sigma
rho_sim = array(0, dim = c(Nsim, dim(Dd)[1], dim(Dd)[2]))
y_sim = array(0, dim = c(Nsim, N.loc, T.sim))
set.seed(2022)
beta.sim = rnorm(N.loc, mean(slope_est[,t]), sd(slope_est[,t]))
b0.sim = rnorm(N.loc, mean(scale_est[,t]), sd(scale_est[,t]))
for(nsim in 1:Nsim)
{
  rho_t = 1 + beta_cor * Dd + matrix(rnorm(N.loc,0, sigma_cor),dim(Dd)[1], dim(Dd)[2])
  diag(rho_t) = 1
  rho_t[] = pmean(rho_t, matrix(rho_t, nrow(rho_t), byrow=TRUE))
  rho_1 = nearPD(rho_t,corr = TRUE, keepDiag = TRUE, maxit = 1000)
  rho_sim[nsim,,] = as.matrix(rho_1$mat)
  err_sim = mvrnorm(n = T.sim, rep(0,N.loc), rho_sim[nsim,,])
  for(tsim in 1:T.sim)
  {
    y_sim[nsim,,tsim] = exp(b0.sim + beta.sim * T.seq[tsim] + err_sim[tsim,])
  }
}

source(paste0(code_path,"gwr.R"))
d_bis = seq(1,5,by = 0.5)
w2 = gwr.lm(Dd,d_bis[4])
w0 = diag(1, N.loc)
fixed_beta = matrix(c(1,1,0), nrow = 3, ncol = 1)
slope_lm = slope_gwr = matrix(0, nrow = Nsim, ncol = N.loc)
b0_lm = b0_gwr = matrix(0, nrow = Nsim, ncol = N.loc)
slope_se_lm = slope_se_gwr = matrix(0, nrow = Nsim, ncol = N.loc)
pred.lm = y_sim

for(nsim in 1:Nsim)
{
  x = list()
  for(i in 1:N.loc)
  {
    x[[i]] = cbind(rep(1, T.sim), T.seq)
  }
  
  y = y_sim[nsim,,]
  
  for(i in 1:N.loc)
  {
    w_i = w2[i,]
    #m_mle2 = bbmle::mle2(gwr_gamma,optimfun = "BFGS",
    #                     start = list(b0 = b0.sim[nsim],b1 = beta.sim[nsim],
    #                                  logshape = 0))
    m_mle = glm(y[i,] ~ T.seq, family=Gamma(link = "log"))
    #slope_gwr[nsim,i] = bbmle::coef(m_mle2)[2]
    #b0_gwr[nsim,i] = bbmle::coef(m_mle2)[1]
    #slope_se_gwr[nsim,i] = sqrt(diag(bbmle::vcov(m_mle2)))[2]
    slope_lm[nsim,i] = m_mle$coefficients[2]
    slope_se_lm[nsim,i] = summary(m_mle)$coefficients[2,2]
    b0_lm[nsim,i] = 1/m_mle$coefficients[1]
    pred.lm[nsim,i,] = exp(predict(m_mle))
  }
}

## MSE of beta
mse.lm = apply(slope_lm, 2, function(x) mean((x - beta.sim)^2))
mse.gwr = apply(slope_gwr, 2, function(x) mean((x - beta.sim)^2))

## MSE of predicted values
pred.gwr = matrix(0, nrow = N.loc, ncol = T.sim)
y.mse.lm = y.mse.gwr = rep(0, Nsim)
for(nsim in 1:Nsim)
{
  for(i in 1:N.loc)
  {
    #pred.lm[i,] = exp(b0_lm[nsim, i] + slope_lm[nsim, i] * T.seq)
    pred.gwr[i,] = exp(b0_gwr[nsim, i] + slope_gwr[nsim, i] * T.seq)
  }
  y.mse.lm[nsim] = mean((pred.lm[nsim,,] - y_sim[nsim,,])^2)
  y.mse.gwr[nsim] = mean((pred.gwr - y_sim[nsim,,])^2)
}

## Zscore comparison
zscore.lm = apply(slope_lm/slope_se_lm, 2, mean)
zscore.gwr = apply(slope_gwr/slope_se_gwr, 2, mean)

sim_dat = list(y_sim = y_sim, b0.sim = b0.sim, beta.sim = beta.sim,
               slope_lm = slope_lm, slope_se_lm = slope_se_lm,
               slope_gwr = slope_gwr, slope_se_gwr = slope_se_gwr,
               b0_lm = b0_lm, b0_gwr = b0_gwr,
               mse.lm = mse.lm, mse.gwr = mse.gwr,
               y.mse.lm = y.mse.lm, y.mse.gwr = y.mse.gwr,
               zscore.lm = zscore.lm, zscore.gwr = zscore.gwr)

#save(sim_dat, file = paste0("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/simulation/sim_temp_dat_N",N,"_Nsim",Nsim,"_T",T.sim,".rdata"))

## Comparison figures
load(paste0("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/simulation/sim_pr_dat_N",N,"_Nsim",Nsim,"_T",T.sim,".rdata"))
pr_sim = sim_dat
load(paste0("/Users/wenyilin/Dropbox/UCSD/Thesis/3.Precipitation/Code/results/simulation/sim_temp_dat_N",N,"_Nsim",Nsim,"_T",T.sim,".rdata"))
temp_sim = sim_dat

### correlation plot of one simulation
sim.cor = as.vector(rho_sim[1,,])
sim.dist = as.vector(Dd)
plot(sim.dist, sim.cor,pch=".", cex=0.8,
     main = "", xlab = "distance",ylab = "correlation")

### coefficient plot (precipitation)
autoimage(points$Var1, points$Var2,
          pr_sim$beta.sim,
          interp.args = list(no.X = 200, no.Y = 200),
          ylab = "",xlab = "",
          #main = "Simulated Slopes",
          #outer.title = "Estimated slopes (C/decade) of four seasons",
          #mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          #legend.axis.args = list(cex.axis=1),
          lratio = 0.15,
          col=cmap(200),
          zlim = c(-0.2,0.1),
          legend = "vertical"
          )

### coefficient plot (temperature)
autoimage(points$Var1, points$Var2,
          temp_sim$beta.sim,
          interp.args = list(no.X = 200, no.Y = 200),
          ylab = "",xlab = "",
          #main = "Simulated Slopes",
          #outer.title = "Estimated slopes (C/decade) of four seasons",
          #mtext.args = list(cex = 1),
          #text = as,
          #text.args = list(cex=1),
          #legend.axis.args = list(cex.axis=1),
          lratio = 0.15,
          col=cmap(200),
          zlim = c(0,0.6),
          legend = "vertical"
)

### boxplot of slope parameters
mse.slope.res = cbind(apply(temp_sim$slope_lm, 2, function(x) mean((x - temp_sim$beta.sim)^2)),
                 apply(temp_sim$slope_gwr, 2, function(x) mean((x - temp_sim$beta.sim)^2)),
                 apply(pr_sim$slope_lm, 2, function(x) mean((x - pr_sim$beta.sim)^2)),
                 apply(pr_sim$slope_gwr, 2, function(x) mean((x - pr_sim$beta.sim)^2)))

b=boxplot(mse.slope.res, col=c(alpha("darkred",0.6),alpha("darkred",0.6),
                       alpha("darkblue",0.6),alpha("darkblue",0.6)),
        pch=16,cex=0.8,frame=FALSE, xaxt = "n", ylim = c(0,0.06), ylab = "MSE",
        names=c(expression(paste("OLS",~beta[1])),
                         expression(paste("GWR",~beta[1])),
                         expression(paste("GR",~beta[1])),
                         expression(paste("GWGR",~beta[1]))))
axis(side = 1, at = seq_along(b$names), labels = b$names,tick = FALSE)
legend("topright",legend = c("temperature","precipitation"), pch=16,
       col=c(alpha("darkred",0.6),alpha("darkblue",0.6)))

### boxplot of predicted MSE
mse.y.res = cbind(temp_sim$y.mse.lm,temp_sim$y.mse.gwr,
                  pr_sim$y.mse.lm,pr_sim$y.mse.gwr)

par(mfrow = c(1,2))
par(mar=c(4,4,2,0))
b = boxplot(mse.y.res[,1:2], col=c(alpha("darkred",0.6),alpha("darkred",0.6)),
        pch=16,cex=0.8,frame=FALSE, ylim = c(0,2),xaxt = "n", yaxt = "n", ylab = "MSE",
        names=c("OLS","GWR"))
axis(side = 1, at = seq_along(b$names), labels = b$names,tick = FALSE)
axis(side = 2, at = seq(0,2,length = 5), labels = seq(0,2,length = 5), col.axis="darkred")

par(mar=c(4,0,2,2))
b = boxplot(mse.y.res[,3:4], col=c(alpha("darkblue",0.6),alpha("darkblue",0.6)),
        pch=16,cex=0.8,frame=F, ylim = c(0,15),xaxt = "n", yaxt = "n",
        names=c("GR","GWGR"))
axis(side = 1, at = seq_along(b$names), labels = b$names, tick = FALSE)
axis(side = 4, at = seq(0,15,length = 6), labels = seq(0,15,length = 6), col.axis="blue")
legend("topright",legend = c("temperature","precipitation"), pch=16,
       col=c(alpha("darkred",0.6),alpha("darkblue",0.6)))

### Zscore comparison
col1 = rgb(1,0,0,1/4)
col2 = rgb(0,1,0,1/4)
zscore.lm = apply(temp_sim$slope_lm/temp_sim$slope_se_lm, 2, mean)
zscore.gwr = apply(temp_sim$slope_gwr/temp_sim$slope_se_gwr, 2, mean)
tmp.z = cbind(zscore.lm,zscore.gwr)
hist(tmp.z[,1],xlab = "Z-score",main = "",
     probability = TRUE, ylim=c(0,1.5), breaks = 20,
     col=col1,xlim = c(1,6))
#hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
hist(tmp.z[,2],col=col2, add = T,breaks = 20,
     probability = TRUE)
legend("topright",legend = c("OLS","GWR"),
       fill = c(col1,col2),
       bty = 'n',cex = 1,
       border = NA)

zscore.lm = apply(pr_sim$slope_lm/pr_sim$slope_se_lm, 2, mean)
zscore.gwr = apply(pr_sim$slope_gwr/pr_sim$slope_se_gwr, 2, mean)
tmp.z = cbind(zscore.lm,zscore.gwr)
hist(tmp.z[,1],xlab = "Z-score",main = "",
     probability = TRUE, ylim=c(0,1.5), breaks = 20,
     col=col1,xlim = c(-4,1))
#hist(z_score$slope_m1[,2],col=rgb(0,1,1,1/4), add = T, breaks =30)
hist(tmp.z[,2],col=col2, add = T,breaks = 20,
     probability = TRUE)
legend("topright",legend = c("GR","GWGR"),
       fill = c(col1,col2),
       bty = 'n',cex = 1,
       border = NA)


# weights matrix
w <- exp(-p * Dd)
w = gwr.lm(Dd,3)
Ww = spam::chol(w)

# errors
z <- t(Ww) %*% rnorm(N,0,1) 

# plot
df <- data.frame(x = x.coord, y = y.coord, z = z)
require(ggplot2)
ggplot(df, aes(x = x, y = y, col = z)) +
  geom_point() +
  scale_colour_gradient(low="red", high="white")
