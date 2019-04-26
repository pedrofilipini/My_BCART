###############################
###### Generating Data ########
###############################

set.seed(2019)

n <- 1000
x1 <- rnorm(n,5)
x2 <- rnorm(n,4)
x3 <- rnorm(n,3)

y <- x1+x2+x3+rnorm(n)



####################
#Auxiliar Functions#
####################

.cp_quantile = function(x, num=100, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))
  
  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }
  
  return(ret)
}

findfinal <- function(x){
  for(i in 1:length(x)){
    if((x[[i]][2]==0)&&(x[[i]][3]==0))  
      x[[i]][5] <- 1
    else{
      x[[i]][5] <- 0
    }
  }
  return(x)
}

is.root <- function(x){
  if(length(x)==1){
    return(1)
  }else{
    return(0)
  }
}

findid <- function(id, x){
  for(i in 1:length(x)){
    if(x[[i]][1]==id){
      return(i)
    }
  }
  return(0)
}

findparent <- function(x){
  for(i in 1:length(x)){
    aux1 <- findid((x[[i]][1]*2+1),x)
    aux2 <- findid((x[[i]][1]*2),x)
    if((aux1!=0)&&(aux2!=0))  
      if((x[[aux1]][2]==0)&&
         (x[[aux1]][3]==0)&&
         (x[[aux2]][2]==0)&&
         (x[[aux2]][3]==0)){
        x[[i]][5] <- 1
      }else{
        x[[i]][5] <- 0
      }
    else{
      x[[i]][5] <- 0
    }
  }
  return(x)
}

loglike <- function(y_suff, sigma_y, sigma_mu, mu_bar){
  (-nrow(y_suff)/2)*(2*pi*sigma_y^2) + 
    1/2*(sigma_y^2/(sigma_y^2 + nrow(y_suff)*sigma_mu^2)) -
    sum((y_suff$y-mean(y_suff$y))^2)/(2*sigma_y^2) - 
    1/2*(nrow(y_suff)*((mean(y_suff$y)-mu_bar)^2))/(sigma_y^2 + nrow(y_suff)*sigma_mu^2)
}

rss <- function(y_suff, tr){
  sum((y_suff$y-tr[4])^2)
}

depth <- function(id){
  aux <- 0
  while(id!=1){
    aux <- aux+1
    id <- floor(id/2)
  }
  return(aux)
}

split <- function(alpha, beta, d){
  alpha*(1+d)^(-beta)
}

can.split <- function(tr, y_suff){
  x <- findfinal(tr) #find the terminal nodes
  
  index <- numeric(0) #find the terminal nodes indexes
  
  for(i in 1:length(x)){
    if(x[[i]][5]==1){ #if node is terminal
      index <- c(index,i) #save index in the vector
    }
  }
  
  if(sum(sapply(y_suff[index],nrow)<15)==length(sapply(y_suff[index],nrow))){
    return(0)
  }else{
    return(1)
  }
}

mu_node_draw <- function(sigma_y, sigma_mu, mu_bar, ysuff){
  mu_new <- (sigma_y^2*mu_bar + sigma_mu^2*sum(ysuff$y))/(sigma_y^2+sigma_mu^2*nrow(ysuff))
  sigma_new <- sqrt((sigma_y^2*sigma_mu^2)/(sigma_y^2+sigma_mu^2*nrow(ysuff)))
  draw <- rnorm(1, mean = mu_new, sd = sigma_new)
  return(draw)
}


###########################
#End of Auxiliar Functions#
###########################




#################
#Hyperparameters#
#################
p <- 3 #number of covariates
c <- 100 #number of cuts
cuts <- lapply(1:ncol(cbind(x1,x2,x3)), function(i) .cp_quantile(cbind(x1,x2,x3)[,i], num = c)) #the cuts
x_matrix <- data.frame(x1,x2,x3,y) #matrix of values

pgrow <- 0.5 #P(GROW)
pprune <- 0.5 #P(PRUNE)

sigma_y <- sd(y)  #Initial guess for the MCMC

sigma_mu <- sigma_y/2  #hyperparameter of \mu_i
mu_bar <- mean(y) #hyperparameter of \mu_i

nu <- 3 #hyperparameter of \sigma_i
sigquant <- 0.90
qchi = qchisq(1.0-sigquant,nu)
lambda = (sd(y)*sd(y)*qchi)/nu #lambda parameter for sigma prior

alpha <- 0.95 #hyperparameter of P(Split)
beta <- 2 #hyperparameter of P(Split)

npost <- 12000 #size of posterior
burn <- 2000
thin <- 10
treedraws <- vector("list", length = npost)

########################
#End of Hyperparameters#
########################




#Only nodes with, at least, 15 observations will be splitted (functon can.split)
#Minimum of observations at terminal nodes is 5 (ok_cut variable)





#######################
######## MCMC #########
#######################
set.seed(2019)
begin <- proc.time()
for(j in 1:npost){

if(j==1){
  print("Starting MCMC:")
}
  
sigma_mu <- sigma_y[j]/2  #hyperparameter of \mu_i

if(j==1){ #If it is the first iteration of the MCMC
  tr <- list(numeric(0)) #Start the List
  tr[[1]] <- c(1,0,0,mean(y)) #Start the Root Tree
  
  y_suff <- list(data.frame(0)) #Start the Y List
  y_suff[[1]] <- x_matrix #Start the first list
  
  n_terminal <- numeric(0) #vector of terminal nodes counter
  accept_ratio <- numeric(0) #acceptance ratio counter
}


k <- runif(1) #select one of the proposals randomly

if((k<=pgrow || is.root(tr)==1) && (can.split(tr, y_suff)==1)){ #if GROW is selected
  
  ok_cut <- 0
  
  while(ok_cut==0){
  
  y_suff_new <- y_suff
  trnew <- tr #creating the new tree
  x <- findfinal(trnew) #find the terminal nodes
  
  b <- 0 #terminal nodes counter
  index <- numeric(0) #find the terminal nodes indexes
  
  for(i in 1:length(x)){
  
      b = b + x[[i]][5] #counter of terminal nodes
    
      if(x[[i]][5]==1){ #if node is terminal
      index <- c(index,i) #save index in the vector
      }
  }
  
  if(length(index)==1){
    base <- index #sample from the available indexes
  }else{
    base <- sample(index,1)  #sample from the available indexes
  }
  
  while(nrow(y_suff[[base]])<15){ #assuring enough data to work with
    base <- sample(index,1)  #sample from the available indexes
  }
    
  y_suff_split <- y_suff[[base]] #creating new list of Y
  
  if(p==1){
    trnew[[base]][2] <- p
  }else{
    trnew[[base]][2] <- sample(1:p, 1) #sample from available covariates
  }
  
  if(c==1){
    trnew[[base]][3] <- c
  }else{
    base_aux <- base #auxiliar variable
    n_j_adj <- cuts[[(trnew[[base_aux]][2])]] #possible cuts
    
    while(base_aux!=1){ #While not in root
      par <- base_aux%%2==0 #verify if it is left or right
      base_aux <- floor(base_aux/2) #Climb the tree
      if(trnew[[base_aux]][2]==trnew[[base]][2]){ #if is the same variable selected
        if(par==1){
          n_j_adj<=trnew[[base_aux]][3] #updating possible cuts
        }else if(par==0){
          n_j_adj>trnew[[base_aux]][3]  #updating possible cuts
        }
      }
    }
    
    cmin <- which(min(n_j_adj)==cuts[[(trnew[[base]][2])]])
    cmax <- which(max(n_j_adj)==cuts[[(trnew[[base]][2])]])
    n_j_adj <- cmax-cmin+1
    
    trnew[[base]][3] <- sample(cmin:cmax, 1) #sample from available cuts
  }
  
  var1 <- 0 #set terminal node config
  var2 <- 0 #set terminal node config
  cut1 <- 0 #set terminal node config
  cut2 <- 0 #set terminal node config
  y1 <- 0 #Just for having the space to sample later
  y2 <- 0 #Just for having the space to sample later
  id1 <- x[[base]][1]*2  #set terminal node new id
  id2 <- x[[base]][1]*2 + 1  #set terminal node new id
  
  node1 <- c(id1, var1, cut1, y1) #new node (left)
  node2 <- c(id2, var2, cut2, y2) #new node (right)
  
  trnew[[(length(tr)+1)]] <- node1 #plugging it on the new tree
  trnew[[(length(tr)+2)]] <- node2 #plugging it on the new tree
  
  rule <- cuts[[(trnew[[base]][2])]][((trnew[[base]][3]))] #find specific rule
  
  y_suff_split$rule <- y_suff_split[,(trnew[[base]][2])]<rule #set the rule
  
  if((nrow(y_suff_split)<15)||(sum(y_suff_split$rule)<5)||(sum(!y_suff_split$rule)<5)||(sum(y_suff_split$rule)==nrow(y_suff_split))){
    ok_cut <- 0
  }else{
    ok_cut <- 1
  }
  
  }
  y1 <- y_suff_split[y_suff_split$rule==TRUE,] #splitting Y
  y2 <- y_suff_split[y_suff_split$rule==FALSE,] #splitting Y
  
  y_suff_new[[(length(y_suff)+1)]] <- y1 #plugging it on the new tree
  y_suff_new[[(length(y_suff)+2)]] <- y2 #plugging it on the new tree
  
  
  trnew[[(length(tr)+1)]][4] <- mu_node_draw(sigma_y[j], sigma_mu, mu_bar, y_suff_new[[(length(tr)+1)]]) #Sample of posterior \mu_i
  trnew[[(length(tr)+2)]][4] <- mu_node_draw(sigma_y[j], sigma_mu, mu_bar, y_suff_new[[(length(tr)+2)]]) #Sample of posterior \mu_i
  
  
  ########################
  #Calculating the ratios#
  ########################
  
  z <- findparent(trnew) #find possible pruning rules
  w <- 0 #pruning possibilites counter
  
  for(i in 1:length(z)){
    w = w + z[[i]][5]  #pruning possibilites counter
  }
  
  p1 <- log((pprune/pgrow)*b*p*length(n_j_adj)/w)
  
  #Calculating the likelihood
  e1 <- loglike(y_suff[[base]], sigma_y[j], sigma_mu, mu_bar)
  e2 <- loglike(y_suff_new[[(findid(trnew[[base]][1]*2,trnew))]], sigma_y[j], sigma_mu, mu_bar)
  e3 <- loglike(y_suff_new[[((findid((trnew[[base]][1]*2+1),trnew)))]], sigma_y[j], sigma_mu, mu_bar)
  
  p2 <- (e2+e3)-e1
  
  #Calculating Tree prior
  p3 <- log(1-(split(alpha, beta, depth(trnew[[((findid((trnew[[base]][1]*2),trnew)))]][1])))) +
          log(1-(split(alpha, beta, depth(trnew[[(((findid((trnew[[base]][1]*2+1),trnew))))]][1])))) +
            log((split(alpha, beta, depth(trnew[[(base)]][1])))) +
              log(p*length(n_j_adj)) -
                log(1-(split(alpha, beta, depth(trnew[[(base)]][1]))))
  
  r <- exp(p1+p2+p3)
  #print(r)
  
}else if(k<=(pgrow+pprune) || (can.split(tr, y_suff)==0)){ #if PRUNE is selected
  trnew <- tr
  x <- findparent(trnew) #find possible pruning rules
  
  w <- 0 #pruning possibilites counter
  index <- numeric(0) #pruning possibilites index
  
  for(i in 1:length(x)){
    
    w = w + x[[i]][5]  #pruning possibilites counter
    
    if(x[[i]][5]==1){
      index <- c(index,i)  #pruning possibilites vector of indexes
    }
  }
  
  if(length(index)==1){
    base <- index
  }else{
    base <- sample(index,1)  #sample from pruning possibilites indexes
  }
  trnew[[base]][2] <- 0 #set prunned node as terminal
  trnew[[base]][3] <- 0 #set prunned node as terminal
  
  new_final_id <- x[[base]][1]
  ind1 <- findid(new_final_id*2, x) #find indexes of old terminal nodes
  ind2 <- findid((new_final_id*2+1), x) #find indexes of old terminal nodes
  
  trnew <- trnew[-c(ind1,ind2)] #remove old terminal nodes
  y_suff_new <- y_suff[-c(ind1,ind2)] #remove old terminal nodes

  ########################
  #Calculating the ratios#
  ########################
  
  z <- findfinal(trnew) #find the terminal nodes
  b <- 0 #terminal nodes counter

  for(i in 1:length(z)){
    b = b + z[[i]][5] #counter of terminal nodes
  }
  
  base_aux <- base
  n_j_adj <- cuts[[(tr[[base_aux]][2])]]
    
  while(base_aux!=1){ #while is not a root
    par <- base_aux%%2==0 #verify if it is left or right
    base_aux <- floor(base_aux/2) #Climb the tree
    if(tr[[base_aux]][2]==tr[[base]][2]){
      if(par==1){
        n_j_adj<=tr[[base_aux]][3] #updating possible cuts
      }else if(par==0){
        n_j_adj>tr[[base_aux]][3] #updating possible cuts
      }
    }
  }
  
  cmin <- which(min(n_j_adj)==cuts[[(tr[[base]][2])]])
  cmax <- which(max(n_j_adj)==cuts[[(tr[[base]][2])]])
  n_j_adj <- cmax-cmin+1
  
  p1 <- log((pgrow/pprune)*w/((b-1)*p*length(n_j_adj)))
  
  #Calculating the likelihood
  e1 <- loglike(y_suff_new[[base]], sigma_y[j], sigma_mu, mu_bar)
  e2 <- loglike(y_suff[[(((findid((tr[[base]][1]*2),tr))))]], sigma_y[j], sigma_mu, mu_bar)
  e3 <- loglike(y_suff[[(((findid((tr[[base]][1]*2+1),tr))))]], sigma_y[j], sigma_mu, mu_bar)
  
  p2 <- e1-(e2+e3)
  
  #Calculating Tree prior
  p3 <- log(1-(split(alpha, beta, depth(trnew[[(base)]][1])))) -
          log(1-(split(alpha, beta, depth(tr[[(((findid((tr[[base]][1]*2),tr))))]][1])))) -
            log(1-(split(alpha, beta, depth(tr[[(((findid((tr[[base]][1]*2+1),tr))))]][1])))) -
              log((split(alpha, beta, depth(trnew[[(base)]][1])))) -
                log(p*length(n_j_adj))
                
  
  r <- exp(p1+p2+p3)
  #print(r)
}

#Selecting new tree
u <- runif(1)

if(u<=r){
  tr <- trnew
  y_suff <- y_suff_new
  
  accept_ratio[j] <-1 
}else{
  accept_ratio[j] <-0 
}

bnode <- findfinal(tr) #find the terminal nodes
bterm <- 0 #terminal nodes counter
for(i in 1:length(bnode)){
  bterm = bterm + bnode[[i]][5] #counter of terminal nodes
}

n_terminal[j] <- bterm


#Drawing sigmas
res <- 0
for(i in 1:length(bnode)){
  if(bnode[[i]][5]==1){
    res <- res + rss(y_suff[[i]],bnode[[i]])
  }
}



sigma_y[j+1] <- sqrt((res + nu*lambda)/rchisq(1,df=(n+nu)))


#Saving the tree
treedraws[[j]] <- tr


if(j %% 100==0){
  print(paste0("Post: ",j," out of ", npost))
}

}
proc.time()-begin

newdata <- cbind(x1,x2,x3)

predict_mybcart <- function(cuts, treedraws, newdata){
    
    #Function to test if it is final node
    
    isleaf <- function(id,trees){
      nodes <- length(trees) #number of nodes on tree i
      local <- 0 #id location
      
      for(q in 1:nodes){
        if(trees[[q]][1]==id){ #walking through the tree
          local <- q #local where the id is
        }
      }
      
      if(local==0){
        return("Error: Could not find node") #Wrong argument: Node does not exist
      }
      
      if((trees[[local]][2]==0) & (trees[[local]][3]==0)){ #if it is final node
        return(1)
      }else{
        return(0)
      }
      
    }
    

    #Verify if id exists
    isid <- function(id,trees){
      nodes <- length(trees) #number of nodes on tree
      local <- 0 #local of the id
      
      for(q in 1:nodes){
        if(trees[[q]][1]==id){ #walks through the tree
          local <- q #local where id is
        }
      }
      
      if(local==0){ #If id does not exist
        return(0) 
      }else{
        return(1)
      }
    }
    
    #Verifica a regra
    isleft <- function(id,treedraws,pred,cuts){
      nodes <- length(treedraws) #number of nodes on tree
      local <- 0 #location of the id
      
      for(q in 1:nodes){
        if(treedraws[[q]][1]==id){ #walks through the tree
          local <- q #place where the id is
        }
      }
      
      if(local==0){
        return("Error: Could not find node") #Some error in the arguments
      }
      
      if(pred[(treedraws[[local]][2])]<cuts[[(treedraws[[local]][2])]][(treedraws[[local]][3])]){ #verify rule
        return(1)
      }else{
        return(0)
      }
      
    }
    
    treeprev <- function(id,trees){
      nodes <- length(trees) #n?mero de n?s na ?rvore
      local <- 0 #local que o id est?
      
      for(q in 1:nodes){
        if(trees[[q]][1]==id){ #caminha pela ?rvore
          local <- q #lugar onde est? o id em quest?o
        }
      }
      
      if(local==0){
        return("Error: Could not find node") #caso tenha colocado argumento errado
      }
      
      return(trees[[local]][4])
    }
    
    compareNodes <- function(x,trees,cuts){
      
      id <- 1 #id do in?cio da ?rvore
      
      while(isleaf(id,trees)==0){
        l <- isleft(id,trees,x,cuts) #verifica a regra
        if(l==1){
          id = 2*id #manda pro n? da esquerda
        }else if(l==0){
          id = 2*id +1 #manda pro n? da direita
        }
      }
      
      return(treeprev(id,trees))
    }
    
    Raux <- matrix(0, ncol = length(treedraws), nrow = nrow(newdata))
    for(h in 1:length(treedraws)){
      for(g in 1:nrow(newdata)){
        Raux[g,h] <- compareNodes(newdata[g,],treedraws[[h]],cuts)
      }
      print(paste0("Post: ",h))
    }
    
    return(Raux)
}

preds <- predict_mybcart(cuts,treedraws[seq(burn+thin,npost, by = thin)],newdata)

preds_means <- rowMeans(as.data.frame(preds))

#################################
######### Some Plots... #########
#################################

par(mfrow=c(2,3))

#Acceptance Ratio Plot
plot(cumsum(accept_ratio)/c(1:npost), type = "l",
     ylab = "Acceptance Ratio",
     xlab = "Posterior Draws")

#Histogram of Posterior Terminal Nodes
hist(n_terminal,
     main = NULL,
     breaks = 30,
     freq = F,
     xlab = "Number of Terminal Nodes")

#Plot of \sigma_i
plot(sigma_y, type = "l",
     ylab = expression(sigma),
     xlab = "Number of Posterior Draws")

#Plot of \sigma_i
plot(sigma_y[seq(burn+thin,npost, by = thin)], type = "l",
     ylab = expression(sigma),
     xlab = "Number of Posterior Draws")
acf(sigma_y[seq(2000,5000, by = 10)],
    main = "")

plot(y,preds_means,
     xlim = c(0,20),
     ylim = c(0,20),
     ylab = "Estimated Y",
     xlab = "True Y")
abline(0,1)

par(mfrow=c(1,1))

##################################
######### STILL MISSING: #########
##################################

#Increase the performance of the model

