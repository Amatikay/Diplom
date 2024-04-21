Generate_Network_t0 <- function(size = 40, sparsity = 0.7) {
	network <- matrix(0, nrow=size, ncol=size)
  
	for (i in 1:size) {
		for (j in 1:size) {
			if (runif(1) > sparsity) {
				network[i, j] <- 1
			}
		}
	}
  
	return(network)
}


Network_function_1 <- function(i, Network) { # outdegree (density)
	net <- sum(Network[i, ])
	return(net)
}


Network_function_3 <- function(i, Network) { # transitive triads (transTrip)
	net <- 0
	for (j in 1:length(Network[i, ])) {
		for(h in 1:length(Network[i, ])){
			if(i!=j & j != h){
				net <- net + (Network[i,j] * Network[i,h] * Network[h,j]) 
			}
		}
	}
	return(net)
}

Generate_Possible_Ministep_Matrix_vector <- function(i, Network) {
	Networks <- list()
  
	for (j in 1:length(Network[i, ])) {
		if (i != j) {
			TMP_Net <- Network
			TMP_Net[i, j] <- 1
			Networks[[j]] <- TMP_Net    
		}
	}
	return(Networks)
}

Ministep <- function(Network, i, j) {
	Network[i, j] <- 1
	return(Network)
}


Distribution_actors <- function(i, Network, beta) {
	Util <- numeric( length(Network[i, ]) )
	total_util <- 0
  
	for (h in 1:length(Network[i, ])) {
		if (i != h) {
			total_util <- total_util + exp(Utility(Network, beta, i, h))
		}
	}
  
	Util_j <- numeric(length(Network[i, ]))
  
	for (j in 1:length(Network[i, ])) {
		if (j != i) {
			Util_j[j] <- exp(Utility(Network, beta, i, j))
            # print(Util_j[j])
            # print(total_util)
			Util[j] <- Util_j[j] / total_util
		}
	}
	Util[which(is.na(Util))] <- 0
	return(Util)
}

Utility <- function(Network, beta, i, j) {
	util <- 0
	util <- util + (beta[1] * Network_function_1(i, Network)) + (beta[2] * Network_function_3(i, Network))
	return(util)
}

Simulator <- function(Network0, t = 0, LambdaV, total_Lambda, T, beta) {
	Network <- Network0 
	while (t < T) {
		deltaT <- rexp(1, rate=total_Lambda)
		i <- sample(1:length(Network0[1, ]), size=1, prob=(LambdaV / total_Lambda))
      
		prob_vector <- 0
		prob_vector <- Distribution_actors(i, Network, beta)
        # print(sum(prob_vector))
		j <- sample(1:length(Network0[i,]), size=1, prob=prob_vector)
		Network[i, j] <- 1
		t <- t + deltaT

	}
	return(Network)
}

size <- 30
Network0 <- Generate_Network_t0(size = size)
LambdaV <- runif(n = size, min = 0.1, max = 1) 
total_Lambda <- sum(LambdaV)
total_Lambda
beta <- c(runif(n = 1, min = -2, max = 2), runif(n = 1, min = -2, max = 2)) 
# beta <- c(1.3,1.6)
T <- 3

library(RSiena)

for (count in 1:2){
	Network0 <- Generate_Network_t0(size = size)
	Network1 <- Simulator(Network0, t = 0, LambdaV, total_Lambda, T, beta)
	Network2 <- Simulator(Network1, t = 0, LambdaV, total_Lambda, T, beta)
	###
	Networks <- array(c(Network0, Network1, Network2), dim = c(size, size, 3))
	Networks <- sienaDependent(Networks, nodeSet="Actors")
	
	mydata <- sienaDataCreate(Networks)
	
	myeff <- getEffects(mydata)
	
	myalgorithm <- sienaAlgorithmCreate(projname = 'Тест симулятора сети')
	
	ans <- siena07( myalgorithm, data = mydata, effects = myeff, clusterType="FORK", useCluster=TRUE, nbrNodes=12)
	###
	totalLambda_file <- file("./data/totalLambda.txt", open = "a")
	write(total_Lambda ,totalLambda_file)
	rate_file <- file("./data/rate.txt", open = "a")
	write(c(ans$rate), rate_file)
	
	beta_file <- file("./data/beta.txt", open = "a")
	write(beta, beta_file)
	theta_file <- file("./data/teta.txt", open = "a")
	write(c(ans$theta), theta_file)
}

