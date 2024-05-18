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


Network_function_2 <- function(i, Network) { # recip
	net <- 0
    for(j in 1:length(Network[1,])){
        if(i != j){
            net <- net + Network[i,j]*Network[j,i]
        }
    }
	return(net)
}

Network_function_3 <- function(i, Network) { # transitive triplets (transTrip)
	net <- 0
	for (j in 1:length(Network[i, ])) {
		for(h in 1:length(Network[, j])){
			if(i!=j & j!=h & i!=h){
                
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
            if ( TMP_Net[i, j] == 0){ TMP_Net[i, j] <- 1 }
			else if ( TMP_Net[i, j] == 1){ TMP_Net[i, j] <- 0 }
			Networks[[j]] <- TMP_Net    
		}
        else{
            Networks[[j]] <- Network
        }
        
	}
	return(Networks)
}




Distribution_actors <- function(i, Network, beta, Possible_Ministep_Matrix_vector) {
	Util <- numeric(length(Network[i, ]))
	total_util <- 0
  
	for (h in 1:length(Network[i, ])) {
		total_util <- total_util + exp(Utility(Network, beta, i, h, Possible_Ministep_Matrix_vector))
	}


    
	Util_j <- numeric(length(Network[i, ]))

	for (j in 1:length(Network[i, ])) {
			Util_j[j] <- exp(Utility(Network, beta, i, j, Possible_Ministep_Matrix_vector))
			Util[j] <- Util_j[j] / total_util
	}
    
	return(Util)
}

Utility <- function(Network, beta, i, j, Possible_Ministep_Matrix_vector) {
    NetIJ <- Possible_Ministep_Matrix_vector
	util <- (beta[1] * Network_function_1(i, NetIJ[[j]])) + (beta[2] * Network_function_2(i, NetIJ[[j]])) #+ (beta[3] * Network_function_3(i, NetIJ[[j]]))
	return(util)
}

Simulator <- function(Network0, t = 0, LambdaV, total_Lambda, T, beta) {
	Network <- Network0 
	while (t < T) {
		deltaT <- rexp(1, rate=total_Lambda)
        
		# i <- sample(1:length(Network0[1,]), size=1, prob=(LambdaV / total_Lambda))
        i <- sample(1:length(Network0[1,]), size=1)

        Possible_Ministep_Matrix_vector <- Generate_Possible_Ministep_Matrix_vector(i, Network)
		prob_vector <- 0
        
		prob_vector <- Distribution_actors(i, Network, beta, Possible_Ministep_Matrix_vector)
        
		j <- sample(1:length(Network0[i, ]), size=1, prob=prob_vector)

        if (j != i ){
            if ( Network[i, j] == 0){ Network[i, j] <- 1 }
    		else if ( Network[i, j] == 1){ Network[i, j] <- 0 }    
        }
        
		t <- t + deltaT
	}
	return(Network)
}


size <- 50
# Network0 <- Generate_Network_t0(size = size, sparsity = .75)
LambdaV <- rep(.8, size)
total_Lambda <- sum(LambdaV)
# beta <- rexp(n = 2)
beta <- c(1.5, 1.8, 2.5)
T <- 2.5

library(RSiena)

for (count in 1:50){
	Network0 <- Generate_Network_t0(size = size)
	Network1 <- Simulator(Network0, t = 2.5, LambdaV, total_Lambda, 5, beta)
	Network2 <- Simulator(Network1, t = 5, LambdaV, total_Lambda, 7.5, beta)
	###
	Networks <- array(c(Network0, Network1, Network2), dim = c(size, size, 3))
	Networks <- sienaDependent(Networks, nodeSet="Actors",  allowOnly=FALSE) 
	
	mydata <- sienaDataCreate(Networks)
	
	myeff <- getEffects(mydata)
    myeff <- includeEffects(myeff, recip , include = TRUE)
    myeff <- includeEffects(myeff, density , include = TRUE)

	myalgorithm <- sienaAlgorithmCreate(projname = 'Тест симулятора сети')
	
	ans <- siena07( myalgorithm, data = mydata, effects = myeff, clusterType="FORK", useCluster=TRUE, nbrNodes=12, batch = TRUE, silent = TRUE)
	###
	# totalLambda_file <- file("../data/totalLambda.txt", open = "a")
	# write(total_Lambda ,totalLambda_file)
	rate_file <- file("../data/rate.txt", open = "a")
	write(c(ans$rate), rate_file)
	
	# beta_file <- file("../data/beta.txt", open = "a")
	# write(beta, beta_file)
	theta_file <- file("../data/teta.txt", open = "a")
	write(c(ans$theta), theta_file)
	# lambda_mean_V <- file("../data/lambda_meam.txt", open = "a")
	# write(mean(LambdaV)*T, lambda_mean_V)

}

