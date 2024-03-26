library(igraph)
library(compiler)




####################
# KURAMOTO
####################
KURAMOTO_r <- function(network, size, dt=5e-2, sd.meas.noise=0.1, sd.dyn.noise=0., sigma=1){
    if(is.null(igraph::E(network)$weight)) stop("No edge weights specified.")
    
    #prepare time series and parameters vectors
    S <- c()    #this is the network state, consisting of 3*N variables
    M <- 1      #number of variables for the model
    idxs <- list()

    W <- igraph::get.adjacency(network, attr="weight", sparse=F)
    Nodes <- length(igraph::V(network))

    omegas <- rnorm(Nodes, 0, 1/Nodes)

    #print("INIT")

    for(m in 1:Nodes){
        #state <- c(X = 0)
        #n-state equivalent:
        S <- c(S, runif(1, 0, 2*pi))
    }

    for(m in 1:M){
        idxs[[m]] <- (1:Nodes-1)*M + m
    }

    #   MODEL
    #KURAMOTO:  d theta_i/dt = w_i + s*sum_j W_ij * sin(theta_j - theta_i) =
    #           w_i + cos(theta_i) * sum_j W_ij * sin(theta_j) - sin(theta_i) * sum_j W_ij * cos(theta_j)
                # also the variant with 1/N or 1/k_i as prefactor for the interaction
                #see https://arxiv.org/pdf/1511.07139.pdf
    model <- function(t, S, parameters) {
        dS <- c()

        ###################################
        # self term
        ###################################

        dS[idxs[[1]]] <- omegas[idxs[[1]]]

        ###################################
        # interaction term
        ###################################

        #strength <- rowSums(W)
        coupling <- sigma/degree(network)
        #coupling <- sigma/Nodes

        dS[idxs[[1]]] <- dS[idxs[[1]]] + coupling * ( cos(S[ idxs[[1]] ]) * (W %*% sin(S[ idxs[[1]] ])) - sin(S[ idxs[[1]] ]) * (W %*% cos(S[ idxs[[1]] ])) )

        ###################################
        # dynamic noise term
        ###################################

        if(sd.dyn.noise>0){
           dS <- dS + rnorm(length(dS), 0, sd.dyn.noise)
        }
        list(dS)
    }

    #generate time series
    times <- seq(0, size*dt, by = dt)
    library(deSolve)

    #print("DSOLVE")

    multi <- ode(y = S, times = times, func = model, parms = parameters)
    #print(paste("DEBUG", dim(multi)))

    #transform in the MultiTS format
    #get only the x variable of the 3D state
    x <- list()
    for(m in 1:Nodes){
        #skip initial condition and column 1 which is time
        x[[m]] <- multi[2:nrow(multi), 1 + (m-1)*M + 1]
    }

    #add measurement noise
    if(sd.meas.noise>0){
        for(m in 1:Nodes){
            x[[m]] <- x[[m]] + rnorm(size,0,sd.meas.noise)
        }
    }

    return(x)
}
KURAMOTO <- cmpfun(KURAMOTO_r)



####################
# OTHER FUNCTIONS
####################


normalize_vec <- function(v_centr){
     return((v_centr - min(v_centr))/( max(v_centr) - min(v_centr) ))
}

vec2pal <- function(v_centr, mypal){
     val_idxs <- 1 + floor((length(mypal)-1)*normalize_vec(v_centr))
     return(mypal[val_idxs])
}

plot.MultiTS <- function(MultiTS){
    library(ggplot2)

    dat <- data.frame()

    tseq <- 1:length(MultiTS[[1]])

    for(m in 1:length(MultiTS)){
        dat <- rbind(dat, data.frame(node=m, time=tseq, value=MultiTS[[m]]))
    }

    return( 
        ggplot(dat, aes(time, sin(value), group=node, color=node)) + theme_bw() + theme(panel.grid=element_blank()) +
        geom_line(alpha=0.3) + 
        scale_color_viridis_c() +
        xlab("Time") + ylab("sin(theta)")
        )
}

MultiTS2Matrix <- function(MultiTS){
    N <- length(MultiTS)
    M <- length(MultiTS[[1]])

    A <- matrix(0, nrow=N, ncol=M)

    for(n in 1:N){
        A[n,] <- MultiTS[[n]]
    }

    return(A)
}


