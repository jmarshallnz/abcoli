source("multigroup.R")

total_procs <- 4

# theta is:
# 1 gamma (recovery rate - same for all groups)
# 2 q (inactivation of bacteria - same for all groups)
# 3 p (transfer rate from local to global environment - same for all groups)
# 4 z_calves (indirect transmission via local environment)
# 5 z_adults (indirect transmission via local environment)
# 6 s_calves (indirect transmission via global environment)
# 7 s_adults (indirect transmission via global environment)
mean_theta <- c(10, 10, 10, 10, 10);
var_theta  <- c(1, 1, 1, 1, 1);

num_theta <- length(mean_theta)

sample_prior <- function()
{
  	return(rgamma(length(mean_theta), scale=mean_theta/var_theta, shape=var_theta));
}

prior <- function(theta)
{
  	return(prod(dgamma(theta, scale=mean_theta/var_theta, shape=var_theta)));
}

sample_kernel <- function(theta, jump_theta=0.1)
{
	# can return negatives, but we don't care as prior() will take care of it
  	return(rnorm(length(theta), mean=theta, sd=mean_theta*jump_theta))
}

kernel <- function(theta, theta_hat, jump_theta=0.1)
{
	p <- 0;
	for (i in 1:length(theta))
		p <- p + dnorm(theta[i], mean=theta_hat[,i], sd=mean_theta[i]*jump_theta, log=T);
	return(exp(p));
}

simulate_data <- function(data, theta, num_sims, threshold)
{
	sum(rho < threshold) / num_sims;
}

#
# Creates output file if it doesn't exist, returns the number of rows in the output file if it does exist
#
create_output_file <- function(out_file)
{
	ignore <- 0;
	if (file.exists(out_file))
	{
		t <- read.table(out_file, header=T);
		ignore <- nrow(t);
	}
	else
	{
		t <- sample_prior();
		for (i in 1:length(t))
		{
			heading <- sprintf("theta%02d", i);
			cat(heading, "\t", sep="", file=out_file, append=T);
		}
		cat("weight\tnum_under\t", sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("a1%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("a2%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("a3%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("a4%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("a5%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("c1%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("c2%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("c3%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:36) cat(sprintf("c4%02d\t", i), sep="", file=out_file, append=T);
		for (i in 1:35) cat(sprintf("c5%02d\t", i), sep="", file=out_file, append=T);
		cat(sprintf("c536\n", i), sep="", file=out_file, append=T);
	}
	return(ignore);
}

#
# Get initial population of particles from our prior.
#  num_particles	number of particles to generate
#  data				our data
#  threshold		the threshold used for the approximate likelihood
#  out				the output filename to use
#  seed				the output seed to use (for parallel processing)
#  num_sims			the number of simulations to compute the approximate likelihood
#
get_prior_particles <- function(num_particles, summary_stats, threshold, out = "particle1_", seed = 1, num_sims = 5)
{
	out_file <- sprintf("%s%0d.txt", out, seed);
	count_file <- sprintf("%s_count_%0d.txt", out, seed);

	ignore <- create_output_file(out_file);

	set.seed(seed*5000 + ignore);

	n <- 1;
	total <- 1;

	while (n + ignore <= num_particles)
	{
		# Step 1: Sample from our prior
		theta2 <- sample_prior();

		# Step 2: Run our simulations to generate our datasets
		out <- multigroup(theta2, data, summary_stats, 3*threshold, num_sims);
		f <- sum(out$rho < threshold) / num_sims;

		if (mean(out$rho) < 3*threshold && f > 0)
		{
			# Step 3: Generate weights
			cat(theta2, f, sum(out$rho < threshold), out$prev, "\n", sep="\t", file=out_file, append=T);
			cat("Up to particle", n + ignore, "accept/reject is", n/total, "mean(rho) is", mean(out$rho), "\n")
			cat(n+ignore, file=count_file);
			n <- n + 1;
		}
		if (count_population_fast(1) >= num_particles)
			break;
		total <- total + 1;
	}
}

#
# Get population of refined particles that match our model
#  num_particles	number of particles to generate
#  data				our data
#  prev_particles	previous particles from a prior population
#  weights			the weights associated with the previous particles
#  threshold		the threshold used for the approximate likelihood
#  transition_sd	the sd of the transition kernel (multiplier)
#  out				the output filename to use
#  seed				the output seed to use (for parallel processing)
#  num_sims			the number of simulations to compute the approximate likelihood
#
get_refined_particles <- function(num_particles, data, prev_particles, weights, summary_stats, threshold, transition_sd, out = "particle", population = 2, seed = 1, num_sims = 5)
{
	num_prev_particles <- nrow(prev_particles);

	out_file <- sprintf("%s%0d.txt", out, seed);
	count_file <- sprintf("%s_count_%0d.txt", out, seed);

	ignore <- create_output_file(out_file);

	set.seed(seed*5000 + ignore);

	n <- 1;
	total <- 1;

	while (TRUE)
	{
		# Step 5.1: Sample from our previous population and weights
		theta1 <- prev_particles[sample(num_prev_particles, 1, prob=weights, replace=TRUE),];

		# Step 5.2: Perturb based on our markov kernel
		theta2 <- sample_kernel(theta1, transition_sd);

		# Step 5.3: Check prior probability
		p <- prior(theta2);
		if (p > 0)
		{
			# Step 5.4: Run our simulations to generate our datasets
			out <- multigroup(theta2, data, summary_stats, 3*threshold, num_sims);
			f <- sum(out$rho < threshold) / num_sims;
			if (mean(out$rho) < 3*threshold && f > 0)
			{
				# Step 5.5: Update particle
				weight <- f * p / sum(weights*kernel(theta2, prev_particles, transition_sd));
				cat(theta2, weight, sum(out$rho < threshold), out$prev, "\n", sep="\t", file=out_file, append=T);
				cat("Up to particle", n + ignore, "accept/reject is", n/total, "mean(rho) is", mean(out$rho), "\n")
				cat(n+ignore, file=count_file);
				n <- n + 1;
			}
		}
		total <- total + 1;
		if (count_population_fast(population) >= num_particles)
			break;
	}
}

#
# Load our previous particles in and normalize.
#  prefix			filename of previous particles
#
load_particles <- function(prefix="particle")
{
	# read in the previous particles
	particles <- NULL
	if (file.exists(sprintf("%s%d.txt", prefix, 1)))
		particles <- rbind(particles, read.table(file=sprintf("%s%d.txt", prefix, 1), header=T))
	for (i in 2:total_procs)
	{
		if (file.exists(sprintf("%s%d.txt", prefix, i)))
		{
			f <- read.table(file=sprintf("%s%d.txt", prefix, i), header=T)
			particles <- rbind(particles, f)
		}
	}
	# normalize the weights
	particles$weight <- particles$weight / sum(particles$weight)
	return(particles)
}

#
# Generate particles for a given population number and threshold.
#  population		which population we're up to (set to the previous filename found)
#  threshold		the threshold to use for the new population
#  num_particles	the number of particles to generate
#  proc				the processor number, for parallel processing
#
generate_particles <- function(population, summary_stats, threshold, num_particles, proc = 1)
{
	transition_sd <- c(0.2, 0.1, 0.05)
	tau <- transition_sd[min(population, length(transition_sd))];
	particles <- load_particles(sprintf("particle%d_", population));
	get_refined_particles(num_particles, data, as.matrix(particles[,1:num_theta]), particles$weight, summary_stats, threshold, tau, sprintf("particle%d_", population+1), population + 1, proc)
}

#
# Sleep until the given population reaches the given number of particles.  Useful for parallel processing
#  population		the population to wait for
#  num_particles	the number of particles to wait for
#
wait_for_population <- function(population, num_particles)
{
	particles <- load_particles(sprintf("particle%d_", population));
	while (nrow(particles) < num_particles)
	{
		# wait for some time
		Sys.sleep(5)
		particles <- load_particles(sprintf("particle%d_", population));	
	}
}

count_population <- function(population)
{
	particles <- load_particles(sprintf("particle%d_", population));
	return(nrow(particles))
}


count_population_fast <- function(population)
{
	# read in the previous particles
	num_particles <- 0
	if (file.exists(sprintf("particle%d__count_%0d.txt", population, 1)))
		num_particles <- num_particles + scan(file=sprintf("particle%d__count_%0d.txt", population, 1), quiet=T)
	for (i in 2:total_procs)
	{
		if (file.exists(sprintf("particle%d__count_%0d.txt", population, i)))
		{
		num_particles <- num_particles + scan(file=sprintf("particle%d__count_%0d.txt", population, i), quiet=T)
		}
	}

	return(num_particles)
}

plot_posterior <- function(population, col="grey60")
{
	particles <- load_particles(sprintf("particle%d_", population));
	pdf(file=sprintf("posterior%d.pdf", population), width=8, height=11);
	par(mfrow=c(3,3))
	for (i in 1:num_theta)
	{
		plot(density(particles[,i]), main="", ylab="", xlab=names(particles)[i], col=col)
	}
	dev.off()
}

plot_correlation <- function(population, col="grey60")
{
	particles <- load_particles(sprintf("particle%d_", population));
	pdf(file=sprintf("correlations%d.pdf", population), width=20, height=20);

        pairs(particles[,1:num_theta])
	dev.off()
}

plot_sims <- function(population, num_sims = 5)
{
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	pdf(file=sprintf("median_fit%d.pdf", population), width=8, height=11);
	param_med <- apply(particles[,1:num_theta], 2, median);
        param_mode <- rep(0,num_theta);
        for (i in 1:num_theta)
 	{
          d <- density(particles[,i])
          param_mode[i] <- d$x[which.max(d$y)]
        }
	plot_multigroup(param_mode, data, num_sims);
	dev.off();
}


plot_test <- function(population, val, num_sims = 5)
{
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	pdf(file=sprintf("median_fit_test%d.pdf", population), width=8, height=11);
	for (i in val)
	{
		plot_multigroup(as.matrix(particles[i,1:num_theta]), data, num_sims);
	}
	dev.off();
}

calc_dist <- function(population, num_simulations=5)
{
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	# repeat our data
	month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
	nyears <- 3
	date_range <- matrix(0, nyears*12, 2);
	date_range[1,] <- c(0, month_days[1]);
	for (i in 2:(nyears*12))
		date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);
	mth_range <- 1:36
	d <- rowMeans(date_range[mth_range,])

        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)

	dist_a <- NULL
	dist_c <- NULL
        peri_a <- NULL
        peri_c <- NULL
	part   <- NULL
	for (j in 1:nrow(particles))
	{

  # run a bunch of simulations
  get_parameters(as.matrix(particles[j,1:num_theta]));
  # load the initial conditions
  get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

  confidence <- 0.95

  p <- matrix(0, num_simulations, nrow(rep_data)*2)
  for (i in 1:num_simulations)
  {
    cat("running simulation", i, "of", num_simulations, "\n");

    get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

    prev <- run_model(params$end_time, temp);

    month_prev <- interpolate_prev(date_range, prev);

    p[i,] <- c(month_prev[,1], month_prev[,2])

    m_c  <- p[i,((0)*nrow(rep_data)+1):(1*nrow(rep_data))];
    m_a  <- p[i,((1)*nrow(rep_data)+1):(2*nrow(rep_data))];

#		prev_c <- matrix(0,nrow(rep_data), 5)
#		prev_a <- matrix(0,nrow(rep_data), 5)
#		for (i in 1:5)
#		{
#			prev_c[,i] <- as.numeric(particles[j,(num_theta+2 + 2*(i-1)*nrow(rep_data) + 1):(num_theta+2 + (2*i-1)*nrow(rep_data))])
#			prev_a[,i] <- as.numeric(particles[j,(num_theta+2 + (2*i-1)*nrow(rep_data) + 1):(num_theta+2 + (2*i)*nrow(rep_data))])
#		}
#		m_c  <- apply(prev_c, 1, mean, na.rm = TRUE);
#		m_a  <- apply(prev_a, 1, mean, na.rm = TRUE);
		dist_c <- c(dist_c,sum((m_c - rep_data[,1])^2*rep_data[,3])/sum(rep_data[,3]))
		dist_a <- c(dist_a,sum((m_a - rep_data[,2])^2*rep_data[,4])/sum(rep_data[,4]))
		p_c <- 0
		p_a <- 0
		for (l in 1:(nyears-1))
		{
			for (k in (l+1):nyears)
			{
				p_c <- p_c + sum((m_c[(k*12-11):(k*12)] - m_c[(l*12-11):(l*12)])^2*rep_data[1:12,3])/sum(rep_data[1:12,3])
				p_a <- p_a + sum((m_a[(k*12-11):(k*12)] - m_a[(l*12-11):(l*12)])^2*rep_data[1:12,4])/sum(rep_data[1:12,4])
			}
		}
		peri_c <- c(peri_c, p_c)
		peri_a <- c(peri_a, p_a)
	}
	part <- rbind(part, p)
}


	dist <- cbind(dist_c, peri_c, dist_a, peri_a, part)
	return(dist)
}


find_summary_stats <- function(population, npowers = 1, trace = 0)
{
	particles <- load_particles(sprintf("particle%d_", population));
	particles <- unique(particles)
	data <- load_data();
	# repeat our data
	month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
	nyears <- 3
	date_range <- matrix(0, nyears*12, 2);
	date_range[1,] <- c(0, month_days[1]);
	for (i in 2:(nyears*12))
		date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);
	mth_range <- 1:36
	d <- rowMeans(date_range[mth_range,])

        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)

	data_set <- c(rep_data[,1],rep_data[,2]);

	num_simulations <- 5;

	# assemble the data
	theta <- matrix(0, nrow(particles)*num_simulations, num_theta);
	sims  <- matrix(0, nrow(particles)*num_simulations, nrow(rep_data)*2);

	for (j in 1:nrow(particles))
	{
		for (i in 1:num_simulations)
		{
			row <- num_simulations*(j-1) + i;
			theta[row,] <- as.numeric(particles[j,1:num_theta]);
			sims[row,]  <- as.numeric(particles[j,(num_theta+2 + (2*i-2)*nrow(rep_data) + 1):(num_theta+2 + (2*i)*nrow(rep_data))])
		}
	}

	summary_stats <- matrix(0, ncol(sims)*npowers, num_theta)
	summaries <- list()
	for (j in 1:num_theta)
	{
		# construct a dataframe for particle linear models
		s <- NULL;
		for (i in 1:npowers)
			s <- cbind(s, sims^i)
		s <- data.frame(s)
		if (j == 2)
		{
			d <- data.frame(param = log(theta[,j]) - log(theta[,j+2]))
		} else if (j == 3) {
			d <- data.frame(param = log(theta[,j]) - log(theta[,j+2]) - log(theta[,1]))
		} else {
			d <- data.frame(param = log(theta[,j]))
		}
		for (i in 1:npowers)
			d <- cbind(d,s[,c(1:3,13:15,25:27,37:72)+(i-1)*72])

		m <- lm(formula(d), data=d)

		# simplify the model
		m <- step(m, direction="both", trace=trace)

		# terms left in the model
		terms <- attr(m$terms,"term.labels")
		index <- rep(0,length(terms))
		for (i in 1:length(terms))
			index[i] <- which(names(s)==terms[i])

		summary_stats[index,j] <- m$coef[-1]

		# do some plotting of the model against the results
#		reg <- colSums(m$coef[-1] * t(s), na.rm=T) + m$coef[1]
#		plot(reg, theta[,j] - reg)
#		lines(c(-1e10,1e10), c(0,0))

		# test fit of data
		data <- NULL;
		for (i in 1:npowers)
			data <- c(data, data_set^i)
		data <- data.frame(t(data))

		names(data) <- names(s)
#		print(data)

		summaries[[j]] <- summary(m)
		print(summary(m))
		print(exp(predict(m, newdata=data, interval="confidence")))

		print(exp(sum(data * summary_stats[,j],na.rm=T)+m$coef[1]))
	}
	for (j in 1:num_theta)
	{
		print(summaries[[j]])
	}
	return(summary_stats)
}


# check our summary stats to see how much data matches

check_summary_stats <- function(population, summary_stats)
{
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	# repeat our data
	month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
	nyears <- 3
	date_range <- matrix(0, nyears*12, 2);
	date_range[1,] <- c(0, month_days[1]);
	for (i in 2:(nyears*12))
		date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);
	mth_range <- 1:36
	d <- rowMeans(date_range[mth_range,])

        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)

	npowers <- nrow(summary_stats) / (nrow(rep_data)*2);

	data_set <- NULL
	for (power in 1:npowers)
		data_set <- c(data_set,rep_data[,1]^power,rep_data[,2]^power);

	# compute the summary stat for the data
	ss_data <- colSums(data_set * summary_stats);

	# grab out our simulation data
	ss_diff <- matrix(0, nrow(particles), num_simulations)
	for (j in 1:nrow(particles))
	{
		for (i in 1:num_simulations)
		{
			theta <- as.numeric(particles[j,1:num_theta]);
			sims  <- as.numeric(particles[j,(num_theta+2 + (2*i-2)*nrow(rep_data) + 1):(num_theta+2 + (2*i)*nrow(rep_data))])
			sim_data <- NULL;
			for (power in 1:npowers)
				sim_data <- c(sim_data, sims^power)
			# compute the summary stat
			diff <- colSums(sim_data * summary_stats) - ss_data;
			ss_diff[j,i] <- sum(diff[c(1:5,6:7)]^2)
		}
	}
	# now order by the ss_diff's and plot the first K
	
	o <- order(rowSums(ss_diff))
	pdf(file=sprintf("ss_fit_test%d.pdf", population), width=8, height=11);
	par(mfrow=c(2,1))
	for (i in 1:min(nrow(particles),200))
	{
		for (an in 1:2)
		{
			max_prev <- 0.5
	  		plot(d, rep_data[mth_range,an], ylim=c(0,max_prev), type="l", xlab="", ylab="Prevelance", main="", xaxt="n", lwd=2)
			for (j in 1:num_simulations)
			{
				sims  <- as.numeric(particles[o[i],(num_theta+2 + (2*j-3+an)*nrow(rep_data) + 1):(num_theta+2 + (2*j-2+an)*nrow(rep_data))])
				lines(d, sims, col=j)
			}
			mtext(o[i], 3, line=0,at=d[18])
			for (j in 1:ncol(ss_diff))
				mtext(ss_diff[o[i],j],4,line=0,at=max_prev*(j)/7, las=1, col=j)
			mtext(sum(ss_diff[o[i],])/5,4,line=0,at=max_prev*(6)/7, las=1)
			axis(1, at=c(date_range[,1], date_range[36,2]), labels=FALSE)
			for (j in 1:36)
			{
				months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
				mtext(months[(j+5) %% 12 + 1], side=1, at=d[j], cex=0.8);
			}
		}
	}
	dev.off()
	print(particles[o, 1:num_theta])
}

update_weights <- function(population, summary_stats, threshold)
{
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	# repeat our data
	month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
	nyears <- 3
	date_range <- matrix(0, nyears*12, 2);
	date_range[1,] <- c(0, month_days[1]);
	for (i in 2:(nyears*12))
		date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);
	mth_range <- 1:36
	d <- rowMeans(date_range[mth_range,])

        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)

	npowers <- nrow(summary_stats) / (nrow(rep_data)*2);

	data_set <- NULL
	for (power in 1:npowers)
		data_set <- c(data_set,rep_data[,1]^power,rep_data[,2]^power);

	# compute the summary stat for the data
	ss_data <- colSums(data_set * summary_stats);

	# grab out our simulation data
 	weights <- rep(0,nrow(particles))
	for (j in 1:nrow(particles))
	{
		ss_diff <- rep(10000, num_simulations)
		for (i in 1:num_simulations)
		{
			theta <- as.numeric(particles[j,1:num_theta]);
			sims  <- as.numeric(particles[j,(num_theta+2 + (2*i-2)*nrow(rep_data) + 1):(num_theta+2 + (2*i)*nrow(rep_data))])
			sim_data <- NULL;
			for (power in 1:npowers)
				sim_data <- c(sim_data, sims^power)
			# compute the summary stat
			diff <- colSums(sim_data * summary_stats) - ss_data;
			ss_diff[i] <- sum(diff[c(1,4:5)]^2)
		}
		weights[j] <- sum(ss_diff < threshold) / num_simulations;

	}
	print(cbind(weights, particles$num_under/5))
	return(weights)
}


# load in our particles from run 1 (assumption: the ABC didn't over-refine them)
# do a linear model of our predictors versus the parameters
find_summaries <- function(param, summary_stats, npowers = 1)
{
  # construct a dataframe for particle linear models
  s <- NULL;
  for (i in 1:npowers)
    s <- cbind(s, summary_stats^i)
  s <- data.frame(s)
  d <- s;
  d$param <- param
  m <- lm(param ~ ., data=d)
  # do some plotting of the model against the results
  reg <- colSums(m$coef[-1] * t(s), na.rm=T) + m$coef[1]
  plot(reg, param - reg)
  lines(c(min(param),max(param)), c(0,0))
  return(m)
}

plot_dist <- function(population, num_simulations=5)
{
	pdf(file=sprintf("median_fit_test%d.pdf", population), width=8, height=11);
	par(mfrow=c(2,1))
	particles <- load_particles(sprintf("particle%d_", population));
	data <- load_data();
	# repeat our data
	month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
	nyears <- 3
	date_range <- matrix(0, nyears*12, 2);
	date_range[1,] <- c(0, month_days[1]);
	for (i in 2:(nyears*12))
		date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);
	mth_range <- 1:36
	d <- rowMeans(date_range[mth_range,])

        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)

	for (j in 1:nrow(particles))
	{

  # run a bunch of simulations
  get_parameters(as.matrix(particles[j,1:num_theta]));
  # load the initial conditions
  get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

  confidence <- 0.95

  p <- matrix(0, num_simulations, nrow(rep_data)*2)
  for (loop in 1:num_simulations)
  {
    cat("running simulation", loop, "of", num_simulations, "\n");

    get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

    prev <- run_model(params$end_time, temp);

    month_prev <- interpolate_prev(date_range, prev);

    p[loop,] <- c(month_prev[,1], month_prev[,2])
  }

		for (an in 1:2)
		{
			prev <- matrix(0, 5, nrow(rep_data))
			dist <- rep(0, 5)
		        peri <- rep(0, 5)
		    	for (i in 1:5)
			{
				prev[i,] <- as.numeric(particles[j,(num_theta+2 + (2*i-3+an)*nrow(rep_data) + 1):(num_theta+2 + (2*i-2+an)*nrow(rep_data))])
				dist[i] <- sum((prev[i,] - rep_data[,an])^2*rep_data[,an+2])/sum(rep_data[,an+2])
				for (l in 1:(nyears-1))
				{
					for (k in (l+1):nyears)
					{
						peri[i] <- peri[i] + sum((prev[i,(k*12-11):(k*12)] - prev[i,(l*12-11):(l*12)])^2*rep_data[1:12,an+2])/sum(rep_data[1:12,an+2])
					}
				}
			}
			max_prev <- max(rep_data[,an], prev)
  			plot(d, rep_data[mth_range,an], ylim=c(0,max_prev), type="l", xlab="", ylab="Prevelance", main="", xaxt="n", lwd=2)
			for (i in 1:5)
			{
				lines(d, prev[i,], col=i+1, lty="dashed")
				mtext(dist[i]*100,4,line=0,at=max_prev*(i+0.2)/7, las=1, col=i+1)
				mtext(peri[i]*100,4,line=0,at=max_prev*(i-0.2)/7, las=1, col=i+1)
    			}
    m  <- apply(p[1:num_simulations,((an-1)*nrow(rep_data)+1):(an*nrow(rep_data))], 2, median, na.rm = TRUE);
    lq <- apply(p[1:num_simulations,((an-1)*nrow(rep_data)+1):(an*nrow(rep_data))], 2, quantile, (1 - confidence)/2, na.rm = TRUE);
    uq <- apply(p[1:num_simulations,((an-1)*nrow(rep_data)+1):(an*nrow(rep_data))], 2, quantile, (1 + confidence)/2, na.rm = TRUE);
    lines(d, m, col="grey60")
    lines(d, lq, col="grey60", lty="dashed")
    lines(d, uq, col="grey60", lty="dashed")
			dist <- sum((m - rep_data[,an])^2*rep_data[,an+2])/sum(rep_data[,an+2])
			peri <- 0;
			for (l in 1:(nyears-1))
			{
				for (k in (l+1):nyears)
				{
					peri <- peri + sum((m[(k*12-11):(k*12)] - m[(l*12-11):(l*12)])^2*rep_data[1:12,an+2])/sum(rep_data[1:12,an+2])
				}
			}
			mtext(dist*100,4,line=0,at=max_prev*(6+0.2)/7, las=1, col="grey60")
			mtext(peri*100,4,line=0,at=max_prev*(6-0.2)/7, las=1, col="grey60")
			mtext(j, 3, line=0,at=d[18])
			axis(1, at=c(date_range[,1], date_range[36,2]), labels=FALSE)
			for (i in 1:36)
			{
				months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
				mtext(months[(i+5) %% 12 + 1], side=1, at=d[i], cex=0.8);
			}
		}
	}
	dev.off()
}

calc_accuracy <- function(population)
{
	particles <- load_particles(sprintf("particle%d_", population));
	particles[order(particles$num_under),1:(num_theta+2)]
#	data <- load_data();
	# repeat our data
	nyears <- 3
        rep_data <- NULL
        for (i in 1:nyears)
          rep_data <- rbind(rep_data, data)
	dist_a <- matrix(0, nrow(particles), 5)
	dist_c <- matrix(0, nrow(particles), 5)
        peri_a <- matrix(0, nrow(particles), 5)
        peri_c <- matrix(0, nrow(particles), 5)
	for (j in 1:nrow(particles))
	{
	    for (i in 1:5)
		{
			prev_c <- particles[j,(num_theta+2 + 2*(i-1)*nrow(rep_data) + 1):(num_theta+2 + (2*i-1)*nrow(rep_data))]
			prev_a <- particles[j,(num_theta+2 + (2*i-1)*nrow(rep_data) + 1):(num_theta+2 + (2*i)*nrow(rep_data))]
			# find difference with rep_data
        	        dist_c[j,i] <-   sum((prev_c - rep_data[,1])^2*rep_data[,3])/sum(rep_data[,3])
        	        dist_a[j,i] <- 4*sum((prev_a - rep_data[,2])^2*rep_data[,4])/sum(rep_data[,4])
			for (k in 2:nyears)
			{
				peri_a[j,i] <- peri_a[j,i] +   3*sum((prev_a[(k*12-11):(k*12)] - prev_a[1:12])^2*rep_data[1:12,4])/sum(rep_data[1:12,4])
				peri_c[j,i] <- peri_c[j,i] + 4*3*sum((prev_a[(k*12-11):(k*12)] - prev_a[1:12])^2*rep_data[1:12,3])/sum(rep_data[1:12,3])
			}
		}
	}
	dist <- sqrt(dist_c + dist_a + peri_c + peri_a)
	o <- order(apply(dist,1,min))
        print(cbind(particles[o,1:num_theta], apply(dist,1,min)[o], apply(dist, 1, mean)[o], apply(dist, 1, max)[o]))
	o <- order(apply(dist,1,mean))
	print(cbind(particles[o,1:num_theta], apply(dist,1,min)[o], apply(dist, 1, mean)[o], apply(dist, 1, max)[o]))
	o <- order(apply(dist,1,max))
	print(cbind(particles[o,1:num_theta], apply(dist,1,min)[o], apply(dist, 1, mean)[o], apply(dist, 1, max)[o]))
}

data <- load_data();
summary_stats <- read.table("summary_stats.txt", header=T)

threshold <- c(150,70,50,30,20)
num_particles <- c(5000,5000,5000,5000)
get_prior_particles(5000, summary_stats, threshold[1], "particle1_", processor);
for (i in 2:length(threshold))
{
	generate_particles(i-1, summary_stats, threshold[i], num_particles[i], processor)
}
