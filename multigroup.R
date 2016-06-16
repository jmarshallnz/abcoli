S <- NULL
I <- NULL
E <- NULL
G <- NULL
params <- NULL
fast_params <- NULL

get_next_event <- function(S, I, E, G, params)
{
  returned_data <- .C('get_exp_event', num_groups = as.integer(length(S[1,])),
                                   num_sub_groups = as.integer(length(S[,1])),
                                                S = as.double(S),
                                                I = as.double(I),
                                                E = as.double(E),
                                                G = as.double(G),
                                           params = as.double(fast_params),
                                            event = double(5))

  return(returned_data$event)
}

rdiscrete <- function(values)
{
  u <- runif(1)*sum(values)

#  cat("matching ", u, " from cumm:\n")
  cumm <- 0;
  n <-length(values);
  for(j in 1:n)
  {
    cumm <- cumm + values[j];
#    cat(cumm, " ")
    if (u < cumm)
      return(j);
  }
  return(n)
}

get_next_event_R <- function(S, I, E, G, params)
{
  num_groups = length(S[1,]);

  s <- colSums(S);
  i <- colSums(I);

  rates <- matrix(0, num_groups, 11)
  rates[, 1] <- params$death      * s;     # death of susceptible
  rates[, 2] <- params$death      * i;     # death of infected
  rates[, 3] <- params$culling    * s;     # culling of susceptible
  rates[, 4] <- params$culling    * i;     # culling of infected
  rates[, 5] <- params$inf_direct * s * i; # direct infection
  rates[, 6] <- params$inf_local  * s * E; # indirect infection via consumption from local env
  rates[, 7] <- params$inf_local  * i * E; # consumption from local env by infected
  rates[, 8] <- params$inf_global * s * G; # indirect infection via consumption from global env
  rates[, 9] <- params$inf_global * i * G; # consumption from global env by infected
  rates[,10] <- params$recovery   * i;     # recovery by infected
  rates[,11] <- params$fakebirth  * (s+i); # lactating events (non-seasonal)

  total_rate = sum(rates);

#  cat("rates:\n")
#  for (i in 1:num_groups)
#    cat(rates[i,],"\n")
#  cat("\n")

  # generate a random exponential for time to next event, followed by a uniform for event type

  event <- rexp(1, total_rate);
  event_num <- rdiscrete(t(rates))
  event[2] <- ((event_num-1) %% 11) + 1;
  event[3] <- ((event_num-1) %/% 11) + 1;
#  cat("\nevent_num=",event_num,"at ", event[2],",", event[3],"\n")
  event[4] <- rdiscrete(S[,event[3]])
  event[5] <- rdiscrete(I[,event[3]])
  return(event)
}

add_to_sorted_queue <- function(queue, event)
{
  if (is.null(queue))
    return(matrix(event, ncol=5))
  queue <- matrix(queue, ncol=5)
  n <- length(queue[,1])
  if (n > 1)
  {
    if (event[1] < queue[1,1])
    {
      return(rbind(event, queue))
    }
    for (i in 2:n)
    {
      if (event[1] < queue[i,1])
      {
        return(rbind(queue[1:(i-1),], event, queue[i:n,]))
      }
    }
  }
  return(rbind(queue, event))
}

queue_exponential_event <- function(event_type, time, event_queue)
{
  # exponential events require rates to be computed, and time is given from now
  new_event <- get_next_event(S, I, E, G, params);
  new_event[1] <- time + new_event[1];

  # remove any old exponential events from the event queue
  if (!is.null(event_queue))
  {
    event_queue <- matrix(event_queue, ncol=5);
    event_queue <- event_queue[event_queue[,2] >= 12]
  }
  return(add_to_sorted_queue(event_queue, new_event))
}

get_annual_event <- function(time)
{
  new_event <- time;
  new_event[2] <- 16;
  new_event[3] <- 0;
  new_event[4] <- 0;
  new_event[5] <- 0;
  return(new_event);
}

get_birth_time <- function(time)
{
  start_of_this_year = floor(time/365)*365;
  if (params$seasonal)
  {
    return(start_of_this_year + rnorm(1, params$birth[1], params$birth[2]));
  }
  else
  {
    return(start_of_this_year + runif(1, 0, 365));
  }
}

get_birth_event <- function(time, group)
{
  new_event    <- time;
#  cat("birth event", new_event, "\n");

  new_event[2] <- 14;    # birth
  new_event[3] <- 0;     # group
  new_event[4] <- group; # subgroup (susceptibles)
  new_event[5] <- 0;     # subgroup (infecteds)

  return(new_event);
}

get_heifer_event <- function(time)
{
  new_event    <- time + rnorm(1, params$maturation_heifer[1], params$maturation_heifer[2]);
#  cat("heifer event", new_event, "\n");
  new_event[2] <- 15;    # heifer merge
  new_event[3] <- 0;     # group
  new_event[4] <- 0;     # subgroup (susceptibles)
  new_event[5] <- 0;     # subgroup (infecteds)
  return(new_event);
}

get_lactating_event <- function(time)
{
  new_event    <- time;
#  cat("queued new lactating event at time", new_event, "\n");
  new_event[2] <- 12;    # lactation starts
  new_event[3] <- 0;     # group
  new_event[4] <- 0;     # subgroup (susceptibles)
  new_event[5] <- 0;     # subgroup (infecteds)
  return(new_event);
}

get_dryoff_event <- function(time)
{
  new_event    <- time + rnorm(1, params$lactation[1], params$lactation[2]);
#  cat("queued new dryoff event at time", new_event, "\n");

  new_event[2] <- 13;    # lactation ends
  new_event[3] <- 0;     # group
  new_event[4] <- 0;     # subgroup (susceptibles)
  new_event[5] <- 0;     # subgroup (infecteds)
  return(new_event);
}

add_replacement <- function(time, cow, event_queue)
{
  birth <- get_birth_event(time, (cow-1) %% length(S[,1]) + 1)
  event_queue <- add_to_sorted_queue(event_queue, birth);

  # and a heifer event
  event_queue <- add_to_sorted_queue(event_queue, get_heifer_event(birth[1]));

  # and queue the lactation event and drying off event based on this
  event_queue <- add_to_sorted_queue(event_queue, get_lactating_event(birth[1]));
  event_queue <- add_to_sorted_queue(event_queue, get_dryoff_event(birth[1]));
  return(event_queue);
}

alter_event_queue <- function(event, time, event_queue)
{
  # alters the queue as a result of an event - doesn't necessarily queue up
  # the specific event that has occured

  if (event == 16)
  {
    # yearly event has occured - we now queue up the birth events (which will occur really soon)
    # and the end of lactation events based on these

    adults <- S[,3:4] + I[,3:4];

    # pick num_ghosts random numbers from num_adults - allows possible hereditary super shedding
    ghost_group <- length(S[1,]);  # ghosts are always the last group
    num_ghosts <- sum(S[,ghost_group] + I[,ghost_group]);

 #   cat("num_ghosts = ", num_ghosts, "num_lact = ", sum(adults[,1]), "num_dry = ", sum(adults[,2]),"\n");

    for (i in 1:num_ghosts)
    {
      birth_time <- get_birth_time(time)
      cow <- rdiscrete(adults);
      # replacements - add birth, heifer, lactation + drying off events for each replacement
      if (params$seasonal)
      {
        event_queue <- add_replacement(birth_time, cow, event_queue)
      }
      adults[cow] <- adults[cow] - 1;
    }
    num_adults <- sum(adults);
    for (i in 1:num_adults) # assume the rest have calves and queue the lactation event
    {
      birth_time <- get_birth_time(time);
      event_queue <- add_to_sorted_queue(event_queue, get_lactating_event(birth_time));
      event_queue <- add_to_sorted_queue(event_queue, get_dryoff_event(birth_time));
    }

    # destroy our ghost cows
    S[,ghost_group] <<- I[,ghost_group] <<- 0;

    # add the next seasonal event
    event_queue <- add_to_sorted_queue(event_queue, get_annual_event(time+365));
  }
  else if (event <= 4 && params$seasonal == 0) # death or culling - replace
  { # add an immediate birth event (plus lactating/dry events that follow it)
    cow <- rdiscrete(S[,3:4] + I[,3:4]);
#cat("death occurred -> queuing birth at time", time, "\n")
    event_queue <- add_replacement(time, cow, event_queue);
  }

  # resample the exponential event
  return(queue_exponential_event(0, time, event_queue));
}

perform_event <- function(time, event, group, S_group, I_group)
{
  ghosts <- length(S[1,]); # ghost group is always the last

  if      (event == 1 || event == 3) # death or culling of susceptible
  {
    S[S_group,group]  <<- S[S_group,group] - 1;
    S[S_group,ghosts] <<- S[S_group,ghosts] + 1;
  }
  else if (event == 2 || event == 4) # death or culling of infected
  {
    I[I_group,group]  <<- I[I_group,group] - 1;
    I[I_group,ghosts] <<- I[I_group,ghosts] + 1;
  }
  else if (event == 5) # direct infection
  {
    I[S_group,group] <<- I[S_group,group] + 1;
    S[S_group,group] <<- S[S_group,group] - 1;
    # TODO: sample shedding profile
  }
  else if (event == 6) # indirect infection, local environment
  {
    I[S_group,group] <<- I[S_group,group] + 1;
    S[S_group,group] <<- S[S_group,group] - 1;
    E[group] <<- max(E[group] - 1, 0);
    # TODO: sample shedding profile
  }
  else if (event == 7) # consumption by infected, local environment
  {
    E[group] <<- max(E[group] - 1, 0);
    # TODO: extend recovery time?
  }
  else if (event == 8) # indirect infection, global environment
  {
    I[S_group,group] <<- I[S_group,group] + 1;
    S[S_group,group] <<- S[S_group,group] - 1;
    G <<- max(G - 1, 0);
    # TODO: sample shedding profile
  }
  else if (event == 9) # consumption by infected, global environment
  {
    G <<- max(G - 1, 0);
  }
  else if (event == 10) # recovery
  {
    S[I_group,group] <<- S[I_group,group] + 1;
    I[I_group,group] <<- I[I_group,group] - 1;
    # TODO: sample shedding profile
  }
  else if (event == 11 || event == 12) # lactating
  {
    # move one cow to lactating...
    adults <- cbind(S[,4],I[,4]);
    if (sum(adults))
    {
#      if (event == 11)
#      {
#        adults_L <- cbind(S[,3],I[,3]);
#        cat("lactating with", sum(adults), "adults at time", time, "(in other group", sum(adults_L), "adults)\n");
#      }
      group <- rdiscrete(adults);
      if (group < 3)
      {
        S[group,4] <<- S[group,4] - 1;
        S[group,3] <<- S[group,3] + 1;
      }
      else
      {
        I[group-2,4] <<- I[group-2,4] - 1;
        I[group-2,3] <<- I[group-2,3] + 1;
      }
    }
  }
  else if (event == 13) # dry off
  {
    # move one cow to dry...
    adults <- cbind(S[,3],I[,3]);
#    cat("dryoff", adults, "\n");
    if (sum(adults))
    {
#      adults_L <- cbind(S[,4],I[,4]);
#      cat("dryoff with", sum(adults), "adults at time", time, "(in other group", sum(adults_L), "adults)\n");
      group <- rdiscrete(adults);
      if (group < 3)
      {
        S[group,4] <<- S[group,4] + 1;
        S[group,3] <<- S[group,3] - 1;
      }
      else
      {
        I[group-2,4] <<- I[group-2,4] + 1;
        I[group-2,3] <<- I[group-2,3] - 1;
      }
    }
  }
  else if (event == 14) # birth
  {
    u1 = runif(1);

    prob_infected = sum(I[,3:4]) / sum(I[,3:4]+S[,3:4]) * params$vertical;
    if (u1 < prob_infected)
    {
      I[S_group,1] <<- I[S_group,1] + 1;
    }
    else
    {
      S[S_group,1] <<- S[S_group,1] + 1;
    }
  }
  else if (event == 15) # heifer rejoin
  {
    # move a calf to being a heifer, and a heifer to the main herd...
    heifers <- cbind(S[,2],I[,2]);
#    cat("heifer to lact event", heifers, "\n");
    if (sum(heifers))
    {
      group <- rdiscrete(heifers);
      if (group < 3)
      {
        S[group,2] <<- S[group,2] - 1;
        S[group,3] <<- S[group,3] + 1;
      }
      else
      {
        I[group-2,2] <<- I[group-2,2] - 1;
        I[group-2,3] <<- I[group-2,3] + 1;
      }
    }
    calves <- cbind(S[,1],I[,1]);
#    cat("calves to heifer event", calves, "\n");
    if (sum(calves))
    {
      group <- rdiscrete(calves);
      if (group < 3)
      {
        S[group,1] <<- S[group,1] - 1;
        S[group,2] <<- S[group,2] + 1;
      }
      else
      {
        I[group-2,1] <<- I[group-2,1] - 1;
        I[group-2,2] <<- I[group-2,2] + 1;
      }
    }
  }
}

simulate_shedding <- function(t)
{
  # TODO: Shedding profiles
  eta <- params$shedding;
  q   <- params$env_local_death;
  p   <- params$env_transfer;
  qg  <- params$env_global_death;

#  cat("eta =", eta, "q =", q, "p =", p, "qg =", qg, "\n")

  # single farm we have the ODEs
  # E'(t) = eta * I - (q + p)*E
  # G'(t) = sum(p*E) - qG * G

  epq <- exp(-(p+q)*t);
  b   <- colSums(eta * I);
  bpq <- b/(q + p);
  eg  <- exp(-qg*t);
#  cat("epq=", epq, ",   b=", b, ",   bpq=", bpq, ",   eg=", eg, "\n");

  G  <<- sum(p*(bpq*(1 - eg)/qg + (E - bpq)*(epq - eg)/(qg - (q+p)))) + G*eg;
  E  <<- bpq + epq*(E - bpq);
}


run_animal_model <- function(end_time)
{
  current_time <- 0;
  if (params$seasonal)
  {
    current_time <- params$birth[1]-4*params$birth[2]; # 4 sd before mean
  }
  if (current_time < 0)
  {
    cat("Warning - possible birth times before January 1, model will be unpredictable\n");
  }
  current_time <- current_time - 365;

  # add an initial drying event to kick off the seasonal component
  event_queue <- NULL;

  num_groups <- length(S[1,]);
  {
    # take a copy of the ghost group...
    ghost_group <- num_groups;  # ghosts are always the last group
    ghost_s <- S[,ghost_group];
    ghost_i  <- I[,ghost_group];
    event_queue <- alter_event_queue(16, current_time, event_queue);
    S[,ghost_group] <<- ghost_s;
    I[,ghost_group] <<- ghost_i;

    # ignore any events from the previous year
    event_queue <- event_queue[event_queue[,1] >= 0]
  }

  current_time <- 0;
  ignore_time <- 365*2 + 31 + 28 + 31 + 30 + 31 + 30; # ignore 2 1/2 years of data for settling

  # and queue up the first exponential event
  event_queue <- queue_exponential_event(0, current_time, event_queue);

#  num_events <- 0;
  event_totals <- rep(0,16*num_groups);
  last_dump_time <- ignore_time;
  
  chunk_size <- 1000;
  p_size <- 0;
  p <- matrix(0, chunk_size, 5);
  while (current_time < ignore_time+end_time)
  {
    event <- event_queue[1,];

    # perform deterministic events (environment)
    simulate_shedding(event[1] - current_time);

    # perform our event
    perform_event(event[1], event[2], event[3], event[4], event[5]);

    # remove it from our queue
    event_queue <- event_queue[-1,];

    # update our time
    current_time <- event[1];

#    num_events <- num_events + 1;
#    if (num_events > 999)
#    {
#      cat("e =", event, ",S =", S, ",I =", I, ",E =", E, ",G =", G, "\n");
#      cat("eventqueue = ", event_queue, "\n")
#      num_events <- 0;
#    }

    # increment the event total
    code <- (event[2]-1)*num_groups+event[3];
    if (event[3] == 0)
      code <- code + 1;
    event_totals[code] <- event_totals[code] + 1;

    # and queue up any consequences of this
    event_queue <- alter_event_queue(event[2], current_time, event_queue);

    # output the current state of the model
    if (current_time - last_dump_time > 1.5)
    {
      if (p_size >= nrow(p))
      	p <- rbind(p, matrix(0, chunk_size, 5))
      p_size <- p_size + 1;
      totals <- colSums(I[,1:4] + S[,1:4])
      p[p_size,] <- c(event[1] - ignore_time, colSums(I[,1:4] + S[,1:4]))
#      cat(file=working, event, colSums(S), colSums(I), E, G, "\n", append=TRUE);
      last_dump_time = current_time;
    }
#    cat("e =", event, ",S =", S, ",I =", I, ",E =", E, ",G =", G, "\n");
  }
  p <- p[1:p_size,]
  # output our totals
  return(p);
}

run_model <- function(end_time, working)
{
  # initialise our data storage
#  if (file.exists(working))
#  {
#    file.remove(working)
#  }

  current_time <- 0;
  if (params$seasonal)
  {
    current_time <- params$birth[1]-4*params$birth[2]; # 4 sd before mean
  }
  if (current_time < 0)
  {
    cat("Warning - possible birth times before January 1, model will be unpredictable\n");
  }
  current_time <- current_time - 365;

  # add an initial drying event to kick off the seasonal component
  event_queue <- NULL;

  num_groups <- length(S[1,]);
#  if (params$seasonal)
  {
    # take a copy of the ghost group...
    ghost_group <- num_groups;  # ghosts are always the last group
    ghost_s <- S[,ghost_group];
    ghost_i  <- I[,ghost_group];
    event_queue <- alter_event_queue(16, current_time, event_queue);
    S[,ghost_group] <<- ghost_s;
    I[,ghost_group] <<- ghost_i;

    # ignore any events from the previous year
    event_queue <- event_queue[event_queue[,1] >= 0]
  }
#  else
#  { # we have to simulate a bunch of dry offs as our lactating group is large
#    # and similarly do a number of births equal to number of heifers (aka expected deaths)
#    num_heifers <- sum(S[,2]+I[,2]);
#    for (i in 1:num_heifers)
#    { # add a death event so that a new birth occurs (this also queues up lactating and dry off)
#      time <- runif(1,-365,0)
#      event_queue <- alter_event_queue(1, time, event_queue)
#    }
#    num_lactating <- sum(S[,3]+I[,3]);
#    for (i in 1:(0.5*num_lactating))
#    { # add a dry off event
#      time <- runif(1,-365,0)
#      event_queue <- alter_event_queue(11, time, event_queue)
#    }
#    # ignore any events from the previous year
#    event_queue <- event_queue[event_queue[,1] >= 0]
#  }

  current_time <- 0;
  ignore_time <- 365*2 + 31 + 28 + 31 + 30 + 31 + 30; # ignore 2 1/2 years of data for settling

  # and queue up the first exponential event
  event_queue <- queue_exponential_event(0, current_time, event_queue);

#  num_events <- 0;
  event_totals <- rep(0,16*num_groups);
  last_dump_time <- ignore_time;
  
  chunk_size <- 1000;
  p_size <- 0;
  p <- matrix(0, chunk_size, 3);
  while (current_time < ignore_time+end_time)
  {
    event <- event_queue[1,];

    # perform deterministic events (environment)
    simulate_shedding(event[1] - current_time);

    # perform our event
    perform_event(event[1], event[2], event[3], event[4], event[5]);

    # remove it from our queue
    event_queue <- event_queue[-1,];

    # update our time
    current_time <- event[1];

#    num_events <- num_events + 1;
#    if (num_events > 999)
#    {
#      cat("e =", event, ",S =", S, ",I =", I, ",E =", E, ",G =", G, "\n");
#      cat("eventqueue = ", event_queue, "\n")
#      num_events <- 0;
#    }

    # increment the event total
    code <- (event[2]-1)*num_groups+event[3];
    if (event[3] == 0)
      code <- code + 1;
    event_totals[code] <- event_totals[code] + 1;

    # and queue up any consequences of this
    event_queue <- alter_event_queue(event[2], current_time, event_queue);

    # output the current state of the model
    if (current_time - last_dump_time > 1.5)
    {
      if (p_size >= nrow(p))
      	p <- rbind(p, matrix(0, chunk_size, 3))
      p_size <- p_size + 1;
      n_calves <- sum(I[,1] + S[,1]);
      n_adults <- sum(I[,2:4] + S[,2:4]);
      p_calves <- 0;
      p_adults <- 0;
      if (n_calves > 0)
      	p_calves <- sum(I[,1]) / n_calves;
      if (n_adults > 0)
      	p_adults <- sum(I[,2:4]) / n_adults;
      p[p_size,] <- c(event[1] - ignore_time, p_calves, p_adults);
#      cat(file=working, event, colSums(S), colSums(I), E, G, "\n", append=TRUE);
      last_dump_time = current_time;
    }
#    cat("e =", event, ",S =", S, ",I =", I, ",E =", E, ",G =", G, "\n");
  }
  p <- p[1:p_size,]
  # output our totals
#  deaths <- event_totals[(0*num_groups+1):(1*num_groups)] +
#            event_totals[(1*num_groups+1):(2*num_groups)]
#  culls  <- event_totals[(2*num_groups+1):(3*num_groups)] +
#            event_totals[(3*num_groups+1):(4*num_groups)]
#  direct   <- event_totals[(4*num_groups+1):(5*num_groups)];
#  indirect <- event_totals[(5*num_groups+1):(6*num_groups)];
#  cross    <- event_totals[(7*num_groups+1):(8*num_groups)];
#  recovery <- event_totals[(9*num_groups+1):(10*num_groups)];
#  fakebirths <- event_totals[(10*num_groups+1):(11*num_groups)];

#  cat("deaths:", deaths, "\n");
#  cat("cullings:", culls, "\n");
#  cat("direct transmission:", direct, "\n");
#  cat("indirect transmission:", indirect, "\n");
#  cat("cross contamination:", cross, "\n");
#  cat("recovery:", recovery, "\n");
#  cat("fake births:", fakebirths, "\n");
  return(p);
}

convert_to_number <- function(d)
{
  value <- rep(0, length(d))
  for (i in 1:length(d))
  {
    value[i] <- as.real(as.character(d[1,i]));
  }
  return(value);
}

get_inits <- function(init_a, init_p, hsp, hsf)
{
  S[1,] <<- init_a;
  I[1,] <<- round(init_a * init_p);
  S[1,] <<- S[1,] - I[1,];

  # expand eta for differential shedding
  eta   <- t(matrix(params$avg_shedding, ncol(S), nrow(S)));

  # if we have X% high shedders, then the total average shedding pressure is
  scale = (hsp*hsf + 1 - hsp);
  eta <- eta / scale;
  eta[2,] <- eta[2,] * hsf;

  S[2,] <<- round(S[1,] * hsp);
  S[1,] <<- S[1,] - S[2,];
  I[2,] <<- round(I[1,] * hsp);
  I[1,] <<- I[1,] - I[2,];

  E <<- rep(0, 5);
  G <<- 0;

  params$shedding <<- eta;

#  cat("S[1,]=", S[1,], "\n");
#  cat("S[2,]=", S[2,], "\n");
#  cat("I[1,]=", I[1,], "\n");
#  cat("I[2,]=", I[2,], "\n");
#  cat("params$shedding[1,]=", params$shedding[1,], "\n");
#  cat("params$shedding[2,]=", params$shedding[2,], "\n");
}

get_parameters <- function(theta = NULL, full_path="none")
{
  S <<- matrix(0, 2, 5);
  I <<- matrix(0, 2, 5);

  if (is.null(theta))
  {
	  # theta is:
	  # 1 beta_calves
	  # 2 beta_adults (beta_heifers = sqrt(beta_adults*beta_calves))
	  # 3 gamma (assumed the same for all animals)
	  # 4 q (assumed the same for all groups)
	  # 5 p (assumed the same for all groups)
	  # 6 z_calves
	  # 7 z_adults
	  # 8 s_calves
	  # 9 s_adults
  	theta <- c(10, 1, 6.833, 1.5, 5, 8, 4, 8, 4)
  }

  # all rates are per day unless otherwise stated
  rho   <- 0.46;                                              # vertical transmission probability

  m     <- c(0, 0, 0, 0.0004);                                # per animal cull rate
  b     <- c(0.000137, 0.000023, 0.0004, 0.0004);             # per animal death rate
#  beta  <- c(theta[1], theta[2], theta[2], theta[2])*0.0001;  # direct transmission (per susceptible per infected)
  beta <- c(0,0,0,0)
  gamma <- rep(6.833, 4)*0.01;                             # per animal recovery rate

  k     <- c(1e5, 1e4, 1e3, 1e3);                             # concentration in faeces     (cfu/g)
  f     <- c(6, 20, 26, 26);                                  # per animal defaecation rate (kg/day)
  iu    <- 100;                                               # size of an infectious unit

  q     <- rep(1.5, 4)*0.01;                             # per iu death rate in local environments
  qG    <- 1.5*0.01;                                     # per iu death rate in global environment
  p     <- rep(theta[1], 4)*0.01;                             # movement of bugs from local to global (per iu)

  z     <- c(theta[2], rep(theta[3], 3)) * 1e-12;             # consumption rate local env (per animal per iu)
  s     <- c(theta[4], rep(theta[5], 3)) * 1e-12;             # consumption rate global env (per animal per iu)

  hsp   <- 0;
  hsf   <- 1;

  seasonal <- 1;
  fb    <- c(0, 0, 0, 0);

  # lactation starts september finishes may
  # calving: Mid July -> End Autust
  # may 15 (end of lact), Aug 20 (start of lact), August1 (birth start), october 20 (heifers join)
  birth  <- c(31+28+31+30+31+30+5, 10);
  lact   <- c(300,3);
  heif   <- c(61,3);
  adult  <- c(365+61,3);

  #        1year + J +  A +  S +  O +  N +  D +  J +  F +  M +  A +  M
              # + 31 + 31 + 30 + 31 + 30 + 31 + 31 + 28 + 31 + 30 + 31;
  end_t  <- 365*3

  init_a <- c(0, 55, 0, 305, 30);
  init_p <- c(0.4,0.2,0.1,0.1,0.2);

  # read in the file...
  if (full_path != "none")
  {
    dat <- read.xls(full_path, colNames = FALSE);
    init_a <- convert_to_number(dat[ 3,2:6]);
    init_p <- convert_to_number(dat[ 4,2:6]);

    m      <- convert_to_number(dat[11,2:5]);
    b      <- convert_to_number(dat[12,2:5]);
    beta   <- convert_to_number(dat[14,2:5]);
    gamma  <- convert_to_number(dat[16,2:5]);
    k      <- convert_to_number(dat[18,2:5]);
    f      <- convert_to_number(dat[19,2:5]);
    q      <- convert_to_number(dat[21,2:5]);
    p      <- convert_to_number(dat[22,2:5]);
    z      <- convert_to_number(dat[24,2:5]);
    s      <- convert_to_number(dat[25,2:5]);

    birth  <- convert_to_number(dat[32,2:5])[1:2];
    lact   <- convert_to_number(dat[34,2:5])[1:2];
    heif   <- convert_to_number(dat[35,2:5])[1:2];
    adult  <- convert_to_number(dat[36,2:5])[1:2];

    rho    <- convert_to_number(dat[43,2:5])[1];
    iu     <- convert_to_number(dat[44,2:5])[1];
    qG     <- convert_to_number(dat[45,2:5])[1];
    hsp    <- convert_to_number(dat[46,2:5])[1];
    hsf    <- convert_to_number(dat[47,2:5])[1];
    end_t  <- convert_to_number(dat[48,2:5])[1];
  }

  eta <- 1000*k*f/iu;

  ghosts <- length(S[1,])

  if (seasonal == 0)
  {
    # I'm not afraid of no ghosts...
    num_ghosts <- init_a[ghosts]
    num_adults <- sum(init_a[3:4])
    init_a[4] <- 25; #round(0.15*(num_adults+num_ghosts));
    init_a[3] <- num_adults+num_ghosts-init_a[4]; # most lact?
    init_a[5] <- 50;
    cat(init_a, "\n")
#end_t <- 10*365

#    num_animals <- sum(init_a)
#    lt_heif <- init_a[2];
#    lt_calves <- round(heif[1]/365*lt_heif);
#    lt_adults <- num_animals - lt_heif - lt_calves
#    lt_lact <- round(lact[1]/365*lt_adults)
#    long_term_animals <- c(lt_calves, lt_heif, lt_lact, lt_adults - lt_lact)
#    cat("long term animals = ", long_term_animals, "\n")
#    init_a <- c(long_term_animals, 0)
#    num_births <- (m+b)*long_term_animals*365
#    cat("num_births = ", num_births, "sum=", sum(num_births), "needed=", (sum(long_term_animals[3:4])-sum(num_births)), "\n")
#    fbr <- (sum(long_term_animals[3:4])-sum(num_births))/sum(long_term_animals[3:4])
#    cat("fbr = ", fbr/365, "\n")
#    fb <- 0.5*c(0, 0, fbr, fbr)/365
#    num_births <- fb*long_term_animals
#    cat("num_births = ", num_births*365, "sum=", sum(num_births)*365, "\n")
#    init_p <- rep(0, length(init_p))
  }

  # make the parameters the correct length for ghosts
  m[ghosts] <- b[ghosts] <- beta[ghosts] <- gamma[ghosts] <- eta[ghosts] <- z[ghosts] <- s[ghosts] <- fb[ghosts] <- 0;
  q[ghosts] <- p[ghosts] <- 1;

  # combine all these into a single parameter dataframe
  params <<- list(vertical=rho,
                  culling=m,
                  death=b,
                  inf_direct=beta,
                  recovery=gamma,
                  avg_shedding=eta,
                  env_local_death=q,
                  env_global_death=qG,
                  env_transfer=p,
                  inf_local=z,
                  inf_global=s,
                  birth=birth,
                  lactation=lact,
                  maturation_heifer=heif,
                  maturation_adult=adult,
                  high_shedder_proportion=hsp,
                  high_shedder_factor=hsf,
                  init_animals=init_a,
                  init_prevalence=init_p,
                  shedding=matrix(0,nrow(S),ncol(S)),
                  end_time=end_t,
                  seasonal=seasonal,
                  fakebirth=fb);

  # and create the "fast" version of these

  n <- length(S[1,]);
  fast_params <<- matrix(0, 7, n);
  fast_params[1,] <<- params$death[1:n];
  fast_params[2,] <<- params$culling[1:n];
  fast_params[3,] <<- params$inf_direct[1:n];
  fast_params[4,] <<- params$inf_local[1:n];
  fast_params[5,] <<- params$inf_global[1:n];
  fast_params[6,] <<- params$recovery[1:n];
  fast_params[7,] <<- params$fakebirth[1:n];
}

interpolate_output <- function(times, num_vars, working)
{
  # read in our output data
  temp <- scan(working);
  dat <- matrix(temp, nrow=(num_vars+5));
  dat_values <- t(dat[6:(num_vars+5),]);
  dat_times  <- dat[1,];

  values <- matrix(0, length(times), num_vars);

  for (i in 1:num_vars)
  {
    # interpolate dat_values[,i] over dat_times on times
    values[,i] <- approx(dat_times, dat_values[,i], times)$y;
  }
  return(values);
}

interpolate_animals <- function(date_range, prev)
{
  times <- seq(from=0.5,to=params$end_time-0.5,by=2);

  # read in our output data
  dat_times  <- prev[,1];
  dat_values <- prev[,2:ncol(prev)];

  values <- matrix(0, length(times), ncol(dat_values));
  for (i in 1:ncol(dat_values))
  {
    # interpolate dat_values[,i] over dat_times on times
    values[,i] <- approx(dat_times, dat_values[,i], times)$y;
  }

  return(values);
}

interpolate_prev <- function(date_range, prev)
{
  times <- seq(from=0.5,to=params$end_time-0.5,by=2);

  # read in our output data
  dat_times  <- prev[,1];
  dat_values <- prev[,2:3];

  values <- matrix(0, length(times), 2);
  for (i in 1:ncol(dat_values))
  {
    # interpolate dat_values[,i] over dat_times on times
    values[,i] <- approx(dat_times, dat_values[,i], times)$y;
  }
  
  # now accumulate up to months - want to do this over multiple runs
  months <- matrix(0, nrow(date_range), 2);
  for (i in 1:nrow(date_range))
  	months[i,] <- colMeans(values[times >= date_range[i,1] & times < date_range[i,2],], na.rm = T)

  # set calves to 0?  TODO: is this sensible?
  for (i in 1:(nrow(date_range)/12))
  	months[(i-1)*12 + 4:12,1] <- 0;

  return(months);
}

analyse <- function(times, values, confidence, ylim, ylabel, plot_legend=FALSE)
{
  cols <- c("blue", "green", "cyan", "red", "purple")
  names <- c("Calves", "Heifers", "Lactating", "Dry")
  m  <- apply(values, 1, mean);
  lq <- apply(values, 1, quantile, (1 - confidence)/2, na.rm = TRUE);
  uq <- apply(values, 1, quantile, (1 + confidence)/2, na.rm = TRUE);

  m  <- matrix(m,  nrow=length(times))
  lq <- matrix(lq, nrow=length(times))
  uq <- matrix(uq, nrow=length(times))

  # now plot them
  xlim <- c(0, ceiling(tail(times,1)/365))
  if (plot_legend)
  {
    plot(times / 365, m[,1], ylim=ylim, xlim=xlim, xaxp=c(0,10,10), xlab="Year", ylab=ylabel, type="l", col=cols[1], log="y");
  }
  else
  {
    plot(times / 365, m[,1], ylim=ylim, xlim=xlim, xaxp=c(0,10,10), xlab="Year", ylab=ylabel, type="l", col=cols[1]);
  }
  for (i in 2:ncol(m))
    lines(times / 365, m[,i], col=cols[i]);
  for (i in 1:ncol(m))
    lines(times / 365, lq[,i], col=cols[i], lty=3);
  for (i in 1:ncol(m))
    lines(times / 365, uq[,i], col=cols[i], lty=3);
  if (plot_legend)
  {
    legend(xlim[1],ylim[2],names,col=cols, lty=1)
  }
}

multigroup <- function(theta=NULL, data=NULL, summary_stats=NULL, threshold=NULL, num_simulations=5)
{
#  full_path <- sprintf("%s/%s", getwd(), file)
  get_parameters(theta);
  # load the initial conditions
  get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

  # output values
  num_groups <- length(S[1,]);
  num_subgroups <- 1;
  num_vars <- num_groups * (num_subgroups*2 + 1) + 1; # S + I + E + G
  month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
  nyears <- 3
  date_range <- matrix(0, nyears*12, 2);
  date_range[1,] <- c(0, month_days[1]);
  for (i in 2:(nyears*12))
    date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);

  if (is.null(data))
  {
  	data <- load_data("EColiSeasonalPrevalence.txt")
  }
  rep_data <- NULL;
  for (i in 1:nyears)
    rep_data <- rbind(rep_data, data);

  totals <- rep_data[,3:4];
  data <- rep_data[,1:2];

  npowers <- nrow(summary_stats) / (nrow(rep_data)*2);

  data_set <- NULL
  for (power in 1:npowers)
	data_set <- c(data_set,rep_data[,1]^power,rep_data[,2]^power);

  ss_data <- colSums(data_set * summary_stats);

  dist <- rep(10000, num_simulations)
  prev_months <- matrix(0, nyears*12, 10);
  for (loop in 1:num_simulations)
  {
    cat("running simulation", loop, "of", num_simulations);

    get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

    prev <- run_model(params$end_time, temp);

    prev <- interpolate_prev(date_range, prev);

    sim_data <- NULL;
    for (power in 1:npowers)
      sim_data <- c(sim_data, prev[,1]^power, prev[,2]^power)

    # compute the summary stat
    diff <- colSums(sim_data * summary_stats) - ss_data;
    dist[loop] <- sum(diff[c(1,2,3,4,5,6,7)]^2)

    cat(" rho =", dist[loop], "\n");
    if (!is.null(threshold) && dist[loop] > threshold)
    	break;
    if (sum(prev == 0) == length(prev))
    {
    	cat("breaking due to zeros\n")
	dist <- rep(10000, num_simulations);
    	break;
    }
    prev_months[,(2*loop-1):(2*loop)] <- prev;
 
 #   for (i in 1:nrow(date_range))
 #   {
 #     cat(i+6, prev[i,], data[i,], (prev - data)[i,], ((prev - data)^2)[i,],"\n");
 #   }
 #   cat("ss = ", dist, "\n");
  }
  return(list(rho=dist, prev=prev_months))
}

plot_animals <- function(num_simulations=5, confidence = 0.9, col="grey60")
{
#  full_path <- sprintf("%s/%s", getwd(), file)
  get_parameters(c(0.1,0,0,0,0));
  # load the initial conditions
  get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

  params$end_time <<- 2*365

  # output values
  num_groups <- length(S[1,]);
  num_subgroups <- 1;
  num_vars <- num_groups * (num_subgroups*2 + 1) + 1; # S + I + E + G
  month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
  nyears <- 2
  date_range <- matrix(0, nyears*12, 2);
  date_range[1,] <- c(0, month_days[1]);
  for (i in 2:(nyears*12))
    date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);

  times <- seq(from=0.5,to=params$end_time-0.5,by=2);
  num_times <- length(times)

  p <- matrix(0, num_simulations, num_times*4)
  for (loop in 1:num_simulations)
  {
    cat("running simulation", loop, "of", num_simulations, "\n");

    get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

    prev <- run_animal_model(params$end_time);

    num_an <- interpolate_animals(date_range, prev);

    p[loop,] <- c(num_an[,1], num_an[,2], num_an[,3], num_an[,4])
  }
  mth_range <- 1:24
  d <- rowMeans(date_range[mth_range,])

#  par(mfrow=c(2,1))
#    pdf(sprintf("out%d.pdf", loop), width=8, height=5.5)
#  	plot(d, data[mth_range,loop], ylim=c(0,ylim), type="l", xlab="", ylab="Number of Animals", main="", xaxt="n")
#  	ylim <- (3-loop)*0.3
    ylim <- max(p, na.rm=T)
    plot(NULL, ylim=c(0,ylim), xlim=c(times[1], times[num_times]), type="l", xlab="", ylab="Prevelance", main="", xaxt="n")
    for (loop in 1:4)
    {
      m  <- apply(p[1:num_simulations,((loop-1)*num_times+1):(loop*num_times)], 2, median, na.rm = TRUE);
      lq <- apply(p[1:num_simulations,((loop-1)*num_times+1):(loop*num_times)], 2, quantile, (1 - confidence)/2, na.rm = TRUE);
      uq <- apply(p[1:num_simulations,((loop-1)*num_times+1):(loop*num_times)], 2, quantile, (1 + confidence)/2, na.rm = TRUE);
#    cat("m=",length(m), "data=", length(data[,loop]), "d=", length(d),"\n");
      lines(times, m, lty=loop, lwd=2)
      lines(times, lq, col="grey60", lty=loop)
      lines(times, uq, col="grey60", lty=loop)

    axis(1, at=c(date_range[,1], date_range[nrow(date_range),2]), labels=FALSE)
    for (i in 1:24)
    {
      months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
      mtext(months[(i+5) %% 12 + 1], side=1, at=d[i], cex=0.8);
    }
  }
  return(p)
}

plot_multigroup <- function(theta=NULL, data=NULL, num_simulations=5, confidence = c(0.6,0.9), col="black")
{
#  full_path <- sprintf("%s/%s", getwd(), file)
  get_parameters(theta);
  # load the initial conditions
  get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

  # output values
  num_groups <- length(S[1,]);
  num_subgroups <- 1;
  num_vars <- num_groups * (num_subgroups*2 + 1) + 1; # S + I + E + G
  month_days <- c(31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30); # starts with July
  nyears <- 3
  date_range <- matrix(0, nyears*12, 2);
  date_range[1,] <- c(0, month_days[1]);
  for (i in 2:(nyears*12))
    date_range[i,] <- c(date_range[i-1,2], date_range[i-1,2] + month_days[(i-1) %% 12 + 1]);

  if (is.null(data))
  {
  	data <- load_data("EColiSeasonalPrevalence.txt")
  }
  rep_data <- NULL;
  for (i in 1:nyears)
    rep_data <- rbind(rep_data, data);

  totals <- rep_data[,3:4];
  data <- rep_data[,1:2];

  p <- matrix(0, num_simulations, nrow(data)*2)
  for (loop in 1:num_simulations)
  {
    cat("running simulation", loop, "of", num_simulations, "\n");

    get_inits(params$init_animals, params$init_prevalence, params$high_shedder_proportion, params$high_shedder_factor);

    prev <- run_model(params$end_time, temp);

    month_prev <- interpolate_prev(date_range, prev);

    p[loop,] <- c(month_prev[,1], month_prev[,2])
  }
  mth_range <- 1:36
  d <- rowMeans(date_range[mth_range,])

  par(mfrow=c(2,1))
  for (loop in 1:2)
  {
#    pdf(sprintf("out%d.pdf", loop), width=8, height=5.5)
#  	plot(d, data[mth_range,loop], ylim=c(0,ylim), type="l", xlab="", ylab="Prevelance", main="", xaxt="n")
#  	ylim <- (3-loop)*0.3
  	ylim <- max(data[,loop],p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))])
  	plot(d, data[mth_range,loop], ylim=c(0,ylim), type="l", xlab="", ylab="Prevelance", main="", xaxt="n")
#    for (i in 1:num_simulations)
#    {
#      lines(d, p[i,((loop-1)*nrow(data)+1):(loop*nrow(data))], lty="dashed")
#    }
    m  <- apply(p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))], 2, median, na.rm = TRUE);
    lq <- apply(p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))], 2, quantile, (1 - confidence[1])/2, na.rm = TRUE);
    uq <- apply(p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))], 2, quantile, (1 + confidence[1])/2, na.rm = TRUE);
    llq <- apply(p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))], 2, quantile, (1 - confidence[2])/2, na.rm = TRUE);
    uuq <- apply(p[1:num_simulations,((loop-1)*nrow(data)+1):(loop*nrow(data))], 2, quantile, (1 + confidence[2])/2, na.rm = TRUE);
#    cat("m=",length(m), "data=", length(data[,loop]), "d=", length(d),"\n");
    lines(d, m, col=col)
    lines(d, lq, col=col, lty="dashed")
    lines(d, uq, col=col, lty="dashed")
    lines(d, llq, col=col, lty="dotted")
    lines(d, uuq, col=col, lty="dotted")
#  	for (sim in 1:num_simulations)
#  	{
#  		lines(d, p[sim,((loop-1)*nrow(data)+1):(loop*nrow(data))], lty="dashed", col="grey80")
#  	}
  	lines(d, data[,loop], type="l", lwd=2, col="black")	
    axis(1, at=c(date_range[,1], date_range[36,2]), labels=FALSE)
    for (i in 1:36)
    {
      months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
      mtext(months[(i+5) %% 12 + 1], side=1, at=d[i], cex=0.8);
    }
  }
}


load_data <- function(file = "EColiSeasonalPrevalence.txt")
{
 	# load our data in
  	data <- matrix(0, 12, 4)
  	d <- read.table(file, header=T)
  	for (t in 7:18)
  	{
  		ra <- d[d$Month == t & d$Type == "Adults",]
  		rc <- d[d$Month == t & d$Type == "Calves",]
  		ra2 <- d[d$Month == t+12 & d$Type == "Adults",]
  		rc2 <- d[d$Month == t+12 & d$Type == "Calves",]
  		pa <- pc <- ta <- tc <- 0;
  		if (nrow(ra) > 0)
  		{
  			pa <- pa + ra$Positive;
  			ta <- ta + ra$Total;
  		}
  		if (nrow(ra2) > 0)
  		{
  			pa <- pa + ra2$Positive;
  			ta <- ta + ra2$Total;
  		}
  		if (nrow(rc) > 0)
  		{
  			pc <- pc + rc$Positive;
  			tc <- tc + rc$Total;
  		}
  		if (nrow(rc2) > 0)
  		{
  			pc <- pc + rc2$Positive;
  			tc <- tc + rc2$Total;
  		}
  		if (ta > 0) { pa <- pa / ta }
  		if (tc > 0) { pc <- pc / tc }
  		data[t-6,] <- c(pc, pa, tc, ta)
  	}
  	return(data)
}

#dyn.load("exp_samp.so")
#require(xlsReadWrite)
