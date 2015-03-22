UpdateStatusFun <- function(alpha, gamma, current.state){
  # updates environmental state between steps in Markov population projection model
  #
  # Args:
  #
  # alpha = pathogen introduction rate
  # gamma = pathogen fade-out rate
  # current state = character indicating disease status in population at year t
  #
  # Returns:
  #
  # current.state.new = list containing character vector of length 1 which label of 
  # new state ("healthy", "spillover", "infected")
  current.state.new <- rep(NA, 1)
  if(current.state == "healthy"){
    gets.infected <- rbinom(1, 1, alpha)
    current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "spillover")
  }
  #  else if(current.state == "infected"){
  else if(current.state == "spillover"){
    gets.infected <- rbinom(1, 1, alpha)
    fade.out <- rbinom(1, 1, gamma)
    current.state.new[1] <- ifelse(gets.infected == 1, "spillover", ifelse((gets.infected == 0 & fade.out == 1), "healthy", "infected"))
  } else if(current.state == "infected"){
    gets.infected <- rbinom(1, 1, alpha)
    fade.out <- rbinom(1, 1, gamma)
    current.state.new[1] <- ifelse(gets.infected == 1, "spillover", ifelse((gets.infected == 0 & fade.out == 1), "healthy", "infected"))
  }
  return(list(current.state.new = current.state.new))
}