library(WheresCroc)
# help("runWheresCroc")

# Probability of all holes based on readings ONLY
probHoles <- function(readingsList, probs) {
  # create empty vector
  prob_holes = matrix(0, nrow = 1, ncol = 40)
  for (i in seq_along(probs$salinity[,1])) {
    # for every lake
    prob_s = dnorm(readingsList$s, probs$salinity[i, 1], probs$salinity[i, 2]) # probability for s value
    prob_p = dnorm(readingsList$p, probs$phosphate[i, 1], probs$phosphate[i, 2]) # probability for p value
    prob_n = dnorm(readingsList$n, probs$nitrogen[i, 1], probs$nitrogen[i, 2]) # probability for n value
    prob = prob_s*prob_p*prob_n # total probability
    prob_holes[1, i] = prob # append it to the vector
  }
  return(prob_holes)  # return the vector containing all the probabilities.
}

# if a tourist die, their location will become negative,
# the turn after, their location will become NA

Death <- function(positions) {
  # Checks if a tourist was killed and returns where, if no kill, return 0
  for (i in range(1:2)){
    if (is.na(positions[[i]])) {
      next
    } else if (positions[[i]] < 0) {
      return(abs(positions[[i]]))
    }}
  return(0)
}

normalize <- function(vector) {
  total = sum(vector)
  return(1/total * vector)
}


visited <- function(position, previous_probs) {
  # Checks which waterholes were visited and sets them to 0
  for (i in 1:2) {
    if (is.na(position[i])) {
      next
    } 
    else {
      previous_probs[position[i]] = 0
    }
  }
  return(previous_probs)
}

transition <- function(newProbs, edges) {
  for (hole_index in seq_along(newProbs)) {
    adjacent = getOptions(hole_index, edges)
    prob = 1/length(adjacent)
    for (adj in adjacent) {
      newProbs[adj] = newProbs[adj] + newProbs[hole_index]*prob
    }
  }
  return(newProbs)
}

BFS <- function(ranger, edges, goal) {
  queue = c(ranger)
  expanded = c(ranger)
  parents = replicate(40, 0)
  parents[ranger] = -1
  while (length(queue) != 0) {
    curr_node = head(queue, n=1)
    queue = setdiff(queue, c(curr_node))
    adjacent = getOptions(curr_node, edges)
    adjacent = setdiff(adjacent, c(curr_node))
    adjacent = setdiff(adjacent, expanded)
    for (adj in adjacent) {
      if (!(adj %in% expanded)) {
        queue = c(queue, adj)
        parents[adj] = curr_node
        expanded = c(expanded, c(adj))
      }
    }
  }
  path = numeric()
  curr_node = goal
  while (curr_node != -1) {
    if(parents[curr_node] != -1) {
      path = c(c(curr_node), path)
    }
    curr_node = parents[curr_node]
  }
  return(path)
}

# Hidden Markov model
HMM <- function(positions, gameStatus, readingsList, probs, prevProbs, tMatrix, edges) {
  newProbs = matrix(0, nrow = 1, ncol = 40)
  death = Death(positions)
  if (death != 0) {
    newProbs[death] = 1
    newProbs = transition(newProbs, edges)
    
  } else if (is.null(gameStatus)) {
    newProbs = probHoles(readingsList, probs)
    newProbs = normalize(newProbs)
    newProbs = visited(positions, newProbs)
    newProbs = transition(newProbs, edges)
    
  } else {
    emission = probHoles(readingsList, probs)
    emission = normalize(emission)
    newProbs = prevProbs * emission
    newProbs = visited(positions, newProbs)
    newProbs = transition(newProbs, edges)
  }
  return(newProbs)
}

# About the function makeMoves (below):
# First argument: a list, made up of 1) a list called "moves" where we return where we want to go, and
# 2) a list called "mem" where we store information we want to remember from turn to turn.
# Second argument: a vector giving the salinity, phosphate and nitrogen reading from Croc's sensor at his current location.
# Third argument: a vector giving the position of the two tourists (elements one and two) and our own (third element).
# Fourth argument: a two-column matrix giving the edges paths between waterholes (edges) present.
# The numbers are from and to numbers for the waterholes. All edges can be crossed both ways.
# Fifth argument: a list of three matrices giving the mean and standard deviation of readings for
# salinity, phosphate and nitrogen respectively at each waterhole.
# Our function returns the first argument passed with an updated moves vector and some changes to the "mem" field we wish to keep.

# We make two moves per turn. Searching is considered a move.
# We receive readings each turn from Croc.

# readings gives us a list with [1] = salinity, 2 = [phosphate], 3 = [nitrogen], at the lake where Croc is.
# positions returns the positions of the people, [1] and [2], backpacker, [3] = ranger
# probs returns 3 matrices, $salinity, $phosphate, $nitrogen. each 40x2, where [1,1] = mean of chem in lake 1, [1,2] = std. dv of chem in lake 1
# edges is a matrix, [71, 2], which gives different combinations of how you can move.
# move info is for us to use
myFunction <- function(moveInfo, readings, positions, edges, probs) {
  ranger = positions[[3]]
  # Evaluate game status to know if we must start with creating the emission vector
  gameStatus = moveInfo$moves
  readingsList = list(s = readings[1], p = readings[2], n = readings[3])
  prevProbs = moveInfo$mem$holeProbs
  holeProbs = HMM(positions, gameStatus, readingsList, probs, prevProbs, tMatrix, edges)
  goal = which.max(holeProbs)
  moves = BFS(ranger, edges, goal)
  if (length(moves) >= 2) {
    moveInfo$moves = c(moves[1], moves[2])
  } 
  if (length(moves) == 1) {
    moveInfo$moves = c(moves[1], 0)
  } 
  if (length(moves) == 0) {
    moveInfo$moves = c(0,0)
  }
  moveInfo$mem$holeProbs = holeProbs
  return(moveInfo)
}

########### from package wheres croc!#############
# REALLY useful for the transformation
getOptions=function(point,edges) {
  return (c(edges[which(edges[,1]==point),2],edges[which(edges[,2]==point),1],point))
}
