# store the first M_0 rows of the matrix
initial_matrix <- matrix[1:M_0,]

# find all infections
traj_values <- c()
index_values <- c()
time_values <- c()
person_values <- c()
date_values <- c()
positive_rows <- 1
for(row_ind in 1:nrow(initial_matrix)){
  if(any(initial_matrix[row_ind,] != 0)){
    traj <- unlist(initial_matrix[row_ind, which(initial_matrix[row_ind,] != 0)])
    traj_values <- append(traj_values, traj)
    index_values <- append(index_values, rep(positive_rows, length(traj)))
    time_values <- append(time_values, 1:length(traj))
    person_values <- append(person_values, rep(row_ind, length(traj)))
    date_values <- append(date_values, which(initial_matrix[row_ind,] != 0))
    
    positive_rows <- positive_rows + 1
  }
}

init_trajectories <- data.frame(index_values, time_values, traj_values, person_values,
                           date_values)
