l_0 <- 5
remaining_matrix <- matrix[(M_0 + 1):length(matrix),]
test_times <- seq(from = 1, to = ncol(remaining_matrix), by = l_0)

traj_values <- c()
index_values <- c()
time_values <- c()
person_values <- c()
date_values <- c()
positive_rows <- 1

for(row_ind in 1:nrow(remaining_matrix)){
  if(any(remaining_matrix[row_ind, test_times] > 0)){
    row <- remaining_matrix[row_ind,test_times]
    traj_indices <- which(row > 0)
    if(max(traj_indices) > length(row) - 2){
      traj_indices <- c(min(traj_indices) - 1, traj_indices)
    }
    else if(min(traj_indices) < 4){
      traj_indices <- c(traj_indices, max(traj_indices) + 1)
    }
    else{
      traj_indices <- c(min(traj_indices) - 1, traj_indices, max(traj_indices) + 1) 
    }
    traj <- unlist(row[traj_indices])
    traj_values <- append(traj_values, traj)
    index_values <- append(index_values, rep(positive_rows, length(traj)))
    time_values <- append(time_values, (seq(1, length(traj)) - 1) * l_0)
    person_values <- append(person_values, rep(row_ind, length(traj)))
    date_values <- append(date_values, traj_indices * l_0)
    
    positive_rows <- positive_rows + 1
  }
}

partial_trajectories <- data.frame(index_values, time_values, traj_values, person_values,
                                  date_values)
