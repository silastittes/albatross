#' Calculate Midpoint
#'
#' For a given [0,1] cutoff value, this function calculates the midpoint between within group clustering and out-of-group exclusion.
#' It's primary purpose is internal use for other functions.
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param cut_off A [0, 1] bounded scalar that specifies cutoff. Defaults to 0.
#' @param group A character string specifying group of interest. Default if first group in matrix.
#' @return returns [0, 1] bounded midpoint
#' @export
#'
#' @examples
#' set.seed(123)
#' groups <- 5
#' n_groups <- 5
#' n_obs <- groups * n_groups
#' group_names <- rep(letters[1:groups], each = 5)
#' example_dist <- sample(1:100, size = 25) %>%
#'   dist %>%
#' as.matrix
#' example_dist <- example_dist/max(example_dist)
#' colnames(example_dist) <- group_names
#' mid_score(example_dist, cut_off = 0.2, group_names[1])

mid_score <- function(dist_mat, cut_off = 0, group = colnames(dist_mat)[1]){
  if(max(dist_mat) > 1){
    warning(
      paste0(
        "The input matrix was not [0,1] bounded. Dividing all values by max(",
        deparse(substitute(dist_mat)),
        ")."
        )
      )
    dist_mat <- dist_mat/max(dist_mat)
  }
  in_cells <- which(colnames(dist_mat) == group)
  out_cells <- which(colnames(dist_mat) != group)
  cutoff_test <- dist_mat <= cut_off
  intest <- mean(cutoff_test[in_cells, in_cells], na.rm = T)
  outtest <- 1 - mean(cutoff_test[out_cells, in_cells], na.rm = T)
  return( (intest+outtest)/2 )
}


#' Optimization function
#'
#' Used internally by opt_mid
#' @export

opt_cut <- function(params, dist_mat, group){
  cut_off <- params[1]
  delta <- params[2]
  in_cells <- which(colnames(dist_mat) == group)
  out_cells <- which(colnames(dist_mat) != group)
  cutoff_test <- dist_mat <= cut_off
  intest <- mean(cutoff_test[in_cells, in_cells], na.rm = T)
  outtest <- 1 - mean(cutoff_test[out_cells, in_cells], na.rm = T)
  return(log(cut_off^exp(delta) + (intest - outtest)^2))
}



#' Optimize mid point
#'
#' The primary function - estimates the optimal cutoff and related statistics for evaluating the groups of the input distance matrix.
#' @import tidyverse
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param group A character string specifying group of interest. Default if first group in matrix.
#' @param params A vector containing intitial values of cut_off and delta (respectively) passed optim function for parameter optimization.
#' @return tibble containing:
#' focal_group - the input group
#' cut_off - the optimized cutoff value that most effectively clusters entities with the same label while excluding entities with alternative labels.
#' delta - the optimized cutoff shrinkage parameter that ensures smaller cutoff values are favored when a range of cutoff values results in a similar mid_point value
#' optim_value - the optimized value of (cut_off ^ delta + (prop grouped - prop exluded)^2), which determines parameter values.
#' mid_point - the compromise between excluding all alternative entities from the focal group while capturing all the members of the focal group, as derived from the optimized cutoff value. A value of 1 means the focal group is perfectly separated from all other groups.
#' converge - convergence dianostic from optim function. If converge != 0, try different initial values. Increasing delta is a good first choice.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Length %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' opt_mid(sepal_dist, "setosa", params = c(0, 1))

opt_mid <- function(dist_mat, group = colnames(dist_mat)[1], params = c(0.0, 1.0)){

  if(max(dist_mat) > 1){
    warning(
      paste0(
        "The input matrix was not [0,1] bounded. Dividing all values by max(",
        deparse(substitute(dist_mat)),
        ")."
      )
    )
    dist_mat <- dist_mat/max(dist_mat)
  }

  opt_val <- optim(par = params,
                   fn = opt_cut,
                   dist_mat = dist_mat,
                   group = group,
                   method = "Nelder-Mead",
                   hessian = F)

  mid_point = mid_score(dist_mat = dist_mat, cut_off = opt_val$par[1], group = group)

  if(opt_val$convergence != 0){
    warning(
      paste(
        "optim() did not converge, try changing initial values of cut_off and delta from their current values of, ",
        params[1], params[2])
    )
  }

  return(data_frame(focal_group = group,
                    cut_off = opt_val$par[1],
                    delta = opt_val$par[2],
                    optim_value = opt_val$value,
                    mid_point = mid_point,
                    converge = opt_val$convergence))
}


#' Runs opt_mid over multiple groups
#'
#' @import tidyverse
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param groups A character vector specifying groups of interest. Defualt is all groups.
#' @param params A vector containing intitial values of cut_off and delta (respectively) passed optim function for parameter optimization.
#' @return tibble containing a row for each input group with columns:
#' focal_group - the input group
#' cut_off - the optimized cutoff value that most effectively clusters entities with the same label while excluding entities with alternative labels.
#' delta - the optimized cutoff shrinkage parameter that ensures smaller cutoff values are favored when a range of cutoff values results in a similar mid_point value
#' optim_value - the optimized value of (cut_off ^ delta + (prop grouped - prop exluded)^2), which determines parameter values.
#' mid_point - the compromise between excluding all alternative entities from the focal group while capturing all the members of the focal group, as derived from the optimized cutoff value. A value of 1 means the focal group is perfectly separated from all other groups.
#' converge - convergence dianostic from optim function. If converge != 0, try different initial values. Increasing delta is a good first choice.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Width %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' opt_mid_multi(dist_mat = sepal_dist, groups = c("setosa", "virginica"), params = c(0.4, 10))


opt_mid_multi <- function(dist_mat, groups = unique(colnames(dist_mat)), params = c(0.0, 1.0)){

  if(max(dist_mat) > 1){
    warning(
      paste0(
        "The input matrix was not [0,1] bounded. Dividing all values by max(",
        deparse(substitute(dist_mat)),
        ")."
      )
    )
    dist_mat <- dist_mat/max(dist_mat)
  }

  optim_df <- groups %>%
    map_df(~ {
      opt_mid(dist_mat = dist_mat, group = .x, params = params)
    })

  if(sum(optim_df$converge) > 0){
    warning(
      paste(
        "optim() did not converge, try changing initial values of cut_off and delta from their current values of, ",
        params[1], params[2])
    )
  }


  return(optim_df)
}


#' Permutation of opt_mid_multi for non parametric significance testing
#'
#' @import tidyverse
#' @import parallel
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param groups A character vector specifying groups of interest. Defualt is all groups.
#' @param params A vector containing intitial values of cut_off and delta (respectively) passed optim function for parameter optimization.
#' @param n_permutes A positive integer from the number of random label permutations to be conducted.
#' @return tibble containing a row for each input group with columns:
#' focal_group - the input group
#' cut_off - the optimized cutoff value that most effectively clusters entities with the same label while excluding entities with alternative labels.
#' delta - the optimized cutoff shrinkage parameter that ensures smaller cutoff values are favored when a range of cutoff values results in a similar mid_point value
#' optim_value - the optimized value of (cut_off ^ delta + (prop grouped - prop exluded)^2), which determines parameter values.
#' mid_point - the compromise between excluding all alternative entities from the focal group while capturing all the members of the focal group, as derived from the optimized cutoff value. A value of 1 means the focal group is perfectly separated from all other groups.
#' converge - convergence dianostic from optim function. If converge != 0, try different initial values. Increasing delta is a good first choice.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Width %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' permute_fit(dist_mat = sepal_dist, params = c(0.1, 10))

permute_fit <- function(dist_mat, groups = unique(colnames(dist_mat)), params = c(0.0, 1.0), n_permutes = 10){

  if(max(dist_mat) > 1){
    warning(
      paste0(
        "The input matrix was not [0,1] bounded. Dividing all values by max(",
        deparse(substitute(dist_mat)),
        ")."
      )
    )
    dist_mat <- dist_mat/max(dist_mat)
  }

  permute_df <- bind_rows(
    1:n_permutes %>%
      mclapply(function(x){
        permute_dists <- dist_mat
        permute_names <- sample(colnames(dist_mat), replace = F)
        colnames(permute_dists) <- permute_names
        rownames(permute_dists) <- permute_names

        x_permute_df <- groups %>%
          map_df(function(focal){
            opt_mid(group = focal, dist_mat = permute_dists, params)
          })
        x_permute_df %>% mutate(iteration = rep(x, n()))
      }))


  if(sum(permute_df$converge) > 0){
    warning(
      paste(
        "optim() did not converge for some permutations, try changing initial values of cut_off and delta from their current values of, ",
        params[1], params[2]
        )
    )
  }

  return(permute_df)
}

#' Calculate mid value statistic for each group in the distance matrix over a grid of cut offs 
#' 
#' @import tidyverse
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param cut_offs An ordered (low to high) sequence of cut_offs to calculate mid values.
#' @param group A character vector specifying groups of interest. Defualt is all groups.
#' @param res The number of digits to round mid values to. Helpful when plateaus occur.
#' 
#' @return A tibble containing the following columns:
#' cut_offs - [0,1] cut off values input by user
#' mids - the mid values for each cut off
#' group - the group id taken from the column names of the input distance matrix
#' @export
#' 
#' @examples 
#' sepal_dist <- iris$Sepal.Width %>% dist() %>% as.matrix()
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' cut_df <- man_multi(dist_mat = sepal_dist, res = 2)
#' ggplot(cut_df, aes(cut_offs, cut_df$mids, colour = group)) + 
#' geom_line()


man_multi <- function(
  dist_mat, 
  cut_offs = seq(0, 1, length.out = 100), 
  group = unique(colnames(dist_mat)), 
  res = 3){
  
  group %>% map_df(function(g){
    mids <- cut_offs %>% 
      map_dbl(~ mid_score(dist_mat = dist_mat, cut_off = .x, group = g)) %>%
      round(res)
    tibble(cut_offs = cut_offs, mids = mids, group = rep(g, length(cut_offs)))
  })
}

#' Create useful summary statistics from the output of man_multi
#'
#' @import tidyverse
#' @param cut_df The output of man_multi() (?man_multi() for further details)
#' @return A tibble contain a row for each group id returned by man_multi():
#' cut_off - The miniumum cut off value of the the max mid value
#' plateau - The the numer cut off values that with the same max mid value
#' mid - The max mid value (1.0 means perfectly grouped)
#' mid_sum - The sum of the mid values across all considered cut_offs, which is useful for quantifying how isolated the groups are
#' rel_sum - The mid_sum scores, scaled by the maximum mid_sum. 
#' @export
#' 
#' @examples 
#' sepal_dist <- iris$Sepal.Width %>% dist() %>% as.matrix()
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' cut_df <- man_multi(dist_mat = sepal_dist, res = 2)
#' summarise_multi(cut_df)

summarise_multi <- function(cut_df){
  cut_df %>%
    group_by(group) %>%
    summarise(cut_off = cut_offs[which.max(mids)],
              plateau = sum(mids == max(mids)),
              mid =  max(mids),
              mid_sum = sum(mids)) %>%
    mutate(rel_sum = mid_sum/max(mid_sum))
}


# summarise_multi <- function(cut_df){
#   cut_df %>% 
#     group_by(group) %>%
#     summarise(cut_off = ifelse( abs(min(mids) - 0.5) < abs(max(mids) - 0.5),
#                                cut_offs[which.max(mids)], 
#                                cut_offs[which.min(mids)]),
#               plateau = ifelse( abs(min(mids) - 0.5) < abs(max(mids) - 0.5), 
#                                sum(mids == max(mids)),
#                                sum(mids == min(mids))),
#               mid = ifelse( abs(min(mids) - 0.5) < abs(max(mids) - 0.5), 
#                            max(mids), min(mids)),
#               mid_sum = ifelse(abs(min(mids) - 0.5) < abs(max(mids) - 0.5), 
#                                sum(mids), sum(1 - mids))) %>%
#     mutate(rel_sum = mid_sum/max(mid_sum))
# }
