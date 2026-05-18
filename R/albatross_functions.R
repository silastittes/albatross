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

mid_score <- function(dist_mat, cut_off = 0, group = colnames(dist_mat)[1], scale = F){
  if(scale){
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


#' Find optimal cutoff and mid_point for a single group
#'
#' The primary function - finds the exact cutoff value that maximises the mid_point
#' statistic (balanced accuracy) for a given group by sweeping over all observed
#' pairwise distances.  When multiple cutoffs yield the same mid_point the smallest
#' is returned, favouring tight groupings.
#' @import tidyverse
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param group A character string specifying group of interest. Default is first group in matrix.
#' @return tibble containing:
#' focal_group - the input group
#' cut_off - the smallest cutoff value that maximises mid_point.
#' mid_point - the maximum balanced accuracy: (proportion of within-group pairs <= cut_off + proportion of out-group pairs > cut_off) / 2. A value of 1 means the focal group is perfectly separated from all other groups.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Length %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' opt_mid(sepal_dist, "setosa")

opt_mid <- function(dist_mat, group = colnames(dist_mat)[1]){

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

  in_cells  <- which(colnames(dist_mat) == group)
  out_cells <- which(colnames(dist_mat) != group)

  in_dists  <- as.vector(dist_mat[in_cells,  in_cells])
  out_dists <- as.vector(dist_mat[out_cells, in_cells])

  n_in  <- length(in_dists)
  n_out <- length(out_dists)

  all_dists <- c(in_dists, out_dists)
  is_in     <- c(rep(TRUE, n_in), rep(FALSE, n_out))
  ord       <- order(all_dists)
  all_dists <- all_dists[ord]
  is_in     <- is_in[ord]

  intest   <- 0.0
  outtest  <- 1.0
  best_mid <- 0.5
  best_cut <- 0.0

  for(i in seq_along(all_dists)){
    if(is_in[i]){
      intest  <- intest  + 1 / n_in
    } else {
      outtest <- outtest - 1 / n_out
    }
    mid <- (intest + outtest) / 2
    if(mid > best_mid){
      best_mid <- mid
      best_cut <- all_dists[i]
    }
  }

  tibble(focal_group = group,
         cut_off     = best_cut,
         mid_point   = best_mid)
}


#' Runs opt_mid over multiple groups
#'
#' @import tidyverse
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param groups A character vector specifying groups of interest. Default is all groups.
#' @return tibble containing a row for each input group with columns:
#' focal_group - the input group
#' cut_off - the smallest cutoff value that maximizes mid_point.
#' mid_point - the maximum balanced accuracy for the group. A value of 1 means the focal group is perfectly separated from all other groups.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Width %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' opt_mid_multi(dist_mat = sepal_dist, groups = c("setosa", "virginica"))


opt_mid_multi <- function(dist_mat, groups = unique(colnames(dist_mat))){

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

  groups %>%
    map_df(~ opt_mid(dist_mat = dist_mat, group = .x))
}


#' Permutation of opt_mid_multi for non-parametric significance testing
#'
#' @import tidyverse
#' @import parallel
#' @param dist_mat A [0, 1] bounded square distance matrix with column names matching groups.
#' @param groups A character vector specifying groups of interest. Default is all groups.
#' @param n_permutes A positive integer for the number of random label permutations to be conducted.
#' @return tibble containing a row per group per permutation with columns:
#' focal_group - the input group
#' cut_off - the smallest cutoff value that maximizes mid_point for this permutation.
#' mid_point - the maximum balanced accuracy for this permutation.
#' iteration - the permutation index.
#' @export
#'
#' @examples
#' set.seed(123)
#' sepal_dist <- iris$Sepal.Width %>% dist %>% as.matrix
#' sepal_dist <- sepal_dist/max(sepal_dist)
#' colnames(sepal_dist) <- iris$Species
#' permute_fit(dist_mat = sepal_dist)

permute_fit <- function(dist_mat, groups = unique(colnames(dist_mat)), n_permutes = 10){

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

  bind_rows(
    1:n_permutes %>%
      mclapply(function(x){
        permute_dists <- dist_mat
        colnames(permute_dists) <- sample(colnames(dist_mat), replace = F)

        groups %>%
          map_df(~ opt_mid(group = .x, dist_mat = permute_dists)) %>%
          mutate(iteration = x)
      }))
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
  res = 3,
  scale = F){
  
  group %>% map_df(function(g){
    mids <- cut_offs %>% 
      map_dbl(~ mid_score(dist_mat = dist_mat, cut_off = .x, group = g, scale = scale)) %>%
      round(res)
    tibble(cut_offs = cut_offs, mids = mids, group = rep(g, length(cut_offs)))
  })
}

#' Create useful summary statistics from the output of man_multi
#'
#' @import tidyverse
#' @param cut_df The output of man_multi() (?man_multi() for further details)
#' @return A tibble contain a row for each group id returned by man_multi():
#' cut_off - The minimum cut off value of the the max mid value
#' plateau - The the number cut off values that with the same max mid value
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
  d_cut <- cut_df$cut_offs[2] - cut_df$cut_offs[1]
  cut_df %>%
    group_by(group) %>%
    summarise(cut_off = cut_offs[which.max(mids)],
              plateau = sum(mids == max(mids)),
              mid =  max(mids),
              mid_sum = sum(mids)*d_cut) %>%
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
