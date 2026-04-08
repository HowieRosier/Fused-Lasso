# R/00_prepare_data.R
#
# Run this ONCE before step1/step2.
# Reads raw cohort + PRS, merges them, and does a 60/20/20 train/test/vali split.

library(dplyr)

# ---------- merge cohort + PRS ----------
data <- read.csv("~/data/relax_prevent_cohort_time_prs.csv")
prs  <- read.delim("~/prs_top100.tsv")

# merge Data
prs$person_id <- prs$IID
data <- data %>%
  inner_join(prs %>%
               mutate(prs = PRS_STD), by = "person_id")
data <- data %>%
  select(-PRS_STD)

write.csv(data, file = "~/data/relax_prevent_cohort_time_prs.csv", row.names = FALSE)


# ---------- Train / test / vali split ----------
data <- read.csv("~/data/relax_prevent_cohort_time_prs.csv")

set.seed(123)
n <- nrow(data)
train_idx <- sample(seq_len(n), size = floor(0.60 * n))

train_data <- data[train_idx, ]
other_data <- data[-train_idx, ]
n1 <- nrow(other_data)
test_idx  <- sample(seq_len(n1), size = floor(0.50 * n1))
test_data <- other_data[test_idx, ]
vail_data <- other_data[-test_idx, ]

write.csv(train_data, file = "~/data/relax_prevent_cohort_time_prs_train.csv", row.names = FALSE)
write.csv(other_data, file = "~/data/relax_prevent_cohort_time_prs_test.csv",  row.names = FALSE)
write.csv(vail_data,  file = "~/data/relax_prevent_cohort_time_prs_vali.csv",  row.names = FALSE)

# Check sizes
cat("Train N:", nrow(train_data), "\n")
cat("Test  N:", nrow(test_data),  "\n")

# Check event proportions
cat("Train event rate:", mean(train_data$cvd), "\n")
cat("Test  event rate:", mean(test_data$cvd),  "\n")
