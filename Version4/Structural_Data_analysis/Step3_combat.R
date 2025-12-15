# 基于多中心合并后的数据进行站点矫正
# 还没有完全自动化，需要自己一个一个输入指标

# Load necessary libraries
library(ComBatFamily)
library(dplyr)
library(mgcv)

# Define paths and parameters for structural data harmonization
project_path <- '/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version4/Structural_data/' # nolint
data_path <- paste(project_path, sep = "")

metric <- 'ci'  # can also be 'gv', 'ct', or 'sa' # nolint: quotes_linter.
# we only want to harmonize 'artifact' version
qc_version <- 'pass'  # can be 'noqc', 'artifact' # nolint: quotes_linter.
pfactor_include <- 'no' # whether to include p-factor in covariates # nolint: quotes_linter.

# Construct the input filename and read the data
if (pfactor_include == 'yes'){
  dtype <- sprintf('%s_%s_pfactor_filter', metric, qc_version)
}
if (pfactor_include == 'no'){
  dtype <- sprintf('%s_%s', metric, qc_version)
}

combined_data <- read.csv(paste(data_path, sprintf('hbn_df_%s.tsv', dtype),  # combined_pnc_hbn_df_%s.tsv
                                sep = ""),
                          sep = '\t')

# Extract covariates and features
age_vec <- combined_data$age
sex_vec <- as.factor(combined_data$sex)
euler_vec <- combined_data$euler
if (pfactor_include == 'yes'){
pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
}

if (metric == 'ct') {
  structural_data <- combined_data[, c(2:402)] #[, c(1:68, 71, 72)]
} else if (metric == 'sa') {
  structural_data <- combined_data[, c(2:402)]
} else {
  structural_data <- combined_data[, c(2:402)]
}
# else if (metric == 'sv') {
#  structural_data <- combined_data[, c(1:15)]
# }
# structural_data = combined_data[0:73] # assumes features are in first 400 columns # nolint: line_length_linter.

# Prepare covariate dataframe 纳入年龄和性别等因素的意义。在矫正站点效应时，避免删除上述因素与脑影像特征的关联。使得更精准的矫正站点效应相关的误差
# Harmonize using covfam with GAM (including nonlinear age effects)

if (pfactor_include == 'yes'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(euler_vec),
                        as.numeric(pfactor_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            euler = ...4,
                            pfactor = ...5)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=structural_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + euler + pfactor)
}

if (pfactor_include == 'no'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(euler_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            euler = ...4)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=structural_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + euler)
}

# -------------------------- 关键修改：合并 participant_id --------------------------
# 1. 提取covar_df中的participant_id（与矫正数据的行顺序完全一致，确保匹配）
participant_ids <- covar_df %>% select(participant_id)
# 2. 将participant_id与矫正后的数据合并（ID列放在第一列，方便查看）
data.harmonized_covbat <- bind_cols(participant_ids, data.frame(data.harmonized$dat.covbat))
# ----------------------------------------------------------------------------------

# Save the harmonized functional data
outputPath <- paste(data_path, sprintf('hbn_df_%s_harmonized.tsv', dtype), #combined_pnc_'
                    sep = "")
write.csv(data.harmonized_covbat, outputPath, row.names = FALSE, sep = '\t')