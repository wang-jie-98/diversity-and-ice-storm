# 1. readying --------------------------------------------------------------------------------------
# __1.1 loading packages ---------------------------------------------------------------------------
pacman::p_load(openxlsx, tidyverse, ggrepel, ggtree, ggimage, RColorBrewer, cowplot, rtrees, lme4, lmerTest, glmmTMB, emmeans, performance, MuMIn, piecewiseSEM)

# __1.2 loading functions ---------------------------------------------------------------------------
# function 01: z-transform for data
scale_z <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# function 02: exported the GLMM results
save_glmm <- function(index, count) {
  model_list <- str_c("model_", index, "_", str_pad(1:count, 2, side = "left", pad = "0"), sep = "")
  model_save <- NULL
  for (i in 1:n_distinct(model_list)) {
    model <- summary(get(model_list[i]))
    call <- as.character(model$call)
    mult_comp <- pairs(emmeans(get(model_list[i]), as.formula(str_c("~ ", as.character(model$call$formula[[3]][2]))), adjust = "tukey")) |> as.data.frame()
    ID <- c(1, nrow(model$coefficients$cond), length(as.character(model$varcor)[1]), nrow(mult_comp))
    model_save <- bind_rows(model_save, data.frame(index = index, 
                                                   formula = c(str_c(call[1], "(formula = ", call[2], ", data = ", call[3], ", family = ", call[4], ")"), rep(NA, max(ID) - 1)), 
                                                   fix_name = c(row.names(model$coefficients$cond), rep(NA, max(ID) - ID[2])), 
                                                   fix_estimate = c(model$coefficients$cond[, 1], rep(NA, max(ID) - ID[2])), 
                                                   fix_se = c(model$coefficients$cond[, 2], rep(NA, max(ID) - ID[2])), 
                                                   fix_z_value = c(model$coefficients$cond[, 3], rep(NA, max(ID) - ID[2])), 
                                                   fix_p_value = c(model$coefficients$cond[, 4], rep(NA, max(ID) - ID[2])), 
                                                   R2_conditional = c(r2(get(model_list[i]))[[1]], rep(NA, max(ID) - 1)), 
                                                   R2_marginal = c(r2(get(model_list[i]))[[2]], rep(NA, max(ID) - 1))
    ))
  }
  row.names(model_save) <- NULL
  return(model_save)
}

# function 03: bootstrap for slope
boot_slope <- function(x1, x2, y, n, data) {
  set.seed(1234); slope <- data.frame(slope1 = as.numeric(), slope2 = as.numeric())
  for(i in 1:n) {
    data_1 <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data_1 <- data[data_1, ]
    data_2 <- data_1[c(x1, x2, "alti_log", y, "n_total", "site")]
    names(data_2) <- c("x1", "x2", "alti_log", "y", "n_total", "site")
    model <- glmmTMB(cbind(y, n_total - y) ~ x1 + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_2, family = binomial)
    std_slope <- effectsize::standardize_parameters(model) |> tibble()
    slope[i, 1] <- std_slope[2, 2]
    model <- glmmTMB(cbind(y, n_total - y) ~ x2 + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_2, family = binomial)
    std_slope <- effectsize::standardize_parameters(model) |> tibble()
    slope[i, 2] <- std_slope[2, 2]
  }
  return(slope)
}

# function 04: exported the rlm results
save_rlm <- function(index, count) {
  model_list <- str_c("model_", index, "_", str_pad(1:count, 2, side = "left", pad = "0"), sep = "")
  model_save <- NULL
  for (i in 1:n_distinct(model_list)) {
    model <- summary(get(model_list[i]))
    p_value <- 2 * pt(-abs(model$coefficients[, 3]), model$df[2])
    call <- as.character(model$call)
    ID <- c(1, nrow(model$coefficients))
    data_1 <- data[index] |> as.vector()
    model_save <- bind_rows(model_save, data.frame(index = index, 
                                                   formula = c(str_c(call[1], "(formula = ", call[2], ", data = ", call[3], ")"), rep(NA, max(ID) - 1)), 
                                                   fix_name = c(row.names(model$coefficients), rep(NA, max(ID) - ID[2])), 
                                                   fix_estimate = c(model$coefficients[, 1], rep(NA, max(ID) - ID[2])), 
                                                   fix_se = c(model$coefficients[, 2], rep(NA, max(ID) - ID[2])), 
                                                   t_value = c(model$coefficients[, 3], rep(NA, max(ID) - ID[2])), 
                                                   t_p_value = c(p_value, rep(NA, max(ID) - ID[2])), 
                                                   Std_Coefficient = c(NA, model$coefficients[2, 1] * sd(data_1[[1]], na.rm = TRUE))
    ))
  }
  row.names(model_save) <- NULL
  return(model_save)
}

# __1.3 loading datasets ---------------------------------------------------------------------------
load(file = "data/data_P1.rdata")
load(file = "data/data_P2_Q1.rdata")
load(file = "data/data_P2_Q2.rdata")
load(file = "data/data_P3.rdata")

# 2. Part 1: Diversity resistance to freezing rain (FR) --------------------------------------------
data_P1 <- data_P1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))
# total damage      (DR_total):  n_damaged/n_total
# uprooting         (DR_uproot): n_uproot/n_total
# clear-bole broken (DR_below):  n_below/n_total
# crown broken      (DR_up):     n_up/n_total
# branch broken     (DR_branch): n_branch/n_total

# __2.1 Q1: Can taxonomic diversity resist the damage caused by FR? --------------------------------
model_save <- data.frame()

# ____2.1.1 Test 1: TD_SR --------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_TD_SR_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_01); r2(model_TD_SR_01); drop1(model_TD_SR_01, test = "Chi")

model_TD_SR_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_TD_SR_02); r2(model_TD_SR_02); drop1(model_TD_SR_02, test = "Chi")

model_TD_SR_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_03); r2(model_TD_SR_03); drop1(model_TD_SR_03, test = "Chi")

model_TD_SR_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_TD_SR_04); r2(model_TD_SR_04); drop1(model_TD_SR_04, test = "Chi")

model_TD_SR_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + TD_SR:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_05); r2(model_TD_SR_05); drop1(model_TD_SR_05, test = "Chi")

model_TD_SR_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_TD_SR_06); r2(model_TD_SR_06); drop1(model_TD_SR_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("TD_SR", 6))

# save model --
# save(model_save, file = "save/1_statistical result/model_Part1_Q1_data_20250703.rdata")
rm(model_save, model_TD_SR_01, model_TD_SR_02, model_TD_SR_03, model_TD_SR_04, model_TD_SR_05, model_TD_SR_06)

# __2.2 Q2: Can structural diversity resist the damage caused by FR? -------------------------------
model_save <- data.frame()

# ____2.2.1 Test 1: CV_TTH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_TTH_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_01); r2(model_CV_TTH_01); drop1(model_CV_TTH_01, test = "Chi")

model_CV_TTH_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_TTH_02); r2(model_CV_TTH_02); drop1(model_CV_TTH_02, test = "Chi")

model_CV_TTH_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_03); r2(model_CV_TTH_03); drop1(model_CV_TTH_03, test = "Chi")

model_CV_TTH_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_TTH_04); r2(model_CV_TTH_04); drop1(model_CV_TTH_04, test = "Chi")

model_CV_TTH_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + CV_TTH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_05); r2(model_CV_TTH_05); drop1(model_CV_TTH_05, test = "Chi")

model_CV_TTH_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_TTH_06); r2(model_CV_TTH_06); drop1(model_CV_TTH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_TTH", 6))

# ____2.2.2 Test 2: CV_DBH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_DBH_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_01); r2(model_CV_DBH_01); drop1(model_CV_DBH_01, test = "Chi")

model_CV_DBH_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_DBH_02); r2(model_CV_DBH_02); drop1(model_CV_DBH_02, test = "Chi")

model_CV_DBH_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_03); r2(model_CV_DBH_03); drop1(model_CV_DBH_03, test = "Chi")

model_CV_DBH_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_DBH_04); r2(model_CV_DBH_04); drop1(model_CV_DBH_04, test = "Chi")

model_CV_DBH_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + CV_DBH:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_05); r2(model_CV_DBH_05); drop1(model_CV_DBH_05, test = "Chi")

model_CV_DBH_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_DBH_06); r2(model_CV_DBH_06); drop1(model_CV_DBH_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_DBH", 6))

# ____2.2.3 Test 3: CV_NMB -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(2) DR_uproot ------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(3) DR_below -------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(4) DR_up ----------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# ______(5) DR_branch ------------------------------------------------------------------------------
model_CV_NMB_01 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_01); r2(model_CV_NMB_01); drop1(model_CV_NMB_01, test = "Chi")

model_CV_NMB_02 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model_CV_NMB_02); r2(model_CV_NMB_02); drop1(model_CV_NMB_02, test = "Chi")

model_CV_NMB_03 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_03); r2(model_CV_NMB_03); drop1(model_CV_NMB_03, test = "Chi")

model_CV_NMB_04 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model_CV_NMB_04); r2(model_CV_NMB_04); drop1(model_CV_NMB_04, test = "Chi")

model_CV_NMB_05 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + CV_NMB:alti_log + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_05); r2(model_CV_NMB_05); drop1(model_CV_NMB_05, test = "Chi")

model_CV_NMB_06 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model_CV_NMB_06); r2(model_CV_NMB_06); drop1(model_CV_NMB_06, test = "Chi")

model_save <- bind_rows(model_save, save_glmm("CV_NMB", 6))

# save model --
# save(model_save, file = "save/1_statistical result/model_Part1_Q2_data_20250703.rdata")
rm(model_save, 
   model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, model_CV_TTH_06, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, model_CV_DBH_06, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, model_CV_NMB_06)

# __2.3 Q3: Do taxonomic and structural diversity differ in their resistance to FR? ----------------
# check colinearity --
model_save <- data.frame()

# ____2.3.1 Test 1: TD_SR vs CV_TTH ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_damaged", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_uproot", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_below", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_up", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_TTH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_TTH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_TTH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_TTH", "n_branch", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_TTH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ____2.3.2 Test 2: TD_SR vs CV_DBH ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_damaged", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_uproot", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_below", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_up", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_DBH") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_DBH_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_DBH_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_DBH", "n_branch", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_DBH_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ____2.3.3 Test 3: TD_SR vs CV_NMB ----------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_damaged", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_total_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_total", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(2) DR_uproot ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_uproot", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_uproot_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_uproot", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(3) DR_below -------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_below, n_total - n_below) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_below, n_total - n_below) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_below", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_below_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_below", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(4) DR_up ----------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_up, n_total - n_up) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_up, n_total - n_up) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_up", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_up_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_up", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# ______(5) DR_branch ------------------------------------------------------------------------------
# (1) standard regression coefficient --
model_1 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
model_2 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)

model_3 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
model_4 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)

model_5 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
model_6 <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)

slope <- bind_rows(effectsize::standardize_parameters(model_1) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_2) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "ESF"), 
                   effectsize::standardize_parameters(model_3) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_4) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "LSF"), 
                   effectsize::standardize_parameters(model_5) |> tibble() |> filter(Parameter == "TD_SR") |> mutate(type = "OF"), 
                   effectsize::standardize_parameters(model_6) |> tibble() |> filter(Parameter == "CV_NMB") |> mutate(type = "OF")); slope

# (2) bootstrap method --
data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_P1 |> filter(forest_age == "ESF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_NMB_ESF.rdata")
slope_test_1 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_P1 |> filter(forest_age == "LSF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_NMB_LSF.rdata")
slope_test_2 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

data_slope <- boot_slope("TD_SR", "CV_NMB", "n_branch", 1000, data_P1 |> filter(forest_age == "OF"))
save(data_slope, file = "save/2_bootstrap slope/DR_branch_CV_NMB_OF.rdata")
slope_test_3 <- wilcox.test(data_slope$slope1, data_slope$slope2, paired = TRUE)

model_save <- bind_rows(model_save, bind_cols(slope |> mutate(Y = "DR_branch", .before = Parameter), boot_slope_p = c(slope_test_1$p.value, NA, slope_test_2$p.value, NA, slope_test_3$p.value, NA)))

# save model --
# save(model_save, file = "save/1_statistical result/model_Part1_Q3_data_20250703.rdata")
rm(data_slope, slope, slope_test_1, slope_test_2, slope_test_3, model_save, model_1, model_2, model_3, model_4, model_5, model_6)

# 3. Part 2: Species-level structural diversity in resistance to FR --------------------------------
# __3.1 Q1: Does the structural diversity of species resist FR? ------------------------------------
data <- data_P2_Q1 |> bind_rows() |> group_by(species_CN) |> filter(DR_total != 0) |> mutate(n_plot = n()) |> filter(n_plot > 5) |> ungroup()
model_save <- data.frame()
list_sp <- data |> distinct(species_CN)

# ____3.1.1 Test 1: CV_TTH -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_TTH)
  model_CV_TTH_01 <- MASS::rlm(DR_total ~ CV_TTH, data_1)
  model_CV_TTH_02 <- MASS::rlm(DR_uproot ~ CV_TTH, data_1)
  model_CV_TTH_03 <- MASS::rlm(DR_below ~ CV_TTH, data_1)
  model_CV_TTH_04 <- MASS::rlm(DR_up ~ CV_TTH, data_1)
  model_CV_TTH_05 <- MASS::rlm(DR_branch ~ CV_TTH, data_1)
  model_save <- bind_rows(model_save, save_rlm("CV_TTH", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# ____3.1.2 Test 2: CV_DBH -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_DBH)
  model_CV_DBH_01 <- MASS::rlm(DR_total ~ CV_DBH, data_1)
  model_CV_DBH_02 <- MASS::rlm(DR_uproot ~ CV_DBH, data_1)
  model_CV_DBH_03 <- MASS::rlm(DR_below ~ CV_DBH, data_1)
  model_CV_DBH_04 <- MASS::rlm(DR_up ~ CV_DBH, data_1)
  model_CV_DBH_05 <- MASS::rlm(DR_branch ~ CV_DBH, data_1)
  model_save <- bind_rows(model_save, save_rlm("CV_DBH", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# ____3.1.3 Test 3: CV_NMB -------------------------------------------------------------------------
for (i in 1:nrow(list_sp)) {
  data_1 <- data |> filter(species_CN == as.character(list_sp[i, 1])) |> drop_na(CV_NMB)
  model_CV_NMB_01 <- MASS::rlm(DR_total ~ CV_NMB, data_1)
  model_CV_NMB_02 <- MASS::rlm(DR_uproot ~ CV_NMB, data_1)
  model_CV_NMB_03 <- MASS::rlm(DR_below ~ CV_NMB, data_1)
  model_CV_NMB_04 <- MASS::rlm(DR_up ~ CV_NMB, data_1)
  model_CV_NMB_05 <- MASS::rlm(DR_branch ~ CV_NMB, data_1)
  model_save <- bind_rows(model_save, save_rlm("CV_NMB", 5) |> mutate(species = as.character(list_sp[i, 1]), .after = index))
}

# save model --
# save(model_save, file = "save/1_statistical result/model_Part2_Q1_data_20250703.rdata")
rm(model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, 
   i, data, data_1, list_sp)

# __3.2 Q2: What species traits determine the response to FR? --------------------------------------
model_save <- data.frame()

# ____3.2.1 Test 1: CV_TTH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_TTH"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_TTH"))

# ______(3) DR_below -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_TTH"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_TTH"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_TTH") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_TTH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_TTH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_TTH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_TTH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_TTH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_TTH"))

# ____3.2.2 Test 2: CV_DBH -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_DBH"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_DBH"))

# ______(3) DR_below -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_DBH"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_DBH"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_DBH") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_DBH_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_DBH_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_DBH_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_DBH_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_DBH_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_DBH"))

# ____3.2.3 Test 3: CV_NMB -------------------------------------------------------------------------
# ______(1) DR_total -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_total") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_total", index = "CV_NMB"))

# ______(2) DR_uproot ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_uproot") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_uproot", index = "CV_NMB"))

# ______(3) DR_below -------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_below") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_below", index = "CV_NMB"))

# ______(4) DR_up ----------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_up") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_up", index = "CV_NMB"))

# ______(5) DR_branch ------------------------------------------------------------------------------
data <- data_P2_Q2 |> filter(index == "CV_NMB") |> filter(damage_type == "DR_branch") |> mutate(across(c(woody_density:NB), ~ scale_z(.x)))
model_CV_NMB_01 <- lm(Std_Coefficient ~ NB, data)
model_CV_NMB_02 <- lm(Std_Coefficient ~ woody_density, data)
model_CV_NMB_03 <- lm(Std_Coefficient ~ P50, data)
model_CV_NMB_04 <- lm(Std_Coefficient ~ rdmax, data)
model_CV_NMB_05 <- lm(Std_Coefficient ~ height, data)

model_all <- lm(Std_Coefficient ~ NB + woody_density + P50 + rdmax + height, data)
model_set <- dredge(model_all, trace = 2, options(na.action = "na.fail"))
model_avg <- model.avg(model_set, delta < 4) |> summary()
model_save <- bind_rows(model_save, model_avg$coefmat.full |> as.data.frame() |> mutate(Y = "DR_branch", index = "CV_NMB"))

# save model --
# model_save <- model_save |> rownames_to_column("fix_name") |> mutate(fix_name = str_remove(fix_name, "\\...[0-9]*"))
# save(model_save, file = "save/1_statistical result/model_Part2_Q2_data_20250703.rdata")
rm(model_CV_TTH_01, model_CV_TTH_02, model_CV_TTH_03, model_CV_TTH_04, model_CV_TTH_05, 
   model_CV_DBH_01, model_CV_DBH_02, model_CV_DBH_03, model_CV_DBH_04, model_CV_DBH_05, 
   model_CV_NMB_01, model_CV_NMB_02, model_CV_NMB_03, model_CV_NMB_04, model_CV_NMB_05, 
   model_all, model_set, model_avg, model_save, data)

# 4. Part 3: The potential mechanism of structural diversity in resisting FR -----------------------
# __4.1 structural equation model ------------------------------------------------------------------
model_save <- data.frame()

# ____4.1.1 Test 1: DR_total -----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_P3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> select(site, forest_age, alti_log, CWM_WD:CWM_H, DR_total:CV_NMB) |> 
  mutate(across(c(alti_log:CV_NMB), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_total ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_total ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Structural diversity --
model <- glmmTMB(DR_total ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_total ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Functional identity --
model <- glmmTMB(DR_total ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_total ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$forest_age, site = data$site, DR = data$DR_total, alti1 = data$alti1, altiq = data$altiq, comp = data$TD_SR, stru = data$stru, func = data$func)

# part model within ESF --
data_sem1 <- data_sem |> filter(stan == "ESF")
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ altiq + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan == "LSF")
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan == "OF")
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_total"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))
rm(data, data_P1, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.2 Test 2: DR_uproot ----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_P3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> select(site, forest_age, alti_log, CWM_WD:CWM_H, DR_total:CV_NMB) |> 
  mutate(across(c(alti_log:CV_NMB), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_uproot ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_uproot ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Structural diversity --
model <- glmmTMB(DR_uproot ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_uproot ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Functional identity --
model <- glmmTMB(DR_uproot ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_uproot ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$forest_age, site = data$site, DR = data$DR_uproot, alti1 = data$alti1, altiq = data$altiq, comp = data$TD_SR, stru = data$stru, func = data$func)

# part model within ESF --
data_sem1 <- data_sem |> filter(stan == "ESF")
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan == "LSF")
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan == "OF")
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_uproot"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))
rm(data, data_P1, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.3 Test 3: DR_below -----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_P3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> select(site, forest_age, alti_log, CWM_WD:CWM_H, DR_total:CV_NMB) |> 
  mutate(across(c(alti_log:CV_NMB), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_below ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_below ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Structural diversity --
model <- glmmTMB(DR_below ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_below ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Functional identity --
model <- glmmTMB(DR_below ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_below ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$forest_age, site = data$site, DR = data$DR_below, alti1 = data$alti1, altiq = data$altiq, comp = data$TD_SR, stru = data$stru, func = data$func)

# part model within ESF --
data_sem1 <- data_sem |> filter(stan == "ESF")
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ altiq + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan == "LSF")
model_sem_2 <- psem(glmmTMB(DR ~ altiq + stru + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan == "OF")
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_below"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))
rm(data, data_P1, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.4 Test 4: DR_up --------------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_P3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> select(site, forest_age, alti_log, CWM_WD:CWM_H, DR_total:CV_NMB) |> 
  mutate(across(c(alti_log:CV_NMB), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_up ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_up ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Structural diversity --
model <- glmmTMB(DR_up ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_up ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Functional identity --
model <- glmmTMB(DR_up ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_up ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$forest_age, site = data$site, DR = data$DR_up, alti1 = data$alti1, altiq = data$altiq, comp = data$TD_SR, stru = data$stru, func = data$func)

# part model within ESF --
data_sem1 <- data_sem |> filter(stan == "ESF")
model_sem_1 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan == "LSF")
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan == "OF")
model_sem_3 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ altiq + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_up"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))
rm(data, data_P1, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_1, model_sem_2, model_sem_3)

# ____4.1.5 Test 5: DR_branch ----------------------------------------------------------------------
# ______(1) build the dataset ----------------------------------------------------------------------
data <- data_P3 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude)) |> select(site, forest_age, alti_log, CWM_WD:CWM_H, DR_total:CV_NMB) |> 
  mutate(across(c(alti_log:CV_NMB), ~ scale_z(.x)))

# (1) construct composite variable Altitude --
model <- glmmTMB(DR_branch ~ poly(alti_log, 2, raw = TRUE) + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$alti1 <- data$alti_log
data$altiq <- -(summary(model)$coefficients$cond[2, 1]*data$alti_log + 
                 summary(model)$coefficients$cond[3, 1]*data$alti_log*data$alti_log)

model <- glmmTMB(DR_branch ~ altiq + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (2) construct composite variable Structural diversity --
model <- glmmTMB(DR_branch ~ CV_TTH + CV_DBH + CV_NMB + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$stru <- -(summary(model)$coefficients$cond[2, 1]*data$CV_TTH + 
                summary(model)$coefficients$cond[3, 1]*data$CV_DBH + 
                summary(model)$coefficients$cond[4, 1]*data$CV_NMB)

model <- glmmTMB(DR_branch ~ stru + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# (3) construct composite variable Functional identity --
model <- glmmTMB(DR_branch ~ CWM_WD + CWM_P50 + CWM_RDmax + CWM_H + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

data$func <- (summary(model)$coefficients$cond[2, 1]*data$CWM_WD +
                summary(model)$coefficients$cond[3, 1]*data$CWM_P50 +
                summary(model)$coefficients$cond[4, 1]*data$CWM_RDmax +
                summary(model)$coefficients$cond[5, 1]*data$CWM_H)

model <- glmmTMB(DR_branch ~ func + (1 | site), data)
summary(model); effectsize::standardize_parameters(model)

# ______(2) implement the SEM ----------------------------------------------------------------------
data_sem <- data.frame(stan = data$forest_age, site = data$site, DR = data$DR_branch, alti1 = data$alti1, altiq = data$altiq, comp = data$TD_SR, stru = data$stru, func = data$func)

# part model within ESF --
data_sem1 <- data_sem |> filter(stan == "ESF")
model_sem_1 <- psem(glmmTMB(DR ~ comp + stru + func + (1 | site), data_sem1), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem1), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem1), 
                    func %~~% stru); summary(model_sem_1)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_1)$coefficients[, -9] |> mutate(type = "ESF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_1)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_1)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_1)$coefficients) - nrow(summary(model_sem_1)$R2)))) |> select(!del)))

# part model within LSF --
data_sem2 <- data_sem |> filter(stan == "LSF")
model_sem_2 <- psem(glmmTMB(DR ~ altiq + comp + stru + func + (1 | site), data_sem2), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem2), 
                    glmmTMB(func ~ altiq + comp + (1 | site), data_sem2), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem2), 
                    func %~~% stru); summary(model_sem_2)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_2)$coefficients[, -9] |> mutate(type = "LSF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_2)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_2)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_2)$coefficients) - nrow(summary(model_sem_2)$R2)))) |> select(!del)))

# part model within OF --
data_sem3 <- data_sem |> filter(stan == "OF")
model_sem_3 <- psem(glmmTMB(DR ~ comp + stru + func + (1 | site), data_sem3), 
                    glmmTMB(stru ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(func ~ alti1 + comp + (1 | site), data_sem3), 
                    glmmTMB(comp ~ alti1 + (1 | site), data_sem3), 
                    func %~~% stru); summary(model_sem_3)

model_save <- bind_rows(model_save, bind_cols(summary(model_sem_3)$coefficients[, -9] |> mutate(type = "OF", DR = "DR_branch"), 
                                              bind_rows(summary(model_sem_3)$ChiSq, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$Cstat, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - 1))) |> select(!del), 
                                              bind_rows(summary(model_sem_3)$R2, data.frame(del = rep(NA, nrow(summary(model_sem_3)$coefficients) - nrow(summary(model_sem_3)$R2)))) |> select(!del)))

model_save <- model_save |> set_names("Response_1", "Predictor", "Estimate", "Std_Error", "DF", "Crit_Value", "P_Value_1", "Std_Estimate", "Type", "DR", "Chisq", "df_1", "P_Value_2", "Fisher_C", "df_2", "P_Value_3", "Response_2", "family", "link", "method", "Marginal", "Conditional")
# save(model_save, file = "save/1_statistical result/model_Part3_Q1_data_20250703.rdata")
rm(data, data_P1, data_P3, data_sem, data_sem1, data_sem2, data_sem3, model, model_sem_1, model_sem_2, model_sem_3, model_save)
