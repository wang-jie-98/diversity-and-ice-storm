# 00. readying -------------------------------------------------------------------------------------
# loading packages --
pacman::p_load(tidyverse, ggtree, ggimage, rtrees, RColorBrewer, cowplot, glmmTMB, performance)
# pacman::p_load(openxlsx, tidyverse, ggtree, ggspatial, ggimage, RColorBrewer, cowplot, rtrees, rgl, lme4, glmmTMB, performance)

# 01. Figure. 01 -----------------------------------------------------------------------------------
# use PPT for drawing --

# 02. Figure. 02 -----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_P1.rdata")
data_P1 <- data_P1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))

# __2.1 Fig. 2A ------------------------------------------------------------------------------------
# TD_SR vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(TD_SR = seq(1, 5, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_01 <- ggplot() + 
  geom_ribbon(data = data, aes(x = TD_SR, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = TD_SR, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(1, 3, 5), labels = c(2, 8, 32)) +
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = "SR", y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" = 0.023"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.681"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# __2.2 Fig. 2B ------------------------------------------------------------------------------------
# CV_TTH vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_02 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_TTH, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_TTH, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(18, 40, 62)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = expression(CV[H]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.565"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# __2.3 Fig. 2C ------------------------------------------------------------------------------------
# CV_DBH vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_03 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_DBH, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_DBH, y = fit, color = type, linetype = type)) + 
  scale_x_continuous(breaks = c(35, 85, 135)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid")) + 
  labs(x = expression(CV[DBH]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" = 0.330"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

# __2.4 Fig. 2D ------------------------------------------------------------------------------------
# CV_NMB vs DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF"))

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data <- bind_rows(data, effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF"))

data <- data |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_04 <- ggplot() + 
  geom_ribbon(data = data, aes(x = CV_NMB, y = fit, ymin = lower, ymax = upper, fill = type), alpha = 0.2) + 
  geom_line(data = data, aes(x = CV_NMB, y = fit, color = type)) + 
  scale_x_continuous(breaks = c(50, 150, 250)) + 
  scale_y_continuous(limits = c(5, 55), breaks = c(10, 30, 50)) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  labs(x = expression(CV[NMB]*" (%)"), y = NULL, fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.94), label = expression("OF:   "*italic(P)*" < 0.001"),  size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.86), label = expression("LSF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") +
  annotate(geom = "text", x = ggpp::as_npc(0.50), y = ggpp::as_npc(0.78), label = expression("ESF: "*italic(P)*" < 0.001"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") +
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_04

# __2.5 Fig. 2E-G ----------------------------------------------------------------------------------
load(file = "save/2_bootstrap slope/DR_total_CV_TTH_ESF.rdata"); data_1 <- data_slope
load(file = "save/2_bootstrap slope/DR_total_CV_DBH_ESF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope3"))
load(file = "save/2_bootstrap slope/DR_total_CV_NMB_ESF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope4"))

load(file = "save/2_bootstrap slope/DR_total_CV_TTH_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> set_names("slope5", "slope6"))
load(file = "save/2_bootstrap slope/DR_total_CV_DBH_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope7"))
load(file = "save/2_bootstrap slope/DR_total_CV_NMB_LSF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope8"))

load(file = "save/2_bootstrap slope/DR_total_CV_TTH_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> set_names("slope9", "slope10"))
load(file = "save/2_bootstrap slope/DR_total_CV_DBH_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope11"))
load(file = "save/2_bootstrap slope/DR_total_CV_NMB_OF.rdata"); data_1 <- bind_cols(data_1, data_slope |> select(slope2) |> set_names("slope12"))

data_1 <- data_1 |> pivot_longer(cols = slope1:slope12, names_to = "slope", values_to = "value") |> mutate(n = rep(4:1, 3000)) |> 
  mutate(X = case_when(slope %in% c("slope1", "slope5", "slope9") ~ "TD_SR", slope %in% c("slope2", "slope6", "slope10") ~ "CV_TTH", slope %in% c("slope3", "slope7", "slope11") ~ "CV_DBH", TRUE ~ "CV_NMB")) |> 
  mutate(forest_age = case_when(slope %in% c("slope1", "slope2", "slope3", "slope4") ~ "young", slope %in% c("slope5", "slope6", "slope7", "slope8") ~ "middle", TRUE ~ "mature")) |> 
  mutate(x = case_when((forest_age == "young") & (n == 1) ~ 1, (forest_age == "young") & (n == 2) ~ 2, (forest_age == "young") & (n == 3) ~ 3, (forest_age == "young") & (n == 4) ~ 4, 
                       (forest_age == "middle") & (n == 1) ~ 6, (forest_age == "middle") & (n == 2) ~ 7, (forest_age == "middle") & (n == 3) ~ 8, (forest_age == "middle") & (n == 4) ~ 9, 
                       (forest_age == "mature") & (n == 1) ~ 11, (forest_age == "mature") & (n == 2) ~ 12, (forest_age == "mature") & (n == 3) ~ 13, (forest_age == "mature") & (n == 4) ~ 14))

load(file = "save/1_statistical result/model_Part1_Q3_data_20250703.rdata"); data_2 <- model_save |> tibble(); rm(model_save)
data_2 <- data_2 |> filter(Y == "DR_total") |> mutate(n = case_when(Parameter == "TD_SR" ~ 4, Parameter == "CV_TTH" ~ 3, Parameter == "CV_DBH" ~ 2, TRUE ~ 1)) |> 
  mutate(forest_age = case_when(str_detect(type, "ESF") ~ "young", str_detect(type, "LSF") ~ "middle", TRUE ~ "mature")) |> 
  mutate(x = case_when((forest_age == "young") & (n == 1) ~ 1, (forest_age == "young") & (n == 2) ~ 2, (forest_age == "young") & (n == 3) ~ 3, (forest_age == "young") & (n == 4) ~ 4, 
                       (forest_age == "middle") & (n == 1) ~ 6, (forest_age == "middle") & (n == 2) ~ 7, (forest_age == "middle") & (n == 3) ~ 8, (forest_age == "middle") & (n == 4) ~ 9, 
                       (forest_age == "mature") & (n == 1) ~ 11, (forest_age == "mature") & (n == 2) ~ 12, (forest_age == "mature") & (n == 3) ~ 13, (forest_age == "mature") & (n == 4) ~ 14))

plot_05 <- ggplot() + 
  gghalves::geom_half_violin(data = data_1, aes(x = x, y = value, group = x, fill = X, color = X), side = "r", alpha = 0.3, linewidth = 0.25) + 
  geom_errorbar(data = data_2, aes(x = x, ymin = CI_low, ymax = CI_high, color = Parameter), alpha = 0.5, width = 0, linewidth = 1, show.legend = FALSE) + 
  geom_point(data = data_2, aes(x = x, y = Std_Coefficient, color = Parameter), alpha = 0.5, shape = 16, size = 1.5, show.legend = FALSE) + 
  
  geom_rect(aes(xmin = 4.7, xmax = 5.3, ymin = -0.5, ymax = 0.5), fill = "#8cc269", alpha = 0.5) + 
  geom_segment(aes(x = 0.5, xend = 4.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 5, y = 0, label = "ESF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 1, y = 0.40, xend = 4, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 2.5, y = 0.45, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  
  geom_rect(aes(xmin = 9.7, xmax = 10.3, ymin = -0.5, ymax = 0.5), fill = "#1ba784", alpha = 0.5) + 
  geom_segment(aes(x = 5.3, xend = 9.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 10, y = 0, label = "LSF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 6, y = 0.40, xend = 9, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 7.5, y = 0.45, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 7, y = 0.32, xend = 9, yend = 0.32), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 8, y = 0.37, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 8, y = 0.24, xend = 9, yend = 0.24), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 8.5, y = 0.29, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  
  geom_rect(aes(xmin = 14.7, xmax = 15.3, ymin = -0.5, ymax = 0.5), fill = "#1a6840", alpha = 0.5) + 
  geom_segment(aes(x = 10.3, xend = 14.7, y = 0, yend = 0), linetype = 5, linewidth = 0.4, alpha = 0.7) + 
  annotate("text", x = 15, y = 0, label = "OF", color = "#000000", size = (9*0.35), family = "serif") +
  geom_segment(aes(x = 11, y = 0.40, xend = 14, yend = 0.40), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 12.5, y = 0.45, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 12, y = 0.32, xend = 14, yend = 0.32), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 13, y = 0.37, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  geom_segment(aes(x = 13, y = 0.24, xend = 14, yend = 0.24), color = "#000000", alpha = 0.5, linewidth = 0.5) + 
  annotate("text", x = 13.5, y = 0.29, label = "***", color = "#000000", size = (10*0.35), family = "serif", angle = 90) +
  
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, 15.3), breaks = c(14, 13, 12, 11, 9, 8, 7, 6, 4, 3, 2, 1), labels = c("SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]), "SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]), "SR", expression(CV[H]), expression(CV[DBH]), expression(CV[NMB]))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 0.5), breaks = c(-0.4, 0, 0.4), labels = c("−0.4", "0", "0.4")) + 
  scale_fill_manual(values = c("#43b244", "#43b244", "#43b244", "#1781b5")) + 
  scale_color_manual(values = c("#43b244", "#43b244", "#43b244", "#1781b5")) + 
  labs(x = NULL, y = "Standardized slopes") + 
  coord_flip() + 
  theme_bw(base_family = "serif") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        strip.text = element_blank(),  
        legend.position = "none"); plot_05

# __2.6 layout -------------------------------------------------------------------------------------
plot_06 <- ggplot() + scale_x_continuous(expand = c(0, 0), limits = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0), panel.grid = element_blank(), panel.background = element_rect(fill = "#ffffff")) + 
  annotate(geom = "text", x = 0.3, y = 5.5, label = "Total damage ratio (%)", size = (10*0.35), family = "serif", angle = 90); plot_06

plot_grid(plot_grid(plot_06, plot_06, nrow = 2), 
  plot_grid(plot_01, plot_02, nrow = 2) + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.98), label = "(a)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.48), label = "(b)", size = (10*0.35), family = "serif"), 
  plot_grid(plot_03, plot_04, nrow = 2) + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.98), label = "(c)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.05), y = ggpp::as_npc(0.48), label = "(d)", size = (10*0.35), family = "serif"), 
  plot_05, rel_widths = c(0.6, 5.7, 5.7, 6.0), nrow = 1) + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.975), label = "(e)", size = (10*0.35), family = "serif") + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.670), label = "(f)", size = (10*0.35), family = "serif") + 
  annotate(geom = "text", x = ggpp::as_npc(0.68), y = ggpp::as_npc(0.363), label = "(g)", size = (10*0.35), family = "serif")
# ggsave(file = "save/fig_2_20260829.tiff", width = 16, height = 11, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

rm(plot_01, plot_02, plot_03, plot_04, plot_05, plot_06, data, data_1, data_2, data_P1, data_slope, model)

# 03. Figure. 03 -----------------------------------------------------------------------------------
# __3.1 Fig. 3A ------------------------------------------------------------------------------------
# data processing --
load(file = "data/data_P2_Q1.rdata")
load(file = "save/1_statistical result/model_Part2_Q1_data_20250703.rdata")

data <- model_save |> filter(!is.na(Std_Coefficient)) |> mutate(damage_type = rep(c("DR_total", "DR_uproot", "DR_below", "DR_up", "DR_branch"), 105)) |> filter(str_detect(fix_name, "CV")) |> 
  select(index, species, damage_type, fix_estimate, Std_Coefficient) |> left_join(bind_rows(data_P2_Q1) |> distinct(species_CN, species_LN), by = c("species" = "species_CN"))

sp_list <- sp_list_df(sp_list = unique(data$species_LN), taxon = "plant")
set.seed(1234); tree_plant <- get_tree(sp_list = sp_list, taxon = "plant", scenario = "random_below_basal", show_grafted = FALSE)
data.frame(x = tree_plant$tip.label, y = 1:35)

data <- data |> mutate(species_LN = case_when(species_LN == "Corylus_ferox_var._thibetica" ~ "Corylus_ferox", species_LN == "Quercus_aliena_var._acutiserrata" ~ "Quercus_aliena", species_LN == "Acer_pictum_subsp._mono" ~ "Acer_pictum", species_LN == "Cornus_kousa_subsp._chinensis" ~ "Cornus_kousa", TRUE ~ species_LN)) |> 
  left_join(data.frame(species_LN = c("Carpinus turczaninowii", "Carpinus fargesiana", "Corylus ferox", "Corylus chinensis", "Betula albosinensis", "Betula luminifera", "Juglans mandshurica", "Platycarya strobilacea", "Quercus serrata", "Quercus aliena", "Quercus variabilis", "Castanea seguinii", "Castanea henryi", "Fagus engleriana", "Sorbus folgneri", "Sorbus alnifolia", "Crataegus wilsonii", "Prunus conradinae", "Prunus padus", "Populus lasiocarpa", "Populus davidiana", "Salix wallichiana", "Rhus potaninii", "Rhus chinensis", "Toxicodendron vernicifluum", "Acer davidii", "Acer pictum", "Cornus macrophylla", "Cornus kousa", "Cornus controversa", "Fraxinus chinensis", "Diospyros lotus", "Litsea ichangensis", "Lindera obtusiloba", "Pinus armandii"), 
                       ID = c(35:1), niche_width = c(0.5693252, 0.5067005, 0.4052766, 0.8033091, 0.7300763, 0.6806774, 0.7154775, 0.7134962, 0.4579634, 0.5269582, 0.5795512, 0.6247113, 0.539175, 0.4583443, 0.4789034, 0.4000346, 0.5557265, 0.5170674, 0.6401926, 0.5663812, 0.6025009, 0.6959552, 0.7151405, 0.7137618, 0.4638151, 0.7080873, 0.6139691, 0.516735, 0.3869133, 0.6158935, 0.4840941, 0.6366028, 0.5825325, 0.4570297, 0.7790024)) |> mutate(species_LN = str_replace_all(species_LN, " ", "_"))) |> arrange(ID)

# phylogenetic heatmap --
plot_01 <- ggplot(tree_plant, branch.length = "none") + 
  geom_tree(linewidth = 0.25) + 
  geom_text(data = data |> distinct(species_LN, ID) |> mutate(species_LN = str_replace_all(species_LN, "_", " ")), aes(x = 11, y = ID, label = species_LN), size = (7*0.35), hjust = "left", fontface = "italic", family = "serif") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_TTH"), aes(x = 30, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_TTH"), aes(x = 30, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_DBH"), aes(x = 33, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_DBH"), aes(x = 33, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_NMB"), aes(x = 36, y = ID, color = Std_Coefficient), shape = 15, size = (10*0.35), stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_total") |> filter(index == "CV_NMB"), aes(x = 36, y = ID, label = ifelse(Std_Coefficient > 0, "+", "−")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 30.0, y = 0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 33.0, y = 0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 36.0, y = 0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 40, hjust = 1) + 
  annotate(geom = "text", x = 3, y = -0.5, label = "Slopes", size = (9*0.35), family = "serif") + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "RdBu")[10:5])(99), breaks = c(-0.15, -0.05, 0.05), labels = c("−0.15", "−0.05", "0.05")) + 
  scale_x_continuous(limits = c(-0.4, 38)) + 
  scale_y_continuous(limits = c(-1.4, 36)) + 
  theme_bw(base_family = "serif") + 
  labs(x = NULL, y = NULL, color = NULL) + 
  theme(legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        legend.direction = "horizontal", 
        legend.position = c(0.40, 0.06),
        legend.background = element_blank(),
        legend.text = element_text(size = 9, margin = margin(t = 1, unit = "pt")),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.margin = margin(t = -22, r = 0, b = -16, l = -9)); plot_01

# __3.2 Fig. 3B ------------------------------------------------------------------------------------
# data processing --
load(file = "save/1_statistical result/model_Part2_Q2_data_20250703.rdata")
data <- data |> filter(fix_name != "(Intercept)") |> filter(index %in% c("CV_TTH", "CV_DBH", "CV_NMB")) |> filter(Y == "DR_total") |> 
  mutate(fix_name = case_when(fix_name == "NB" ~ "1_NB", fix_name == "woody_density" ~ "2_WD", fix_name == "P50" ~ "3_P50", fix_name == "rdmax" ~ "4_RDmax", TRUE ~ "5_Hmax")) |> arrange(index, fix_name) |> 
  mutate(ID = c(14, 16, 18, 20, 22, 26, 28, 30, 32, 34, 2, 4, 6, 8, 10)) |> set_names("fix_name", "estimate", "se", "se_adjusted", "z_palue", "p_value", "Y", "index", "ID")

# forest plot --
plot_02 <- ggplot(data, aes(x = estimate, y = ID)) + 
  geom_rect(aes(xmin = -0.024, xmax = 0.044, ymin =  0.5, ymax = 11.5), fill = "#8cc269", alpha = 0.008) + 
  geom_rect(aes(xmin = -0.024, xmax = 0.044, ymin = 12.5, ymax = 23.5), fill = "#1ba784", alpha = 0.008) + 
  geom_rect(aes(xmin = -0.024, xmax = 0.044, ymin = 24.5, ymax = 35.5), fill = "#1a6840", alpha = 0.008) + 
  
  geom_errorbarh(aes(xmax = estimate + 1.96 * se, xmin = estimate - 1.96 * se, color = index), height = 0, linewidth = 3, alpha = 0.5) + 
  geom_point(aes(color = index), size = 2.5, alpha = 0.7, shape = 16) + 
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 35.5), linetype = 2, linewidth = 0.4) + 

  annotate(geom = "text", x = -0.023, y = 10.8, label = expression(CV[H]), size = (10*0.35), family = "serif", fontface = "bold", color = "#8cc269", hjust = "left") + 
  annotate(geom = "text", x = -0.015, y = 2, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = -0.009, y = 2, label = expression(""^"*"), size = (10*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 4, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 6, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 8, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 10, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  
  annotate(geom = "text", x = -0.023, y = 22.8, label = expression(CV[DBH]), size = (10*0.35), family = "serif", fontface = "bold", color = "#1ba784", hjust = "left") + 
  annotate(geom = "text", x = -0.015, y = 14, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = -0.009, y = 14, label = expression(""^"*"), size = (10*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 16, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 18, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 20, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 22, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  
  annotate(geom = "text", x = -0.023, y = 34.8, label = expression(CV[NMB]), size = (10*0.35), family = "serif", fontface = "bold", color = "#1a6840", hjust = "left") + 
  annotate(geom = "text", x = -0.015, y = 26, label = "NB", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = -0.009, y = 26, label = expression(""^"**"), size = (10*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 28, label = "WD", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 30, label = "P50", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 32, label = "RDmax", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = 0.025, y = 34, label = "Hmax", size = (9*0.35), family = "serif", hjust = "left") + 
  scale_x_continuous(limits = c(-0.024, 0.044), breaks = c(-0.02, 0, 0.04), labels = c("−0.02", "0", "0.04")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 36)) + 
  scale_fill_manual(values = c("#1ba784", "#1a6840", "#8cc269")) + 
  scale_color_manual(values = c("#1ba784", "#1a6840", "#8cc269")) + 
  labs(x = "Standardized effect sizes", y = NULL) + 
  theme_classic(base_family = "serif") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        axis.line.y  = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "none", 
        plot.margin = margin(t = 0, r = 2, b = 2, l = 0)); plot_02

# __3.3 layout -------------------------------------------------------------------------------------
plot_grid(plot_01, plot_02, nrow = 1, rel_widths = c(0.52, 0.48)) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.98), label = "(a)", size = (10*0.35), family = "serif") + 
  annotate(geom = "text", x = ggpp::as_npc(0.51), y = ggpp::as_npc(0.98), label = "(b)", size = (10*0.35), family = "serif")
# ggsave(file = "save/fig_3_20260829.tiff", width = 12, height = 12, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

rm(data, data_P2_Q1, model_save, sp_list, tree_plant, plot_01, plot_02)

# 04. Figure. 04 -----------------------------------------------------------------------------------
# __4.1 Fig. 4A ------------------------------------------------------------------------------------
# use PPT for drawing --

# __4.2 Fig. 4B ------------------------------------------------------------------------------------
# use PPT for drawing --

# __4.3 Fig. 4C ------------------------------------------------------------------------------------
# use PPT for drawing --

# __4.4 Fig. 4D ------------------------------------------------------------------------------------
plot_01 <- data.frame(n = c(1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15), 
                      type = c("young", "young", "young", "young", "middle", "middle", "middle", "middle", "mature", "mature", "mature", "mature"), 
                      effect = c(-0.448, -0.003, -0.168, 0.082, -0.482, -0.121, -0.219, 0.153, -0.168, -0.199, -0.188, 0.111)) |> 
  ggplot(aes(x = n, y = effect, fill = type)) + 
  geom_col(position = "stack", width = 0.9, alpha = 0.5) + 
  geom_hline(yintercept = 0, linewidth = 0.2) + 
  
  geom_vline(xintercept = 4, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 
  geom_vline(xintercept = 8, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 
  geom_vline(xintercept = 12, linewidth = 0.3, linetype = 5, color = "#7f7f7f") + 

  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.20, ymax = -0.24), fill = "#c5e0b4") + 
  annotate(geom = "text", x = 13.6, y = -0.22, label = "ESF", size = (10*0.35), family = "serif", hjust = "left") + 
  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.28, ymax = -0.32), fill = "#8dd3c1") + 
  annotate(geom = "text", x = 13.6, y = -0.30, label = "LSF", size = (10*0.35), family = "serif", hjust = "left") + 
  geom_rect(aes(xmin = 12.5, xmax = 13.3, ymin = -0.36, ymax = -0.40), fill = "#8cb39f") + 
  annotate(geom = "text", x = 13.6, y = -0.38, label = "OF", size = (10*0.35), family = "serif", hjust = "left") + 
  
  scale_fill_manual(values = c("#1a6840", "#1ba784", "#8cc269"), labels = c("OF", "LSF", "ESF")) + 
  scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15), labels = c("", "Elevation", "", "", "Species\nrichness", "", "", "Structural\ndiversity", "", "", "Structural\nidentity", "")) + 
  scale_y_continuous(limits = c(-0.49, 0.2), breaks = c(-0.4, -0.2, 0, 0.2), labels = c("−0.4", "−0.2", "0", "0.2")) + 
  labs(x = NULL, y = "Total effects", fill = NULL) + 
  theme_bw(base_family = "serif") + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)); plot_01
# ggsave(file = "save/fig_4_20260829.tiff", width = 7.6, height = 6.4, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

rm(plot_01)

# 05. Figure. S1 -----------------------------------------------------------------------------------
png_1 <- "data/icon/Uprooting.png"
png_2 <- "data/icon/Clear-bole broken.png"
png_3 <- "data/icon/Crown broken.png"
png_4 <- "data/icon/Branch broken.png"
load(file = "data/data_P1.rdata")

data_1 <- data_P1 |> filter(forest_age == "ESF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_1$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_1

data_2 <- data_P1 |> filter(forest_age == "LSF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_2$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_2

data_3 <- data_P1 |> filter(forest_age == "OF") |> mutate(n_uproot = sum(n_uproot), n_below = sum(n_below), n_up = sum(n_up), n_branch = sum(n_branch), n_well = sum(n_total - n_damaged)) |> slice_head(n = 1) |> 
  select(n_uproot:n_branch, n_well) |> t() |> data.frame() |> rownames_to_column() |> set_names("damage_type", "n") |> mutate(ratio = round(n/sum(n)*100, 2)) |> select(!n)
data_3$damage_type <- c("1_Uprooting", "2_Clear-bole broken", "3_Crown broken", "4_Branch broken", "5_Undamaged"); data_3

data <- bind_rows(data_1 |> mutate(type = "Early secondary-growth forest"), data_2 |> mutate(type = "Late secondary-growth forest"), data_3 |> mutate(type = "Old-growth forest")) |> mutate(ratio = ifelse(damage_type == "5_Undamaged", 100 - ratio, ratio))

ggplot(data, aes(x = damage_type, y = ratio, fill = type)) + 
  geom_col(position = "dodge", width = 0.8, alpha = 0.5) + 
  scale_x_discrete(labels = c("Uprooting", "Clear-bole\nbroken", "Crown\nbroken", "Branch\nbroken", "Total")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 27), breaks = c(5, 15, 25)) +
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  geom_image(aes(x = 1, y = 24), image = png_1, size = 0.19) + 
  geom_image(aes(x = 2, y = 24), image = png_2, size = 0.19) + 
  geom_image(aes(x = 3, y = 24), image = png_3, size = 0.19) + 
  geom_image(aes(x = 4, y = 24), image = png_4, size = 0.19) + 
  labs(x = NULL, y = "Proportion (%)", fill = "Forest stand categories") + 
  theme_bw(base_family = "serif") + 
  theme(panel.grid = element_blank(), 
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.5, "cm"), 
        legend.position = c(0.297, 0.6), 
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2), 
        axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10))
# ggsave(file = "save/fig_S1_20260829.tiff", width = 10, height = 8, units = "cm", dpi = 300, limitsize = FALSE)

rm(data, data_1, data_2, data_3, data_P1, png_1, png_2, png_3, png_4)

# 06. Figure. S2 -----------------------------------------------------------------------------------
# use PPT for drawing --

# 07. Figure. S3 -----------------------------------------------------------------------------------
ggplot() + scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + scale_y_continuous(expand = c(0, 0), limits = c(0.8, 4.8)) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0), panel.grid = element_blank(), panel.background = element_rect(fill = "#ffffff")) + 
  
  geom_rect(aes(xmin = 0.5, xmax = 1.9, ymin = 4.1, ymax = 4.7), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 1.2, y = 4.4, label = "Elevation", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 3.1, xmax = 4.5, ymin = 4.1, ymax = 4.7), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 3.8, y = 4.4, label = "Species\nrichness", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 0.5, xmax = 1.9, ymin = 2.5, ymax = 3.1), fill = "#7f7f7f", alpha = 0.5) + 
  geom_polygon(aes(x = c(0.5, 0.7, 0.5), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(0.5, 0.7, 0.5), y = c(2.8, 3.1, 3.1)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(1.7, 1.9, 1.9), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(1.7, 1.9, 1.9), y = c(3.1, 3.1, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_text(aes(x = 1.2, y = 2.8, label = "Structural\ndiversity", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 3.1, xmax = 4.5, ymin = 2.5, ymax = 3.1), fill = "#7f7f7f", alpha = 0.5) + 
  geom_polygon(aes(x = c(3.1, 3.3, 3.1), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(3.1, 3.3, 3.1), y = c(2.8, 3.1, 3.1)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(4.3, 4.5, 4.5), y = c(2.5, 2.5, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_polygon(aes(x = c(4.3, 4.5, 4.5), y = c(3.1, 3.1, 2.8)), fill = "#ffffff", color = "#ffffff") + 
  geom_text(aes(x = 3.8, y = 2.8, label = "Structural\nidentity", family = "serif"), size = (9*0.35)) + 
  
  geom_rect(aes(xmin = 1.8, xmax = 3.2, ymin = 0.9, ymax = 1.5), fill = "#7f7f7f", alpha = 0.5) + 
  geom_text(aes(x = 2.5, y = 1.2, label = "Damage\nratio", family = "serif"), size = (9*0.35)) + 
  
  # 01 line
  geom_segment(aes(x = 0.2, y = 4.4, xend = 0.5, yend = 4.4), linewidth = 0.2) + 
  geom_segment(aes(x = 0.2, y = 4.4, xend = 0.2, yend = 1.2), linewidth = 0.2) + 
  geom_segment(aes(x = 0.2, y = 1.2, xend = 1.8, yend = 1.2), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 0.2, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 0.2, y = 3.75, label = expression(italic("1")), size = (9*0.35), family = "serif") + 
  
  # 02 line
  geom_segment(aes(x = 1.9, y = 4.4, xend = 3.1, yend = 4.4), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.5, y0 = 4.4, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.5, y = 4.4, label = expression(italic("2")), size = (9*0.35), family = "serif") + 
  
  # 03 line
  geom_segment(aes(x = 1.2, y = 4.1, xend = 1.2, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 1.2, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 1.2, y = 3.75, label = expression(italic("3")), size = (9*0.35), family = "serif") + 
  
  # 04 line
  geom_segment(aes(x = 1.2, y = 4.1, xend = 3.8, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.11, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.11, y = 3.75, label = expression(italic("4")), size = (9*0.35), family = "serif") + 
  
  # 05 line
  geom_segment(aes(x = 3.8, y = 4.1, xend = 1.2, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 2.89, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.89, y = 3.75, label = expression(italic("5")), size = (9*0.35), family = "serif") + 
  
  # 06 line
  geom_segment(aes(x = 3.8, y = 4.1, xend = 3.8, yend = 3.1), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 3.8, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 3.8, y = 3.75, label = expression(italic("6")), size = (9*0.35), family = "serif") + 
  
  # 07 line
  geom_segment(aes(x = 4.8, y = 4.4, xend = 4.5, yend = 4.4), linewidth = 0.2) + 
  geom_segment(aes(x = 4.8, y = 4.4, xend = 4.8, yend = 1.2), linewidth = 0.2) + 
  geom_segment(aes(x = 4.8, y = 1.2, xend = 3.2, yend = 1.2), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 4.8, y0 = 3.75, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 4.8, y = 3.75, label = expression(italic("7")), size = (9*0.35), family = "serif") + 
  
  # 08 line
  geom_segment(aes(x = 1.2, y = 2.5, xend = 2.5, yend = 1.5), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 1.85, y0 = 2.0, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 1.85, y = 2.0, label = expression(italic("8")), size = (9*0.35), family = "serif") + 
  
  # 09 line
  geom_segment(aes(x = 3.8, y = 2.5, xend = 2.5, yend = 1.5), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.2) + 
  ggforce::geom_circle(aes(x0 = 3.15, y0 = 2.0, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 3.15, y = 2.0, label = expression(italic("9")), size = (9*0.35), family = "serif") + 
  
  # 10 line
  geom_curve(aes(x = 1.9, y = 2.8, xend = 3.1, yend = 2.8), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.25, curvature = 0.3) + 
  geom_curve(aes(x = 3.1, y = 2.8, xend = 1.9, yend = 2.8), arrow = arrow(angle = 15, length = unit(0.1, "cm"), type = "closed"), linewidth = 0.25, curvature = -0.3) + 
  ggforce::geom_circle(aes(x0 = 2.5, y0 = 2.6, r = 0.15), fill = "#ffffff", alpha = 0.8, linewidth = 0.1) + 
  annotate(geom = "text", x = 2.5, y = 2.6, label = expression(italic("10")), size = (9*0.35), family = "serif")

# ggsave(file = "save/fig_S3_20260829.tiff", width = 8, height = 6, units = "cm", dpi = 300, limitsize = FALSE)

# 08. Figure. S4 -----------------------------------------------------------------------------------
load(file = "data/data_P1.rdata")
data_P1 <- data_P1 |> mutate(TD_SR = log2(TD_SR), alti_log = log2(altitude))

# OF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 4.25, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(16.3, 62.3, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29, 139, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "OF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "OF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(41.9, 184, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "OF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_01 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 4.4), breaks = c(0, 2.2, 4.4)) +
  scale_y_continuous(limits = c(14, 40), breaks = c(14, 27, 40)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Old-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression("−0.166"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression("−0.175"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression("−0.195"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.201"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 4), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# LSF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1, 4.39, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(17.2, 56.8, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(31.5, 119, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "LSF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "LSF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(39.7, 252, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "LSF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_02 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 4.4), breaks = c(0, 2.2, 4.4)) +
  scale_y_continuous(limits = c(10, 40), breaks = c(10, 25, 40)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Late secondary-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression("−0.076"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression("−0.141"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression("−0.179"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.201"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 4), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# ESF --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ TD_SR + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(TD_SR), max = max(TD_SR))
data1 <- effects::allEffects(model, xlevels = list(TD_SR = seq(1.58, 4.75, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data1 <- data1 |> mutate(fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_TTH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_TTH), max = max(CV_TTH))
data2 <- effects::allEffects(model, xlevels = list(CV_TTH = seq(19.1, 55.8, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data2 <- data2 |> mutate(CV_TTH = CV_TTH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_DBH + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_DBH), max = max(CV_DBH))
data3 <- effects::allEffects(model, xlevels = list(CV_DBH = seq(29.4, 126, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data3 <- data3 |> mutate(CV_DBH = CV_DBH/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ CV_NMB + poly(alti_log, 2, raw = TRUE) + (1 | site), data = data_P1 |> filter(forest_age == "ESF"), family = binomial)
summary(model); r2(model); data_P1 |> filter(forest_age == "ESF") |> summarise(min = min(CV_NMB), max = max(CV_NMB))
data4 <- effects::allEffects(model, xlevels = list(CV_NMB = seq(46.7, 248, 0.01)))[[1]] |> as.data.frame() |> mutate(type = "ESF")
data4 <- data4 |> mutate(CV_NMB = CV_NMB/100, fit = 100*fit, lower = 100*lower, upper = 100*upper)

plot_03 <- ggplot() + 
  geom_line(data = data1, aes(x = TD_SR, y = fit, color = type, linetype = type), color = "#f43e06") + 
  geom_line(data = data2, aes(x = CV_TTH, y = fit, color = type, linetype = type), color = "#8cc269") + 
  geom_line(data = data3, aes(x = CV_DBH, y = fit, color = type, linetype = type), color = "#1ba784") + 
  geom_line(data = data4, aes(x = CV_NMB, y = fit, color = type, linetype = type), color = "#1a6840") + 
  scale_x_continuous(limits = c(0, 5), breaks = c(0, 2.5, 5)) +
  scale_y_continuous(limits = c(9, 25), breaks = c(9, 17, 25)) + 
  labs(x = "Diversity metrics", y = "Total damage ratio (%)", fill = NULL, color = NULL) + 
  annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.94), label = "Early secondary-growth forest", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.73), y = ggpp::as_npc(0.94), label = "Standardized slopes", size = (9*0.35), family = "serif", hjust = "left") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.87), label = expression(SR*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.87), label = expression(" 0.019"), size = (9*0.35), family = "serif", hjust = "left", color = "#f43e06") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.80), label = expression(CV[H]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.80), label = expression(" 0.023"), size = (9*0.35), family = "serif", hjust = "left", color = "#8cc269") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.73), label = expression(CV[DBH]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.73), label = expression(" 0.044"), size = (9*0.35), family = "serif", hjust = "left", color = "#1ba784") + 
  annotate(geom = "text", x = ggpp::as_npc(0.75), y = ggpp::as_npc(0.66), label = expression(CV[NMB]*":"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  annotate(geom = "text", x = ggpp::as_npc(0.85), y = ggpp::as_npc(0.66), label = expression("−0.190"), size = (9*0.35), family = "serif", hjust = "left", color = "#1a6840") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        plot.margin = margin(t = 2, r = 2, b = 2, l = 4), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

plot_grid(plot_01, plot_02, plot_03, nrow = 3) +  
    annotate(geom = "text", x = ggpp::as_npc(0.022), y = ggpp::as_npc(0.990), label = "(a)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.022), y = ggpp::as_npc(0.655), label = "(b)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.022), y = ggpp::as_npc(0.322), label = "(c)", size = (10*0.35), family = "serif")
# ggsave(file = "save/fig_S4_20260829.tiff", width = 12, height = 18, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

rm(plot_01, plot_02, plot_03, data1, data2, data3, data4, data_P1, model)

# 09. Figure. S5 -----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_P2_Q1.rdata")
load(file = "save/1_statistical result/model_Part2_Q1_data_20250703.rdata")

data <- model_save |> filter(!is.na(Std_Coefficient)) |> mutate(damage_type = rep(c("DR_total", "DR_uproot", "DR_below", "DR_up", "DR_branch"), 105)) |> filter(str_detect(fix_name, "CV")) |> 
  select(index, species, damage_type, fix_estimate, Std_Coefficient) |> left_join(bind_rows(data_P2_Q1) |> distinct(species_CN, species_LN), by = c("species" = "species_CN"))

sp_list <- sp_list_df(sp_list = unique(data$species_LN), taxon = "plant")
set.seed(1234); tree_plant <- get_tree(sp_list = sp_list, taxon = "plant", scenario = "random_below_basal", show_grafted = FALSE)
data.frame(x = tree_plant$tip.label, y = 1:35)

data <- data |> mutate(species_LN = case_when(species_LN == "Corylus_ferox_var._thibetica" ~ "Corylus_ferox", species_LN == "Quercus_aliena_var._acutiserrata" ~ "Quercus_aliena", species_LN == "Acer_pictum_subsp._mono" ~ "Acer_pictum", species_LN == "Cornus_kousa_subsp._chinensis" ~ "Cornus_kousa", TRUE ~ species_LN)) |> 
  left_join(data.frame(species_LN = c("Carpinus turczaninowii", "Carpinus fargesiana", "Corylus ferox", "Corylus chinensis", "Betula albosinensis", "Betula luminifera", "Juglans mandshurica", "Platycarya strobilacea", "Quercus serrata", "Quercus aliena", "Quercus variabilis", "Castanea seguinii", "Castanea henryi", "Fagus engleriana", "Sorbus folgneri", "Sorbus alnifolia", "Crataegus wilsonii", "Prunus conradinae", "Prunus padus", "Populus lasiocarpa", "Populus davidiana", "Salix wallichiana", "Rhus potaninii", "Rhus chinensis", "Toxicodendron vernicifluum", "Acer davidii", "Acer pictum", "Cornus macrophylla", "Cornus kousa", "Cornus controversa", "Fraxinus chinensis", "Diospyros lotus", "Litsea ichangensis", "Lindera obtusiloba", "Pinus armandii"), 
                       ID = c(35:1), niche_width = c(0.5693252, 0.5067005, 0.4052766, 0.8033091, 0.7300763, 0.6806774, 0.7154775, 0.7134962, 0.4579634, 0.5269582, 0.5795512, 0.6247113, 0.539175, 0.4583443, 0.4789034, 0.4000346, 0.5557265, 0.5170674, 0.6401926, 0.5663812, 0.6025009, 0.6959552, 0.7151405, 0.7137618, 0.4638151, 0.7080873, 0.6139691, 0.516735, 0.3869133, 0.6158935, 0.4840941, 0.6366028, 0.5825325, 0.4570297, 0.7790024)) |> mutate(species_LN = str_replace_all(species_LN, " ", "_"))) |> arrange(ID)
range(data$Std_Coefficient)

# phylogenetic heatmap --
ggplot(tree_plant, branch.length = "none") + 
  geom_tree() + 
  geom_text(data = data |> distinct(species_LN, ID) |> mutate(species_LN = str_replace_all(species_LN, "_", " ")), aes(x = 11, y = ID, label = species_LN), size = (7*0.35), hjust = "left", fontface = "italic", family = "serif") + 
  annotate(geom = "text", x = 1.7, y = 36.4, label = "Slopes", size = (9*0.35), family = "serif", fontface = "bold") + 
  
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_TTH"), aes(x = 22, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_TTH"), aes(x = 22, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_DBH"), aes(x = 24, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_DBH"), aes(x = 24, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_NMB"), aes(x = 26, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_uproot") |> filter(index == "CV_NMB"), aes(x = 26, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 24, y = 36.8, label = "Uprooting", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 21, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 23, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 25, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_TTH"), aes(x = 29, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_TTH"), aes(x = 29, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_DBH"), aes(x = 31, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_DBH"), aes(x = 31, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_NMB"), aes(x = 33, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_below") |> filter(index == "CV_NMB"), aes(x = 33, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 31, y = 36.8, label = "Clear-bole", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 31, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 28, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 30, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 32, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_TTH"), aes(x = 36, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_TTH"), aes(x = 36, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_DBH"), aes(x = 38, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_DBH"), aes(x = 38, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_NMB"), aes(x = 40, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_up") |> filter(index == "CV_NMB"), aes(x = 40, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 38, y = 36.8, label = "Crown", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 38, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 35, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 37, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 39, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_TTH"), aes(x = 43, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_TTH"), aes(x = 43, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_DBH"), aes(x = 45, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_DBH"), aes(x = 45, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  geom_point(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_NMB"), aes(x = 47, y = ID, color = Std_Coefficient), shape = 15, size = 3.8, stroke = 0) + 
  geom_text(data = data |> filter(damage_type == "DR_branch") |> filter(index == "CV_NMB"), aes(x = 47, y = ID, label = case_when(Std_Coefficient > 0 ~ "+", Std_Coefficient < 0 ~ "−", TRUE ~ "")), size = (10*0.35), color = "#000000") + 
  annotate(geom = "text", x = 45, y = 36.8, label = "Branch", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 45, y = 36.0, label = "broken", size = (9*0.35), family = "serif") + 
  annotate(geom = "text", x = 42, y = -0.5, label = expression(CV["H    "]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 44, y = -0.5, label = expression(CV[DBH]), size = (9*0.35), family = "serif", angle = 45) + 
  annotate(geom = "text", x = 46, y = -0.5, label = expression(CV[NMB]), size = (9*0.35), family = "serif", angle = 45) + 
  
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "RdBu")[10:2])(99), limits = c(-0.2, 0.42), breaks = c(-0.2, 0, 0.2, 0.4), labels = c("−0.2", " 0", " 0.2", " 0.4")) + 
  scale_y_continuous(limits = c(-2, 38)) + 
  theme_bw(base_family = "serif") + 
  labs(x = NULL, y = NULL, color = NULL) + 
  theme(legend.key.width = unit(0.2, "cm"), 
        legend.key.height = unit(0.5, "cm"), 
        legend.position = c(0.08, 0.80), 
        legend.background = element_blank(),
        legend.text = element_text(size = 9, margin = margin(l = 2, unit = "pt")),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.margin = margin(t = -25, r = -10, b = -25, l = -15))
# ggsave(file = "save/fig_S5_20260829.tiff", width = 14, height = 14, units = "cm", dpi = 300, limitsize = FALSE)

rm(data, data_P2_Q1, model_save, sp_list, tree_plant)

# 10. Figure. S6 -----------------------------------------------------------------------------------
# data processing --
load(file = "data/data_P1.rdata")

# DR_total --
model <- glmmTMB(cbind(n_damaged, n_total - n_damaged) ~ forest_age + (1 | site), data_P1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_01 <- ggplot(data_P1, aes(x = forest_age, y = DR_total*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Total damage", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (10*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_01

# DR_uproot --
model <- glmmTMB(cbind(n_uproot, n_total - n_uproot) ~ forest_age + (1 | site), data_P1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_02 <- ggplot(data_P1, aes(x = forest_age, y = DR_uproot*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Uprooting", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (10*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "**", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_02

# DR_below --
model <- glmmTMB(cbind(n_below, n_total - n_below) ~ forest_age + (1 | site), data_P1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_03 <- ggplot(data_P1, aes(x = forest_age, y = DR_below*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Clear-bole broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "*", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_03

# DR_up --
model <- glmmTMB(cbind(n_up, n_total - n_up) ~ forest_age + (1 | site), data_P1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_04 <- ggplot(data_P1, aes(x = forest_age, y = DR_up*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Crown broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "***", size = (10*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_04

# DR_branch --
model <- glmmTMB(cbind(n_branch, n_total - n_branch) ~ forest_age + (1 | site), data_P1, family = binomial)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_05 <- ggplot(data_P1, aes(x = forest_age, y = DR_branch*100, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.5, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  labs(x = NULL, y = "Damage ratio (%)") + 
  annotate(geom = "text", x = ggpp::as_npc(0.5), y = ggpp::as_npc(0.95), label = "Branch broken", size = (9*0.35), family = "serif", hjust = "left") +
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "**", size = (10*0.35), family = "serif") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.81), xend = 2, yend = ggpp::as_npc(0.81)), color = "#4d4d4d", linewidth = 0.45) + 
  annotate(geom = "text", x = 1.5, y = ggpp::as_npc(0.82), label = "***", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_05

# TD_SR --
model <- glmmTMB(TD_SR ~ forest_age + (1 | site), data_P1)
summary(model); r2(model); emmeans::emmeans(model, ~ forest_age) |> pairs()

plot_06 <- ggplot(data_P1, aes(x = forest_age, y = TD_SR, color = forest_age, fill = forest_age)) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8, shape = 16) + 
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), color = "#000000", size = 0.8, linewidth = 0.5, alpha = 0.7, shape = 16) +
  geom_point(stat = "summary", fun = "mean", fun.args = list(mult = 1), size = 2.5, alpha = 0.7, shape = 16) + 
  scale_fill_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_color_manual(values = c("#8cc269", "#1ba784", "#1a6840")) + 
  scale_x_discrete(labels = c("ESF", "LSF", "OF")) + 
  scale_y_continuous(limits = c(2, 32), breaks = c(2, 17, 32), labels = c("2", "17", "32")) + 
  labs(x = NULL, y = "Species richness") + 
  geom_segment(aes(x = 1, y = ggpp::as_npc(0.87), xend = 3, yend = ggpp::as_npc(0.87)), color = "#4d4d4d", linewidth = 0.4) + 
  annotate(geom = "text", x = 2, y = ggpp::as_npc(0.88), label = "**", size = (10*0.35), family = "serif") + 
  theme_bw(base_family = "serif") + 
  theme(axis.text = element_text(size = 9, color = "#000000"), 
        axis.title = element_text(size = 10), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none"); plot_06

plot_grid(plot_01, plot_02, plot_03, plot_04, plot_05, plot_06, nrow = 3) + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.990), label = "(a)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.656), label = "(c)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.02), y = ggpp::as_npc(0.323), label = "(e)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.990), label = "(b)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.656), label = "(d)", size = (10*0.35), family = "serif") + 
    annotate(geom = "text", x = ggpp::as_npc(0.52), y = ggpp::as_npc(0.323), label = "(f)", size = (10*0.35), family = "serif")
ggsave(file = "save/fig_S6_20260829.tiff", width = 12, height = 18, units = "cm", dpi = 300, limitsize = FALSE, bg = "#ffffff")

rm(data_P1, model, plot_01, plot_02, plot_03, plot_04, plot_05, plot_06)
