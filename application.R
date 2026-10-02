library(tidyverse)
library(patchwork)
library(gt)
library(gtsummary)
library(ggh4x)
library(scales)
library(clusrank)

# create directories
out <- "Output"
dir.create("Output", recursive = TRUE, showWarnings = FALSE)
dir.create("Images", showWarnings = FALSE)

source("age_grade_method.R")

# read in age records
df <- read.csv("Data/age_grades_25.csv")

# data checks
stopifnot(
  all(c("Event", "Sex", "Age", "Best_Time_Sec", "Standard_Sec") %in%
        names(df)),
  all(complete.cases(df)),
  all(df$Sex %in% c("F", "M")),
  all(df$Event %in% c("5 km", "10 km", "H Mar", "Mar")),
  all(is.finite(df$Age)), all(df$Age == floor(df$Age)),
  all(df$Age >= 13 & df$Age <= 80),
  all(is.finite(df$Best_Time_Sec)), all(df$Best_Time_Sec > 0),
  all(is.finite(df$Standard_Sec)), all(df$Standard_Sec > 0),
  !anyDuplicated(df[c("Event", "Sex", "Age")])
)

# code events
df$Event[df$Event == "H Mar"] <- "Half Marathon"
df$Event[df$Event == "Mar"] <- "Marathon"
df$Sex <- ifelse(df$Sex=="F", "Female", "Male")
df$Event <- factor(df$Event, levels = c("5 km", "10 km", "Half Marathon", 
                                        "Marathon"))
df <- df %>%
  arrange(Sex, Event, Age) 

# report missing ages without filling in unobserved records.
df %>%
  group_by(Event, Sex) %>%
  summarise(
    n = n(), min_age = min(Age), max_age = max(Age),
    missing_ages = paste(setdiff(13:80, Age), collapse = ";"),
    .groups = "drop"
  )

# set parameters
eta <- 1000
rho <- 30
convex_after <- 35
lambdas <- 10^seq(-4, 4, length.out = 50)
knots <- seq(13,80,4)

# store results on influence
influence <- list()
fit_checks <- list()

# run method for all sex-event combinations and store 
new_df <- data.frame(Sex = character(), 
                     Event = character(),
                     Age = numeric(),
                     New_Std = numeric())

for (sex in c("Female", "Male")){
  for (event in c("5 km", "10 km", "Half Marathon", "Marathon")){
    # get subset of data
    df_filter <- df %>%
      filter(Event == event, Sex == sex)
    x <- df_filter$Age
    y <- df_filter$Best_Time_Sec
    
    # find optimal smoothing parameter and get final model
    spline_spec <- make_spline_spec(x, knots, degree = 3)
    cv_fit <- cv_asymmetric_spline(x = x, y = y, lambdas = lambdas,
                                   knots = knots, rho = rho, eta = eta,
                                   degree = 3, k = 5,
                                   convex_after = convex_after)
    best_lambda <- cv_fit$best_lambda
    mod_fit <- fit_asymmetric_spline(x = x, y = y, lambda =  best_lambda,
                                     spline_spec = spline_spec, 
                                     rho = rho, eta = eta,
                                     convex_after = convex_after)
    
    # predict to find full curve and add to results
    new_std <- predict_asymmetric_spline(mod_fit, newx=c(min(x):max(x)))
   
    # use leave-one-out to look at influence for each point
    loo <- sapply(seq_along(x), function(i) {
      deleted_fit <- fit_asymmetric_spline(
        x[-i], y[-i], lambda = best_lambda, spline_spec = spline_spec,
        rho = rho, eta = eta, convex_after = convex_after
      )
      predict_asymmetric_spline(deleted_fit, min(x):max(x))
    })
    influence[[length(influence) + 1]] <- data.frame(
      Sex = sex, Event = event, Age = min(x):max(x),
      New_Std = new_std, lower = apply(loo, 1, min),
      upper = apply(loo, 1, max)
    )
    fit_checks[[length(fit_checks) + 1]] <- data.frame(
      Sex = sex, Event = event, lambda = best_lambda,
      iterations = mod_fit$iterations
    )
    

    new_df <- rbind(new_df,
                    data.frame(Sex = sex,
                               Event = event,
                               Age = c(min(x):max(x)),
                               New_Std = new_std))
  }
}

# add standards
df <- left_join(df, new_df, by = c("Sex", "Event", "Age"))
df$Event <- factor(df$Event,
  levels = c("5 km", "10 km", "Half Marathon", "Marathon"))

# combine the influence results and keep events in the same order
influence <- bind_rows(influence) %>%
  mutate(Event = factor(Event, levels = levels(df$Event)))

# plot influence
ggplot(influence, aes(x = Age)) +
  geom_ribbon(aes(ymin = lower/60, ymax = upper/60), fill = "grey75") +
  geom_line(aes(y = New_Std/60), color = "red") +
  facet_grid(Event ~ Sex, scales = "free_y") +
  labs(x = "Age (years)", y = "Standard time (minutes)",
       caption = "Delete-one envelope at fixed lambda; not a 95% CI") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("Images/delete_one_envelopes.png", width = 8, height = 9,
       bg = "white")

# sensitivity to the two penalties, one at a time
sensitivity_settings <- data.frame(
  Penalty = rep(c("rho", "eta"), each = 3),
  Level = rep(c("Low", "Current", "High"), 2),
  rho = c(1, rho, 1000, rep(rho, 3)),
  eta = c(rep(eta, 3), 0, eta, 1000000)
)
sensitivity <- list()
sensitivity_checks <- list()

for (sex in c("Female", "Male")) {
  for (event in levels(df$Event)) {
    df_filter <- df %>%
      filter(Event == event, Sex == sex)
    x <- df_filter$Age
    y <- df_filter$Best_Time_Sec
    spline_spec <- make_spline_spec(x, knots, degree = 3)
    best_lambda <- bind_rows(fit_checks) %>%
      filter(Event == event, Sex == sex) %>%
      pull(lambda)

    for (j in seq_len(nrow(sensitivity_settings))) {
      setting <- sensitivity_settings[j, ]
      sensitivity_fit <- fit_asymmetric_spline(
        x, y, lambda = best_lambda, rho = setting$rho,
        eta = setting$eta, spline_spec = spline_spec,
        convex_after = convex_after, max_iter = 2000
      )
      if (sensitivity_fit$iterations >= 2000) {
        stop("Sensitivity fit reached iteration limit: ",
             sex, " / ", event, " / ", setting$Penalty, " / ",
             setting$Level)
      }
      standard <- predict_asymmetric_spline(sensitivity_fit, x)
      sensitivity[[length(sensitivity) + 1]] <- data.frame(
        Sex = sex, Event = event, Age = x,
        Penalty = setting$Penalty, Level = setting$Level,
        rho = setting$rho, eta = setting$eta,
        Standard = standard, Baseline = df_filter$New_Std,
        Change = 100*(standard/df_filter$New_Std - 1)
      )
      sensitivity_checks[[length(sensitivity_checks) + 1]] <-
        data.frame(Sex = sex, Event = event,
                   Penalty = setting$Penalty, Level = setting$Level,
                   rho = setting$rho, eta = setting$eta,
                   lambda = best_lambda,
                   iterations = sensitivity_fit$iterations)
    }
  }
}
sensitivity <- bind_rows(sensitivity) %>%
  mutate(Event = factor(Event, levels = levels(df$Event)),
         Penalty = factor(Penalty, levels = c("rho", "eta")),
         Level = factor(Level, levels = c("Low", "Current", "High")))
write.csv(sensitivity, file.path(out, "penalty_sensitivity.csv"),
          row.names = FALSE)
write.csv(bind_rows(sensitivity_checks),
          file.path(out, "penalty_sensitivity_checks.csv"),
          row.names = FALSE)

# relative changes 
sensitivity_caption <- paste0(
  "rho: 1, ", rho, ", 1000 (eta = ", eta, "); ",
  "eta: 0, ", eta, ", 1000000 (rho = ", rho, ").\n",
  "Lambda fixed at the current fit's CV choice in each sex-event group."
)
ggplot(sensitivity, aes(Age, Change, color = Level,
                       linetype = Level)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(linewidth = 0.7) +
  facet_grid(Event ~ Penalty + Sex, scales = "free_y") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  labs(x = "Age (years)", y = "Change in standard time (%)",
       color = "Penalty setting", linetype = "Penalty setting",
       caption = sensitivity_caption) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("Images/penalty_sensitivity.png", width = 13, height = 10,
       bg = "white")

# show the fitted standards 
ggplot(filter(sensitivity, Event == "5 km"),
       aes(Age, Standard/60, color = Level,
                       linetype = Level)) +
  geom_point(data = filter(df, Event == "5 km"),
             aes(Age, Best_Time_Sec/60),
             inherit.aes = FALSE, size = 0.5, alpha = 0.5) +
  geom_line(linewidth = 0.7) +
  facet_grid(Penalty ~ Sex) +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  labs(title = "5 km: sensitivity to penalty settings",
       x = "Age (years)", y = "Standard time (minutes)",
       color = "Penalty setting", linetype = "Penalty setting",
       caption = sensitivity_caption) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("Images/penalty_sensitivity_standards.png", width = 10,
       height = 8, bg = "white")

# reformat data with both standards
df_long <- df %>%
  pivot_longer(
    cols = c(Standard_Sec, New_Std),
    names_to = "standard_type",
    values_to = "standard_time"
  ) %>%
  mutate(
    standard_type = recode(
      standard_type,
      Standard_Sec = "Existing Standard",
      New_Std = "Proposed Standard"
    )
  )

# check events ordered consistently
df$Event <- factor(df$Event,
                   levels = c("5 km", "10 km", "Half Marathon", "Marathon"))
df_long$Event <- factor(df_long$Event,
                        levels = levels(df$Event))

# plot existing 
ggplot(df, aes(Age, Standard_Sec/60)) +
  geom_line(color = "red") +
  geom_point(data = df, aes(Age, Best_Time_Sec/60), inherit.aes = FALSE,
             size = 0.6) +
  facet_grid(Event ~ Sex, scales = "free_y") +
  labs(y = "Time (minutes)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"))
ggsave(paste0(getwd(),"/Images/current_standards.png"),
       width = 10, height = 7, bg = "white")

# plot existing vs proposed
ggplot(df, aes(x = Age, y = Best_Time_Sec/60)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(
    data = df_long,
    aes(y = standard_time/60, color = standard_type),
    linewidth = 1
  ) +
  facet_grid(
    Event ~ Sex,
    scales = "free_y"
  ) +
  labs(
    x = "Age (years)",
    y = "Time (minutes)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/new_standards.png"),
       width = 10, height = 7, bg = "white")

# statistical comparison for existing standards on single-age best

df$Standard_AG <- 100 * (df$Standard_Sec/df$Best_Time_Sec)
df$New_AG <- 100*df$New_Std/df$Best_Time_Sec

ttest_results <- df %>%
  select(Event, Sex, Age, Standard_AG) %>%
  pivot_wider(
    names_from = Sex,
    values_from = Standard_AG
  ) %>%
  group_by(Event) %>%
  summarise(
    n_pairs = sum(!is.na(Female - Male)),
    mean_female = mean(Female, na.rm = TRUE),
    mean_male = mean(Male, na.rm = TRUE),
    mean_difference = mean(Female - Male, na.rm = TRUE),
    ci_lower = t.test(Female, Male, paired = TRUE)$conf.int[1],
    ci_upper = t.test(Female, Male, paired = TRUE)$conf.int[2],
    p_value = t.test(
      Female,
      Male,
      paired = TRUE,
      alternative = "two.sided"
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Mean diff = ", round(mean_difference, 2),
      " pp\n95% CI [", round(ci_lower, 2), ", ",
      round(ci_upper, 2), "]",
      "\npaired t-test p (two-sided) ",
      format.pval(p_value, digits = 3, eps = 0.001)
    )
  )

write.csv(ttest_results, file.path(out, "records_existing.csv"),
          row.names = FALSE)

# annotations
plot_labels <- df %>%
  group_by(Event) %>%
  summarise(
    x_pos = max(Standard_AG, na.rm = TRUE),
    y_pos = Inf,
    .groups = "drop"
  ) %>%
  left_join(ttest_results, by = "Event")

# 2x2 histogram grid
ggplot(df, aes(x = Standard_AG, fill = Sex)) +
  geom_histogram(
    alpha = 0.6,
    bins = 30,
    position = "identity"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  geom_text(
    data = plot_labels,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.2,
    size = 3
  ) +
  facet_wrap(~ Event, nrow = 2, ncol = 2, scales = "free_y") +
  labs(
    x = "Age grade (%)",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/stat_analysis_curr.png"),
       width = 10, height = 7, bg = "white")

# statistical comparison of existing standards marathon data

df_mar <- df %>%
  filter(Event == "Marathon") %>%
  select(Age, Sex, Standard_Sec, New_Std) 

# marathon data
df_23 <- read.csv("Data/data_23.csv")
df_23 <- df_23 %>%
  filter(is.finite(Finishers), Finishers >= 200) %>%
  filter(Sex %in% c("Female", "Male"),
         is.finite(Age), Age == floor(Age),
         is.finite(Time), Time > 0, !is.na(Race), trimws(Race) != "")

df_23 <- df_23 %>%
  left_join(df_mar, by = c("Age", "Sex")) %>%
  filter(is.finite(Standard_Sec), is.finite(New_Std)) %>%
  mutate(Standard_AG = 100*Standard_Sec/(Time*60),
         New_AG = 100*New_Std/(Time*60))

# participation trends
participation <- df_23 %>%
  count(Race, Sex) %>%
  group_by(Sex) %>%
  mutate(share = n/sum(n)) %>%
  ungroup()
participation_wide <- participation %>%
  select(Race, Sex, share) %>%
  pivot_wider(names_from = Sex, values_from = share, values_fill = 0)
write.csv(participation, file.path(out, "participation.csv"),
          row.names = FALSE)
p_race <- ggplot(participation_wide, aes(x = Male, y = Female)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  coord_equal() +
  labs(x = "Share of male finishers", y = "Share of female finishers",
       title = "Participation across races") +
  theme_minimal()
p_age <- ggplot(df_23, aes(x = Age, fill = Sex)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5,
                 boundary = 0, position = "identity", alpha = 0.6) +
  labs(x = "Age (years)", y = "Density within sex",
       title = "Participation across ages") +
  theme_minimal() +
  theme(legend.position = "bottom")
p_race + p_age
ggsave("Images/marathon_participation.png", width = 10, height = 5,
       bg = "white")

# age group winners - signed ranks use Rosner-Glynn-Lee (RGL)
# all-finisher rank sums use Datta-Satten (DS) to retain 
# unequal-size races.
cluster_rank_sum <- function(value, race, sex) {
  stopifnot(length(value) == length(race), length(value) == length(sex),
            all(is.finite(value)), !anyNA(race),
            all(sex %in% c("Female", "Male")))
  id <- as.integer(factor(race))
  size <- tabulate(id)
  n_races <- length(size)
  stopifnot(n_races >= 2)
  female <- as.numeric(sex == "Female")
  proportion <- as.vector(rowsum(female, id))/size
  stopifnot(any(female == 1), any(female == 0))

  # sum the race-specific midpoint CDFs at each observed value.
  order_value <- order(value)
  ties <- rle(value[order_value])$lengths
  tie_id <- rep(seq_along(ties), ties)
  mass <- as.vector(rowsum(1/size[id[order_value]], tie_id))
  total_cdf <- numeric(length(value))
  total_cdf[order_value] <- rep(cumsum(mass) - mass/2, ties)
  within_cdf <- ave(value, id, FUN = function(x) {
    (rank(x, ties.method = "average") - 0.5)/length(x)
  })
  other_cdf <- total_cdf - within_cdf
  pooled_cdf <- (rank(value, ties.method = "average") - 0.5)/
    length(value)

  # Cluster contributions to the DS statistic and variance estimate.
  score <- sum(female/size[id]*(1 + other_cdf))/(n_races + 1)
  expected <- sum(proportion)/2
  contribution <- as.vector(rowsum(
    ((n_races - 1)*female -
       (sum(proportion) - proportion[id]))*pooled_cdf/
      (size[id]*(n_races + 1)), id
  ))
  expected_contribution <- n_races/(2*(n_races + 1))*
    (proportion - mean(proportion))
  variance <- sum((contribution - expected_contribution)^2)
  if (!is.finite(variance) || variance <= 0) {
    stop("Clustered rank-sum variance is zero or undefined")
  }
  z <- (score - expected)/sqrt(variance)
  structure(list(
    statistic = c(Z = z), p.value = 2*pnorm(-abs(z)),
    alternative = "two.sided", nobs = length(value), nclus = n_races,
    method = "Clustered rank-sum test (Datta-Satten)",
    data.name = "value by sex, clustered by race"
  ), class = "htest")
}

# bootstrap the descriptive mean difference by sampling whole races.
race_mean_ci <- function(value, race, sex = NULL,
                         B = 1999, seed = 20261002) {
  stopifnot(length(value) == length(race), all(is.finite(value)),
            !anyNA(race))
  ids <- as.integer(factor(race))
  if (is.null(sex)) {
    totals <- rowsum(cbind(value, rep(0, length(value))), ids)
    counts <- rowsum(cbind(rep(1, length(value)),
                          rep(1, length(value))), ids)
  } else {
    stopifnot(length(sex) == length(value),
              all(sex %in% c("Female", "Male")))
    female <- as.numeric(sex == "Female")
    male <- as.numeric(sex == "Male")
    totals <- rowsum(cbind(value*female, value*male), ids)
    counts <- rowsum(cbind(female, male), ids)
  }
  # Each sampled race contributes all its observations, including repeats.
  estimate <- function(index) {
    total <- colSums(totals[index, , drop = FALSE])
    count <- colSums(counts[index, , drop = FALSE])
    if (any(count == 0)) return(NA_real_)
    total[1]/count[1] - total[2]/count[2]
  }
  n_races <- nrow(totals)
  stopifnot(n_races >= 2)
  set.seed(seed)
  draws <- replicate(B, estimate(sample.int(n_races, n_races,
                                            replace = TRUE)))
  if (any(!is.finite(draws))) {
    stop("A bootstrap sample lacks one sex; inspect race coverage")
  }
  data.frame(
    mean_difference = estimate(seq_len(n_races)),
    ci_lower = unname(quantile(draws, 0.025)),
    ci_upper = unname(quantile(draws, 0.975)),
    ci_method = "Whole-race bootstrap of mean, percentile 95%",
    n = length(value), n_races = n_races, bootstrap_draws = B
  )
}

# report results
rank_summary <- function(test, interval) {
  stopifnot(test$alternative == "two.sided")
  z <- unname(test$statistic)
  stopifnot(length(z) == 1, is.finite(z))
  cbind(interval, data.frame(
    statistic = z, p_value = test$p.value,
    log_p_value = log(2) + pnorm(-abs(z), log.p = TRUE),
    p_display = format.pval(test$p.value, digits = 3, eps = 0.001),
    test = test$method, alternative = "two.sided"
  ))
}

winners <- df_23 %>%
  mutate(
    AgeGroup = paste0(
      5 * floor(Age / 5), "-",
      5 * floor(Age / 5) + 4
    )
  ) %>%
  group_by(Race, AgeGroup, Sex) %>%
  slice_min(Time, n = 1, with_ties = FALSE) %>%
  ungroup()

paired_dat <- winners %>%
  select(Race, AgeGroup, Sex, Standard_AG) %>%
  pivot_wider(
    names_from = Sex,
    values_from = Standard_AG
  ) %>%
  drop_na(Female, Male)

# retain paired age-group differences and cluster by race.
paired_dat <- paired_dat %>%
  mutate(Difference = Female - Male)
winner_test <- clusWilcox.test(
  paired_dat$Difference, cluster = as.integer(factor(paired_dat$Race)),
  paired = TRUE, method = "rgl", exact = FALSE,
  alternative = "two.sided"
)
winner_result <- rank_summary(winner_test,
  race_mean_ci(paired_dat$Difference, paired_dat$Race))

# all finishers are unpaired across sexes: use clustered rank sums.
overall_dat <- df_23
overall_test <- cluster_rank_sum(
  overall_dat$Standard_AG, overall_dat$Race, overall_dat$Sex
)
overall_result <- rank_summary(overall_test,
  race_mean_ci(overall_dat$Standard_AG, overall_dat$Race,
               overall_dat$Sex))
write.csv(bind_rows(
  mutate(winner_result, Analysis = "Matched winners"),
  mutate(overall_result, Analysis = "All finishers")
), file.path(out, "marathon_existing.csv"), row.names = FALSE)

winner_label <- sprintf(
  paste0("Mean gap = %.2f pp\nBootstrap 95%% CI [%.2f, %.2f]",
         "\nClustered signed-rank p %s"),
  winner_result$mean_difference,
  winner_result$ci_lower, winner_result$ci_upper,
  winner_result$p_display
)
overall_label <- sprintf(
  paste0("Mean gap = %.2f pp\nBootstrap 95%% CI [%.2f, %.2f]",
         "\nClustered rank-sum p %s"),
  overall_result$mean_difference,
  overall_result$ci_lower, overall_result$ci_upper,
  overall_result$p_display
)

# Plot the same matched winners that enter the paired rank test.
winners <- winners %>%
  semi_join(paired_dat, by = c("Race", "AgeGroup"))

p1 <- ggplot(winners,
             aes(Standard_AG, fill = Sex)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 30
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = winner_label,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5
  ) +
  labs(
    title = "Matched age-group winners",
    subtitle = "Two-sided rank test; CI describes the mean gap",
    x = "Age grade (%)",
    y = "Density within sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(df_23,
             aes(Standard_AG, fill = Sex)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 40
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = overall_label,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5
  ) +
  labs(
    title = "All finishers",
    subtitle = "Two-sided rank test; CI describes the mean gap",
    x = "Age grade (%)",
    y = "Density within sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p1 + p2
ggsave(paste0(getwd(),"/Images/marathon_current.png"),
       width = 10, height = 7, bg = "white")

# statistical comparison of proposed standards on single age bests

ttest_results <- df %>%
  select(Event, Sex, Age, New_AG) %>%
  pivot_wider(
    names_from = Sex,
    values_from = New_AG
  ) %>%
  group_by(Event) %>%
  summarise(
    n_pairs = sum(!is.na(Female - Male)),
    mean_female = mean(Female, na.rm = TRUE),
    mean_male = mean(Male, na.rm = TRUE),
    mean_difference = mean(Female - Male, na.rm = TRUE),
    ci_lower = t.test(Female, Male, paired = TRUE)$conf.int[1],
    ci_upper = t.test(Female, Male, paired = TRUE)$conf.int[2],
    p_value = t.test(
      Female,
      Male,
      paired = TRUE,
      alternative = "two.sided"
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Mean diff = ", round(mean_difference, 2),
      " pp\n95% CI [", round(ci_lower, 2), ", ",
      round(ci_upper, 2), "]",
      "\npaired t-test p (two-sided) ",
      format.pval(p_value, digits = 3, eps = 0.001)
    )
  )

write.csv(ttest_results, file.path(out, "records_proposed.csv"),
          row.names = FALSE)

# Annotation locations
plot_labels <- df %>%
  group_by(Event) %>%
  summarise(
    x_pos = max(New_AG, na.rm = TRUE),
    y_pos = Inf,
    .groups = "drop"
  ) %>%
  left_join(ttest_results, by = "Event")

# 2x2 histogram grid
ggplot(df, aes(x = New_AG, fill = Sex)) +
  geom_histogram(
    alpha = 0.6,
    bins = 30,
    position = "identity"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  geom_text(
    data = plot_labels,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.2,
    size = 3
  ) +
  facet_wrap(~ Event, nrow = 2, ncol = 2, scales = "free_y") +
  labs(
    x = "Proposed age grade (%)",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/stat_analysis_new.png"),
       width = 10, height = 7, bg = "white")

# statistical comparison of proposed standards on marathon data

winners <- df_23 %>%
  mutate(
    AgeGroup = paste0(
      5 * floor(Age / 5), "-",
      5 * floor(Age / 5) + 4
    )
  ) %>%
  group_by(Race, AgeGroup, Sex) %>%
  slice_min(Time, n = 1, with_ties = FALSE) %>%
  ungroup()

paired_dat <- winners %>%
  select(Race, AgeGroup, Sex, New_AG) %>%
  pivot_wider(
    names_from = Sex,
    values_from = New_AG
  ) %>%
  drop_na(Female, Male)

# Retain paired age-group differences and cluster by race.
paired_dat <- paired_dat %>%
  mutate(Difference = Female - Male)
winner_test <- clusWilcox.test(
  paired_dat$Difference, cluster = as.integer(factor(paired_dat$Race)),
  paired = TRUE, method = "rgl", exact = FALSE,
  alternative = "two.sided"
)
winner_result <- rank_summary(winner_test,
  race_mean_ci(paired_dat$Difference, paired_dat$Race))

# All finishers are unpaired across sexes: use clustered rank sums.
overall_dat <- df_23
overall_test <- cluster_rank_sum(
  overall_dat$New_AG, overall_dat$Race, overall_dat$Sex
)
overall_result <- rank_summary(overall_test,
  race_mean_ci(overall_dat$New_AG, overall_dat$Race,
               overall_dat$Sex))
write.csv(bind_rows(
  mutate(winner_result, Analysis = "Matched winners"),
  mutate(overall_result, Analysis = "All finishers")
), file.path(out, "marathon_proposed.csv"), row.names = FALSE)

winner_label <- sprintf(
  paste0("Mean gap = %.2f pp\nBootstrap 95%% CI [%.2f, %.2f]",
         "\nClustered signed-rank p %s"),
  winner_result$mean_difference,
  winner_result$ci_lower, winner_result$ci_upper,
  winner_result$p_display
)
overall_label <- sprintf(
  paste0("Mean gap = %.2f pp\nBootstrap 95%% CI [%.2f, %.2f]",
         "\nClustered rank-sum p %s"),
  overall_result$mean_difference,
  overall_result$ci_lower, overall_result$ci_upper,
  overall_result$p_display
)

# Plot the same matched winners that enter the paired rank test.
winners <- winners %>%
  semi_join(paired_dat, by = c("Race", "AgeGroup"))

p1 <- ggplot(winners,
             aes(New_AG, fill = Sex)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 30
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = winner_label,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5
  ) +
  labs(
    title = "Matched age-group winners",
    subtitle = "Two-sided rank test; CI describes the mean gap",
    x = "Age grade (%)",
    y = "Density within sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(df_23,
             aes(New_AG, fill = Sex)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 40
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = overall_label,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5
  ) +
  labs(
    title = "All finishers",
    subtitle = "Two-sided rank test; CI describes the mean gap",
    x = "Age grade (%)",
    y = "Density within sex"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p1 + p2
ggsave(paste0(getwd(),"/Images/marathon_new.png"),
       width = 10, height = 7, bg = "white")



# compare standards

df_23 <- df_23 %>%
  mutate(Grade_Change = New_AG - Standard_AG)
winners <- winners %>%
  mutate(Grade_Change = New_AG - Standard_AG)
paired_change <- winners %>%
  select(Race, AgeGroup, Sex, Grade_Change) %>%
  pivot_wider(names_from = Sex, values_from = Grade_Change) %>%
  drop_na(Female, Male) %>%
  mutate(Difference = Female - Male)

change_winner_test <- clusWilcox.test(
  paired_change$Difference,
  cluster = as.integer(factor(paired_change$Race)),
  paired = TRUE, method = "rgl", exact = FALSE,
  alternative = "two.sided"
)
change_winner_result <- rank_summary(change_winner_test,
  race_mean_ci(paired_change$Difference, paired_change$Race))

change_overall_test <- cluster_rank_sum(
  df_23$Grade_Change, df_23$Race, df_23$Sex
)
change_overall_result <- rank_summary(change_overall_test,
  race_mean_ci(df_23$Grade_Change, df_23$Race, df_23$Sex))
method_comparison <- bind_rows(
  mutate(change_winner_result, Analysis = "Matched winners"),
  mutate(change_overall_result, Analysis = "All finishers")
) %>%
  rename(change_in_mean_gap = mean_difference)

# attach the two baseline gaps for an interpretable method comparison.
existing_results <- read.csv(file.path(out, "marathon_existing.csv"))
proposed_results <- read.csv(file.path(out, "marathon_proposed.csv"))
method_comparison <- method_comparison %>%
  left_join(existing_results %>%
    select(Analysis, existing_gap = mean_difference), by = "Analysis") %>%
  left_join(proposed_results %>%
    select(Analysis, proposed_gap = mean_difference), by = "Analysis")
stopifnot(all(abs(method_comparison$change_in_mean_gap -
  (method_comparison$proposed_gap - method_comparison$existing_gap)) <
  1e-10))
write.csv(method_comparison, file.path(out, "method_comparison.csv"),
          row.names = FALSE)
