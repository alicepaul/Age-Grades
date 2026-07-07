library(tidyverse)
library(patchwork)
library(gt)
library(gtsummary)
library(ggh4x)
library(scales)

source("age_grade_method.R")

# read in age records
df <- read.csv("Data/age_grades_25.csv")  
df$Event[df$Event == "H Mar"] <- "Half Marathon"
df$Event[df$Event == "Mar"] <- "Marathon"
df$Sex <- ifelse(df$Sex=="F", "Female", "Male")
df$Event <- factor(df$Event, levels = c("5 km", "10 km", "Half Marathon", 
                                        "Marathon"))
df <- df %>%
  arrange(Sex, Event, Age) 

# set parameters
eta <- 1000
rho <- 50
lambdas <- 10^seq(-4, 4, length.out = 50)
knots <- seq(13,80,4)

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
    cv_fit <- cv_asymmetric_spline(x = x, y = y, lambdas = lambdas, knots,
                                   rho = rho, eta = eta, degree = 3, k = 5)
    best_lambda <- cv_fit$best_lambda
    mod_fit <- fit_asymmetric_spline(x = x, y = y, lambda =  best_lambda,
                                     spline_spec = spline_spec, 
                                     rho = rho, eta = eta)
    
    # predict to find full curve and add to results
    new_std <- predict_asymmetric_spline(mod_fit, newx=c(min(x):max(x)))
    new_df <- rbind(new_df,
                    data.frame(Sex = sex,
                               Event = event,
                               Age = c(min(x):max(x)),
                               New_Std = new_std))
  }
}
df <- left_join(df, new_df)

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
ggplot(df, aes(x = Age, y = Best_Time_Sec)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(
    data = df,
    aes(y = Standard_Sec),
    #data = df_long,
    #aes(y = standard_time, color = standard_type),
    linewidth = 1
  ) +
  facet_grid(
    Event ~ Sex,
    scales = "free_y"
  ) +
  facetted_pos_scales(
    y = list(
      Event == "5 km" ~ scale_y_continuous(
        labels = \(x) {
          mins <- floor(x / 60)
          secs <- floor((x %% 60))
          sprintf("%d:%02d", mins, secs)
        }
      ),
      Event == "10 km" ~ scale_y_continuous(
        labels = \(x) {
          mins <- floor(x / 60)
          secs <- floor((x %% 60))
          sprintf("%d:%02d", mins, secs)
        }
      ),
      Event == "Half Marathon" ~ scale_y_continuous(
        labels = \(x) {
          hrs <- floor(x / 3600)
          mins <- floor((x %% 3600) / 60)
          sprintf("%d:%02d", hrs, mins)
        }
      ),
      Event == "Marathon" ~ scale_y_continuous(
        labels = \(x) {
          hrs <- floor(x / 3600)
          mins <- floor((x %% 3600) / 60)
          sprintf("%d:%02d", hrs, mins)
        }
      )
    )
  ) +
  labs(
    x = "Age",
    y = "Time",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/current_standards.png"), bg="white")

# plot existing vs proposed
ggplot(df, aes(x = Age, y = Best_Time_Sec)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(
    data = df_long,
    aes(y = standard_time, color = standard_type),
    linewidth = 1
  ) +
  facet_grid(
    Event ~ Sex,
    scales = "free_y"
  ) +
  facetted_pos_scales(
    y = list(
      Event == "5 km" ~ scale_y_continuous(
        labels = \(x) {
          mins <- floor(x / 60)
          secs <- floor((x %% 60))
          sprintf("%d:%02d", mins, secs)
        }
      ),
      Event == "10 km" ~ scale_y_continuous(
        labels = \(x) {
          mins <- floor(x / 60)
          secs <- floor((x %% 60))
          sprintf("%d:%02d", mins, secs)
        }
      ),
      Event == "Half Marathon" ~ scale_y_continuous(
        labels = \(x) {
          hrs <- floor(x / 3600)
          mins <- floor((x %% 3600) / 60)
          sprintf("%d:%02d", hrs, mins)
        }
      ),
      Event == "Marathon" ~ scale_y_continuous(
        labels = \(x) {
          hrs <- floor(x / 3600)
          mins <- floor((x %% 3600) / 60)
          sprintf("%d:%02d", hrs, mins)
        }
      )
    )
  ) +
  labs(
    x = "Age",
    y = "Time",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/new_standards.png"), bg="white")

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
    p_value = t.test(
      Female,
      Male,
      paired = TRUE,
      alternative = "less"
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Mean diff = ", round(mean_difference, 2),
      "\npaired t-test p = ", signif(p_value, 3)
    )
  )

# Annotation locations
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
    x = "Age Grade",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/stat_analysis_curr.png"), bg="white")

# statistical comparison of existing standards marathon data

df_mar <- df %>%
  filter(Event == "Marathon") %>%
  select(Age, Sex, Standard_Sec, New_Std) 

df_23 <- read.csv("Data/data_23.csv") %>%
  left_join(df_mar) %>%
  select(c("Race", "Finishers", "Age", "Sex", "Time", 
           "Standard_Sec", "New_Std")) %>%
  filter(Finishers >= 200) %>%
  mutate(Standard_AG = 100*Standard_Sec/(Time*60),
         New_AG = 100*New_Std/(Time*60))

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

winner_test <- t.test(
  paired_dat$Female,
  paired_dat$Male,
  paired = TRUE,
  alternative = "less"
)

overall_dat <- df_23 %>%
  filter(
    Sex %in% c("Female", "Male"),
    is.finite(Standard_AG)
  )

overall_test <- t.test(
  Standard_AG ~ Sex,
  data = overall_dat,
  alternative = "less"
)

winner_label <- sprintf(
  "Paired t-test\nn = %d\nMean diff = %.2f\np = %.3g",
  nrow(paired_dat),
  mean(paired_dat$Female - paired_dat$Male),
  winner_test$p.value
)

overall_label <- sprintf(
  "Welch t-test\nFemale = %.2f\nMale = %.2f\np = %.3g",
  mean(df_23$Standard_AG[df_23$Sex == "Female"], na.rm=TRUE),
  mean(df_23$Standard_AG[df_23$Sex == "Male"], na.rm =TRUE),
  overall_test$p.value
)

p1 <- ggplot(winners,
             aes(Standard_AG, fill = Sex)) +
  geom_histogram(
    position = "identity",
    alpha = 0.6,
    bins = 30
  ) +
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
    title = "Age Group Winners",
    x = "Age Grade",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(df_23,
             aes(Standard_AG, fill = Sex)) +
  geom_histogram(
    position = "identity",
    alpha = 0.6,
    bins = 40
  ) +
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
    title = "All Finishers",
    x = "Age Grade",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p1 + p2
ggsave(paste0(getwd(),"/Images/marathon_current.png"), bg="white")

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
    p_value = t.test(
      Female,
      Male,
      paired = TRUE,
      alternative = "less"
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Mean diff = ", round(mean_difference, 2),
      "\npaired t-test p = ", signif(p_value, 3)
    )
  )

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
    x = "Age Grade (Proposed)",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(paste0(getwd(),"/Images/stat_analysis_new.png"), bg="white")

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

winner_test <- t.test(
  paired_dat$Female,
  paired_dat$Male,
  paired = TRUE,
  alternative = "less"
)

overall_dat <- df_23 %>%
  filter(
    Sex %in% c("Female", "Male"),
    is.finite(New_AG)
  )

overall_test <- t.test(
  New_AG ~ Sex,
  data = overall_dat,
  alternative = "less"
)

winner_label <- sprintf(
  "Paired t-test\nn = %d\nMean diff = %.2f\np = %.3g",
  nrow(paired_dat),
  mean(paired_dat$Female - paired_dat$Male),
  winner_test$p.value
)

overall_label <- sprintf(
  "Welch t-test\nFemale = %.2f\nMale = %.2f\np = %.3g",
  mean(df_23$New_AG[df_23$Sex == "Female"], na.rm=TRUE),
  mean(df_23$New_AG[df_23$Sex == "Male"], na.rm =TRUE),
  overall_test$p.value
)

p1 <- ggplot(winners,
             aes(New_AG, fill = Sex)) +
  geom_histogram(
    position = "identity",
    alpha = 0.6,
    bins = 30
  ) +
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
    title = "Age Group Winners",
    x = "Age Grade",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(df_23,
             aes(New_AG, fill = Sex)) +
  geom_histogram(
    position = "identity",
    alpha = 0.6,
    bins = 40
  ) +
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
    title = "All Finishers",
    x = "Age Grade",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p1 + p2
ggsave(paste0(getwd(),"/Images/marathon_new.png"), bg="white")
