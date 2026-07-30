# Load data in R:
# Clear environment
rm(list = ls(all = FALSE))
library(tidyverse)

spatial_dir <- "output/spatial_profiles"

spatial_files <- list.files(
  path = spatial_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


spatial_profiles = read.csv(file="output/spatial_profiles.csv")

spatial_profiles$times = as.factor(round(spatial_profiles$time_days))

#Soluble C
ggplot(spatial_profiles, aes(x=radius_mm, y=(C), color=times))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)


ggplot(spatial_profiles, aes(x=radius_mm, y=log10(B), color=times))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)

ggplot(spatial_profiles, aes(x=radius_mm, y=log10(F_i), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)

p_glue =ggplot(spatial_profiles, aes(x=radius_mm, y=log10(F_i+0.5*E), color=factor(times)))+
  geom_line()+
  geom_hline(yintercept = log10(0.019), linetype = "dashed")+
  facet_grid(factor(diam_mm)~.)+
  labs(
    title = "Glue Distribution",
    subtitle = "F_i + E/2",
    x = "r (mm)", y = "Glue",
    color = "Days"
  ) +
  theme_minimal(base_size = 12)

print(p_glue)

stop()

ggplot(spatial_profiles, aes(x=radius_mm, y=log10(F_m), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)


ggplot(spatial_profiles, aes(x=radius_mm, y=log10(F_n), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)


ggplot(spatial_profiles, aes(x=radius_mm, y=log10(E), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)


ggplot(spatial_profiles, aes(x=radius_mm, y=log10(M), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)

ggplot(spatial_profiles, aes(x=radius_mm, y=(O), color=factor(time_days)))+
  geom_line()+
  facet_grid(factor(diam_mm)~.)



summary_data = read.csv(file="output/summary.csv")

ggplot(summary_data, aes(x=time_days, y=r_agg_mm, color=factor(r_0_mm)))+
  geom_point()+geom_line()


ggplot(summary_data, aes(x=time_days, y=POM_permille, color=factor(r_0_mm)))+
  geom_point()+geom_line()


summary_data[summary_data$diam_mm==1.25,]

plot(CO2_permille~time_days,summary_data[summary_data$diam_mm==1.25,])
abline(v=10)
abline(v=29)

ggplot(summary_data, aes(x = time_days, y = CO2_permille, color=factor(diam_mm))) +
  geom_line() +
  labs(x = "Time (days)", y = "CO@_permille",
       title = "Bacterial distribution at day 365")




