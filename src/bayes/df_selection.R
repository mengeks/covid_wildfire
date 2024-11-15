library(abind)
library(tidyr)
library(splines)
library(lme4)
library(pscl)
library(rjags)
library(ggpubr)
library(here)

remove(list = ls())

### Data Loading
source(here("src/Utilities.R"))
source(here("src/bayes/model.R"))
source(here("src/bayes/bayes_fun.R"))

dff <- load.data()
dff$FIPS <- as.numeric(as.character(dff$FIPS))
dff$pm_counter <- dff$pm25
dff$pm_counter[dff$pm_wildfire != 0] <- dff$pm25_history[dff$pm_wildfire != 0]
dff$FIPS <- as.numeric(as.character(dff$FIPS))

### Optimal DF
models = c("Constrained", "Unconstrained")
causes = c("cases", "deaths")
df.combo = list(c(3,1), c(4,1), c(5,2), c(6,2), c(7,2),c(8,2), c(9,3), c(10,3)) 

### execute modelling 
result <- c()
for (model in models) {
  for (cause in causes) {
    for (idf.combo in df.combo) {
      # pm_model takes in the whole dataset and specified combo, returns a point estimate and the confidence interval UB and LB
      glmer_result <- pm_model(dff, lags=0:28, df.date=idf.combo[1], df.tmmx=idf.combo[2], df.rmax=idf.combo[2], cause = cause, model = model)
      out <- c(coefs = glmer_result, cause = cause, combo = paste("(", paste(idf.combo[1], idf.combo[2], idf.combo[2], sep = ","), ")", sep = ""), model = model)
      result <- rbind(result, out)
      
    }
  }
}

result = data.frame(result)
write.csv(result, 
          file = here("output/bayes/df_selection.csv"), 
          row.names = FALSE)

### Sensitivity of DF
glmer_result <- read.csv(here("output/bayes/df_selection.csv"))
ford <- glmer_result$combo[1:(nrow(glmer_result)/4)]
glmer_result$combo <- factor(glmer_result$combo, levels = ford)
glmer_result$cause <- stringr::str_to_title(glmer_result$cause)

df_select <- 
  ggplot(aes(x = combo, y = coefs.pm, group = cause, colour = cause), 
         data = subset(glmer_result, model == "Constrained")) +
  geom_errorbar(aes(ymin = coefs.pm.low, ymax = coefs.pm.high), width = .2, size = 0.5, position = position_dodge(0.25)) +
  geom_line(position = position_dodge(0.25), size = 1) +
  geom_point(position = position_dodge(0.25), size = 2) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(-10, 15) +
  theme_bw() +
  labs(x = "DF Combination", y = "Cumulative Effect", color = "Event")

pdf("~/Dropbox/Projects/Wildfires/Output/bayes/df_selection.pdf", width = 5, height = 5)  
df_select
dev.off()
