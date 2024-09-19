# libraries
library(tidyverse)

metabolite_mdd <- read.table(
  "./Output/MDF.tsv",
  sep = "\t",
  header = TRUE
)
metabolite_mdd[, 66]

# running the models
length <- 60
results <- data.frame(
  metabolite = 1:length,
  mdd_beta = 1:length,
  mdd_SE = 1:length,
  mdd_P = 1:length,
  n = 1:length
)

for (i in 6:65) {
  metabolite <- colnames(metabolite_mdd)[i]
  print(i)
  print(metabolite)
  mod <- lm(
    metabolite_mdd[
      ,
      metabolite
    ] ~
      as.factor(MDD_status) + as.factor(batch_indicator),
    data = metabolite_mdd
  )
  results[
    i - 5,
    1
  ] <- metabolite
  results[
    i - 5,
    2
  ] <- coef(summary(mod))[2, 1]
  results[
    i - 5,
    3
  ] <- coef(summary(mod))[2, 2]
  results[
    i - 5,
    4
  ] <- coef(summary(mod))[2, 4]
  results[
    i - 5,
    5
  ] <- nobs(mod)
}

results$FDR_P <- p.adjust(
  results[
    ,
    "mdd_P"
  ],
  method = "BH"
)
results <- results[order(results$FDR_P), ]

write.csv(results,
  "./Output/Results_MDD.csv",
  row.names = F
) # change to appropriate file
