# autocorrelation in disease presence

compd.data <- read.csv("./Data/compiled_data_summary_130919.csv", header = T)
compd.data <- subset(compd.data, is.na(CLASS) == FALSE & is.na(year) == F)

disease.set <- c("ADULTS", "ALL_AGE", "ALL_AGE_SUSP", "LAMBS", "LAMBS_SUSP")
  
compd.data$lag1match <- compd.data$lag2match <- rep(NA, dim(compd.data)[1])
compd.data$lag3match <- compd.data$lag4match <- rep(NA, dim(compd.data)[1])
compd.data$lag5match <- compd.data$lag6match <- rep(NA, dim(compd.data)[1])
compd.data$lag7match <- compd.data$lag8match <- rep(NA, dim(compd.data)[1])
compd.data$lag9match <- compd.data$lag10match <- rep(NA, dim(compd.data)[1])

for(i in 1:(dim(compd.data)[1] - 10)){
  if(compd.data$year[i + 1] == compd.data$year[i] + 1)
    {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 1] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 1] %in% disease.set))
      {
      compd.data$lag1match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 1] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 1] %in% disease.set))
      {
      compd.data$lag1match[i] <- 0
    }
  } else {
    compd.data$lag1match[i] <- NA
  }
  
  if(compd.data$year[i + 2] == compd.data$year[i] + 2)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 2] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 2] %in% disease.set))
    {
      compd.data$lag2match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 2] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 2] %in% disease.set))
    {
      compd.data$lag2match[i] <- 0
    }
  } else {
    compd.data$lag2match[i] <- NA
  }
  
  if(compd.data$year[i + 3] == compd.data$year[i] + 3)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 3] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 3] %in% disease.set))
    {
      compd.data$lag3match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 3] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 3] %in% disease.set))
    {
      compd.data$lag3match[i] <- 0
    }
  } else {
    compd.data$lag3match[i] <- NA
  }
  
  if(compd.data$year[i + 4] == compd.data$year[i] + 4)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 4] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 4] %in% disease.set))
    {
      compd.data$lag4match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 4] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 4] %in% disease.set))
    {
      compd.data$lag4match[i] <- 0
    }
  } else {
    compd.data$lag4match[i] <- NA
  }
  
  if(compd.data$year[i + 5] == compd.data$year[i] + 5)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 5] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 5] %in% disease.set))
    {
      compd.data$lag5match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 5] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 5] %in% disease.set))
    {
      compd.data$lag5match[i] <- 0
    }
  } else {
    compd.data$lag5match[i] <- NA
  }
  
  if(compd.data$year[i + 6] == compd.data$year[i] + 6)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 6] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 6] %in% disease.set))
    {
      compd.data$lag6match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 6] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 6] %in% disease.set))
    {
      compd.data$lag6match[i] <- 0
    }
  } else {
    compd.data$lag6match[i] <- NA
  }
  
  if(compd.data$year[i + 7] == compd.data$year[i] + 7)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 7] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 7] %in% disease.set))
    {
      compd.data$lag7match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 7] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 7] %in% disease.set))
    {
      compd.data$lag7match[i] <- 0
    }
  } else {
    compd.data$lag7match[i] <- NA
  }
  
  if(compd.data$year[i + 8] == compd.data$year[i] + 8)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 8] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 8] %in% disease.set))
    {
      compd.data$lag8match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 8] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 8] %in% disease.set))
    {
      compd.data$lag8match[i] <- 0
    }
  } else {
    compd.data$lag8match[i] <- NA
  }
  
  if(compd.data$year[i + 9] == compd.data$year[i] + 9)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 9] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 9] %in% disease.set))
    {
      compd.data$lag9match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 9] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 9] %in% disease.set))
    {
      compd.data$lag9match[i] <- 0
    }
  } else {
    compd.data$lag9match[i] <- NA
  }
  
  if(compd.data$year[i + 10] == compd.data$year[i] + 10)
  {
    if((compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 10] == "HEALTHY") | 
         (compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 10] %in% disease.set))
    {
      compd.data$lag10match[i] <- 1
    } else if((compd.data$CLASS[i] %in% disease.set & compd.data$CLASS[i + 10] == "HEALTHY") |
                (compd.data$CLASS[i] == "HEALTHY" & compd.data$CLASS[i + 10] %in% disease.set))
    {
      compd.data$lag10match[i] <- 0
    }
  } else {
    compd.data$lag10match[i] <- NA
  }
  
  print(i)
}

# lag 1 model
require(lme4)
lag1.glmer <- glmer(lag1match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag2.glmer <- glmer(lag2match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag3.glmer <- glmer(lag3match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag4.glmer <- glmer(lag4match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag5.glmer <- glmer(lag5match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag6.glmer <- glmer(lag6match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag7.glmer <- glmer(lag7match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag8.glmer <- glmer(lag8match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag9.glmer <- glmer(lag9match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)
lag10.glmer <- glmer(lag10match ~ 1 + (1 | Pop), family = "binomial", data = compd.data)

lag1.est <- summary(lag1.glmer)$coef[1]
lag2.est <- summary(lag2.glmer)$coef[1]
lag3.est <- summary(lag3.glmer)$coef[1]
lag4.est <- summary(lag4.glmer)$coef[1]
lag5.est <- summary(lag5.glmer)$coef[1]
lag6.est <- summary(lag6.glmer)$coef[1]
lag7.est <- summary(lag7.glmer)$coef[1]
lag8.est <- summary(lag8.glmer)$coef[1]
lag9.est <- summary(lag9.glmer)$coef[1]
lag10.est <- summary(lag10.glmer)$coef[1]

lag1.se <- 1.96 * summary(lag1.glmer)$coef[2]
lag2.se <- 1.96 * summary(lag2.glmer)$coef[2]
lag3.se <- 1.96 * summary(lag3.glmer)$coef[2]
lag4.se <- 1.96 * summary(lag4.glmer)$coef[2]
lag5.se <- 1.96 * summary(lag5.glmer)$coef[2]
lag6.se <- 1.96 * summary(lag6.glmer)$coef[2]
lag7.se <- 1.96 * summary(lag7.glmer)$coef[2]
lag8.se <- 1.96 * summary(lag8.glmer)$coef[2]
lag9.se <- 1.96 * summary(lag9.glmer)$coef[2]
lag10.se <- 1.96 * summary(lag10.glmer)$coef[2]

acf.vector <- c(lag1.est, lag2.est, lag3.est, lag4.est, lag5.est, lag6.est,
                lag7.est, lag8.est, lag9.est, lag10.est)
acf.se.vector <- c(lag1.se, lag2.se, lag3.se, lag4.se, lag5.se, lag6.se,
                   lag7.se, lag8.se, lag9.se, lag10.se)
lag.vector <- seq(1:10)

plot(exp(acf.vector) ~ lag.vector, ylim = c(0, 6), xlab = "Lag", ylab = "Odds ratio of same status to different status")
for(i in 1:10){
  segments(x0 = i, x1 = i, y0 = exp(acf.vector[i] - acf.se.vector[i]), y1 = exp(acf.vector[i] + acf.se.vector[i]))
}
abline(h = 1, lty = 2, lwd = 2)
