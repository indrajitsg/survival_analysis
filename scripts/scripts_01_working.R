### Survival Analysis in R ###
# https://webr.sh

# Load libraries
library(survival)
library(ggfortify)
library(SMPracticals)
library(flexsurv)

setwd("C:/Stuff/DeepCortex Dropbox/Indrajit Sen Gupta/Training/survival_analysis")

# Load dataset
prost <- read.table("prostate_cancer.txt")
colnames(prost) <- c("patient", "treatment", "time", "status",
                    "age", "sh", "size", "index")

head(prost)
dim(prost)

prost_sv <- Surv(time=prost$time, event=prost$status)
prost_sv

# Simple Kaplan Meier plot
mykm1 <- survfit(prost_sv ~ 1, data = prost)
plot(mykm1, mark.time = TRUE)
summary(mykm1)

# Stratified KM plot
mykm2 <- survfit(prost_sv ~ treatment, data = prost)
## alternative 2 line chart
plot(mykm2, mark.time = TRUE)

# Advanced plot
autoplot(mykm2)

# Log Rank Test
survdiff(prost_sv ~ treatment, data=prost)


# Exercise Kaplan Meier
head(motorette)

plot(motorette$x, motorette$y)

motor_sv <- Surv(time=motorette$y, event=motorette$cens)
motor_km1 <- survfit(motor_sv ~ 1, data=motorette)
plot(motor_km1, mark.time = TRUE)

# Create a binary variable
motorette$temp_binary <- ifelse(motorette$x < 180, "low", "high")

motor_km2 <- survfit(motor_sv ~ temp_binary, data=motorette)
plot(motor_km2, mark.time=TRUE)

# Perform log-rank test
survdiff(motor_sv ~ temp_binary, data=motorette)

# Accelerated Failure Time Model
aft1 <- survreg(prost_sv ~ treatment + age + sh + size + index, 
                 data = prost, 
                 dist = "weibull")

summary(aft1)

# Cox Proportional Hazards Model
cox <- coxph(prost_sv ~ treatment + age + sh + size + index,  data = prost)
summary(cox)

# Getting a plot of the model
autoplot(survfit(cox))

# Aalens additive regression
aalen <- aareg(prost_sv ~ treatment + age + sh+ size + index,
			   data = prost)

summary(aalen)
autoplot(aalen)

# Parametric Models
weib_mod <- flexsurvreg(prost_sv ~ treatment + age + sh+ size + index,
			data = prost, dist = "weibull")

exp_mod <- flexsurvreg(prost_sv ~ treatment + age + sh+ size + index,
			data = prost, dist = "exp")

plot(weib_mod)

lines(exp_mod, col="blue")

# ------------------------------------
# Exercise: Fit Cox PH Model to udca1
# ------------------------------------
data(udca)
head(udca1)

# Create a survival object
udca_sv <- Surv(time=udca1$futime, event=udca1$status)
udca_sv

udca_cox <- coxph(udca_sv ~ trt + stage + bili + riskscore,  data = udca1)
summary(udca_cox)
autoplot(survfit(udca_cox))
AIC(udca_cox)

udca_cox1 <- coxph(udca_sv ~ trt + bili + riskscore,  data = udca1)
summary(udca_cox1)
AIC(udca_cox1)

udca_cox2 <- coxph(udca_sv ~ trt + riskscore,  data = udca1)
summary(udca_cox2)
AIC(udca_cox2)

udca_km <- survfit(udca_sv ~ trt, data=udca1)
autoplot(udca_km)


