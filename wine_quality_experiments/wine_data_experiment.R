require(ggplot2, gridExtra, mlr, parallelMap)

red <- read.table("winequality/winequality-red.csv", header=TRUE, sep=";")
white <- read.table("winequality/winequality-white.csv", header=TRUE, sep=";")

red$colour <- "red"; white$colour <- "white"

wines <- rbind(red, white)
wines$colour <- as.factor(wines$colour)
wines <- wines[sample(nrow(wines)),]  # makes plots look a bit nicer

summary(wines)

plot_colours <- ifelse(as.matrix(wines$colour)=="white", "grey", "red")
plot_colours <- adjustcolor(plot_colours, alpha.f=0.3)

pairs(wines[, 1:4], col=plot_colours, pch=16)
pairs(wines[, 5:8], col=plot_colours, pch=16)
pairs(wines[, 9:12], col=plot_colours, pch=16)

# sadly aes won't take colnames(wines)[i] properly so have select manually
density_plot1 <- ggplot(wines, aes(fixed.acidity, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[1])
density_plot2 <- ggplot(wines, aes(volatile.acidity, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[2])
density_plot3 <- ggplot(wines, aes(citric.acid, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[3])
density_plot4 <- ggplot(wines, aes(residual.sugar, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[4])
density_plot5 <- ggplot(wines, aes(chlorides, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[5])
density_plot6 <- ggplot(wines, aes(free.sulfur.dioxide, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[6])
density_plot7 <- ggplot(wines, aes(total.sulfur.dioxide, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[7])
density_plot8 <- ggplot(wines, aes(density, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[8])
density_plot9 <- ggplot(wines, aes(pH, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[9])
density_plot10 <- ggplot(wines, aes(sulphates, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[10])
density_plot11 <- ggplot(wines, aes(alcohol, colour=colour)) + 
  geom_density() + theme_classic() + xlab(colnames(wines)[11])
density_plot12 <- ggplot(wines, aes(quality, colour=colour)) +      # replace with a histogram
  geom_histogram() + theme_classic() + xlab(colnames(wines)[12])

# TODO: sort out colours
density_grid <- grid.arrange(density_plot1, density_plot2, density_plot3, 
                             density_plot4, density_plot5, density_plot6,
                             density_plot7, density_plot8, density_plot9,
                             density_plot10, density_plot11, density_plot12, ncol=3)

cor(wines[, -13])

# don't see any reason to do any transformations yet

##################################################################################################
##################################################################################################
###                                                                                            ###  
###                             Benchmark Experiment (i)                                       ###                                    
###                                                                                            ###
##################################################################################################
##################################################################################################

# use parallelisation for benchmarking to speed it up
parallelStartSocket(3, level="mlr.benchmark")

# make the classif & regr tasks
task.c <- makeClassifTask(id="classif", wines, "quality")
task.c.col <- makeClassifTask(id="classif.col", wines[c(12, 13)], "quality")
task.c.chem <- makeClassifTask(id="classif.chem", wines[-13], "quality")

task.r <- makeRegrTask(id="regr", wines, "quality")
task.r.col <- makeRegrTask(id="regr.col", wines[c(12, 13)], "quality")
task.r.chem <- makeRegrTask(id="regr.chem", wines[-13], "quality")


# make the inner & outer validation objects; CV for tuning, bootstrap for outer. 
inner <- makeResampleDesc("CV", iters=3)
outer <- makeResampleDesc("Bootstrap", iters=100) # should do more iters but ends up being slow
control <- makeTuneControlGrid() 


# param tuning; use the inner validation set-up for selected parameter tuning.
ps <- makeParamSet(
  makeDiscreteParam("ntree", values=seq(10, 50, 10))
)
rf.c <- makeTuneWrapper("classif.randomForest", resampling=inner, par.set=ps, 
                        control=control, show.info=FALSE)
rf.r <- makeTuneWrapper("regr.randomForest", resampling=inner, par.set=ps, 
                        control=control, show.info=FALSE)

ps <- makeParamSet(
  makeDiscreteParam("C", values=seq(0.5, 1.5, 0.5))
)
svm.c <- makeTuneWrapper("classif.ksvm", resampling=inner, par.set=ps, 
                         control=control, show.info=FALSE)
svm.r <- makeTuneWrapper("regr.ksvm", resampling=inner, par.set=ps, 
                         control=control, show.info=FALSE)

ps <- makeParamSet(
  makeDiscreteParam("size", values=1:3)
)
nnet.c <- makeTuneWrapper("classif.nnet", resampling=inner, par.set=ps, 
                          control=control, show.info=FALSE)
nnet.r <- makeTuneWrapper("regr.nnet", resampling=inner, par.set=ps, 
                          control=control, show.info=FALSE)

# make some baseline learners; we will use featureless for both, i.e. average over the labels
baseline.c <- makeLearner("classif.featureless")
baseline.r <- makeLearner("regr.featureless")

# group all our learners into lists for benchmarking
learners.c <- list(rf.c, svm.c, nnet.c, baseline.c)
learners.r <- list(rf.r, svm.r, nnet.r, baseline.r)

# choose our validation metrics; we use MMCE for classification and MSE, MAE for regression
metrics.c <- list(mmce)
metrics.r <- list(mse, mae)

# run the benchmark experiments
set.seed(100)
bmresults.c <- benchmark(learners=learners.c, tasks=list(task.c, task.c.col, task.c.chem), 
                         resamplings=outer, measures=metrics.c)
bmresults.r <- benchmark(learners=learners.r, tasks=list(task.r, task.r.col, task.r.chem), 
                         resamplings=outer, measures=metrics.r)

parallelStop() # we aren't doing any more parallelisation so turn it off

# store the validation samples for the two experiments as dataframes
bootstrap_results.c <- getBMRPerformances(bmresults.c, as.df=TRUE)
bootstrap_results.r <- getBMRPerformances(bmresults.r, as.df=TRUE)

# for comparison we can do Wilcoxon signed-rank tests for differences between validation samples
# and calculate metric standard errors for confidence intervals. 

# Wilcoxon signed-rank test for RF and SVM
wilcox.test(subset(bootstrap_results.c, learner.id=="classif.randomForest.tuned")$mmce,
            subset(bootstrap_results.c, learner.id=="classif.ksvm.tuned")$mmce)

# bootstrap estimate for MMCE 
mean(subset(bootstrap_results.c, (learner.id=="classif.randomForest.tuned" &
                                  task.id=="classif.chem"))$mmce)

# bootstrap estimate for MMCE standard error 
sd(subset(bootstrap_results.c, (learner.id=="classif.randomForest.tuned" &
                                task.id=="classif"))$mmce)

# we might also want to check the regularity of our bootstrap validation samples to check 
# inference is legitimate; we want symmetric looking distributions
hist(subset(bootstrap_results.c, learner.id=="classif.randomForest.tuned")$mmce)





















