library(vegan)
library(dplyr)

# import data
spe <- read.csv("week-11/DoubsSpe.csv", row.names = 1)
env <- read.csv("week-11/DoubsEnv.csv", row.names = 1)

# remove site w no species (8)
sitesums <- rowSums(spe)
which(sitesums == 0)

spe <- spe[-which(sitesums == 0), ]
env <- env[-which(sitesums == 0), ] %>%
  select(-das)

pen2 <- rep("very_steep", nrow(env))
pen2[env$pen <= quantile(env$pen)[4]] <- "steep"
pen2[env$pen <= quantile(env$pen)[3]] <- "moderate"
pen2[env$pen <= quantile(env$pen)[2]] <- "low"
pen2 <- factor(pen2, levels = c("low", "moderate", "steep", "very_steep"))
table(pen2)

env2 <- env
env2$pen <- pen2

# Do Hellinger transformation (sqrt of the row normalized data)
# Use vegan::decostand

spe.hel <- decostand(spe, 'hellinger')

# Run an rda using vegan::rda

spe.rda <- rda(spe.hel ~ ., env2)

# unadjusted and adjusted r2
R2 <- RsquareAdj(spe.rda)$r.squared
R2adj <- RsquareAdj(spe.rda)$adj.r.squared

# plot w f scores
plot(spe.rda, scaling = 1, main = "Triplot RDA spe.hel ~ env2 - scaling 1 - w a scores")
spe.sc <- scores(spe.rda, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc[,1], spe.sc[,2], length = 0, lty = 1, col = 'red')

# plot w z scores
plot(spe.rda, scaling = 1, display = c("sp", "lc", "cn"), main = "Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length = 0, lty = 1, col = "red")

anova(spe.rda, step = 1000)
anova(spe.rda, by = "axis", step = 1000)


# Reduce the number of variables for the most parsimonious model.

set.seed(11)
step.forward <- ordiR2step(rda(spe.hel ~ 1, data = env2), scope = formula(spe.rda), 
                           R2scope = F, direction = "forward", pstep = 1000)

# Most parsimonious is Call: spe.hel ~ alt + oxy + dbo 

spe.rda.pars <- rda(spe.hel ~ alt + oxy + dbo, env2)

R2p <- RsquareAdj(spe.rda.pars)$r.squared
R2adjp <- RsquareAdj(spe.rda.pars)$adj.r.squared

# plot w f scores
plot(spe.rda.pars, scaling = 1, main = "Triplot RDA spe.hel ~ alt + oxy + dbo, env2 - scaling 1 - w a scores")
spe.sc.p <- scores(spe.rda.pars, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc.p[,1], spe.sc.p[,2], length = 0, lty = 1, col = 'red')

# plot w z scores
plot(spe.rda.pars, scaling = 1, display = c("sp", "lc", "cn"), main = "Triplot RDA spe.hel ~ alt + oxy + dbo, env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc.p[, 1], spe.sc.p[, 2], length = 0, lty = 1, col = "red")

anova(spe.rda.pars, step = 1000)
anova(spe.rda.pars, by = "axis", step = 1000)


# Partial RDA

partial.alt <- rda(spe.hel ~ alt + Condition(oxy + dbo + dur), data = env2)
anova(partial.alt, step = 1000)

# Variance partitioning
??varpart
spe.part <- varpart(spe.hel, ~alt, ~oxy, ~dur, ~dbo, data = env2)
spe.part

plot(spe.part, digits = 2)
