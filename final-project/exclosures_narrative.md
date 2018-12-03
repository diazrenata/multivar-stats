Everything On The Exclosures!
================
Renata Diaz
12/3/2018

``` r
knitr::opts_chunk$set(echo = TRUE)

library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(dplyr)
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-3

``` r
library(ca)
```

### Summer PCoA

``` r
summer_plants_e <- read.csv('data/summer-exclosure-plants-adjusted.csv',
                          stringsAsFactors = F)

summer_plants_e_wis <- vegan::wisconsin(summer_plants_e[,2:ncol(summer_plants_e)])

summer_dist_mat <- vegdist(summer_plants_e_wis, 'bray')
  
summer_pcoa <- cmdscale(summer_dist_mat, k = nrow(summer_plants_e_wis) - 1, eig = T)
```

    ## Warning in cmdscale(summer_dist_mat, k = nrow(summer_plants_e_wis) - 1, :
    ## only 28 of the first 32 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants_e_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  Summer_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   1.8890336 0.16969596 0.1696960
    ##  [2,]   1.3745622 0.12347988 0.2931758
    ##  [3,]   1.1238464 0.10095754 0.3941334
    ##  [4,]   0.8499738 0.07635497 0.4704884
    ##  [5,]   0.7678629 0.06897878 0.5394671
    ##  [6,]   0.7059653 0.06341838 0.6028855
    ##  [7,]   0.6175524 0.05547606 0.6583616
    ##  [8,]   0.4830310 0.04339172 0.7017533
    ##  [9,]   0.4394322 0.03947514 0.7412284
    ## [10,]   0.3879450 0.03484993 0.7760784
    ## [11,]   0.3597845 0.03232022 0.8083986
    ## [12,]   0.3065980 0.02754236 0.8359410
    ## [13,]   0.2828934 0.02541292 0.8613539
    ## [14,]   0.2746660 0.02467384 0.8860277
    ## [15,]   0.2355847 0.02116308 0.9071908

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](exclosures_narrative_files/figure-markdown_github/summer%20x%20pcoa-1.png)

``` r
ordiplot(scores(summer_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Summer plant PCoA")
```

    ## species scores not available

``` r
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(summer_pcoa$points[, 1:2], summer_plants_e_wis)
text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
```

![](exclosures_narrative_files/figure-markdown_github/summer%20x%20pcoa-2.png)

### Winter PCoA

``` r
winter_plants_e <- read.csv('data/winter-exclosure-plants-adjusted.csv',
                          stringsAsFactors = F)

winter_plants_e_wis <- vegan::wisconsin(winter_plants_e[,2:ncol(winter_plants_e)])

winter_dist_mat <- vegdist(winter_plants_e_wis, 'bray')
  
winter_pcoa <- cmdscale(winter_dist_mat, k = nrow(winter_plants_e_wis) - 1, eig = T)
```

    ## Warning in cmdscale(winter_dist_mat, k = nrow(winter_plants_e_wis) - 1, :
    ## only 25 of the first 31 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants_e_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  winter_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   1.9888596 0.20199046 0.2019905
    ##  [2,]   1.3231717 0.13438257 0.3363730
    ##  [3,]   1.0325658 0.10486835 0.4412414
    ##  [4,]   0.8944273 0.09083888 0.5320803
    ##  [5,]   0.6576370 0.06679023 0.5988705
    ##  [6,]   0.5556789 0.05643527 0.6553058
    ##  [7,]   0.4919643 0.04996436 0.7052701
    ##  [8,]   0.4397650 0.04466295 0.7499331
    ##  [9,]   0.3879896 0.03940459 0.7893377
    ## [10,]   0.3298701 0.03350192 0.8228396
    ## [11,]   0.3105875 0.03154356 0.8543832
    ## [12,]   0.2773584 0.02816878 0.8825519
    ## [13,]   0.2364860 0.02401775 0.9065697
    ## [14,]   0.1894455 0.01924026 0.9258099
    ## [15,]   0.1786617 0.01814505 0.9439550

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](exclosures_narrative_files/figure-markdown_github/winter%20x%20pcoa-1.png)

``` r
ordiplot(scores(winter_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "winter plant PCoA")
```

    ## species scores not available

``` r
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(winter_pcoa$points[, 1:2], winter_plants_e_wis)
text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
```

![](exclosures_narrative_files/figure-markdown_github/winter%20x%20pcoa-2.png)

``` r
write.csv(species_pc, 'models/winter_exclosures_species_scores.csv', row.names = T)
```

``` r
rodents <- read.csv('data/exclosure-rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('models/ex_summer_pcoa_vals_3.csv', 
                        stringsAsFactors = F)

winter_axes <- read.csv('models/ex_winter_pcoa_vals_3.csv',
                        stringsAsFactors = F)


pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

rodents <- filter(rodents, year %in% pred_vals$year) %>%
  select(-year)

pred_vals_noy <- select(pred_vals, -year)
pred_vals_y <- select(pred_vals, year)

rodents_hel <- decostand(rodents, 'hellinger')

rodents_prda <- rda(rodents_hel ~ . + Condition(pred_vals_y$year), pred_vals_noy)

R2 <- RsquareAdj(rodents_prda)$r.squared
R2adj <- RsquareAdj(rodents_prda)$adj.r.squared

R2
```

    ## [1] 0.3623587

``` r
R2adj
```

    ## [1] NA

``` r
anova(rodents_prda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + winterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)    
    ## Model     6 0.068943 5.2442  0.001 ***
    ## Residual 21 0.046013                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rodents_prda, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + winterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.048636 22.1971  0.001 ***
    ## RDA2      1 0.010603  4.8391  0.011 *  
    ## RDA3      1 0.006234  2.8454  0.135    
    ## RDA4      1 0.002426  1.1073  0.877    
    ## RDA5      1 0.000573  0.2614  1.000    
    ## RDA6      1 0.000471  0.2150  0.994    
    ## Residual 21 0.046013                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Find the most parsimonious model...

``` r
set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_prda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                   R2.adjusted
    ## + winterPCoAxis_1  0.47598660
    ## + SummerPCoAxis_1  0.21921918
    ## + winterPCoAxis_2  0.06372072
    ## + SummerPCoAxis_2  0.03540181
    ## + SummerPCoAxis_3  0.02831296
    ## <none>             0.00000000
    ## + winterPCoAxis_3 -0.02073203
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + winterPCoAxis_1  1 -64.935 26.434  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.4759866 
    ## Call: rodents_hel ~ winterPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + winterPCoAxis_2   0.5563157
    ## + SummerPCoAxis_1   0.5361416
    ## + SummerPCoAxis_3   0.5175246
    ## + SummerPCoAxis_2   0.4987705
    ## <none>              0.4759866
    ## + winterPCoAxis_3   0.4746232
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + winterPCoAxis_2  1 -68.855 5.8884  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5563157 
    ## Call: rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_1   0.5998437
    ## + SummerPCoAxis_3   0.5859096
    ## + SummerPCoAxis_2   0.5732289
    ## + winterPCoAxis_3   0.5586413
    ## <none>              0.5563157
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + SummerPCoAxis_1  1 -70.987 3.8282  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5998437 
    ## Call: rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_3   0.6333053
    ## + SummerPCoAxis_2   0.6072131
    ## + winterPCoAxis_3   0.6008931
    ## <none>              0.5998437
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + SummerPCoAxis_3  1 -72.703 3.2813  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6333053 
    ## Call: rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 +      SummerPCoAxis_3 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_2   0.6424874
    ## + winterPCoAxis_3   0.6362336
    ## <none>              0.6333053
    ## 
    ##                   Df     AIC      F Pr(>F)
    ## + SummerPCoAxis_2  1 -72.673 1.6164  0.122

``` r
# Most parsimonious is Call: rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 

rod_prda_pars <- rda(rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 + Condition(pred_vals_y$year), pred_vals_noy)

pR2p <- RsquareAdj(rod_prda_pars)$r.squared
pR2adjp <- RsquareAdj(rod_prda_pars)$adj.r.squared

pR2p
```

    ## [1] 0.282135

``` r
pR2adjp
```

    ## [1] NA

``` r
anova(rod_prda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)    
    ## Model     3 0.053679 7.0082  0.001 ***
    ## Residual 24 0.061276                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rod_prda_pars, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ winterPCoAxis_1 + winterPCoAxis_2 + SummerPCoAxis_1 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.044688 17.5028  0.001 ***
    ## RDA2      1 0.008148  3.1912  0.017 *  
    ## RDA3      1 0.000844  0.3306  0.962    
    ## Residual 24 0.061276                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rod_ppart <- varpart(rodents_hel, ~ winterPCoAxis_1, ~  winterPCoAxis_2, ~SummerPCoAxis_1, ~ pred_vals_y$year, data = pred_vals_noy)
rod_ppart
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~winterPCoAxis_1,
    ## ~winterPCoAxis_2, ~SummerPCoAxis_1, ~pred_vals_y$year, data =
    ## pred_vals_noy)
    ## 
    ## Explanatory tables:
    ## X1:  ~winterPCoAxis_1
    ## X2:  ~winterPCoAxis_2
    ## X3:  ~SummerPCoAxis_1
    ## X4:  ~pred_vals_y$year 
    ## 
    ## No. of explanatory tables: 4 
    ## Total variation (SS): 5.3273 
    ##             Variance: 0.19026 
    ## No. of observations: 29 
    ## 
    ## Partition table:
    ##                             Df R.square Adj.R.square Testable
    ## [aeghklno] = X1              1  0.49470      0.47599     TRUE
    ## [befiklmo] = X2              1  0.09716      0.06372     TRUE
    ## [cfgjlmno] = X3              1  0.24710      0.21922     TRUE
    ## [dhijkmno] = X4              1  0.39580      0.37342     TRUE
    ## [abefghiklmno] = X1+X2       2  0.58801      0.55632     TRUE
    ## [acefghjklmno] = X1+X3       2  0.56927      0.53614     TRUE
    ## [adeghijklmno] = X1+X4       2  0.60729      0.57709     TRUE
    ## [bcefgijklmno] = X2+X3       2  0.29229      0.23785     TRUE
    ## [bdefhijklmno] = X2+X4       2  0.46171      0.42030     TRUE
    ## [cdfghijklmno] = X3+X4       2  0.43692      0.39360     TRUE
    ## [abcefghijklmno] = X1+X2+X3  3  0.64272      0.59984     TRUE
    ## [abdefghijklmno] = X1+X2+X4  3  0.62285      0.57759     TRUE
    ## [acdefghijklmno] = X1+X3+X4  3  0.66387      0.62353     TRUE
    ## [bcdefghijklmno] = X2+X3+X4  3  0.50302      0.44338     TRUE
    ## [abcdefghijklmno] = All      4  0.67794      0.62426     TRUE
    ## Individual fractions                                         
    ## [a] = X1 | X2+X3+X4          1               0.18088     TRUE
    ## [b] = X2 | X1+X3+X4          1               0.00073     TRUE
    ## [c] = X3 | X1+X2+X4          1               0.04667     TRUE
    ## [d] = X4 | X1+X2+X3          1               0.02442     TRUE
    ## [e]                          0               0.04906    FALSE
    ## [f]                          0              -0.00022    FALSE
    ## [g]                          0              -0.02359    FALSE
    ## [h]                          0               0.18112    FALSE
    ## [i]                          0               0.06298    FALSE
    ## [j]                          0              -0.00314    FALSE
    ## [k]                          0              -0.09413    FALSE
    ## [l]                          0              -0.00268    FALSE
    ## [m]                          0               0.01685    FALSE
    ## [n]                          0               0.15419    FALSE
    ## [o]                          0               0.03115    FALSE
    ## [p] = Residuals              0               0.37574    FALSE
    ## Controlling 2 tables X                                       
    ## [ae] = X1 | X3+X4            1               0.22993     TRUE
    ## [ag] = X1 | X2+X4            1               0.15729     TRUE
    ## [ah] = X1 | X2+X3            1               0.36199     TRUE
    ## [be] = X2 | X3+X4            1               0.04978     TRUE
    ## [bf] = X2 | X1+X4            1               0.00051     TRUE
    ## [bi] = X2 | X1+X3            1               0.06370     TRUE
    ## [cf] = X3 | X1+X4            1               0.04645     TRUE
    ## [cg] = X3 | X2+X4            1               0.02308     TRUE
    ## [cj] = X3 | X1+X2            1               0.04353     TRUE
    ## [dh] = X4 | X2+X3            1               0.20553     TRUE
    ## [di] = X4 | X1+X3            1               0.08739     TRUE
    ## [dj] = X4 | X1+X2            1               0.02128     TRUE
    ## Controlling 1 table X                                        
    ## [aghn] = X1 | X2             1               0.49259     TRUE
    ## [aehk] = X1 | X3             1               0.31692     TRUE
    ## [aegl] = X1 | X4             1               0.20366     TRUE
    ## [bfim] = X2 | X1             1               0.08033     TRUE
    ## [beik] = X2 | X3             1               0.01863     TRUE
    ## [befl] = X2 | X4             1               0.04688     TRUE
    ## [cfjm] = X3 | X1             1               0.06016     TRUE
    ## [cgjn] = X3 | X2             1               0.17413     TRUE
    ## [cfgl] = X3 | X4             1               0.02018     TRUE
    ## [dijm] = X4 | X1             1               0.10110     TRUE
    ## [dhjn] = X4 | X2             1               0.35658     TRUE
    ## [dhik] = X4 | X3             1               0.17438     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_ppart, digits = 2)
```

![](exclosures_narrative_files/figure-markdown_github/partial%20variance%20partitioning-1.png)
