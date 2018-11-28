Multivar stats final project
================
Renata Diaz
11/23/2018

Project plan
------------

1.  Get data
    1.  load Portal rodent and plant data
    2.  summarize by year
    3.  standardize according to trapping/survey effort

2.  Reduce dimensionality of plant data
    1.  try LDA
        1.  If base LDA doesn't work, try removing 5% most/least common species

    2.  pcoa
    3.  separate winter/summer censuses a priori?

3.  Distance based redundancy analysis
    1.  try to use plant summary axes to predict rodent community
    2.  (?) Use LDA rodent topics as ind variable?

4.  Variance partitioning to isolate effects of year, plants on rodents

Get data
--------

### Plant data

-   Extracted summer & winter plant censuses for all years
-   Standardized plant abundances according to sampling effort
-   Kept seasons separate

### Rodent data

-   Rodent censuses on control plots for all years, restricted to granivores
-   Standardized according to sampling effort (per census period)
-   Summed across all months in each calendar year

PCoA
----

-   Ran separate principle coordinates analyses on winter and summer datasets.

### Summer PCoA

``` r
summer_plants <- read.csv('data/summer-plants-adjusted.csv',
                          stringsAsFactors = F)

summer_dist_mat <- vegdist(summer_plants[,2:ncol(summer_plants)], 'bray')
  
summer_pcoa <- cmdscale(summer_dist_mat, k = nrow(summer_plants) - 1, eig = T)
```

    ## Warning in cmdscale(summer_dist_mat, k = nrow(summer_plants) - 1, eig = T):
    ## only 24 of the first 31 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  Summer_PCoA_Table
```

    ##         eigenvalues       propVar    cumVar
    ##  [1,]  2.339598e+00  2.125538e-01 0.2125538
    ##  [2,]  1.904286e+00  1.730054e-01 0.3855593
    ##  [3,]  1.285763e+00  1.168123e-01 0.5023716
    ##  [4,]  1.123691e+00  1.020880e-01 0.6044596
    ##  [5,]  7.959522e-01  7.231271e-02 0.6767723
    ##  [6,]  5.829654e-01  5.296274e-02 0.7297350
    ##  [7,]  4.955673e-01  4.502257e-02 0.7747576
    ##  [8,]  4.403899e-01  4.000967e-02 0.8147673
    ##  [9,]  4.093077e-01  3.718584e-02 0.8519531
    ## [10,]  3.597252e-01  3.268124e-02 0.8846344
    ## [11,]  2.711839e-01  2.463721e-02 0.9092716
    ## [12,]  2.573015e-01  2.337598e-02 0.9326476
    ## [13,]  1.754937e-01  1.594370e-02 0.9485913
    ## [14,]  1.540241e-01  1.399318e-02 0.9625844
    ## [15,]  1.257192e-01  1.142166e-02 0.9740061
    ## [16,]  9.808451e-02  8.911033e-03 0.9829171
    ## [17,]  8.454061e-02  7.680563e-03 0.9905977
    ## [18,]  7.117380e-02  6.466180e-03 0.9970639
    ## [19,]  6.361291e-02  5.779269e-03 1.0028431
    ## [20,]  5.101497e-02  4.634739e-03 1.0074779
    ## [21,]  4.444792e-02  4.038119e-03 1.0115160
    ## [22,]  2.341441e-02  2.127212e-03 1.0136432
    ## [23,]  1.595917e-02  1.449899e-03 1.0150931
    ## [24,]  5.869408e-03  5.332390e-04 1.0156263
    ## [25,] -3.330669e-16 -3.025932e-17 1.0156263
    ## [26,] -5.337132e-04 -4.848815e-05 1.0155779
    ## [27,] -6.190269e-03 -5.623895e-04 1.0150155
    ## [28,] -2.755041e-02 -2.502971e-03 1.0125125
    ## [29,] -3.537574e-02 -3.213906e-03 1.0092986
    ## [30,] -4.316576e-02 -3.921634e-03 1.0053770
    ## [31,] -5.918462e-02 -5.376956e-03 1.0000000

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](narrative_files/figure-markdown_github/summer%20pcoa-1.png)

There seems to be an inflection point in the scree plot around axis 3.

Moving forward, keeping the first 3 axes as predictor variables for the rodent community.

### Winter PCoA

``` r
winter_plants <- read.csv('data/winter-plants-adjusted.csv',
                          stringsAsFactors = F)

winter_dist_mat <- vegdist(winter_plants[,2:ncol(winter_plants)], 'bray')
  
winter_pcoa <- cmdscale(winter_dist_mat, k = nrow(winter_plants) - 1, eig = T)
```

    ## Warning in cmdscale(winter_dist_mat, k = nrow(winter_plants) - 1, eig = T):
    ## only 22 of the first 30 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  winter_PCoA_Table
```

    ##         eigenvalues       propVar    cumVar
    ##  [1,]  2.952057e+00  2.644179e-01 0.2644179
    ##  [2,]  1.412192e+00  1.264910e-01 0.3909089
    ##  [3,]  1.151511e+00  1.031417e-01 0.4940506
    ##  [4,]  9.806955e-01  8.784160e-02 0.5818922
    ##  [5,]  8.641039e-01  7.739841e-02 0.6592906
    ##  [6,]  5.451243e-01  4.882717e-02 0.7081177
    ##  [7,]  5.196633e-01  4.654661e-02 0.7546644
    ##  [8,]  4.824972e-01  4.321762e-02 0.7978820
    ##  [9,]  4.492552e-01  4.024011e-02 0.8381221
    ## [10,]  4.191044e-01  3.753948e-02 0.8756616
    ## [11,]  3.301053e-01  2.956777e-02 0.9052293
    ## [12,]  2.762754e-01  2.474618e-02 0.9299755
    ## [13,]  2.119215e-01  1.898196e-02 0.9489575
    ## [14,]  1.506384e-01  1.349279e-02 0.9624503
    ## [15,]  1.455869e-01  1.304032e-02 0.9754906
    ## [16,]  1.253993e-01  1.123210e-02 0.9867227
    ## [17,]  7.765375e-02  6.955502e-03 0.9936782
    ## [18,]  7.231843e-02  6.477613e-03 1.0001558
    ## [19,]  6.857242e-02  6.142080e-03 1.0062979
    ## [20,]  4.531939e-02  4.059290e-03 1.0103572
    ## [21,]  2.197016e-02  1.967883e-03 1.0123251
    ## [22,]  1.877693e-02  1.681863e-03 1.0140069
    ## [23,] -4.857226e-17 -4.350652e-18 1.0140069
    ## [24,] -2.456597e-04 -2.200391e-05 1.0139849
    ## [25,] -3.498485e-03 -3.133618e-04 1.0136716
    ## [26,] -1.026284e-02 -9.192502e-04 1.0127523
    ## [27,] -2.000302e-02 -1.791685e-03 1.0109606
    ## [28,] -2.913202e-02 -2.609376e-03 1.0083513
    ## [29,] -3.530142e-02 -3.161973e-03 1.0051893
    ## [30,] -5.793504e-02 -5.189283e-03 1.0000000

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](narrative_files/figure-markdown_github/winter%20pcoa-1.png)

There seems to be an inflection point in the scree plot around axis 2 or 3. Since stopping at 2 would only capture 39% of variation, going to go for 3.

RDA
---

Redundancy analysis, using combined winter and summer axes and year to predict the rodent community.

Restricted to years with both a winter & summer census (n = 27).

``` r
rodents <- read.csv('data/rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('models/Summer_pcoa_vals_3.csv', 
                        stringsAsFactors = F)

winter_axes <- read.csv('models/Winter_pcoa_vals_3.csv',
                        stringsAsFactors = F)


pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

rodents <- filter(rodents, year %in% pred_vals$year) %>%
  select(-year)

rodents_hel <- decostand(rodents, 'hellinger')

rodents_rda <- rda(rodents_hel ~ ., pred_vals)

R2 <- RsquareAdj(rodents_rda)$r.squared
R2adj <- RsquareAdj(rodents_rda)$adj.r.squared

R2
```

    ## [1] 0.7553613

``` r
R2adj
```

    ## [1] 0.6652313

``` r
anova(rodents_rda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ year + WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     7  0.13620 8.3808  0.001 ***
    ## Residual 19  0.04411                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rodents_rda, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ year + WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3, data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.100080 43.1083  0.001 ***
    ## RDA2      1 0.018872  8.1288  0.001 ***
    ## RDA3      1 0.009484  4.0851  0.035 *  
    ## RDA4      1 0.004300  1.8520  0.550    
    ## RDA5      1 0.002483  1.0697  0.896    
    ## RDA6      1 0.000846  0.3645  1.000    
    ## RDA7      1 0.000133  0.0573  1.000    
    ## Residual 19 0.044110                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Find the most parsimonious model...

``` r
set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_rda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_1  0.44246414
    ## + year             0.39300101
    ## + SummerPCoAxis_2  0.10562176
    ## + SummerPCoAxis_3  0.03158667
    ## + WinterPCoAxis_2  0.01695581
    ## + WinterPCoAxis_3  0.01672135
    ## + SummerPCoAxis_1  0.01071350
    ## <none>             0.00000000
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_1  1 -60.105 21.634  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.4424641 
    ## Call: rodents_hel ~ WinterPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + year              0.5009630
    ## + SummerPCoAxis_2   0.4985225
    ## + WinterPCoAxis_3   0.4865958
    ## + WinterPCoAxis_2   0.4822580
    ## + SummerPCoAxis_1   0.4682159
    ## + SummerPCoAxis_3   0.4544798
    ## <none>              0.4424641
    ## 
    ##        Df     AIC      F Pr(>F)   
    ## + year  1 -62.201 3.9306  0.004 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.500963 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + year 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_2   0.5600700
    ## + WinterPCoAxis_2   0.5553529
    ## + WinterPCoAxis_3   0.5508614
    ## + SummerPCoAxis_1   0.5341521
    ## + SummerPCoAxis_3   0.5158944
    ## <none>              0.5009630
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + SummerPCoAxis_2  1 -64.753 4.2245  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.56007 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_2   0.6132639
    ## + SummerPCoAxis_1   0.5974728
    ## + WinterPCoAxis_3   0.5777751
    ## + SummerPCoAxis_3   0.5670713
    ## <none>              0.5600700
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_2  1 -67.433 4.1636  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6132639 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_1   0.6490379
    ## + WinterPCoAxis_3   0.6324017
    ## + SummerPCoAxis_3   0.6197349
    ## <none>              0.6132639
    ## 
    ##                   Df    AIC      F Pr(>F)   
    ## + SummerPCoAxis_1  1 -69.31 3.2425  0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6490379 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 +      SummerPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_3   0.6610786
    ## + SummerPCoAxis_3   0.6529810
    ## <none>              0.6490379
    ## 
    ##                   Df    AIC      F Pr(>F)  
    ## + WinterPCoAxis_3  1 -69.57 1.7461  0.066 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Most parsimonious is Call: rodents_hel ~ Winter 1 + year + Summer 2 + Winter 2 

rod_rda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 + SummerPCoAxis_1 , pred_vals)

R2p <- RsquareAdj(rod_rda_pars)$r.squared
R2adjp <- RsquareAdj(rod_rda_pars)$adj.r.squared

R2p
```

    ## [1] 0.7165306

``` r
R2adjp
```

    ## [1] 0.6490379

``` r
anova(rod_rda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 + SummerPCoAxis_1, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     5 0.129196 10.616  0.001 ***
    ## Residual 21 0.051112                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rod_rda_pars, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 + SummerPCoAxis_1, data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.099422 40.8491  0.001 ***
    ## RDA2      1 0.016304  6.6986  0.001 ***
    ## RDA3      1 0.008738  3.5900  0.007 ** 
    ## RDA4      1 0.003504  1.4399  0.431    
    ## RDA5      1 0.001228  0.5046  0.871    
    ## Residual 21 0.051112                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Variance partitioning
---------------------

``` r
# Variance partitioning
rod_part <- varpart(rodents_hel, ~WinterPCoAxis_1, ~SummerPCoAxis_2, ~WinterPCoAxis_2, ~SummerPCoAxis_1, data = pred_vals)
rod_part
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~WinterPCoAxis_1,
    ## ~SummerPCoAxis_2, ~WinterPCoAxis_2, ~SummerPCoAxis_1, data =
    ## pred_vals)
    ## 
    ## Explanatory tables:
    ## X1:  ~WinterPCoAxis_1
    ## X2:  ~SummerPCoAxis_2
    ## X3:  ~WinterPCoAxis_2
    ## X4:  ~SummerPCoAxis_1 
    ## 
    ## No. of explanatory tables: 4 
    ## Total variation (SS): 4.688 
    ##             Variance: 0.18031 
    ## No. of observations: 27 
    ## 
    ## Partition table:
    ##                             Df R.square Adj.R.square Testable
    ## [aeghklno] = X1              1  0.46391      0.44246     TRUE
    ## [befiklmo] = X2              1  0.14002      0.10562     TRUE
    ## [cfgjlmno] = X3              1  0.05477      0.01696     TRUE
    ## [dhijkmno] = X4              1  0.04876      0.01071     TRUE
    ## [abefghiklmno] = X1+X2       2  0.53710      0.49852     TRUE
    ## [acefghjklmno] = X1+X3       2  0.52208      0.48226     TRUE
    ## [adeghijklmno] = X1+X4       2  0.50912      0.46822     TRUE
    ## [bcefgijklmno] = X2+X3       2  0.18513      0.11722     TRUE
    ## [bdefhijklmno] = X2+X4       2  0.18556      0.11769     TRUE
    ## [cdfghijklmno] = X3+X4       2  0.12360      0.05057     TRUE
    ## [abcefghijklmno] = X1+X2+X3  3  0.58843      0.53475     TRUE
    ## [abdefghijklmno] = X1+X2+X4  3  0.58215      0.52765     TRUE
    ## [acdefghijklmno] = X1+X3+X4  3  0.56569      0.50904     TRUE
    ## [bcdefghijklmno] = X2+X3+X4  3  0.24539      0.14697     TRUE
    ## [abcdefghijklmno] = All      4  0.63178      0.56483     TRUE
    ## Individual fractions                                         
    ## [a] = X1 | X2+X3+X4          1               0.41786     TRUE
    ## [b] = X2 | X1+X3+X4          1               0.05579     TRUE
    ## [c] = X3 | X1+X2+X4          1               0.03718     TRUE
    ## [d] = X4 | X1+X2+X3          1               0.03008     TRUE
    ## [e]                          0               0.04061    FALSE
    ## [f]                          0               0.00364    FALSE
    ## [g]                          0              -0.00791    FALSE
    ## [h]                          0              -0.00033    FALSE
    ## [i]                          0              -0.00330    FALSE
    ## [j]                          0              -0.00095    FALSE
    ## [k]                          0               0.00716    FALSE
    ## [l]                          0               0.00694    FALSE
    ## [m]                          0              -0.00008    FALSE
    ## [n]                          0              -0.01672    FALSE
    ## [o]                          0              -0.00515    FALSE
    ## [p] = Residuals              0               0.43517    FALSE
    ## Controlling 2 tables X                                       
    ## [ae] = X1 | X3+X4            1               0.45847     TRUE
    ## [ag] = X1 | X2+X4            1               0.40995     TRUE
    ## [ah] = X1 | X2+X3            1               0.41753     TRUE
    ## [be] = X2 | X3+X4            1               0.09640     TRUE
    ## [bf] = X2 | X1+X4            1               0.05943     TRUE
    ## [bi] = X2 | X1+X3            1               0.05249     TRUE
    ## [cf] = X3 | X1+X4            1               0.04083     TRUE
    ## [cg] = X3 | X2+X4            1               0.02927     TRUE
    ## [cj] = X3 | X1+X2            1               0.03623     TRUE
    ## [dh] = X4 | X2+X3            1               0.02975     TRUE
    ## [di] = X4 | X1+X3            1               0.02678     TRUE
    ## [dj] = X4 | X1+X2            1               0.02912     TRUE
    ## Controlling 1 table X                                        
    ## [aghn] = X1 | X2             1               0.39290     TRUE
    ## [aehk] = X1 | X3             1               0.46530     TRUE
    ## [aegl] = X1 | X4             1               0.45750     TRUE
    ## [bfim] = X2 | X1             1               0.05606     TRUE
    ## [beik] = X2 | X3             1               0.10027     TRUE
    ## [befl] = X2 | X4             1               0.10698     TRUE
    ## [cfjm] = X3 | X1             1               0.03979     TRUE
    ## [cgjn] = X3 | X2             1               0.01160     TRUE
    ## [cfgl] = X3 | X4             1               0.03986     TRUE
    ## [dijm] = X4 | X1             1               0.02575     TRUE
    ## [dhjn] = X4 | X2             1               0.01207     TRUE
    ## [dhik] = X4 | X3             1               0.03361     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_part, digits = 2)
```

![](narrative_files/figure-markdown_github/var%20part-1.png)

WinterPCoAxis\_1 has the largest chunk (.4)

*look more at composition of pc axes!*
