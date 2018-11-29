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
  
  Summer_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.3395983 0.21255384 0.2125538
    ##  [2,]   1.9042858 0.17300545 0.3855593
    ##  [3,]   1.2857635 0.11681234 0.5023716
    ##  [4,]   1.1236910 0.10208796 0.6044596
    ##  [5,]   0.7959522 0.07231271 0.6767723
    ##  [6,]   0.5829654 0.05296274 0.7297350
    ##  [7,]   0.4955673 0.04502257 0.7747576
    ##  [8,]   0.4403899 0.04000967 0.8147673
    ##  [9,]   0.4093077 0.03718584 0.8519531
    ## [10,]   0.3597252 0.03268124 0.8846344
    ## [11,]   0.2711839 0.02463721 0.9092716
    ## [12,]   0.2573015 0.02337598 0.9326476
    ## [13,]   0.1754937 0.01594370 0.9485913
    ## [14,]   0.1540241 0.01399318 0.9625844
    ## [15,]   0.1257192 0.01142166 0.9740061

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
  
  winter_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.9520570 0.26441785 0.2644179
    ##  [2,]   1.4121916 0.12649101 0.3909089
    ##  [3,]   1.1515115 0.10314171 0.4940506
    ##  [4,]   0.9806955 0.08784160 0.5818922
    ##  [5,]   0.8641039 0.07739841 0.6592906
    ##  [6,]   0.5451243 0.04882717 0.7081177
    ##  [7,]   0.5196633 0.04654661 0.7546644
    ##  [8,]   0.4824972 0.04321762 0.7978820
    ##  [9,]   0.4492552 0.04024011 0.8381221
    ## [10,]   0.4191044 0.03753948 0.8756616
    ## [11,]   0.3301053 0.02956777 0.9052293
    ## [12,]   0.2762754 0.02474618 0.9299755
    ## [13,]   0.2119215 0.01898196 0.9489575
    ## [14,]   0.1506384 0.01349279 0.9624503
    ## [15,]   0.1455869 0.01304032 0.9754906

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](narrative_files/figure-markdown_github/winter%20pcoa-1.png)

``` r
ordiplot(scores(winter_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Winter plant PCoA")
```

    ## species scores not available

``` r
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(winter_pcoa$points[, 1:2], winter_plants)
text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
```

![](narrative_files/figure-markdown_github/winter%20pcoa-2.png)

There seems to be an inflection point in the scree plot around axis 2 or 3. Since stopping at 2 would only capture 39% of variation, going to go for 3.

Look at axis composition:

``` r
winter_pc <- wascores(winter_pcoa$points[, 1:3], winter_plants)
winter_pc
```

    ##                   [,1]          [,2]          [,3]
    ## year       0.001314352 -0.0003197344 -3.911284e-05
    ## ambr.arte  0.347620933  0.0932001182  8.205587e-02
    ## amsi.tess  0.272179935  0.0488011383  1.840448e-02
    ## astr.allo -0.302042493  0.0563240035 -4.348513e-02
    ## astr.nutt  0.196783526  0.1216213400  9.559205e-02
    ## bail.mult -0.349730258  0.0318042729 -2.816683e-03
    ## caly.wrig -0.150063910  0.1281348827  2.534284e-01
    ## chae.stev  0.181462957  0.1576475321  3.968224e-02
    ## chen.frem  0.174015633  0.0614183015  1.595803e-01
    ## chor.tene  0.032819160  0.2142319042  3.119201e-01
    ## cryp.cras -0.246254363  0.0243128057  6.354981e-03
    ## cryp.micr -0.269819330  0.0827256901  2.017683e-02
    ## dale.brac -0.249153532  0.1455530740 -1.910574e-01
    ## desc.obtu  0.335960238 -0.0216744148  1.205310e-01
    ## desc.pinn -0.215060705  0.0627305380 -1.569288e-03
    ## dith.wisl -0.067638093 -0.0015408730  3.257677e-01
    ## eria.diff -0.310584946 -0.0087005003 -3.233826e-02
    ## erig.dive -0.351736851  0.0293889879 -5.656956e-02
    ## erio.aber -0.286000486  0.0781337951 -4.249522e-02
    ## erod.cicu  0.329081955  0.1694267594 -1.794562e-01
    ## erod.texa -0.211061332  0.1524584540  1.077483e-02
    ## esch.mexi -0.387655101  0.1050356672  3.774531e-02
    ## gili.sinu -0.335395608  0.0963853845 -2.423223e-03
    ## hapl.grac -0.419917554  0.1104684240 -1.334116e-01
    ## laen.coul  0.321967405 -0.1595238544  1.667011e-01
    ## lapp.redo  0.289865941  0.1170394765  1.215048e-01
    ## lepi.lasi -0.033133773  0.0448629565  8.973910e-02
    ## lesq.gord  0.242210203  0.0664921068  4.515143e-02
    ## lina.bige  0.252852158  0.1381147402  3.824457e-02
    ## lupi.brev  0.280237438 -0.0947265093  2.252263e-01
    ## lupi.conc -0.244075704  0.1321220311  3.496476e-02
    ## mala.fend -0.338415787  0.0326247952 -2.791676e-02
    ## micr.lene -0.137252626  0.1095524615  1.329083e-01
    ## oeno.prim  0.144649512  0.1065147716  1.599352e-01
    ## pani.hirt  0.154727172  0.0033173117 -1.206497e-01
    ## pect.recu  0.340758521 -0.0175003860  1.457291e-01
    ## phac.ariz  0.016187050  0.0888375434  1.376843e-01
    ## plag.ariz -0.129963420  0.1448658978  2.847516e-01
    ## plan.purs -0.126021019  0.1388997393  3.780076e-02
    ## schi.barb  0.149839276  0.1744337745  2.165952e-01
    ## sisy.irio  0.345676373  0.1505196255 -3.373036e-01
    ## sper.echi -0.211737492  0.1237111798  7.744768e-02
    ## step.exig -0.343451856  0.1039092616  1.075226e-02
    ## vulp.octo -0.245704194  0.1418332636  5.945074e-02

``` r
  rm(list=ls())
```

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
    ## RDA3      1 0.009484  4.0851  0.026 *  
    ## RDA4      1 0.004300  1.8520  0.556    
    ## RDA5      1 0.002483  1.0697  0.898    
    ## RDA6      1 0.000846  0.3645  0.999    
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

### Years - rodents comparison

``` r
plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
abline(v = 1990, col = 'red')
```

![](narrative_files/figure-markdown_github/plot%20winter%20pcoa1%20v%20year-1.png)
