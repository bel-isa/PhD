### Number of clusters for clustering condition averages silhouette coeff ###

> num.clus(data3)
2  avg:  0.06 clus:  -0.01 0.25 
3  avg:  0.07 clus:  0.02 0.22 0 
4  avg:  0.09 clus:  0.06 0.18 0 0 
5  avg:  0.11 clus:  0 0.15 0.12 0 0 
6  avg:  0.13 clus:  0 0.23 0.08 0.14 0 0 
7  avg:  0.11 clus:  0 0.2 0 0.13 0 0 0 
8  avg:  0.09 clus:  0 0.19 0 0 0 0.12 0 0 
9  avg:  0.06 clus:  0 0.18 0 0 0 0 0 0.05 0 
10  avg:  0.07 clus:  0 0.17 0 0 0 0 0 0.1 0 0 

> num.clus(t(data3))
2  avg:  0.18 clus:  0.12 0.22 
3  avg:  0.16 clus:  0.1 0.21 0.18 
4  avg:  0.16 clus:  0.1 0.37 0.17 0.09 
5  avg:  0.17 clus:  0.22 0.22 0.02 0.09 0.26 
6  avg:  0.2 clus:  0.21 0.22 0.02 0.24 0.15 0.35 
7  avg:  0.26 clus:  0.9 0.35 0.06 0.14 0.23 0.12 0.35 
8  avg:  0.22 clus:  0.9 0.34 0.05 0.14 0.27 0.13 0.33 0.08 
9  avg:  0.22 clus:  0.9 0.34 0.13 0.12 0.25 0 0.12 0.33 0.08 
10  avg:  0.23 clus:  0.9 0.33 0.3 0.12 0.24 0 0 0.12 0.33 0.08

> Samp.clu1 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==1,]
> rownames(Samp.clu1)
[1] "Syn exp1 HC" "Syn exp3 HC"
> Samp.clu2 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==2,]
> rownames(Samp.clu2)
[1] "Syn exp1 3h LC"    "Syn exp3 3h LC"    "Ara exp4 light HC" "Ara exp7 4h HC"   
> Samp.clu3 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==3,]
> rownames(Samp.clu3)
[1] "Syn exp1 24h LC"      "Syn exp2 24h LC"      "Syn exp3 24h LC"     
[4] "Ara exp4 light 1d LC" "Ara exp7 6h NC"      
> Samp.clu4 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==4,]
> rownames(Samp.clu4)
[1] "Syn exp2 HC"    "Ara exp6 2h VC" "Ara exp6 4h VC" "Ara exp6 6h VC"
> Samp.clu5 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==5,]
> rownames(Samp.clu5)
[1] "Syn exp2 3h LC"       "Ara exp4 light 5d LC" "Ara exp6 6h NC"      
[4] "Ara exp6 8h NC"       "Ara exp6 10h NC"      "Ara exp7 4h NC"      
> Samp.clu6 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==6,]
> rownames(Samp.clu6)
[1] "Ara exp4 light 3d LC" "Ara exp6 2h NC"       "Ara exp6 4h NC"      
[4] "Ara exp7 2h NC"       "Ara exp7 2h HC"      
> Samp.clu7 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==7,]
> rownames(Samp.clu7)
[1] "Ara exp6 8h VC"  "Ara exp6 10h VC" "Ara exp7 8h NC"  "Ara exp7 10h NC"
[5] "Ara exp7 6h HC"  "Ara exp7 8h HC"  "Ara exp7 10h HC"

### Number of clusters for clustering condition averages affinity propagation ###

APResult object

Number of samples     =  33 
Number of iterations  =  153 
Input preference      =  -31.34694 
Sum of similarities   =  -282.3255 
Sum of preferences    =  -219.4286 
Net similarity        =  -501.7541 
Number of clusters    =  7 

Exemplars:
  Syn exp1 HC Syn exp1 3h LC Syn exp3 24h LC Ara exp6 4h NC Ara exp6 10h NC 
Ara exp6 2h VC Ara exp6 10h VC
Clusters:
  Cluster 1, exemplar Syn exp1 HC:
  Syn exp1 HC Syn exp3 HC
Cluster 2, exemplar Syn exp1 3h LC:
  Syn exp1 3h LC Syn exp3 3h LC Ara exp4 light HC Ara exp7 4h HC
Cluster 3, exemplar Syn exp3 24h LC:
  Syn exp1 24h LC Syn exp2 24h LC Syn exp3 24h LC Ara exp4 light 1d LC Ara exp7 6h NC
Cluster 4, exemplar Ara exp6 4h NC:
  Ara exp4 light 3d LC Ara exp6 2h NC Ara exp6 4h NC Ara exp7 2h NC Ara exp7 2h HC
Cluster 5, exemplar Ara exp6 10h NC:
  Syn exp2 3h LC Ara exp4 light 5d LC Ara exp6 6h NC Ara exp6 8h NC Ara exp6 10h NC 
Ara exp7 4h NC
Cluster 6, exemplar Ara exp6 2h VC:
  Syn exp2 HC Ara exp6 2h VC Ara exp6 4h VC Ara exp6 6h VC
Cluster 7, exemplar Ara exp6 10h VC:
  Ara exp6 8h VC Ara exp6 10h VC Ara exp7 8h NC Ara exp7 10h NC Ara exp7 6h HC 
Ara exp7 8h HC Ara exp7 10h HC

### CV correlation: regression analysis ###

> summary(reg)

Call:
  lm(formula = CV.S ~ CV.A - 1)

Residuals:
  Min       1Q   Median       3Q      Max 
-0.36659 -0.13521  0.02736  0.21851  0.92517 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
CV.A   0.6437     0.1134   5.674 2.52e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3358 on 33 degrees of freedom
Multiple R-squared:  0.4938,  Adjusted R-squared:  0.4785 
F-statistic:  32.2 on 1 and 33 DF,  p-value: 2.519e-06

> summary(impr)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-1.1290 -0.4065  0.0826  0.2689  0.6531  2.7650 

> summary(reg1)

Call:
  lm(formula = CV.A[data4] ~ CV.S[data4] - 1)

Residuals:
  Min       1Q   Median       3Q      Max 
-0.36611 -0.05088  0.09488  0.27211  0.34412 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
CV.S[data4]   0.8664     0.1029    8.42 1.76e-08 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2112 on 23 degrees of freedom
Multiple R-squared:  0.755,  Adjusted R-squared:  0.7444 
F-statistic: 70.89 on 1 and 23 DF,  p-value: 1.765e-08

> summary(reg2)

Call:
  lm(formula = CV.A[data5] ~ CV.S[data5] - 1)

Residuals:
  Min       1Q   Median       3Q      Max 
-0.25740 -0.04956  0.05442  0.11575  0.18252 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
CV.S[data5]  0.79459    0.06817   11.66 1.36e-08 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1223 on 14 degrees of freedom
Multiple R-squared:  0.9066,  Adjusted R-squared:  0.8999 
F-statistic: 135.9 on 1 and 14 DF,  p-value: 1.357e-08

### Number of clusters for clustering rv_matrix affinity propagation ###

> d.apclus

APResult object

Number of samples     =  45 
Number of iterations  =  152 
Input preference      =  -3.462455 
Sum of similarities   =  -38.57479 
Sum of preferences    =  -38.087 
Net similarity        =  -76.66179 
Number of clusters    =  11 

Exemplars:
  Syn exp1 HC Syn exp1 24h LC Ara exp4 light HC Ara exp5 dark HC Ara exp6 10h NC 
Ara exp6 4h VC Ara exp6 8h VC Ara exp7 10h NC Ara exp7 2h HC Ara exp7 12h HC 
Ara exp7 24h HC
Clusters:
  Cluster 1, exemplar Syn exp1 HC:
  Syn exp1 HC Syn exp3 HC
Cluster 2, exemplar Syn exp1 24h LC:
  Syn exp1 24h LC Syn exp2 HC Syn exp3 24h LC
Cluster 3, exemplar Ara exp4 light HC:
  Syn exp1 3h LC Syn exp3 3h LC Ara exp4 light HC Ara exp5 dark 1d LC
Cluster 4, exemplar Ara exp5 dark HC:
  Ara exp5 dark HC Ara exp5 dark 5d LC Ara exp6 2h NC Ara exp6 6h VC Ara exp7 4h HC
Cluster 5, exemplar Ara exp6 10h NC:
  Syn exp2 3h LC Ara exp6 6h NC Ara exp6 8h NC Ara exp6 10h NC
Cluster 6, exemplar Ara exp6 4h VC:
  Syn exp2 24h LC Ara exp4 light 3d LC Ara exp4 light 5d LC Ara exp6 12h NC 
Ara exp6 4h VC
Cluster 7, exemplar Ara exp6 8h VC:
  Ara exp5 3d LC Ara exp6 2h VC Ara exp6 8h VC Ara exp6 10h VC Ara exp7 6h HC 
Ara exp7 8h HC
Cluster 8, exemplar Ara exp7 10h NC:
  Ara exp7 8h NC Ara exp7 10h NC Ara exp7 10h HC
Cluster 9, exemplar Ara exp7 2h HC:
  Ara exp4 light 1d LC Ara exp7 2h NC Ara exp7 2h HC
Cluster 10, exemplar Ara exp7 12h HC:
  Ara exp6 4h NC Ara exp6 12h VC Ara exp7 4h NC Ara exp7 6h NC Ara exp7 12h NC 
Ara exp7 12h HC
Cluster 11, exemplar Ara exp7 24h HC:
  Ara exp6 24h NC Ara exp6 24h VC Ara exp7 24h NC Ara exp7 24h HC

### Rv_matrix1 <- no dark samples ###

> num.clus(rv_matrix1)
2  avg:  0.18 clus:  0.12 0.41 
3  avg:  0.16 clus:  0.09 0.09 0.4 
4  avg:  0.19 clus:  0.06 0.19 0.22 0.38 
5  avg:  0.22 clus:  0.38 0.18 0.2 0.1 0.37 
6  avg:  0.26 clus:  0.95 0.16 0.31 0.18 0.1 0.37 
7  avg:  0.27 clus:  0.95 0.41 0.47 0.17 0.01 0.07 0.37 
8  avg:  0.26 clus:  0.95 0.41 0.44 0.01 0.11 0.13 0.18 0.34 
9  avg:  0.26 clus:  0.95 0.41 0.44 0.01 0.11 0.13 0.18 0.38 0.24 
10  avg:  0.26 clus:  0.95 0.41 0.45 0.22 -0.03 0.08 0.24 0.09 0.38 0.24 

> cl1 <- PAM1[PAM1[,ncol(PAM1)]==1,]
> rownames(cl1)
[1] "Syn exp1 HC" "Syn exp3 HC"
> cl2 <- PAM1[PAM1[,ncol(PAM1)]==2,]
> rownames(cl2)
[1] "Syn exp1 3h LC"    "Syn exp3 3h LC"    "Ara exp4 light HC"
> cl3 <- PAM1[PAM1[,ncol(PAM1)]==3,]
> rownames(cl3)
[1] "Syn exp1 24h LC" "Syn exp2 HC"     "Syn exp3 24h LC"
> cl4 <- PAM1[PAM1[,ncol(PAM1)]==4,]
> rownames(cl4)
[1] "Syn exp2 3h LC"       "Ara exp4 light 5d LC" "Ara exp6 6h NC"      
[4] "Ara exp6 8h NC"       "Ara exp6 10h NC"      "Ara exp6 2h VC"      
[7] "Ara exp6 4h VC"      
> cl5 <- PAM1[PAM1[,ncol(PAM1)]==5,]
> rownames(cl5)
[1] "Syn exp2 24h LC"      "Ara exp4 light 1d LC" "Ara exp6 6h VC"      
[4] "Ara exp7 4h HC"      
> cl6 <- PAM1[PAM1[,ncol(PAM1)]==6,]
> rownames(cl6)
[1] "Ara exp4 light 3d LC" "Ara exp6 2h NC"       "Ara exp6 4h NC"      
[4] "Ara exp7 2h NC"       "Ara exp7 4h NC"       "Ara exp7 6h NC"      
[7] "Ara exp7 2h HC"      
> cl7 <- PAM1[PAM1[,ncol(PAM1)]==7,]
> rownames(cl7)
[1] "Ara exp6 8h VC"  "Ara exp6 10h VC" "Ara exp7 8h NC"  "Ara exp7 10h NC"
[5] "Ara exp7 6h HC"  "Ara exp7 8h HC"  "Ara exp7 10h HC"

> d.apclus

APResult object

Number of samples     =  33 
Number of iterations  =  193 
Input preference      =  -3.145012 
Sum of similarities   =  -31.17184 
Sum of preferences    =  -22.01509 
Net similarity        =  -53.18693 
Number of clusters    =  7 

Exemplars:
  Syn exp1 HC Syn exp1 3h LC Syn exp2 HC Ara exp6 4h NC Ara exp6 6h NC Ara exp6 10h VC 
Ara exp7 4h HC
Clusters:
  Cluster 1, exemplar Syn exp1 HC:
  Syn exp1 HC Syn exp3 HC
Cluster 2, exemplar Syn exp1 3h LC:
  Syn exp1 3h LC Syn exp3 3h LC Ara exp4 light HC
Cluster 3, exemplar Syn exp2 HC:
  Syn exp1 24h LC Syn exp2 HC Syn exp2 3h LC Syn exp2 24h LC Syn exp3 24h LC 
Ara exp6 4h VC
Cluster 4, exemplar Ara exp6 4h NC:
  Ara exp4 light 3d LC Ara exp6 2h NC Ara exp6 4h NC Ara exp7 2h NC Ara exp7 6h NC 
Ara exp7 2h HC
Cluster 5, exemplar Ara exp6 6h NC:
  Ara exp4 light 5d LC Ara exp6 6h NC Ara exp6 8h NC Ara exp6 10h NC Ara exp6 2h VC 
Ara exp7 4h NC
Cluster 6, exemplar Ara exp6 10h VC:
  Ara exp6 8h VC Ara exp6 10h VC Ara exp7 8h NC Ara exp7 10h NC Ara exp7 6h HC 
Ara exp7 8h HC Ara exp7 10h HC
Cluster 7, exemplar Ara exp7 4h HC:
  Ara exp4 light 1d LC Ara exp6 6h VC Ara exp7 4h HC