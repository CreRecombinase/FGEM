LabMeeting
========================================================
author: Nicholas Knoblauch
date: 10/10/2016
autosize: true

Features
========================================================

-Gene Ontology
  - 5898 features
-ExAC
  - 19 features
-Brainspan Expression
  - 32 features (16 regions x Sex, averaged over developmental stages)


Since last time
========================================================

  -Fixed a bug wherein the likelihood ratio was computed incorrectly
  -Added new features (more ExAC)

Significant Features
========================================================


```
  feature          Beta  Intercept     Chisq         pval prior_mean class
1   pNull -7.784540e+05  -2.221050 192.67979 8.268230e-44 0.01418624  ExAC
2    pRec -8.643961e+01  -2.011946 175.10086 5.690921e-40 0.01424332  ExAC
3   lof_z  4.566989e-01  -5.483995 150.45611 1.378051e-34 0.01734152  ExAC
4     pLI  7.104705e+01 -72.977466 148.88915 3.032298e-34 0.01304573  ExAC
5   mis_z  6.037565e-01  -5.138497  87.56902 8.138964e-21 0.01579879  ExAC
6 exp_lof  9.877854e-03  -4.074519  36.26839 1.719296e-09 0.02034761  ExAC
          bh.p
1 1.570964e-42
2 5.406375e-39
3 8.727655e-34
4 1.440342e-33
5 3.092806e-20
6 5.444437e-09
```


```r
group_by(results,class) %>% summarise(n_sig=sum(bh.p<0.05))
```

```
# A tibble: 3 x 2
       class n_sig
       <chr> <int>
1       ExAC    16
2 Expression     0
3         GO     0
```

Overall Model
========================================================

Jointly estimate the prior, using significant features



```
Error in inner_join(datadf, annodf) : object 'datadf' not found
```
