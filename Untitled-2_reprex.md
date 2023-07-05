This template demonstrates many of the bells and whistles of the `reprex::reprex_document()` output format. The YAML sets many options to non-default values, such as using `#;-)` as the comment in front of output.

## Code style

Since `style` is `TRUE`, this difficult-to-read code (look at the `.Rmd` source file) will be restyled according to the Tidyverse style guide when it’s rendered. Whitespace rationing is not in effect!

``` r
x <- 1
```

<details style="margin-bottom:10px;">
<summary>
Standard output and standard error
</summary>

``` sh
-- nothing to show --
```

</details>
<details style="margin-bottom:10px;">
<summary>
Session info
</summary>

``` r
sessioninfo::session_info()
#;-) ─ Session info ───────────────────────────────────────────────────────────────
#;-)  setting  value
#;-)  version  R version 4.3.0 (2023-04-21 ucrt)
#;-)  os       Windows 11 x64 (build 23493)
#;-)  system   x86_64, mingw32
#;-)  ui       RTerm
#;-)  language (EN)
#;-)  collate  English_United States.utf8
#;-)  ctype    English_United States.utf8
#;-)  tz       Europe/Amsterdam
#;-)  date     2023-07-04
#;-)  pandoc   3.1.1 @ C:/Users/Eriba/AppData/Local/Pandoc/ (via rmarkdown)
#;-) 
#;-) ─ Packages ───────────────────────────────────────────────────────────────────
#;-)  package     * version date (UTC) lib source
#;-)  cli           3.6.1   2023-03-23 [1] CRAN (R 4.3.0)
#;-)  digest        0.6.31  2022-12-11 [1] CRAN (R 4.3.0)
#;-)  evaluate      0.21    2023-05-05 [1] CRAN (R 4.3.0)
#;-)  fastmap       1.1.1   2023-02-24 [1] CRAN (R 4.3.0)
#;-)  fs            1.6.2   2023-04-25 [1] CRAN (R 4.3.0)
#;-)  glue          1.6.2   2022-02-24 [1] CRAN (R 4.3.0)
#;-)  htmltools     0.5.5   2023-03-23 [1] CRAN (R 4.3.0)
#;-)  knitr         1.43    2023-05-25 [1] CRAN (R 4.3.0)
#;-)  lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.3.0)
#;-)  magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
#;-)  purrr         1.0.1   2023-01-10 [1] CRAN (R 4.3.0)
#;-)  R.cache       0.16.0  2022-07-21 [1] CRAN (R 4.3.0)
#;-)  R.methodsS3   1.8.2   2022-06-13 [1] CRAN (R 4.3.0)
#;-)  R.oo          1.25.0  2022-06-12 [1] CRAN (R 4.3.0)
#;-)  R.utils       2.12.2  2022-11-11 [1] CRAN (R 4.3.0)
#;-)  reprex        2.0.2   2022-08-17 [1] CRAN (R 4.3.0)
#;-)  rlang         1.1.1   2023-04-28 [1] CRAN (R 4.3.0)
#;-)  rmarkdown     2.22    2023-06-01 [1] CRAN (R 4.3.0)
#;-)  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
#;-)  styler        1.10.1  2023-06-05 [1] CRAN (R 4.3.1)
#;-)  vctrs         0.6.3   2023-06-14 [1] CRAN (R 4.3.1)
#;-)  withr         2.5.0   2022-03-03 [1] CRAN (R 4.3.0)
#;-)  xfun          0.39    2023-04-20 [1] CRAN (R 4.3.0)
#;-)  yaml          2.3.7   2023-01-23 [1] CRAN (R 4.3.0)
#;-) 
#;-)  [1] C:/Program Files/R/R-4.3.0/library
#;-) 
#;-) ──────────────────────────────────────────────────────────────────────────────
```

</details>
