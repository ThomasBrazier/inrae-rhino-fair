
package_list = c("devtools",
                 "MASS",
                 "fitdistrplus",
                 "gamlss.dist",
                 "ggplot2",
                 "ggpubr",
                 "gnorm",
                 "adegenet",
                 "coda",
                 "geosphere",
                 "hierfstat",
                 "ade4",
                 "fitdistrplus",
                 "moments",
                 "extraDistr",
                 "fdrtool",
                 "sp",
                 "distr")

for (p in 1:length(package_list)) {
  if (!(package_list[p] %in% installed.packages())) {
    install.packages(package_list[p])
  }
  library(package_list[p], character.only = TRUE)
}

