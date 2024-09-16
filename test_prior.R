lambda0 <- 0.5
compl_rate <- 1
lambda1 <- 0.5
mu0 <- 0.0
mu1 <- 0.0

v1 <- pbdABC::test_prior(num_samples = 1000,
                         means = c(1 / lambda0, 1 / lambda1, 10000, 10000, 1 / compl_rate),
                         use_inv = FALSE)

v2 <- pbdABC::test_prior(num_samples = 1000,
                         means = c(1 / lambda0, 1 / lambda1, 10000, 10000, 1 / compl_rate),
                         use_inv = TRUE)





require(tidyverse)
colnames(v1) <- c("compl_rate", "density")
v1 <- as_tibble(v1)
v1$dens2 <- dexp(x = v1$compl_rate, rate = compl_rate)


colnames(v2) <- c("compl_rate", "density")
v2 <- as_tibble(v2)
v2$den2 <- dexp(x = 1 / v2$compl_rate, rate = compl_rate)

v1$type <- "inv"
v2$type <- "normal"

v3 <- rbind(v1, v2)



v3 %>%
  gather(key = "statistic", value = "val", -c(type, density)) %>%
  ggplot(aes(x = val, col = type)) +
    geom_density() +
    facet_wrap(~statistic, scales = "free") +
    scale_x_log10()

ggplot(v1, aes(x = compl_rate, y = density)) +
    geom_point() +
    scale_x_log10() +
    ggtitle("normal")

ggplot(v2, aes(x = compl_rate, y = density)) +
  geom_point() +
  scale_x_log10() +
  ggtitle("inverse")
