lambda0 <- 0.5
compl_rate <- 1e5
lambda1 <- 0.5
mu0 <- 0.0
mu1 <- 0.0


v2 <- pbdABC::test_prior(num_samples = 100000,
                         means = c(1 / lambda0, 1 / lambda1, 10000, 10000, compl_rate),
                         use_inv = FALSE)

v1 <- pbdABC::test_prior(num_samples = 100000,
                         means = c(1 / lambda0, 1 / lambda1, 10000, 10000, compl_rate),
                         use_inv = TRUE)



require(tidyverse)
colnames(v1) <- c("l0", "l1", "m0", "m1", "compl_rate")
v1 <- as_tibble(v1)

colnames(v2) <- c("l0", "l1", "m0", "m1", "compl_rate")
v2 <- as_tibble(v2)

v1$type <- "inv"
v2$type <- "normal"

v3 <- rbind(v1, v2)



v3 %>%
  gather(key = "statistic", value = "val", -c(type)) %>%
  ggplot(aes(x = val, col = type)) +
    geom_density() +
    facet_wrap(~statistic, scales = "free") +
    scale_x_log10()

v3 %>%
  ggplot(aes(x = compl_rate, col = type)) +
    geom_density()


colMeans(v1[, 1:5])
colMeans(v2[, 1:5])
