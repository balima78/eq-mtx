library(eq.mtx)
library(histoc)

library(tictoc)
tic()
# executa todas as possibilidades para 'eqm'
res_eqm <- donor_recipient_pairs_v2(df.donors = dns,
                                    df.candidates = cds,
                                    df.abs = ant,
                                    algorithm = eqm,
                                    n = 0,
                                    q2 = 60, q3 = 100,
                                    check.validity = FALSE)
# executa todas as possibilidades para 'lima'
res_lima <- donor_recipient_pairs_v2(df.donors = dns,
                                     df.candidates = cds,
                                     df.abs = ant,
                                     algorithm = lima,
                                     n = 0,
                                     q2 = 60, q3 = 100,
                                     check.validity = FALSE)
toc()


dns$nrow <- 1:nrow(dns)
# resultados multi iterações para 'eqm'
tic()
all.statistics.eqm <- list()
set.seed(1)
for (. in 1:100){
  used.candidates <- NULL
  current.iteration.statistics <- NULL
  shuffled_donors <- sample(dns$nrow)

  for (j in 1:length(shuffled_donors)){
    tmp <- res_eqm[[shuffled_donors[j]]][
      !ID %in% used.candidates,][
        1:2,]

    current.iteration.statistics <- data.table::rbindlist(list(current.iteration.statistics,
                                                               tmp))

    used.candidates <- c(used.candidates, tmp$ID)
  }

  all.statistics.eqm <- append(all.statistics.eqm, list(current.iteration.statistics))
}
toc()

# resultados multi iterações para 'lima'
tic()
all.statistics.lima <- list()
set.seed(1)
for (. in 1:100){
  used.candidates <- NULL
  current.iteration.statistics <- NULL
  shuffled_donors <- sample(dns$nrow)

  for (j in 1:length(shuffled_donors)){
    tmp <- res_lima[[shuffled_donors[j]]][
      !ID %in% used.candidates,][
        1:2,]

    current.iteration.statistics <- data.table::rbindlist(list(current.iteration.statistics,
                                                               tmp))

    used.candidates <- c(used.candidates, tmp$ID)
  }

  all.statistics.lima <- append(all.statistics.lima, list(current.iteration.statistics))
}
toc()



bind_rows(all.statistics.eqm, .id = 'it') %>%
  # as_tibble() %>%
  # nest(data = !it) %>%
  group_by(it) %>%
  nest() %>%
  ungroup() %>%
  mutate(age_eqm = map_dbl(data, ~mean(.x$age)),
         mmHLA_eqm = map_dbl(data, ~mean(.x$mmHLA)),
         mmHLA0_eqm = map_dbl(data, ~sum(.x$mmHLA == 0)),
         mmHLA1_eqm = map_dbl(data, ~sum(.x$mmHLA == 1)),
         mmHLA2_eqm = map_dbl(data, ~sum(.x$mmHLA == 2)),
         mmHLA3_eqm = map_dbl(data, ~sum(.x$mmHLA == 3)),
         mmHLA4_eqm = map_dbl(data, ~sum(.x$mmHLA == 4)),
         mmHLA5_eqm = map_dbl(data, ~sum(.x$mmHLA == 5)),
         mmHLA6_eqm = map_dbl(data, ~sum(.x$mmHLA == 6)),
         dialysis_eqm = map_dbl(data, ~mean(.x$dialysis)),
         cPRA_eqm = map_dbl(data, ~mean(.x$cPRA)),
         HI_eqm = map_dbl(data, ~sum(.x$HI)),
         SP_eqm = map_dbl(data, ~sum(.x$SP))
         ) %>%
  rename(data.eqm = data)


bind_rows(all.statistics.lima, .id = 'it') %>%
  group_by(it) %>%
  nest() %>%
  ungroup() %>%
  mutate(age_lima = map_dbl(data, ~mean(.x$age)),
         mmHLA_lima = map_dbl(data, ~mean(.x$mmHLA)),
         mmHLA0_lima = map_dbl(data, ~sum(.x$mmHLA == 0)),
         mmHLA1_lima = map_dbl(data, ~sum(.x$mmHLA == 1)),
         mmHLA2_lima = map_dbl(data, ~sum(.x$mmHLA == 2)),
         mmHLA3_lima = map_dbl(data, ~sum(.x$mmHLA == 3)),
         mmHLA4_lima = map_dbl(data, ~sum(.x$mmHLA == 4)),
         mmHLA5_lima = map_dbl(data, ~sum(.x$mmHLA == 5)),
         mmHLA6_lima = map_dbl(data, ~sum(.x$mmHLA == 6)),
         dialysis_lima = map_dbl(data, ~mean(.x$dialysis)),
         cPRA_lima = map_dbl(data, ~mean(.x$cPRA)),
         HI_lima = map_dbl(data, ~sum(.x$HI))
  ) %>%
  rename(data.lima = data)



# ID bg A1 A2 B1 B2 DR1 DR2 mmA mmB mmDR mmHLA age donor_age dialysis cPRA HI ptsEQM SP urgent
#
# AM = ifelse(SP == 0 & cPRA >= 85, 1, 0)
