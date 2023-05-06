library(tidyverse)
library(simK)


cds <- candidates_df(n = 2000) %>%
  select(ID, bg, A1, A2, B1, B2, DR1, DR2, age, dialysis, cPRA, urgent)
dns <- donors_df(n = 50)
ant <- Abs_df(candidates =cds)


library(histoc)

library(tictoc)
tic()
several(iteration.number = 10, df.donors = donors_df(n=50), df.candidates = cds,
        df.abs = Abs_df(cds), algorithm = lima, n = 0, seed.number = 123,
        check.validity = TRUE)
toc()


tic()
res_pairs <- donor_recipient_pairs(df.donors = dns,
                                  df.candidates = cds,
                                  df.abs = ant,
                                  algorithm = lima,
                                  n = 0,
                                  check.validity = TRUE
                                  , q2 = 60
                                  , q3 = 100
                                  , cPRA1 = 50
                                  , cPRA2 = 85)
toc()

lima(iso = TRUE
     , dABO = "A"
     , dA = c("29","80"), dB = c("14","39"), dDR = c("3","11")
     , donor.age = 36
     , df.abs = ant
     , data = cds
     , n = 0
     , q2 = 60
     , q3 = 100
     , cPRA1 = 50
     , cPRA2 = 85
     , check.validity = TRUE)


dns$nrow <- 1:nrow(dns)



purrr::map(all.statistics, ~table(.x$SP)) |> bind_rows()

res_pairs[[1]] %>%
  rowwise() %>%
  mutate(txsc = txscore(recipient.age = age,
                        recipient.dialysis = dialysis,
                        donor.age = donor_age,
                        mmHLA_A = mmA,
                        mmHLA_B = mmB,
                        mmHLA_DR = mmDR
                        )$prob5y)

txscore(recipient.age = 59,
        recipient.dialysis = 33,
        donor.age = 53,
        mmHLA_A = 1,
        mmHLA_B = 2,
        mmHLA_DR = 0)



lima(iso = TRUE
                 , dABO = "O"
                 , dA = c("1","2"), dB = c("15","44"), dDR = c("1","4")
                 , donor.age = 60
                 , df.abs = ant
                 , data = cds
                 , n = 0
                 , q2 = 60
                 , q3 = 100
                 , cPRA1 = 50
                 , cPRA2 = 85
                 , check.validity = TRUE)

library(eq.mtx)
library(histoc)
# teste de execução do algoritmo 'eqm' 1 vez
eqm(iso = TRUE, dABO = "O", dA = c("1", "2"),
    dB = c("15", "44"), dDR = c("1", "4"),
    donor.age = 60, df.abs = ant, data = cds, n = 2,
    q2 = 60, q3 = 80, uj.matx = uj_matx(), check.validity = TRUE)
# teste de execução do algoritmo 'lima' 1 vez
lima(iso = TRUE, dABO = "O", dA = c("1", "2"),
     dB = c("15", "44"), dDR = c("1", "4"),
     donor.age = 60, df.abs = ant, data = cds, n = 2,
     q2 = 60, q3 = 80
     , cPRA1 = 50
     , cPRA2 = 85
     , check.validity = TRUE)


several(
  iteration.number = 100,
  df.donors = tar_read(dnrs),
  df.candidates = tar_read(cndts),
  df.abs = tar_read(antbs),
  algorithm = eqm,
  n = 0,
  seed.number = 123,
  check.validity = TRUE,
  q2 = 60,
  q3 = 80,
  uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.1))


tar_read(cndts)$dialysis |> quantile(probs=seq(0, 1, 0.1))
