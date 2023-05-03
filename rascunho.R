library(tidyverse)
library(simK)


cds <- simK::candidates_df(n = 2000) %>%
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

tic()
all.statistics <- list()

for (. in 1:100){
  used.candidates <- NULL
  current.iteration.statistics <- NULL
  shuffled_donors <- sample(dns$nrow)

  for (j in 1:length(shuffled_donors)){
    tmp <- res_pairs[[shuffled_donors[j]]][
      !ID %in% used.candidates,][
        1:2,]

    current.iteration.statistics <- data.table::rbindlist(list(current.iteration.statistics,
                                                               tmp))

    used.candidates <- c(used.candidates, tmp$ID)
  }

  all.statistics <- append(all.statistics, list(current.iteration.statistics))
}
toc()


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
library(histoc)
data("ant")
eq.mtx::eqm(
iso = iso,
dABO = dABO,
dA = dA,
dB = dB,
dDR = dDR,
donor.age = donor.age,
df.abs = ant,
data = candidates,
n = 6,
q2 = q2,
q3 = q3,
uj.matx = eq.mtx::uj_matx()
)

