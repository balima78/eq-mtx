library(tidyverse)
library(simK)

tab_res %>% count(model, SP_n) %>% as.data.frame()


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

#library(eq.mtx)
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


histoc::several(
  iteration.number = 10,
  df.donors = tar_read(dnrs),
  df.candidates = tar_read(cndts),
  df.abs = tar_read(antbs),
  algorithm = eqm,
  n = 0,
  seed.number = 123,
  check.validity = TRUE,
  q2 = 60,
  q3 = 80,
  uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.5))


tar_read(cndts)$dialysis |> quantile(probs=seq(0, 1, 0.1))

library(tidyverse)
library(ggridges)
tar_read(cndts) %>% ggplot(aes(cPRA)) + geom_histogram(binwidth=20)


  ggplot(aes(x = cPRA, y = urgent)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  )

tar_read(antbs) %>% filter(ID == 'K87') %>% arrange(abs) %>% as.data.frame()

tar_read(cndts) %>% filter(cPRA > 90,ID == 'K87')


tab_res <- res_modelos(age_avg) %>% .$data %>%
  left_join(res_modelos_extra(variavel = ageDiff_avg) %>% .$data) %>%
  left_join(res_modelos(mmHLA_avg) %>% .$data) %>%
  left_join(res_modelos(variavel = mmHLA0_n) %>% .$data) %>%
  left_join(res_modelos_extra(variavel = mmHLA02_n) %>% .$data) %>%
  left_join(res_modelos(variavel = mmHLA6_n) %>% .$data) %>%
  left_join(res_modelos(variavel = dialysis_avg) %>% .$data) %>%
  left_join(res_modelos(variavel = cPRA_avg) %>% .$data) %>%
  left_join(res_modelos(variavel = HI_n) %>% .$data) %>%
  left_join(res_modelos(variavel = SP_n) %>% .$data) %>%
  left_join(res_modelos_extra(variavel = txs_avg) %>% .$data)

tab_res %>%
  select(-it) %>%
  filter(model %in% c('Lima', 'ET', 'eqm0101', 'eqm_am_0101', 'eqm_et_0101')) %>%
  gtsummary::tbl_summary(by = model,
                         statistic = c(age_avg, ageDiff_avg, mmHLA_avg,
                                       mmHLA0_n, mmHLA02_n,
                                       mmHLA6_n, dialysis_avg,cPRA_avg,
                                       HI_n, SP_n, txs_avg
                                       ) ~ "{mean} (+/- {sd})",
                         digits = c(age_avg, ageDiff_avg, mmHLA_avg,
                                    mmHLA0_n, mmHLA02_n,
                                    mmHLA6_n, dialysis_avg,cPRA_avg,
                                    HI_n, SP_n, txs_avg) ~ 1,
                         label = list(age_avg ~ "recipients' age",
                                      ageDiff_avg ~ "age differences |donor-recipient|",
                                      mmHLA_avg ~ "number of HLA mm",
                                      mmHLA0_n ~ "recipients' with 0 HLA mm",
                                      mmHLA02_n ~ "recipients' with 0 to 2 HLA mm",
                                      mmHLA6_n ~ "recipients' with 6 HLA mm",
                                      dialysis_avg ~ "time on dialysis",
                                      cPRA_avg ~ "cPRA",
                                      HI_n ~ "Hipersensitized recipients",
                                      SP_n ~ "transplants with Senior Program",
                                      txs_avg ~ "TxScore"),
                         ) #%>%
  #gtsummary::add_p() %>%
  # gtsummary::add_overall()

tar_read(cndts) %>%
  filter(cPRA == 85)


lsr::cohensD(tar_read(res_et_extra)$ageDiff_avg, tar_read(res_lima_extra)$ageDiff_avg)
effsize::cohen.d(tar_read(res_et_extra)$ageDiff_avg, tar_read(res_lima_extra)$ageDiff_avg)


teste <- tar_read(res_et_extra) %>% select(-data_extra) %>% mutate(model = 'ET') %>%
  bind_rows(tar_read(res_lima_extra) %>% select(-data_extra) %>% mutate(model = 'Lima')) %>%
  mutate(mod = ifelse(model == 'ET', 1, 2))

library(psych)
res <- psych::cohen.d(as.data.frame(teste[c(-1,-5)]), 'mod')

res$p


tab_et_lima <- tar_read(res_et) %>%
  left_join(tar_read(res_et_extra)) %>%
  select(-it, -starts_with('data')) %>%
  mutate(model = 'ET') %>%
  bind_rows(tar_read(res_lima) %>%
              left_join(tar_read(res_lima_extra)) %>%
              select(-it, -starts_with('data')) %>%
              mutate(model = 'Lima')) %>%
  mutate(model = ifelse(model == 'ET', 1, 2))


res_et_lima <- psych::cohen.d(as.data.frame(tab_et_lima), 'model')

res_et_lima$p

x1<-rnorm(50,2,0.25)
b1<-boot::boot(x1,function(u,i) mean(u[i]),R=10)
boot::boot.ci(b1,type=c("norm","basic","perc"))


corr.fun <- function(data, idx)
{
  df <- data[idx, ]

  # Find the spearman correlation between
  # the 3rd and 4th columns of dataset
  c(cor(df[, 3], df[, 4], method = 'spearman'))
}


corr.fun(iris)
