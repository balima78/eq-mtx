library(targets)
# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(data_summary) to view the results.
library(tarchetypes)

# run in parallel
library(future)
library(future.callr)
plan(callr)

# funçoes adhoc
dados_extra <- function(dados){
  purrr::map(dados,
      ~.x %>%
        dplyr::rowwise() %>%
        dplyr::mutate(txs = txscore(recipient.age = age,
                             recipient.race = "White", recipient.causeESRD = "Other",
                             recipient.dialysis = dialysis,
                             recipient.diabetes = FALSE, recipient.coronary = FALSE, recipient.albumin = 1.5, recipient.hemoglobin = 10,
                             donor.age = donor_age,
                             donor.diabetes = "Absence", donor.ECD = FALSE,
                             mmHLA_A = mmA, mmHLA_B = mmB, mmHLA_DR = mmDR)$prob5y,
               ageDiff = abs(donor_age - age),
               mmHLA02 = mmHLA <=2) %>%
        dplyr::ungroup())
}

# Set target-specific options such as packages:
tar_option_set(packages = c("histoc","simK"),
               memory = 'transient',
               garbage_collection = TRUE)

# End this file with a list of target objects.
list(
  # # dados v1 #########
  # tar_target(
  #   dnrs,
  #   simK::donors_df_v1(n = 280,
  #                   replace = TRUE,
  #                   origin = 'PT',
  #                   probs = c(0.4658, 0.0343, 0.077, 0.4229),
  #                   lower=18, upper=75,
  #                   mean.age = 60, sd.age = 12,
  #                   uk = FALSE,
  #                   seed.number = 3)
  # ),
  # tar_target(
  #   cndts,
  #   simK::candidates_df_v1(n = 2000,
  #                       replace = TRUE,
  #                       origin = 'PT',
  #                       probs.abo = c(0.43, 0.03, 0.08, 0.46),
  #                       probs.cpra = c(0.7, 0.1, 0.1, 0.1),
  #                       lower=18, upper=75,
  #                       mean.age = 45, sd.age = 15,
  #                       prob.dm = 0.12,
  #                       prob.urgent = 0,
  #                       uk = FALSE,
  #                       seed.number = 3) |>
  #     dplyr::select(ID, bg, A1, A2, B1, B2, DR1, DR2, age, dialysis, cPRA, urgent)
  # ),
  # tar_target(
  #   antbs,
  #   simK::Abs_df_v1(candidates = cndts,
  #                origin = 'PT',
  #                seed.number = 3)
  # ),
  # dados #########
  tar_target(
    dnrs,
    simK::donors_df(n = 280,
                    replace = TRUE,
                    origin = 'PT',
                    probs = c(0.4658, 0.0343, 0.077, 0.4229),
                    lower=18, upper=75,
                    mean.age = 60, sd.age = 20,
                    uk = FALSE,
                    seed.number = 44)
  ),
  tar_target(
    cndts,
    simK::candidates_df(n = 2000,
                        replace = TRUE,
                        origin = 'PT',
                        probs.abo = c(0.44, 0.022, 0.042, 0.496),
                        probs.cpra = c(0.7, 0.1, 0.1, 0.1),
                        lower=18, upper=75,
                        mean.age = 45, sd.age = 15,
                        prob.dm = 0.12,
                        prob.urgent = 0,
                        uk = FALSE,
                        seed.number = 3) |>
      dplyr::mutate(dialysis = purrr::map2_dbl(.x = bg,
                                               .y = hiper,
                                               ~simK::dial(hiper = .y, bg = .x,
                                                           mean.dial1 = 25,
                                                           mean.dial2 = 60,
                                                           mean.dial3 = 75,
                                                           sd.dial = 20,
                                                           seed.number = NA)))|>
      dplyr::select(ID, bg, A1, A2, B1, B2, DR1, DR2, age, dialysis, cPRA, urgent)
  ),
  tar_target(
    antbs,
    simK::Abs_df(candidates = cndts,
                 origin = 'PT',
                 seed.number = 3)
  ),
  tar_target(it.numb, 1000),
  tar_target(Q2, 50),
  tar_target(Q3, 70),
  # Lima ##########
  tar_target(res_lima,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = lima,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3)
             ),
  ## extra ####
  tar_target(res_lima_extra,
             res_lima %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                           ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                             ~sum(.x$mmHLA02))
                             )
  ),
  # ET #########
  tar_target(res_et,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = et,
               n = 0,
               seed.number = 123,
               check.validity = TRUE
               )
  ),
  ## extra ####
  tar_target(res_et_extra,
             res_et %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  # EQM ###############
  tar_target(res_eqm_c01r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.1))
             ),
  tar_target(res_eqm_c01r01_extra,
             res_eqm_c01r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c01r02,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.2))
             ),
  tar_target(res_eqm_c01r02_extra,
             res_eqm_c01r02 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c01r03,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.3))
             ),
  tar_target(res_eqm_c01r03_extra,
             res_eqm_c01r03 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c01r04,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.4))
             ),
  tar_target(res_eqm_c01r04_extra,
             res_eqm_c01r04 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c01r05,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.5))
             ),
  tar_target(res_eqm_c01r05_extra,
             res_eqm_c01r05 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c02r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.2, ratio.just = 0.1))
             ),
  tar_target(res_eqm_c02r01_extra,
             res_eqm_c02r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c03r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.3, ratio.just = 0.1))
             ),
  tar_target(res_eqm_c03r01_extra,
             res_eqm_c03r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c04r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.4, ratio.just = 0.1))
             ),
  tar_target(res_eqm_c04r01_extra,
             res_eqm_c04r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_c05r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               uj.matx = uj_matx(ratio.util = 0.5, ratio.just = 0.1))
             ),
  tar_target(res_eqm_c05r01_extra,
             res_eqm_c05r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                               ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                   ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                 ~sum(.x$mmHLA02))
               )
  ),
  # EQM com priorizações do ET ##############
  tar_target(res_eqm_et_c01r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.1))
  ),
  tar_target(res_eqm_et_c01r01_extra,
             res_eqm_et_c01r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c01r02,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.2))
  ),
  tar_target(res_eqm_et_c01r02_extra,
             res_eqm_et_c01r02 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c01r03,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.3))
  ),
  tar_target(res_eqm_et_c01r03_extra,
             res_eqm_et_c01r03 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c01r04,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.4))
  ),
  tar_target(res_eqm_et_c01r04_extra,
             res_eqm_et_c01r04 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c01r05,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.5))
  ),
  tar_target(res_eqm_et_c01r05_extra,
             res_eqm_et_c01r05 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c02r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.2, ratio.just = 0.1))
  ),
  tar_target(res_eqm_et_c02r01_extra,
             res_eqm_et_c02r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c03r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.3, ratio.just = 0.1))
  ),
  tar_target(res_eqm_et_c03r01_extra,
             res_eqm_et_c03r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c04r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.4, ratio.just = 0.1))
  ),
  tar_target(res_eqm_et_c04r01_extra,
             res_eqm_et_c04r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  ),
  tar_target(res_eqm_et_c05r01,
             histoc::several(
               iteration.number = it.numb,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = Q2,
               q3 = Q3,
               SP = TRUE,
               AM = TRUE,
               mm000 = TRUE,
               uj.matx = uj_matx(ratio.util = 0.5, ratio.just = 0.1))
  ),
  tar_target(res_eqm_et_c05r01_extra,
             res_eqm_et_c05r01 %>%
               dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
               dplyr::select(it, data_extra) %>%
               dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
                                                      ~mean(.x$txs)),
                             ageDiff_avg = purrr::map_dbl(data_extra,
                                                          ~mean(.x$ageDiff)),
                             mmHLA02_n = purrr::map_dbl(data_extra,
                                                        ~sum(.x$mmHLA02))
               )
  )#,

  # # EQM-AM (com priorização do AM program) ##############
  # tar_target(res_eqm_am_c01r01,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.1))
  # ),
  # tar_target(res_eqm_am_c01r01_extra,
  #            res_eqm_am_c01r01 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c01r02,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.2))
  # ),
  # tar_target(res_eqm_am_c01r02_extra,
  #            res_eqm_am_c01r02 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c01r03,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.3))
  # ),
  # tar_target(res_eqm_am_c01r03_extra,
  #            res_eqm_am_c01r03 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c01r04,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.4))
  # ),
  # tar_target(res_eqm_am_c01r04_extra,
  #            res_eqm_am_c01r04 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c01r05,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.5))
  # ),
  # tar_target(res_eqm_am_c01r05_extra,
  #            res_eqm_am_c01r05 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c02r01,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.2, ratio.just = 0.1))
  # ),
  # tar_target(res_eqm_am_c02r01_extra,
  #            res_eqm_am_c02r01 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c03r01,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.3, ratio.just = 0.1))
  # ),
  # tar_target(res_eqm_am_c03r01_extra,
  #            res_eqm_am_c03r01 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c04r01,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.4, ratio.just = 0.1))
  # ),
  # tar_target(res_eqm_am_c04r01_extra,
  #            res_eqm_am_c04r01 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # ),
  # tar_target(res_eqm_am_c05r01,
  #            histoc::several(
  #              iteration.number = it.numb,
  #              df.donors = dnrs,
  #              df.candidates = cndts,
  #              df.abs = antbs,
  #              algorithm = eqm,
  #              n = 0,
  #              seed.number = 123,
  #              check.validity = TRUE,
  #              q2 = Q2,
  #              q3 = Q3,
  #              SP = FALSE,
  #              AM = TRUE,
  #              mm000 = FALSE,
  #              uj.matx = uj_matx(ratio.util = 0.5, ratio.just = 0.1))
  # ),
  # tar_target(res_eqm_am_c05r01_extra,
  #            res_eqm_am_c05r01 %>%
  #              dplyr::mutate(data_extra = dados_extra(dados = data)) %>%
  #              dplyr::select(it, data_extra) %>%
  #              dplyr::mutate(txs_avg = purrr::map_dbl(data_extra,
  #                                                     ~mean(.x$txs)),
  #                            ageDiff_avg = purrr::map_dbl(data_extra,
  #                                                         ~mean(.x$ageDiff)),
  #                            mmHLA02_n = purrr::map_dbl(data_extra,
  #                                                       ~sum(.x$mmHLA02))
  #              )
  # )
  ## Report #############
  ,tar_render(report_eqm,
              "ReportEQM.Rmd") # Here is our call to tar_render()

)
