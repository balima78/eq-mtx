library(targets)
# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(data_summary) to view the results.
library(tarchetypes)

# Set target-specific options such as packages:
tar_option_set(packages = c("histoc","eq.mtx"),
               memory = 'transient',
               garbage_collection = TRUE)

# End this file with a list of target objects.
list(
  tar_target(
    dnrs,
    data("dns")
  ),
  tar_target(
    cndts,
    data("cds")
  ),
  tar_target(
    antbs,
    data("ant")
  ),
  tar_target(res_lima,
             several(
               iteration.number = 100,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = lima,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = 60,
               q3 = 80)
             ),
  tar_target(res_eqm_c01r01,
             several(
               iteration.number = 100,
               df.donors = dnrs,
               df.candidates = cndts,
               df.abs = antbs,
               algorithm = eqm,
               n = 0,
               seed.number = 123,
               check.validity = TRUE,
               q2 = 60,
               q3 = 80,
               uj.matx = uj_matx(ratio.util = 0.1, ratio.just = 0.1))
  ) # Call your custom functions.
)
