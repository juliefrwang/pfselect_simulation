#!/bin/bash

Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation.R equal TRUE
Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation.R not-equal TRUE
# Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation-exclude-AMR.R equal TRUE TRUE
# Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation-exclude-AMR.R equal FALSE TRUE
# Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation-exclude-AMR.R not_equal TRUE TRUE
# Rscript ~/Project/Knockoff/debug_trials/Control-var-simulation-exclude-AMR.R not_equal FALSE TRUE
