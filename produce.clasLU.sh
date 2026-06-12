#!/usr/bin/env bash
set -euo pipefail

num=500000
seed=0

install/bin/clas-stringspinner-slurm $num $seed ~/j/bihadro/out/lund.sss.prod2.clasLU.l \
  --config clas12 \
  --pol-type LU \
  --cut-inclusive 11,211,-211 \
  --cut-lepton-theta 2,60 \
  --cut-pion-multiplicity 4 \
  --save-hipo \
  --set 'StringSpinner:GLGT=20.0'

install/bin/clas-stringspinner-slurm $num $seed ~/j/bihadro/out/lund.sss.prod2.clasLU.t \
  --config clas12 \
  --pol-type LU \
  --cut-inclusive 11,211,-211 \
  --cut-lepton-theta 2,60 \
  --cut-pion-multiplicity 4 \
  --save-hipo \
  --set 'StringSpinner:GLGT=0.0'
