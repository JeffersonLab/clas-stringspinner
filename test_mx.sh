#!/usr/bin/env bash
set -e
ninja -C build
build/clas-stringspinner \
  --patch-boost none \
  --seed 29877 --trig 10000 --cut-inclusive 11,211,-211 \
  |grep M_X| awk '{print $2}' > output.before_patch.dat
build/clas-stringspinner \
  --patch-boost beam \
  --seed 29877 --trig 10000 --cut-inclusive 11,211,-211 \
  |grep M_X| awk '{print $2}' > output.after_patch.dat
paste output.before_patch.dat output.after_patch.dat > output.both.dat
root -b -q draw_mx_2.C'("output.both.dat")'
nsxiv output.both.dat.png
