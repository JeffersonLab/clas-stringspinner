#!/usr/bin/env bash
if [ $# -ne 1 ]; then
  echo "USAGE: $0 [EXECUTABLE]" >&2
  exit 1
fi

$1 --help |\
  perl -pe 's/(--.*?)([^a-zA-Z-])/<font color="#ff00ff">$1<\/font>$2/' |\
  perl -pe 's/(^[^\s].*)/<font color="#ffffff">$1<\/font>/' |\
  sed 's;https://github.com.*;<a href="&">&</a>;g' \
  > tmp_usage

date > tmp_date

cat $(dirname $0)/usage.html |\
  sed $'/@usage@/{r tmp_usage\nd}' |\
  sed $'/@date@/{r tmp_date\nd}'

rm -f tmp_usage tmp_date
