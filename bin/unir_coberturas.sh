#!/bin/bash
set -euo pipefail

# Se asume que este script se procesa donde están las 8 carpetas
(
  echo -e "POS\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8";
  paste \
    <(awk 'NR > 1 { print $3 }' ./S1/COBERTURA/S1_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S2/COBERTURA/S2_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S3/COBERTURA/S3_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S4/COBERTURA/S4_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S5/COBERTURA/S5_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S6/COBERTURA/S6_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S7/COBERTURA/S7_cobertura) \
    <(awk 'NR > 1 { print $3 }' ./S8/COBERTURA/S8_cobertura) \
  | awk '{ print NR "\t" $0 }'
) > coberturas_finales.tsv