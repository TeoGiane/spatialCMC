#!/bin/sh

# Set working directory
cd "$(dirname "$0")"
echo "BASH Current directory: $PWD\n"

# Getting raw data
wget -P $PWD/raw/ https://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip
wget -P $PWD/raw/ https://www.istat.it/wp-content/uploads/2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_5marzo.zip
wget -P $PWD/raw/ https://demo.istat.it/data/posas/POSAS_2021_it_Comuni.zip
wget -P $PWD/raw/ https://www1.finanze.gov.it/finanze/analisi_stat/public/v_4_0_0/contenuti/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.zip

# Un-zipping files
unzip -d $PWD/raw $PWD/raw/Limiti01012020.zip
unzip -d $PWD/raw $PWD/raw/Dataset-decessi-comunali-giornalieri-e-tracciato-record_5marzo.zip
unzip -d $PWD/raw $PWD/raw/POSAS_2021_it_Comuni.zip
unzip -d $PWD/raw $PWD/raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.zip

# Running script to generate data
echo "\n"
Rscript --vanilla src/build_full_dataset.R

# Clean-up raw folder
rm -rf $PWD/raw/
