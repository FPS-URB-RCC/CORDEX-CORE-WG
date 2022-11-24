#!/bin/bash
# Download required data

DEST_FOLDER="."

if ! test -f ${DEST_FOLDER}/koppen_1901-2010.tsv ; then
  wget http://hanschen.org/uploads/koppen/data.zip
  unzip data.zip
  mv data/koppen_1901-2010.tsv ${DEST_FOLDER}/
  mv data/README.txt ${DEST_FOLDER}/koppen_README.txt
  rm data.zip
fi

if ! test -f ${DEST_FOLDER}/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg ; then
  wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_FUA_UCDB2015_GLOBE_R2019A/V1-0/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.zip
  unzip GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.zip
  mv GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg ${DEST_FOLDER}/
  mv GHSL_FUA_2019.pdf ${DEST_FOLDER}/
  rm GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.zip
fi

if ! test -f ${DEST_FOLDER}/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg ; then
  wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip
  unzip GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip
  mv GHS_STAT_UCDB2015MT_GLOBE_R2019A/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg ${DEST_FOLDER}/
  mv GHS_STAT_UCDB2015MT_GLOBE_R2019A/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0_web.pdf ${DEST_FOLDER}/
  rm GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip
  rm -r GHS_STAT_UCDB2015MT_GLOBE_R2019A/
fi

