#!/bin/bash
FileIn="$1"
JSON="$2"
o2-analysis-timestamp -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-trackselection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-track-propagation -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-bc-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-tracks-extra-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-event-selection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-finder -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-fje-filter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-emcal-correction-task -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-event-selection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-deriveddata-producer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-full-jet-spectra-pp -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
