#!/bin/bash

wget -P $(dirname $0) http://portals.broadinstitute.org/collaboration/giant/images/e/eb/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz
tar -xzf $(dirname $0)/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz
rm $(dirname $0)/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz
gzip $(dirname $0)/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt
