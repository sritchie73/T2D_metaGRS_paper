#!/bin/bash

wget -P $(dirname $0) http://portals.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz 
tar -xzf $(dirname $0)/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
rm $(dirname $0)/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
gzip $(dirname $0)/GIANT_2015_WCadjBMI_COMBINED_EUR.txt
