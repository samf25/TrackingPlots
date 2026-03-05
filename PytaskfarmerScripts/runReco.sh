#!/bin/bash

source /global/cfs/cdirs/m5197/sferrar2/TrackingPaper/mucoll-benchmarks/k4MuCPlayground/setup_digireco.sh /global/cfs/cdirs/m5197/sferrar2/TrackingPaper/mucoll-benchmarks/ MAIA_v0
cd /global/cfs/cdirs/m5197/sferrar2/TrackingPaper/TemplateWorkspace/
source setup.sh

#Gen-parameters
PARTICLE=${2:-"mu-"}
PT_LOW=${3:-"10"}
PT_HIGH=${4:-"5000"}
THETA_LOW=${5:-"10"}
THETA_HIGH=${6:-"170"}

INDEX=$1

INPUT_FILE="/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/digi_${PARTICLE}_pt${PT_LOW}-${PT_HIGH}_theta${THETA_LOW}-${THETA_HIGH}_${INDEX}.edm4hep.root"
OUTPUT_FILE="/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/reco_${PARTICLE}_pt${PT_LOW}-${PT_HIGH}_theta${THETA_LOW}-${THETA_HIGH}_${INDEX}.edm4hep.root"


k4run /global/cfs/cdirs/m5197/sferrar2/TrackingPaper/mucoll-benchmarks/reconstruction/reco_steer.py \
        --IOSvc.Input $INPUT_FILE \
        --IOSvc.Output $OUTPUT_FILE \
        --TrackingThreads 5 \
        --useMT \
        -n 100 \
        --numThreads 10