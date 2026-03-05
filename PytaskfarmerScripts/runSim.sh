#!/bin/bash

source /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/k4MuCPlayground/setup_digireco.sh /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/ MAIA_v0

#Gen-parameters
PARTICLE=${2:-"mu-"}
PT_LOW=${3:-"10"}
PT_HIGH=${4:-"5000"}
THETA_LOW=${5:-"40"}
THETA_HIGH=${6:-"140"}

num_events=10000

INPUT_FILE="/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/Mu_pgun_barrel/gen_${PARTICLE}_pt${PT_LOW}-${PT_HIGH}_theta${THETA_LOW}-${THETA_HIGH}.edm4hep.root"
OUTPUT_FILE="/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/Mu_pgun_barrel/sim_${PARTICLE}_pt${PT_LOW}-${PT_HIGH}_theta${THETA_LOW}-${THETA_HIGH}.edm4hep.root"


ddsim --steeringFile /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/simulation/steer_baseline.py \
        --inputFiles $INPUT_FILE \
        --outputFile $OUTPUT_FILE \
        --numberOfEvents ${num_events}
