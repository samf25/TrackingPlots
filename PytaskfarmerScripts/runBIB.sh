#!/bin/bash

source /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/k4MuCPlayground/setup_digireco.sh /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/ MuSIC_v2

PARTICLE=$1
NUM=$2

ddsim --steeringFile /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/simulation/steer_baseline.py \
        --inputFiles /global/cfs/cdirs/m5197/data/bib/FLUKA_to_SLCIO/${PARTICLE}/summary${NUM}_DET_IP.slcio \
        --outputFile /global/cfs/cdirs/m5197/data/bib/10TeV_MuSIC_edm4hep/${PARTICLE}/summary${NUM}_DET_IP.edm4hep.root \
        --numberOfEvents 1