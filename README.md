# TrackingPlots

This directory contains scripts and utilities for generating tracking performance plots and NTuple data produced by the MuColl analysis workflow.

## Overview
TrackingPlots provides ROOT-based NTupleWriter and plotting scripts to extract, visualize, and analyze tracking data. These tools are used to assess tracking efficiency, fake rate, resolution, hit density, and other key metrics for detector and reconstruction studies.

## Typical Workflow
1. **NTuple Generation**: Use the NTupleWriter scripts (included in this repository) to extract data from EDM4hep ROOT files into NTuples:
    - `WriteTracksMT.C`: Multithreaded track data extraction
    - `WriteSeedsMT.C`: Multithreaded seed data extraction
    - `WriteHitsMT.C`: Multithreaded hit data extraction
    - `WriteAllMT.C`: Runs all NTupleWriter modules at once
2. **Plotting**: Run the plotting scripts in this directory to produce summary plots from the NTuple files:
    - `PlotTracks.C`: Plots tracking efficiency, fake rate, resolution, and track quality metrics.
    - `PlotSeeds.C`: Plots seed distributions, layer occupancy, and seed resolution.
    - `PlotHits.C`: Plots hit density and event-level hit statistics.
    - `PlotAll.C`: Runs all plot scripts for a given NTuple directory and outputs all plots to a specified folder.

## Usage Example
```bash
# Generate NTuples (from NTupleWriter directory)
root -l -q 'WriteAllMT.C("/path/to/reco_", "./ntuples", 4)'

# Make plots (from TrackingPlots directory)
root -l -b -q 'PlotAll.C("./ntuples", "./plots/")'
```

## Output
Plots are saved as ROOT files in the specified output directory. These include:
- Efficiency vs. $p_T$ and $\theta$
- Fake rate
- Track parameter resolutions
- Hit density maps
- Seed and layer distributions

## Requirements
- ROOT (>=6.32)
- EDM4hep libraries