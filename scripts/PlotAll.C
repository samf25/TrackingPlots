#include "TSystem.h"
#include <iostream>
#include "PlotTracks.C"
#include "PlotSeeds.C"
#include "PlotHits.C"

void PlotAll(const char* inputDir, const char* outputDir) {
    
    if (!inputDir || !outputDir) {
        printf("Usage: PlotAll(\"<ntuple_directory>\", \"<output_directory>\")\n");
        printf("Example: PlotAll(\"./ntuple_output\", \"./plots\")\n");
        return;
    }
    
    // Create output directory if it doesn't exist
    gSystem->mkdir(outputDir, true);
    
    std::string inDir = inputDir;
    if (inDir.back() != '/') inDir += '/';
    
    std::string outDir = outputDir;
    if (outDir.back() != '/') outDir += '/';
    
    // Run each plotter
    printf("\n========================================\n");
    printf("Running PlotTracks...\n");
    printf("========================================\n");
    PlotTracks((inDir + "tracks_ntuple.root").c_str(), (outDir + "tracks_plots.root").c_str());
    
    printf("\n========================================\n");
    printf("Running PlotSeeds...\n");
    printf("========================================\n");
    PlotSeeds((inDir + "seeds_ntuple.root").c_str(), (outDir + "seeds_plots.root").c_str());
    
    printf("\n========================================\n");
    printf("Running PlotHits...\n");
    printf("========================================\n");
    PlotHits((inDir + "hits_ntuple.root").c_str(), (outDir + "hits_plots.root").c_str());
    
    printf("\n========================================\n");
    printf("All plotters complete!\n");
    printf("Output files:\n");
    printf("  %stracks_plots.root\n", outDir.c_str());
    printf("  %sseeds_plots.root\n", outDir.c_str());
    printf("  %shits_plots.root\n", outDir.c_str());
    printf("========================================\n");
}
