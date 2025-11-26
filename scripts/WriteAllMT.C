#include "WriteTracksMT.C"
#include "WriteSeedsMT.C"
#include "WriteHitsMT.C"

void WriteAllMT(const char* inputFilePrefix, const char* outputDir = ".", int nThreads = 4,
                const char* trackBranch = "SiTrack", bool isSubsetCollection = true,
                const char* seedBranch = "SeedTracks") {
    
    if (!inputFilePrefix) {
        printf("Usage: WriteAllMT(\"<path/prefix_>\", \"output_dir\", nThreads, \"trackBranch\", isSubset, \"seedBranch\")\n");
        printf("  inputFilePrefix: Path and filename prefix (e.g., \"data/reco_\")\n");
        printf("  outputDir: Directory for output files (default: \".\")\n");
        printf("  nThreads: Number of parallel threads (default: 4)\n");
        printf("  trackBranch: Track collection name (default: \"SiTracks\")\n");
        printf("  isSubsetCollection: Whether tracks are subset collection (default: true)\n");
        printf("  seedBranch: Seed collection name (default: \"SeedTracks\")\n");
        return;
    }
    
    TString trackOutput = TString::Format("%s/tracks_ntuple.root", outputDir);
    TString seedOutput = TString::Format("%s/seeds_ntuple.root", outputDir);
    TString hitOutput = TString::Format("%s/hits_ntuple.root", outputDir);
    
    printf("============================================\n");
    printf("Running multithreaded NTuple writers with %d threads\n", nThreads);
    printf("============================================\n\n");
    
    printf("=== Writing Tracks ===\n");
    WriteTracksMT(inputFilePrefix, trackOutput.Data(), trackBranch, isSubsetCollection, nThreads);
    printf("\n");
    
    printf("=== Writing Seeds ===\n");
    WriteSeedsMT(inputFilePrefix, seedOutput.Data(), seedBranch, nThreads);
    printf("\n");
    
    printf("=== Writing Hits ===\n");
    WriteHitsMT(inputFilePrefix, hitOutput.Data(), nThreads);
    printf("\n");
    
    printf("============================================\n");
    printf("All NTuples written to %s/\n", outputDir);
    printf("  - tracks_ntuple.root\n");
    printf("  - seeds_ntuple.root\n");
    printf("  - hits_ntuple.root\n");
    printf("============================================\n");
}
