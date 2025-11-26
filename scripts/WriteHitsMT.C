#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TROOT.h"
#include <glob.h>
#include <vector>
#include <string>
#include <iostream>
#include <mutex>
#include <thread>
#include <atomic>
#include <array>

// Load EDM4hep
R__LOAD_LIBRARY(libedm4hep)
#include "edm4hep/TrackerHitPlaneData.h"
#include "edm4hep/Vector3d.h"

static std::vector<std::string> GetMatchingFiles_HitsMT(const char* pathPrefix) {
    std::vector<std::string> files;
    std::string pattern = std::string(pathPrefix) + "*.root";
    
    glob_t globResult;
    if (glob(pattern.c_str(), GLOB_TILDE, NULL, &globResult) == 0) {
        for (size_t i = 0; i < globResult.gl_pathc; ++i) {
            files.push_back(std::string(globResult.gl_pathv[i]));
        }
    }
    globfree(&globResult);
    
    return files;
}

struct HitFileResult {
    std::vector<std::array<float, 3>> hit_data;    // theta, layer, isBarrel
    std::vector<std::array<float, 4>> event_data;  // nHits_VXD, nHits_IT, nHits_OT, nHits_total
    Long64_t entries = 0;
};

static const char* branchNames_MT[] = {"VXDBarrelHits", "VXDEndcapHits", "ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits"};
static const int layerMapping_MT[] = {0, 0, 1, 1, 2, 2}; // VXD=0, IT=1, OT=2

HitFileResult ProcessHitFile(const std::string& filename) {
    HitFileResult result;
    
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) return result;
    
    TTree* tree = (TTree*)file->Get("events");
    if (!tree) { file->Close(); return result; }

    std::vector<std::vector<edm4hep::TrackerHitPlaneData>*> hitCollections(6, nullptr);
    for (int i = 0; i < 6; i++) {
        tree->SetBranchAddress(branchNames_MT[i], &hitCollections[i]);
    }
    
    Long64_t nEntries = tree->GetEntries();
    result.entries = nEntries;
    
    // Pre-allocate
    result.hit_data.reserve(nEntries * 100);
    result.event_data.reserve(nEntries);
    
    for (Long64_t evt = 0; evt < nEntries; evt++) {
        tree->GetEntry(evt);
        
        int nHits[3] = {0, 0, 0}; // VXD, IT, OT
        
        for (int collIdx = 0; collIdx < 6; collIdx++) {
            if (!hitCollections[collIdx]) continue;
            
            int layer = layerMapping_MT[collIdx];
            float isBarrel = (collIdx % 2 == 0) ? 1.0f : 0.0f;
            
            for (const auto& hit : *hitCollections[collIdx]) {
                edm4hep::Vector3d position = hit.position;
                float theta = atan2(sqrt(position.x*position.x + position.y*position.y), position.z);
                
                result.hit_data.push_back({theta, (float)layer, isBarrel});
                nHits[layer]++;
            }
        }
        
        result.event_data.push_back({(float)nHits[0], (float)nHits[1], (float)nHits[2], 
                                     (float)(nHits[0]+nHits[1]+nHits[2])});
    }
    
    file->Close();
    return result;
}

void WriteHitsMT(const char* inputFilePrefix, const char* outputFile, int nThreads = 4) {
    
    if (!inputFilePrefix) {
        printf("Usage: WriteHitsMT(\"<path/prefix_>\", \"output.root\", nThreads)\n");
        return;
    }
    
    ROOT::EnableThreadSafety();
    
    std::vector<std::string> files = GetMatchingFiles_HitsMT(inputFilePrefix);
    if (files.empty()) {
        printf("Error: No files found matching pattern %s*.root\n", inputFilePrefix);
        return;
    }
    
    printf("Found %zu files to process with %d threads\n", files.size(), nThreads);
    
    std::vector<HitFileResult> results(files.size());
    std::atomic<int> filesProcessed(0);
    std::atomic<int> nextFileIdx(0);
    
    auto worker = [&]() {
        while (true) {
            int idx = nextFileIdx.fetch_add(1);
            if (idx >= (int)files.size()) break;
            
            results[idx] = ProcessHitFile(files[idx]);
            int done = filesProcessed.fetch_add(1) + 1;
            printf("Processed %d/%zu files: %s (%lld events)\n", done, files.size(), files[idx].c_str(), results[idx].entries);
        }
    };
    
    std::vector<std::thread> threads;
    int actualThreads = std::min(nThreads, (int)files.size());
    for (int t = 0; t < actualThreads; t++) threads.emplace_back(worker);
    for (auto& t : threads) t.join();
    
    printf("\nMerging results and writing to %s...\n", outputFile);
    
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TNtuple* ntuple_hits = new TNtuple("hits", "Tracker Hits", "theta:layer:isBarrel");
    TNtuple* ntuple_events = new TNtuple("events_summary", "Event Summary", "nHits_VXD:nHits_IT:nHits_OT:nHits_total");
    
    Long64_t totalEntries = 0;
    for (const auto& result : results) {
        totalEntries += result.entries;
        for (const auto& d : result.hit_data) ntuple_hits->Fill(d[0], d[1], d[2]);
        for (const auto& d : result.event_data) ntuple_events->Fill(d[0], d[1], d[2], d[3]);
    }
    
    printf("Total events processed: %lld\n", totalEntries);
    
    outFile->Write();
    outFile->Close();
    delete outFile;
    
    printf("Done!\n");
}
