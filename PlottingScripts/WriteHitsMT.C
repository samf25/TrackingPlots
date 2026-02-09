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

static const char* branchNames_MT[] = {"VXDBarrelHits", "VXDEndcapHits", "ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits"};
static const int layerMapping_MT[] = {0, 0, 1, 1, 2, 2}; // VXD=0, IT=1, OT=2

void ProcessAndWriteHitFile(const std::string& filename, TNtuple* ntuple_hits, TNtuple* ntuple_events, 
                            std::mutex& writeMutex, std::atomic<Long64_t>& totalEntries) {
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) return;
    
    TTree* tree = (TTree*)file->Get("events");
    if (!tree) { file->Close(); return; }

    std::vector<std::vector<edm4hep::TrackerHitPlaneData>*> hitCollections(6, nullptr);
    for (int i = 0; i < 6; i++) {
        tree->SetBranchAddress(branchNames_MT[i], &hitCollections[i]);
    }
    
    Long64_t nEntries = tree->GetEntries();
    totalEntries += nEntries;
    
    // Local buffers for this file
    std::vector<std::array<float, 3>> hit_data;
    std::vector<std::array<float, 4>> event_data;
    hit_data.reserve(nEntries * 100);
    event_data.reserve(nEntries);
    
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
                
                hit_data.push_back({theta, (float)layer, isBarrel});
                nHits[layer]++;
            }
        }
        
        event_data.push_back({(float)nHits[0], (float)nHits[1], (float)nHits[2], 
                             (float)(nHits[0]+nHits[1]+nHits[2])});
    }
    
    file->Close();
    
    // Write to TNtuples with mutex protection
    {
        std::lock_guard<std::mutex> lock(writeMutex);
        for (const auto& d : hit_data) ntuple_hits->Fill(d[0], d[1], d[2]);
        for (const auto& d : event_data) ntuple_events->Fill(d[0], d[1], d[2], d[3]);
    }
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
    
    // Create output file and TNtuples upfront
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TNtuple* ntuple_hits = new TNtuple("hits", "Tracker Hits", "theta:layer:isBarrel");
    TNtuple* ntuple_events = new TNtuple("events_summary", "Event Summary", "nHits_VXD:nHits_IT:nHits_OT:nHits_total");
    
    std::mutex writeMutex;
    std::atomic<int> filesProcessed(0);
    std::atomic<int> nextFileIdx(0);
    std::atomic<Long64_t> totalEntries(0);
    
    auto worker = [&]() {
        while (true) {
            int idx = nextFileIdx.fetch_add(1);
            if (idx >= (int)files.size()) break;
            
            ProcessAndWriteHitFile(files[idx], ntuple_hits, ntuple_events, writeMutex, totalEntries);
            int done = filesProcessed.fetch_add(1) + 1;
            printf("Processed %d/%zu files: %s\n", done, files.size(), files[idx].c_str());
        }
    };
    
    std::vector<std::thread> threads;
    int actualThreads = std::min(nThreads, (int)files.size());
    for (int t = 0; t < actualThreads; t++) threads.emplace_back(worker);
    for (auto& t : threads) t.join();
    
    printf("\nTotal events processed: %lld\n", totalEntries.load());
    printf("Writing to %s...\n", outputFile);
    
    outFile->Write();
    outFile->Close();
    delete outFile;
    
    printf("Done!\n");
}
