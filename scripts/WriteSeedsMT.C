#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "DDSegmentation/BitFieldCoder.h"
#include <glob.h>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <atomic>

// Load EDM4hep
R__LOAD_LIBRARY(libedm4hep)
R__LOAD_LIBRARY(libpodio)
R__LOAD_LIBRARY(libpodioRootIO)
#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHitPlaneData.h"
#include "podio/ObjectID.h"
#include "podio/utilities/RootHelpers.h"

static std::vector<std::string> GetMatchingFiles_SeedsMT(const char* pathPrefix) {
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

struct SeedFileResult {
    std::vector<std::array<float, 2>> seed_data;       // theta, isMatched
    std::vector<std::array<float, 10>> matched_data;   // resolution fields
    std::vector<std::array<float, 3>> layer_data;      // layer, isBarrel, isMatched
    std::vector<std::array<float, 4>> event_data;      // nSeeds, nMatched, nUnmatched, avgSeedsPerMCP
    Long64_t entries = 0;
};

static std::vector<std::string> trackerHitCollections_MT = {
    "ITBarrelHits", "ITEndcapHits", "VXDBarrelHits", 
    "VXDEndcapHits", "OTBarrelHits", "OTEndcapHits"
};
static std::vector<std::string> simTrackerHitCollections_MT = {
    "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", 
    "VertexBarrelCollection", "VertexEndcapCollection", 
    "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"
};

SeedFileResult ProcessSeedFile(const std::string& filename, const char* branchName) {
    SeedFileResult result;
    dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder("system:5,side:-2,layer:6,module:11,sensor:8");
    
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) return result;
    
    TTree* tree = (TTree*)file->Get("events");
    if (!tree) { file->Close(); return result; }
    
    // Get collection ID mapping
    TTree* meta = static_cast<TTree*>(file->Get("podio_metadata"));
    std::vector<podio::root_utils::CollectionWriteInfo>* infos = nullptr;
    std::unordered_map<std::string, uint32_t> name2id;
    if (meta) {
        meta->SetBranchAddress("events___CollectionTypeInfo", &infos);
        meta->GetEntry(0);
        if (infos) {
            for (const auto& ci : *infos) name2id[ci.name] = ci.collectionID;
        }
    }
    
    std::unordered_map<uint32_t, size_t> idToSlot;
    for (size_t slot = 0; slot < trackerHitCollections_MT.size(); ++slot) {
        auto it = name2id.find(trackerHitCollections_MT[slot]);
        if (it != name2id.end()) idToSlot[it->second] = slot;
    }

    std::vector<edm4hep::MCParticleData>* mcParticles = nullptr;
    std::vector<edm4hep::TrackData>* tracks = nullptr;
    std::vector<edm4hep::TrackState>* trackStates = nullptr;
    std::vector<podio::ObjectID>* seedTrackerHits = nullptr;
    std::vector<std::vector<edm4hep::TrackerHitPlaneData>*> trackerHits(trackerHitCollections_MT.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> toRelations(trackerHitCollections_MT.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> fromRelations(trackerHitCollections_MT.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> simParticles(trackerHitCollections_MT.size(), nullptr);

    tree->SetBranchAddress("MCParticles", &mcParticles);
    tree->SetBranchAddress(branchName, &tracks);
    tree->SetBranchAddress(("_" + std::string(branchName) + "_trackStates").c_str(), &trackStates);
    tree->SetBranchAddress(("_" + std::string(branchName) + "_trackerHits").c_str(), &seedTrackerHits);
    for (size_t i = 0; i < trackerHitCollections_MT.size(); ++i) {
        tree->SetBranchAddress((trackerHitCollections_MT[i]).c_str(), &trackerHits[i]);
        tree->SetBranchAddress(("_" + trackerHitCollections_MT[i] + "Relations_to").c_str(), &toRelations[i]);
        tree->SetBranchAddress(("_" + trackerHitCollections_MT[i] + "Relations_from").c_str(), &fromRelations[i]);
        tree->SetBranchAddress(("_" + simTrackerHitCollections_MT[i] + "_particle").c_str(), &simParticles[i]);
    }

    Long64_t nEntries = tree->GetEntries();
    result.entries = nEntries;
    
    for (Long64_t evt = 0; evt < nEntries; evt++) {
        tree->GetEntry(evt);

        int matchedSeeds = 0, unmatchedSeeds = 0, totalSeeds = 0;
        std::unordered_map<unsigned int, int> mcpMatchCount;

        if (!tracks || !trackStates) continue;
        
        for (const auto& trk : *tracks) {
            totalSeeds++;
            const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
            double tanLambda = firstTrackState.tanLambda;
            float theta = atan2(1.0, tanLambda);

            std::vector<int> matchedMCPIndices;
            for (unsigned int hitIdx = trk.trackerHits_begin; hitIdx < trk.trackerHits_end; ++hitIdx) {
                if (!seedTrackerHits || hitIdx >= seedTrackerHits->size()) continue;
                const auto& hitID = (*seedTrackerHits)[hitIdx];
                
                auto itSlot = idToSlot.find(hitID.collectionID);
                if (itSlot == idToSlot.end()) continue;
                const size_t collIdx = itSlot->second;

                const auto* fromVec = fromRelations[collIdx];
                const auto* toVec = toRelations[collIdx];
                const auto* mcVec = simParticles[collIdx];
                if (!fromVec || !toVec || !mcVec) continue;
                
                for (size_t relIdx = 0; relIdx < fromVec->size(); ++relIdx) {
                    const auto& fromOID = fromVec->at(relIdx);
                    if (fromOID.index == hitID.index && fromOID.collectionID == hitID.collectionID) {
                        unsigned int mcIdx = toVec->at(relIdx).index;
                        if (mcIdx < mcVec->size()) matchedMCPIndices.push_back(mcVec->at(mcIdx).index);
                        break;
                    }
                }
            }

            bool isMatched = false;
            unsigned int matchedMCPIndex = 0;
            if (matchedMCPIndices.size() >= 3) {
                std::unordered_map<unsigned int, int> mcpCount;
                for (unsigned int mcpIdx : matchedMCPIndices) mcpCount[mcpIdx]++;
                for (const auto& pair : mcpCount) {
                    if (pair.second >= 2) {
                        isMatched = true;
                        matchedMCPIndex = pair.first;
                        mcpMatchCount[matchedMCPIndex]++;
                        break;
                    }
                }
            }

            result.seed_data.push_back({theta, isMatched ? 1.0f : 0.0f});

            // Layer info
            for (unsigned int hitIdx = trk.trackerHits_begin; hitIdx < trk.trackerHits_end; ++hitIdx) {
                if (!seedTrackerHits || hitIdx >= seedTrackerHits->size()) continue;
                const auto& hitID = (*seedTrackerHits)[hitIdx];
                auto itSlot = idToSlot.find(hitID.collectionID);
                if (itSlot == idToSlot.end()) continue;
                const size_t collIdx = itSlot->second;
                const auto* hitVec = trackerHits[collIdx];
                if (!hitVec || static_cast<size_t>(hitID.index) >= hitVec->size()) continue;
                
                int layer = bitFieldCoder.get(hitVec->at(hitID.index).cellID, "layer");
                float isBarrel = (collIdx % 2 == 0) ? 1.0f : 0.0f;
                result.layer_data.push_back({(float)layer, isBarrel, isMatched ? 1.0f : 0.0f});
            }

            if (isMatched) {
                matchedSeeds++;
                const auto& matchedMCP = (*mcParticles)[matchedMCPIndex];
                const edm4hep::Vector3d& mom = matchedMCP.momentum;
                float true_pt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                float true_theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                float true_d0 = std::sqrt(matchedMCP.vertex.x * matchedMCP.vertex.x + matchedMCP.vertex.y * matchedMCP.vertex.y);
                float true_z0 = matchedMCP.vertex.z;
                
                float reco_pt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
                float reco_d0 = firstTrackState.D0;
                float reco_z0 = firstTrackState.Z0;
                
                float true_q_over_pt = matchedMCP.charge / true_pt;
                float reco_q_over_pt = matchedMCP.charge / reco_pt;
                float res_q_over_pt = (reco_q_over_pt - true_q_over_pt) / true_q_over_pt;
                float res_d0 = reco_d0 - true_d0;
                float res_z0 = reco_z0 - true_z0;
                
                result.matched_data.push_back({true_pt, true_theta, reco_pt, reco_d0, reco_z0, 
                                              true_d0, true_z0, res_q_over_pt, res_d0, res_z0});
            } else {
                unmatchedSeeds++;
            }
        }

        float avgSeedsPerMCP = 0;
        if (!mcpMatchCount.empty()) {
            float total = 0;
            for (const auto& pair : mcpMatchCount) total += pair.second;
            avgSeedsPerMCP = total / mcpMatchCount.size();
        }
        
        result.event_data.push_back({(float)totalSeeds, (float)matchedSeeds, (float)unmatchedSeeds, avgSeedsPerMCP});
    }
    
    file->Close();
    return result;
}

void WriteSeedsMT(const char* inputFilePrefix, const char* outputFile, const char* branchName = "SeedTracks", int nThreads = 4) {
    
    if (!inputFilePrefix) {
        printf("Usage: WriteSeedsMT(\"<path/prefix_>\", \"output.root\", branchName, nThreads)\n");
        return;
    }
    
    ROOT::EnableThreadSafety();
    
    std::vector<std::string> files = GetMatchingFiles_SeedsMT(inputFilePrefix);
    if (files.empty()) {
        printf("Error: No files found matching pattern %s*.root\n", inputFilePrefix);
        return;
    }
    
    printf("Found %zu files to process with %d threads\n", files.size(), nThreads);
    
    std::vector<SeedFileResult> results(files.size());
    std::atomic<int> filesProcessed(0);
    std::atomic<int> nextFileIdx(0);
    
    auto worker = [&]() {
        while (true) {
            int idx = nextFileIdx.fetch_add(1);
            if (idx >= (int)files.size()) break;
            
            results[idx] = ProcessSeedFile(files[idx], branchName);
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
    TNtuple* ntuple_seeds = new TNtuple("seeds", "All Seeds", "theta:isMatched");
    TNtuple* ntuple_matched = new TNtuple("matched_seeds", "Matched Seeds",
        "true_pt:true_theta:reco_pt:reco_d0:reco_z0:true_d0:true_z0:res_q_over_pt:res_d0:res_z0");
    TNtuple* ntuple_layers = new TNtuple("seed_layers", "Seed Hit Layers", "layer:isBarrel:isMatched");
    TNtuple* ntuple_events = new TNtuple("events_summary", "Event Summary", "nSeeds:nMatched:nUnmatched:avgSeedsPerMCP");
    
    Long64_t totalEntries = 0;
    for (const auto& result : results) {
        totalEntries += result.entries;
        for (const auto& d : result.seed_data) ntuple_seeds->Fill(d[0], d[1]);
        for (const auto& d : result.matched_data) ntuple_matched->Fill(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9]);
        for (const auto& d : result.layer_data) ntuple_layers->Fill(d[0], d[1], d[2]);
        for (const auto& d : result.event_data) ntuple_events->Fill(d[0], d[1], d[2], d[3]);
    }
    
    printf("Total events processed: %lld\n", totalEntries);
    
    outFile->Write();
    outFile->Close();
    delete outFile;
    
    printf("Done!\n");
}
