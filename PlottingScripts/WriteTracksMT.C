#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TROOT.h"
#include <glob.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
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
#include "edm4hep/utils/bit_utils.h"
#include "podio/ObjectID.h"

static const int BITCreatedInSimulation = 30;
static const int BITDecayedInTracker = 27;

static std::vector<std::string> GetMatchingFiles_TracksMT(const char* pathPrefix) {
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

// Data structures to hold results from each file - REMOVED, now using incremental writing

void ProcessAndWriteTrackFile(const std::string& filename, const char* branchName, bool isSubsetCollection,
                              TNtuple* ntuple_truth, TNtuple* ntuple_tracks, TNtuple* ntuple_matched, TNtuple* ntuple_events,
                              std::mutex& writeMutex, std::atomic<Long64_t>& totalEntries) {
    
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) return;
    
    TTree* tree = (TTree*)file->Get("events");
    if (!tree) {
        file->Close();
        return;
    }
    
    // Set up branches
    std::vector<edm4hep::TrackState>* trackStates = nullptr;
    std::vector<edm4hep::MCParticleData>* mcParticles = nullptr;
    std::vector<edm4hep::TrackData>* tracks = nullptr;
    std::vector<podio::ObjectID>* toRelations = nullptr;
    std::vector<podio::ObjectID>* fromRelations = nullptr;
    std::vector<podio::ObjectID>* subsetTrackIndices = nullptr;
    
    tree->SetBranchAddress("MCParticles", &mcParticles);
    if (isSubsetCollection) {
        tree->SetBranchAddress("AllTracks", &tracks);
        tree->SetBranchAddress("_AllTracks_trackStates", &trackStates);
        tree->SetBranchAddress((std::string(branchName)+"s_objIdx").c_str(), &subsetTrackIndices);
    } else {
        tree->SetBranchAddress(branchName, &tracks);
        tree->SetBranchAddress(("_" + std::string(branchName) + "_trackStates").c_str(), &trackStates);
    }
    tree->SetBranchAddress(("_" + std::string(branchName) + "Relations_to").c_str(), &toRelations);
    tree->SetBranchAddress(("_" + std::string(branchName) + "Relations_from").c_str(), &fromRelations);

    Long64_t nEntries = tree->GetEntries();
    totalEntries += nEntries;
    
    // Local buffers for this file
    std::vector<std::array<float, 3>> truth_data;
    std::vector<std::array<float, 5>> track_data;
    std::vector<std::array<float, 11>> matched_data;
    std::vector<std::array<float, 4>> event_data;
    truth_data.reserve(nEntries * 2);
    track_data.reserve(nEntries * 10);
    matched_data.reserve(nEntries * 5);
    event_data.reserve(nEntries);
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        int nTruths = 0, nMatched = 0, nFake = 0;
        
        // Select acceptable MC particles
        std::vector<edm4hep::MCParticleData> mcpSet;
        if (mcParticles) {
            for (const auto& mcp : *mcParticles) {
                if (mcp.generatorStatus != 1) continue;
                if (mcp.charge == 0) continue;
                if (edm4hep::utils::checkBit(mcp.simulatorStatus, BITCreatedInSimulation)) continue;
                if (edm4hep::utils::checkBit(mcp.simulatorStatus, BITDecayedInTracker)) continue;

                const edm4hep::Vector3d& mom = mcp.momentum;
                double pt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                double theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                
                mcpSet.push_back(mcp);
                truth_data.push_back({(float)pt, (float)theta, (float)mcp.charge});
                nTruths++;
            }
        }
        
        // Process tracks
        std::vector<edm4hep::TrackData> trkSet;
        std::vector<int> trkIndices;
        if (tracks && trackStates) {
            if (isSubsetCollection && subsetTrackIndices) {
                for (const auto& idxObj : *subsetTrackIndices) {
                    const auto& trk = (*tracks)[idxObj.index];
                    trkSet.push_back(trk);
                    trkIndices.push_back(idxObj.index);
                }
            } else {
                for (size_t t = 0; t < tracks->size(); t++) {
                    trkSet.push_back((*tracks)[t]);
                    trkIndices.push_back(t);
                }
            }
        }
        
        std::vector<bool> trackMatched(trkSet.size(), false);
        
        // Process MC Relations
        if (toRelations && fromRelations && tracks && trackStates) {
            for (size_t r = 0; r < toRelations->size(); ++r) {
                const edm4hep::MCParticleData& mcpObj = mcParticles->at((*toRelations)[r].index);
                const edm4hep::TrackData& trkObj = tracks->at((*fromRelations)[r].index);
                
                auto itMC = std::find_if(mcpSet.begin(), mcpSet.end(), 
                    [&mcpObj](const edm4hep::MCParticleData& obj) { 
                        return (obj.momentum.x == mcpObj.momentum.x && 
                               obj.momentum.y == mcpObj.momentum.y && 
                               obj.momentum.z == mcpObj.momentum.z &&
                               obj.charge == mcpObj.charge);
                    });
                
                if (itMC != mcpSet.end()) {
                    int trkIdx = (*fromRelations)[r].index;
                    auto itTrk = std::find(trkIndices.begin(), trkIndices.end(), trkIdx);
                    
                    if (itTrk != trkIndices.end()) {
                        size_t setIdx = itTrk - trkIndices.begin();
                        if (!trackMatched[setIdx]) {
                            trackMatched[setIdx] = true;
                            
                            const edm4hep::Vector3d& mom = mcpObj.momentum;
                            float true_pt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                            float true_theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                            float true_d0 = std::sqrt(mcpObj.vertex.x * mcpObj.vertex.x + mcpObj.vertex.y * mcpObj.vertex.y);
                            float true_z0 = mcpObj.vertex.z;
                            
                            const auto& firstTrackState = (*trackStates)[trkObj.trackStates_begin];
                            float reco_pt = fabs(0.3 * 5.0 / firstTrackState.omega / 1000);
                            float reco_d0 = firstTrackState.D0;
                            float reco_z0 = firstTrackState.Z0;
                            float nHits = trkObj.trackerHits_end - trkObj.trackerHits_begin;
                            float nHoles = trkObj.Nholes;
                            float chi2ndof = (trkObj.ndf > 0) ? trkObj.chi2 / trkObj.ndf : -1;
                            
                            matched_data.push_back({true_pt, true_theta, (float)mcpObj.charge, reco_pt, 
                                                   reco_d0, reco_z0, true_d0, true_z0, nHits, nHoles, chi2ndof});
                            nMatched++;
                        }
                    }
                }
            }
        }
        
        // Fill track data
        for (size_t t = 0; t < trkSet.size(); t++) {
            const auto& trk = trkSet[t];
            const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
            float trackPt = fabs(0.3 * 5.0 / firstTrackState.omega / 1000);
            float nHits = trk.trackerHits_end - trk.trackerHits_begin;
            float nHoles = trk.Nholes;
            float chi2ndof = (trk.ndf > 0) ? trk.chi2 / trk.ndf : -1;
            float isReal = trackMatched[t] ? 1.0f : 0.0f;
            
            track_data.push_back({trackPt, nHits, nHoles, chi2ndof, isReal});
            if (!trackMatched[t]) nFake++;
        }
        
        event_data.push_back({(float)trkSet.size(), (float)nTruths, (float)nMatched, (float)nFake});
    }
    
    file->Close();
    
    // Write to TNtuples with mutex protection
    {
        std::lock_guard<std::mutex> lock(writeMutex);
        for (const auto& d : truth_data) ntuple_truth->Fill(d[0], d[1], d[2]);
        for (const auto& d : track_data) ntuple_tracks->Fill(d[0], d[1], d[2], d[3], d[4]);
        for (const auto& d : matched_data) ntuple_matched->Fill(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10]);
        for (const auto& d : event_data) ntuple_events->Fill(d[0], d[1], d[2], d[3]);
    }
}

void WriteTracksMT(const char* inputFilePrefix, const char* outputFile, const char* branchName = "SiTrack", bool isSubsetCollection = true, int nThreads = 4) {
    
    if (!inputFilePrefix) {
        printf("Usage: WriteTracksMT(\"<path/prefix_>\", \"output.root\", branchName, isSubsetCollection, nThreads)\n");
        return;
    }
    
    ROOT::EnableThreadSafety();
    
    std::vector<std::string> files = GetMatchingFiles_TracksMT(inputFilePrefix);
    if (files.empty()) {
        printf("Error: No files found matching pattern %s*.root\n", inputFilePrefix);
        return;
    }
    
    printf("Found %zu files to process with %d threads\n", files.size(), nThreads);
    
    // Create output file and TNtuples upfront
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TNtuple* ntuple_truth = new TNtuple("truth", "MC Truth Particles", "pt:theta:charge");
    TNtuple* ntuple_tracks = new TNtuple("tracks", "All Tracks", "pt:nHits:nHoles:chi2ndof:isReal");
    TNtuple* ntuple_matched = new TNtuple("matched", "Matched Tracks", 
        "true_pt:true_theta:true_charge:reco_pt:reco_d0:reco_z0:true_d0:true_z0:nHits:nHoles:chi2ndof");
    TNtuple* ntuple_events = new TNtuple("events_summary", "Event Summary", "nTracks:nTruths:nMatched:nFake");
    
    std::mutex writeMutex;
    std::atomic<int> filesProcessed(0);
    std::atomic<int> nextFileIdx(0);
    std::atomic<Long64_t> totalEntries(0);
    
    // Worker function
    auto worker = [&]() {
        while (true) {
            int idx = nextFileIdx.fetch_add(1);
            if (idx >= (int)files.size()) break;
            
            ProcessAndWriteTrackFile(files[idx], branchName, isSubsetCollection, ntuple_truth, ntuple_tracks, ntuple_matched, ntuple_events, writeMutex, totalEntries);
            int done = filesProcessed.fetch_add(1) + 1;
            printf("Processed %d/%zu files: %s\n", done, files.size(), files[idx].c_str());
        }
    };
    
    // Launch threads
    std::vector<std::thread> threads;
    int actualThreads = std::min(nThreads, (int)files.size());
    for (int t = 0; t < actualThreads; t++) {
        threads.emplace_back(worker);
    }
    
    // Wait for completion
    for (auto& t : threads) {
        t.join();
    }
    
    printf("\nTotal events processed: %lld\n", totalEntries.load());
    printf("Writing to %s...\n", outputFile);
    
    outFile->Write();
    outFile->Close();
    delete outFile;
    
    printf("Done!\n");
}
