#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

// Load EDM4hep
R__LOAD_LIBRARY(libpodio)
R__LOAD_LIBRARY(libpodioDict)
R__LOAD_LIBRARY(libpodioRootIO)
R__LOAD_LIBRARY(libedm4hep)
R__LOAD_LIBRARY(libedm4hepDict)
#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/utils/bit_utils.h"
#include "podio/ObjectID.h"

static const int BITCreatedInSimulation = 30;
static const int BITDecayedInTracker = 27;

void SummarizeNoTrackEvents(const char* inputFile = "MC/Mu_pgun/NOBIB_reco_mu-_pt10-5000_theta10-170.edm4hep.root",
                            const char* outputFile = "MC/Mu_pgun/NOBIB_ntuples/no_track_summary.root",
                            const char* trackBranch = "SiTrack",
                            const char* seedBranch = "SeedTracks") {
    TFile* file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open input file: " << inputFile << std::endl;
        return;
    }

    TTree* tree = static_cast<TTree*>(file->Get("events"));
    if (!tree) {
        std::cerr << "Could not find events tree in " << inputFile << std::endl;
        file->Close();
        return;
    }

    // Branch pointers
    std::vector<edm4hep::MCParticleData>* mcParticles = nullptr;
    std::vector<edm4hep::TrackData>* tracks = nullptr;
    std::vector<edm4hep::TrackData>* allTracks = nullptr;
    std::vector<edm4hep::TrackData>* seeds = nullptr;
    std::vector<edm4hep::TrackState>* trackStates = nullptr;
    std::vector<podio::ObjectID>* subsetTrackIndices = nullptr;

    tree->SetBranchAddress("MCParticles", &mcParticles);

    std::string subsetIdxName = std::string(trackBranch) + "s_objIdx";
    bool hasSubset = (tree->GetBranch("AllTracks") && tree->GetBranch(subsetIdxName.c_str()));
    bool hasDirect = tree->GetBranch(trackBranch);

    if (hasSubset) {
        tree->SetBranchAddress("AllTracks", &allTracks);
        tree->SetBranchAddress("_AllTracks_trackStates", &trackStates);
        tree->SetBranchAddress(subsetIdxName.c_str(), &subsetTrackIndices);
    } else if (hasDirect) {
        tree->SetBranchAddress(trackBranch, &tracks);
        std::string tsName = "_" + std::string(trackBranch) + "_trackStates";
        if (tree->GetBranch(tsName.c_str())) tree->SetBranchAddress(tsName.c_str(), &trackStates);
    } else {
        std::cerr << "No track branch found for " << trackBranch << std::endl;
    }

    if (seedBranch && tree->GetBranch(seedBranch)) {
        tree->SetBranchAddress(seedBranch, &seeds);
    }

    // Histograms for zero-track events
    TH1F* hTheta = new TH1F("h_mc_theta", "MCParticles (no reconstructed tracks);#theta [rad];Entries", 80, 0.0, 3.2);
    TH1F* hPt = new TH1F("h_mc_pt", "MCParticles (no reconstructed tracks);p_{T} [GeV];Entries", 100, 0.0, 5500.0);
    TH2F* hPtVsTheta = new TH2F("h_mc_pt_vs_theta", "MCParticles (no reconstructed tracks);#theta [rad];p_{T} [GeV]", 60, 0.0, 3.2, 100, 0.0, 5500.0);
    TH1F* hSeedCounts = new TH1F("h_seed_counts", "SeedTracks per event (no reconstructed tracks);N seeds;Events", 11, -0.5, 10.5);
    TH2F* hSeedsVsTheta = new TH2F("h_seeds_vs_theta", "SeedTracks vs MC #theta (no reconstructed tracks);#theta [rad];N seeds", 60, 0.0, 3.2, 11, -0.5, 10.5);
    TH2F* hSeedsVsPt = new TH2F("h_seeds_vs_pt", "SeedTracks vs MC p_{T} (no reconstructed tracks);p_{T} [GeV];N seeds", 100, 0.0, 5500.0, 11, -0.5, 10.5);
    TH3F* hSeedsVsPtTheta = new TH3F("h_seeds_vs_pt_theta", "SeedTracks vs MC p_{T} vs #theta (no reconstructed tracks);p_{T} [GeV];#theta [rad];N seeds", 80, 0.0, 5500.0, 48, 0.0, 3.2, 11, -0.5, 10.5);
    TH1F* hMCPCounts = new TH1F("h_mc_counts", "Generator-status=1 charged MCParticles per event (no reconstructed tracks);N MCParticles;Events", 20, -0.5, 19.5);

    Long64_t nEntries = tree->GetEntries();
    Long64_t nNoTrackEvents = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        int nTracks = 0;
        if (hasSubset && subsetTrackIndices) {
            nTracks = subsetTrackIndices->size();
        } else if (tracks) {
            nTracks = tracks->size();
        } else if (allTracks) {
            nTracks = allTracks->size();
        }

        if (nTracks > 0) continue;

        nNoTrackEvents++;

        int nSeedsRaw = seeds ? static_cast<int>(seeds->size()) : 0;
        int nSeeds = std::min(nSeedsRaw, 10); // cap at 10 for plotting
        hSeedCounts->Fill(static_cast<double>(nSeeds));

        int mcpCount = 0;
        if (mcParticles) {
            for (const auto& mcp : *mcParticles) {
                if (mcp.generatorStatus != 1) continue;
                if (mcp.charge == 0) continue;
                if (edm4hep::utils::checkBit(mcp.simulatorStatus, BITCreatedInSimulation)) continue;
                if (edm4hep::utils::checkBit(mcp.simulatorStatus, BITDecayedInTracker)) continue;

                const edm4hep::Vector3d& mom = mcp.momentum;
                double pt = std::sqrt(mom.x * mom.x + mom.y * mom.y);
                double theta = std::atan2(std::sqrt(mom.x * mom.x + mom.y * mom.y), mom.z);

                hPt->Fill(pt);
                hTheta->Fill(theta);
                hPtVsTheta->Fill(theta, pt);
                hSeedsVsTheta->Fill(theta, nSeeds);
                hSeedsVsPt->Fill(pt, nSeeds);
                hSeedsVsPtTheta->Fill(pt, theta, nSeeds);
                mcpCount++;
            }
        }
        hMCPCounts->Fill(static_cast<double>(mcpCount));
    }

    TFile outFile(outputFile, "RECREATE");
    hTheta->Write();
    hPt->Write();
    hPtVsTheta->Write();
    hSeedsVsTheta->Write();
    hSeedsVsPt->Write();
    hSeedsVsPtTheta->Write();
    hSeedCounts->Write();
    hMCPCounts->Write();
    outFile.Close();

    std::cout << "Processed " << nEntries << " events; " << nNoTrackEvents << " had zero "
              << trackBranch << " tracks. Histograms saved to " << outputFile << std::endl;

    file->Close();
}
