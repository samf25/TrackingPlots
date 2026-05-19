// ============================================================
//  RunAnalysis.C  –  Steering file for NTuple writing
//
//  Configure all runtime options in an external .conf file and
//  pass it to RunAnalysis(configPath).
//
//  Run with default config in current working directory:
//    root -l -q RunAnalysis.C
//
//  Run with an explicit config file:
//    root -l -q 'RunAnalysis.C("RunAnalysis.conf")'
// ============================================================

#include "WriteAllMT.C"   // pulls in all sub-writers + SelectionConfig.h

#include <limits>
#include <string>

#include "TEnv.h"

struct AnalysisConfig {
    std::string inputFilePrefix;
    std::string outputDir;
    int nThreads = 4;
    int nThreadsTracks = -1;
    int nThreadsSeeds = -1;
    int nThreadsHits = -1;
    std::string trackBranch;
    bool isSubsetCollection = true;
    std::string seedBranch;
    EventSelectionConfig evtSel;
    TrackSelectionConfig trkSel;
};

AnalysisConfig DefaultAnalysisConfig() {
    AnalysisConfig cfg;

    cfg.inputFilePrefix = "/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/Mu_pgun_MAIA_v0/reco_mu-_p1-5000_theta15-165_";
    cfg.outputDir = "/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/ntuples_Mu_pgun_MAIA_v0_ENDCAP";
    cfg.nThreads = 50;
    cfg.nThreadsTracks = -1;
    cfg.nThreadsSeeds = -1;
    cfg.nThreadsHits = -1;
    cfg.trackBranch = "SiTrack";
    cfg.isSubsetCollection = true;
    cfg.seedBranch = "SeedTracks";

    cfg.evtSel.absEtaMin = 0.7f;
    cfg.evtSel.absEtaMax = std::numeric_limits<float>::max();
    cfg.evtSel.thetaMin  = 0.0f;
    cfg.evtSel.thetaMax  = (float)M_PI;
    cfg.evtSel.ptMin     = 0.0f;
    cfg.evtSel.ptMax     = std::numeric_limits<float>::max();

    cfg.trkSel.ptMin     = 0.5f;
    cfg.trkSel.ptMax     = std::numeric_limits<float>::max();
    cfg.trkSel.absEtaMin = 0.0f;
    cfg.trkSel.absEtaMax = std::numeric_limits<float>::max();
    cfg.trkSel.thetaMin  = 0.0f;
    cfg.trkSel.thetaMax  = (float)M_PI;
    cfg.trkSel.phiMin    = -(float)M_PI;
    cfg.trkSel.phiMax    =  (float)M_PI;
    cfg.trkSel.d0Min     = -std::numeric_limits<float>::max();
    cfg.trkSel.d0Max     =  std::numeric_limits<float>::max();
    cfg.trkSel.z0Min     = -std::numeric_limits<float>::max();
    cfg.trkSel.z0Max     =  std::numeric_limits<float>::max();
    cfg.trkSel.chi2Min   = 0.0f;
    cfg.trkSel.chi2Max   = std::numeric_limits<float>::max();
    cfg.trkSel.nHitsMin  = 0;
    cfg.trkSel.nHolesMax = std::numeric_limits<int>::max();

    return cfg;
}

bool LoadAnalysisConfig(const char* configPath, AnalysisConfig& cfg) {
    if (!configPath || !configPath[0]) {
        printf("ERROR: Empty config path.\n");
        return false;
    }

    TEnv env;
    if (env.ReadFile(configPath, kEnvLocal) != 0) {
        printf("ERROR: Failed to read config file: %s\n", configPath);
        return false;
    }

    cfg.inputFilePrefix = env.GetValue("io.inputFilePrefix", cfg.inputFilePrefix.c_str());
    cfg.outputDir = env.GetValue("io.outputDir", cfg.outputDir.c_str());
    cfg.nThreads = env.GetValue("io.nThreads", cfg.nThreads);
    cfg.nThreadsTracks = env.GetValue("io.nThreadsTracks", cfg.nThreadsTracks);
    cfg.nThreadsSeeds = env.GetValue("io.nThreadsSeeds", cfg.nThreadsSeeds);
    cfg.nThreadsHits = env.GetValue("io.nThreadsHits", cfg.nThreadsHits);
    cfg.trackBranch = env.GetValue("io.trackBranch", cfg.trackBranch.c_str());
    cfg.isSubsetCollection = env.GetValue("io.isSubsetCollection", cfg.isSubsetCollection ? 1 : 0) != 0;
    cfg.seedBranch = env.GetValue("io.seedBranch", cfg.seedBranch.c_str());

    cfg.evtSel.ptMin = env.GetValue("event.ptMin", cfg.evtSel.ptMin);
    cfg.evtSel.ptMax = env.GetValue("event.ptMax", cfg.evtSel.ptMax);
    cfg.evtSel.thetaMin = env.GetValue("event.thetaMin", cfg.evtSel.thetaMin);
    cfg.evtSel.thetaMax = env.GetValue("event.thetaMax", cfg.evtSel.thetaMax);
    cfg.evtSel.absEtaMin = env.GetValue("event.absEtaMin", cfg.evtSel.absEtaMin);
    cfg.evtSel.absEtaMax = env.GetValue("event.absEtaMax", cfg.evtSel.absEtaMax);

    cfg.trkSel.ptMin = env.GetValue("track.ptMin", cfg.trkSel.ptMin);
    cfg.trkSel.ptMax = env.GetValue("track.ptMax", cfg.trkSel.ptMax);
    cfg.trkSel.thetaMin = env.GetValue("track.thetaMin", cfg.trkSel.thetaMin);
    cfg.trkSel.thetaMax = env.GetValue("track.thetaMax", cfg.trkSel.thetaMax);
    cfg.trkSel.absEtaMin = env.GetValue("track.absEtaMin", cfg.trkSel.absEtaMin);
    cfg.trkSel.absEtaMax = env.GetValue("track.absEtaMax", cfg.trkSel.absEtaMax);
    cfg.trkSel.phiMin = env.GetValue("track.phiMin", cfg.trkSel.phiMin);
    cfg.trkSel.phiMax = env.GetValue("track.phiMax", cfg.trkSel.phiMax);
    cfg.trkSel.d0Min = env.GetValue("track.d0Min", cfg.trkSel.d0Min);
    cfg.trkSel.d0Max = env.GetValue("track.d0Max", cfg.trkSel.d0Max);
    cfg.trkSel.z0Min = env.GetValue("track.z0Min", cfg.trkSel.z0Min);
    cfg.trkSel.z0Max = env.GetValue("track.z0Max", cfg.trkSel.z0Max);
    cfg.trkSel.chi2Min = env.GetValue("track.chi2Min", cfg.trkSel.chi2Min);
    cfg.trkSel.chi2Max = env.GetValue("track.chi2Max", cfg.trkSel.chi2Max);
    cfg.trkSel.nHitsMin = env.GetValue("track.nHitsMin", cfg.trkSel.nHitsMin);
    cfg.trkSel.nHolesMax = env.GetValue("track.nHolesMax", cfg.trkSel.nHolesMax);

    if (cfg.inputFilePrefix.empty()) {
        printf("ERROR: io.inputFilePrefix is empty in %s\n", configPath);
        return false;
    }

    return true;
}

void PrintLoadedConfig(const AnalysisConfig& cfg, const char* configPath) {
    printf("============================================\n");
    printf("Loaded analysis config: %s\n", configPath);
    printf("============================================\n");
    printf("io.inputFilePrefix     = %s\n", cfg.inputFilePrefix.c_str());
    printf("io.outputDir           = %s\n", cfg.outputDir.c_str());
    printf("io.nThreads            = %d\n", cfg.nThreads);
    printf("io.nThreadsTracks      = %d\n", cfg.nThreadsTracks);
    printf("io.nThreadsSeeds       = %d\n", cfg.nThreadsSeeds);
    printf("io.nThreadsHits        = %d\n", cfg.nThreadsHits);
    printf("io.trackBranch         = %s\n", cfg.trackBranch.c_str());
    printf("io.isSubsetCollection  = %d\n", cfg.isSubsetCollection ? 1 : 0);
    printf("io.seedBranch          = %s\n", cfg.seedBranch.c_str());

    printf("event.ptMin            = %g\n", cfg.evtSel.ptMin);
    printf("event.ptMax            = %g\n", cfg.evtSel.ptMax);
    printf("event.thetaMin         = %g\n", cfg.evtSel.thetaMin);
    printf("event.thetaMax         = %g\n", cfg.evtSel.thetaMax);
    printf("event.absEtaMin        = %g\n", cfg.evtSel.absEtaMin);
    printf("event.absEtaMax        = %g\n", cfg.evtSel.absEtaMax);

    printf("track.ptMin            = %g\n", cfg.trkSel.ptMin);
    printf("track.ptMax            = %g\n", cfg.trkSel.ptMax);
    printf("track.thetaMin         = %g\n", cfg.trkSel.thetaMin);
    printf("track.thetaMax         = %g\n", cfg.trkSel.thetaMax);
    printf("track.absEtaMin        = %g\n", cfg.trkSel.absEtaMin);
    printf("track.absEtaMax        = %g\n", cfg.trkSel.absEtaMax);
    printf("track.phiMin           = %g\n", cfg.trkSel.phiMin);
    printf("track.phiMax           = %g\n", cfg.trkSel.phiMax);
    printf("track.d0Min            = %g\n", cfg.trkSel.d0Min);
    printf("track.d0Max            = %g\n", cfg.trkSel.d0Max);
    printf("track.z0Min            = %g\n", cfg.trkSel.z0Min);
    printf("track.z0Max            = %g\n", cfg.trkSel.z0Max);
    printf("track.chi2Min          = %g\n", cfg.trkSel.chi2Min);
    printf("track.chi2Max          = %g\n", cfg.trkSel.chi2Max);
    printf("track.nHitsMin         = %d\n", cfg.trkSel.nHitsMin);
    printf("track.nHolesMax        = %d\n", cfg.trkSel.nHolesMax);
    printf("============================================\n\n");
}

void RunAnalysis(const char* configPath = "RunAnalysis.conf") {
    AnalysisConfig cfg = DefaultAnalysisConfig();
    if (!LoadAnalysisConfig(configPath, cfg)) {
        return;
    }

    PrintLoadedConfig(cfg, configPath);

    WriteAllMT(cfg.inputFilePrefix.c_str(), cfg.outputDir.c_str(), cfg.nThreads,
               cfg.trackBranch.c_str(), cfg.isSubsetCollection, cfg.seedBranch.c_str(),
               cfg.evtSel, cfg.trkSel,
               cfg.nThreadsTracks, cfg.nThreadsSeeds, cfg.nThreadsHits);
}
