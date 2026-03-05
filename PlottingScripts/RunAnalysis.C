// ============================================================
//  RunAnalysis.C  –  Steering file for NTuple writing
//
//  Edit the sections below to configure:
//    ① Input / output paths and processing options
//    ② Event-level selection (filter events by MC particle)
//    ③ Track-level selection (filter reconstructed tracks)
//
//  Run with:
//    root -l -q RunAnalysis.C
//  or interactively:
//    root -l RunAnalysis.C
// ============================================================

#include "WriteAllMT.C"   // pulls in all sub-writers + SelectionConfig.h

void RunAnalysis() {

    // ──────────────────────────────────────────────────────────
    //  ① INPUT / OUTPUT  — edit these paths
    // ──────────────────────────────────────────────────────────

    // Path prefix that glob-expands to your reco files.
    // All files matching  <inputFilePrefix>*.root  are processed.
    const char* inputFilePrefix = "/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/MC/Mu_pgun/reco_mu-_pt10-5000_theta10-170_3";

    // Directory where tracks_ntuple.root, seeds_ntuple.root,
    // and hits_ntuple.root will be written.
    const char* outputDir = "/global/cfs/cdirs/m5197/sferrar2/TrackingPaper/test";

    // Number of parallel threads.
    int nThreads = 11;

    // Name of the final track collection in the EDM4hep file.
    const char* trackBranch = "SiTrack";

    // Is the track collection a subset (index) collection?
    bool isSubsetCollection = true;

    // Name of the seed-track collection.
    const char* seedBranch = "SeedTracks";


    // ──────────────────────────────────────────────────────────
    //  ② EVENT SELECTION
    //
    //  An event is kept if at least one primary MC particle
    //  (generatorStatus==1, charged, stable through tracker)
    //  satisfies ALL of the pT, theta, and eta windows.
    //
    //  Use etaMin/etaMax OR thetaMin/thetaMax for angular cuts
    //  — applying both simultaneously will double-cut the same
    //  angle and may have unintended effects.
    //
    //  Preset examples:
    //    Full acceptance (no cut):  all fields at defaults
    //    Barrel |η|<0.8:            etaMin=-0.8, etaMax=0.8
    //    Forward region η>1.5:      etaMin=1.5,  etaMax=FLT_MAX
    //    High-pT barrel:            etaMin=-0.8, etaMax=0.8, ptMin=5
    // ──────────────────────────────────────────────────────────

    EventSelectionConfig evtSel;

    // --- Pseudorapidity window (symmetric about η=0) ----------
    // Barrel only (|η| < 0.8 corresponds to θ ∈ [0.72, 2.42] rad)
    evtSel.etaMin = -std::numeric_limits<float>::max(); // no cut
    evtSel.etaMax =  std::numeric_limits<float>::max(); // no cut

    // --- Polar-angle window [radians] -------------------------
    // Use this instead of etaMin/etaMax if you prefer theta.
    evtSel.thetaMin = 0.0f;        // rad  (0   = no lower cut)
    evtSel.thetaMax = (float)M_PI; // rad  (π   = no upper cut)

    // --- Transverse-momentum window [GeV] ---------------------
    evtSel.ptMin = 4000.0f;          // GeV  (0   = no lower cut)
    evtSel.ptMax = std::numeric_limits<float>::max(); // GeV (= no upper cut)


    // ──────────────────────────────────────────────────────────
    //  ③ TRACK SELECTION
    //
    //  Applied to every reconstructed track.  A track must pass
    //  ALL enabled cuts to appear in the output ntuples.
    //  Defaults keep everything.
    //
    //  As with the event selection, use etaMin/etaMax OR
    //  thetaMin/thetaMax for angular cuts, not both.
    // ──────────────────────────────────────────────────────────

    TrackSelectionConfig trkSel;

    // --- Transverse momentum [GeV] ----------------------------
    trkSel.ptMin = 0.0f;
    trkSel.ptMax = std::numeric_limits<float>::max();

    // --- Pseudorapidity (recommended for barrel/endcap) -------
    // Barrel |η| < 0.8:  etaMin = -0.8,  etaMax = 0.8
    trkSel.etaMin = -std::numeric_limits<float>::max(); // no cut
    trkSel.etaMax =  std::numeric_limits<float>::max(); // no cut

    // --- Polar angle [rad] (use instead of eta if preferred) --
    trkSel.thetaMin = 0.0f;
    trkSel.thetaMax = (float)M_PI;

    // --- Azimuthal angle (phi from track state) [rad] ---------
    trkSel.phiMin = -(float)M_PI;
    trkSel.phiMax =  (float)M_PI;

    // --- Transverse impact parameter D0 [mm] ------------------
    trkSel.d0Min = -std::numeric_limits<float>::max();
    trkSel.d0Max =  std::numeric_limits<float>::max();

    // --- Longitudinal impact parameter Z0 [mm] ----------------
    trkSel.z0Min = -std::numeric_limits<float>::max();
    trkSel.z0Max =  std::numeric_limits<float>::max();

    // --- Chi-squared per degree of freedom -------------------
    trkSel.chi2Max = std::numeric_limits<float>::max();

    // --- Minimum number of tracker hits ----------------------
    trkSel.nHitsMin = 0;

    // --- Maximum number of holes (missing expected hits) ------
    trkSel.nHolesMax = std::numeric_limits<int>::max();


    // ──────────────────────────────────────────────────────────
    //  Run
    // ──────────────────────────────────────────────────────────
    WriteAllMT(inputFilePrefix, outputDir, nThreads,
               trackBranch, isSubsetCollection, seedBranch,
               evtSel, trkSel);
}
