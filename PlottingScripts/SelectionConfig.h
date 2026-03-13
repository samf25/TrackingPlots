#pragma once
// ============================================================
//  SelectionConfig.h
//
//  Shared configuration structs and inline helpers for
//  event-level and track-level selections applied during
//  NTuple writing.
//
//  EventSelectionConfig: keep an event if at least one
//    primary MC particle (generatorStatus==1, charged, not
//    created in simulation, not decayed in tracker) falls
//    within the requested pT, theta, AND |eta| windows.
//
//  TrackSelectionConfig: keep a reconstructed track if all
//    specified kinematic and quality cuts are satisfied.
//
//  Theta / eta relationship:
//    eta = -ln( tan(theta/2) )
//    theta = 2 * atan( exp(-eta) )
//
//  Endcap selection example (|eta| > 0.7):
//    cfg.absEtaMin = 0.7;
//    (thetaMin / thetaMax can be left at defaults)
//
//  NOTE: theta and |eta| cuts are both applied. Set only one
//  set of fields to avoid double-counting the same cut.
// ============================================================

#include <cmath>
#include <limits>

// ─── Simulator-status bit indices (EDM4hep / LCIO convention) ─
static const int BITCreatedInSimulation = 30;
static const int BITDecayedInTracker    = 27;

// ─── Event-level selection ────────────────────────────────────
struct EventSelectionConfig {
    // Retain the event if at least one acceptable primary MC
    // particle satisfies ALL of the pT, theta, and |eta|
    // requirements simultaneously.
    // Defaults cover the full phase space (no cut).

    // --- pT window [GeV] --------------------------------------
    float ptMin    =  0.0f;
    float ptMax    =  std::numeric_limits<float>::max();

    // --- Polar-angle window [rad] (0 = beam, π/2 = barrel) ----
    // Use thetaMin/thetaMax OR absEtaMin/absEtaMax, not both.
    float thetaMin =  0.0f;
    float thetaMax =  (float)M_PI;

    // --- Absolute pseudorapidity window ----------------------
    // |η| = |-ln(tan(θ/2))|.  Endcap |η|>0.7 → absEtaMin=0.7
    // Barrel |η|<0.8 → absEtaMax=0.8.  Default: no cut.
    float absEtaMin =  0.0f;
    float absEtaMax =  std::numeric_limits<float>::max();
};

// ─── Track-level selection ────────────────────────────────────
struct TrackSelectionConfig {
    // All cut variables are applied simultaneously.
    // Defaults pass everything.

    // --- Transverse momentum [GeV] ----------------------------
    float ptMin     =  0.0f;
    float ptMax     =  std::numeric_limits<float>::max();

    // --- Polar angle [rad] ------------------------------------
    // Use thetaMin/thetaMax OR absEtaMin/absEtaMax, not both.
    float thetaMin  =  0.0f;
    float thetaMax  =  (float)M_PI;

    // --- Absolute pseudorapidity ------------------------------
    // |η| = |-ln(tan(θ/2))|.  Endcap |η|>0.7 → absEtaMin=0.7
    // Barrel |η|<0.8 → absEtaMax=0.8.  Default: no cut.
    float absEtaMin =  0.0f;
    float absEtaMax =  std::numeric_limits<float>::max();

    // --- Azimuthal angle [rad] --------------------------------
    float phiMin    = -(float)M_PI;
    float phiMax    =  (float)M_PI;

    // --- Impact parameters [mm] -------------------------------
    float d0Min     = -std::numeric_limits<float>::max();
    float d0Max     =  std::numeric_limits<float>::max();

    float z0Min     = -std::numeric_limits<float>::max();
    float z0Max     =  std::numeric_limits<float>::max();

    // --- Track quality ----------------------------------------
    float chi2Max   =  std::numeric_limits<float>::max();  // chi2/ndf
    int   nHitsMin  =  0;
    int   nHolesMax =  std::numeric_limits<int>::max();
};

// ─── Inline helpers ───────────────────────────────────────────

/// Convert polar angle theta [rad] to pseudorapidity eta.
/// theta = 0 or π gives ±infinity; safe for finite values.
inline float ThetaToEta(float theta) {
    // Clamp to avoid exact 0 or π which give ±inf
    const float eps = 1e-6f;
    float t = std::max(eps, std::min((float)M_PI - eps, theta));
    return -std::log(std::tan(t * 0.5f));
}

/// Returns true if a primary MC particle with the given pT and
/// theta satisfies the event-level selection (pT, theta, |eta|
/// windows are all applied simultaneously).
inline bool MCPassesEventSelection(float pt, float theta,
                                   const EventSelectionConfig& cfg) {
    if (pt    < cfg.ptMin    || pt    > cfg.ptMax)    return false;
    if (theta < cfg.thetaMin || theta > cfg.thetaMax) return false;
    float absEta = std::abs(ThetaToEta(theta));
    if (absEta < cfg.absEtaMin || absEta > cfg.absEtaMax) return false;
    return true;
}

/// Returns true if a reconstructed track passes all track cuts.
/// phi     – azimuthal angle of momentum from track state [rad]
/// theta   – polar angle of momentum (from tanLambda) [rad]
/// pt      – transverse momentum [GeV]
/// d0, z0  – impact parameters [mm]
/// chi2ndof – chi2 per degree of freedom (pass -1 to skip)
/// nHits, nHoles – from TrackData
inline bool TrackPassesSelection(float pt, float theta, float phi,
                                 float d0, float z0, float chi2ndof,
                                 int nHits, int nHoles,
                                 const TrackSelectionConfig& cfg) {
    if (pt    < cfg.ptMin    || pt    > cfg.ptMax)    return false;
    if (theta < cfg.thetaMin || theta > cfg.thetaMax) return false;
    float absEta = std::abs(ThetaToEta(theta));
    if (absEta < cfg.absEtaMin || absEta > cfg.absEtaMax) return false;
    if (phi   < cfg.phiMin   || phi   > cfg.phiMax)   return false;
    if (d0    < cfg.d0Min    || d0    > cfg.d0Max)     return false;
    if (z0    < cfg.z0Min    || z0    > cfg.z0Max)     return false;
    if (chi2ndof >= 0 && chi2ndof > cfg.chi2Max)       return false;
    if (nHits  < cfg.nHitsMin)                         return false;
    if (nHoles > cfg.nHolesMax)                        return false;
    return true;
}

/// Returns true when any event-level cut is active (i.e. not
/// all fields are at their full-phase-space defaults).  Used to
/// skip the MC scan entirely when no event selection is needed.
inline bool EventSelectionIsActive(const EventSelectionConfig& cfg) {
    const float fmax = std::numeric_limits<float>::max();
    return !(cfg.ptMin    == 0.0f  && cfg.ptMax    == fmax &&
             cfg.thetaMin == 0.0f  && cfg.thetaMax == (float)M_PI &&
             cfg.absEtaMin == 0.0f && cfg.absEtaMax == fmax);
}
