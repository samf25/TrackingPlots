#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>

// Include the plot style header
#include "SimplePlotTemplate.h"

void PlotTracks(const char* inputFile, const char* outputFile) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    
    if (!inputFile || !outputFile) {
        printf("Usage: PlotTracks(\"tracks.root\", \"output_plots.root\")\n");
        return;
    }
    
    // Open input file with ntuples
    TFile* inFile = TFile::Open(inputFile);
    if (!inFile || inFile->IsZombie()) {
        printf("Error: Could not open input file %s\n", inputFile);
        return;
    }
    
    // Get the ntuples (names match what WriteTracks.C creates)
    TNtuple* truth_ntuple = (TNtuple*)inFile->Get("truth");
    TNtuple* track_ntuple = (TNtuple*)inFile->Get("tracks");
    TNtuple* matched_ntuple = (TNtuple*)inFile->Get("matched");
    
    if (!truth_ntuple || !track_ntuple || !matched_ntuple) {
        printf("Error: Could not find required ntuples in file\n");
        printf("  Looking for: truth, tracks, matched\n");
        inFile->Close();
        return;
    }
    
    printf("Loaded ntuples from %s\n", inputFile);
    printf("  truth: %lld entries\n", truth_ntuple->GetEntries());
    printf("  tracks: %lld entries\n", track_ntuple->GetEntries());
    printf("  matched: %lld entries\n", matched_ntuple->GetEntries());
    
    // Define binning
    const int nPtBins = 20;
    const double ptMin = 0.5;
    const double ptMax = 110.0;
    
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;
    
    // === Create histograms from ntuples ===
    
    // Efficiency histograms
    TH1D* h_allTruths_pt = new TH1D("allTruths_pt", "", nPtBins, ptMin, ptMax);
    TH1D* h_realTruths_pt = new TH1D("realTruths_pt", "", nPtBins, ptMin, ptMax);
    TH1D* h_allTruths_theta = new TH1D("allTruths_theta", "", nThetaBins, thetaMin, thetaMax);
    TH1D* h_realTruths_theta = new TH1D("realTruths_theta", "", nThetaBins, thetaMin, thetaMax);
    
    // Fake rate histograms
    TH1D* h_allTracks = new TH1D("allTracks", "", nPtBins, ptMin, ptMax);
    TH1D* h_realTracks = new TH1D("realTracks", "", nPtBins, ptMin, ptMax);
    TH1D* h_fakeTracks = new TH1D("fakeTracks", "", nPtBins, ptMin, ptMax);
    
    // Track quality histograms
    TH1D* h_allTracks_nHits = new TH1D("allTracks_nHits", "", 20, 0, 20);
    TH1D* h_allTracks_nHoles = new TH1D("allTracks_nHoles", "", 10, 0, 10);
    TH1D* h_allTracks_chi2ndof = new TH1D("allTracks_chi2ndof", "", 50, 0, 10);
    TH1D* h_realTracks_nHits = new TH1D("realTracks_nHits", "", 20, 0, 20);
    TH1D* h_realTracks_nHoles = new TH1D("realTracks_nHoles", "", 10, 0, 10);
    TH1D* h_realTracks_chi2ndof = new TH1D("realTracks_chi2ndof", "", 50, 0, 10);
    
    // Resolution histograms
    TH1D* h_resolutions_q_over_pt = new TH1D("resolutions_q/pt", "", 100, -0.5, 0.5);
    TH1D* h_resolutions_d0 = new TH1D("resolutions_d0", "", 100, -1., 1.);
    TH1D* h_resolutions_z0 = new TH1D("resolutions_z0", "", 100, -1., 1.);
    
    // Event-level histogram
    TH1D* h_numberOfTracks = new TH1D("numberOfTracks", "", 100, 0, 1000);
    
    // Fill efficiency histograms from truth ntuple
    // Truth ntuple has: pt:theta:charge
    // We'll get matched truth info from matched ntuple instead
    Float_t t_pt, t_theta, t_charge;
    truth_ntuple->SetBranchAddress("pt", &t_pt);
    truth_ntuple->SetBranchAddress("theta", &t_theta);
    truth_ntuple->SetBranchAddress("charge", &t_charge);
    
    // All truths go into denominator
    for (Long64_t i = 0; i < truth_ntuple->GetEntries(); i++) {
        truth_ntuple->GetEntry(i);
        h_allTruths_pt->Fill(t_pt);
        h_allTruths_theta->Fill(t_theta);
    }
    
    // Fill track histograms from track ntuple
    // Track ntuple has: pt:nHits:nHoles:chi2ndof:isReal
    Float_t tr_pt, tr_nHits, tr_nHoles, tr_chi2ndof, tr_isReal;
    track_ntuple->SetBranchAddress("pt", &tr_pt);
    track_ntuple->SetBranchAddress("nHits", &tr_nHits);
    track_ntuple->SetBranchAddress("nHoles", &tr_nHoles);
    track_ntuple->SetBranchAddress("chi2ndof", &tr_chi2ndof);
    track_ntuple->SetBranchAddress("isReal", &tr_isReal);
    
    for (Long64_t i = 0; i < track_ntuple->GetEntries(); i++) {
        track_ntuple->GetEntry(i);
        h_allTracks->Fill(tr_pt);
        h_allTracks_nHits->Fill(tr_nHits);
        h_allTracks_nHoles->Fill(tr_nHoles);
        h_allTracks_chi2ndof->Fill(tr_chi2ndof);
        
        if (tr_isReal > 0.5) {
            h_realTracks->Fill(tr_pt);
            h_realTracks_nHits->Fill(tr_nHits);
            h_realTracks_nHoles->Fill(tr_nHoles);
            h_realTracks_chi2ndof->Fill(tr_chi2ndof);
        } else {
            h_fakeTracks->Fill(tr_pt);
        }
    }
    
    // Fill resolution histograms from matched ntuple
    // Matched ntuple has: true_pt:true_theta:true_charge:reco_pt:reco_d0:reco_z0:true_d0:true_z0:nHits:nHoles:chi2ndof
    Float_t m_true_pt, m_true_theta, m_true_charge, m_reco_pt, m_reco_d0, m_reco_z0, m_true_d0, m_true_z0;
    matched_ntuple->SetBranchAddress("true_pt", &m_true_pt);
    matched_ntuple->SetBranchAddress("true_theta", &m_true_theta);
    matched_ntuple->SetBranchAddress("true_charge", &m_true_charge);
    matched_ntuple->SetBranchAddress("reco_pt", &m_reco_pt);
    matched_ntuple->SetBranchAddress("reco_d0", &m_reco_d0);
    matched_ntuple->SetBranchAddress("reco_z0", &m_reco_z0);
    matched_ntuple->SetBranchAddress("true_d0", &m_true_d0);
    matched_ntuple->SetBranchAddress("true_z0", &m_true_z0);
    
    for (Long64_t i = 0; i < matched_ntuple->GetEntries(); i++) {
        matched_ntuple->GetEntry(i);
        // Fill numerator for efficiency (matched truths)
        h_realTruths_pt->Fill(m_true_pt);
        h_realTruths_theta->Fill(m_true_theta);
        // Calculate resolutions
        double q_over_pt_res = (m_true_charge/m_reco_pt - m_true_charge/m_true_pt) / (m_true_charge/m_true_pt);
        double d0_res = m_reco_d0 - m_true_d0;
        double z0_res = m_reco_z0 - m_true_z0;
        h_resolutions_q_over_pt->Fill(q_over_pt_res);
        h_resolutions_d0->Fill(d0_res);
        h_resolutions_z0->Fill(z0_res);
    }
    
    // Fill event-level histogram from events_summary ntuple
    TNtuple* events_ntuple = (TNtuple*)inFile->Get("events_summary");
    if (events_ntuple) {
        Float_t e_nTracks;
        events_ntuple->SetBranchAddress("nTracks", &e_nTracks);
        for (Long64_t i = 0; i < events_ntuple->GetEntries(); i++) {
            events_ntuple->GetEntry(i);
            h_numberOfTracks->Fill(e_nTracks);
        }
    }
    
    // === Create plots ===
    
    // Create efficiency plots
    TEfficiency* eff_pt = new TEfficiency(*h_realTruths_pt, *h_allTruths_pt);
    TEfficiency* eff_theta = new TEfficiency(*h_realTruths_theta, *h_allTruths_theta);
    TEfficiency* fake_rate = new TEfficiency(*h_fakeTracks, *h_allTracks);
    
    // Create output file
    TFile* outFile = new TFile(outputFile, "RECREATE");
    
    // === Efficiency vs pT ===
    TDirectory* effDir = outFile->mkdir("Efficiency");
    effDir->cd();
    
    TCanvas* c_eff_pt = MuCollStyle::CreateCanvas("c_eff_pt", "Tracking Efficiency vs p_{T}");
    eff_pt->SetTitle(";p_{T} [GeV];Efficiency");
    eff_pt->Draw("APE");
    MuCollStyle::AddStandardLabels(c_eff_pt, "10 TeV");
    c_eff_pt->Write();
    
    // === Efficiency vs theta ===
    TCanvas* c_eff_theta = MuCollStyle::CreateCanvas("c_eff_theta", "Tracking Efficiency vs #theta");
    eff_theta->SetTitle(";#theta [rad];Efficiency");
    eff_theta->Draw("APE");
    MuCollStyle::AddStandardLabels(c_eff_theta, "10 TeV");
    c_eff_theta->Write();
    
    // === Fake rate vs pT ===
    TDirectory* fakeDir = outFile->mkdir("FakeRate");
    fakeDir->cd();
    
    TCanvas* c_fake = MuCollStyle::CreateCanvas("c_fake_rate", "Fake Rate vs p_{T}");
    fake_rate->SetTitle(";p_{T} [GeV];Fake Rate");
    fake_rate->Draw("APE");
    MuCollStyle::AddStandardLabels(c_fake, "10 TeV");
    c_fake->Write();
    
    // === Resolution plots ===
    TDirectory* resDir = outFile->mkdir("Resolutions");
    resDir->cd();
    
    TCanvas* c_res_qpt = MuCollStyle::CreateCanvas("c_res_qpt", "q/p_{T} Resolution");
    MuCollStyle::StyleHist(h_resolutions_q_over_pt, MuCollStyle::GetColor(0));
    h_resolutions_q_over_pt->GetXaxis()->SetTitle("#Delta(q/p_{T}) / (q/p_{T})");
    h_resolutions_q_over_pt->GetYaxis()->SetTitle("Entries");
    h_resolutions_q_over_pt->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_qpt, "10 TeV");
    c_res_qpt->Write();
    
    TCanvas* c_res_d0 = MuCollStyle::CreateCanvas("c_res_d0", "d_{0} Resolution");
    MuCollStyle::StyleHist(h_resolutions_d0, MuCollStyle::GetColor(1));
    h_resolutions_d0->GetXaxis()->SetTitle("#Delta d_{0} [mm]");
    h_resolutions_d0->GetYaxis()->SetTitle("Entries");
    h_resolutions_d0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_d0, "10 TeV");
    c_res_d0->Write();
    
    TCanvas* c_res_z0 = MuCollStyle::CreateCanvas("c_res_z0", "z_{0} Resolution");
    MuCollStyle::StyleHist(h_resolutions_z0, MuCollStyle::GetColor(2));
    h_resolutions_z0->GetXaxis()->SetTitle("#Delta z_{0} [mm]");
    h_resolutions_z0->GetYaxis()->SetTitle("Entries");
    h_resolutions_z0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_z0, "10 TeV");
    c_res_z0->Write();
    
    // === Track quality plots ===
    TDirectory* qualDir = outFile->mkdir("TrackQuality");
    qualDir->cd();
    
    TCanvas* c_nHits = MuCollStyle::CreateCanvas("c_nHits", "Number of Hits per Track");
    MuCollStyle::StyleHist(h_allTracks_nHits, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_nHits, MuCollStyle::GetColor(1));
    h_allTracks_nHits->GetXaxis()->SetTitle("Number of Hits");
    h_allTracks_nHits->GetYaxis()->SetTitle("Entries");
    h_allTracks_nHits->Draw("HIST");
    h_realTracks_nHits->Draw("HIST SAME");
    TLegend* leg_nHits = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg_nHits->AddEntry(h_allTracks_nHits, "All Tracks", "l");
    leg_nHits->AddEntry(h_realTracks_nHits, "Real Tracks", "l");
    leg_nHits->Draw();
    MuCollStyle::AddStandardLabels(c_nHits, "10 TeV");
    c_nHits->Write();
    
    TCanvas* c_nHoles = MuCollStyle::CreateCanvas("c_nHoles", "Number of Holes per Track");
    MuCollStyle::StyleHist(h_allTracks_nHoles, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_nHoles, MuCollStyle::GetColor(1));
    h_allTracks_nHoles->GetXaxis()->SetTitle("Number of Holes");
    h_allTracks_nHoles->GetYaxis()->SetTitle("Entries");
    h_allTracks_nHoles->Draw("HIST");
    h_realTracks_nHoles->Draw("HIST SAME");
    TLegend* leg_nHoles = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg_nHoles->AddEntry(h_allTracks_nHoles, "All Tracks", "l");
    leg_nHoles->AddEntry(h_realTracks_nHoles, "Real Tracks", "l");
    leg_nHoles->Draw();
    MuCollStyle::AddStandardLabels(c_nHoles, "10 TeV");
    c_nHoles->Write();
    
    TCanvas* c_chi2 = MuCollStyle::CreateCanvas("c_chi2ndof", "#chi^{2}/ndof");
    MuCollStyle::StyleHist(h_allTracks_chi2ndof, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_chi2ndof, MuCollStyle::GetColor(1));
    h_allTracks_chi2ndof->GetXaxis()->SetTitle("#chi^{2}/ndof");
    h_allTracks_chi2ndof->GetYaxis()->SetTitle("Entries");
    h_allTracks_chi2ndof->Draw("HIST");
    h_realTracks_chi2ndof->Draw("HIST SAME");
    TLegend* leg_chi2 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg_chi2->AddEntry(h_allTracks_chi2ndof, "All Tracks", "l");
    leg_chi2->AddEntry(h_realTracks_chi2ndof, "Real Tracks", "l");
    leg_chi2->Draw();
    MuCollStyle::AddStandardLabels(c_chi2, "10 TeV");
    c_chi2->Write();
    
    // Number of tracks per event
    TCanvas* c_nTracks = MuCollStyle::CreateCanvas("c_numberOfTracks", "Tracks per Event");
    MuCollStyle::StyleHist(h_numberOfTracks, MuCollStyle::GetColor(0));
    h_numberOfTracks->GetXaxis()->SetTitle("Number of Tracks");
    h_numberOfTracks->GetYaxis()->SetTitle("Events");
    h_numberOfTracks->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_nTracks, "10 TeV");
    c_nTracks->Write();
    
    // Write raw objects
    outFile->cd();
    eff_pt->Write("eff_pt");
    eff_theta->Write("eff_theta");
    fake_rate->Write("fake_rate");
    
    outFile->Write();
    outFile->Close();
    inFile->Close();
    
    printf("Plots written to %s\n", outputFile);
}
