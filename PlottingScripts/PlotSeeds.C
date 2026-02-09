#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <algorithm>

// Include the plot style header
#include "SimplePlotTemplate.h"

void PlotSeeds(const char* inputFile, const char* outputFile) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    
    if (!inputFile || !outputFile) {
        printf("Usage: PlotSeeds(\"seeds.root\", \"output_plots.root\")\n");
        return;
    }
    
    // Open input file with ntuples
    TFile* inFile = TFile::Open(inputFile);
    if (!inFile || inFile->IsZombie()) {
        printf("Error: Could not open input file %s\n", inputFile);
        return;
    }
    
    // Get the ntuples (names match what WriteSeeds.C creates)
    TNtuple* seed_ntuple = (TNtuple*)inFile->Get("seeds");
    TNtuple* layer_ntuple = (TNtuple*)inFile->Get("seed_layers");
    TNtuple* matched_seed_ntuple = (TNtuple*)inFile->Get("matched_seeds");
    
    if (!seed_ntuple || !layer_ntuple || !matched_seed_ntuple) {
        printf("Error: Could not find required ntuples in file\n");
        printf("  Looking for: seeds, seed_layers, matched_seeds\n");
        inFile->Close();
        return;
    }
    
    printf("Loaded ntuples from %s\n", inputFile);
    printf("  seeds: %lld entries\n", seed_ntuple->GetEntries());
    printf("  seed_layers: %lld entries\n", layer_ntuple->GetEntries());
    printf("  matched_seeds: %lld entries\n", matched_seed_ntuple->GetEntries());
    
    // Define binning
    const int nPtBins = 20;
    const double ptMin = 0.5;
    const double ptMax = 110.0;
    
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;
    
    // === Create histograms ===
    
    // Seed theta histogram
    TH1D* h_seed_theta = new TH1D("seed_theta", "", nThetaBins, thetaMin, thetaMax);
    
    // Layer histograms - separate barrel and endcap, all/matched/unmatched
    TH1D* h_seed_layer_barrel = new TH1D("seed_layer_barrel", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap = new TH1D("seed_layer_endcap", "", 15, 0, 15);
    TH1D* h_seed_layer_barrel_matched = new TH1D("seed_layer_barrel_matched", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap_matched = new TH1D("seed_layer_endcap_matched", "", 15, 0, 15);
    TH1D* h_seed_layer_barrel_unmatched = new TH1D("seed_layer_barrel_unmatched", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap_unmatched = new TH1D("seed_layer_endcap_unmatched", "", 15, 0, 15);
    
    // Event-level histograms
    std::vector<double> seed_number_vec;
    std::vector<double> seed_matched_vec;
    std::vector<double> seed_unmatched_vec;
    TH1D* h_seeds_per_MCP = new TH1D("seeds_per_MCP", "", 10, 0, 10);
    
    // Resolution histograms
    TH1D* h_realSeeds_pt = new TH1D("realSeeds_pt", "", nPtBins, ptMin, ptMax);
    TH1D* h_resolutions_q_over_pt = new TH1D("seed_resolutions_q/pt", "", 100, -0.5, 0.5);
    TH1D* h_resolutions_d0 = new TH1D("seed_resolutions_d0", "", 100, -1., 1.);
    TH1D* h_resolutions_z0 = new TH1D("seed_resolutions_z0", "", 100, -1., 1.);
    
    // === Fill histograms from ntuples ===
    
    // Seed ntuple has: theta:isMatched
    Float_t s_theta, s_isMatched;
    seed_ntuple->SetBranchAddress("theta", &s_theta);
    seed_ntuple->SetBranchAddress("isMatched", &s_isMatched);
    
    for (Long64_t i = 0; i < seed_ntuple->GetEntries(); i++) {
        seed_ntuple->GetEntry(i);
        h_seed_theta->Fill(s_theta);
    }
    
    // Layer ntuple has: layer:isBarrel:isMatched
    Float_t h_layer, h_isBarrel, h_layerIsMatched;
    layer_ntuple->SetBranchAddress("layer", &h_layer);
    layer_ntuple->SetBranchAddress("isBarrel", &h_isBarrel);
    layer_ntuple->SetBranchAddress("isMatched", &h_layerIsMatched);
    
    for (Long64_t i = 0; i < layer_ntuple->GetEntries(); i++) {
        layer_ntuple->GetEntry(i);
        
        if (h_isBarrel > 0.5) {
            h_seed_layer_barrel->Fill(h_layer);
            if (h_layerIsMatched > 0.5) {
                h_seed_layer_barrel_matched->Fill(h_layer);
            } else {
                h_seed_layer_barrel_unmatched->Fill(h_layer);
            }
        } else {
            h_seed_layer_endcap->Fill(h_layer);
            if (h_layerIsMatched > 0.5) {
                h_seed_layer_endcap_matched->Fill(h_layer);
            } else {
                h_seed_layer_endcap_unmatched->Fill(h_layer);
            }
        }
    }
    
    // Events summary ntuple has: nSeeds:nMatched:nUnmatched:avgSeedsPerMCP
    TNtuple* events_ntuple = (TNtuple*)inFile->Get("events_summary");
    if (events_ntuple) {
        Float_t e_nSeeds, e_nMatched, e_nUnmatched, e_avgSeedsPerMCP;
        events_ntuple->SetBranchAddress("nSeeds", &e_nSeeds);
        events_ntuple->SetBranchAddress("nMatched", &e_nMatched);
        events_ntuple->SetBranchAddress("nUnmatched", &e_nUnmatched);
        events_ntuple->SetBranchAddress("avgSeedsPerMCP", &e_avgSeedsPerMCP);
        
        for (Long64_t i = 0; i < events_ntuple->GetEntries(); i++) {
            events_ntuple->GetEntry(i);
            seed_number_vec.push_back(e_nSeeds);
            seed_matched_vec.push_back(e_nMatched);
            seed_unmatched_vec.push_back(e_nUnmatched);
            if (e_avgSeedsPerMCP > 0) h_seeds_per_MCP->Fill(e_avgSeedsPerMCP);
        }
    }
    
    // Matched seed ntuple has: true_pt:true_theta:reco_pt:reco_d0:reco_z0:true_d0:true_z0:res_q_over_pt:res_d0:res_z0
    Float_t m_true_pt, m_true_theta, m_reco_pt, m_reco_d0, m_reco_z0, m_true_d0, m_true_z0;
    Float_t m_res_q_over_pt, m_res_d0, m_res_z0;
    matched_seed_ntuple->SetBranchAddress("true_pt", &m_true_pt);
    matched_seed_ntuple->SetBranchAddress("true_theta", &m_true_theta);
    matched_seed_ntuple->SetBranchAddress("reco_pt", &m_reco_pt);
    matched_seed_ntuple->SetBranchAddress("reco_d0", &m_reco_d0);
    matched_seed_ntuple->SetBranchAddress("reco_z0", &m_reco_z0);
    matched_seed_ntuple->SetBranchAddress("true_d0", &m_true_d0);
    matched_seed_ntuple->SetBranchAddress("true_z0", &m_true_z0);
    matched_seed_ntuple->SetBranchAddress("res_q_over_pt", &m_res_q_over_pt);
    matched_seed_ntuple->SetBranchAddress("res_d0", &m_res_d0);
    matched_seed_ntuple->SetBranchAddress("res_z0", &m_res_z0);
    
    for (Long64_t i = 0; i < matched_seed_ntuple->GetEntries(); i++) {
        matched_seed_ntuple->GetEntry(i);
        // Fill pt histograms for matched seeds
        h_realSeeds_pt->Fill(m_reco_pt);
        h_resolutions_q_over_pt->Fill(m_res_q_over_pt);
        h_resolutions_d0->Fill(m_res_d0);
        h_resolutions_z0->Fill(m_res_z0);
    }
    
    // === Create plots ===
    
    // Create output file
    TFile* outFile = new TFile(outputFile, "RECREATE");
    
    // === Seed info plots ===
    TDirectory* seedDir = outFile->mkdir("SeedInfo");
    seedDir->cd();
    
    // Helper: adjust canvas margins for right axis
    auto adjust_margins = [](TCanvas* c) {
        c->SetLeftMargin(0.10);  // Reduce left margin
        c->SetRightMargin(0.18); // Increase right margin for right axis
        c->SetTicky(0);          // Disable right-side tick marks (we have a custom axis)
    };
    
    // Seed theta distribution
    TCanvas* c_theta = MuCollStyle::CreateCanvas("c_seed_theta", "Seed #theta Distribution");
    MuCollStyle::StyleHist(h_seed_theta, MuCollStyle::GetColor(0));
    h_seed_theta->GetXaxis()->SetTitle("#theta [rad]");
    h_seed_theta->GetYaxis()->SetTitle("Seeds");
    h_seed_theta->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_theta, "10 TeV");
    c_theta->Write();
    
    // Seeds per event - with matched on right axis
    // Determine histogram range from data
    double maxSeeds = 10;
    if (!seed_number_vec.empty()) {
        maxSeeds = *std::max_element(seed_number_vec.begin(), seed_number_vec.end()) * 1.1;
        if (maxSeeds < 10) maxSeeds = 10;
    }
    TH1D* h_seed_number = new TH1D("seed_number", "", 50, 0, maxSeeds);
    TH1D* h_seed_matched = new TH1D("seed_matched", "", 50, 0, maxSeeds);
    TH1D* h_seed_unmatched = new TH1D("seed_unmatched", "", 50, 0, maxSeeds);

    for (size_t i = 0; i < seed_number_vec.size(); ++i) {
        h_seed_number->Fill(seed_number_vec[i]);
    }
    for (size_t i = 0; i < seed_matched_vec.size(); ++i) {
        h_seed_matched->Fill(seed_matched_vec[i]);
    }
    for (size_t i = 0; i < seed_unmatched_vec.size(); ++i) {
        h_seed_unmatched->Fill(seed_unmatched_vec[i]);
    }
    double scale_factor = 1.0;
    if (h_seed_number->GetMaximum() > 0 && h_seed_matched->GetMaximum() > 0) {
        scale_factor = h_seed_number->GetMaximum() / h_seed_matched->GetMaximum();
        if (scale_factor < 5.0) scale_factor = 5.0;
    }
    TH1D* h_seed_matched_scaled = (TH1D*)h_seed_matched->Clone("seed_matched_scaled");
    h_seed_matched_scaled->Scale(scale_factor);
    
    TCanvas* c_nSeeds = MuCollStyle::CreateCanvas("c_nSeeds", "Seeds per Event");
    adjust_margins(c_nSeeds);
    MuCollStyle::StyleHist(h_seed_number, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_seed_matched_scaled, MuCollStyle::GetColor(1));
    MuCollStyle::StyleHist(h_seed_unmatched, MuCollStyle::GetColor(2));
    h_seed_unmatched->SetLineStyle(2);  // Dotted line for unmatched
    h_seed_number->GetXaxis()->SetTitle("Number of Seeds");
    h_seed_number->GetYaxis()->SetTitle("Events (Total/Unmatched)");
    h_seed_number->Draw("HIST");
    h_seed_unmatched->Draw("HIST SAME");
    h_seed_matched_scaled->Draw("HIST SAME");
    TLegend* leg_nSeeds = MuCollStyle::CreateLegend(0.45, 0.65, 0.78, 0.88);
    leg_nSeeds->AddEntry(h_seed_number, "Total Seeds", "l");
    leg_nSeeds->AddEntry(h_seed_matched_scaled, "Matched Seeds", "l");
    leg_nSeeds->AddEntry(h_seed_unmatched, "Unmatched Seeds", "l");
    leg_nSeeds->Draw();
    MuCollStyle::AddStandardLabels(c_nSeeds, "10 TeV");
    c_nSeeds->Update();
    double left_min = 0, left_max = h_seed_number->GetMaximum();
    double right_min = 0, right_max = h_seed_matched->GetMaximum();
    TGaxis* axis_nSeeds = new TGaxis(
        h_seed_number->GetXaxis()->GetXmax(), left_min,
        h_seed_number->GetXaxis()->GetXmax(), left_max,
        right_min, right_max, 510, "+L"
    );
    axis_nSeeds->SetTitle("Events (Matched)");
    axis_nSeeds->SetTitleOffset(1.2);
    axis_nSeeds->SetLineColor(kBlack);
    axis_nSeeds->SetLabelColor(kBlack);
    axis_nSeeds->SetTitleColor(kBlack);
    axis_nSeeds->SetLabelFont(h_seed_number->GetYaxis()->GetLabelFont());
    axis_nSeeds->SetLabelSize(h_seed_number->GetYaxis()->GetLabelSize());
    axis_nSeeds->SetTitleFont(h_seed_number->GetYaxis()->GetTitleFont());
    axis_nSeeds->SetTitleSize(h_seed_number->GetYaxis()->GetTitleSize());
    axis_nSeeds->Draw();
    c_nSeeds->Write();
    
    // Seeds per MCP
    TCanvas* c_perMCP = MuCollStyle::CreateCanvas("c_seeds_per_MCP", "Seeds per MC Particle");
    MuCollStyle::StyleHist(h_seeds_per_MCP, MuCollStyle::GetColor(0));
    h_seeds_per_MCP->GetXaxis()->SetTitle("Average Seeds per MCP");
    h_seeds_per_MCP->GetYaxis()->SetTitle("Events");
    h_seeds_per_MCP->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_perMCP, "10 TeV");
    c_perMCP->Write();
    
    // === Layer distribution plots ===
    TDirectory* layerDir = outFile->mkdir("LayerDistribution");
    layerDir->cd();
    
    // Barrel layers - with matched on right axis
    scale_factor = 1.0;
    if (h_seed_layer_barrel->GetMaximum() > 0 && h_seed_layer_barrel_matched->GetMaximum() > 0) {
        scale_factor = h_seed_layer_barrel->GetMaximum() / h_seed_layer_barrel_matched->GetMaximum();
        if (scale_factor < 5.0) scale_factor = 5.0;
    }
    TH1D* h_seed_layer_barrel_matched_scaled = (TH1D*)h_seed_layer_barrel_matched->Clone("seed_layer_barrel_matched_scaled");
    h_seed_layer_barrel_matched_scaled->Scale(scale_factor);
    
    TCanvas* c_barrel = MuCollStyle::CreateCanvas("c_barrel_layers", "Barrel Layer Distribution");
    adjust_margins(c_barrel);
    MuCollStyle::StyleHist(h_seed_layer_barrel, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_seed_layer_barrel_matched_scaled, MuCollStyle::GetColor(1));
    MuCollStyle::StyleHist(h_seed_layer_barrel_unmatched, MuCollStyle::GetColor(2));
    h_seed_layer_barrel_unmatched->SetLineStyle(2);  // Dotted line for unmatched
    h_seed_layer_barrel->GetXaxis()->SetTitle("Layer");
    h_seed_layer_barrel->GetYaxis()->SetTitle("Hits (All/Unmatched)");
    h_seed_layer_barrel->Draw("HIST");
    h_seed_layer_barrel_unmatched->Draw("HIST SAME");
    h_seed_layer_barrel_matched_scaled->Draw("HIST SAME");
    TLegend* leg_barrel = MuCollStyle::CreateLegend(0.45, 0.65, 0.78, 0.88);
    leg_barrel->AddEntry(h_seed_layer_barrel, "All", "l");
    leg_barrel->AddEntry(h_seed_layer_barrel_matched_scaled, "Matched", "l");
    leg_barrel->AddEntry(h_seed_layer_barrel_unmatched, "Unmatched", "l");
    leg_barrel->Draw();
    MuCollStyle::AddStandardLabels(c_barrel, "10 TeV");
    c_barrel->Update();
    left_max = h_seed_layer_barrel->GetMaximum();
    right_max = h_seed_layer_barrel_matched->GetMaximum();
    TGaxis* axis_barrel = new TGaxis(
        h_seed_layer_barrel->GetXaxis()->GetXmax(), left_min,
        h_seed_layer_barrel->GetXaxis()->GetXmax(), left_max,
        right_min, right_max, 510, "+L"
    );
    axis_barrel->SetTitle("Hits (Matched)");
    axis_barrel->SetTitleOffset(1.2);
    axis_barrel->SetLineColor(kBlack);
    axis_barrel->SetLabelColor(kBlack);
    axis_barrel->SetTitleColor(kBlack);
    axis_barrel->SetLabelFont(h_seed_layer_barrel->GetYaxis()->GetLabelFont());
    axis_barrel->SetLabelSize(h_seed_layer_barrel->GetYaxis()->GetLabelSize());
    axis_barrel->SetTitleFont(h_seed_layer_barrel->GetYaxis()->GetTitleFont());
    axis_barrel->SetTitleSize(h_seed_layer_barrel->GetYaxis()->GetTitleSize());
    axis_barrel->Draw();
    c_barrel->Write();
    
    // Endcap layers - with matched on right axis
    scale_factor = 1.0;
    if (h_seed_layer_endcap->GetMaximum() > 0 && h_seed_layer_endcap_matched->GetMaximum() > 0) {
        scale_factor = h_seed_layer_endcap->GetMaximum() / h_seed_layer_endcap_matched->GetMaximum();
        if (scale_factor < 5.0) scale_factor = 5.0;
    }
    TH1D* h_seed_layer_endcap_matched_scaled = (TH1D*)h_seed_layer_endcap_matched->Clone("seed_layer_endcap_matched_scaled");
    h_seed_layer_endcap_matched_scaled->Scale(scale_factor);
    
    TCanvas* c_endcap = MuCollStyle::CreateCanvas("c_endcap_layers", "Endcap Layer Distribution");
    adjust_margins(c_endcap);
    MuCollStyle::StyleHist(h_seed_layer_endcap, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_seed_layer_endcap_matched_scaled, MuCollStyle::GetColor(1));
    MuCollStyle::StyleHist(h_seed_layer_endcap_unmatched, MuCollStyle::GetColor(2));
    h_seed_layer_endcap_unmatched->SetLineStyle(2);  // Dotted line for unmatched
    h_seed_layer_endcap->GetXaxis()->SetTitle("Layer");
    h_seed_layer_endcap->GetYaxis()->SetTitle("Hits (All/Unmatched)");
    h_seed_layer_endcap->Draw("HIST");
    h_seed_layer_endcap_unmatched->Draw("HIST SAME");
    h_seed_layer_endcap_matched_scaled->Draw("HIST SAME");
    TLegend* leg_endcap = MuCollStyle::CreateLegend(0.45, 0.65, 0.78, 0.88);
    leg_endcap->AddEntry(h_seed_layer_endcap, "All", "l");
    leg_endcap->AddEntry(h_seed_layer_endcap_matched_scaled, "Matched", "l");
    leg_endcap->AddEntry(h_seed_layer_endcap_unmatched, "Unmatched", "l");
    leg_endcap->Draw();
    MuCollStyle::AddStandardLabels(c_endcap, "10 TeV");
    c_endcap->Update();
    left_max = h_seed_layer_endcap->GetMaximum();
    right_max = h_seed_layer_endcap_matched->GetMaximum();
    TGaxis* axis_endcap = new TGaxis(
        h_seed_layer_endcap->GetXaxis()->GetXmax(), left_min,
        h_seed_layer_endcap->GetXaxis()->GetXmax(), left_max,
        right_min, right_max, 510, "+L"
    );
    axis_endcap->SetTitle("Hits (Matched)");
    axis_endcap->SetTitleOffset(1.2);
    axis_endcap->SetLineColor(kBlack);
    axis_endcap->SetLabelColor(kBlack);
    axis_endcap->SetTitleColor(kBlack);
    axis_endcap->SetLabelFont(h_seed_layer_endcap->GetYaxis()->GetLabelFont());
    axis_endcap->SetLabelSize(h_seed_layer_endcap->GetYaxis()->GetLabelSize());
    axis_endcap->SetTitleFont(h_seed_layer_endcap->GetYaxis()->GetTitleFont());
    axis_endcap->SetTitleSize(h_seed_layer_endcap->GetYaxis()->GetTitleSize());
    axis_endcap->Draw();
    c_endcap->Write();
    
    // === Resolution plots ===
    TDirectory* resDir = outFile->mkdir("Resolutions");
    resDir->cd();
    
    TCanvas* c_res_qpt = MuCollStyle::CreateCanvas("c_seed_res_qpt", "Seed q/p_{T} Resolution");
    MuCollStyle::StyleHist(h_resolutions_q_over_pt, MuCollStyle::GetColor(0));
    h_resolutions_q_over_pt->GetXaxis()->SetTitle("#Delta(q/p_{T}) / (q/p_{T})");
    h_resolutions_q_over_pt->GetYaxis()->SetTitle("Entries");
    h_resolutions_q_over_pt->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_qpt, "10 TeV");
    c_res_qpt->Write();
    
    TCanvas* c_res_d0 = MuCollStyle::CreateCanvas("c_seed_res_d0", "Seed d_{0} Resolution");
    MuCollStyle::StyleHist(h_resolutions_d0, MuCollStyle::GetColor(1));
    h_resolutions_d0->GetXaxis()->SetTitle("#Delta d_{0} [mm]");
    h_resolutions_d0->GetYaxis()->SetTitle("Entries");
    h_resolutions_d0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_d0, "10 TeV");
    c_res_d0->Write();
    
    TCanvas* c_res_z0 = MuCollStyle::CreateCanvas("c_seed_res_z0", "Seed z_{0} Resolution");
    MuCollStyle::StyleHist(h_resolutions_z0, MuCollStyle::GetColor(2));
    h_resolutions_z0->GetXaxis()->SetTitle("#Delta z_{0} [mm]");
    h_resolutions_z0->GetYaxis()->SetTitle("Entries");
    h_resolutions_z0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_z0, "10 TeV");
    c_res_z0->Write();
    
    outFile->Write();
    outFile->Close();
    inFile->Close();
    
    printf("Plots written to %s\n", outputFile);
}
