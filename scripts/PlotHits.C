#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <cmath>

// Include the plot style header
#include "SimplePlotTemplate.h"

void PlotHits(const char* inputFile, const char* outputFile) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    
    if (!inputFile || !outputFile) {
        printf("Usage: PlotHits(\"hits.root\", \"output_plots.root\")\n");
        return;
    }
    
    // Open input file with ntuples
    TFile* inFile = TFile::Open(inputFile);
    if (!inFile || inFile->IsZombie()) {
        printf("Error: Could not open input file %s\n", inputFile);
        return;
    }
    
    // Get the ntuples
    TNtuple* hits_ntuple = (TNtuple*)inFile->Get("hits");
    TNtuple* events_ntuple = (TNtuple*)inFile->Get("events_summary");
    
    if (!hits_ntuple || !events_ntuple) {
        printf("Error: Could not find required ntuples in file\n");
        inFile->Close();
        return;
    }
    
    printf("Loaded ntuples from %s\n", inputFile);
    printf("  hits: %lld entries\n", hits_ntuple->GetEntries());
    printf("  events_summary: %lld entries\n", events_ntuple->GetEntries());
    
    // Define theta bins
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;
    const double thetaBinWidth = (thetaMax - thetaMin) / nThetaBins;
    
    // Layer names
    const char* layerNames[] = {"VXD", "IT", "OT"};
    const char* layerLabels[] = {"Vertex Detector", "Inner Tracker", "Outer Tracker"};
    
    // Create histograms for hit density per layer
    TH1D* h_density[3];
    for (int layer = 0; layer < 3; layer++) {
        h_density[layer] = new TH1D(Form("h_%s_density", layerNames[layer]), "", nThetaBins, thetaMin, thetaMax);
    }
    
    // Hit count arrays
    int hitCounts[3][nThetaBins] = {0};
    
    // Count total events
    Long64_t totalEvents = events_ntuple->GetEntries();
    
    // Fill hit counts from ntuple
    Float_t h_theta, h_layer, h_isBarrel;
    hits_ntuple->SetBranchAddress("theta", &h_theta);
    hits_ntuple->SetBranchAddress("layer", &h_layer);
    hits_ntuple->SetBranchAddress("isBarrel", &h_isBarrel);
    
    for (Long64_t i = 0; i < hits_ntuple->GetEntries(); i++) {
        hits_ntuple->GetEntry(i);
        
        int layer = (int)h_layer;
        int thetaBin = (int)((h_theta - thetaMin) / thetaBinWidth);
        
        if (layer >= 0 && layer < 3 && thetaBin >= 0 && thetaBin < nThetaBins) {
            hitCounts[layer][thetaBin]++;
        }
    }
    
    // Calculate hit densities
    for (int layer = 0; layer < 3; layer++) {
        for (int bin = 0; bin < nThetaBins; bin++) {
            double thetaLow = thetaMin + bin * thetaBinWidth;
            double thetaHigh = thetaMin + (bin + 1) * thetaBinWidth;
            double angularArea = 2 * 3.14159 * (cos(thetaLow) - cos(thetaHigh));
            
            double density = (double)hitCounts[layer][bin] / (angularArea * totalEvents);
            h_density[layer]->SetBinContent(bin + 1, density);
        }
    }
    
    // === Create plots ===
    
    // Create output file
    TFile* outFile = new TFile(outputFile, "RECREATE");
    
    // Directory for hit density plots
    TDirectory* hitDir = outFile->mkdir("HitDensity");
    hitDir->cd();
    
    // Individual layer plots
    for (int layer = 0; layer < 3; layer++) {
        TCanvas* c = MuCollStyle::CreateCanvas(Form("c_%s", layerNames[layer]), 
                                                Form("%s Hit Density", layerLabels[layer]));
        
        MuCollStyle::StyleHist(h_density[layer], MuCollStyle::GetColor(layer));
        h_density[layer]->GetXaxis()->SetTitle("#theta [rad]");
        h_density[layer]->GetYaxis()->SetTitle("Hit Density [hits/sr/event]");
        h_density[layer]->Draw("PE");
        
        MuCollStyle::AddStandardLabels(c, "10 TeV");
        
        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(h_density[layer], layerLabels[layer], "pe");
        leg->Draw();
        
        c->Write();
    }
    
    // Combined plot
    TCanvas* c_combined = MuCollStyle::CreateCanvas("c_combined", "Hit Density - All Layers");
    
    // Find maximum for y-axis range
    double maxVal = 0;
    for (int layer = 0; layer < 3; layer++) {
        double thisMax = h_density[layer]->GetMaximum();
        if (thisMax > maxVal) maxVal = thisMax;
    }
    
    h_density[0]->SetMaximum(maxVal * 1.2);
    h_density[0]->GetXaxis()->SetTitle("#theta [rad]");
    h_density[0]->GetYaxis()->SetTitle("Hit Density [hits/sr/event]");
    h_density[0]->Draw("PE");
    h_density[1]->Draw("PE SAME");
    h_density[2]->Draw("PE SAME");
    
    TLegend* leg_combined = new TLegend(0.65, 0.65, 0.88, 0.88);
    for (int layer = 0; layer < 3; layer++) {
        leg_combined->AddEntry(h_density[layer], layerLabels[layer], "pe");
    }
    leg_combined->Draw();
    
    MuCollStyle::AddStandardLabels(c_combined, "10 TeV");
    c_combined->Write();
    
    // === Event summary plots ===
    TDirectory* evtDir = outFile->mkdir("EventSummary");
    evtDir->cd();
    
    // Create histograms for event summary
    TH1D* h_nHits_total = new TH1D("h_nHits_total", "", 100, 0, 10000);
    TH1D* h_nHits_VXD = new TH1D("h_nHits_VXD", "", 100, 0, 5000);
    TH1D* h_nHits_IT = new TH1D("h_nHits_IT", "", 100, 0, 3000);
    TH1D* h_nHits_OT = new TH1D("h_nHits_OT", "", 100, 0, 2000);
    
    Float_t e_nHits_VXD, e_nHits_IT, e_nHits_OT, e_nHits_total;
    events_ntuple->SetBranchAddress("nHits_VXD", &e_nHits_VXD);
    events_ntuple->SetBranchAddress("nHits_IT", &e_nHits_IT);
    events_ntuple->SetBranchAddress("nHits_OT", &e_nHits_OT);
    events_ntuple->SetBranchAddress("nHits_total", &e_nHits_total);
    
    for (Long64_t i = 0; i < events_ntuple->GetEntries(); i++) {
        events_ntuple->GetEntry(i);
        h_nHits_total->Fill(e_nHits_total);
        h_nHits_VXD->Fill(e_nHits_VXD);
        h_nHits_IT->Fill(e_nHits_IT);
        h_nHits_OT->Fill(e_nHits_OT);
    }
    
    TCanvas* c_total = MuCollStyle::CreateCanvas("c_total_hits", "Total Hits per Event");
    MuCollStyle::StyleHist(h_nHits_total, MuCollStyle::GetColor(0));
    h_nHits_total->GetXaxis()->SetTitle("Hits per Event");
    h_nHits_total->GetYaxis()->SetTitle("Events");
    h_nHits_total->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_total, "10 TeV");
    c_total->Write();
    
    TCanvas* c_layers = MuCollStyle::CreateCanvas("c_layer_hits", "Hits per Event by Layer");
    MuCollStyle::StyleHist(h_nHits_VXD, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_nHits_IT, MuCollStyle::GetColor(1));
    MuCollStyle::StyleHist(h_nHits_OT, MuCollStyle::GetColor(2));
    
    // Find max
    double maxHits = h_nHits_VXD->GetMaximum();
    if (h_nHits_IT->GetMaximum() > maxHits) maxHits = h_nHits_IT->GetMaximum();
    if (h_nHits_OT->GetMaximum() > maxHits) maxHits = h_nHits_OT->GetMaximum();
    
    h_nHits_VXD->SetMaximum(maxHits * 1.2);
    h_nHits_VXD->GetXaxis()->SetTitle("Hits per Event");
    h_nHits_VXD->GetYaxis()->SetTitle("Events");
    h_nHits_VXD->Draw("HIST");
    h_nHits_IT->Draw("HIST SAME");
    h_nHits_OT->Draw("HIST SAME");
    
    TLegend* leg_layers = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg_layers->AddEntry(h_nHits_VXD, "Vertex Detector", "l");
    leg_layers->AddEntry(h_nHits_IT, "Inner Tracker", "l");
    leg_layers->AddEntry(h_nHits_OT, "Outer Tracker", "l");
    leg_layers->Draw();
    MuCollStyle::AddStandardLabels(c_layers, "10 TeV");
    c_layers->Write();
    
    outFile->Write();
    outFile->Close();
    inFile->Close();
    
    printf("Plots written to %s\n", outputFile);
}
