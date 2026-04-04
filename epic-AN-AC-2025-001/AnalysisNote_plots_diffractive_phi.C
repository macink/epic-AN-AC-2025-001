#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void plot_e_theta()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_theta_MC = (TH1D*) file->Get("h_theta_MC"); 
    TH1D* h_theta_REC_EEMC = (TH1D*) file->Get("h_theta_REC_EEMC");
    TH2D* h_theta_response_EEMC_beforeCuts = (TH2D*) file->Get("h_theta_response_EEMC_beforeCuts"); 
    TH2D* h_theta_resolution_beforeCuts = (TH2D*) file->Get("h_theta_resolution");

    TCanvas* c6 = new TCanvas("c6","c6",1,1,1600,600);
    c6->Divide(3,1,0.01,0.01);
    c6->cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    gStyle->SetOptStat(0);
    h_theta_MC->GetXaxis()->SetTitleOffset(1.2);
    h_theta_MC->GetYaxis()->SetTitleOffset(1.8);
    h_theta_MC->GetXaxis()->SetLabelSize(0.04);  
    h_theta_MC->GetYaxis()->SetLabelSize(0.04);
    h_theta_MC->GetXaxis()->SetTitleSize(0.04);  
    h_theta_MC->GetYaxis()->SetTitleSize(0.04);
	h_theta_MC->GetYaxis()->SetTitle("counts");
	h_theta_MC->GetXaxis()->SetTitle("#theta");
    h_theta_MC->SetLineColor(kBlack);
    h_theta_MC->Draw();
    h_theta_REC_EEMC->SetMarkerColor(kBlue);
    h_theta_REC_EEMC->SetMarkerStyle(24);
    h_theta_REC_EEMC->Draw("PEsame");
    TLegend *w6 = new TLegend(0.25,0.7,0.5,0.85);
    w6->AddEntry(h_theta_MC, "MC", "L");
    w6->AddEntry(h_theta_REC_EEMC, "RECO", "P");
    w6->SetBorderSize(0);   
    w6->SetFillStyle(0);
    w6->SetTextSize(0.05);
	w6->Draw("same");

    c6->cd(2);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_theta_resolution_beforeCuts->GetXaxis()->SetNdivisions(505);
    h_theta_resolution_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_theta_resolution_beforeCuts->GetYaxis()->SetTitleOffset(2.1);
    h_theta_resolution_beforeCuts->GetXaxis()->SetRangeUser(2.7,3.1);
    h_theta_resolution_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_theta_resolution_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_theta_resolution_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_theta_resolution_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_theta_resolution_beforeCuts->GetYaxis()->SetRangeUser(-0.02,0.02);
    h_theta_resolution_beforeCuts->Draw("colzsame");

    c6->cd(3);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetNdivisions(505);
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_theta_response_EEMC_beforeCuts->GetYaxis()->SetTitleOffset(2.1);
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetRangeUser(2.7,3.1);
    h_theta_response_EEMC_beforeCuts->GetYaxis()->SetRangeUser(2.7,3.1);
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_theta_response_EEMC_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_theta_response_EEMC_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_theta_response_EEMC_beforeCuts->GetYaxis()->SetTitle("#theta_{MC}");
    h_theta_response_EEMC_beforeCuts->GetXaxis()->SetTitle("#theta_{RECO}");
    h_theta_response_EEMC_beforeCuts->Draw("colzsame");

    c6->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title6 = new TLatex();
    title6->SetNDC(); 
    title6->SetTextSize(0.05);
    title6->SetTextAlign(22);  
    title6->DrawLatex(0.5, 0.97, "#theta Truth vs Reco");  

    c6->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_4 = new TLatex();
    title6_4->SetNDC(); 
    title6_4->SetTextSize(0.05);
    title6_4->SetTextAlign(22);  
    title6_4->DrawLatex(0.5, 0.97, "#theta Resolution");  

    c6->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title6_3 = new TLatex();
    title6_3->SetNDC(); 
    title6_3->SetTextSize(0.05);
    title6_3->SetTextAlign(22);  
    title6_3->DrawLatex(0.5, 0.97, "#theta Response");  
	
    c6->Print("./figures/analysisNote_plots/plot_theta.pdf");
}

void e_theta_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_theta_MC_after = (TH1D*) file->Get("h_theta_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_theta_migration_beforeCuts = (TH2D*) file2->Get("h_theta_migration_beforeCuts"); 
    TH1D* h_theta_MC = (TH1D*) file2->Get("h_theta_MC");
    TH1D* h_theta_REC_EEMC = (TH1D*) file2->Get("h_theta_REC_EEMC");

    TCanvas* c19 = new TCanvas("c19","c19",1,1,1400,800);
    c19->Divide(3,2,0.01,0.01);
    c19->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetLogz(1);
    h_theta_migration_beforeCuts->GetXaxis()->SetTitleOffset(1.1);
    h_theta_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_theta_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_theta_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_theta_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_theta_migration_beforeCuts->GetXaxis()->SetRangeUser(2.7,3.07);
    h_theta_migration_beforeCuts->GetYaxis()->SetRangeUser(2.7,3.07);
    h_theta_migration_beforeCuts->Draw("colzsame");

    int nx_theta = h_theta_migration_beforeCuts->GetNbinsX();
    int ny_theta = h_theta_migration_beforeCuts->GetNbinsY();
    TH1D* h_theta_purity = new TH1D("h_theta_purity",";#theta_{reco bin};Purity",nx_theta,0,3.14);
    for (int j=1; j<=ny_theta; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_theta; ++i) recoColSum += h_theta_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_theta_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_theta_purity->SetBinContent(j, diag / recoColSum);
    }
    c19->cd(2);
    gStyle->SetOptStat(0);
    h_theta_purity->GetXaxis()->SetRangeUser(2.2,3.14);
    h_theta_purity->GetXaxis()->SetLabelSize(0.04);  
    h_theta_purity->GetYaxis()->SetLabelSize(0.04);
    h_theta_purity->GetXaxis()->SetTitleSize(0.04);  
    h_theta_purity->GetYaxis()->SetTitleSize(0.04);
    h_theta_purity->GetXaxis()->SetTitle("#theta_{reco bin}");
    h_theta_purity->GetYaxis()->SetTitle("Purity");
    h_theta_purity->Draw("same");

    double binWidth_thetaPurity = h_theta_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label19_2 = new TLatex(0.32,0.71, Form("bin width=%.3f",binWidth_thetaPurity));
    label19_2->SetNDC();
    label19_2->SetTextSize(0.04);
    label19_2->Draw("same");

    TH1D* h_theta_stability = new TH1D("h_theta_stability",";#theta_{true bin};stability",nx_theta,0,3.14);
    for (int i=1; i<=nx_theta; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_theta; ++j) trueRowSum += h_theta_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_theta_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_theta_stability->SetBinContent(i, diag / trueRowSum);
    }
    c19->cd(3);
    gStyle->SetOptStat(0);
    h_theta_stability->GetXaxis()->SetLabelSize(0.04);  
    h_theta_stability->GetYaxis()->SetLabelSize(0.04);
    h_theta_stability->GetXaxis()->SetTitleSize(0.04);  
    h_theta_stability->GetYaxis()->SetTitleSize(0.04);
    h_theta_stability->GetXaxis()->SetRangeUser(2.2,3.14);
    h_theta_stability->GetYaxis()->SetTitle("Stability");
    h_theta_stability->Draw("same");

    double binWidth_thetaStability = h_theta_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label19_3 = new TLatex(0.32,0.71, Form("bin width=%.3f",binWidth_thetaStability));
    label19_3->SetNDC();
    label19_3->SetTextSize(0.04);
    label19_3->Draw("same");
    
    int nb_theta = h_theta_MC->GetNbinsX();
    TH1D* h_theta_efficiency = new TH1D("h_theta_efficiency", ";#theta;Efficiency", nb_theta, 0, 3.14);
    TH1D* h_theta_acceptance = new TH1D("h_theta_acceptance", ";#theta;Acceptance", nb_theta, 0, 3.14);
    TH1D* h_theta_corrected = new TH1D("h_theta_corrected", ";#theta;Acceptance Corrected", nb_theta, 0, 3.14);
    for (int i=1; i<=nb_theta; ++i) 
    {
        double Ntrue = h_theta_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_theta; ++j) Naccepted += h_theta_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_theta_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_theta_REC_EEMC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_theta_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_theta_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_theta_efficiency->SetBinContent(i, efficiency);
    }
    c19->cd(4);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetLogz(1);
    h_theta_acceptance->GetXaxis()->SetRangeUser(2.2,3.14);
    h_theta_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_theta_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_theta_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_theta_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_theta_acceptance->GetXaxis()->SetTitle("#theta");
    h_theta_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_theta_acceptance->Draw("same");

    c19->cd(5);
    gStyle->SetOptStat(0);
    h_theta_efficiency->GetXaxis()->SetRangeUser(2.2,3.14);
    h_theta_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_theta_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_theta_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_theta_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_theta_efficiency->GetXaxis()->SetTitle("#theta");
    h_theta_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_theta_efficiency->Draw("same");

    c19->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_theta_corrected->GetXaxis()->SetRangeUser(2.2,3.14);
    h_theta_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_theta_corrected->GetYaxis()->SetLabelSize(0.04);
    h_theta_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_theta_corrected->GetYaxis()->SetTitleSize(0.04);
    h_theta_corrected->GetYaxis()->SetTitle("Corrected");
    h_theta_corrected->Draw("same");

    c19->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_1 = new TLatex();
    title19_1->SetNDC(); 
    title19_1->SetTextSize(0.05);
    title19_1->SetTextAlign(22);  
    title19_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c19->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_3 = new TLatex();
    title19_3->SetNDC(); 
    title19_3->SetTextSize(0.05);
    title19_3->SetTextAlign(22);  
    title19_3->DrawLatex(0.5, 0.97, "Purity");  

    c19->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_4 = new TLatex();
    title19_4->SetNDC(); 
    title19_4->SetTextSize(0.05);
    title19_4->SetTextAlign(22);  
    title19_4->DrawLatex(0.5, 0.97, "Stability");

    c19->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_5 = new TLatex();
    title19_5->SetNDC(); 
    title19_5->SetTextSize(0.05);
    title19_5->SetTextAlign(22);  
    title19_5->DrawLatex(0.5, 0.97, "Acceptance");

    c19->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_2 = new TLatex();
    title19_2->SetNDC(); 
    title19_2->SetTextSize(0.05);
    title19_2->SetTextAlign(22);  
    title19_2->DrawLatex(0.5, 0.97, "Efficiency");  

    c19->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title19_6 = new TLatex();
    title19_6->SetNDC(); 
    title19_6->SetTextSize(0.05);
    title19_6->SetTextAlign(22);  
    title19_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");

    c19->Print("./figures/analysisNote_plots/plot_theta_migrationCorrectedPurityAcceptance.pdf");
}

void plot_e_energy()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_energy_MC = (TH1D*) file->Get("h_energy_MC"); 
    TH1D* h_e_energy_beforeCuts = (TH1D*) file->Get("h_e_energy_beforeCuts");
    TH2D* h_energy_res_EEMC_beforeCuts = (TH2D*) file->Get("h_energy_res_EEMC_beforeCuts");
    TH2D* h_energy_response_EEMC_beforeCuts = (TH2D*) file->Get("h_energy_response_EEMC_beforeCuts");

    TCanvas* c3 = new TCanvas("c3","c3",1,1,1600,600);
    c3->Divide(3,1,0.01,0.01);
    c3->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    h_energy_MC->GetXaxis()->SetTitleOffset(1.2);
    h_energy_MC->GetYaxis()->SetTitleOffset(1.5);
    h_energy_MC->GetXaxis()->SetLabelSize(0.04);  
    h_energy_MC->GetYaxis()->SetLabelSize(0.04);
    h_energy_MC->GetXaxis()->SetTitleSize(0.04);  
    h_energy_MC->GetYaxis()->SetTitleSize(0.04);
    h_energy_MC->GetYaxis()->SetTitle("counts");
    h_energy_MC->GetXaxis()->SetTitle("E [GeV]");
    h_energy_MC->SetLineColor(kBlack);
    h_energy_MC->Draw();
    h_e_energy_beforeCuts->SetMarkerStyle(24);
    h_e_energy_beforeCuts->SetMarkerColor(kBlue);
    h_e_energy_beforeCuts->Draw("PEsame");
    TLegend *w3 = new TLegend(0.65,0.7,0.8,0.85);
    w3->AddEntry(h_energy_MC, " MC", "L");
    w3->AddEntry(h_e_energy_beforeCuts, " RECO", "P");
    w3->SetBorderSize(0);   
    w3->SetFillStyle(0);
    w3->SetTextSize(0.05);
	w3->Draw("same");

    c3->cd(2);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_energy_res_EEMC_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_energy_res_EEMC_beforeCuts->GetYaxis()->SetTitleOffset(1.8);
    h_energy_res_EEMC_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_energy_res_EEMC_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_energy_res_EEMC_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_energy_res_EEMC_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_energy_res_EEMC_beforeCuts->GetXaxis()->SetRangeUser(6,11);
    h_energy_res_EEMC_beforeCuts->GetYaxis()->SetRangeUser(-0.2,0.2);
    h_energy_res_EEMC_beforeCuts->GetYaxis()->SetTitle("(E_{MC}-E_{RECO})/E_{MC}");
    h_energy_res_EEMC_beforeCuts->GetXaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_res_EEMC_beforeCuts->Draw("colzsame");

    c3->cd(3);
    gPad->SetLeftMargin(0.18);
    gPad->SetLogz(1);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_energy_response_EEMC_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_energy_response_EEMC_beforeCuts->GetYaxis()->SetTitleOffset(1.5);
    h_energy_response_EEMC_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_energy_response_EEMC_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_energy_response_EEMC_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_energy_response_EEMC_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_energy_response_EEMC_beforeCuts->GetXaxis()->SetRangeUser(6,11);
    h_energy_response_EEMC_beforeCuts->GetYaxis()->SetRangeUser(6,11);
    h_energy_response_EEMC_beforeCuts->GetYaxis()->SetTitle("E_{MC} [GeV]");
    h_energy_response_EEMC_beforeCuts->GetXaxis()->SetTitle("E_{RECO} [GeV]");
    h_energy_response_EEMC_beforeCuts->Draw("colzsame");

    c3->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title3 = new TLatex();
    title3->SetNDC(); 
    title3->SetTextSize(0.05);
    title3->SetTextAlign(22);  
    title3->DrawLatex(0.5, 0.97, "e' Energy Truth vs. Reco");  

    c3->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_2 = new TLatex();
    title3_2->SetNDC(); 
    title3_2->SetTextSize(0.05);
    title3_2->SetTextAlign(22);  
    title3_2->DrawLatex(0.5, 0.97, "e' Energy Resolution");

    c3->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title3_3 = new TLatex();
    title3_3->SetNDC(); 
    title3_3->SetTextSize(0.05);
    title3_3->SetTextAlign(22);  
    title3_3->DrawLatex(0.5, 0.97, "e' Energy Response");  

    c3->Print("./figures/analysisNote_plots/plot_energy.pdf");
}

void e_energy_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_energy_MC_after = (TH1D*) file->Get("h_energy_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_energy_migration_beforeCuts = (TH2D*) file2->Get("h_energy_migration_beforeCuts"); 
    TH1D* h_energy_MC = (TH1D*) file2->Get("h_energy_MC");
    TH1D* h_e_energy_beforeCuts = (TH1D*) file2->Get("h_e_energy_beforeCuts");

    TCanvas* c18 = new TCanvas("c18","c18",1,1,1400,800);
    c18->Divide(3,2,0.01,0.01);
    c18->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.18);
    gPad->SetLogz(1);
    h_energy_migration_beforeCuts->GetXaxis()->SetTitleOffset(1.5);
    h_energy_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_energy_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_energy_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_energy_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_energy_migration_beforeCuts->GetXaxis()->SetRangeUser(6,11);
    h_energy_migration_beforeCuts->GetYaxis()->SetRangeUser(6,11);
    h_energy_migration_beforeCuts->Draw("colzsame");

    int nx_energy = h_energy_migration_beforeCuts->GetNbinsX();
    int ny_energy = h_energy_migration_beforeCuts->GetNbinsY();
    TH1D* h_energy_purity = new TH1D("h_energy_purity",";E_{reco bin} [GeV];Purity",nx_energy,0,20);
    // purity for reco bin j: fraction of events in reco bin j that came from same true bin j
    for (int j=1; j<=ny_energy; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_energy; ++i) recoColSum += h_energy_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_energy_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_energy_purity->SetBinContent(j, diag / recoColSum);
    }
    c18->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    h_energy_purity->GetXaxis()->SetTitleOffset(1.5);
    h_energy_purity->GetXaxis()->SetRangeUser(0,12);
    h_energy_purity->GetXaxis()->SetLabelSize(0.04);  
    h_energy_purity->GetYaxis()->SetLabelSize(0.04);
    h_energy_purity->GetXaxis()->SetTitleSize(0.04);  
    h_energy_purity->GetYaxis()->SetTitleSize(0.04);
    h_energy_purity->GetXaxis()->SetTitle("E_{reco bin} [GeV]");
    h_energy_purity->GetYaxis()->SetTitle("Purity");
    h_energy_purity->Draw("same");

    double binWidth_energyPurtiy = h_energy_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label18_2 = new TLatex(0.22,0.81, Form("bin width=%.1f GeV",binWidth_energyPurtiy));
    label18_2->SetNDC();
    label18_2->SetTextSize(0.04);
    label18_2->Draw("same");

    TH1D* h_energy_stability = new TH1D("h_energy_stability",";E_{true bin} [GeV];stability",nx_energy,0,20);
    for (int i=1; i<=nx_energy; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_energy; ++j) trueRowSum += h_energy_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_energy_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_energy_stability->SetBinContent(i, diag / trueRowSum);
    }
    c18->cd(3);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    h_energy_stability->GetXaxis()->SetTitleOffset(1.5);
    h_energy_stability->GetXaxis()->SetRangeUser(0,12);
    h_energy_stability->GetXaxis()->SetLabelSize(0.04);  
    h_energy_stability->GetYaxis()->SetLabelSize(0.04);
    h_energy_stability->GetXaxis()->SetTitleSize(0.04);  
    h_energy_stability->GetYaxis()->SetTitleSize(0.04);
    h_energy_stability->GetYaxis()->SetTitle("Stability");
    h_energy_stability->Draw("same");

    double binWidth_energyStability = h_energy_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label18_3 = new TLatex(0.25,0.81, Form("bin width=%.1f GeV",binWidth_energyStability));
    label18_3->SetNDC();
    label18_3->SetTextSize(0.04);
    label18_3->Draw("same");

    int nb_energy = h_energy_MC->GetNbinsX();
    TH1D* h_energy_efficiency = new TH1D("h_energy_efficiency", ";E [GeV];Efficiency", nb_energy, 0, 20);
    TH1D* h_energy_acceptance = new TH1D("h_energy_acceptance", ";E [GeV];Acceptance", nb_energy, 0, 20);
    TH1D* h_energy_corrected = new TH1D("h_energy_corrected", ";E [GeV];Acceptance Corrected", nb_energy, 0, 20);
    for (int i=1; i<=nb_energy; ++i) 
    {
        double Ntrue = h_energy_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_energy; ++j) Naccepted += h_energy_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_energy_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_e_energy_beforeCuts->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_energy_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_energy_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_energy_efficiency->SetBinContent(i, efficiency);
    }
    c18->cd(4);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    h_energy_acceptance->GetXaxis()->SetTitleOffset(1.);
    h_energy_acceptance->GetXaxis()->SetRangeUser(0,12);
    h_energy_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_energy_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_energy_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_energy_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_energy_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_energy_acceptance->Draw("same");

    c18->cd(5);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    h_energy_efficiency->GetXaxis()->SetTitleOffset(1.5);
    h_energy_efficiency->GetXaxis()->SetRangeUser(0,12);
    h_energy_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_energy_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_energy_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_energy_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_energy_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_energy_efficiency->GetXaxis()->SetTitle("E [GeV]");
    h_energy_efficiency->Draw("same");

    c18->cd(6);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);
    h_energy_corrected->GetXaxis()->SetTitleOffset(1.5);
    h_energy_corrected->GetXaxis()->SetRangeUser(0,12);
    h_energy_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_energy_corrected->GetYaxis()->SetLabelSize(0.04);
    h_energy_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_energy_corrected->GetYaxis()->SetTitleSize(0.04);
    h_energy_corrected->GetYaxis()->SetTitle("Corrected");
    h_energy_corrected->Draw("same");

    c18->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_1 = new TLatex();
    title18_1->SetNDC(); 
    title18_1->SetTextSize(0.05);
    title18_1->SetTextAlign(22);  
    title18_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c18->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_3 = new TLatex();
    title18_3->SetNDC(); 
    title18_3->SetTextSize(0.05);
    title18_3->SetTextAlign(22);  
    title18_3->DrawLatex(0.5, 0.97, "Purity");  

    c18->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_4 = new TLatex();
    title18_4->SetNDC(); 
    title18_4->SetTextSize(0.05);
    title18_4->SetTextAlign(22);  
    title18_4->DrawLatex(0.5, 0.97, "Stability");

    c18->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_5 = new TLatex();
    title18_5->SetNDC(); 
    title18_5->SetTextSize(0.05);
    title18_5->SetTextAlign(22);  
    title18_5->DrawLatex(0.5, 0.97, "Acceptance");

    c18->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_2 = new TLatex();
    title18_2->SetNDC(); 
    title18_2->SetTextSize(0.05);
    title18_2->SetTextAlign(22);  
    title18_2->DrawLatex(0.5, 0.97, "Efficiency");  

    c18->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title18_6 = new TLatex();
    title18_6->SetNDC(); 
    title18_6->SetTextSize(0.05);
    title18_6->SetTextAlign(22);  
    title18_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");
        
    c18->Print("./figures/analysisNote_plots/plot_energy_migrationCorrectedPurityAcceptance.pdf");
}

void plot_Q2()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_Q2_e = (TH1D*) file->Get("h_Q2_e"); 
    TH1D* h_Q2REC_e_EEMC = (TH1D*) file->Get("h_Q2_afterCut"); 
    TH2D* h_Q2_res_beforecut = (TH2D*) file->Get("h_Q2_res_beforecut"); 
    TH2D* h_Q2_response_beforecut = (TH2D*) file->Get("h_Q2_response_beforecut"); 

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1800,600);
    c1->Divide(3,1,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0); 
    h_Q2_e->SetStats(0);
    h_Q2_e->GetXaxis()->SetTitleOffset(1.2);
    h_Q2_e->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_e->GetXaxis()->SetRangeUser(1,10);
    h_Q2_e->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_e->GetYaxis()->SetLabelSize(0.04);
    h_Q2_e->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_e->GetYaxis()->SetTitleSize(0.04);
	h_Q2_e->GetYaxis()->SetTitle("counts");
	h_Q2_e->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
	h_Q2_e->SetLineColor(kBlack);
    h_Q2_e->Draw();
    h_Q2REC_e_EEMC->SetMarkerStyle(24);
    h_Q2REC_e_EEMC->SetMarkerColor(kBlue);
    h_Q2REC_e_EEMC->Draw("PEsame");
    TLegend *w1 = new TLegend(0.65,0.7,0.85,0.85);
	w1->AddEntry(h_Q2_e, " MC", "L");
    w1->AddEntry(h_Q2REC_e_EEMC, " RECO", "P");
    w1->SetBorderSize(0);   
    w1->SetFillStyle(0);
    w1->SetTextSize(0.05);
	w1->Draw("same");

    c1->cd(2);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_res_beforecut->GetXaxis()->SetTitleOffset(1.2);
    h_Q2_res_beforecut->GetYaxis()->SetTitleOffset(1.8);
    h_Q2_res_beforecut->GetXaxis()->SetRangeUser(1,10);
    h_Q2_res_beforecut->GetYaxis()->SetRangeUser(-0.3,0.3);
    h_Q2_res_beforecut->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_res_beforecut->GetYaxis()->SetLabelSize(0.04);
    h_Q2_res_beforecut->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_res_beforecut->GetYaxis()->SetTitleSize(0.04);
    h_Q2_res_beforecut->GetYaxis()->SetTitle("(Q^{2}_{MC}-Q^{2}_{RECO})/Q^{2}_{MC}");
    h_Q2_res_beforecut->GetXaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_res_beforecut->Draw("colzsame");

    c1->cd(3);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_response_beforecut->GetXaxis()->SetTitleOffset(1.2);
    h_Q2_response_beforecut->GetXaxis()->SetRangeUser(1,10);
    h_Q2_response_beforecut->GetYaxis()->SetRangeUser(1,10);
    h_Q2_response_beforecut->GetYaxis()->SetTitleOffset(1.5);
    h_Q2_response_beforecut->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_response_beforecut->GetYaxis()->SetLabelSize(0.04);
    h_Q2_response_beforecut->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_response_beforecut->GetYaxis()->SetTitleSize(0.04);
    h_Q2_response_beforecut->GetYaxis()->SetTitle("Q^{2}_{MC} [GeV/c]^{2}");
    h_Q2_response_beforecut->GetXaxis()->SetTitle("Q^{2}_{RECO} [GeV/c]^{2}");
    h_Q2_response_beforecut->Draw("colzsame");

    c1->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1 = new TLatex();
    title1_1->SetNDC(); 
    title1_1->SetTextSize(0.05);
    title1_1->SetTextAlign(22);  
    title1_1->DrawLatex(0.5, 0.97, "Q^{2} Truth vs Reco");  

    c1->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_2 = new TLatex();
    title1_2->SetNDC(); 
    title1_2->SetTextSize(0.05);
    title1_2->SetTextAlign(22);  
    title1_2->DrawLatex(0.5, 0.97, "Q^{2} Resolution");

    c1->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_3 = new TLatex();
    title1_3->SetNDC(); 
    title1_3->SetTextSize(0.05);
    title1_3->SetTextAlign(22);  
    title1_3->DrawLatex(0.5, 0.97, "Q^{2} Response");

    c1->Print("./figures/analysisNote_plots/plot_Q2.pdf");
}

void Q2_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_Q2_MC_after = (TH1D*) file->Get("h_Q2_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_Q2_migration_beforeCuts = (TH2D*) file2->Get("h_Q2_migration");
    TH1D* h_Q2REC_e_EEMC = (TH1D*) file2->Get("h_Q2REC_e_EEMC");
    TH1D* h_Q2_e = (TH1D*) file2->Get("h_Q2_e");

    TCanvas* c23 = new TCanvas("c23","c23",1,1,1400,800);
    c23->Divide(3,2,0.01,0.01);
    c23->cd(1);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_Q2_migration_beforeCuts->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_Q2_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_Q2_migration_beforeCuts->GetYaxis()->SetRangeUser(1,10);
    h_Q2_migration_beforeCuts->GetXaxis()->SetRangeUser(1,10);
    h_Q2_migration_beforeCuts->Draw("colzsame");

    int nx_Q2 = h_Q2_migration_beforeCuts->GetNbinsX();
    int ny_Q2 = h_Q2_migration_beforeCuts->GetNbinsY();
    TH1D* h_Q2_purity = new TH1D("h_Q2_purity",";Q^{2}_{reco bin} [GeV/c]^{2};Purity",nx_Q2,0,10);
    for (int j=1; j<=ny_Q2; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_Q2; ++i) recoColSum += h_Q2_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_Q2_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_Q2_purity->SetBinContent(j, diag / recoColSum);
    }
    c23->cd(2);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_purity->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_purity->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_purity->GetYaxis()->SetLabelSize(0.04);
    h_Q2_purity->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_purity->GetYaxis()->SetTitleSize(0.04);
    h_Q2_purity->GetXaxis()->SetTitle("Q^{2}_{reco bin} [GeV/c]^{2}");
    h_Q2_purity->GetYaxis()->SetTitle("Purity");
    h_Q2_purity->Draw("same");

    double binWidth_Q2Purity = h_Q2_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label23_2 = new TLatex(0.32,0.81, Form("bin width=%.1f GeV^{2}/c^{2}",binWidth_Q2Purity));
    label23_2->SetNDC();
    label23_2->SetTextSize(0.04);
    label23_2->Draw("same");

    TH1D* h_Q2_stability = new TH1D("h_Q2_stability",";Q^{2}_{true bin};Stability",nx_Q2,0,10);
    for (int i=1; i<=nx_Q2; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_Q2; ++j) trueRowSum += h_Q2_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_Q2_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_Q2_stability->SetBinContent(i, diag / trueRowSum);
    }
    c23->cd(3);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_stability->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_stability->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_stability->GetYaxis()->SetLabelSize(0.04);
    h_Q2_stability->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_stability->GetYaxis()->SetTitleSize(0.04);
    h_Q2_stability->GetXaxis()->SetTitle("Q^{2}_{true bin} [GeV/c]^{2}");
    h_Q2_stability->GetYaxis()->SetTitle("Stability");
    h_Q2_stability->Draw("same");

    double binWidth_Q2Stability = h_Q2_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label23_3 = new TLatex(0.32,0.81, Form("bin width=%.1f GeV^{2}/c^{2}",binWidth_Q2Stability));
    label23_3->SetNDC();
    label23_3->SetTextSize(0.04);
    label23_3->Draw("same");
    
    int nb_Q2 = h_Q2_e->GetNbinsX();
    TH1D* h_Q2_efficiency = new TH1D("h_Q2_efficiency", ";Q^{2} [GeV/c]^{2};Efficiency", nb_Q2, 0, 10);
    TH1D* h_Q2_acceptance = new TH1D("h_Q2_acceptance", ";Q^{2} [GeV/c]^{2};Acceptance", nb_Q2, 0, 10);
    TH1D* h_Q2_corrected = new TH1D("h_Q2_corrected", ";Q^{2} [GeV/c]^{2};Acceptance Corrected", nb_Q2, 0, 10);
    for (int i=1; i<=nb_Q2; ++i) 
    {
        double Ntrue = h_Q2_e->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_Q2; ++j) Naccepted += h_Q2_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_Q2_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_Q2REC_e_EEMC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_Q2_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_Q2_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_Q2_efficiency->SetBinContent(i, efficiency);
    }
    c23->cd(4);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_acceptance->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_Q2_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_Q2_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_Q2_acceptance->Draw("same");

    c23->cd(5);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_efficiency->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_Q2_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_Q2_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_Q2_efficiency->GetXaxis()->SetTitle("Q^{2} [GeV/c]^{2}");
    h_Q2_efficiency->Draw("same");
    
    c23->cd(6);
    gPad->SetLogy(1);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_Q2_corrected->GetXaxis()->SetTitleOffset(1.5);
    h_Q2_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_Q2_corrected->GetYaxis()->SetLabelSize(0.04);
    h_Q2_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_Q2_corrected->GetYaxis()->SetTitleSize(0.04);
    h_Q2_corrected->GetYaxis()->SetTitle("Corrected");
    h_Q2_corrected->Draw("same");

    c23->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_1 = new TLatex();
    title23_1->SetNDC(); 
    title23_1->SetTextSize(0.05);
    title23_1->SetTextAlign(22);  
    title23_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c23->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_3 = new TLatex();
    title23_3->SetNDC(); 
    title23_3->SetTextSize(0.05);
    title23_3->SetTextAlign(22);  
    title23_3->DrawLatex(0.5, 0.97, "Purity");  

    c23->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_4 = new TLatex();
    title23_4->SetNDC(); 
    title23_4->SetTextSize(0.05);
    title23_4->SetTextAlign(22);  
    title23_4->DrawLatex(0.5, 0.97, "Stability");

    c23->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_5 = new TLatex();
    title23_5->SetNDC(); 
    title23_5->SetTextSize(0.05);
    title23_5->SetTextAlign(22);  
    title23_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c23->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_2 = new TLatex();
    title23_2->SetNDC(); 
    title23_2->SetTextSize(0.05);
    title23_2->SetTextAlign(22);  
    title23_2->DrawLatex(0.5, 0.97, "Efficiency");  

    c23->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title23_6 = new TLatex();
    title23_6->SetNDC(); 
    title23_6->SetTextSize(0.05);
    title23_6->SetTextAlign(22);  
    title23_6->DrawLatex(0.5, 0.97, "Acceptance Corrected"); 
    
    c23->Print("./figures/analysisNote_plots/plot_Q2_migrationCorrectedPurityAcceptance.pdf");
}

void plot_y()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_y_e = (TH1D*) file->Get("h_y_e");
    TH1D* h_yREC_e_EEMC = (TH1D*) file->Get("h_yREC_e_EEMC"); 
    TH2D* h_y_res_beforeCuts = (TH2D*) file->Get("h_y_res_beforeCuts"); 
    TH2D* h_y_response_beforeCuts = (TH2D*) file->Get("h_y_response_beforeCuts");

    TCanvas* c2 = new TCanvas("c2","c2",1,1,1800,600);
    c2->Divide(3,1,0.01,0.01);
    c2->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
	h_y_e->GetYaxis()->SetTitle("counts");
	h_y_e->GetXaxis()->SetTitle("y");
	h_y_e->SetLineColor(kBlack);
    h_y_e->GetXaxis()->SetLabelSize(0.04);  
    h_y_e->GetYaxis()->SetLabelSize(0.04);
    h_y_e->GetXaxis()->SetTitleSize(0.04);  
    h_y_e->GetYaxis()->SetTitleSize(0.04);
    h_y_e->Draw();
    h_yREC_e_EEMC->SetMarkerStyle(24);
    h_yREC_e_EEMC->SetMarkerColor(kBlue);
    h_yREC_e_EEMC->Draw("PEsame");
    TLegend *w2 = new TLegend(0.65,0.7,0.85,0.85);
	w2->AddEntry(h_y_e, " MC", "L");
    w2->AddEntry(h_yREC_e_EEMC, " RECO", "P");
    w2->SetBorderSize(0);   
    w2->SetFillStyle(0);
    w2->SetTextSize(0.05);
	w2->Draw("same");

    c2->cd(2);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gStyle->SetPalette(kBird); 
    gStyle->SetPadRightMargin(0.15); 
    gStyle->SetLabelSize(0.03,"Z");  
    gStyle->SetTitleSize(0.03,"Z");  
    h_y_res_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_y_res_beforeCuts->GetYaxis()->SetTitleOffset(1.5);
    h_y_res_beforeCuts->GetXaxis()->SetRangeUser(0,0.8);
    h_y_res_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_y_res_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_y_res_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_y_res_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_y_res_beforeCuts->GetXaxis()->SetRangeUser(0,0.4);
	h_y_res_beforeCuts->GetYaxis()->SetTitle("(y_{MC}-y_{RECO})/y_{MC}");
	h_y_res_beforeCuts->GetXaxis()->SetTitle("y_{MC}");
    h_y_res_beforeCuts->Draw("colzsame");
    gPad->Update(); 

    c2->cd(3);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.15);
    gStyle->SetPalette(kBird); 
    gStyle->SetPadRightMargin(0.15); 
    gStyle->SetLabelSize(0.03,"Z"); 
    gStyle->SetTitleSize(0.03,"Z"); 
    h_y_response_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_y_response_beforeCuts->GetYaxis()->SetTitleOffset(1.5);
    h_y_response_beforeCuts->GetXaxis()->SetRangeUser(0,0.4);
    h_y_response_beforeCuts->GetYaxis()->SetRangeUser(0,0.4);
    h_y_response_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_y_response_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_y_response_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_y_response_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_y_response_beforeCuts->GetYaxis()->SetTitle("y_{MC}");
    h_y_response_beforeCuts->GetXaxis()->SetTitle("y_{RECO}");
    h_y_response_beforeCuts->Draw("colzsame");
    gPad->Update(); 

    c2->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title2 = new TLatex();
    title2->SetNDC(); 
    title2->SetTextSize(0.05);
    title2->SetTextAlign(22);  
    title2->DrawLatex(0.5, 0.97, "y Truth vs Reco");  

    c2->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_2 = new TLatex();
    title2_2->SetNDC(); 
    title2_2->SetTextSize(0.05);
    title2_2->SetTextAlign(22);  
    title2_2->DrawLatex(0.5, 0.97, "y Resolution");

    c2->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title2_3 = new TLatex();
    title2_3->SetNDC(); 
    title2_3->SetTextSize(0.05);
    title2_3->SetTextAlign(22);  
    title2_3->DrawLatex(0.5, 0.97, "y Response");
	
    c2->Print("./figures/analysisNote_plots/plot_y.pdf");
}

void y_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_y_MC_after = (TH1D*) file->Get("h_y_MC_after");
    
    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_y_migration_beforeCuts = (TH2D*) file2->Get("h_y_migration_beforeCuts"); 
    TH1D* h_yREC_e_EEMC = (TH1D*) file2->Get("h_yREC_e_EEMC");
    TH1D* h_y_e = (TH1D*) file2->Get("h_y_e");

    TCanvas* c24 = new TCanvas("c24","c24",1,1,1400,800);
    c24->Divide(3,2,0.01,0.01);
    c24->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogz(1);
    h_y_migration_beforeCuts->GetXaxis()->SetTitleOffset(1.5);
    h_y_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_y_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_y_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_y_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_y_migration_beforeCuts->GetXaxis()->SetRangeUser(0,0.4);
    h_y_migration_beforeCuts->GetYaxis()->SetRangeUser(0,0.4);
    h_y_migration_beforeCuts->Draw("colzsame");

    int nx_y = h_y_migration_beforeCuts->GetNbinsX();
    int ny_y = h_y_migration_beforeCuts->GetNbinsY();
    TH1D* h_y_purity = new TH1D("h_y_purity", ";y_{reco bin};Purity", nx_y, 0.01, 0.85);
    for (int j=1; j<=ny_y; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_y; ++i) recoColSum += h_y_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_y_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_y_purity->SetBinContent(j, diag / recoColSum);
    }
    c24->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_y_purity->GetXaxis()->SetTitleOffset(1.2);
    h_y_purity->GetYaxis()->SetTitleOffset(1.5);
    h_y_purity->GetXaxis()->SetLabelSize(0.04);  
    h_y_purity->GetYaxis()->SetLabelSize(0.04);
    h_y_purity->GetXaxis()->SetTitleSize(0.04);  
    h_y_purity->GetYaxis()->SetTitleSize(0.04);
    h_y_purity->GetXaxis()->SetTitle("y_{reco bin}");
    h_y_purity->GetYaxis()->SetTitle("Purity");
    h_y_purity->Draw("same");

    double binWidth_yPurity = h_y_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label24_2 = new TLatex(0.5,0.81, Form("bin width=%.4f",binWidth_yPurity));
    label24_2->SetNDC();
    label24_2->SetTextSize(0.04);
    label24_2->Draw("same");

    TH1D* h_y_stability = new TH1D("h_y_stability",";y_{true bin};stability",nx_y,0.01,0.85);
    for (int i=1; i<=nx_y; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_y; ++j) trueRowSum += h_y_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_y_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_y_stability->SetBinContent(i, diag / trueRowSum);
    }
    c24->cd(3);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    h_y_stability->GetXaxis()->SetTitleOffset(1.2);
    h_y_stability->GetYaxis()->SetTitleOffset(1.5);
    h_y_stability->GetXaxis()->SetLabelSize(0.04);  
    h_y_stability->GetYaxis()->SetLabelSize(0.04);
    h_y_stability->GetXaxis()->SetTitleSize(0.04);  
    h_y_stability->GetYaxis()->SetTitleSize(0.04);
    h_y_stability->GetYaxis()->SetTitle("Stability");
    h_y_stability->Draw("same");

    double binWidth_yStability = h_y_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label24_3 = new TLatex(0.5,0.81, Form("bin width=%.4f",binWidth_yStability));
    label24_3->SetNDC();
    label24_3->SetTextSize(0.04);
    label24_3->Draw("same");
    
    int nb_y = h_y_e->GetNbinsX();
    TH1D* h_y_efficiency = new TH1D("h_y_efficiency", ";y;Efficiency", nb_y, 0.01, 0.85);
    TH1D* h_y_acceptance = new TH1D("h_y_acceptance", ";y;Acceptance", nb_y, 0.01, 0.85);
    TH1D* h_y_corrected = new TH1D("h_y_corrected", ";y;Acceptance Corrected", nb_y, 0.01, 0.85);
    for (int i=1; i<=nb_y; ++i) 
    {
        double Ntrue = h_y_e->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_y; ++j) Naccepted += h_y_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_y_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_yREC_e_EEMC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_y_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_y_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_y_efficiency->SetBinContent(i, efficiency);
    }
    c24->cd(4);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    h_y_acceptance->GetXaxis()->SetTitleOffset(1.2);
    h_y_acceptance->GetYaxis()->SetTitleOffset(1.5);
    h_y_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_y_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_y_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_y_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_y_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_y_acceptance->Draw("same");

    c24->cd(5);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    h_y_efficiency->GetXaxis()->SetTitleOffset(1.2);
    h_y_efficiency->GetYaxis()->SetTitleOffset(1.8);
    h_y_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_y_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_y_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_y_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_y_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_y_efficiency->GetXaxis()->SetTitle("y");
    h_y_efficiency->Draw("same");
    
    c24->cd(6);
    gPad->SetLogy(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    h_y_corrected->GetXaxis()->SetTitleOffset(1.2);
    h_y_corrected->GetYaxis()->SetTitleOffset(1.5);
    h_y_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_y_corrected->GetYaxis()->SetLabelSize(0.04);
    h_y_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_y_corrected->GetYaxis()->SetTitleSize(0.04);
    h_y_corrected->GetYaxis()->SetTitle("Corrected");
    h_y_corrected->Draw("same");

    c24->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_1 = new TLatex();
    title24_1->SetNDC(); 
    title24_1->SetTextSize(0.05);
    title24_1->SetTextAlign(22);  
    title24_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c24->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_3 = new TLatex();
    title24_3->SetNDC(); 
    title24_3->SetTextSize(0.05);
    title24_3->SetTextAlign(22);  
    title24_3->DrawLatex(0.5, 0.97, "Purity");  

    c24->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_4 = new TLatex();
    title24_4->SetNDC(); 
    title24_4->SetTextSize(0.05);
    title24_4->SetTextAlign(22);  
    title24_4->DrawLatex(0.5, 0.97, "Stability");

    c24->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_5 = new TLatex();
    title24_5->SetNDC(); 
    title24_5->SetTextSize(0.05);
    title24_5->SetTextAlign(22);  
    title24_5->DrawLatex(0.5, 0.97, "Acceptance");

    c24->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_2 = new TLatex();
    title24_2->SetNDC(); 
    title24_2->SetTextSize(0.05);
    title24_2->SetTextAlign(22);  
    title24_2->DrawLatex(0.5, 0.97, "Efficiency");  


    c24->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title24_6 = new TLatex();
    title24_6->SetNDC(); 
    title24_6->SetTextSize(0.05);
    title24_6->SetTextAlign(22);  
    title24_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");
    
    c24->Print("./figures/analysisNote_plots/plot_y_migrationCorrectedPurityAcceptance.pdf");
}

void plot_x()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_x_res_cut_beforeCuts = (TH2D*) file->Get("h_x_res_cut_beforeCuts");
    TH2D* h_x_response_cut_beforeCuts = (TH2D*) file->Get("h_x_response_cut_beforeCuts");
    TH1D* h_x = (TH1D*) file->Get("h_x");
    TH1D* h_x_REC = (TH1D*) file->Get("h_x_REC");

    TCanvas* c1_123 = new TCanvas("c1_123","c1_123",1,1,1600,600);
    c1_123->Divide(3,1,0.01,0.01);
    c1_123->cd(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0); 
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    h_x->GetXaxis()->SetTitleOffset(1.2);
    h_x->GetYaxis()->SetTitleOffset(1.5);
    h_x->GetXaxis()->SetLabelSize(0.04);  
    h_x->GetYaxis()->SetLabelSize(0.04);
    h_x->GetXaxis()->SetTitleSize(0.04);  
    h_x->GetYaxis()->SetTitleSize(0.04);
	h_x->GetYaxis()->SetTitle("counts");
	h_x->GetXaxis()->SetTitle("x");
    h_x->GetXaxis()->SetRangeUser(0,0.3);
	h_x->SetLineColor(kBlack);
    h_x->Draw();
    h_x_REC->SetMarkerStyle(24);
    h_x_REC->SetMarkerColor(kBlue);
    h_x_REC->Draw("PEsame");
    TLegend *w1_123 = new TLegend(0.65,0.7,0.75,0.85);
	w1_123->AddEntry(h_x, " MC", "L");
    w1_123->AddEntry(h_x_REC, " RECO", "P");
    w1_123->SetBorderSize(0);   
    w1_123->SetFillStyle(0);
    w1_123->SetTextSize(0.05);
	w1_123->Draw("same");

    c1_123->cd(2);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0); 
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.12);
    h_x_res_cut_beforeCuts->GetXaxis()->SetNdivisions(505);
    h_x_res_cut_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_x_res_cut_beforeCuts->GetYaxis()->SetTitleOffset(1.7);
    h_x_res_cut_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_x_res_cut_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_x_res_cut_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_x_res_cut_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_x_res_cut_beforeCuts->GetYaxis()->SetRangeUser(-0.06,0.06);
    h_x_res_cut_beforeCuts->GetXaxis()->SetRangeUser(0,0.1);
    h_x_res_cut_beforeCuts->GetYaxis()->SetTitle("(x_{MC}-x_{RECO})/x_{MC}");
    h_x_res_cut_beforeCuts->GetXaxis()->SetTitle("x_{MC}");
    h_x_res_cut_beforeCuts->Draw("colzsame");

    c1_123->cd(3);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0); 
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.12);
    h_x_response_cut_beforeCuts->GetXaxis()->SetNdivisions(505);
    h_x_response_cut_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_x_response_cut_beforeCuts->GetYaxis()->SetTitleOffset(1.8);
    h_x_response_cut_beforeCuts->GetXaxis()->SetRangeUser(0,0.1);
    h_x_response_cut_beforeCuts->GetYaxis()->SetRangeUser(0,0.1);
    h_x_response_cut_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_x_response_cut_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_x_response_cut_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_x_response_cut_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_x_response_cut_beforeCuts->GetYaxis()->SetTitle("x_{MC}");
    h_x_response_cut_beforeCuts->GetXaxis()->SetTitle("x_{RECO}");
    h_x_response_cut_beforeCuts->Draw("colzsame");

    c1_123->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1_1_123 = new TLatex();
    title1_1_123->SetNDC(); 
    title1_1_123->SetTextSize(0.05);
    title1_1_123->SetTextAlign(22);  
    title1_1_123->DrawLatex(0.5, 0.97, "x Truth vs Reco");  

    c1_123->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_2_123 = new TLatex();
    title1_2_123->SetNDC(); 
    title1_2_123->SetTextSize(0.05);
    title1_2_123->SetTextAlign(22);  
    title1_2_123->DrawLatex(0.5, 0.97, "x Resolution");

    c1_123->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title1_3_123 = new TLatex();
    title1_3_123->SetNDC(); 
    title1_3_123->SetTextSize(0.05);
    title1_3_123->SetTextAlign(22);  
    title1_3_123->DrawLatex(0.5, 0.97, "x Response");

    c1_123->Print("./figures/analysisNote_plots/plot_x.pdf");
}

void x_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_x_MC_after = (TH1D*) file->Get("h_x_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_x_migration_beforeCuts = (TH2D*) file2->Get("h_x_migration_beforeCuts");
    TH1D* h_x_REC = (TH1D*) file2->Get("h_x_REC");
    TH1D* h_x = (TH1D*) file2->Get("h_x");

    TCanvas* c22 = new TCanvas("c22","c22",1,1,1400,800);
    c22->Divide(3,2,0.01,0.01);
    c22->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_x_migration_beforeCuts->GetXaxis()->SetNdivisions(505);
    h_x_migration_beforeCuts->GetYaxis()->SetTitleOffset(1.5);
    h_x_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_x_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_x_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_x_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_x_migration_beforeCuts->GetYaxis()->SetRangeUser(0,0.2);
    h_x_migration_beforeCuts->GetXaxis()->SetRangeUser(0,0.1);
    h_x_migration_beforeCuts->Draw("colzsame");   

    int nx_x = h_x_migration_beforeCuts->GetNbinsX();
    int ny_x = h_x_migration_beforeCuts->GetNbinsY();
    TH1D* h_x_purity = new TH1D("h_x_purity",";x_{reco bin};Purity",nx_x,0,0.5);
    for (int j=1; j<=ny_x; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_x; ++i) recoColSum += h_x_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_x_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_x_purity->SetBinContent(j, diag / recoColSum);
    }
    c22->cd(2);
    //gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_x_purity->GetXaxis()->SetNdivisions(505);
    h_x_purity->GetXaxis()->SetRangeUser(0,0.04);
    h_x_purity->GetXaxis()->SetLabelSize(0.04);  
    h_x_purity->GetYaxis()->SetLabelSize(0.04);
    h_x_purity->GetXaxis()->SetTitleSize(0.04);  
    h_x_purity->GetYaxis()->SetTitleSize(0.04);
    h_x_purity->GetXaxis()->SetTitle("x_{reco bin}");
    h_x_purity->GetYaxis()->SetTitle("Purity");
    h_x_purity->Draw("same");

    double binWidth_xPurity = h_x_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label22_2 = new TLatex(0.32,0.81, Form("bin width=%.3f",binWidth_xPurity));
    label22_2->SetNDC();
    label22_2->SetTextSize(0.04);
    label22_2->Draw("same");

    TH1D* h_x_stability = new TH1D("h_x_stability",";x_{true bin};Stability",nx_x,0,0.5);
    for (int i=1; i<=nx_x; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_x; ++j) trueRowSum += h_x_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_x_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_x_stability->SetBinContent(i, diag / trueRowSum);
    }
    c22->cd(3);
    //gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_x_stability->GetXaxis()->SetNdivisions(505);
    h_x_stability->GetXaxis()->SetRangeUser(0,0.2);
    h_x_stability->GetYaxis()->SetTitle("Stability");
    h_x_stability->GetXaxis()->SetLabelSize(0.04);  
    h_x_stability->GetYaxis()->SetLabelSize(0.04);
    h_x_stability->GetXaxis()->SetTitleSize(0.04);  
    h_x_stability->GetYaxis()->SetTitleSize(0.04);
    h_x_stability->Draw("same");

    double binWidth_xStability = h_x_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label22_3 = new TLatex(0.32,0.81, Form("bin width=%.3f",binWidth_xStability));
    label22_3->SetNDC();
    label22_3->SetTextSize(0.04);
    label22_3->Draw("same");
    
    int nb_x = h_x->GetNbinsX();
    TH1D* h_x_efficiency = new TH1D("h_x_efficiency", ";x;Efficiency", nb_x, 0, 0.5);
    TH1D* h_x_acceptance = new TH1D("h_x_acceptance", ";x;Acceptance", nb_x, 0, 0.5);
    TH1D* h_x_corrected = new TH1D("h_x_corrected", ";x;Acceptance Corrected", nb_x, 0, 0.5);
    for (int i=1; i<=nb_x; ++i) 
    {
        double Ntrue = h_x->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_x; ++j) Naccepted += h_x_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_x_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_x_REC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_x_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_x_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_x_efficiency->SetBinContent(i, efficiency);
    }

    c22->cd(4);
    //gPad->SetLogx(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_x_acceptance->GetXaxis()->SetNdivisions(505);
    h_x_acceptance->GetXaxis()->SetRangeUser(0,0.2);
    h_x_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_x_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_x_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_x_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_x_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_x_acceptance->Draw("same");

    c22->cd(5);
    //gPad->SetLogx(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_x_efficiency->GetXaxis()->SetNdivisions(505);
    h_x_efficiency->GetXaxis()->SetRangeUser(0,0.2);
    h_x_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_x_efficiency->GetXaxis()->SetTitle("x");
    h_x_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_x_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_x_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_x_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_x_efficiency->Draw("same");
    
    c22->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_x_corrected->GetXaxis()->SetNdivisions(505);
    h_x_corrected->GetXaxis()->SetRangeUser(0,0.2);
    h_x_corrected->GetYaxis()->SetTitle("Corrected");
    h_x_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_x_corrected->GetYaxis()->SetLabelSize(0.04);
    h_x_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_x_corrected->GetYaxis()->SetTitleSize(0.04);
    h_x_corrected->Draw("same");

    c22->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_1 = new TLatex();
    title22_1->SetNDC(); 
    title22_1->SetTextSize(0.05);
    title22_1->SetTextAlign(22);  
    title22_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c22->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_3 = new TLatex();
    title22_3->SetNDC(); 
    title22_3->SetTextSize(0.05);
    title22_3->SetTextAlign(22);  
    title22_3->DrawLatex(0.5, 0.97, "Purity");  

    c22->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_4 = new TLatex();
    title22_4->SetNDC(); 
    title22_4->SetTextSize(0.05);
    title22_4->SetTextAlign(22);  
    title22_4->DrawLatex(0.5, 0.97, "Stability");

    c22->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_5 = new TLatex();
    title22_5->SetNDC(); 
    title22_5->SetTextSize(0.05);
    title22_5->SetTextAlign(22);  
    title22_5->DrawLatex(0.5, 0.97, "Acceptance");

    c22->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_2 = new TLatex();
    title22_2->SetNDC(); 
    title22_2->SetTextSize(0.05);
    title22_2->SetTextAlign(22);  
    title22_2->DrawLatex(0.5, 0.97, "Efficiency");  


    c22->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title22_6 = new TLatex();
    title22_6->SetNDC(); 
    title22_6->SetTextSize(0.05);
    title22_6->SetTextAlign(22);  
    title22_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");

    c22->Print("./figures/analysisNote_plots/plot_x_migrationCorrectedPurityAcceptance.pdf");
}

void plot_VM_pT()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_VM_pt_MC = (TH1D*) file->Get("h_VM_pt_MC");
    TH1D* h_VM_pt_REC = (TH1D*) file->Get("h_VM_pt_REC");
    TH2D* h_VM_pt_response = (TH2D*) file->Get("h_VM_pt_response");
    TH2D* h_VM_pt_resolution = (TH2D*) file->Get("h_VM_pt_resolution");

    TCanvas* c8 = new TCanvas("c8","c8",1,1,1600,600);
    c8->Divide(3,1,0.01,0.01);
    c8->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h_VM_pt_MC->GetXaxis()->SetTitleOffset(1.2);
    h_VM_pt_MC->GetYaxis()->SetTitleOffset(1.5);
    h_VM_pt_MC->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_MC->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_MC->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_MC->GetYaxis()->SetTitleSize(0.04);
	h_VM_pt_MC->GetYaxis()->SetTitle("counts");
	h_VM_pt_MC->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
	h_VM_pt_MC->SetLineColor(kBlack);
    h_VM_pt_MC->Draw();
    h_VM_pt_REC->SetMarkerColor(kBlue);
    h_VM_pt_REC->SetMarkerStyle(24);
    h_VM_pt_REC->Draw("PEsame");
    TLegend *w8 = new TLegend(0.65,0.7,0.75,0.8);
	w8->AddEntry(h_VM_pt_MC, " MC", "L");
    w8->AddEntry(h_VM_pt_REC, " RECO", "P");
    w8->SetBorderSize(0);   
    w8->SetFillStyle(0);
    w8->SetTextSize(0.05);
	w8->Draw("same");

    c8->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_VM_pt_resolution->GetYaxis()->SetRangeUser(-0.1,0.1);
    h_VM_pt_resolution->GetXaxis()->SetRangeUser(0.5,3.1);
    h_VM_pt_resolution->GetXaxis()->SetTitleOffset(1.2);
    h_VM_pt_resolution->GetYaxis()->SetTitleOffset(2);
    h_VM_pt_resolution->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_resolution->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_resolution->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_resolution->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_resolution->GetXaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
    h_VM_pt_resolution->Draw("colzsame");

    c8->cd(3);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.12);
    gStyle->SetPalette(kBird); 
    gStyle->SetPadRightMargin(0.15); 
    gStyle->SetLabelSize(0.03,"Z");  
    gStyle->SetTitleSize(0.03,"Z"); 
    gStyle->SetOptStat(0);
    h_VM_pt_response->GetXaxis()->SetRangeUser(0.5,3.1);
    h_VM_pt_response->GetYaxis()->SetRangeUser(0.5,3.1);
    h_VM_pt_response->GetXaxis()->SetTitleOffset(1.2);
    h_VM_pt_response->GetYaxis()->SetTitleOffset(1.6);
    h_VM_pt_response->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_response->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_response->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_response->GetYaxis()->SetTitleSize(0.04);
	h_VM_pt_response->GetYaxis()->SetTitle("p_{T,VM,MC} [GeV/c]");
	h_VM_pt_response->GetXaxis()->SetTitle("p_{T,VM,RECO} [GeV/c]");
    h_VM_pt_response->Draw("colz");
    gPad->Update(); 
    
    c8->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title8 = new TLatex();
    title8->SetNDC(); 
    title8->SetTextSize(0.05);
    title8->SetTextAlign(22);  
    title8->DrawLatex(0.5, 0.97, "p_{T} Truth vs Reco (VM)");  

    c8->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_3 = new TLatex();
    title8_3->SetNDC(); 
    title8_3->SetTextSize(0.05);
    title8_3->SetTextAlign(22);  
    title8_3->DrawLatex(0.5, 0.97, "p_{T} Resolution (VM)"); 
    
    c8->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title8_2 = new TLatex();
    title8_2->SetNDC(); 
    title8_2->SetTextSize(0.05);
    title8_2->SetTextAlign(22);  
    title8_2->DrawLatex(0.5, 0.97, "p_{T} Response (VM)");

    c8->Print("./figures/analysisNote_plots/plot_VM_pt.pdf");
}

void VM_pT_QA()
{     
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_VM_pt_MC_after = (TH1D*) file->Get("h_VM_pt_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_VM_pt_migration = (TH2D*) file2->Get("h_VM_pt_migration"); 
    TH1D* h_VM_pt_REC = (TH1D*) file2->Get("h_VM_pt_REC");
    TH1D* h_VM_pt_MC = (TH1D*) file2->Get("h_VM_pt_MC");

    TCanvas* c20 = new TCanvas("c20","c20",1,1,1400,800);
    c20->Divide(3,2,0.01,0.01);
    c20->cd(1);   
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogz(1);
    h_VM_pt_migration->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_migration->GetYaxis()->SetLabelSize(0.04); 
    h_VM_pt_migration->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_migration->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_migration->GetXaxis()->SetRangeUser(0.5,3.1);
    h_VM_pt_migration->GetYaxis()->SetRangeUser(0.5,3.1);
    h_VM_pt_migration->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_migration->GetYaxis()->SetTitleOffset(1.2); 
    h_VM_pt_migration->Draw("colzsame");

    int nx_VMpt = h_VM_pt_migration->GetNbinsX();
    int ny_VMpt = h_VM_pt_migration->GetNbinsY();
    TH1D* h_VM_pt_purity = new TH1D("h_VM_pt_purity",";p_{T,VM,reco bin} [GeV/c];Purity",nx_VMpt,0,10);
    for (int j=1; j<=ny_VMpt; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_VMpt; ++i) recoColSum += h_VM_pt_migration->GetBinContent(i,j);
        double diag = h_VM_pt_migration->GetBinContent(j,j);
        if (recoColSum>0) h_VM_pt_purity->SetBinContent(j, diag / recoColSum);
    }
    c20->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_pt_purity->GetXaxis()->SetRangeUser(0,3.3);
    h_VM_pt_purity->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_purity->GetYaxis()->SetLabelSize(0.04); 
    h_VM_pt_purity->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_purity->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_purity->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_purity->GetYaxis()->SetTitleOffset(1.2);
    h_VM_pt_purity->GetXaxis()->SetTitle("p_{T,VM,reco bin} [GeV/c]");
    h_VM_pt_purity->GetYaxis()->SetTitle("Purity");
    h_VM_pt_purity->Draw("same");

    double binWidth_VMptPurity = h_VM_pt_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label20_2 = new TLatex(0.35,0.21, Form("bin width=%.2f GeV/c",binWidth_VMptPurity));
    label20_2->SetNDC();
    label20_2->SetTextSize(0.04);
    label20_2->Draw("same");

    TH1D* h_VM_pt_stability = new TH1D("h_VM_pt_stability",";p_{T,VM,true bin} [GeV/c];Stability",nx_VMpt,0,10);
    for (int i=1; i<=nx_VMpt; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_VMpt; ++j) trueRowSum += h_VM_pt_migration->GetBinContent(i,j);
        double diag = h_VM_pt_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_VM_pt_stability->SetBinContent(i, diag / trueRowSum);
    }
    c20->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_pt_stability->GetXaxis()->SetRangeUser(0,3.3);
    h_VM_pt_stability->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_stability->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_stability->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_stability->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_stability->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_stability->GetYaxis()->SetTitleOffset(1.2);
    h_VM_pt_stability->Draw("same");

    double binWidth_VMptStability = h_VM_pt_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label20_3 = new TLatex(0.35,0.21, Form("bin width=%.2f GeV/c",binWidth_VMptStability));
    label20_3->SetNDC();
    label20_3->SetTextSize(0.04);
    label20_3->Draw("same");
    
    int nb_VMpt = h_VM_pt_MC->GetNbinsX();
    TH1D* h_VM_pt_efficiency = new TH1D("h_VM_pt_efficiency", ";p_{T,VM} [GeV/c];Efficiency", nb_VMpt, 0, 10);
    TH1D* h_VM_pt_acceptance = new TH1D("h_VM_pt_acceptance", ";p_{T,VM} [GeV/c];Acceptance", nb_VMpt, 0, 10);
    TH1D* h_VM_pt_corrected = new TH1D("h_VM_pt_corrected", ";p_{T,VM} [GeV/c];Acceptance Corrected", nb_VMpt, 0, 10);
    for (int i=1; i<=nb_VMpt; ++i) 
    {
        double Ntrue = h_VM_pt_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_VMpt; ++j) Naccepted += h_VM_pt_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_VM_pt_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_VM_pt_REC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_VM_pt_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_VM_pt_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_VM_pt_efficiency->SetBinContent(i, efficiency);
    }
    c20->cd(4);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_pt_acceptance->GetXaxis()->SetRangeUser(0,3.3);
    h_VM_pt_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_acceptance->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_acceptance->GetYaxis()->SetTitleOffset(1.2);
    h_VM_pt_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_VM_pt_acceptance->Draw("same");

    c20->cd(5);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_pt_efficiency->GetXaxis()->SetRangeUser(0,3.3);
    h_VM_pt_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_efficiency->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_efficiency->GetYaxis()->SetTitleOffset(1.2);
    h_VM_pt_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_VM_pt_efficiency->GetXaxis()->SetTitle("p_{T,VM} [GeV/c]");
    h_VM_pt_efficiency->Draw("same");
    
    c20->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_pt_corrected->GetXaxis()->SetRangeUser(0,3.3);
    h_VM_pt_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_VM_pt_corrected->GetYaxis()->SetLabelSize(0.04);
    h_VM_pt_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_VM_pt_corrected->GetYaxis()->SetTitleSize(0.04);
    h_VM_pt_corrected->GetXaxis()->SetTitleOffset(1.2); 
    h_VM_pt_corrected->GetYaxis()->SetTitleOffset(1.2);
    h_VM_pt_corrected->GetYaxis()->SetTitle("Corrected");
    h_VM_pt_corrected->Draw("same");

    c20->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_1 = new TLatex();
    title20_1->SetNDC(); 
    title20_1->SetTextSize(0.05);
    title20_1->SetTextAlign(22);  
    title20_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c20->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_3 = new TLatex();
    title20_3->SetNDC(); 
    title20_3->SetTextSize(0.05);
    title20_3->SetTextAlign(22);  
    title20_3->DrawLatex(0.5, 0.97, "Purity");  

    c20->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_4 = new TLatex();
    title20_4->SetNDC(); 
    title20_4->SetTextSize(0.05);
    title20_4->SetTextAlign(22);  
    title20_4->DrawLatex(0.5, 0.97, "Stability");

    c20->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_5 = new TLatex();
    title20_5->SetNDC(); 
    title20_5->SetTextSize(0.05);
    title20_5->SetTextAlign(22);  
    title20_5->DrawLatex(0.5, 0.97, "Acceptance");

    c20->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_2 = new TLatex();
    title20_2->SetNDC(); 
    title20_2->SetTextSize(0.05);
    title20_2->SetTextAlign(22);  
    title20_2->DrawLatex(0.5, 0.97, "Efficiency");  


    c20->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title20_6 = new TLatex();
    title20_6->SetNDC(); 
    title20_6->SetTextSize(0.05);
    title20_6->SetTextAlign(22);  
    title20_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");
    
    c20->Print("./figures/analysisNote_plots/plot_VM_migrationCorrectedPurityAcceptance.pdf");
}

void plot_VM_Epz()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_VM_Epz_MC = (TH1D*) file->Get("h_VM_Epz_MC");
    TH1D* h_VM_Epz_REC = (TH1D*) file->Get("h_VM_Epz_REC");
    TH2D* h_VM_Epz_res = (TH2D*) file->Get("h_VM_Epz_res");
    TH2D* h_VM_Epz_response = (TH2D*) file->Get("h_VM_Epz_response");

    TCanvas* c7 = new TCanvas("c7","c7",1,1,1600,600);
    c7->Divide(3,1,0.01,0.01);
    c7->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.15);
    h_VM_Epz_MC->GetXaxis()->SetTitleOffset(1.2);
    h_VM_Epz_MC->GetYaxis()->SetTitleOffset(1.8);
	h_VM_Epz_MC->GetYaxis()->SetTitle("counts");
	h_VM_Epz_MC->GetXaxis()->SetTitle("(E-p_{z})_{VM} [GeV]");
	h_VM_Epz_MC->SetLineColor(kBlack);
    h_VM_Epz_MC->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_MC->GetYaxis()->SetLabelSize(0.04);
    h_VM_Epz_MC->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_MC->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_MC->Draw();
    h_VM_Epz_REC->SetMarkerStyle(24);
    h_VM_Epz_REC->SetMarkerColor(kBlue);
    h_VM_Epz_REC->Draw("PEsame");
    TLegend *w7 = new TLegend(0.65,0.75,0.75,0.85);
	w7->AddEntry(h_VM_Epz_MC, "MC", "L");
    w7->AddEntry(h_VM_Epz_REC, "RECO", "P");
    w7->SetBorderSize(0);   
    w7->SetFillStyle(0);
    w7->SetTextSize(0.05);
	w7->Draw("same");

    c7->cd(2);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gStyle->SetOptStat(0);
    h_VM_Epz_res->GetYaxis()->SetRangeUser(-0.1,0.1);
    h_VM_Epz_res->GetXaxis()->SetRangeUser(0,7);
    h_VM_Epz_res->GetXaxis()->SetTitleOffset(1.2);
    h_VM_Epz_res->GetYaxis()->SetTitleOffset(2.2);
    h_VM_Epz_res->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_res->GetYaxis()->SetLabelSize(0.04);
    h_VM_Epz_res->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_res->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_res->GetXaxis()->SetTitle("(E-p_{z})_{VM,MC} [GeV]");
    h_VM_Epz_res->GetYaxis()->SetTitle("[(E-p_{z})_{VM,RECO} - (E-p_{z})_{VM,MC}]/(E-p_{z})_{VM,MC}");
    h_VM_Epz_res->Draw("colzsame");

    c7->cd(3);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gStyle->SetOptStat(0);
    h_VM_Epz_response->GetXaxis()->SetTitleOffset(1.2);
    h_VM_Epz_response->GetYaxis()->SetTitleOffset(1.8);
    h_VM_Epz_response->GetXaxis()->SetRangeUser(0,7);
    h_VM_Epz_response->GetYaxis()->SetRangeUser(0,7);
    h_VM_Epz_response->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_response->GetYaxis()->SetLabelSize(0.04);
    h_VM_Epz_response->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_response->GetYaxis()->SetTitleSize(0.04);
	h_VM_Epz_response->GetYaxis()->SetTitle("(E-p_{z})_{VM,MC} [GeV]");
	h_VM_Epz_response->GetXaxis()->SetTitle("(E-p_{z})_{VM,RECO} [GeV]");
    h_VM_Epz_response->Draw("colzsame");

    c7->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title7 = new TLatex();
    title7->SetNDC(); 
    title7->SetTextSize(0.05);
    title7->SetTextAlign(22);  
    title7->DrawLatex(0.5, 0.97, "E-p_{z} Truth vs Reco (VM)");  

    c7->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_3 = new TLatex();
    title7_3->SetNDC(); 
    title7_3->SetTextSize(0.05);
    title7_3->SetTextAlign(22);  
    title7_3->DrawLatex(0.5, 0.97, "E-p_{z} Resolution (VM)");  

    c7->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title7_2 = new TLatex();
    title7_2->SetNDC(); 
    title7_2->SetTextSize(0.05);
    title7_2->SetTextAlign(22);  
    title7_2->DrawLatex(0.5, 0.97, "E-p_{z} Response (VM)");

    c7->Print("./figures/analysisNote_plots/plot_VM_Epz.pdf");
}

void VM_Epz_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_VM_Epz_MC_after = (TH1D*) file->Get("h_VM_Epz_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_VM_Epz_migration = (TH2D*) file2->Get("h_VM_Epz_migration");
    TH1D* h_VM_Epz_REC = (TH1D*) file2->Get("h_VM_Epz_REC");
    TH1D* h_VM_Epz_MC = (TH1D*) file2->Get("h_VM_Epz_MC");

    TCanvas* c21 = new TCanvas("c21","c21",1,1,1400,800);
    c21->Divide(3,2,0.01,0.01);
    c21->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetLogz(1);
    h_VM_Epz_migration->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_migration->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_migration->GetXaxis()->SetRangeUser(0,7);
    h_VM_Epz_migration->GetYaxis()->SetRangeUser(0,7);
    h_VM_Epz_migration->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_migration->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_migration->GetXaxis()->SetTitleOffset(1.5);  
    h_VM_Epz_migration->GetXaxis()->SetTitle("(E-p_{z})_{VM,MC} [GeV]");
    h_VM_Epz_migration->GetYaxis()->SetTitle("(E-p_{z})_{VM,RECO} [GeV]");
    h_VM_Epz_migration->Draw("colzsame");

    int nx_VMEpz = h_VM_Epz_migration->GetNbinsX();
    int ny_VMEpz = h_VM_Epz_migration->GetNbinsY();
    TH1D* h_VM_Epz_purity = new TH1D("h_VM_Epz_purity",";VM (E-p_{z})_{VM,reco bin};Purity",nx_VMEpz,0,20);
    for (int j=1; j<=ny_VMEpz; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_VMEpz; ++i) recoColSum += h_VM_Epz_migration->GetBinContent(i,j);
        double diag = h_VM_Epz_migration->GetBinContent(j,j);
        if (recoColSum>0) h_VM_Epz_purity->SetBinContent(j, diag / recoColSum);
    }
    c21->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_Epz_purity->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_purity->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_purity->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_purity->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_purity->GetXaxis()->SetTitleOffset(1.5); 
    h_VM_Epz_purity->GetXaxis()->SetTitle("(E-p_{z})_{VM,reco bin} [GeV]");
    h_VM_Epz_purity->GetYaxis()->SetTitle("Purity");
    h_VM_Epz_purity->Draw("same");

    double binWidth_VMEpzPurity = h_VM_Epz_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label21_2 = new TLatex(0.5,0.81, Form("bin width=%.1f GeV",binWidth_VMEpzPurity));
    label21_2->SetNDC();
    label21_2->SetTextSize(0.04);
    label21_2->Draw("same");


    TH1D* h_VM_Epz_stability = new TH1D("h_VM_Epz_stability",";VM (E-p_{z})_{VM,true bin};Stability",nx_VMEpz,0,20);
    for (int i=1; i<=nx_VMEpz; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_VMEpz; ++j) trueRowSum += h_VM_Epz_migration->GetBinContent(i,j);
        double diag = h_VM_Epz_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_VM_Epz_stability->SetBinContent(i, diag / trueRowSum);
    }
    c21->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_Epz_stability->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_stability->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_stability->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_stability->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_stability->GetXaxis()->SetTitleOffset(1.5); 
    h_VM_Epz_stability->GetXaxis()->SetTitle("(E-p_{z})_{VM,true bin} [GeV]");
    h_VM_Epz_stability->GetYaxis()->SetTitle("Stability");
    h_VM_Epz_stability->Draw("same");

    double binWidth_VMEpzStability = h_VM_Epz_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label21_3 = new TLatex(0.5,0.81, Form("bin width=%.1f GeV",binWidth_VMEpzStability));
    label21_3->SetNDC();
    label21_3->SetTextSize(0.04);
    label21_3->Draw("same");
    
    int nb_VMEpz = h_VM_Epz_MC->GetNbinsX();
    TH1D* h_VM_Epz_efficiency = new TH1D("h_VM_Epz_efficiency", ";(E-p_{z})_{VM} [GeV];Efficiency", nb_VMEpz, 0, 20);
    TH1D* h_VM_Epz_acceptance = new TH1D("h_VM_Epz_acceptance", ";(E-p_{z})_{VM} [GeV];Acceptance", nb_VMEpz, 0, 20);
    TH1D* h_VM_Epz_corrected = new TH1D("h_VM_Epz_corrected", ";(E-p_{z})_{VM} [GeV];Acceptance Corrected", nb_VMEpz, 0, 20);
    for (int i=1; i<=nb_VMEpz; ++i) 
    {
        double Ntrue = h_VM_Epz_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_VMEpz; ++j) Naccepted += h_VM_Epz_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_VM_Epz_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_VM_Epz_REC->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_VM_Epz_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_VM_Epz_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_VM_Epz_efficiency->SetBinContent(i, efficiency);
    }
    c21->cd(4);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_Epz_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_acceptance->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_acceptance->GetXaxis()->SetTitleOffset(1.5); 
    h_VM_Epz_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_VM_Epz_acceptance->GetXaxis()->SetTitle("(E-p_{z})_{VM} [GeV]");
    h_VM_Epz_acceptance->Draw("same");

    c21->cd(5);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_Epz_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_efficiency->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_efficiency->GetXaxis()->SetTitleOffset(1.5); 
    h_VM_Epz_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_VM_Epz_efficiency->GetXaxis()->SetTitle("(E-p_{z})_{VM} [GeV]");
    h_VM_Epz_efficiency->Draw("same");
    
    c21->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    h_VM_Epz_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_VM_Epz_corrected->GetYaxis()->SetLabelSize(0.04); 
    h_VM_Epz_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_VM_Epz_corrected->GetYaxis()->SetTitleSize(0.04);
    h_VM_Epz_corrected->GetXaxis()->SetTitleOffset(1.5); 
    h_VM_Epz_corrected->GetYaxis()->SetTitle("Corrected");
    h_VM_Epz_corrected->GetXaxis()->SetTitle("(E-p_{z})_{VM} [GeV]");
    h_VM_Epz_corrected->Draw("same");

    c21->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_1 = new TLatex();
    title21_1->SetNDC(); 
    title21_1->SetTextSize(0.05);
    title21_1->SetTextAlign(22);  
    title21_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c21->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_3 = new TLatex();
    title21_3->SetNDC(); 
    title21_3->SetTextSize(0.05);
    title21_3->SetTextAlign(22);  
    title21_3->DrawLatex(0.5, 0.97, "Purity");  

    c21->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_4 = new TLatex();
    title21_4->SetNDC(); 
    title21_4->SetTextSize(0.05);
    title21_4->SetTextAlign(22);  
    title21_4->DrawLatex(0.5, 0.97, "Stability");

    c21->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_5 = new TLatex();
    title21_5->SetNDC(); 
    title21_5->SetTextSize(0.05);
    title21_5->SetTextAlign(22);  
    title21_5->DrawLatex(0.5, 0.97, "Acceptance");

    c21->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_2 = new TLatex();
    title21_2->SetNDC(); 
    title21_2->SetTextSize(0.05);
    title21_2->SetTextAlign(22);  
    title21_2->DrawLatex(0.5, 0.97, "Efficiency");  


    c21->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title21_6 = new TLatex();
    title21_6->SetNDC(); 
    title21_6->SetTextSize(0.05);
    title21_6->SetTextAlign(22);  
    title21_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");
    
    c21->Print("./figures/analysisNote_plots/plot_VM_Epz_migrationCorrectedPurityAcceptance.pdf");
}

void plot_t_2d_wRES()
{
    TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_REC_2d_wRES_cut = (TH1D*) phi_t_file->Get("h_t_REC_2d_wRES_cut");

    TCanvas* c16 = new TCanvas("c16","c16",1,1,800,800);
    c16->Divide(1,1,0.01,0.01);
    c16->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.11);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz(1);
    h_t_REC_2d_wRES_cut->GetYaxis()->SetLabelSize(0.03); 
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitleOffset(1.6);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetTitleOffset(1.4);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetLabelSize(0.03);
    h_t_REC_2d_wRES_cut->GetXaxis()->SetTitle("#Delta_{x} = #sqrt{|t|_{x}} [GeV/c]");
    h_t_REC_2d_wRES_cut->GetYaxis()->SetTitle("#Delta_{y} = #sqrt{|t|_{y}} = #sqrt{|t|_{#hat{n}}} [GeV/c]");
    gPad->SetLogz(1);
    h_t_REC_2d_wRES_cut->Draw("colzsame");

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    double length = 0.14;
    double omega  = M_PI/12; 
    double x_end  = length * sin(omega);
    double y_end  = length * cos(omega);

    TArrow *arrow = new TArrow(0, 0, x_end, y_end, 0.02, "|>");
    arrow->SetLineColor(kBlack);
    arrow->SetLineWidth(2);
    arrow->Draw();

    double radius = 0.1;
    TEllipse *arc = new TEllipse(0, 0, radius, radius, 75, 90); 
    arc->SetFillStyle(0);
    arc->SetLineColor(kBlack);
    arc->SetLineWidth(2);
    arc->Draw();

    TLatex* r421111 = new TLatex(0.14, 0.56, "#omega_{max}"); 
	r421111->SetNDC();
	r421111->SetTextSize(30);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");
    gPad->Update();

    gStyle->SetOptStat(0);
    gPad->SetLogz(1);

    c16->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title16 = new TLatex();
    title16->SetNDC(); 
    title16->SetTextSize(0.05);
    title16->SetTextAlign(22);  
    title16->DrawLatex(0.5, 0.97, "2D |t| Distribution with Detector Resolution");  

    c16->Print("./figures/analysisNote_plots/plot_t_2D_withRes.pdf");
}

void plot_e_Epz()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_Epz_MC = (TH1D*) file->Get("h_Epz_MC");
    TH1D* h_Epz_beforeCut = (TH1D*) file->Get("h_Epz_beforeCut");
    TH2D* h_Epz_response_beforeCuts = (TH2D*) file->Get("h_Epz_response_beforeCuts"); 
    TH2D* h_Epz_res_beforeCuts = (TH2D*) file->Get("h_Epz_res_beforeCuts"); 

    TCanvas* c11 = new TCanvas("c11","c11",1,1,1600,600);
    c11->Divide(3,1,0.01,0.01);
    c11->cd(1);
    //gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetTopMargin(0.18);  
    gStyle->SetOptStat(0);
    h_Epz_MC->GetXaxis()->SetTitleOffset(1.2);
    h_Epz_MC->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_MC->GetXaxis()->SetRangeUser(15,25);
    h_Epz_MC->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_MC->GetYaxis()->SetLabelSize(0.04);
    h_Epz_MC->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_MC->GetYaxis()->SetTitleSize(0.04);
	h_Epz_MC->GetYaxis()->SetTitle("counts");
	h_Epz_MC->GetXaxis()->SetTitle("(E - p_{z}) [GeV]");
    h_Epz_MC->SetLineColor(kBlack);
    h_Epz_MC->Draw();
    h_Epz_beforeCut->SetMarkerColor(kBlue);
    h_Epz_beforeCut->SetMarkerStyle(24);
    h_Epz_beforeCut->Draw("PEsame");
    TLegend *w10 = new TLegend(0.65,0.7,0.75,0.8);
	w10->AddEntry(h_Epz_MC, " MC", "L");
    w10->AddEntry(h_Epz_beforeCut, " RECO", "P");
    w10->SetBorderSize(0);   
    w10->SetFillStyle(0);
    w10->SetTextSize(0.05);
	w10->Draw("same");

    c11->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    //gPad->SetBottomMargin(0.15);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_Epz_res_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_Epz_res_beforeCuts->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_res_beforeCuts->GetXaxis()->SetRangeUser(19.45,20.1);
    h_Epz_res_beforeCuts->GetYaxis()->SetRangeUser(-0.3,0.3);
    h_Epz_res_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_res_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_Epz_res_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_res_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_Epz_res_beforeCuts->Draw("colzsame");

    c11->cd(3);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    //gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_response_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_Epz_response_beforeCuts->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_response_beforeCuts->GetYaxis()->SetRangeUser(19.45,20.1);
    h_Epz_response_beforeCuts->GetXaxis()->SetRangeUser(14,26);
    h_Epz_response_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_response_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_Epz_response_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_response_beforeCuts->GetYaxis()->SetTitleSize(0.04);
	h_Epz_response_beforeCuts->GetYaxis()->SetTitle("(E_{MC}-p_{z,MC}) [GeV]");
	h_Epz_response_beforeCuts->GetXaxis()->SetTitle("(E_{EMCal} - p_{z,trk}) [GeV]");
    h_Epz_response_beforeCuts->Draw("colzsame");

    c11->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title11 = new TLatex();
    title11->SetNDC(); 
    title11->SetTextSize(0.05);
    title11->SetTextAlign(22);  
    title11->DrawLatex(0.6, 0.97, "(E-p_{z}) Truth vs. Reco (e'+HFS)");  

    c11->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_3 = new TLatex();
    title11_3->SetNDC(); 
    title11_3->SetTextSize(0.05);
    title11_3->SetTextAlign(22);  
    title11_3->DrawLatex(0.5, 0.97, "(E-p_{z}) Resolution (e'+HFS)");
    
    c11->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title11_2 = new TLatex();
    title11_2->SetNDC(); 
    title11_2->SetTextSize(0.05);
    title11_2->SetTextAlign(22);  
    title11_2->DrawLatex(0.5, 0.97, "(E-p_{z}) Response (e'+HFS)");
   
    c11->Print("./figures/analysisNote_plots/plot_Epz.pdf");
}

void e_Epz_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_Epz_MC_after = (TH1D*) file->Get("h_Epz_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_Epz_migration_noHFS = (TH2D*) file2->Get("h_Epz_migration_noHFS"); 
    TH1D* h_Epz_MC = (TH1D*) file2->Get("h_Epz_MC");
    TH1D* h_e_Epz_noHFS = (TH1D*) file2->Get("h_e_Epz_noHFS");

    TCanvas* c25 = new TCanvas("c25","c25",1,1,1400,800);
    c25->Divide(3,2,0.01,0.01);
    c25->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogz(1);
    h_Epz_migration_noHFS->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_migration_noHFS->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_migration_noHFS->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_migration_noHFS->GetYaxis()->SetLabelSize(0.04);
    h_Epz_migration_noHFS->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_migration_noHFS->GetYaxis()->SetTitleSize(0.04);
    h_Epz_migration_noHFS->GetXaxis()->SetRangeUser(19,20.5);
    h_Epz_migration_noHFS->GetYaxis()->SetRangeUser(14,26);
    h_Epz_migration_noHFS->GetXaxis()->SetTitle("(E-p_{z})_{MC}");
    h_Epz_migration_noHFS->GetYaxis()->SetTitle("(E-p_{z})_{RECO}");
    h_Epz_migration_noHFS->Draw("colzsame");

    int nx_Epz = h_Epz_migration_noHFS->GetNbinsX();
    int ny_Epz = h_Epz_migration_noHFS->GetNbinsY();
    TH1D* h_Epz_purity = new TH1D("h_Epz_purity",";(E-p_{z})_{reco bin};Purity",nx_Epz,0,45);
    for (int j=1; j<=ny_Epz; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_Epz; ++i) recoColSum += h_Epz_migration_noHFS->GetBinContent(i,j);

        double diag = h_Epz_migration_noHFS->GetBinContent(j,j);
        if (recoColSum>0) h_Epz_purity->SetBinContent(j, diag / recoColSum);
    }
    c25->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_purity->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_purity->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_purity->GetXaxis()->SetRangeUser(18,22);
    h_Epz_purity->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_purity->GetYaxis()->SetLabelSize(0.04);
    h_Epz_purity->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_purity->GetYaxis()->SetTitleSize(0.04);
    h_Epz_purity->GetXaxis()->SetTitle("(E-p_{z})_{reco bin} [GeV]");
    h_Epz_purity->GetYaxis()->SetTitle("Purity");
    h_Epz_purity->Draw("same");

    double binWidth_EpzPurity = h_Epz_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label25_2 = new TLatex(0.22,0.81, Form("bin width=%.2f GeV",binWidth_EpzPurity));
    label25_2->SetNDC();
    label25_2->SetTextSize(0.04);
    label25_2->Draw("same");

    TH1D* h_Epz_stability = new TH1D("h_Epz_stability",";(E-p_{z})_{true bin};stability",nx_Epz,0,45);
    for (int i=1; i<=nx_Epz; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_Epz; ++j) trueRowSum += h_Epz_migration_noHFS->GetBinContent(i,j);

        double diag = h_Epz_migration_noHFS->GetBinContent(i,i);
        if (trueRowSum>0) h_Epz_stability->SetBinContent(i, diag / trueRowSum);
    }
    c25->cd(3);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_stability->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_stability->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_stability->GetXaxis()->SetRangeUser(18,22);
    h_Epz_stability->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_stability->GetYaxis()->SetLabelSize(0.04);
    h_Epz_stability->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_stability->GetYaxis()->SetTitleSize(0.04);
    h_Epz_stability->GetYaxis()->SetTitle("Stability");
    h_Epz_stability->GetXaxis()->SetTitle("(E-p_{z})_{true bin} [GeV]");
    h_Epz_stability->Draw("same");

    double binWidth_EpzStability = h_Epz_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label25_3 = new TLatex(0.22,0.81, Form("bin width=%.2f GeV",binWidth_EpzStability));
    label25_3->SetNDC();
    label25_3->SetTextSize(0.04);
    label25_3->Draw("same");
    
    int nb_Epz_true = h_Epz_MC->GetNbinsX();
    int nb_Epz_reco = h_e_Epz_noHFS->GetNbinsY();
    TH1D* h_Epz_efficiency = new TH1D("h_Epz_efficiency", ";(E-p_{z})_{true bin};Efficiency", nb_Epz_true, 0, 45);
    TH1D* h_Epz_acceptance = new TH1D("h_Epz_acceptance", ";(E-p_{z}) [GeV];Acceptance", nb_Epz_true, 0, 45);
    TH1D* h_Epz_corrected = new TH1D("h_Epz_corrected", ";(E-p_{z}) [GeV];Acceptance Corrected", nb_Epz_true, 0, 45);
    for (int i=1; i<=nb_Epz_true; ++i) 
    {
        double Ntrue = h_Epz_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_Epz_reco; ++j) Naccepted += h_Epz_migration_noHFS->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_Epz_acceptance->SetBinContent(i, acceptance);

        double Nreco_in_bin_i = h_e_Epz_noHFS->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_Epz_corrected->SetBinContent(i, corrected);

        double Nreco_from_bin_i = h_Epz_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_Epz_efficiency->SetBinContent(i, efficiency);
    }
    c25->cd(4);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_acceptance->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_acceptance->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_acceptance->GetXaxis()->SetRangeUser(18,22);
    h_Epz_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_Epz_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_Epz_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_Epz_acceptance->GetXaxis()->SetTitle("(E-p_{z}) [GeV]");
    h_Epz_acceptance->Draw("same");

    c25->cd(5);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_efficiency->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_efficiency->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_efficiency->GetXaxis()->SetRangeUser(18,22);
    h_Epz_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_Epz_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_Epz_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_Epz_efficiency->GetXaxis()->SetTitle("(E-p_{z}) [GeV]");
    h_Epz_efficiency->Draw("same");
    
    c25->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(0);
    h_Epz_corrected->GetXaxis()->SetTitleOffset(1.5);
    h_Epz_corrected->GetYaxis()->SetTitleOffset(1.8);
    h_Epz_corrected->GetXaxis()->SetRangeUser(18,22);
    h_Epz_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_Epz_corrected->GetYaxis()->SetLabelSize(0.04);
    h_Epz_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_Epz_corrected->GetYaxis()->SetTitleSize(0.04);
    h_Epz_corrected->GetYaxis()->SetTitle("Corrected");
    h_Epz_corrected->GetXaxis()->SetTitle("(E-p_{z}) [GeV]");
    h_Epz_corrected->Draw("same");    

    c25->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_1 = new TLatex();
    title25_1->SetNDC(); 
    title25_1->SetTextSize(0.05);
    title25_1->SetTextAlign(22);  
    title25_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c25->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_3 = new TLatex();
    title25_3->SetNDC(); 
    title25_3->SetTextSize(0.05);
    title25_3->SetTextAlign(22);  
    title25_3->DrawLatex(0.5, 0.97, "Purity");  

    c25->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_4 = new TLatex();
    title25_4->SetNDC(); 
    title25_4->SetTextSize(0.05);
    title25_4->SetTextAlign(22);  
    title25_4->DrawLatex(0.5, 0.97, "Stability"); 

    c25->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_5 = new TLatex();
    title25_5->SetNDC(); 
    title25_5->SetTextSize(0.05);
    title25_5->SetTextAlign(22);  
    title25_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c25->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_2 = new TLatex();
    title25_2->SetNDC(); 
    title25_2->SetTextSize(0.05);
    title25_2->SetTextAlign(22);  
    title25_2->DrawLatex(0.5, 0.97, "Efficiency");  

    c25->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title25_6 = new TLatex();
    title25_6->SetNDC(); 
    title25_6->SetTextSize(0.05);
    title25_6->SetTextAlign(22);  
    title25_6->DrawLatex(0.5, 0.97, "Acceptance Corrected"); 
    
    c25->Print("./figures/analysisNote_plots/plot_Epz_migrationCorrectedPurityAcceptance.pdf");
}

void plot_e_EoP()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_EoverP_MC = (TH1D*) file->Get("h_EoverP_MC"); 
	TH1D* h_EoverP_beforeCut = (TH1D*) file->Get("h_EoverP_beforeCut");
    TH2D* h_EoverP_response_beforeCuts = (TH2D*) file->Get("h_EoverP_response_beforeCuts"); 
    TH2D* h_EoverP_res_beforeCuts = (TH2D*) file->Get("h_EoverP_res_beforeCuts"); 

    TCanvas* c12 = new TCanvas("c12","c12",1,1,1600,600);
    c12->Divide(3,1,0.01,0.01);
    c12->cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    gStyle->SetOptStat(0);
    h_EoverP_MC->GetXaxis()->SetTitleOffset(1.2);
    h_EoverP_MC->GetYaxis()->SetTitleOffset(2.5);
    h_EoverP_MC->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_MC->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_MC->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_MC->GetYaxis()->SetTitleSize(0.04);
	h_EoverP_MC->GetYaxis()->SetTitle("counts");
	h_EoverP_MC->GetXaxis()->SetTitle("E/|p|");
    h_EoverP_MC->SetLineColor(kBlack);
    h_EoverP_MC->Draw();
    h_EoverP_beforeCut->SetMarkerStyle(24);
    h_EoverP_beforeCut->SetMarkerColor(kBlue);
    h_EoverP_beforeCut->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_beforeCut->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_beforeCut->Draw("PEsame");
    TLegend *w11 = new TLegend(0.65,0.75,0.75,0.85);
	w11->AddEntry(h_EoverP_MC, " MC", "L");
	w11->AddEntry(h_EoverP_beforeCut, " RECO", "P");
    w11->SetBorderSize(0);   
    w11->SetFillStyle(0);
    w11->SetTextSize(0.05);
	w11->Draw("same");

    c12->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_EoverP_res_beforeCuts->GetXaxis()->SetRangeUser(0.99,1.03);
    h_EoverP_res_beforeCuts->GetYaxis()->SetRangeUser(-0.21,0.11);
    h_EoverP_res_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_res_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_res_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_res_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_res_beforeCuts->GetYaxis()->SetTitle("E_{EEMC}/|p|_{trk}");
    h_EoverP_res_beforeCuts->GetXaxis()->SetTitle("(E/|p|)_{MC}");
    h_EoverP_res_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_EoverP_res_beforeCuts->GetYaxis()->SetTitleOffset(2);
    h_EoverP_res_beforeCuts->Draw("colzsame");

    c12->cd(3);
    gPad->SetLogz(1);    
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    h_EoverP_response_beforeCuts->GetXaxis()->SetRangeUser(0.85,1.25);
    h_EoverP_response_beforeCuts->GetXaxis()->SetTitleOffset(1.2);
    h_EoverP_response_beforeCuts->GetYaxis()->SetTitleOffset(2);
    h_EoverP_response_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_response_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_response_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_response_beforeCuts->GetYaxis()->SetTitleSize(0.04);
	h_EoverP_response_beforeCuts->GetYaxis()->SetTitle("(E/|p|)_{MC}");
	h_EoverP_response_beforeCuts->GetXaxis()->SetTitle("E_{EEMC}/|p|_{trk}");
	h_EoverP_response_beforeCuts->SetLineColor(kBlack);
    h_EoverP_response_beforeCuts->Draw("colzsame");

    c12->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title1 = new TLatex();
    title1->SetNDC(); 
    title1->SetTextSize(0.05);
    title1->SetTextAlign(22);  
    title1->DrawLatex(0.5, 0.97, "E/|p| Truth vs. Reco");  

    c12->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_3 = new TLatex();
    title12_3->SetNDC(); 
    title12_3->SetTextSize(0.05);
    title12_3->SetTextAlign(22);  
    title12_3->DrawLatex(0.5, 0.97, "E/|p| Resolution");  

    c12->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title12_2 = new TLatex();
    title12_2->SetNDC(); 
    title12_2->SetTextSize(0.05);
    title12_2->SetTextAlign(22);  
    title12_2->DrawLatex(0.5, 0.97, "E/|p| Response");  

    c12->Print("./figures/analysisNote_plots/plot_EoverP.pdf");
}

void e_EoP_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_EoverP_MC_after = (TH1D*) file->Get("h_EoverP_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_EoverP_migration_beforeCuts = (TH2D*) file2->Get("h_EoverP_migration_beforeCuts"); 
    TH1D* h_EoverP_MC = (TH1D*) file2->Get("h_EoverP_MC");
    TH1D* h_EoverP_beforeCut = (TH1D*) file2->Get("h_EoverP_beforeCut");

    TCanvas* c26 = new TCanvas("c26","c26",1,1,1400,800);
    c26->Divide(3,2,0.01,0.01);
    c26->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    h_EoverP_migration_beforeCuts->GetXaxis()->SetTitleOffset(1.5);
    h_EoverP_migration_beforeCuts->GetYaxis()->SetTitleOffset(1.7);
    h_EoverP_migration_beforeCuts->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_migration_beforeCuts->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_migration_beforeCuts->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_migration_beforeCuts->GetXaxis()->SetTitle("(E/|p|)_{MC}");
    h_EoverP_migration_beforeCuts->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_migration_beforeCuts->GetXaxis()->SetRangeUser(0.98,1.025);
    h_EoverP_migration_beforeCuts->GetYaxis()->SetRangeUser(0.98,1.025);
    h_EoverP_migration_beforeCuts->Draw("colzsame");

    int nx_EoverP = h_EoverP_migration_beforeCuts->GetNbinsX();
    int ny_EoverP = h_EoverP_migration_beforeCuts->GetNbinsY();
    TH1D* h_EoverP_purity = new TH1D("h_EoverP_purity",";(E/|p|)_{reco bin};Purity",nx_EoverP,0,2);
    for (int j=1; j<=ny_EoverP; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx_EoverP; ++i) recoColSum += h_EoverP_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_EoverP_migration_beforeCuts->GetBinContent(j,j);
        if (recoColSum>0) h_EoverP_purity->SetBinContent(j, diag / recoColSum);
    }
    c26->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_EoverP_purity->GetXaxis()->SetTitleOffset(1.5);
    h_EoverP_purity->GetYaxis()->SetTitleOffset(1.5);
    h_EoverP_purity->GetXaxis()->SetRangeUser(0.98,1.03);
    h_EoverP_purity->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_purity->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_purity->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_purity->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_purity->GetXaxis()->SetTitle("(E/|p|)_{reco bin}");
    h_EoverP_purity->GetYaxis()->SetTitle("Purity");
    h_EoverP_purity->Draw("same");

    double binWidth_EoverPpurity = h_EoverP_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label26_2 = new TLatex(0.16,0.81, Form("bin width=%.2f",binWidth_EoverPpurity));
    label26_2->SetNDC();
    label26_2->SetTextSize(0.04);
    label26_2->Draw("same");

    TH1D* h_EoverP_stability = new TH1D("h_EoverP_stability",";(E/|p|)_{true bin};stability",nx_EoverP,0,2);
    for (int i=1; i<=nx_EoverP; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny_EoverP; ++j) trueRowSum += h_EoverP_migration_beforeCuts->GetBinContent(i,j);
        double diag = h_EoverP_migration_beforeCuts->GetBinContent(i,i);
        if (trueRowSum>0) h_EoverP_stability->SetBinContent(i, diag / trueRowSum);
    }
    c26->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_EoverP_stability->GetXaxis()->SetTitleOffset(1.5);
    h_EoverP_stability->GetYaxis()->SetTitleOffset(1.7);
    h_EoverP_stability->GetXaxis()->SetRangeUser(0.98,1.03);
    h_EoverP_stability->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_stability->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_stability->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_stability->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_stability->GetYaxis()->SetTitle("Stability");
    h_EoverP_stability->Draw("same");

    double binWidth_EoverPStability = h_EoverP_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label26_3 = new TLatex(0.16,0.81, Form("bin width=%.2f",binWidth_EoverPStability));
    label26_3->SetNDC();
    label26_3->SetTextSize(0.04);
    label26_3->Draw("same");
    
    int nb_EoverP = h_EoverP_MC->GetNbinsX();
    TH1D* h_EoverP_efficiency = new TH1D("h_EoverP_efficiency", ";(E/|p|);Efficiency", nb_EoverP, 0, 2);
    TH1D* h_EoverP_acceptance = new TH1D("h_EoverP_acceptance", ";(E/|p|);Acceptance", nb_EoverP, 0, 2);
    TH1D* h_EoverP_corrected = new TH1D("h_EoverP_corrected", ";(E/|p|);Acceptance Corrected", nb_EoverP, 0, 2);
    for (int i=1; i<=nb_EoverP; ++i) 
    {
        double Ntrue = h_EoverP_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb_EoverP; ++j) Naccepted += h_EoverP_migration_beforeCuts->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_EoverP_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_EoverP_beforeCut->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_EoverP_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_EoverP_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_EoverP_efficiency->SetBinContent(i, efficiency);
    }
    c26->cd(4);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_EoverP_acceptance->GetXaxis()->SetRangeUser(0.98,1.03);
    h_EoverP_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_EoverP_acceptance->Draw("same");

    c26->cd(5);
    gStyle->SetOptStat(0);
    h_EoverP_efficiency->GetXaxis()->SetRangeUser(0.98,1.03);
    h_EoverP_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_EoverP_efficiency->GetXaxis()->SetTitle("E/|p|");
    h_EoverP_efficiency->Draw("same");
    
    c26->cd(6);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);
    h_EoverP_corrected->GetXaxis()->SetRangeUser(0.98,1.03);
    h_EoverP_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_EoverP_corrected->GetYaxis()->SetLabelSize(0.04);
    h_EoverP_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_EoverP_corrected->GetYaxis()->SetTitleSize(0.04);
    h_EoverP_corrected->GetYaxis()->SetTitle("Corrected");
    h_EoverP_corrected->Draw("same");    

    c26->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_1 = new TLatex();
    title26_1->SetNDC(); 
    title26_1->SetTextSize(0.05);
    title26_1->SetTextAlign(22);  
    title26_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c26->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_3 = new TLatex();
    title26_3->SetNDC(); 
    title26_3->SetTextSize(0.05);
    title26_3->SetTextAlign(22);  
    title26_3->DrawLatex(0.5, 0.97, "Purity");  

    c26->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_4 = new TLatex();
    title26_4->SetNDC(); 
    title26_4->SetTextSize(0.05);
    title26_4->SetTextAlign(22);  
    title26_4->DrawLatex(0.5, 0.97, "Stability");

    c26->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_5 = new TLatex();
    title26_5->SetNDC(); 
    title26_5->SetTextSize(0.05);
    title26_5->SetTextAlign(22);  
    title26_5->DrawLatex(0.5, 0.97, "Acceptance"); 

    c26->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_2 = new TLatex();
    title26_2->SetNDC(); 
    title26_2->SetTextSize(0.05);
    title26_2->SetTextAlign(22);  
    title26_2->DrawLatex(0.5, 0.97, "Efficiency");  


    c26->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title26_6 = new TLatex();
    title26_6->SetNDC(); 
    title26_6->SetTextSize(0.05);
    title26_6->SetTextAlign(22);  
    title26_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");
    
    c26->Print("./figures/analysisNote_plots/plot_EoverP_migrationCorrectedPurityAcceptance.pdf");
}

void plot_tdist_allMethods_preTDR()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = (TH1D*)file->Get("h_t_MC");
	TH1D* h_t_REC_L = (TH1D*)file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_proj12 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*)file->Get("h_Nevents");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

    h_t_REC_L->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_L->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_L->SetBinError(i, newError);
    }
	h_t_REC_L->SetMarkerStyle(20); // method L RECO
	h_t_REC_L->SetMarkerColor(kP8Blue);
	h_t_REC_L->Draw("PEsame");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");

	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r142111 = new TLatex(0.18, 0.18, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r142111->SetNDC();
	r142111->SetTextSize(25);
	r142111->SetTextFont(43);
	r142111->SetTextColor(kBlack);
	r142111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.7,0.73,0.85);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC_L, "Sartre "+vm_label+" Method L RECO", "P");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #omega_{max} = #pi/12", "P");
	w7->Draw("same");

	c1->Print("./figures/analysisNote_plots/plot_t_dist_preTDR_proj.pdf");
}

void plot_tdist_manyAngles_preTDR()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = (TH1D*)file->Get("h_t_MC");
	TH1D* h_t_REC_proj12 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_REC_proj2 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi2");
    TH1D* h_t_REC_proj3 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi3");
    TH1D* h_t_REC_proj6 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi6");
    TH1D* h_t_REC_proj24 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi24");
    TH1D* h_phi_sartre_events = (TH1D*)file->Get("h_Nevents");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

	h_t_REC_proj2->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj2->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/2)));
    for (int i = 1; i <= h_t_REC_proj2->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj2->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj2->SetBinError(i, newError);
    }
	h_t_REC_proj2->SetMarkerStyle(42);
	h_t_REC_proj2->SetMarkerColor(kP10Violet);
	h_t_REC_proj2->SetLineColor(kP10Violet);
    h_t_REC_proj2->Draw("PEsame");

    h_t_REC_proj3->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj3->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/3)));
    for (int i = 1; i <= h_t_REC_proj3->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj3->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj3->SetBinError(i, newError);
    }
	h_t_REC_proj3->SetMarkerStyle(4);
	h_t_REC_proj3->SetMarkerColor(kP8Cyan);
	h_t_REC_proj3->SetLineColor(kP8Cyan);
	h_t_REC_proj3->Draw("PEsame");

    h_t_REC_proj6->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj6->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/6)));
    for (int i = 1; i <= h_t_REC_proj6->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj6->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj6->SetBinError(i, newError);
    }
	h_t_REC_proj6->SetMarkerStyle(31);
	h_t_REC_proj6->SetMarkerColor(kP8Gray);
	h_t_REC_proj6->SetLineColor(kP8Gray);
	h_t_REC_proj6->Draw("PEsame");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");

    h_t_REC_proj24->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj24->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/24)));
    for (int i = 1; i <= h_t_REC_proj24->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj24->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj24->SetBinError(i, newError);
    }
	h_t_REC_proj24->SetMarkerStyle(38);
	h_t_REC_proj24->SetMarkerColor(kP10Yellow);
	h_t_REC_proj24->SetLineColor(kP10Yellow);
	h_t_REC_proj24->Draw("PEsame");
	
	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r142111 = new TLatex(0.18, 0.18, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r142111->SetNDC();
	r142111->SetTextSize(25);
	r142111->SetTextFont(43);
	r142111->SetTextColor(kBlack);
	r142111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.6,0.73,0.88);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC_proj2, "Sartre "+vm_label+" RECO #omega_{max} = #pi/2", "P");
    w7->AddEntry(h_t_REC_proj3, "Sartre "+vm_label+" RECO #omega_{max} = #pi/3", "P");
    w7->AddEntry(h_t_REC_proj6, "Sartre "+vm_label+" RECO #omega_{max} = #pi/6", "P");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #omega_{max} = #pi/12", "P");
    w7->AddEntry(h_t_REC_proj24, "Sartre "+vm_label+" RECO #omega_{max} = #pi/24", "P");
	w7->Draw("same");

	c1->Print("./figures/analysisNote_plots/plot_t_dist_preTDR_ALL.pdf");
}

void t_QA()
{
    TFile* file = TFile::Open("z_coherentPhi_MCcuts.root","READ");
    TH1D* h_t_MC_after = (TH1D*) file->Get("h_t_MC_after");

    TFile* file2 = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH2D* h_t_MC = (TH2D*) file2->Get("h_t_MC");
    TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) file2->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_migration = (TH1D*) file2->Get("h_t_migration");

    TCanvas* c17 = new TCanvas("c17","c17",1,1,1400,800);
    c17->Divide(3,2,0.01,0.01);
    c17->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz(1);
    h_t_migration->GetYaxis()->SetTitleOffset(1.8);
    h_t_migration->GetXaxis()->SetTitleOffset(1.5);
    h_t_migration->GetXaxis()->SetLabelSize(0.04);  
    h_t_migration->GetYaxis()->SetLabelSize(0.04);
    h_t_migration->GetXaxis()->SetTitleSize(0.04);  
    h_t_migration->GetYaxis()->SetTitleSize(0.04);
    h_t_migration->GetXaxis()->SetNdivisions(505);
    h_t_migration->Draw("colzsame");

    int nx = h_t_migration->GetNbinsX();
    int ny = h_t_migration->GetNbinsY();
    TH1D* h_t_purity = new TH1D("h_t_purity",";|t|_{reco bin} [GeV/c];Purity",nx,0,0.2);
    for (int j=1; j<=ny; ++j) 
    {
        double recoColSum = 0;
        for (int i=1; i<=nx; ++i) recoColSum += h_t_migration->GetBinContent(i,j);
        double diag = h_t_migration->GetBinContent(j,j);
        if (recoColSum>0) h_t_purity->SetBinContent(j, diag / recoColSum);
    }
    c17->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_t_purity->GetYaxis()->SetTitleOffset(1.8);
    h_t_purity->GetXaxis()->SetTitleOffset(1.5);
    h_t_purity->GetXaxis()->SetLabelSize(0.04);  
    h_t_purity->GetYaxis()->SetLabelSize(0.04);
    h_t_purity->GetXaxis()->SetTitleSize(0.04);  
    h_t_purity->GetYaxis()->SetTitleSize(0.04);
    h_t_purity->GetXaxis()->SetNdivisions(505);
    h_t_purity->GetXaxis()->SetTitle("|t|_{reco bin} [GeV/c]^{2}");
    h_t_purity->GetYaxis()->SetTitle("Purity");
    h_t_purity->Draw("same");

    double binWidth_tpurtiy = h_t_purity->GetXaxis()->GetBinWidth(1);
    TLatex* label17_2 = new TLatex(0.23,0.81,"bin width");
    label17_2->SetNDC();
    label17_2->SetTextSize(0.04);
    label17_2->Draw("same");
    TLatex* label17_22 = new TLatex(0.23,0.75, Form("=%.3f GeV^{2}/c^{2}",binWidth_tpurtiy));
    label17_22->SetNDC();
    label17_22->SetTextSize(0.04);
    label17_22->Draw("same");

    TH1D* h_t_stability = new TH1D("h_t_stability",";|t|_{true bin} [GeV/c]^{2};Stability",nx,0,0.2);
    for (int i=1; i<=nx; ++i) 
    {
        double trueRowSum = 0;
        for (int j=1; j<=ny; ++j) trueRowSum += h_t_migration->GetBinContent(i,j);
        double diag = h_t_migration->GetBinContent(i,i);
        if (trueRowSum>0) h_t_stability->SetBinContent(i, diag / trueRowSum);
    }
    c17->cd(3);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_t_stability->GetYaxis()->SetTitleOffset(1.8);
    h_t_stability->GetXaxis()->SetTitleOffset(1.5);
    h_t_stability->GetXaxis()->SetLabelSize(0.04);  
    h_t_stability->GetYaxis()->SetLabelSize(0.04);
    h_t_stability->GetXaxis()->SetTitleSize(0.04);  
    h_t_stability->GetYaxis()->SetTitleSize(0.04);
    h_t_stability->GetXaxis()->SetNdivisions(505);
    h_t_stability->Draw("same");

    double binWidth_tstability = h_t_stability->GetXaxis()->GetBinWidth(1);
    TLatex* label17_3 = new TLatex(0.23,0.81, "bin width=");
    label17_3->SetNDC();
    label17_3->SetTextSize(0.04);
    label17_3->Draw("same");
    TLatex* label17_32 = new TLatex(0.23,0.75, Form("=%.3f GeV^{2}/c^{2}",binWidth_tstability));
    label17_32->SetNDC();
    label17_32->SetTextSize(0.04);
    label17_32->Draw("same");

    int nb = h_t_MC->GetNbinsX();
    TH1D* h_t_efficiency = new TH1D("h_t_efficiency", ";|t| [GeV/c]^{2};Efficiency", nb, 0, 0.2);
    TH1D* h_t_acceptance = new TH1D("h_t_acceptance",";|t| [GeV/c]^{2};Acceptance",nb,0,0.2);
    TH1D* h_t_corrected = new TH1D("h_t_corrected",";|t| [GeV/c]^{2};Acceptance Corrected",nb,0,0.2);
    for (int i=1; i<=nb; ++i) 
    {
        double Ntrue = h_t_MC->GetBinContent(i);   
        double Naccepted = 0;
        for (int j=1; j<=nb; ++j) Naccepted += h_t_migration->GetBinContent(i,j);
        double acceptance = (Ntrue>0) ? Naccepted / Ntrue : 0;
        h_t_acceptance->SetBinContent(i, acceptance);
        double Nreco_in_bin_i = h_t_REC_wRES_cut_pi12->GetBinContent(i);
        double corrected = (acceptance>0) ? Nreco_in_bin_i / acceptance : 0;
        h_t_corrected->SetBinContent(i, corrected);
        double Nreco_from_bin_i = h_t_MC_after->GetBinContent(i);
        double efficiency = (Ntrue > 0) ? Nreco_from_bin_i / Ntrue : 0;
        h_t_efficiency->SetBinContent(i, efficiency);
    }
    c17->cd(4);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_t_acceptance->GetYaxis()->SetTitleOffset(1.8);
    h_t_acceptance->GetXaxis()->SetTitleOffset(1.5);
    h_t_acceptance->GetXaxis()->SetLabelSize(0.04);  
    h_t_acceptance->GetYaxis()->SetLabelSize(0.04);
    h_t_acceptance->GetXaxis()->SetTitleSize(0.04);  
    h_t_acceptance->GetYaxis()->SetTitleSize(0.04);
    h_t_acceptance->GetXaxis()->SetNdivisions(505);
    h_t_acceptance->GetYaxis()->SetTitle("Acceptance");
    h_t_acceptance->Draw("same");

    c17->cd(5);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_t_efficiency->GetYaxis()->SetTitleOffset(1.8);
    h_t_efficiency->GetXaxis()->SetTitleOffset(1.5);
    h_t_efficiency->GetXaxis()->SetLabelSize(0.04);  
    h_t_efficiency->GetYaxis()->SetLabelSize(0.04);
    h_t_efficiency->GetXaxis()->SetTitleSize(0.04);  
    h_t_efficiency->GetYaxis()->SetTitleSize(0.04);
    h_t_efficiency->GetXaxis()->SetNdivisions(505);
    h_t_efficiency->GetYaxis()->SetTitle("Efficiency");
    h_t_efficiency->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_efficiency->Draw("same");

    c17->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    h_t_corrected->GetYaxis()->SetTitleOffset(1.8);
    h_t_corrected->GetXaxis()->SetTitleOffset(1.5);
    h_t_corrected->GetXaxis()->SetLabelSize(0.04);  
    h_t_corrected->GetYaxis()->SetLabelSize(0.04);
    h_t_corrected->GetXaxis()->SetTitleSize(0.04);  
    h_t_corrected->GetYaxis()->SetTitleSize(0.04);
    h_t_corrected->GetXaxis()->SetNdivisions(505);
    h_t_corrected->GetYaxis()->SetTitle("Corrected");
    h_t_corrected->Draw("same");

    c17->cd(1);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_1 = new TLatex();
    title17_1->SetNDC(); 
    title17_1->SetTextSize(0.05);
    title17_1->SetTextAlign(22);  
    title17_1->DrawLatex(0.5, 0.97, "Bin Migration");  

    c17->cd(2);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_3 = new TLatex();
    title17_3->SetNDC(); 
    title17_3->SetTextSize(0.05);
    title17_3->SetTextAlign(22);  
    title17_3->DrawLatex(0.5, 0.97, "Purity");  

    c17->cd(3);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_4 = new TLatex();
    title17_4->SetNDC(); 
    title17_4->SetTextSize(0.05);
    title17_4->SetTextAlign(22);  
    title17_4->DrawLatex(0.5, 0.97, "Stability");

    c17->cd(4);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_5 = new TLatex();
    title17_5->SetNDC(); 
    title17_5->SetTextSize(0.05);
    title17_5->SetTextAlign(22);  
    title17_5->DrawLatex(0.5, 0.97, "Acceptance");

    c17->cd(5);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_2 = new TLatex();
    title17_2->SetNDC(); 
    title17_2->SetTextSize(0.05);
    title17_2->SetTextAlign(22);  
    title17_2->DrawLatex(0.5, 0.97, "Efficiency");  

    c17->cd(6);  
    gPad->SetTopMargin(0.08);  
    TLatex* title17_6 = new TLatex();
    title17_6->SetNDC(); 
    title17_6->SetTextSize(0.05);
    title17_6->SetTextAlign(22);  
    title17_6->DrawLatex(0.5, 0.97, "Acceptance Corrected");

    c17->Print("./figures/analysisNote_plots/plot_t_migrationCorrectedPurityAcceptance.pdf");
}

void plot_t_resolution()
{
	TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");	
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	TH2D* h_t_res = (TH2D*) file->Get("h_t_res_proj_percent_pi12");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,800);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.13);
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.13);

	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", " #delta t/|t| (resolution) ", 100,0,0.2,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-3, 1000);
	base1->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base1,1.2,1.6);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.8);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.7);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();
	TH2D* h_res = (TH2D*) h_t_res;

	TH1D* h_res_1D = new TH1D("h_res_1D","",100,0,0.2);
	for(int ibin=0;ibin<h_res->GetNbinsX();ibin++)
	{
		TH1D* tmp=h_res->ProjectionY("tmp",ibin+1,ibin+1);
		double sigma = tmp->GetStdDev();
		double sigmaerror = tmp->GetStdDevError();
		h_res_1D->SetBinContent(ibin+1, sigma);
		h_res_1D->SetBinError(ibin+1, sigmaerror);
	}

	h_res_1D->SetMarkerSize(1.6);
	h_res_1D->SetMarkerColor(kBlack);
	h_res_1D->SetLineColor(kBlack);
	h_res_1D->SetMarkerStyle(20);

	h_res_1D->Fit("pol0","RMS0","",0.011,0.019);//first  dip
	h_res_1D->Fit("pol0","RMS0","",0.045,0.053);//second dip
	h_res_1D->Fit("pol0","RMS0","",0.095,0.102);//third  dip
    h_res_1D->Fit("pol0","RMS0","",0.167,0.175);//fourth  dip
	h_res_1D->Draw("Psame");

	TLatex* r42 = new TLatex(0.15, 0.91, "eAu 10x100 GeV");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.8,0.91, "ePIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44_2 = new TLatex(0.17, 0.18, vm_label+" #rightarrow "+daug_label );
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TPad* drawPad = new TPad("pad_etalab_11","pad_etalab_11",0.16,0.53,0.47,0.83);
	drawPad->SetLeftMargin(0.08);
	drawPad->SetRightMargin(0.08);
	drawPad->SetTopMargin(0.0);
	drawPad->SetBottomMargin(0.08);
	drawPad->Draw("same");
	drawPad->SetTicks();
	drawPad->SetLogz(1);
	drawPad->cd();
	TH1D* base2 = makeHist("base2", "", "MC", " resolution ", 100,0,0.2,kBlack);
	base2->GetYaxis()->SetRangeUser(-10, 1.5);
	base2->GetXaxis()->SetTitleColor(kBlack);
	base2->GetXaxis()->SetLabelColor(kBlack);
	base2->GetYaxis()->SetLabelColor(kBlack);
	base2->GetXaxis()->SetTitle("Resolution");
	base2->GetYaxis()->SetTitle("MC");
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base2,3,3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->Draw();

	h_res->Draw("colzsame");
	gPad->Update(); 
    TPaletteAxis *palette = (TPaletteAxis*)h_res->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9); 
    palette->SetX2NDC(0.95); 
    palette->SetY1NDC(0.08); 
    palette->SetY2NDC(0.99); 
    gPad->Modified();
    gPad->Update();

	c2->Print("./figures/analysisNote_plots/plot_t_resolution.pdf");
}

void plot_background_preTDR()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_phi_sartre_events = (TH1F*) phi_t_file->Get("h_Nevents");

	TFile* sartre_rho = TFile::Open("z_coherentRho_noPID_allVeotesCuts.root","READ");
    TH1D* h_t_rho = (TH1D*) sartre_rho->Get("h_t_MC");
    TH1F* h_rho_events = (TH1F*) sartre_rho->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("z_incoherentPhi_allVeotesCuts.root","READ");
    TH1D* h_incoherent_phi_noVetoes = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1F* h_incoherent_events = (TH1F*) incoherent_phi_noVetoes->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("z_DIS_allVeotesCuts.root","READ");
    TH1D* h_DIS_noVetoes = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1F* h_DIS_events = (TH1F*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double sigma_rho = 5418.68; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_rho = h_rho_events->GetEntries();
	double nEvents_DIS = h_DIS_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
 	h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->GetYaxis()->SetRangeUser(1e-2,1e7);
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
    
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("Psame");

    h_t_rho->Rebin(rebin_width);
    h_t_rho->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_rho->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    double simu_rho_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_rho_preTDR = preTDR_lumi/simu_rho_lumi;
    for (int i = 1; i <= h_t_rho->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho->GetBinError(i);
        double newError = binError / sqrt(ratio_rho_preTDR);
        h_t_rho->SetBinError(i, newError);
    }
    h_t_rho->SetMarkerStyle(28);
    h_t_rho->SetMarkerColor(kP8Gray);
    h_t_rho->SetLineColor(kP8Gray);
    h_t_rho->Draw("PEsame");

	h_DIS_noVetoes->Rebin(rebin_width);
    h_DIS_noVetoes->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width);
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes->SetBinError(i, newError);
    }
    h_DIS_noVetoes->SetMarkerStyle(2);
    h_DIS_noVetoes->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes->SetLineColor(kP8Cyan);
    h_DIS_noVetoes->Draw("PEsame");

    h_incoherent_phi_noVetoes->Rebin(rebin_width);
    h_incoherent_phi_noVetoes->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width);
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes->SetMarkerStyle(21);
    h_incoherent_phi_noVetoes->SetMarkerColor(kP8Orange);
    h_incoherent_phi_noVetoes->SetLineColor(kP8Orange);
    h_incoherent_phi_noVetoes->Draw("PEsame");
	
    TLegend *w14_213 = new TLegend(0.6,0.63,0.72,0.88);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #omega_{max}= #pi/12", "P");
	w14_213->AddEntry(h_DIS_noVetoes,"DIS", "P");
    w14_213->AddEntry(h_incoherent_phi_noVetoes,"Incoherent #phi: #omega_{max}= #pi/12", "P");    
    w14_213->AddEntry(h_t_rho,"Coherent #rho: #omega_{max}= #pi/12", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	
	
    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.16, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.21, 0.85, " Simulation 25.10.2/3, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.16, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.16, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_t_backgroundv5_preTDR.pdf");
}

void plot_DIS_noVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* DIS_noVetoes = TFile::Open("z_DIS_allVeotesCuts.root","READ");
    TH1D* h_DIS_noVetoes = (TH1D*) DIS_noVetoes->Get("h_t_MC");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_DIS = h_DIS_events->GetEntries();
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_DIS_noVetoes->Rebin(rebin_width);
    h_DIS_noVetoes->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes->SetBinError(i, newError);
    }
    h_DIS_noVetoes->SetMarkerStyle(2);
    h_DIS_noVetoes->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes->SetLineColor(kP8Cyan);
    h_DIS_noVetoes->Draw("PEsame");
        
	TLegend *w14_213 = new TLegend(0.6,0.74,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_DIS_noVetoes,"DIS", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_t_DIS_noVetoes_preTDR_proj.pdf");
}

void plot_DIS_together_preTDR()
{
    TFile* coherent_phi_allVetoes = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_MC = (TH1D*) coherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) coherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) coherent_phi_allVetoes->Get("h_Nevents");

    TFile* DIS_file = TFile::Open("z_DIS_allVeotesCuts.root","READ");
    TH1D* h_t_DIS_MC = (TH1D*) DIS_file->Get("h_t_MC");
    TH1D* h_t_pi12_afterHFS = (TH1D*) DIS_file->Get("h_t_pi12_afterHFS");
    TH1D* h_t_pi12_afterKaons = (TH1D*) DIS_file->Get("h_t_pi12_afterKaons");
    TH1D* h_t_pi12_afterMass = (TH1D*) DIS_file->Get("h_t_pi12_afterMass");
    TH1D* h_t_pi12_afterOMD = (TH1D*) DIS_file->Get("h_t_pi12_afterOMD");
    TH1D* h_events = (TH1D*) DIS_file->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_DIS = 50184.8; // nb
    double nEvents_DIS = h_events->GetEntries();
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    // fit wedge cut distribution to DIS truth distribution:
    TH1D* h_raw_DIS = (TH1D*) h_t_DIS_MC->Clone("h_raw_DIS");
    h_raw_DIS->Reset();
    h_raw_DIS->Add(h_t_DIS_MC);

    TF1* fDIS = new TF1("fDIS", "[0]*exp(-[1]*x)", 0, 0.2);
    h_raw_DIS->Fit(fDIS, "L");

    // expected number of events
    double lumi_10fb = 1e7; // nb^-1
    double DIS_10fb_to_nb = (sigma_DIS/nEvents_DIS) * lumi_10fb;
    double N_pi12_raw = h_t_REC_wRES_cut_pi12->Integral(); // 38 events total
    double Nexp_pi12_10fb = N_pi12_raw * DIS_10fb_to_nb;

    // sample the fit wedge cut curve
    TH1D* h_t_pi12_pseudo = (TH1D*) h_t_REC_wRES_cut_pi12->Clone("h_t_pi12_pseudo");
    h_t_pi12_pseudo->Reset();
    h_t_pi12_pseudo->Sumw2();

    double tmin = h_t_REC_wRES_cut_pi12->GetXaxis()->GetXmin();
    double tmax = h_t_REC_wRES_cut_pi12->GetXaxis()->GetXmax();
    int Ngen = (int) std::round(Nexp_pi12_10fb);
    cout << Ngen << endl;
    for (int i = 0; i < Ngen; ++i) 
    {
        double t = fDIS->GetRandom(tmin, tmax);
        h_t_pi12_pseudo->Fill(t);
    }

    // put in terms of cross section
    double sigma_pi12 = Nexp_pi12_10fb / lumi_10fb;  // nb
    double scale = sigma_pi12 / h_t_pi12_pseudo->Integral("width");
    h_t_pi12_pseudo->Scale(scale);

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1200,800);
    c1->Divide(1,1,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
    h_t_MC->SetStats(0);
    h_t_MC->GetXaxis()->SetTitle("");
    h_t_MC->SetTitle("DIS Background");
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
    h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherentMC_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherentMC_phi_preTDR = preTDR_lumi/simu_coherentMC_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherentMC_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineStyle(1);
    h_t_MC->SetLineColor(kBlack);
    h_t_MC->Draw("HIST");

    h_t_DIS_MC->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_DIS_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_DIS_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_DIS_MC->SetBinError(i, newError);
    }
    h_t_DIS_MC->SetMarkerStyle(2);
    h_t_DIS_MC->SetMarkerColor(kP10Blue);
    h_t_DIS_MC->SetLineColor(kP10Blue);
    h_t_DIS_MC->Draw("PE same");

    h_t_pi12_afterHFS->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterHFS->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterHFS->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterHFS->SetBinError(i, newError);
    }
    h_t_pi12_afterHFS->SetMarkerStyle(4);
    h_t_pi12_afterHFS->SetMarkerColor(kP10Yellow);
    h_t_pi12_afterHFS->SetLineColor(kP10Yellow);
    h_t_pi12_afterHFS->Draw("PE same");

    h_t_pi12_afterKaons->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterKaons->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterKaons->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterKaons->SetBinError(i, newError);
    }
    h_t_pi12_afterKaons->SetMarkerStyle(22);
    h_t_pi12_afterKaons->SetMarkerColor(kP8Pink);
    h_t_pi12_afterKaons->SetLineColor(kP8Pink);
    h_t_pi12_afterKaons->Draw("PE same");

    h_t_pi12_afterMass->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterMass->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterMass->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterMass->SetBinError(i, newError);
    }
    h_t_pi12_afterMass->SetMarkerStyle(31);
    h_t_pi12_afterMass->SetMarkerColor(kP10Violet);
    h_t_pi12_afterMass->SetLineColor(kP10Violet);
    h_t_pi12_afterMass->Draw("PE same");

    h_t_pi12_afterOMD->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterOMD->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterOMD->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterOMD->SetBinError(i, newError);
    }

    h_t_pi12_afterOMD->SetMarkerStyle(48);
    h_t_pi12_afterOMD->SetMarkerColor(kP10Brown);
    h_t_pi12_afterOMD->SetLineColor(kP10Brown);
    h_t_pi12_afterOMD->Draw("PE same");

    h_t_pi12_pseudo->SetLineColor(kP10Cyan);
    h_t_pi12_pseudo->SetLineStyle(3);
    h_t_pi12_pseudo->Draw("PE same"); // Fit wedge cut distribution

    TLegend *w14_213 = new TLegend(0.58,0.58,0.89,0.89);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");//Form("MC: %.1f",h_t_MC->GetEntries()),"L");
    w14_213->AddEntry(h_t_DIS_MC,"DIS", "P");//Form("DIS truth: %.1f",h_t_DIS_MC->GetEntries()),"P");
    w14_213->AddEntry(h_t_pi12_afterHFS,"After HFS", "P");//Form("After HFS==2: %.1f",h_t_pi12_afterHFS->GetEntries()),"P");
    w14_213->AddEntry(h_t_pi12_afterKaons,"After kaons", "P");//Form("After kaons==2: %.1f",h_t_pi12_afterKaons->GetEntries()),"P");
    w14_213->AddEntry(h_t_pi12_afterMass,"After mass", "P");//Form("After mass: %.1f",h_t_pi12_afterMass->GetEntries()),"P");
    w14_213->AddEntry(h_t_pi12_afterOMD,"After RP+OMD", "P");//Form("After RP+OMD: %.1f",h_t_pi12_afterOMD->GetEntries()),"P");
    w14_213->AddEntry(h_t_pi12_pseudo,"Fit #omega_{max}=#pi/12 to MC", "L");//Form("Fit #omega_{max}=#pi/12 to MC: %.1f",h_t_pi12_pseudo->GetEntries()),"L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");
    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_t_DIS_compareALL_preTDR_proj.pdf");
}

void plot_DIS_allVetoes()
{
    TFile* coherent_phi_allVetoes = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) coherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) coherent_phi_allVetoes->Get("h_Nevents");

    TFile* DIS_file = TFile::Open("z_DIS_allVeotesCuts.root","READ");
    TH1D* h_t_DIS_MC = (TH1D*) DIS_file->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12_DIS = (TH1D*) DIS_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_events = (TH1D*) DIS_file->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_DIS = 50184.8; // nb
    double nEvents_DIS = h_events->GetEntries();
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    // fit wedge cut distribution to DIS truth distribution:
    TH1D* h_raw_DIS = (TH1D*) h_t_DIS_MC->Clone("h_raw_DIS");
    h_raw_DIS->Reset();
    h_raw_DIS->Add(h_t_DIS_MC);

    TF1* fDIS = new TF1("fDIS", "[0]*exp(-[1]*x)", 0, 0.2);
    h_raw_DIS->Fit(fDIS, "L");

    // expected number of events
    double lumi_10fb = 1e7; // nb^-1
    double DIS_10fb_to_nb = (sigma_DIS/nEvents_DIS) * lumi_10fb;
    double N_pi12_raw = h_t_REC_wRES_cut_pi12_DIS->Integral(); // 38 events total
    double Nexp_pi12_10fb = N_pi12_raw * DIS_10fb_to_nb;

    // sample the fit wedge cut curve
    TH1D* h_t_pi12_pseudo = (TH1D*) h_t_REC_wRES_cut_pi12_DIS->Clone("h_t_pi12_pseudo");
    h_t_pi12_pseudo->Reset();
    h_t_pi12_pseudo->Sumw2();

    double tmin = h_t_REC_wRES_cut_pi12_DIS->GetXaxis()->GetXmin();
    double tmax = h_t_REC_wRES_cut_pi12_DIS->GetXaxis()->GetXmax();
    int Ngen = (int) std::round(Nexp_pi12_10fb);
    cout << Ngen << endl;
    for (int i = 0; i < Ngen; ++i) 
    {
        double t = fDIS->GetRandom(tmin, tmax);
        h_t_pi12_pseudo->Fill(t);
    }

    // put in terms of cross section
    double sigma_pi12 = Nexp_pi12_10fb / lumi_10fb;  // nb
    double scale = sigma_pi12 / h_t_pi12_pseudo->Integral("width");
    h_t_pi12_pseudo->Scale(scale);

    TCanvas* c1 = new TCanvas("c1","c1",1,1,1200,800);
    c1->Divide(1,1,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(1);
    h_t_REC_wRES_cut_pi12->SetStats(0);
    h_t_REC_wRES_cut_pi12->GetXaxis()->SetTitle("");
    h_t_REC_wRES_cut_pi12->SetTitle("DIS Background");
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
    h_t_REC_wRES_cut_pi12->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    double simu_coherentMC_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherentMC_phi_preTDR = preTDR_lumi/simu_coherentMC_lumi;
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherentMC_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
    h_t_REC_wRES_cut_pi12->SetLineStyle(1);
    h_t_REC_wRES_cut_pi12->SetLineColor(kBlack);
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetRangeUser(1e-2,1e5);
    h_t_REC_wRES_cut_pi12->Draw("HIST");

    h_t_pi12_pseudo->SetLineColor(kP10Cyan);
    h_t_pi12_pseudo->SetLineStyle(3);
    h_t_pi12_pseudo->Draw("PE same"); // Fit wedge cut distribution

    TLegend *w14_213 = new TLegend(0.58,0.58,0.89,0.89);
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #omega_{max}=#pi/12", "L");
    w14_213->AddEntry(h_t_pi12_pseudo,"Fit #omega_{max}=#pi/12 to MC", "L");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");
    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_t_DIS_preTDR_proj.pdf");
}

void plot_rho_mass_noPID()
{
    TFile* sartre_rho = TFile::Open("z_coherentRho_noPID_allVeotesCuts.root","READ");
    TH1D* h_t_rho_MC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_MC");
    TH1D* h_t_rho_REC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_REC");
    TH1D* h_t_rho_REC_mass_beforeCut = (TH1D*) sartre_rho->Get("h_VM_mass_REC_before");
    
    TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_phi_MC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_MC");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_phi_MC_mass->GetYaxis()->SetRangeUser(1e-1,1e11);
    h_t_phi_MC_mass->GetYaxis()->SetTitle("events");
    h_t_phi_MC_mass->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
    h_t_phi_MC_mass->SetTitle("Invariant Mass Distribution (no PID)");
    h_t_phi_MC_mass->SetLineColor(kBlack);
    h_t_phi_MC_mass->Draw();
    h_t_rho_MC_mass->SetLineColor(kP8Blue);
    h_t_rho_MC_mass->SetLineStyle(7);
    h_t_rho_MC_mass->Draw("same");
    h_t_rho_REC_mass->SetLineColor(kRed);
    h_t_rho_REC_mass->SetLineStyle(3);
    h_t_rho_REC_mass->Draw("same");
    h_t_rho_REC_mass_beforeCut->SetLineColor(kP8Green);
    h_t_rho_REC_mass_beforeCut->SetLineStyle(9);
    h_t_rho_REC_mass_beforeCut->Draw("same");

    TLegend *w14_213 = new TLegend(0.45,0.62,0.75,0.88);
    w14_213->AddEntry(h_t_rho_MC_mass,"Coherent #rho: MC", "L");//Form("#rho MC: %.1f",h_t_rho_MC_mass->GetEntries()),"L");
    w14_213->AddEntry(h_t_rho_REC_mass_beforeCut,"Coherent #rho: REC before mass cut", "L");//Form("#rho reco before mass cut: %.1f",h_t_rho_REC_mass_beforeCut->GetEntries()),"L"); 
    w14_213->AddEntry(h_t_rho_REC_mass,"Coherent #rho: REC after mass cut", "L");//Form("#rho reco after mass cut: %.1f",h_t_rho_REC_mass->GetEntries()),"L");
    w14_213->AddEntry(h_t_phi_MC_mass,"Coherent #phi: MC", "L"); //Form("#phi MC: %.1f",h_t_phi_MC_mass->GetEntries()),"L");   
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_VM_mass_noPID.pdf");
}

void plot_rho_mass_wPID()
{
    TFile* sartre_rho = TFile::Open("z_coherentRho_PID_allVeotesCuts.root","READ");
    TH1D* h_t_rho_MC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_MC");
    TH1D* h_t_rho_REC_mass = (TH1D*) sartre_rho->Get("h_VM_mass_REC");
    TH1D* h_t_rho_REC_mass_beforeCut = (TH1D*) sartre_rho->Get("h_VM_mass_REC_before");
    
    TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_phi_MC_mass = (TH1D*) phi_t_file->Get("h_VM_mass_MC");

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    c14_213->Divide(1,1,0.01,0.01);
    c14_213->cd(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_phi_MC_mass->GetYaxis()->SetRangeUser(1e-1,1e11);
    h_t_phi_MC_mass->GetYaxis()->SetTitle("events");
    h_t_phi_MC_mass->GetXaxis()->SetTitle("VM mass [GeV/c^{2}]");
    h_t_phi_MC_mass->SetTitle("Invariant Mass Distribution (PID)");
    h_t_phi_MC_mass->SetLineColor(kBlack);
    h_t_phi_MC_mass->Draw();
    h_t_rho_MC_mass->SetLineColor(kP8Blue);
    h_t_rho_MC_mass->SetLineStyle(7);
    h_t_rho_MC_mass->Draw("same");
    h_t_rho_REC_mass->SetLineColor(kRed);
    h_t_rho_REC_mass->SetLineStyle(3);
    h_t_rho_REC_mass->Draw("same");
    h_t_rho_REC_mass_beforeCut->SetLineColor(kP8Green);
    h_t_rho_REC_mass_beforeCut->SetLineStyle(9);
    h_t_rho_REC_mass_beforeCut->Draw("same");

    TLegend *w14_213 = new TLegend(0.45,0.62,0.75,0.88);
    w14_213->AddEntry(h_t_rho_MC_mass,"Coherent #rho: MC", "L");//Form("#rho MC: %.1f",h_t_rho_MC_mass->GetEntries()),"L");
    w14_213->AddEntry(h_t_rho_REC_mass_beforeCut,"Coherent #rho: REC before mass cut", "L");//Form("#rho reco before mass cut: %.1f",h_t_rho_REC_mass_beforeCut->GetEntries()),"L");    
    w14_213->AddEntry(h_t_rho_REC_mass,"Coherent #rho: REC after mass cut", "L");//Form("#rho reco after mass cut: %.1f",h_t_rho_REC_mass->GetEntries()),"L");
    w14_213->AddEntry(h_t_phi_MC_mass,"Coherent #phi: MC", "L"); //Form("#phi MC: %.1f",h_t_phi_MC_mass->GetEntries()),"L");   
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_VM_mass_wPID.pdf");
}

void plot_t_rho_preTDR()
{
    TFile* coherent_phi_allVetoes = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_MC = (TH1D*) coherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) coherent_phi_allVetoes->Get("h_Nevents");

    TFile* coherent_rho_allVetoes = TFile::Open("z_coherentRho_PID_allVeotesCuts.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) coherent_rho_allVetoes->Get("h_t_MC");
    TH1D* h_rho_sartre_events = (TH1D*) coherent_rho_allVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_rho = 5418.68; // nb
    double nEvents_rho = h_rho_sartre_events->GetEntries();
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   
	
    TCanvas* c1 = new TCanvas("c1","c1",1,1,1200,800);
    c1->Divide(1,1,0.01,0.01);
    c1->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->SetTitle("#rho Distribution");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherentMC_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherentMC_phi_preTDR = preTDR_lumi/simu_coherentMC_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherentMC_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineStyle(1);
    h_t_MC->SetLineColor(kBlack);
    h_t_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_MC->Draw("HISTsame");

    h_t_rho_MC->Rebin(rebin_width);
    h_t_rho_MC->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    double simu_coherent_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_rho_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_rho_MC->SetBinError(i, newError);
    }
    h_t_rho_MC->SetMarkerStyle(28);
    h_t_rho_MC->SetMarkerColor(kP8Gray);
    h_t_rho_MC->SetLineColor(kP8Gray);
    h_t_rho_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_rho_MC->Draw("PE same");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.7, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.7, " Simulation 25.10.2/3");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.13, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.13, 0.85, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r4 = new TLatex(0.13, 0.8, "10x100 GeV");
	r4->SetNDC();
	r4->SetTextSize(30);
	r4->SetTextFont(43);
	r4->SetTextColor(kBlack);
	r4->Draw("same"); 

    TLatex* r42 = new TLatex(0.15, 0.19, "|t|_{#omega_{max}} normalized by #pi/2/#omega_{max}");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	//r42->Draw("same");

    TLatex* r421111 = new TLatex(0.15, 0.14, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    TLegend *w14_213 = new TLegend(0.5,0.7,0.9,0.89);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_rho_MC,"Coherent #rho: MC", "P");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_t_rho_preTDR_projv2.pdf");
}

void plot_coherent_rho_wPID()
{
    TFile* coherent_phi_allVetoes = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_MC = (TH1D*) coherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) coherent_phi_allVetoes->Get("h_Nevents");

    TFile* coherent_rho_allVetoes = TFile::Open("z_coherentRho_PID_allVeotesCuts.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) coherent_rho_allVetoes->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) coherent_rho_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_after_VMrec = (TH1D*) coherent_rho_allVetoes->Get("h_t_after_VMrec");
    TH1D* h_t_pi12_afterHFS = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterHFS");
    TH1D* h_t_pi12_afterKaons = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterKaons");
    TH1D* h_t_pi12_afterMass = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterMass");
    TH1D* h_t_pi12_afterOMD = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterOMD");
    TH1D* h_rho_sartre_events = (TH1D*) coherent_rho_allVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_rho = 5418.68; // nb
    double nEvents_rho = h_rho_sartre_events->GetEntries();
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   
	
    TCanvas* c1 = new TCanvas("c1","c1",1,1,1200,800);
    c1->Divide(1,1,0.01,0.01);
    c1->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->SetTitle("#rho Distribution with PID");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherentMC_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherentMC_phi_preTDR = preTDR_lumi/simu_coherentMC_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherentMC_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineStyle(1);
    h_t_MC->SetLineColor(kBlack);
    h_t_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_MC->Draw("HISTsame");

    h_t_rho_MC->Rebin(rebin_width);
    h_t_rho_MC->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    double simu_coherent_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_rho_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_rho_MC->SetBinError(i, newError);
    }
    h_t_rho_MC->SetMarkerStyle(28);
    h_t_rho_MC->SetMarkerColor(kP8Gray);
    h_t_rho_MC->SetLineColor(kP8Gray);
    h_t_rho_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_rho_MC->Draw("PE same");

    h_t_after_VMrec->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_after_VMrec->GetNbinsX(); ++i) 
    {
        double binError = h_t_after_VMrec->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_after_VMrec->SetBinError(i, newError);
    }
    h_t_after_VMrec->SetMarkerStyle(25);
    h_t_after_VMrec->SetMarkerColor(kP10Violet);
    h_t_after_VMrec->SetLineColor(kP10Violet);
    h_t_after_VMrec->Draw("PE same");

    h_t_pi12_afterHFS->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterHFS->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterHFS->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterHFS->SetBinError(i, newError);
    }
    h_t_pi12_afterHFS->SetMarkerStyle(4);
    h_t_pi12_afterHFS->SetMarkerColor(kP10Orange);
    h_t_pi12_afterHFS->SetLineColor(kP10Orange);
    h_t_pi12_afterHFS->Draw("PE same");

    h_t_pi12_afterKaons->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterKaons->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterKaons->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterKaons->SetBinError(i, newError);
    }
    h_t_pi12_afterKaons->SetMarkerStyle(22);
    h_t_pi12_afterKaons->SetMarkerColor(kP10Brown);
    h_t_pi12_afterKaons->SetLineColor(kP10Brown);
    h_t_pi12_afterKaons->Draw("PE same");

    h_t_pi12_afterMass->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterMass->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterMass->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterMass->SetBinError(i, newError);
    }
    h_t_pi12_afterMass->SetMarkerStyle(31);
    h_t_pi12_afterMass->SetMarkerColor(kP10Blue);
    h_t_pi12_afterMass->SetLineColor(kP10Blue);
    h_t_pi12_afterMass->Draw("PE same");

    h_t_pi12_afterOMD->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterOMD->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterOMD->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterOMD->SetBinError(i, newError);
    }
    h_t_pi12_afterOMD->SetMarkerStyle(45);
    h_t_pi12_afterOMD->SetMarkerColor(kP10Gray);
    h_t_pi12_afterOMD->SetLineColor(kP10Gray);
    h_t_pi12_afterOMD->Draw("PE same");

    h_t_REC_wRES_cut_pi12->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)*(M_PI/2)/(M_PI/12));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
    h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PE same");


    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.7, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.7, " Simulation 25.10.2/3");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.13, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.13, 0.85, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r4 = new TLatex(0.13, 0.8, "10x100 GeV");
	r4->SetNDC();
	r4->SetTextSize(30);
	r4->SetTextFont(43);
	r4->SetTextColor(kBlack);
	r4->Draw("same"); 

    TLatex* r42 = new TLatex(0.15, 0.19, "|t|_{#omega_{max}} normalized by #pi/2/#omega_{max}");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	//r42->Draw("same");

    TLatex* r421111 = new TLatex(0.15, 0.14, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    TLegend *w14_213 = new TLegend(0.5,0.49,0.9,0.89);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_rho_MC,"Coherent #rho: MC", "P");
    w14_213->AddEntry(h_t_after_VMrec,"After VM rec", "P");
    w14_213->AddEntry(h_t_pi12_afterHFS,"After requiring 2 trks", "P");
    w14_213->AddEntry(h_t_pi12_afterKaons,"After requiring K^{+}K^{-}", "P");
    w14_213->AddEntry(h_t_pi12_afterMass,"After mass selection", "P");
    w14_213->AddEntry(h_t_pi12_afterOMD,"After ZDC+RP+OMD", "P");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"After #omega_{max}=#pi/12", "P");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_coherent_rho_wPID.pdf");
}

void plot_coherent_rho_noPID()
{
    TFile* coherent_phi_allVetoes = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_MC = (TH1D*) coherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) coherent_phi_allVetoes->Get("h_Nevents");

    TFile* coherent_rho_allVetoes = TFile::Open("z_coherentRho_noPID_allVeotesCuts.root","READ");
    TH1D* h_t_rho_MC = (TH1D*) coherent_rho_allVetoes->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) coherent_rho_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_t_after_VMrec = (TH1D*) coherent_rho_allVetoes->Get("h_t_after_VMrec");
    TH1D* h_t_pi12_afterHFS = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterHFS");
    TH1D* h_t_pi12_afterKaons = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterKaons");
    TH1D* h_t_pi12_afterMass = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterMass");
    TH1D* h_t_pi12_afterOMD = (TH1D*) coherent_rho_allVetoes->Get("h_t_pi12_afterOMD");
    TH1D* h_rho_sartre_events = (TH1D*) coherent_rho_allVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_rho = 5418.68; // nb
    double nEvents_rho = h_rho_sartre_events->GetEntries();
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   
	
    TCanvas* c1 = new TCanvas("c1","c1",1,1,1200,800);
    c1->Divide(1,1,0.01,0.01);
    c1->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->SetTitle("#rho Distribution without PID");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherentMC_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherentMC_phi_preTDR = preTDR_lumi/simu_coherentMC_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherentMC_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineStyle(1);
    h_t_MC->SetLineColor(kBlack);
    h_t_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_MC->Draw("HISTsame");

    h_t_rho_MC->Rebin(rebin_width);
    h_t_rho_MC->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)/rebin_width);
    double simu_coherent_lumi = nEvents_rho/sigma_rho; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_rho_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_rho_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_rho_MC->SetBinError(i, newError);
    }
    h_t_rho_MC->SetMarkerStyle(28);
    h_t_rho_MC->SetMarkerColor(kP8Gray);
    h_t_rho_MC->SetLineColor(kP8Gray);
    h_t_rho_MC->GetYaxis()->SetRangeUser(1e-1,1e8);
    h_t_rho_MC->Draw("PE same");

    h_t_after_VMrec->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_after_VMrec->GetNbinsX(); ++i) 
    {
        double binError = h_t_after_VMrec->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_after_VMrec->SetBinError(i, newError);
    }
    h_t_after_VMrec->SetMarkerStyle(25);
    h_t_after_VMrec->SetMarkerColor(kP10Violet);
    h_t_after_VMrec->SetLineColor(kP10Violet);
    h_t_after_VMrec->Draw("PE same");

    h_t_pi12_afterHFS->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterHFS->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterHFS->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterHFS->SetBinError(i, newError);
    }
    h_t_pi12_afterHFS->SetMarkerStyle(4);
    h_t_pi12_afterHFS->SetMarkerColor(kP10Orange);
    h_t_pi12_afterHFS->SetLineColor(kP10Orange);
    h_t_pi12_afterHFS->Draw("PE same");

    h_t_pi12_afterKaons->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterKaons->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterKaons->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterKaons->SetBinError(i, newError);
    }
    h_t_pi12_afterKaons->SetMarkerStyle(22);
    h_t_pi12_afterKaons->SetMarkerColor(kP10Brown);
    h_t_pi12_afterKaons->SetLineColor(kP10Brown);
    h_t_pi12_afterKaons->Draw("PE same");

    h_t_pi12_afterMass->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterMass->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterMass->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterMass->SetBinError(i, newError);
    }
    h_t_pi12_afterMass->SetMarkerStyle(31);
    h_t_pi12_afterMass->SetMarkerColor(kP10Blue);
    h_t_pi12_afterMass->SetLineColor(kP10Blue);
    h_t_pi12_afterMass->Draw("PE same");

    h_t_pi12_afterOMD->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterOMD->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterOMD->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_pi12_afterOMD->SetBinError(i, newError);
    }
    h_t_pi12_afterOMD->SetMarkerStyle(45);
    h_t_pi12_afterOMD->SetMarkerColor(kP10Gray);
    h_t_pi12_afterOMD->SetLineColor(kP10Gray);
    h_t_pi12_afterOMD->Draw("PE same");

    h_t_REC_wRES_cut_pi12->Scale(sigma_rho*(1/nEvents_rho)*(1/dN_dt)*(M_PI/2)/(M_PI/12));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
    h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PE same");


    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.7, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.7, " Simulation 25.10.2/3");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.13, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.13, 0.85, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r4 = new TLatex(0.13, 0.8, "10x100 GeV");
	r4->SetNDC();
	r4->SetTextSize(30);
	r4->SetTextFont(43);
	r4->SetTextColor(kBlack);
	r4->Draw("same"); 

    TLatex* r42 = new TLatex(0.15, 0.19, "|t|_{#omega_{max}} normalized by #pi/2/#omega_{max}");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	//r42->Draw("same");

    TLatex* r421111 = new TLatex(0.15, 0.14, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    TLegend *w14_213 = new TLegend(0.5,0.49,0.9,0.89);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_rho_MC,"Coherent #rho: MC", "P");
    w14_213->AddEntry(h_t_after_VMrec,"After VM rec", "P");
    w14_213->AddEntry(h_t_pi12_afterHFS,"After requiring 2 trks", "P");
    w14_213->AddEntry(h_t_pi12_afterKaons,"After requiring K^{+}K^{-}", "P");
    w14_213->AddEntry(h_t_pi12_afterMass,"After mass selection", "P");
    w14_213->AddEntry(h_t_pi12_afterOMD,"After ZDC+RP+OMD", "P");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"After #omega_{max}=#pi/12", "P");
    w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);   
    w14_213->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_coherent_rho_noPID.pdf");
}

void plot_incoherent_phi_noVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_noVetoes = TFile::Open("z_incoherentPhi_allVeotesCuts.root","READ");
    TH1D* h_incoherent_phi_noVetoes = (TH1D*) incoherent_phi_noVetoes->Get("h_t_MC");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_noVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_incoherent_phi_noVetoes->Rebin(rebin_width);
    h_incoherent_phi_noVetoes->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_noVetoes->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_noVetoes->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_noVetoes->SetBinError(i, newError);
    }
    h_incoherent_phi_noVetoes->SetMarkerStyle(7);
    h_incoherent_phi_noVetoes->SetMarkerColor(kP10Orange);
    h_incoherent_phi_noVetoes->SetLineColor(kP10Orange);
    h_incoherent_phi_noVetoes->Draw("PEsame");
    
	TLegend *w14_213 = new TLegend(0.6,0.7,0.72,0.9);
    w14_213->AddEntry(h_t_MC,"Coherent #phi: MC", "L");
    w14_213->AddEntry(h_incoherent_phi_noVetoes,"Incoherent #phi: MC", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with Incoherent #phi No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_t_incoherent_noVetoes_preTDR.pdf");
}

void plot_incoherent_together_preTDR()
{
    TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allVetoes = TFile::Open("z_incoherentPhi_allVeotesCuts.root","READ");
    TH1D* h_t_incoherent_MC = (TH1D*) incoherent_phi_allVetoes->Get("h_t_MC");
    TH1D* h_t_pi12_afterElectronSelect = (TH1D*) incoherent_phi_allVetoes->Get("h_t_pi12_afterElectronSelect");
    TH1D* h_t_pi12_afterHFS = (TH1D*) incoherent_phi_allVetoes->Get("h_t_pi12_afterHFS");
    TH1D* h_t_pi12_afterMass = (TH1D*) incoherent_phi_allVetoes->Get("h_t_pi12_afterMass");
    TH1D* h_t_pi12_afterOMD = (TH1D*) incoherent_phi_allVetoes->Get("h_t_pi12_afterOMD");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

	TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)]^{2}");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineWidth(1);
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_incoherent_MC->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_t_incoherent_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_incoherent_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_t_incoherent_MC->SetBinError(i, newError);
    }
    h_t_incoherent_MC->SetMarkerStyle(7);
    h_t_incoherent_MC->SetMarkerColor(kP10Orange);
    h_t_incoherent_MC->SetLineColor(kP10Orange);
    h_t_incoherent_MC->Draw("PEsame");

    h_t_pi12_afterElectronSelect->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterElectronSelect->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterElectronSelect->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_t_pi12_afterElectronSelect->SetBinError(i, newError);
    }
    h_t_pi12_afterElectronSelect->SetMarkerStyle(2);
    h_t_pi12_afterElectronSelect->SetMarkerColor(kP10Blue);
    h_t_pi12_afterElectronSelect->SetLineColor(kP10Blue);
    h_t_pi12_afterElectronSelect->Draw("PEsame");

    h_t_pi12_afterHFS->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterHFS->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterHFS->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_t_pi12_afterHFS->SetBinError(i, newError);
    }
    h_t_pi12_afterHFS->SetMarkerStyle(4);
    h_t_pi12_afterHFS->SetMarkerColor(kP10Yellow);
    h_t_pi12_afterHFS->SetLineColor(kP10Yellow);
    h_t_pi12_afterHFS->Draw("PEsame");

    h_t_pi12_afterMass->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterMass->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterMass->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_t_pi12_afterMass->SetBinError(i, newError);
    }
    h_t_pi12_afterMass->SetMarkerStyle(31);
    h_t_pi12_afterMass->SetMarkerColor(kP10Violet);
    h_t_pi12_afterMass->SetLineColor(kP10Violet);
    h_t_pi12_afterMass->Draw("PEsame");

    h_t_pi12_afterOMD->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_pi12_afterOMD->GetNbinsX(); ++i) 
    {
        double binError = h_t_pi12_afterOMD->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_t_pi12_afterOMD->SetBinError(i, newError);
    }
    h_t_pi12_afterOMD->SetMarkerStyle(48);
    h_t_pi12_afterOMD->SetMarkerColor(kP10Brown);
    h_t_pi12_afterOMD->SetLineColor(kP10Brown);
    h_t_pi12_afterOMD->Draw("PEsame");

    TLegend *w14_213 = new TLegend(0.6,0.6,0.72,0.9);
    w14_213->AddEntry(h_t_MC," Coherent #phi: MC", "L");
    w14_213->AddEntry(h_t_incoherent_MC,"Incoh. #phi: MC", "P");
    w14_213->AddEntry(h_t_pi12_afterElectronSelect,"After e' selection", "P");
    w14_213->AddEntry(h_t_pi12_afterHFS,"After requiring 2 trks", "P");
    w14_213->AddEntry(h_t_pi12_afterMass,"After mass selection", "P");
    w14_213->AddEntry(h_t_pi12_afterOMD,"After ZDC+RP+OMD", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
	w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./figures/analysisNote_plots/plot_t_incoherent_compareAll_preTDR_proj.pdf");
}

void plot_incoherent_phi_allVetoes_preTDR()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

	TFile* incoherent_phi_allVetoes = TFile::Open("z_incoherentPhi_allVeotesCuts.root","READ");
    TH1D* h_incoherent_phi_allVetoes_proj = (TH1D*) incoherent_phi_allVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_incoherent_events = (TH1D*) incoherent_phi_allVetoes->Get("h_Nevents");

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_incoherent = 261.95; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_incoherent = h_incoherent_events->GetEntries();
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   
	
	TCanvas* c14_2131 = new TCanvas("c14_2131","c14_2131",1,1,1200,800);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
    h_t_REC_wRES_cut_pi12->GetYaxis()->SetRangeUser(1e-1,1e6);
    h_t_REC_wRES_cut_pi12->SetLineColor(kBlack);
	h_t_REC_wRES_cut_pi12->Draw("HISTsame");

    h_incoherent_phi_allVetoes_proj->Scale(sigma_incoherent*(1/nEvents_incoherent)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_incoherent_lumi = nEvents_incoherent/sigma_incoherent; // nb^-1 
    double ratio_incoherent_preTDR = preTDR_lumi/simu_incoherent_lumi;
    for (int i = 1; i <= h_incoherent_phi_allVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_incoherent_phi_allVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_incoherent_preTDR);
        h_incoherent_phi_allVetoes_proj->SetBinError(i, newError);
    }
    h_incoherent_phi_allVetoes_proj->SetMarkerStyle(7);
    h_incoherent_phi_allVetoes_proj->SetMarkerColor(kP10Orange);
    h_incoherent_phi_allVetoes_proj->SetLineColor(kP10Orange);
    h_incoherent_phi_allVetoes_proj->Draw("PEsame");

	TLegend *w14_2131 = new TLegend(0.6,0.7,0.72,0.9);
    w14_2131->AddEntry(h_t_REC_wRES_cut_pi12,"Coherent #phi: #omega_{max}= #pi/12", "L");    
    w14_2131->AddEntry(h_incoherent_phi_allVetoes_proj,"Incoherent #phi: #omega_{max}= #pi/12", "P");
	w14_2131->SetBorderSize(0);   
    w14_2131->SetFillStyle(0);
    w14_2131->SetTextSize(30);
	w14_2131->SetTextFont(45);
    w14_2131->Draw("same");	

    c14_2131->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_2131 = new TLatex();
    title14_2131->SetNDC(); 
    title14_2131->SetTextSize(0.05);
    title14_2131->SetTextAlign(22);  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.18, 0.13, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_2131->Print("./figures/analysisNote_plots/plot_t_incoherent_allvetoes_preTDR_proj.pdf");
}

void plot_tdist_projOnly_preTDR()
{
    TFile* file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = (TH1D*)file->Get("h_t_MC");
	TH1D* h_t_REC_proj12 = (TH1D*)file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*)file->Get("h_Nevents");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");
	
	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r142111 = new TLatex(0.18, 0.18, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r142111->SetNDC();
	r142111->SetTextSize(25);
	r142111->SetTextFont(43);
	r142111->SetTextColor(kBlack);
	r142111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.75,0.73,0.85);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #omega_{max} = #pi/12", "P");
	w7->Draw("same");

	c1->Print("./figures/analysisNote_plots/plot_t_dist_preTDR_projOnly.pdf");
}

void plot_transform_allMethods()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* hdsigmadt_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* hdsigmadt_REC = (TH1D*) phi_t_file->Get("h_t_REC_EEMC_cut");
	TH1D* hdsigmadt_REC_new = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");

	const int nTrials = 1000; 
	const double hbarc = 0.197;
	const double t_cut = 0.2;
	const double bmin = -12;
	const double bmax = 12;
	const int noOfBins = 300;

	int nbins = hdsigmadt_MC->GetNbinsX();
	TRandom3 randGen(0); // Random seed

	TH1D* hF_b_MC_2d = new TH1D("hF_b_MC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_2d = new TH1D("hF_b_REC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_new_2d = new TH1D("hF_b_REC_new_2d", "", noOfBins, bmin, bmax);

	vector<vector<double>> trials_MC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC_new(noOfBins, vector<double>(nTrials));

	for (int trial = 0; trial < nTrials; ++trial) 
	{
    	for (int j = 1; j <= noOfBins; ++j) 
    	{
        	double b_2d = hF_b_MC_2d->GetBinCenter(j);
        	double prefactor = 1.0 / (2*TMath::Pi());

        	double F_b_MC_2d = 0, F_b_REC_2d = 0, F_b_REC_new_2d = 0;

        	for (int i = 1; i <= nbins; ++i) 
        	{
            	double tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            	double t = hdsigmadt_MC->GetBinCenter(i);
            	double delta = sqrt(fabs(t));

            	// Gaussian sampling
            	double dsigmadt_MC = randGen.Gaus(hdsigmadt_MC->GetBinContent(i), sqrt(60)*hdsigmadt_MC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC = randGen.Gaus(hdsigmadt_REC->GetBinContent(i), sqrt(60)*hdsigmadt_REC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC_new = randGen.Gaus(hdsigmadt_REC_new->GetBinContent(i), sqrt(60)*hdsigmadt_REC_new->GetBinError(i)) / 1e7;

            	double bessel = TMath::BesselJ0(b_2d * delta / hbarc);

            	if (t > t_cut) continue;

            	double amp_MC = dsigmadt_MC > 0 ? sqrt(dsigmadt_MC) : 0;
            	double amp_REC = dsigmadt_REC > 0 ? sqrt(dsigmadt_REC) : 0;
            	double amp_REC_new = dsigmadt_REC_new > 0 ? sqrt(dsigmadt_REC_new) : 0;

            	if (t > 0.014) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.048) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.098) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t > 0.17) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }

            	F_b_MC_2d += amp_MC * bessel * tBinWidth / 2;
            	F_b_REC_2d += amp_REC * bessel * tBinWidth / 2;
            	F_b_REC_new_2d += amp_REC_new * bessel * tBinWidth / 2;
        	}

        	F_b_MC_2d *= prefactor / hbarc;
        	F_b_REC_2d *= prefactor / hbarc;
        	F_b_REC_new_2d *= prefactor / hbarc;

        	trials_MC[j - 1][trial] = F_b_MC_2d;
        	trials_REC[j - 1][trial] = F_b_REC_2d;
        	trials_REC_new[j - 1][trial] = F_b_REC_new_2d;
    	}
	}

	for (int j = 1; j <= noOfBins; ++j) 
	{
    	auto& vec_MC = trials_MC[j - 1];
    	auto& vec_REC = trials_REC[j - 1];
    	auto& vec_REC_new = trials_REC_new[j - 1];

    	double mean_MC = TMath::Mean(nTrials, vec_MC.data());
    	double std_MC = TMath::RMS(nTrials, vec_MC.data());
    	double mean_REC = TMath::Mean(nTrials, vec_REC.data());
    	double std_REC = TMath::RMS(nTrials, vec_REC.data());
    	double mean_REC_new = TMath::Mean(nTrials, vec_REC_new.data());
    	double std_REC_new = TMath::RMS(nTrials, vec_REC_new.data());

    	hF_b_MC_2d->SetBinContent(j, mean_MC);
    	hF_b_MC_2d->SetBinError(j, std_MC);

    	hF_b_REC_2d->SetBinContent(j, mean_REC);
    	hF_b_REC_2d->SetBinError(j, std_REC);

    	hF_b_REC_new_2d->SetBinContent(j, mean_REC_new);
    	hF_b_REC_new_2d->SetBinError(j, std_REC_new);
	}

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);
	gStyle->SetOptStat(0);
        
    hF_b_MC_2d->Scale(1.0 / hF_b_MC_2d->Integral("width"));
    hF_b_MC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC_2d->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC_2d->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_MC_2d->SetLineColor(kBlack);
    hF_b_MC_2d->SetLineWidth(4);
    hF_b_MC_2d->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC_2d->GetYaxis()->SetLabelSize(0.04);

	hF_b_REC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_2d->Scale(1.0 / hF_b_REC_2d->Integral("width"));
    hF_b_REC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_REC_2d->SetMarkerStyle(20); 
    hF_b_REC_2d->SetMarkerColor(kP8Blue);
    hF_b_REC_2d->SetLineColor(kP8Blue);

	hF_b_REC_new_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_new_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_new_2d->Scale(1.0 / hF_b_REC_new_2d->Integral("width"));
    hF_b_REC_new_2d->SetMarkerStyle(30); 
    hF_b_REC_new_2d->SetMarkerColor(kP8Pink);
    hF_b_REC_new_2d->SetLineColor(kP8Pink);
 
    hF_b_REC_2d->GetYaxis()->SetRangeUser(-0.07, 0.12);
    hF_b_REC_2d->Draw();
    hF_b_REC_new_2d->Draw("PEsame"); 
    hF_b_MC_2d->SetMarkerStyle(0);
    hF_b_MC_2d->Draw("Lsame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.22, 0.3, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.28, 0.3, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.22, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.22, 0.25, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLegend *leg = new TLegend(0.68,0.75,0.72,0.85);
    leg->AddEntry(hF_b_MC_2d, " MC", "L");
    leg->AddEntry(hF_b_REC_2d, " Method L", "PE"); 
    leg->AddEntry(hF_b_REC_new_2d, " Projection method", "PE"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(45);       
    leg->SetTextSize(20);    
    leg->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_Fb_2d_allMethods.pdf");

    cout << "All done. Bye." << endl;
}

void plot_transform_wSTAR()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* hdsigmadt_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* hdsigmadt_REC_new = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");

	const int nTrials = 1000; 
	const double hbarc = 0.197;
	const double t_cut = 0.2;
	const double bmin = -12;
	const double bmax = 12;
	const int noOfBins = 300;

	int nbins = hdsigmadt_MC->GetNbinsX();
	TRandom3 randGen(0); // Random seed

	TH1D* hF_b_MC_2d = new TH1D("hF_b_MC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_new_2d = new TH1D("hF_b_REC_new_2d", "", noOfBins, bmin, bmax);

	vector<vector<double>> trials_MC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC_new(noOfBins, vector<double>(nTrials));

	for (int trial = 0; trial < nTrials; ++trial) 
	{
    	for (int j = 1; j <= noOfBins; ++j) 
    	{
        	double b_2d = hF_b_MC_2d->GetBinCenter(j);
        	double prefactor = 1.0 / (2*TMath::Pi());

        	double F_b_MC_2d = 0, F_b_REC_new_2d = 0;

        	for (int i = 1; i <= nbins; ++i) 
        	{
            	double tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            	double t = hdsigmadt_MC->GetBinCenter(i);
            	double delta = sqrt(fabs(t));

            	// Gaussian sampling
            	double dsigmadt_MC = randGen.Gaus(hdsigmadt_MC->GetBinContent(i), sqrt(60)*hdsigmadt_MC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC_new = randGen.Gaus(hdsigmadt_REC_new->GetBinContent(i), sqrt(60)*hdsigmadt_REC_new->GetBinError(i)) / 1e7;

            	double bessel = TMath::BesselJ0(b_2d * delta / hbarc);

            	if (t > t_cut) continue;

            	double amp_MC = dsigmadt_MC > 0 ? sqrt(dsigmadt_MC) : 0;
            	double amp_REC_new = dsigmadt_REC_new > 0 ? sqrt(dsigmadt_REC_new) : 0;

            	if (t > 0.014) { amp_MC *= -1; amp_REC_new *= -1; }
            	if (t > 0.048) { amp_MC *= -1; amp_REC_new *= -1; }
            	if (t > 0.098) { amp_MC *= -1; amp_REC_new *= -1; }
                if (t > 0.17) { amp_MC *= -1; amp_REC_new *= -1; }

            	F_b_MC_2d += amp_MC * bessel * tBinWidth / 2;
            	F_b_REC_new_2d += amp_REC_new * bessel * tBinWidth / 2;
        	}

        	F_b_MC_2d *= prefactor / hbarc;
        	F_b_REC_new_2d *= prefactor / hbarc;

        	trials_MC[j - 1][trial] = F_b_MC_2d;
        	trials_REC_new[j - 1][trial] = F_b_REC_new_2d;
    	}
	}

	for (int j = 1; j <= noOfBins; ++j) 
	{
    	auto& vec_MC = trials_MC[j - 1];
    	auto& vec_REC_new = trials_REC_new[j - 1];

    	double mean_MC = TMath::Mean(nTrials, vec_MC.data());
    	double std_MC = TMath::RMS(nTrials, vec_MC.data());
    	double mean_REC_new = TMath::Mean(nTrials, vec_REC_new.data());
    	double std_REC_new = TMath::RMS(nTrials, vec_REC_new.data());

    	hF_b_MC_2d->SetBinContent(j, mean_MC);
    	hF_b_MC_2d->SetBinError(j, std_MC);

    	hF_b_REC_new_2d->SetBinContent(j, mean_REC_new);
    	hF_b_REC_new_2d->SetBinError(j, std_REC_new);
	}

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);
	gStyle->SetOptStat(0);

    hF_b_MC_2d->Scale(1.0 / hF_b_MC_2d->Integral("width"));
    hF_b_MC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC_2d->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC_2d->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_MC_2d->SetLineColor(kBlack);
    hF_b_MC_2d->SetLineWidth(4);
    hF_b_MC_2d->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC_2d->GetYaxis()->SetLabelSize(0.04);

	hF_b_REC_new_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_new_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_new_2d->Scale(1.0 / hF_b_REC_new_2d->Integral("width"));
    hF_b_REC_new_2d->SetMarkerStyle(30); 
    hF_b_REC_new_2d->SetMarkerColor(kP8Pink);
    hF_b_REC_new_2d->SetLineColor(kP8Pink);

    TFile* fSTAR = TFile::Open("STAR_Fb.root");
    TGraphErrors* gSTAR = (TGraphErrors*) fSTAR->Get("STAR_Fb");

    gSTAR->SetMarkerStyle(20);
    gSTAR->SetMarkerColor(kBlue);
    gSTAR->SetLineColor(kBlue);
    gSTAR->Scale(1.0 / gSTAR->Integral());
 
    hF_b_REC_new_2d->GetYaxis()->SetRangeUser(-0.07, 0.12);
    hF_b_REC_new_2d->Draw("PEsame"); 
    gSTAR->Draw("P SAME");
    hF_b_MC_2d->SetMarkerStyle(0);
    hF_b_MC_2d->Draw("Lsame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.22, 0.3, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.28, 0.3, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.22, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.22, 0.25, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLegend *leg = new TLegend(0.68,0.75,0.72,0.85);
    leg->AddEntry(hF_b_MC_2d, " MC", "L");
    leg->AddEntry(hF_b_REC_new_2d, " Projection method", "PE"); 
    leg->AddEntry(gSTAR, " STAR XnXn", "PE"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(45);       
    leg->SetTextSize(20);    
    leg->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_Fb_2d_wSTAR.pdf");

    cout << "All done. Bye." << endl;
}

void plot_transform_projOnly()
{
	TFile* phi_t_file = TFile::Open("z_coherentPhi_allVeotesCuts.root","READ");
	TH1D* hdsigmadt_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* hdsigmadt_REC_new = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");

	const int nTrials = 1000; 
	const double hbarc = 0.197;
	const double t_cut = 0.2;
	const double bmin = -12;
	const double bmax = 12;
	const int noOfBins = 300;

	int nbins = hdsigmadt_MC->GetNbinsX();
	TRandom3 randGen(0); // Random seed

	TH1D* hF_b_MC_2d = new TH1D("hF_b_MC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_new_2d = new TH1D("hF_b_REC_new_2d", "", noOfBins, bmin, bmax);

	vector<vector<double>> trials_MC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC_new(noOfBins, vector<double>(nTrials));

	for (int trial = 0; trial < nTrials; ++trial) 
	{
    	for (int j = 1; j <= noOfBins; ++j) 
    	{
        	double b_2d = hF_b_MC_2d->GetBinCenter(j);
        	double prefactor = 1.0 / (2*TMath::Pi());

        	double F_b_MC_2d = 0, F_b_REC_2d = 0, F_b_REC_new_2d = 0;

        	for (int i = 1; i <= nbins; ++i) 
        	{
            	double tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            	double t = hdsigmadt_MC->GetBinCenter(i);
            	double delta = sqrt(fabs(t));

            	// Gaussian sampling
            	double dsigmadt_MC = randGen.Gaus(hdsigmadt_MC->GetBinContent(i), sqrt(60)*hdsigmadt_MC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC_new = randGen.Gaus(hdsigmadt_REC_new->GetBinContent(i), sqrt(60)*hdsigmadt_REC_new->GetBinError(i)) / 1e7;

            	double bessel = TMath::BesselJ0(b_2d * delta / hbarc);

            	if (t > t_cut) continue;

            	double amp_MC = dsigmadt_MC > 0 ? sqrt(dsigmadt_MC) : 0;
            	double amp_REC_new = dsigmadt_REC_new > 0 ? sqrt(dsigmadt_REC_new) : 0;

            	if (t > 0.014) { amp_MC *= -1; amp_REC_new *= -1; }
            	if (t > 0.048) { amp_MC *= -1; amp_REC_new *= -1; }
            	if (t > 0.098) { amp_MC *= -1; amp_REC_new *= -1; }
                if (t > 0.17) { amp_MC *= -1; amp_REC_new *= -1; }

            	F_b_MC_2d += amp_MC * bessel * tBinWidth / 2;
            	F_b_REC_new_2d += amp_REC_new * bessel * tBinWidth / 2;
        	}

        	F_b_MC_2d *= prefactor / hbarc;
        	F_b_REC_new_2d *= prefactor / hbarc;

        	trials_MC[j - 1][trial] = F_b_MC_2d;
        	trials_REC_new[j - 1][trial] = F_b_REC_new_2d;
    	}
	}

	for (int j = 1; j <= noOfBins; ++j) 
	{
    	auto& vec_MC = trials_MC[j - 1];
    	auto& vec_REC_new = trials_REC_new[j - 1];

    	double mean_MC = TMath::Mean(nTrials, vec_MC.data());
    	double std_MC = TMath::RMS(nTrials, vec_MC.data());
    	double mean_REC_new = TMath::Mean(nTrials, vec_REC_new.data());
    	double std_REC_new = TMath::RMS(nTrials, vec_REC_new.data());

    	hF_b_MC_2d->SetBinContent(j, mean_MC);
    	hF_b_MC_2d->SetBinError(j, std_MC);

    	hF_b_REC_new_2d->SetBinContent(j, mean_REC_new);
    	hF_b_REC_new_2d->SetBinError(j, std_REC_new);
	}

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);
	gStyle->SetOptStat(0);
        
    hF_b_MC_2d->Scale(1.0 / hF_b_MC_2d->Integral("width"));
    hF_b_MC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC_2d->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC_2d->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_MC_2d->SetLineColor(kBlack);
    hF_b_MC_2d->SetLineWidth(4);
    hF_b_MC_2d->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC_2d->GetYaxis()->SetLabelSize(0.04);

	hF_b_REC_new_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_new_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_new_2d->Scale(1.0 / hF_b_REC_new_2d->Integral("width"));
    hF_b_REC_new_2d->SetMarkerStyle(30); 
    hF_b_REC_new_2d->SetMarkerColor(kP8Pink);
    hF_b_REC_new_2d->SetLineColor(kP8Pink);

    hF_b_REC_new_2d->GetYaxis()->SetRangeUser(-0.07, 0.12);
    hF_b_REC_new_2d->Draw("PEsame"); 
    hF_b_MC_2d->SetMarkerStyle(0);
    hF_b_MC_2d->Draw("Lsame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.22, 0.3, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.28, 0.3, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.22, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.22, 0.25, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLegend *leg = new TLegend(0.68,0.75,0.72,0.85);
    leg->AddEntry(hF_b_MC_2d, " MC", "L");
    leg->AddEntry(hF_b_REC_new_2d, " Projection method", "PE"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(45);       
    leg->SetTextSize(20);    
    leg->Draw("same");

    c1->Print("./figures/analysisNote_plots/plot_Fb_2d_projOnly.pdf");

    cout << "All done. Bye." << endl;
}

