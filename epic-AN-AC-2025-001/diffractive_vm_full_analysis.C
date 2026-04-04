#include "pleaseIncludeMe.h"
#include <TVector3.h>
#include <TVector.h>
#include <cmath>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"

/*-----------------------------------------------------
    Analysis for Exclusive Diffractive VM Production
- Purpose: reconstuct kinematic variables and the t-distribution
    using methods E, L, and projection with wedge cut
- Input: root file/s
- Output: root file
------------------------------------------------------*/

using namespace std;

/*-----------------------------
    For incoherent analysis
------------------------------*/
struct Cluster_EEMC
{
    float x, y, z, energy;
};

struct Cluster_ZDC 
{
    float x, y, z, energy;
};

struct Hit_neutral_ZDC 
{
    float x, y, z, energy;
};

struct Hit_REC_ZDC 
{
    float x, y, z, energy;
};

struct Hit_RP 
{
    float x, y, z, energy;
};

struct Hit_OMD 
{
    float x, y, z, energy;
};

struct Event 
{
    vector<Cluster_EEMC> clusters_eemc;
    vector<Hit_neutral_ZDC> zdc_neutrals;        
    vector<Hit_REC_ZDC> zdc_REC_hits;
    vector<Cluster_ZDC> zdc_clusters;    // combined ECal + HCal clusters
    vector<Hit_RP> hit_rp;
    vector<Hit_OMD> hit_omd;
};

/*----------------------------
    for kaon identification
------------------------------*/
struct KaonCand 
{
    TLorentzVector pKaon;
    int pdg;
};

auto giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    /*---------------------------------------------------------
        Function to calculate method L for t reconstruction
    Input: (four momenta) P_e, P_e', P_A, P_VM
    Returns: double, |t| = -(P_A'corr - P_A)^2 
    ----------------------------------------------------------*/
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	return method_L; 
}

auto giveme_t_new_method(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn, TLorentzVector vmOut)
{
    /*---------------------------------------------------------
        Function to calculate method L for t reconstruction
    Input: four momenta, P_e, P_e', P_A, P_VM
    Returns: four vector, P = sqrt(|t|) = -(P_A'corr - P_A)
    ----------------------------------------------------------*/
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();
	TLorentzVector method_L_4vect = -(a_beam_scattered_corr-aInVec); 
	return method_L_4vect;
}

int diffractive_vm_full_analysis(TString rec_file, TString mode, TString cutMode, TString vetoes, TString pid, TString MCcuts, TString outputfile)
{	
    /*---------------------------------------------------------
        Macro to 
    Input: 
    Returns: 
    ----------------------------------------------------------*/

    // reads in lists of root files	
    TString listname = rec_file;
    TString outname = outputfile;
    TChain *chain = new TChain("events");
    int nfiles = 0;
    char filename[512];
    ifstream *inputstream = new ifstream;
    inputstream->open(listname.Data());
    if(!inputstream)
    {
      printf("[e] Cannot open file list: %s\n", listname.Data());
    }
    cout << "Running analysis mode = " << mode << endl;
    while(inputstream->good())
    {
        inputstream->getline(filename, 512);
        if(inputstream->good())
        {
            TFile *ftmp = TFile::Open(filename, "read");
            if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) 
            {
                printf("[e] Could you open file: %s\n", filename);
            } 
            else
            {
                cout<<"[i] Add "<<nfiles<<"th file: "<<filename<<endl;
                chain->Add(filename);
                nfiles++;
            }
        }
    }
    inputstream->close();
    printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());
    TTreeReader tree_reader(chain);

    /*-----------------------------------
        Define different analysis modes
    ------------------------------------*/
    bool isDIS = false;
    bool isRho = false;
    bool isPhi = false;
    bool useZDC = true;
    bool applyCuts = true; 
    bool applyVetoes = true; 
    bool usePID = true;
    bool applyMCcuts = true;

    // choose "mode" argument
    if (mode == "DIS") 
    {
        isDIS = true;
        useZDC = false;   // DIS doesn't have ZDC yet
        cout << "DIS mode, no ZDC" << endl;
    }
    else if (mode == "rho") 
    {
        isRho = true; // uses pions for MC VM
        cout << "rho mode, use pions for MC" << endl;
    }
    else if (mode == "phi_coh") // coherent phi Sartre files
    {
        isPhi = true; // uses kaons for MC VM 
        cout << "coherent phi mode" << endl;
    }
    else if (mode == "phi_incoh") // incoherent phi BeAGLE files
    {
        isPhi = true; // uses kaons for MC VM
        cout << "incoherent phi mode" << endl;
    }
    else 
    {
        isPhi = true;
        cout << "[w] Unknown mode, defaulting to phi analysis." << endl;
    }

    // choose "cutMode" argument
    if (cutMode == "noCuts") 
    {
        applyCuts = false; // E-pz, E/p, E==0, HFS==2, kaons==2, VMmass
        cout << "no selection cuts applied" << endl;
    }
    else if (cutMode == "allCuts")
    {
        applyCuts = true;
        cout << "all selection cuts applied" << endl;
    }
    else
    {
        applyCuts = true;
        cout << "[w] Unknown cut mode, defaulting to all cuts" << endl;
    }

    // choose "vetoes" argument
    if (vetoes == "noVetoes") 
    {
        applyVetoes = false; // ZDC, OMD, RP
        cout << "no detector vetoes applied" << endl;
    }
    else if (vetoes == "allVetoes")
    {
        applyVetoes = true;
        cout << "all detector vetoes applied" << endl;
    }
    else
    {
        applyVetoes = true;
        cout << "[w] Unknown veto mode, defaulting to all detector vetoes" << endl;
    }

    if (pid == "noPID")
    {
        usePID = false;
        cout << "no PID used" << endl;
    }
    else if (pid == "withPID")
    {
        usePID = true;
        cout << "PID used" << endl;
    }
    else
    {
        usePID = true;
        cout << "[w] Unknown PID selection, defaulting to using PID" << endl;
    }

    if (MCcuts=="yes")
    {
        applyMCcuts = true;
        cout << "Applying MC level selection cuts" << endl;
    }
    else if (MCcuts=="no")
    {
        applyMCcuts = false;
        cout << "No MC level selection cuts" << endl;
    }
    else
    {
        applyMCcuts = false;
        cout << "[w] Unknown MC cuts selection, defaulting to no MC cuts" << endl;
    }


    // MC particle arrays for each MC particle
    TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<double> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

    // Reconstructed EcalEndcapNClusters arrays
    TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
    TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
    TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};
    TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
    TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
    TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};
    TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClusterAssociations.recID"};
    TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClusterAssociations.simID"};

    // Reconstructed Calorimeter track arrays
    TTreeReaderArray<float> emCal_trk_x_array = {tree_reader, "_CalorimeterTrackProjections_points.position.x"};
    TTreeReaderArray<float> emCal_trk_y_array = {tree_reader, "_CalorimeterTrackProjections_points.position.y"};
    TTreeReaderArray<float> emCal_trk_px_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.x"};
    TTreeReaderArray<float> emCal_trk_py_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.y"};
    TTreeReaderArray<float> emCal_trk_pz_array = {tree_reader, "_CalorimeterTrackProjections_points.momentum.y"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
    TTreeReaderArray<int>   reco_pdg_array = {tree_reader, "ReconstructedChargedParticles.PDG"};
    TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};
    TTreeReaderArray<int> rec_id = {tree_reader, "_ReconstructedChargedParticleAssociations_rec.index"};
    TTreeReaderArray<int> sim_id = {tree_reader, "_ReconstructedChargedParticleAssociations_sim.index"};

    // Reconstructed Far Forward ZDC Neutrals (particles)
    TTreeReaderArray<float> zdc_x_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.x"};
    TTreeReaderArray<float> zdc_y_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.y"};
    TTreeReaderArray<float> zdc_z_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.referencePoint.z"};
    TTreeReaderArray<float> zdc_energy_array = {tree_reader, "ReconstructedFarForwardZDCNeutrals.energy"};

    // Reconstructed ECal ZDC Clusters
    TTreeReaderArray<float> Ecal_zdc_x_array = {tree_reader, "EcalFarForwardZDCClusters.position.x"};
    TTreeReaderArray<float> Ecal_zdc_y_array = {tree_reader, "EcalFarForwardZDCClusters.position.y"};
    TTreeReaderArray<float> Ecal_zdc_z_array = {tree_reader, "EcalFarForwardZDCClusters.position.z"};
    TTreeReaderArray<float> Ecal_zdc_energy_array = {tree_reader, "EcalFarForwardZDCClusters.energy"};

    // Reconstructed HCal ZDC Clusters
    TTreeReaderArray<float> Hcal_zdc_x_array = {tree_reader, "HcalFarForwardZDCClusters.position.x"};
    TTreeReaderArray<float> Hcal_zdc_y_array = {tree_reader, "HcalFarForwardZDCClusters.position.y"};
    TTreeReaderArray<float> Hcal_zdc_z_array = {tree_reader, "HcalFarForwardZDCClusters.position.z"};
    TTreeReaderArray<float> Hcal_zdc_energy_array = {tree_reader, "HcalFarForwardZDCClusters.energy"};

    // Reconstructed HCal ZDC hits
    TTreeReaderArray<float> REChits_zdc_x_array = {tree_reader, "HcalFarForwardZDCRecHits.position.x"};
    TTreeReaderArray<float> REChits_zdc_y_array = {tree_reader, "HcalFarForwardZDCRecHits.position.y"};
    TTreeReaderArray<float> REChits_zdc_z_array = {tree_reader, "HcalFarForwardZDCRecHits.position.z"};
    TTreeReaderArray<float> REChits_zdc_energy_array = {tree_reader, "HcalFarForwardZDCRecHits.energy"};

    // Reconstructed HCal ZDC raw hits
    TTreeReaderArray<int> RAWhits_zdc_amp_array = {tree_reader, "HcalFarForwardZDCRawHits.amplitude"};

    // Forward Roman Pot Rec Hits
    TTreeReaderArray<float> rp_x_array = {tree_reader, "ForwardRomanPotRecHits.position.x"};
    TTreeReaderArray<float> rp_y_array = {tree_reader, "ForwardRomanPotRecHits.position.y"};
    TTreeReaderArray<float> rp_z_array = {tree_reader, "ForwardRomanPotRecHits.position.z"};

    // Forward Off Momentum Tracker Rec Hits
    TTreeReaderArray<float> omd_x_array = {tree_reader, "ForwardOffMTrackerRecHits.position.x"};
    TTreeReaderArray<float> omd_y_array = {tree_reader, "ForwardOffMTrackerRecHits.position.y"};
    TTreeReaderArray<float> omd_z_array = {tree_reader, "ForwardOffMTrackerRecHits.position.z"};

    TString output_name_dir = outputfile+"_output.root";
    cout << "Output file = " << output_name_dir << endl;
    TFile* output = new TFile(output_name_dir,"RECREATE");

    TTree* outputTree = new TTree("miniTree", "Tree with structured event data");
    Event event;
    outputTree->Branch("event", &event);

    /*--------------------------------------------------------------------------------------
        - Below are all of the histograms for analysis, they are labelled accordingly
        - "_after" refers the histogram being generated after the event cut, etc
    ----------------------------------------------------------------------------------------*/

    // MC events
    TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{MC} [GeV/c]^{2}",100,0,10);
	TH1D* h_y_e = new TH1D("h_y_e",";y_{MC}",100,0.01,0.85);
    TH1D* h_Q2_MC_after = new TH1D("h_Q2_MC_after",";Q^{2}_{MC} [GeV/c]^{2}",100,0,10);
	TH1D* h_y_MC_after = new TH1D("h_y_MC_after",";y_{MC}",100,0.01,0.85);
	TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} [GeV]",100,0,20);
    TH1D* h_energy_MC_after = new TH1D("h_energy_MC_after",";E_{MC} [GeV]",100,0,20);
	TH1D* h_e_pt_MC = new TH1D("h_e_pt_MC",";p_{T,e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_pt_MC_after = new TH1D("h_e_pt_MC_after",";p_{T,e,MC} [GeV/c]",200,0,20);
	TH1D* h_e_pz_MC = new TH1D("h_e_pz_MC",";p_{z,e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_pz_MC_after = new TH1D("h_e_pz_MC_after",";p_{z,e,MC} [GeV/c]",200,0,20);
	TH1D* h_e_p_MC = new TH1D("h_e_p_MC",";p_{e,MC} [GeV/c]",200,0,20);
    TH1D* h_e_p_MC_after = new TH1D("h_e_p_MC_after",";p_{e,MC} [GeV/c]",200,0,20);
    TH1D* h_eta_MC = new TH1D("h_eta_MC",";#eta_{MC}",100,-4,4);
    TH1D* h_eta_MC_after = new TH1D("h_eta_MC_after",";#eta_{MC}",100,-4,4);
	TH1D* h_phi_MC = new TH1D("h_phi_MC",";#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_MC_after = new TH1D("h_phi_MC_after",";#phi_{MC}",100,-3.14,3.14);
	TH1D* h_theta_MC = new TH1D("h_theta_MC",";#theta_{MC}",100,0,3.14);
    TH1D* h_theta_MC_after = new TH1D("h_theta_MC_after",";#theta_{MC}",100,0,3.14);
    TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";VM_{MC} mass [GeV/c^{2}]",200,0,2);
    TH1D* h_VM_pt_MC = new TH1D("h_VM_pt_MC",";p_{T,VM,MC} [GeV/c]",200,0,10);
    TH1D* h_VM_pt_MC_after = new TH1D("h_VM_pt_MC_after",";p_{T,VM,MC} [GeV/c]",200,0,10);
	TH1D* h_VM_pz_MC = new TH1D("h_VM_pz_MC",";p_{z,VM,MC} [GeV/c]",200,0,20);
	TH1D* h_VM_p_MC = new TH1D("h_VM_p_MC",";p_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_Epz_MC = new TH1D("h_VM_Epz_MC",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_VM_Epz_MC_after = new TH1D("h_VM_Epz_MC_after",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_Epz_MC = new TH1D("h_Epz_MC", ";(E_{MC} - p_{z,MC}) [GeV]",100,15,25);
    TH1D* h_Epz_MC_after = new TH1D("h_Epz_MC_after", ";(E_{MC} - p_{z,MC}) [GeV]",100,15,25);
    TH2D* h_EvsP_MC = new TH2D("h_EvsP_MC",";|p|_{MC} [GeV]; E_{MC} [GeV]",100,0,20,100,0,20);
    TH1D* h_EoverP_MC = new TH1D("h_EoverP_MC",";E_{MC}/|p|_{MC}",100,0,2);
    TH1D* h_EoverP_MC_after = new TH1D("h_EoverP_MC_after",";E_{MC}/|p|_{MC}",100,0,2);
    TH1D* h_t_MC = new TH1D("h_t_MC",";|t|_{MC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_MC_after = new TH1D("h_t_MC_after",";|t|_{MC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_sigma = new TH1D("h_sigma",":#sigma; counts",100,0,20);

    // RECO track
	TH1D* h_energy_REC_trk = new TH1D("h_energy_REC_trk",";E_{trk} [GeV]",100,0,20);
    TH2D* h_energy_res_trk = new TH2D("h_energy_res_trk",";E_{MC} [GeV]; (E_{MC}-E_{trk})/E_{MC}",100,0,20,1000,-1,1);
	TH1D* h_eta_REC_trk = new TH1D("h_eta_REC_trk",";#eta_{trk}",100,-4,4);
	TH1D* h_e_pt_REC_trk = new TH1D("h_e_pt_REC_trk",";p_{T,e,trk} [GeV/c]",200,0,20);
	TH1D* h_e_pz_REC_trk = new TH1D("h_e_pz_REC_trk",";p_{z,e,trk} [GeV/c]",200,0,20);
	TH1D* h_e_p_REC_trk = new TH1D("h_e_p_REC_trk",";p_{e,trk} [GeV/c]",200,0,20);
    TH1D* h_t_REC_trk_cut = new TH1D("h_t_REC_trk_cut",";|t|_{trk} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_trk = new TH2D("h_t_res_trk",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_trk_cut = new TH2D("h_t_res_trk_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{trk})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_response_trk = new TH2D("h_t_response_trk","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_trk_cut = new TH2D("h_t_response_trk_cut","; |t|_{trk} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH1D* h_trk_position_x_REC = new TH1D("h_trk_position_x_REC",";x_{trk} [mm]",80,-800,800);
    TH1D* h_trk_position_y_REC = new TH1D("h_trk_position_y_REC",";y_{trk} [mm]",80,-800,800);
    TH2D* h_XvsY_trk = new TH2D("h_XvsY_trk",";x [mm]; y [mm]",80,-800,800,80,-800,800);

	// RECO Q2
    TH1D* h_Q2REC_e_EEMC = new TH1D("h_Q2REC_e_EEMC",";Q^{2}_{EEMC} [GeV/c]^{2}",100,0,10);
    TH2D* h_Q2_res = new TH2D("h_Q2_res",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,0,10,1000,-1,1);
    TH1D* h_Q2_res2 = new TH1D("h_Q2_res2",";(Q^{2}_{RECO}-Q^{2}_{MC})/Q^{2}_{MC}",100,-0.2,0.2);
    TH2D* h_Q2_response = new TH2D("h_Q2_response",";Q^{2}_{MC} [GeV/c]^{2};Q^{2}_{EEMC} [GeV/c]^{2}",100,0,10,100,0,10);
	TH1D* h_dQ2overQ2_REC = new TH1D("h_dQ2overQ2_REC",";dQ^{2}_{RECO}/Q^{2}_{MC}",100,-2,2);
    TH2D* h_Q2_migration = new TH2D("h_Q2_migration",";Q^{2}_{MC} [GeV/c]^{2};Q^{2}_{RECO} [GeV/c]^{2}",100,0,10,100,0,10);
    TH1D* h_Q2_res_vs_counts = new TH1D("h_Q2_res_vs_counts",";(Q^{2}_{RECO}-Q^{2}_{MC})/Q^{2}_{MC}",100,0,10);
    
    // RECO Bjorken-x
    TH2D* h_x_res_cut = new TH2D("h_x_res_cut",";x_{MC} ; x_{MC}-x_{RECO}/x_{MC}",100,0,0.5,100,-0.5,0.5);
    TH1D* h_x_res_cut2 = new TH1D("h_x_res_cut2",";x_{RECO}-x_{MC}/x_{MC}",100,-0.2,0.2);
    TH2D* h_x_response_cut = new TH2D("h_x_response_cut",";x_{MC};x_{RECO} ",100,0,0.5,100,0,0.5);
    TH2D* h_x_migration = new TH2D("h_x_migration",";x_{MC};x_{RECO} ",100,0,0.5,100,0,0.5);
    TH1D* h_x_beforeCut = new TH1D("h_x_beforeCut",";x",100,0,0.5);
    TH1D* h_x = new TH1D("h_x",";x_{MC}",100,0,0.5);
    TH1D* h_x_MC_after = new TH1D("h_x_MC_after",";x_{MC}",100,0,0.5);
    TH1D* h_x_REC = new TH1D("h_x_REC",";x_{RECO}",100,0,0.5);
    TH1D* h_x_afterCut = new TH1D("h_x_afterCut",";x_{RECO}",100,0,0.5);
    TH1D* h_x_res_vs_counts = new TH1D("h_x_res_vs_counts","; x_{RECO}-x_{MC}/x_{MC}",100,0,0.5);

    // RECO y
	TH1D* h_yREC_e_EEMC = new TH1D("h_yREC_e_EEMC",";y_{EEMC}",100,0.01,0.85);
    TH2D* h_y_res = new TH2D("h_y_res",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH1D* h_y_res2 = new TH1D("h_y_res2",";(y_{e,REC}-y_{e,MC})/y_{e,MC}",100,-0.2,0.2);
    TH2D* h_y_response = new TH2D("h_y_response"," ; y_{MC};y_{EEMC}",100,0.01,0.85,100,0.01,0.85);
    TH1D* h_dyOvery_REC = new TH1D("h_dyOvery_REC",";dy/y",100,-2,2);
    TH2D* h_y_migration = new TH2D("h_y_migration",";y_{MC};y_{RECO}",100,0.01,0.85,100,0.01,0.85);

    // RECO energy
	TH1D* h_energy_REC_EEMC = new TH1D("h_energy_REC_EEMC",";E_{EEMC} [GeV]",100,0,20);
    TH2D* h_energy_res_EEMC = new TH2D("h_energy_res_EEMC",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_res_EEMC_after = new TH2D("h_energy_res_EEMC_after",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_response_EEMC = new TH2D("h_energy_response_EEMC","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH2D* h_energy_response_EEMC_after = new TH2D("h_energy_response_EEMC_after","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH2D* h_energy_migration = new TH2D("h_energy_migration",";E_{MC} [GeV];E_{RECO} [GeV]",100,0,20,100,0,20);
    TH1D* h_energy_REC_EEMC_after = new TH1D("h_energy_REC_EEMC_after",";E_{EEMC} [GeV]",100,0,20);
    TH1D* h_energy_res_counts = new TH1D("h_energy_res_counts",";(E_{RECO}-E_{MC})/E_{MC}",100,-0.2,0.2);
    TH1D* h_e_energy_beforeCuts = new TH1D("h_e_energy_beforeCuts",";E_{EEMC} [GeV]",100,0,20);

    // RECO eta
	TH1D* h_eta_REC_EEMC = new TH1D("h_eta_REC_EEMC",";#eta_{EEMC}",100,-4,4);
    TH1D* h_eta_REC_EEMC_after = new TH1D("h_eta_REC_EEMC_after",";#eta_{EEMC}",100,-4,4);
    TH1D* h_eta_diff = new TH1D("h_eta_diff",";#eta_{EEMC}-#eta_{MC}",100,-4,4);
    TH1D* h_eta_diff_after = new TH1D("h_eta_diff_after",";#eta_{EEMC}-#eta_{MC}",100,-4,4);

    // RECO phi
    TH1D* h_phi_REC_EEMC = new TH1D("h_phi_REC_EEMC",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_phi_REC_EEMC_after = new TH1D("h_phi_REC_EEMC_after",";#phi_{EEMC}",100,-3.14,3.14);
    TH1D* h_phi_diff = new TH1D("h_phi_diff",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_cut = new TH1D("h_phi_diff_cut",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH1D* h_phi_diff_after = new TH1D("h_phi_diff_after",";#phi_{EEMC}-#phi_{MC}",100,-3.14,3.14);
    TH2D* h_phi_resolution = new TH2D("h_phi_resolution",";#phi_{MC}; (#phi_{RECO}-#phi_{MC})/#phi_{MC}",100,-3.14,3.14,1000,-1,1);
    TH1D* h_phi_res_diff = new TH1D("h_phi_res_diff",";#phi_{RECO}-#phi_{MC}",100,-3.14,3.14);

    // RECO theta
    TH1D* h_theta_REC_EEMC = new TH1D("h_theta_REC_EEMC",";#theta_{EEMC}",100,0,3.14);
    TH1D* h_theta_REC_EEMC_after = new TH1D("h_theta_REC_EEMC_after",";#theta_{EEMC}",100,-0,3.14);
    TH2D* h_theta_response_EEMC = new TH2D("h_theta_response_EEMC","; #theta_{EEMC}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_response_EEMC_after = new TH2D("h_theta_response_EEMC_after","; #theta_{EEMC}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH1D* h_theta_diff = new TH1D("h_theta_diff",";#theta_{EEMC}-#theta_{MC}",100,0,3.14);
    TH1D* h_theta_diff_after = new TH1D("h_theta_diff_after",";#theta_{EEMC}-#theta_{MC}",100,0,3.14);
    TH2D* h_theta_migration = new TH2D("h_theta_migration","; #theta_{MC};#theta_{RECO}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_resolution = new TH2D("h_theta_resolution",";#theta_{MC}; (#theta_{MC}-#theta_{RECO})/#theta_{MC}",100,0,3.14,1000,-1,1);
    TH1D* h_theta_res_counts = new TH1D("h_theta_res_counts",";(#theta_{RECO}-#theta_{MC})/#theta_{MC}",100,-0.2,0.2);

    // RECO e' P
	TH1D* h_e_pt_REC_EEMC = new TH1D("h_e_pt_REC_EEMC",";p_{T,e,EEMC} [GeV/c]",200,0,20);
    TH2D* h_e_pt_res = new TH2D("h_e_pt_res",";p_{T,e,MC} [GeV/c]; (p_{T,e,EEMC}-p_{T,e,MC})/p_{T,e,MC}",200,0,20,2000,-1,1);
	TH1D* h_e_pz_REC_EEMC = new TH1D("h_e_pz_REC_EEMC",";p_{z,e,EEMC} [GeV/c]",200,0,20);
    TH2D* h_e_pz_res = new TH2D("h_e_pz_res",";p_{z,e,MC} [GeV/c]; (p_z,e,EEMC}-p_{z,e,MC})/p_{z,e,MC}",200,0,20,2000,-1,1);
    TH2D* h_e_pz_response = new TH2D("h_e_pz_response","; p_{z,e,EEMC} [GeV/c];p_{z,e,MC} [GeV/c]",200,0,20,200,0,2);
	TH1D* h_e_p_REC_EEMC = new TH1D("h_e_p_REC_EEMC",";p_{e,EEMC} [GeV/c]",200,0,20);
	TH2D* h_e_p_res = new TH2D("h_e_p_res",";p_{e,MC} [GeV/c]; (p_{e,EEMC}-p_{e,MC})/p_{e,MC}",200,0,20,2000,-1,1);
    
    // RECO VM 
    TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";VM_{RECO} mass [GeV/c^{2}]",200,0,2);
    TH1D* h_VM_p_REC = new TH1D("h_VM_p_REC",";p_{VM,RECO} [GeV/c]",200,0,20);
    TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T,VM,RECO} [GeV/c]",200,0,10);
    TH2D* h_VM_pt_response = new TH2D("h_VM_pt_response","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,10,200,0,10);
	TH1D* h_VM_pz_REC = new TH1D("h_VM_pz_REC",";p_{z,VM,RECO} [GeV/c]",200,0,20);
    TH2D* h_VM_Epz_res = new TH2D("h_VM_Epz_res",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{MC}-(E-p_{z})_{RECO})/(E-p_{z})_{MC}",200,0,20,200,-1,1);
    TH2D* h_VM_Epz_response = new TH2D("h_VM_Epz_response","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",200,0,20,200,0,20);
    TH1D* h_VM_Epz_REC = new TH1D("h_VM_Epz_REC",";(E_{VM,REC}-p_{z,VM,REC}) [GeV]",200,0,20);
    TH2D* h_VM_Epz_migration = new TH2D("h_VM_Epz_migration",";(E_{VM,MC}-p_{z,VM,MC}) [GeV];(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]",200,0,20,200,0,20);
    TH2D* h_VM_pt_migration = new TH2D("h_VM_pt_migration",";p_{T,VM,MC} [GeV/c];p_{T,VM,RECO} [GeV/c]",200,0,10,200,0,10);
    TH2D* h_VM_pt_resolution = new TH2D("h_VM_pt_resolution",";p_{T,MC} [GeV/c]; (p_{T,MC}-p_{T,RECO})/p_{T,MC}",200,0,10,1000,-1,1);
    TH2D* h_VM_pz_resolution = new TH2D("h_VM_pz_resolution",";p_{z,MC} [GeV/c]; (p_{z,MC}-p_{z,RECO})/p_{z,MC}",200,0,20,1000,-1,1);
    TH2D* h_VM_p_resolution = new TH2D("h_VM_p_resolution",";p_{MC} [GeV/c]; (p_{MC}-p_{RECO})/p_{MC}",200,0,20,1000,-1,1);
    TH1D* h_VM_theta_MC = new TH1D("h_VM_theta_MC",";#theta_{VM MC}",100,0,3.14);
    TH1D* h_VM_theta_REC = new TH1D("h_VM_theta_REC",";#theta_{VM RECO}",100,0,3.14);
    TH1D* h_VM_phi_MC = new TH1D("h_VM_phi_MC",";#phi_{VM MC}",100,0,2*M_PI);
    TH1D* h_VM_phi_REC = new TH1D("h_VM_phi_REC",";#phi_{VM RECO}",100,0,2*M_PI);
    TH1D* h_VM_eta_REC = new TH1D("h_VM_eta_REC",";#eta_{VM}",100,-4,4);
    TH1D* h_VM_pt_res_counts = new TH1D("h_VM_pt_res_counts",";(p_{T,RECO}-p_{T,MC})/p_{T,MC}",100,-0.2,0.2);
    TH1D* h_VM_Epz_res_counts = new TH1D("h_VM_Epz_res_counts",";((E-p_{z})_{RECO}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",100,-0.2,0.2);
    TH1D* h_VM_daughterPlus_mass_MC = new TH1D("h_VM_daughterPlus_mass_MC",";VM_{MC} daughter mass [GeV/c^{2}]",100,0,1);
    TH1D* h_VM_daughterMinus_mass_MC = new TH1D("h_VM_daughterMinus_mass_MC",";VM_{MC} daughter mass [GeV/c^{2}]",100,0,1);
    TH1D* h_VM_daughterPlus_mass_REC = new TH1D("h_VM_daughterPlus_mass_REC",";VM_{RECO} daughter mass [GeV/c^{2}]",100,0,1);
    TH1D* h_VM_daughterMinus_mass_REC = new TH1D("h_VM_daughterMinus_mass_REC",";VM_{RECO} daughter mass [GeV/c^{2}]",100,0,1);
    TH1D* h_VM_daughterPlus_mass_REC_before = new TH1D("h_VM_daughterPlus_mass_REC_before",";VM_{RECO} daughter mass [GeV/c^{2}]",100,0,1);
    TH1D* h_VM_daughterMinus_mass_REC_before = new TH1D("h_VM_daughterMinus_mass_REC_before",";VM_{RECO} daughter mass [GeV/c^{2}]",100,0,1);
        
    // RECO position
	TH1D* h_emHits_position_x_REC = new TH1D("h_emHits_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emHits_position_y_REC = new TH1D("h_emHits_position_y_REC","y [mm]",80,-800,800);
	TH1D* h_emClus_position_x_REC = new TH1D("h_emClus_position_x_REC","x [mm]",80,-800,800);
	TH1D* h_emClus_position_y_REC = new TH1D("h_emClus_position_y_REC","y [mm]",80,-800,800);
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (mm);y (mm)",80,-800,800,80,-800,800);
    TH2D* h_emClus_position_REC_cut = new TH2D("h_emClus_position_REC_cut",";x (mm);y (mm)",80,-800,800,80,-800,800);
    TH2D* h_XvsY_hits = new TH2D("h_XvsY_hits",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH2D* h_XvsY_clus = new TH2D("h_XvsY_clus",";x [mm]; y [mm]",80,-800,800,80,-800,800);
    TH1D* h_Xclus_minus_Xtrk = new TH1D("h_Xclus_minus_Xtrk",";x_{clus}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yclus_minus_Ytrk = new TH1D("h_Yclus_minus_Ytrk",";y_{clus}-y_{trk} [mm]",80,-800,800);
    TH1D* h_Xhits_minus_Xtrk = new TH1D("h_Xhits_minus_Xtrk",";x_{hits}-x_{trk} [mm]",80,-800,800);
    TH1D* h_Yhits_minus_Ytrk = new TH1D("h_Yhits_minus_Ytrk",";y_{hits}-y_{trk} [mm]",80,-800,800);
    
    // RECO E-pz
    TH1D* h_Epz_REC = new TH1D("h_Epz_REC", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,15,25);
    TH2D* h_Epz_res = new TH2D("h_Epz_res",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{MC}-(E-p_{z})_{RECO})/(E-p_{z})_{MC}",100,15,25,100,-1,1);
    TH2D* h_Epz_response = new TH2D("h_Epz_response", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,15,25,100,15,25);
    TH2D* h_Epz_migration = new TH2D("h_Epz_migration",";(E_{MC}-p_{z,MC}) [GeV];(E_{RECO} - p_{z,RECO}) [GeV]",100,15,25,100,15,25);
    TH1D* h_Epz_res_counts = new TH1D("h_Epz_res_counts",";((E-p_{z})_{RECO}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",100,-0.2,0.2);
    TH1D* h_e_Epz_noHFS = new TH1D("h_e_Epz_noHFS", ";(E_{EEMC} - p_{z,EEMC}) [GeV]",100,15,25);
    
    // RECO E/P
    TH1D* h_EcalOverPtrk = new TH1D("h_EcalOverPtrk",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_EtrkOverPcal = new TH1D("h_EtrkOverPcal",";E_{trk}/|p|_{EEMC}",100,0,2);
    TH2D* h_EoverP_res = new TH2D("h_EoverP_res",";(E/p)_{MC}; ((E/p)_{MC}-(E/p)_{RECO})/(E/p)_{MC}",100,0,2,100,-1,1);
    TH2D* h_EoverP_response = new TH2D("h_EoverP_response","; E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,2,100,0,2);
    TH2D* h_EoverP_migration = new TH2D("h_EoverP_migration",";E_{MC}/|p|_{MC};E_{EEMC}/|p|_{trk}",100,0,2,100,0,2);
    TH1D* h_EoverP_res_counts = new TH1D("h_EoverP_res_counts",";((E/p)_{REC}-(E/p)_{MC})/(E/p)_{MC}",100,-0.2,0.2);

    // rho
    TH1D* h_rho_rec_before_mass_cut = new TH1D("h_rho_rec_before_mass_cut",";m_{#rho,REC} [GeV/c]",200,0,2);
    TH1D* h_rho_rec_after_mass_cut = new TH1D("h_rho_rec_after_mass_cut",";m_{#rho,REC} [GeV/c]",200,0,2);
    TH1D* h_rho_mass = new TH1D("h_rho_mass",";m_{#rho,MC} [GeV/c]",200,0,2);
    
    // RECO t
    TH1D* h_t_REC_EEMC = new TH1D("h_t_REC_EEMC",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_EEMC_cut = new TH1D("h_t_REC_EEMC_cut",";|t|_{EEMC} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_res_EEMC_cut_percent = new TH2D("h_t_res_EEMC_cut_percent",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_EEMC = new TH2D("h_t_res_EEMC",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res = new TH2D("h_t_res",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH1D* h_t_res_counts = new TH1D("h_t_res_counts",";(|t|_{RECO}-|t|_{MC}) [GeV/c]^{2}",100,-0.2,0.2);
    TH2D* h_t_res_EEMC_cut = new TH2D("h_t_res_EEMC_cut",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{EEMC}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH2D* h_t_migration = new TH2D("h_t_migration",";|t|_{RECO} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_cut = new TH2D("h_t_response_cut","; |t|_{RECO} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_response_EEMC_cut = new TH2D("h_t_response_EEMC_cut","; |t|_{EEMC} [GeV/c]^{2};|t|_{MC} [GeV/c]^{2}",100,0,0.2,100,0,0.2);
    TH2D* h_t_res_proj_percent_pi12 = new TH2D("h_t_res_proj_percent_pi12",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO})/|t|_{MC}",100,0,0.2,1000,-10,10);
    TH2D* h_t_res_proj_12 = new TH2D("h_t_res_proj_12",";|t|_{MC} [GeV/c]^{2}; (|t|_{MC}-|t|_{RECO}) [GeV/c]^{2}",100,0,0.2,1000,-10,10);
    TH1D* h_t_REC_new_method = new TH1D("h_t_REC_new_method",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES = new TH1D("h_t_REC_wRES",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut = new TH1D("h_t_REC_wRES_cut",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi2 = new TH1D("h_t_REC_wRES_cut_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi3 = new TH1D("h_t_REC_wRES_cut_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi4 = new TH1D("h_t_REC_wRES_cut_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi6 = new TH1D("h_t_REC_wRES_cut_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi9 = new TH1D("h_t_REC_wRES_cut_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi12 = new TH1D("h_t_REC_wRES_cut_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi16 = new TH1D("h_t_REC_wRES_cut_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi20 = new TH1D("h_t_REC_wRES_cut_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wRES_cut_pi24 = new TH1D("h_t_REC_wRES_cut_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH2D* h_t_REC_2d = new TH2D("h_t_REC_2d",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_cut = new TH2D("h_t_REC_2d_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_cut2 = new TH2D("h_t_REC_2d_cut2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES = new TH2D("h_t_REC_2d_wRES",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut = new TH2D("h_t_REC_2d_wRES_cut",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi2 = new TH2D("h_t_REC_2d_wRES_cut_pi2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi3 = new TH2D("h_t_REC_2d_wRES_cut_pi3",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi4 = new TH2D("h_t_REC_2d_wRES_cut_pi4",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi6 = new TH2D("h_t_REC_2d_wRES_cut_pi6",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi9 = new TH2D("h_t_REC_2d_wRES_cut_pi9",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi12 = new TH2D("h_t_REC_2d_wRES_cut_pi12",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi16 = new TH2D("h_t_REC_2d_wRES_cut_pi16",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi20 = new TH2D("h_t_REC_2d_wRES_cut_pi20",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wRES_cut_pi24 = new TH2D("h_t_REC_2d_wRES_cut_pi24",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT = new TH2D("h_t_REC_2d_wCUT",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi2 = new TH2D("h_t_REC_2d_wCUT_pi2",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi3 = new TH2D("h_t_REC_2d_wCUT_pi3",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi4 = new TH2D("h_t_REC_2d_wCUT_pi4",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi6 = new TH2D("h_t_REC_2d_wCUT_pi6",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi9 = new TH2D("h_t_REC_2d_wCUT_pi9",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi12 = new TH2D("h_t_REC_2d_wCUT_pi12",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi16 = new TH2D("h_t_REC_2d_wCUT_pi16",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi20 = new TH2D("h_t_REC_2d_wCUT_pi20",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    TH2D* h_t_REC_2d_wCUT_pi24 = new TH2D("h_t_REC_2d_wCUT_pi24",";#sqrt{|t|_{x}} [GeV/c]; #sqrt{|t|_{y}} [GeV/c]",100,0,0.2,100,0,0.2);
    
        /*-----------------------------------------------------------------------------------------
            - t distribution with cut only (for testing)
            - distributions should all look the same since no resolution has been added
        ------------------------------------------------------------------------------------------*/
    TH1D* h_t_REC_wCUT_pi2 = new TH1D("h_t_REC_wCUT_pi2",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi3 = new TH1D("h_t_REC_wCUT_pi3",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi4 = new TH1D("h_t_REC_wCUT_pi4",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi6 = new TH1D("h_t_REC_wCUT_pi6",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi9 = new TH1D("h_t_REC_wCUT_pi9",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi12 = new TH1D("h_t_REC_wCUT_pi12",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi16 = new TH1D("h_t_REC_wCUT_pi16",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi20 = new TH1D("h_t_REC_wCUT_pi20",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    TH1D* h_t_REC_wCUT_pi24 = new TH1D("h_t_REC_wCUT_pi24",";|t|_{RECO} [GeV/c]^{2}; counts",100,0,0.2);
    
    // Compare cuts and extra histograms for analysis
    TH1D* h_Epz_beforeCut = new TH1D("h_Epz_beforeCut", ";(E_{EEMC} - p_{z,trk}) [GeV]",100,15,25);
    TH1D* h_EoverP_afterCut = new TH1D("h_EoverP_afterCut",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_EoverP_beforeCut = new TH1D("h_EoverP_beforeCut",";E_{EEMC}/|p|_{trk}",100,0,2);
    TH1D* h_Q2_beforeCut = new TH1D("h_Q2_beforeCut",";Q^{2} [GeV/c]^{2}",100,0,10);
    TH1D* h_Q2_afterCut = new TH1D("h_Q2_afterCut",";Q^{2} [GeV/c]^{2}",100,0,10);
    TH1D* h_y_beforeCut = new TH1D("h_y_beforeCut",";y",100,0.01,0.85);
    TH1D* h_y_afterCut = new TH1D("h_y_afterCut",";y",100,0.01,0.85);
    TH1D* h_VM_Epz_MC_before = new TH1D("h_VM_Epz_MC_before",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_VM_pt_MC_before = new TH1D("h_VM_pt_MC_before",";p_{T,VM,MC} [GeV/c]",200,0,10);
    TH1D* h_VM_mass_MC_before = new TH1D("h_VM_mass_MC_before",";m_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_pz_MC_before = new TH1D("h_VM_pz_MC_before",";p_{z,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_p_MC_before = new TH1D("h_VM_p_MC_before",";|p|_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_mass_REC_before = new TH1D("h_VM_mass_REC_before",";m_{VM,MC} [GeV/c]",200,0,2);
    TH1D* h_VM_pt_REC_before = new TH1D("h_VM_pt_REC_before",";p_{T,VM,MC} [GeV/c]",200,0,10);
    TH1D* h_VM_pz_REC_before = new TH1D("h_VM_pz_REC_before",";p_{z,VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_p_REC_before = new TH1D("h_VM_p_REC_before",";|p|_{VM,MC} [GeV/c]",200,0,20);
    TH1D* h_VM_Epz_REC_before = new TH1D("h_VM_Epz_REC_before",";(E_{VM,MC}-p_{z,VM,MC}) [GeV]",200,0,20);
    TH1D* h_solid_angle_REC = new TH1D("h_solid_angle_REC",";#Omega between #theta_{e'} and #phi_{e'}",1000,-0.1,0.1);
    TH1D* h_angle_resolution = new TH1D("h_angle_resolution",";Angular Resolution",1000,-0.1,0.1);
    TH2D* h_phi_response = new TH2D("h_phi_response",";#phi_{RECO}; #phi_{MC}",100,-3.14,3.14,100,-3.14,3.14);
    TH2D* h_theta_response = new TH2D("h_theta_response",";#theta_{RECO}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH1D* h_omega_MC = new TH1D("h_omega_MC",";#Omega_{MC}",100,-3.14,3.14);
    TH1D* h_omega_REC = new TH1D("h_omega_REC",";#Omega_{REC}",100,-3.14,3.14);
	TH2D* h_e_pz_trk_res = new TH2D("h_e_pz_trk_res",";p_{z,e,MC} [GeV/c]; (p_{z,e,trk}-p_{z,e,MC})/p_{z,e,MC}",200,0,20,2000,-1,1);
	TH2D* h_e_p_trk_res = new TH2D("h_e_p_trk_res",";|p|_{e,MC} [GeV/c]; (|p|_{e,trk}-|p|_{e,MC})/|p|_{e,MC}",200,0,20,2000,-1,1);
    TH2D* h_Q2_vs_x_REC = new TH2D("h_Q2_vs_x_REC",";x_{RECO}; Q^{2}_{RECO} [GeV/c]^{2}",100,0,0.5,100,0,10);
    TH2D* h_Q2_vs_x_MC = new TH2D("h_Q2_vs_x_MC",";x_{MC}; Q^{2}_{MC} [GeV/c]^{2}",100,0,0.5,100,0,10);
    TH2D* h_Q2_res_beforecut = new TH2D("h_Q2_res_beforecut",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,0,10,1000,-1,1);
    TH2D* h_Q2_response_beforecut = new TH2D("h_Q2_response_beforecut","; Q^{2}_{EEMC} [GeV/c]^{2};Q^{2}_{MC} [GeV/c]^{2}",100,0,10,100,0,10);     
    TH2D* h_EvsP_REC = new TH2D("h_EvsP_REC",";|p|_{trk} [GeV]; E_{EEMC} [GeV]",100,0,20,100,0,20);

    // detectors
    TH1D* hEEMC_hits = new TH1D("hEEMC_hits", "Number of EEMC clusters per event", 50, 0, 50);
    TH1D* hZDC_hits_neutral  = new TH1D("hZDC_hits_neutral", "ZDC cluster multiplicity;clusters/event;Counts", 50, 0, 50);
    TH1D* hZDC_hits  = new TH1D("hZDC_hits", "ZDC cluster multiplicity;clusters/event;Counts", 50,0,50);
    TH1D* h_zdc_raw_multiplicity = new TH1D("h_zdc_raw_multiplicity", "ZDC Raw Hit Multiplicity;N_{raw hits};Events", 50, 0, 50);
    TH1D* h_zdc_Ecal_multiplicity = new TH1D("h_zdc_Ecal_multiplicity", "ZDC Ecal Hit Multiplicity;N_{raw hits};Events", 50, 0, 50);
    TH1D* h_zdc_Hcal_multiplicity = new TH1D("h_zdc_Hcal_multiplicity", "ZDC Hcal Hit Multiplicity;N_{raw hits};Events", 50, 0, 50);
    TH1D* hRP_hits   = new TH1D("hRP_hits",   "Number of RP hits per event", 50, 0, 50);
    TH1D* hOMD_hits  = new TH1D("hOMD_hits",  "Number of OMD hits per event", 50, 0, 50);
    TH1D* hEEMC_energy = new TH1D("hEEMC_energy", "EEMC Cluster Energy", 50,0,30);
    TH1D* hZDC_energy = new TH1D("hZDC_energy", "ZDC Energy;E [GeV];Counts", 200, 0, 200);
    TH2D* hEEMC_xy = new TH2D("hEEMC_xy", "EEMC hit map;x [mm];y [mm]", 100, -200, 200, 100, -200, 200);
    TH2D* hZDC_xy = new TH2D("hZDC_xy", "ZDC XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hRP_xy = new TH2D("hRP_xy", "RP XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hOMD_xy = new TH2D("hOMD_xy", "OMD XY;x [mm];y [mm]", 200, -1200, -600, 200, -300, 300);
    TH2D* hZDC_hits_xy = new TH2D("hZDC_hits_xy", "ZDC Hit Positions;x [mm];y [mm]",200, -1200, 600, 200, -300, 300);
    TH1D* h_zdc_raw_amplitude = new TH1D("h_zdc_raw_amplitude", "ZDC Raw Hit Amplitude;ADC;Counts", 200, 0, 2000);

    TH1F* h_Nevents_MC = new TH1F("h_Nevents_MC","Total MC Events; ;Count",1, 0, 1);
    TH1F* h_Nevents_inAcceptance = new TH1F("h_Nevents_inAcceptance","Events with >=1 Forward MC Particle; ;Count", 1, 0, 1);
    TH1F* h_Nevents_ZDC = new TH1F("h_Nevents_ZDC", "Events with >=1 ZDC Cluster; ;Count",1, 0, 1);

    TH1F* h_eta_neutrons = new TH1F("h_eta_neutrons","MC Neutron #eta;#eta;Counts",200, -10, 10);
    TH1F* h_multiplicity_neutrons = new TH1F("h_multiplicity_neutrons","MC Neutron Multiplicity;N_{raw hits};Events", 50, 0, 50);
    TH1F* h_multiplicity_neutrons_accept = new TH1F("h_multiplicity_neutrons_accept","MC Neutron Multiplicity;N_{raw hits};Events", 50, 0, 50);
    TH1F* h_eta_ZDC_accept = new TH1F("h_eta_ZDC_accept","ZDC Acceptance;#eta;Counts",200, -10, 10);
    TH1F* h_eta_ZDC_reco = new TH1F("h_eta_ZDC_reco","ZDC Clusters;#eta;Counts",200, -10, 10);
    TH1F* hZDC_hits_eta = new TH1F("hZDC_hits_eta","ZDC Position;#eta;Counts",200, -10, 10);

    TH1D* h_Ntotal = new TH1D("h_Ntotal", "Total Events; ;Count",1, 0, 1);
    TH1D* h_Nvetoed = new TH1D("h_Nvetoed", "Events Vetoed; ;Count",1, 0, 1);
    TH1D* h_NafterVeto = new TH1D("h_NafterVeto", "Events After Veto; ;Count",1, 0, 1);

    // Diagnostic before/after veto histograms 
    TH1D* h_t_incoh_before      = new TH1D("h_t_incoh_before",      "Incoherent t before veto; t; Events", 200, 0, 1);
    TH1D* h_t_incoh_reco      = new TH1D("h_t_incoh_reco",      "Incoherent t reco after vetoes; t; Events", 200, 0, 1);
    TH1D* h_t_incoh_after_ZDC   = new TH1D("h_t_incoh_after_ZDC",   "Incoherent t after ZDC veto; t; Events", 200, 0, 1);
    TH1D* h_t_incoh_after_all   = new TH1D("h_t_incoh_after_all",   "Incoherent t after ALL vetoes; t; Events", 200, 0, 1);

    // Diagnostic distributions
    TH1D* h_ZDC_energy_lowt     = new TH1D("h_ZDC_energy_lowt",     "ZDC energy for low-t incoherent events", 200, 0, 20);
    TH1D* h_ZDC_energy_lowt_rec     = new TH1D("h_ZDC_energy_lowt_rec",     "ZDC energy for low-t reco incoherent events", 200, 0, 20);
    TH1D* h_forward_eta_lowt    = new TH1D("h_forward_eta_lowt",    "Forward eta for low-t incoherent events", 200, -5, 8);
    TH1D* h_forward_eta_lowt_rec    = new TH1D("h_forward_eta_lowt_rec",    "Forward eta for low-t reco incoherent events", 200, -5, 8);
    TH1D* h_HFS_lowt            = new TH1D("h_HFS_lowt",            "HFS magnitude for low-t incoherent events", 200, 0, 20);
    TH1D* h_HFS_lowt_rec            = new TH1D("h_HFS_lowt_rec",            "HFS magnitude for low-t reco incoherent events", 200, 0, 20);
    TH1D* h_eff_ZDC  = (TH1D*) h_t_incoh_after_ZDC->Clone("h_eff_ZDC");
    TH1D* h_eff_all  = (TH1D*) h_t_incoh_after_all->Clone("h_eff_all");

    TH1D* h_zdc_Ecal_eta = new TH1D("h_zdc_Ecal_eta", "ZDC Ecal; #eta",200,-10,10); 
    TH1D* h_zdc_Hcal_eta = new TH1D("h_zdc_Hcal_eta", "ZDC Hcal; #eta",200,-10,10); 
    TH1D* h_zdc_REChits_eta = new TH1D("h_zdc_REChits_eta", "ZDC Rec Hits; #eta",200,-10,10); 
    TH1D* h_zdc_neutrals_eta = new TH1D("h_zdc_neutrals_eta", "ZDC Neutrals; #eta",200,-10,10); 
    TH1D* h_zdc_Ecal_theta = new TH1D("h_zdc_Ecal_theta", "ZDC Ecal; #theta",200,-10,10); 
    TH1D* h_zdc_Hcal_theta = new TH1D("h_zdc_Hcal_theta", "ZDC Hcal; #theta",200,-10,10); 
    TH1D* h_zdc_theta = new TH1D("h_zdc_theta", "ZDC Neutrals; #theta",200,-10,10); 
    TH1D* h_rp_eta = new TH1D("h_rp_eta", "RP; #eta",200,-10,10); 
    TH1D* h_rp_theta = new TH1D("h_rp_theta", "RP; #theta",200,-10,10); 
    TH1D* h_omd_eta = new TH1D("h_omd_eta", "OMD; #eta",200,-10,10); 
    TH1D* h_omd_theta = new TH1D("h_omd_theta", "OMD; #theta",200,-10,10); 
    TH1D* h_EEMC_eta = new TH1D("h_EEMC_eta", "EEMC; #eta",200,-10,10); 
    TH1D* h_EEMC_theta = new TH1D("h_EEMC_theta", "EEMC; #theta",200,-10,10); 

    TH1D* h_eta_REC_had_in_ZDC = new TH1D("h_eta_REC_had_in_ZDC","",200,0,5);
    TH1D* h_eta_REC_had_in_RP = new TH1D("h_eta_REC_had_in_RP","",200,0,5);
    TH1D* h_eta_REC_had_in_OMD = new TH1D("h_eta_REC_had_in_OMD","",200,0,5);
    TH1D* h_eta_REC_had_in_ZDC_trk = new TH1D("h_eta_REC_had_in_ZDC_trk","",200,0,5);
    TH1D* h_eta_REC_had_in_RP_trk = new TH1D("h_eta_REC_had_in_RP_trk","",200,0,5);
    TH1D* h_eta_REC_had_in_OMD_trk = new TH1D("h_eta_REC_had_in_OMD_trk","",200,0,5);

    TH1D* h_eta_REC_had_in_ZDC_e = new TH1D("h_eta_REC_had_in_ZDC_e","",200,0,5);
    TH1D* h_eta_REC_had_in_RP_e = new TH1D("h_eta_REC_had_in_RP_e","",200,0,5);
    TH1D* h_eta_REC_had_in_OMD_e = new TH1D("h_eta_REC_had_in_OMD_e","",200,0,5);

    TH1D* h_t_pi12_notVetoed_hits = new TH1D("h_t_pi12_notVetoed_hits","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_hits = new TH1D("h_t_pi12_vetoed_hits","",100,0,0.2);
    TH1D* h_t_pi12_notVetoed_neut = new TH1D("h_t_pi12_notVetoed_neut","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_neut = new TH1D("h_t_pi12_vetoed_neut","",100,0,0.2);
    TH1D* h_t_pi12_notVetoed_clust = new TH1D("h_t_pi12_notVetoed_clust","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_clust = new TH1D("h_t_pi12_vetoed_clust","",100,0,0.2);
    TH1D* h_t_pi12_notVetoed_RP = new TH1D("h_t_pi12_notVetoed_RP","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_RP = new TH1D("h_t_pi12_vetoed_RP","",100,0,0.2);
    TH1D* h_t_pi12_notVetoed_OMD = new TH1D("h_t_pi12_notVetoed_OMD","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_OMD = new TH1D("h_t_pi12_vetoed_OMD","",100,0,0.2);
    TH1D* h_t_pi12_notVetoed_etaHFS = new TH1D("h_t_pi12_notVetoed_etaHFS","",100,0,0.2);
    TH1D* h_t_pi12_vetoed_etaHFS = new TH1D("h_t_pi12_vetoed_etaHFS","",100,0,0.2);

    TH1I* h_total_hfs_particles = new TH1I("h_total_hfs_particles","HFS multiplicity per event;N_{HFS};Counts",4, 0, 4);
    TH2D* h_MC_neutrons_vs_t_2D = new TH2D("h_MC_neutrons_vs_t_2D",";neutrons; |t|",200,0,5,200,0,0.2);
    TH2D* h_MC_neutrons_vs_t_2Dv2 = new TH2D("h_MC_neutrons_vs_t_2Dv2",";neutrons; |t|",200,0,5,200,0,0.2);
    TH2D* h_REC_neutrons_vs_t_2D = new TH2D("h_REC_neutrons_vs_t_2D",";neutrons_{reco}; neutrons_{MC}",200,0,5,200,0,5);

    TH1D* h_mult_neut_pi12 = new TH1D("h_mult_neut_pi12","REC Neutron Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_zdc_hits_pi12 = new TH1D("h_mult_zdc_hits_pi12","ZDC REC Hits Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_zdc_clusters_pi12 = new TH1D("h_mult_zdc_clusters_pi12","ZDC Cluster Multiplicity in #theta_{max} = #pi/12;N_{clusters};Events", 50, 0, 50);
    TH1D* h_mult_rp_hits_pi12 = new TH1D("h_mult_rp_hits_pi12","RP Hit Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_omd_hits_pi12 = new TH1D("h_mult_omd_hits_pi12","OMD Hit Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_eemc_clusters_pi12 = new TH1D("h_mult_eemc_clusters_pi12","EEMC Cluster Multiplicity in #theta_{max} = #pi/12;N_{clusters};Events", 50, 0, 50);
    TH1D* h_mult_hfs_pi12 = new TH1D("h_mult_hfs_pi12","HFS Particle Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_eta_pi12 = new TH1D("h_mult_eta_pi12","#eta Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_eta_plus_hfs_pi12 = new TH1D("h_mult_eta_plus_hfs_pi12","#eta + HFS Particle Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_electron_pi12 = new TH1D("h_mult_electron_pi12","Electron Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_kaon_pi12 = new TH1D("h_mult_kaon_pi12","Kaons Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1D* h_mult_HFSplusKaon = new TH1D("h_mult_HFSplusKaon","HFS+kaon Particle Multiplicity in #theta_{max} = #pi/12;N_{hits};Events", 50, 0, 50);
    TH1I* h_pass_pi12_flag = new TH1I("h_pass_pi12_flag","Pass flag for #theta_{max} = #pi/12;flag;counts",2, 0, 2);

    TH1I* h_Nevents = new TH1I("h_Nevents", "Total Number of Events", 1, 0, 1);
    TH1I* h_passMass = new TH1I("h_passMass", "Events passing mass selection", 1, 0, 1);
    TH1I* h_passZDC = new TH1I("h_passZDC", "Events passing ZDC veto", 1, 0, 1);
    TH1I* h_passRP = new TH1I("h_passRP", "Events passing RP veto", 1, 0, 1);
    TH1I* h_passOMD = new TH1I("h_passOMD", "Events passing OMD veto", 1, 0, 1);
    TH1I* h_passEtaHFS = new TH1I("h_passEtaHFS", "Events passing #eta+HFS veto", 1, 0, 1);
    TH1I* h_passEventSelection = new TH1I("h_passEventSelection", "Events passing event cuts", 1, 0, 1);
    TH1I* h_passPi12 = new TH1I("h_passPi12", "Events passing #theta_{max} cut", 1, 0, 1);
    TH1I* h_pass_e = new TH1I("h_pass_e", "Events with electron", 1, 0, 1);
    TH1I* h_passHFS = new TH1I("h_passHFS", "Events passing HFS+kaon criteria", 1, 0, 1);
    TH1I* h_totalHFS = new TH1I("h_totalHFS", "Events with 2 HFS particles", 1, 0, 1);
    TH1I* h_nKaons = new TH1I("h_nKaons", "Events with 2 kaons", 1, 0, 1);
    TH1I* h_didntPass_e = new TH1I("h_didntPass_e", "Events Didn't Pass e' Requirement", 1, 0, 1);
    TH1I* h_didntPasstotalHFS = new TH1I("h_didntPasstotalHFS", "Events Didn't Pass HFS==2 Requirement", 1, 0, 1);
    TH1I* h_didntPassnKaons = new TH1I("h_didntPassnKaons", "Events Didn't Pass Kaon==2 Requirement", 1, 0, 1);
    TH1I* h_didntPassHFS = new TH1I("h_didntPassHFS", "Events Didn't Pass HFS==2 plus Kaon==2 Requirement", 1, 0, 1);
    TH1I* h_didntPassEtaHFS = new TH1I("h_didntPassEtaHFS", "Events Didn't Pass #eta Requirement", 1, 0, 1);
    TH1I* h_didntPassMass = new TH1I("h_didntPassMass", "Events Didn't Pass VM Mass Requirement", 1, 0, 1);
    TH1I* h_didntPassZDC = new TH1I("h_didntPassZDC", "Events Didn't Pass ZDC veto", 1, 0, 1);
    TH1I* h_didntPassRP = new TH1I("h_didntPassRP", "Events Didn't Pass RP veto", 1, 0, 1);
    TH1I* h_didntPassOMD = new TH1I("h_didntPassOMD", "Events Didn't Pass OMD veto", 1, 0, 1);

    TH1D* h_t_pi12_afterElectronSelect = new TH1D("h_t_pi12_afterElectronSelect","",100,0,0.2);
    TH1D* h_t_pi12_afterHFS = new TH1D("h_t_pi12_afterHFS","",100,0,0.2);
    TH1D* h_t_pi12_afterKaons = new TH1D("h_t_pi12_afterKaons","",100,0,0.2);
    TH1D* h_t_pi12_afterEta = new TH1D("h_t_pi12_afterEta","",100,0,0.2);
    TH1D* h_t_pi12_afterMass = new TH1D("h_t_pi12_afterMass","",100,0,0.2);
    TH1D* h_t_pi12_afterQ2y = new TH1D("h_t_pi12_afterQ2y","",100,0,0.2);
    TH1D* h_t_pi12_afterZDC = new TH1D("h_t_pi12_afterZDC","",100,0,0.2);
    TH1D* h_t_pi12_afterRP = new TH1D("h_t_pi12_afterRP","",100,0,0.2);
    TH1D* h_t_pi12_afterOMD = new TH1D("h_t_pi12_afterOMD","",100,0,0.2);
    TH1D* h_t_after_VMrec = new TH1D("h_t_after_VMrec","",100,0,0.2);
    TH1D* h_t_after_2tracks = new TH1D("h_t_after_2tracks","",100,0,0.2);

    TH1D* h_theta_L = new TH1D("h_theta_L",";#theta_{max} [rad];Counts", 60, 0, M_PI/2);
    TH1D* h_theta_wedge = new TH1D("h_theta_wedge",";#theta_{max} [rad];Counts", 60, 0, M_PI/2);
    TH1D* h_theta_all = new TH1D("h_theta_all",";#theta_{max} [rad];Counts", 60, 0, M_PI/2);

    TH2D* h_kPlus_p_resolution = new TH2D("h_kPlus_p_resolution",";p_{MC} [GeV/c]; (p_{MC}-p_{RECO})/p_{MC}",100,0,20,1000,-1,1);
    TH2D* h_kMinus_p_resolution = new TH2D("h_kMinus_p_resolution",";p_{MC} [GeV/c]; (p_{MC}-p_{RECO})/p_{MC}",100,0,20,1000,-1,1);
    TH1D* h_kPlus_p = new TH1D("h_kPlus_p",";p_{K^{+}} [GeV/c]",100,0,10);
    TH1D* h_kMinus_p = new TH1D("h_kMinus_p",";p_{K^{-}} [GeV/c];",100,0,10);
    TH1D* h_kMinus_phi = new TH1D("h_kMinus_phi",";#phi_{K^{-}};",100,0,2*M_PI);
    TH1D* h_kPlus_phi = new TH1D("h_kPlus_phi",";#phi_{K^{+}};",100,0,2*M_PI);
    TH1D* h_kMinus_theta = new TH1D("h_kMinus_theta",";#theta_{K^{-}};",100,0,M_PI);
    TH1D* h_kPlus_theta = new TH1D("h_kPlus_theta",";#theta_{K^{+}};",100,0,M_PI);
    TH1D* h_t_phi = new TH1D("h_t_phi",";#phi_{|t|};",100,0,2*M_PI);
    TH1D* h_q_phi = new TH1D("h_q_phi",";#phi_{q};",100,0,2*M_PI);
    TH1D* h_t_theta = new TH1D("h_t_theta",";#theta_{|t|};",100,0,M_PI);
    TH1D* h_q_theta = new TH1D("h_q_theta",";#theta_{q};",100,0,M_PI);
    TH1D* h_t_perp = new TH1D("h_t_perp",";|t|_{perp} [GeV/c]^{2};",100,0,0.2);
    TH1D* h_q_perp = new TH1D("h_q_perp",";q_{perp} [GeV/c];",100,0,0.5);
    TH1D* h_q_x = new TH1D("h_q_x",";q_{x} [GeV/c];",100,0,0.5);
    TH1D* h_q_y = new TH1D("h_q_y",";q_{y} [GeV/c];",100,0,0.5);
    TH1D* h_nHat_phi = new TH1D("h_nHat_phi",";#phi_{#hat{n}};",100,0,2*M_PI);
    TH1D* h_vm_px = new TH1D("h_vm_px",";p_{VM,x};",100,0,2);
    TH1D* h_vm_py = new TH1D("h_vm_py",";p_{VM,y};",100,0,2);
    TH1D* h_vm_pz = new TH1D("h_vm_pz",";p_{VM,z};",100,0,2);
    TH1D* h_qx_phi = new TH1D("h_qx_phi",";#phi_{q_{x}};",100,0,2*M_PI);
    TH1D* h_qy_phi = new TH1D("h_qy_phi",";#phi_{q_{y}};",100,0,2*M_PI);
    TH1D* h_qx_theta = new TH1D("h_qx_theta",";#theta_{q_{x}};",100,0,M_PI);
    TH1D* h_qy_theta = new TH1D("h_qy_theta",";#theta_{q_{y}};",100,0,M_PI);
    TH1D* h_q_mag = new TH1D("h_q_mag","",100,0,0.5);

    TH1D* h_res_t = new TH1D("h_res_t","",100,0,0.2);
    TH1D* h_res_q = new TH1D("h_res_q","",100,0,0.2);
    TH1D* h_res_tperp = new TH1D("h_res_tperp","",100,0,0.2);
    TH1D* h_res_qperp = new TH1D("h_res_qperp","",100,0,0.2);
    TH1D* h_res_eScat = new TH1D("h_res_eScat","",100,0,0.2);
    TH1D* h_res_vm = new TH1D("h_res_vm","",100,0,0.2);
    TH1D* h_res_vm_phi = new TH1D("h_res_vm_phi","",100,0,2*M_PI);
    TH1D* h_res_vm_theta = new TH1D("h_res_vm_theta","",100,0,M_PI);
    TH1D* h_res_vm_mag = new TH1D("h_res_vm_mag","",100,0,10);
    TH1D* h_res_omega = new TH1D("h_res_omega","",100,0,2*M_PI);
    TH1D* h_res_nHat = new TH1D("h_res_nHat","",100,0,0.2);
    TH1D* h_res_qy = new TH1D("h_res_qy","",100,0,0.2);

    TH2D* h_Q2_res_beforeCuts = new TH2D("h_Q2_res_beforeCuts",";Q^{2}_{MC} [GeV/c]^{2}; (Q^{2}_{MC}-Q^{2}_{EEMC})/Q^{2}_{MC}",100,0,10,1000,-1,1);
    TH1D* h_Q2_res2_beforeCuts = new TH1D("h_Q2_res2_beforeCuts",";(Q^{2}_{RECO}-Q^{2}_{MC})/Q^{2}_{MC}",100,-0.2,0.2);
    TH2D* h_Q2_response_beforeCuts = new TH2D("h_Q2_response_beforeCuts",";Q^{2}_{MC} [GeV/c]^{2};Q^{2}_{EEMC} [GeV/c]^{2}",100,0,10,100,0,10);
    TH2D* h_Q2_migration_beforeCuts = new TH2D("h_Q2_migration_beforeCuts",";Q^{2}_{MC} [GeV/c]^{2};Q^{2}_{RECO} [GeV/c]^{2}",100,0,10,100,0,10);
    TH2D* h_Q2_vs_x_REC_beforeCuts = new TH2D("h_Q2_vs_x_REC_beforeCuts",";x_{RECO}; Q^{2}_{RECO} [GeV/c]^{2}",100,0,0.5,100,0,10);
    TH2D* h_x_res_cut_beforeCuts = new TH2D("h_x_res_cut_beforeCuts",";x_{MC} ; x_{MC}-x_{RECO}/x_{MC}",100,0,0.5,100,-0.5,0.5);
    TH1D* h_x_res_cut2_beforeCuts = new TH1D("h_x_res_cut2_beforeCuts",";x_{RECO}-x_{MC}/x_{MC}",100,-0.2,0.2);
    TH2D* h_x_response_cut_beforeCuts = new TH2D("h_x_response_cut_beforeCuts",";x_{MC};x_{RECO} ",100,0,0.5,100,0,0.5);
    TH2D* h_x_migration_beforeCuts = new TH2D("h_x_migration_beforeCuts",";x_{MC};x_{RECO} ",100,0,0.5,100,0,0.5);
    TH2D* h_y_migration_beforeCuts = new TH2D("h_y_migration_beforeCuts",";y_{MC};y_{RECO}",100,0.01,0.85,100,0.01,0.85);
    TH2D* h_y_res_beforeCuts = new TH2D("h_y_res_beforeCuts",";y_{e,MC} ;(y_{e,MC}-y_{e,EEMC})/y_{e,MC}",100,0.01,0.85,1000,-1,1);
    TH1D* h_y_res2_beforeCuts = new TH1D("h_y_res2_beforeCuts",";(y_{e,REC}-y_{e,MC})/y_{e,MC}",100,-0.2,0.2);
    TH2D* h_y_response_beforeCuts = new TH2D("h_y_response_beforeCuts"," ; y_{MC};y_{EEMC}",100,0.01,0.85,100,0.01,0.85);
    TH2D* h_Epz_res_beforeCuts = new TH2D("h_Epz_res_beforeCuts",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{MC}-(E-p_{z})_{RECO})/(E-p_{z})_{MC}",100,15,25,100,-1,1);
    TH2D* h_Epz_response_beforeCuts = new TH2D("h_Epz_response_beforeCuts", ";(E_{EEMC} - p_{z,EEMC}) [GeV];(E_{MC}-p_{z,MC}) [GeV]",100,15,25,100,15,25);
    TH2D* h_Epz_migration_beforeCuts = new TH2D("h_Epz_migration_beforeCuts",";(E_{MC}-p_{z,MC}) [GeV];(E_{RECO} - p_{z,RECO}) [GeV]",100,15,25,100,15,25);
    TH2D* h_Epz_migration_noHFS = new TH2D("h_Epz_migration_noHFS",";(E_{MC}-p_{z,MC}) [GeV];(E_{RECO} - p_{z,RECO}) [GeV]",100,15,25,100,15,25);
    TH1D* h_Epz_res2_beforeCuts = new TH1D("h_Epz_res2_beforeCuts",";((E-p_{z})_{RECO}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",100,-0.2,0.2);
    TH2D* h_EoverP_res_beforeCuts = new TH2D("h_EoverP_res_beforeCuts",";(E/p)_{MC}; ((E/p)_{MC}-(E/p)_{RECO})/(E/p)_{MC}",100,0,2,100,-1,1);
    TH2D* h_EoverP_response_beforeCuts = new TH2D("h_EoverP_response_beforeCuts","; E_{EEMC}/|p|_{trk};E_{MC}/|p|_{MC}",100,0,2,100,0,2);
    TH2D* h_EoverP_migration_beforeCuts = new TH2D("h_EoverP_migration_beforeCuts",";E_{MC}/|p|_{MC};E_{EEMC}/|p|_{trk}",100,0,2,100,0,2);
    TH1D* h_EoverP_res2_beforeCuts = new TH1D("h_EoverP_res2_beforeCuts",";((E/p)_{REC}-(E/p)_{MC})/(E/p)_{MC}",100,-0.2,0.2);
    TH2D* h_energy_res_EEMC_beforeCuts = new TH2D("h_energy_res_EEMC_beforeCuts",";E_{MC} [GeV]; (E_{MC}-E_{EEMC})/E_{MC}",100,0,20,1000,-1,1);
    TH2D* h_energy_response_EEMC_beforeCuts = new TH2D("h_energy_response_EEMC_beforeCuts","; E_{EEMC} [GeV];E_{MC} [GeV]",100,0,20,100,0,20);
    TH2D* h_energy_migration_beforeCuts = new TH2D("h_energy_migration_beforeCuts",";E_{MC} [GeV];E_{RECO} [GeV]",100,0,20,100,0,20);
    TH1D* h_energy_res2_beforeCuts = new TH1D("h_energy_res2_beforeCuts",";(E_{RECO}-E_{MC})/E_{MC}",100,-0.2,0.2);
    TH2D* h_theta_response_EEMC_beforeCuts = new TH2D("h_theta_response_EEMC_beforeCuts","; #theta_{EEMC}; #theta_{MC}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_migration_beforeCuts = new TH2D("h_theta_migration_beforeCuts","; #theta_{MC};#theta_{RECO}",100,0,3.14,100,0,3.14);
    TH2D* h_theta_resolution_beforeCuts = new TH2D("h_theta_resolution_beforeCuts",";#theta_{MC}; (#theta_{MC}-#theta_{RECO})/#theta_{MC}",100,0,3.14,1000,-1,1);
    TH1D* h_theta_res2_beforeCuts = new TH1D("h_theta_res2_beforeCuts",";(#theta_{RECO}-#theta_{MC})/#theta_{MC}",100,-0.2,0.2);
    TH2D* h_VM_pt_response_beforeCuts = new TH2D("h_VM_pt_response_beforeCuts","; p_{T,VM,RECO} [GeV/c];p_{T,VM,MC} [GeV/c]",200,0,10,200,0,10);
    TH2D* h_VM_pt_migration_beforeCuts = new TH2D("h_VM_pt_migration_beforeCuts",";p_{T,VM,MC} [GeV/c];p_{T,VM,RECO} [GeV/c]",200,0,10,200,0,10);
    TH2D* h_VM_pt_resolution_beforeCuts = new TH2D("h_VM_pt_resolution_beforeCuts",";p_{T,MC} [GeV/c]; (p_{T,MC}-p_{T,RECO})/p_{T,MC}",200,0,10,1000,-1,1);
    TH1D* h_VM_pt_res2_beforeCuts = new TH1D("h_VM_pt_res2_beforeCuts",";(p_{T,RECO}-p_{T,MC})/p_{T,MC}",100,-0.2,0.2);
    TH2D* h_VM_Epz_res_beforeCuts = new TH2D("h_VM_Epz_res_beforeCuts",";(E-p_{z})_{MC} [GeV]; ((E-p_{z})_{MC}-(E-p_{z})_{RECO})/(E-p_{z})_{MC}",200,0,20,200,-1,1);
    TH2D* h_VM_Epz_response_beforeCuts = new TH2D("h_VM_Epz_response_beforeCuts","; (E_{VM,RECO}-p_{z,VM,RECO}) [GeV];(E_{VM,MC}-p_{z,VM,MC} [GeV])",200,0,20,200,0,20);
    TH2D* h_VM_Epz_migration_beforeCuts = new TH2D("h_VM_Epz_migration_beforeCuts",";(E_{VM,MC}-p_{z,VM,MC}) [GeV];(E_{VM,RECO}-p_{z,VM,RECO}) [GeV]",200,0,20,200,0,20);
    TH1D* h_VM_Epz_res2_beforeCuts = new TH1D("h_VM_Epz_res2_beforeCuts",";((E-p_{z})_{RECO}-(E-p_{z})_{MC})/(E-p_{z})_{MC}",100,-0.2,0.2);

    double totalMCevents = 0;
    double totalEventsWithForward = 0;
    double totalEventsWithZDC = 0;
    double n_total = 0; 
    double n_vetoed = 0;
    double n_after_veto = 0;
    double n_sel = 0;
    double n_sel_vetoed = 0;
    double nPassMassRange = 0;
    double nPassZDC = 0;
    double nPassRP = 0;
    double nPassOMD = 0;
    double nPassEtaHFS = 0;
    double nPassEventCuts = 0;
    double nPassRec_pi12 = 0;

    long nAll = 0;
    long nAll_before = 0; 
    long nAfterElectron = 0; 
    long nAfterHFS2 = 0; 
    long nAfterKaons2 = 0; 
    long nAfterEta = 0;
    long nAfterIncoherent = 0;
    long nAfterBeam = 0;
    long nAfterVM = 0;
    long nAfterQ2Y = 0;
    long nAfterTruthT = 0;
    long nAfterScatClus = 0;
    long nAfterRadius = 0;
    long nAfterRecVM = 0;
    long nAfterIncoherentRec = 0;
    long nAfterVMrec = 0;
    long nAfter2Tracks = 0;

    int nEpz = 0;
    int nEoP = 0;
    int ny = 0;
    int nQ2 = 0;
    int nMCy = 0;
    int nMCQ2 = 0;
    int nMCincoherent = 0;

    chain->GetEntries();
    while (tree_reader.Next()) 
    {   
        /*---------------
            Event Loop
        ----------------*/

        event.clusters_eemc.clear();
        event.zdc_clusters.clear();
        event.zdc_REC_hits.clear();
        event.zdc_neutrals.clear();
        event.hit_rp.clear();
        event.hit_omd.clear();
        
        nAll++;
    
        int nNeutrons = 0;
        int nNeutrons_accept = 0;
        n_total++;
        bool have_truth = false;
        bool have_reco  = false;

        double t_total = 0;
        double t_total_rec = 0;
        
        h_Nevents->Fill(0.5);
        cout << "Total events processed" << h_Nevents->GetEntries() << endl;

        // Beam particles
    	TLorentzVector ebeam(0,0,0,0);
    	TLorentzVector pbeam(0,0,0,0);
        TLorentzVector abeam(0,0,0,0);
        TLorentzVector vmMC(0,0,0,0);
        TLorentzVector rho_MC(0,0,0,0);
        TLorentzVector vmMC_before(0,0,0,0);
    	TLorentzVector kplusMC(0,0,0,0);
    	TLorentzVector kminusMC(0,0,0,0);
        TLorentzVector piplusMC(0,0,0,0);
    	TLorentzVector piminusMC(0,0,0,0);
        TLorentzVector hfsMC(0,0,0,0);
        TLorentzVector particleMC(0,0,0,0);

    	//MC level
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
        int incoherent=0;

        // loop over all MC particles in the event
    	for(int imc=0;imc<mc_px_array.GetSize();imc++)
        {
            TVector3 mctrk(static_cast<double>(mc_px_array[imc]), 
                static_cast<double>(mc_py_array[imc]), 
                static_cast<double>(mc_pz_array[imc]));
            particleMC.SetVectM(mctrk,MASS_PION); // assume pions;
    		if(mc_genStatus_array[imc]==4) // 4 is Sartre.
            {
                /*--------------------------------------
                    -genStatus = 4 is Sartre
                    -genStatus = 1 is stable particle
                    -PDG 11 = electron
                    -PDG 2212 = proton
                    -PDG 321 = kaon
                    -PDG 211 = pion
                    -PDG 2112 = neutron
                    -All other PDG set as ion
                ---------------------------------------*/
                if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
                else
                {
                    double MASS_A = mc_mass_array[imc];
                    abeam.SetVectM(mctrk,MASS_A);
                }
            } 
            // Select neutrons
            if (mc_pdg_array[imc] == 2112 && (mc_genStatus_array[imc] == 1 || mc_genStatus_array[imc] == 4))
            {
                nNeutrons++;
                double eta = mctrk.Eta();
                h_eta_neutrons->Fill(eta);
                if(eta>4 && eta<4.8)
                {
                    nNeutrons_accept++;
                    h_eta_ZDC_accept->Fill(eta);
                }
            }
            if(mc_genStatus_array[imc]!=1) continue;
    		if(mc_pdg_array[imc]==11 && mctrk.Perp()>maxPt)
            {
                /*---------------------------------------------------
                        electron with highest pT = scattered electron
                    ----------------------------------------------------*/
    			maxPt=mctrk.Perp();
    			mc_elect_index=imc;
    			scatMC.SetVectM(mctrk,mc_mass_array[imc]); 
    		}
            if(isRho)
            {
                //if particle is stable and pion PDG == 211
        		if(mc_pdg_array[imc]==211 
                    && mc_genStatus_array[imc]==1) piplusMC.SetVectM(mctrk,MASS_PION);
        		if(mc_pdg_array[imc]==-211 
                    && mc_genStatus_array[imc]==1) piminusMC.SetVectM(mctrk,MASS_PION);
            }
            else
            {
                // if particle is stable and kaon PDG == 321
        		if(mc_pdg_array[imc]==321 
                    && mc_genStatus_array[imc]==1) kplusMC.SetVectM(mctrk,MASS_KAON);
        		if(mc_pdg_array[imc]==-321 
                    && mc_genStatus_array[imc]==1) kminusMC.SetVectM(mctrk,MASS_KAON);
            }
            if(imc!=mc_elect_index) hfsMC += particleMC;
        }
        if(isRho) 
        {
            vmMC=piplusMC+piminusMC;
        }
        else
        {
            vmMC=kplusMC+kminusMC;
        }
        if(vmMC.E()==0) continue;
        h_multiplicity_neutrons->Fill(nNeutrons);
        h_multiplicity_neutrons_accept->Fill(nNeutrons_accept);

        // checks
        cout<<"A energy: "<<abeam.E()<<" p Energy: "<<pbeam.E()<<" e Energy: "<<ebeam.E()<<endl;
        cout<<"Electron Beam: "<<" px: "<<ebeam.Px()<<" py: "<<ebeam.Py()<<" pz: "<<ebeam.Pz()<< " E: "<<ebeam.E()<< "e beam theta: " << ebeam.Theta() << endl;
		cout<<"p Beam: "<<" px: "<<pbeam.Px()<<" py: "<<pbeam.Py()<<" pz: "<<pbeam.Pz()<<" E: "<<pbeam.E()<< endl;
        cout <<"A beam: "<<" px: "<< abeam.Px()<<" py: "<< abeam.Py()<<" pz: "<< abeam.Pz()<< " E: "<< abeam.E()<< endl;
		cout <<"Scattered Electron: "<<" px: "<<scatMC.Px()<<" py: "<<scatMC.Py()<<" pz: "<< scatMC.Pz()<<" E: "<< scatMC.E()<<endl;
    
    	// protection.
    	if(ebeam.E()==abeam.E() && ebeam.E()==0) 
        {
            nAfterBeam++;
    		cout << "problem with MC incoming beams" << endl;
    		continue;
    	}

        // VM
        double vm_Epz_MC = vmMC.E()-vmMC.Pz();
        double vm_pT_MC = vmMC.Pt();
        double vm_mass_MC = vmMC.M();
        double vm_pz_MC = vmMC.Pz();
        double vm_p_MC = vmMC.P();
        double theta_VM_MC = vmMC.Theta();
        double vm_phi_MC = vmMC.Phi();
    
        double daughterPlus_mass_MC  = isRho ? piplusMC.M()  : kplusMC.M();
        double daughterMinus_mass_MC = isRho ? piminusMC.M() : kminusMC.M();
        h_VM_daughterPlus_mass_MC->Fill(daughterPlus_mass_MC);
        h_VM_daughterMinus_mass_MC->Fill(daughterMinus_mass_MC);
        
        h_VM_theta_MC->Fill(theta_VM_MC);
        h_VM_phi_MC->Fill(vm_phi_MC);
        h_VM_mass_MC->Fill(vm_mass_MC);
        h_VM_pt_MC->Fill(vm_pT_MC);
    	h_VM_pz_MC->Fill(vm_pz_MC);
    	h_VM_p_MC->Fill(vm_p_MC);
        h_VM_Epz_MC->Fill(vm_Epz_MC);

        // electron
        double phi_MC = scatMC.Phi();
        double theta_MC = scatMC.Theta();
        double eta_MC = scatMC.Eta();
        double e_pT_MC = scatMC.Pt();
        double e_pz_MC = scatMC.Pz();
        double e_p_MC = scatMC.P();
        double energy_MC = scatMC.E();
        double EoverP_MC = scatMC.E()/scatMC.P();
        double Epz_MC = (scatMC+hfsMC).E()-(scatMC+hfsMC).Pz();
        
        h_phi_MC->Fill(phi_MC);
        h_theta_MC->Fill(theta_MC);
        h_eta_MC->Fill(eta_MC);
        h_e_pt_MC->Fill(e_pT_MC);
		h_e_pz_MC->Fill(e_pz_MC);
		h_e_p_MC->Fill(e_p_MC);
		h_energy_MC->Fill(energy_MC);
        h_EoverP_MC->Fill(EoverP_MC);
        h_Epz_MC->Fill(Epz_MC);

    	TLorentzVector qbeam = ebeam-scatMC; // p_e - p_e'

        // electron method
        double Q2 = 2*ebeam.E()*energy_MC*(1+TMath::Cos(theta_MC));
        double y = 1-(energy_MC/ebeam.E())*TMath::Sin(theta_MC/2)*TMath::Sin(theta_MC/2);
        double BjorkenX = Q2/(4*ebeam.E()*abeam.E()*y);

        h_Q2_e->Fill(Q2);
        h_y_e->Fill(y);
        h_x->Fill(BjorkenX);
        h_Q2_vs_x_MC->Fill(BjorkenX,Q2);

		// t dist
        double t_MC = 0;
        double method_E = -(qbeam-vmMC).Mag2();  // t = -(p_e - p_e' - p_v)^2
        double t_MC_perp = 0;
        double method_E_perp = -(qbeam.Pt()-vmMC.Pt())*(qbeam.Pt()-vmMC.Pt());  // t = -(p_e - p_e' - p_v)^2
    		
        t_MC = method_E;
        t_MC_perp = method_E_perp;
        h_MC_neutrons_vs_t_2Dv2->Fill(nNeutrons, t_MC);
    	h_t_MC->Fill(method_E);
            
    	/*------------------------------------------------
            projection method
            - Same as method L but takes a theta_max cut
        -------------------------------------------------*/
    	TLorentzVector T = -ebeam + scatMC + vmMC;  // true t (method E)
            
    	// define nHat direction -> direction perpendicular to electron scattering plane
    	TVector3 eScattered_momentum = scatMC.Vect(); 
    	TVector3 e_momentum = ebeam.Vect();
    	TVector3 nHat = e_momentum.Cross(eScattered_momentum); // nHat = p_e x p_e'
    	nHat = nHat.Unit();
        double nHat_MC = nHat.Mag2();
        
        // define projected VM (y-component of t)
    	double pv = vmMC.Vect().Dot(nHat);
    	TVector3 py_vector = nHat*pv;  // VM in nhat direction
    	TLorentzVector py(py_vector.X(),py_vector.Y(),0,0);
        double ty = -(py).Mag2();
        double qy = sqrt(ty);
        
        // define x and z-components
        TLorentzVector px = T - py; 
        TLorentzVector pz(0,0,px.Z(),px.E());
    	TLorentzVector px_true = T - py - pz;
        double tx_true = -(px_true).Mag2();
        double qx_true = sqrt(tx_true);
        double tz = -(pz).Mag2();

        TLorentzVector q_MC_vect(qx_true,qy,sqrt(tz),0);
        double q_MC = q_MC_vect.Mag2();
        double q_MC_perp = q_MC_vect.Pt();
        
        t_total = tx_true+ty+tz; // no resolution added so should still be t_MC      
        have_truth = true;

    	h_t_REC_2d->Fill(qx_true,qy);
        h_t_REC_new_method->Fill(t_total);
       
    	// apply cut with qz and E subtracted out of qx-component (just for testing, distributions should match method E)
    	double theta = atan(fabs(qx_true)/fabs(qy));
        h_theta_all->Fill(theta);
    	if(fabs(theta)<PI/2) 
    	{
    		h_t_REC_wCUT_pi2->Fill(t_total);
    		h_t_REC_2d_wCUT_pi2->Fill(qx_true,qy);
    	}
    	if(fabs(theta)<PI/3) 
    	{
    		h_t_REC_wCUT_pi3->Fill(t_total);
    		h_t_REC_2d_wCUT_pi3->Fill(qx_true,qy);
    	}
    	if(fabs(theta)<PI/4) 
    	{
    		h_t_REC_wCUT_pi4->Fill(t_total);
    		h_t_REC_2d_wCUT_pi4->Fill(qx_true,qy);
    	}
    	if(fabs(theta)<PI/6) 
    	{
    		h_t_REC_wCUT_pi6->Fill(t_total);
    		h_t_REC_2d_wCUT_pi6->Fill(qx_true,qy);
    	}
    	if(fabs(theta)<PI/9) 
    	{
    		h_t_REC_wCUT_pi9->Fill(t_total);
    		h_t_REC_2d_wCUT_pi9->Fill(qx_true,qy);
    	}
    	if(fabs(theta)<PI/12) 
    	{
    		h_t_REC_wCUT_pi12->Fill(t_total);
    		h_t_REC_2d_wCUT_pi12->Fill(qx_true,qy);
            h_MC_neutrons_vs_t_2D->Fill(nNeutrons, t_total);
    	}
        if(fabs(theta)<PI/16) 
    	{
            h_t_REC_wCUT_pi16->Fill(t_total);
    		h_t_REC_2d_wCUT_pi16->Fill(qx_true,qy);
    	}
        if(fabs(theta)<PI/20) 
    	{
            h_t_REC_wCUT_pi20->Fill(t_total);
    		h_t_REC_2d_wCUT_pi20->Fill(qx_true,qy);
    	}
        if(fabs(theta)<PI/24) 
    	{
            h_t_REC_wCUT_pi24->Fill(t_total);
    		h_t_REC_2d_wCUT_pi24->Fill(qx_true,qy);
    	}
        if(applyMCcuts==true)
        {
            // MC level phase space cut
        	if(Q2<1.||Q2>10.) continue;
        	if(y<0.01||y>0.85) continue;
            
            h_phi_MC_after->Fill(phi_MC);
            h_theta_MC_after->Fill(theta_MC);
            h_eta_MC_after->Fill(eta_MC);
            h_e_pt_MC_after->Fill(e_pT_MC);
    		h_e_pz_MC_after->Fill(e_pz_MC);
    		h_e_p_MC_after->Fill(e_p_MC);
    		h_energy_MC_after->Fill(energy_MC);
            h_EoverP_MC_after->Fill(EoverP_MC);
            h_Epz_MC_after->Fill(Epz_MC);
            h_EvsP_MC->Fill(e_p_MC,energy_MC);
            h_Q2_MC_after->Fill(Q2);
        	h_y_MC_after->Fill(y);
            h_x_MC_after->Fill(BjorkenX);
            h_VM_Epz_MC_after->Fill(vm_Epz_MC);
            h_VM_pt_MC_after->Fill(vm_pT_MC);
            h_t_MC_after->Fill(method_E);
            
        }

        // rec level 
        /*--------------------------------------------------------------------------------
            -find highest energy cluster in negative endcap to find scattered electron
        ---------------------------------------------------------------------------------*/
    	double maxEnergy=-99.;
    	double xpos=-999.;
    	double ypos=-999.;
                
        // find highest energy cluster in negative endcap to find scattered electron
    	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++)
        {
        	if(em_energy_array[iclus]>maxEnergy)
            {
                /*---------------------------------------------------------------
                    -check if cluster has higher energy
                    -if it does, define as maxEnergy and define that position
                ----------------------------------------------------------------*/
        		maxEnergy=em_energy_array[iclus];
    			xpos=em_x_array[iclus]; 
    			ypos=em_y_array[iclus];
    		}
    	}
                
    	/*----------------------------------------------------------------------------------
            find highest energy hits (above threshold) and store that position and index
        -----------------------------------------------------------------------------------*/
    	double maxHitEnergy=0.01;  // threshold 10 MeV
    	double xhitpos=-999.;
    	double yhitpos=-999.;
    	int hit_index=-1;
                
        // individual calorimeter cells that registered energy
    	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++)
        {	
    		if(emhits_energy_array[ihit]>maxHitEnergy)
            {
                // find highest energy hits (above threshold) and store that position and index
    			maxHitEnergy=emhits_energy_array[ihit];
			    xhitpos=emhits_x_array[ihit];
			    yhitpos=emhits_y_array[ihit];
			    hit_index=ihit;
    		}
    	}

    	/*----------------------------------------------------------------------------------------------
            -sum over all 3x3 towers around the leading tower
            -agregate nearby hits to form a shower to reconstruct particles energy and impact point
        -----------------------------------------------------------------------------------------------*/
    	double xClus=xhitpos*maxHitEnergy;
    	double yClus=yhitpos*maxHitEnergy;
    	for(int ihit=0;ihit<emhits_energy_array.GetSize();ihit++)
        {
    		double hitenergy=emhits_energy_array[ihit];
    		double x=emhits_x_array[ihit];
    		double y=emhits_y_array[ihit];
            // sum over energies within 70 mm of maxHitEnergy
    		double d=sqrt( (x-xhitpos)*(x-xhitpos) + (y-yhitpos)*(y-yhitpos));
    		if(d<70. && ihit!=hit_index && hitenergy>0.01)  
            {
    			maxHitEnergy+=hitenergy;
    			xClus+=x*hitenergy;
    			yClus+=y*hitenergy;
    		}
    	}
	
    	// weighted average cluster position.
    	xClus = xClus/maxHitEnergy; // built from hits
    	yClus = yClus/maxHitEnergy;
    	double radius=sqrt(xClus*xClus+yClus*yClus);

    	double clusEnergy=1.044*maxHitEnergy; // 4.4% energy calibration.
    	h_energy_REC_EEMC->Fill(clusEnergy);
    	h_emClus_position_REC->Fill(xClus,yClus); 

        double min_distance = 1e6;
        double xtrk = -999.;
        double ytrk = -999.;
                
        // get positions from calorimeter track
        for(int itrk=0; itrk<emCal_trk_x_array.GetSize(); itrk++) 
        {
            double dx = xpos - emCal_trk_x_array[itrk];
            double dy = ypos - emCal_trk_y_array[itrk];
            double dist = sqrt(dx*dx + dy*dy);
            if(dist<min_distance) 
            {
                min_distance = dist;
                xtrk = emCal_trk_x_array[itrk];
                ytrk = emCal_trk_y_array[itrk];
            }
        }
        
        double x_clus_trk_diff = xpos-xtrk; // built from EE cluster
        double x_hit_trk_diff = xClus-xtrk; // built from hits
        double y_clus_trk_diff = ypos-ytrk;
        double y_hit_trk_diff = yClus-ytrk;
        
        // energy resolution
        double res= (energy_MC-clusEnergy)/energy_MC;
        h_energy_res_EEMC->Fill(energy_MC,res);
        h_energy_response_EEMC->Fill(energy_MC,clusEnergy);
        
        /*---------------------------------------------------------------------
                -association of rec level scat' e
                -sim_id = index of simulated MC particle 
                -rec_id = index of reconstructed particle
        		-find the rec_id that matches MC scattered electron (sim_id)
            ----------------------------------------------------------------------*/
        int rec_elect_index=-1;
        for(int i=0;i<sim_id.GetSize();i++)
        {
        	if(sim_id[i]==mc_elect_index) // mc_elect_index found when defining scatMC
            {
        		rec_elect_index = rec_id[i]; 
        	}
        }

        // reconstructed particles
        TLorentzVector scatMCmatchREC(0,0,0,0);
        TLorentzVector scatREC(0,0,0,0);
        TLorentzVector scatClusEREC(0,0,0,0);
        TLorentzVector hfs(0,0,0,0);
        TLorentzVector particle(0,0,0,0);
        TLorentzVector kplusREC(0,0,0,0);
        TLorentzVector kminusREC(0,0,0,0);
        TLorentzVector kplusREC_temp(0,0,0,0);
        TLorentzVector kminusREC_temp(0,0,0,0);
        TLorentzVector piplusREC(0,0,0,0);
        TLorentzVector piminusREC(0,0,0,0);
        TLorentzVector vmREC(0,0,0,0);
        TLorentzVector rho_REC(0,0,0,0);
        TLorentzVector vmREC_before(0,0,0,0);

        double maxP=-1.;
        double scat_electron = 0;
        // track loop
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
            int q = reco_charge_array[itrk];  
        	if(rec_elect_index!=-1 && itrk==rec_elect_index && trk.Mag()>maxP && q<0) 
            {
                // if stable and matches electron index
        		maxP=trk.Mag();
        		scatREC.SetVectM(trk,MASS_ELECTRON);
        		// use emcal energy to define 4 vector
        		double p = sqrt(clusEnergy*clusEnergy- MASS_ELECTRON*MASS_ELECTRON );
                double pt = TMath::Sin(scatREC.Theta())*p;
        		scatClusEREC.SetPtEtaPhiM(pt,scatREC.Eta(),scatREC.Phi(),MASS_ELECTRON);
        	}
        }
        double e_pT_REC_cal = scatClusEREC.Pt();
        double e_pz_REC_cal = scatClusEREC.Pz();
        double e_p_REC_cal = scatClusEREC.P();
        double e_pT_REC_trk = scatREC.Pt();
        double e_pz_REC_trk = scatREC.Pz();
        double e_p_REC_trk = scatREC.P();
    
        double eta_e = scatClusEREC.Eta();
        double eta_REC = scatREC.Eta();
        double phi_REC = scatREC.Phi();
        double theta_REC = scatREC.Theta();
        
        double energy_REC = scatClusEREC.E();
        double EoverP_REC = energy_REC/e_p_REC_trk;
        double Epz_REC = scatClusEREC.E() - scatREC.Pz();
        
        if(applyCuts)
        {
            if(rec_elect_index < 0) {h_didntPass_e->Fill(0.5); nAfterElectron++; continue;} 
            else
            {
                scat_electron++;
                h_pass_e->Fill(0.5);
            }
            if (scatClusEREC.E() == 0) {nAfterScatClus++; continue;}
            if( Epz_REC<15||Epz_REC>25 ) {nEpz++; continue;}
            if( EoverP_REC<0.8||EoverP_REC>1.2 ) {nEoP++; continue;} 
        }
        h_t_pi12_afterElectronSelect->Fill(t_MC); 
        h_EoverP_beforeCut->Fill(EoverP_REC);
        h_e_energy_beforeCuts->Fill(energy_REC);
        h_e_Epz_noHFS->Fill(Epz_REC);
        h_Epz_migration_noHFS->Fill(Epz_MC,Epz_REC);

        // cluster-base DIS kine;
        TLorentzVector qbeamREC=ebeam-scatClusEREC; 
        
        // electron method
        double Q2REC = 2*ebeam.E()*energy_REC*(1+TMath::Cos(theta_REC));
        double yREC = 1-(energy_REC/ebeam.E())*TMath::Sin(theta_REC/2)*TMath::Sin(theta_REC/2);
        double BjorkenX_REC = Q2REC/(4*ebeam.E()*abeam.E()*yREC);

        h_Q2REC_e_EEMC->Fill(Q2REC);
        h_Q2_res_beforecut->Fill(Q2,(Q2-Q2REC)/Q2);
        h_Q2_response_beforecut->Fill(Q2REC,Q2);
        h_yREC_e_EEMC->Fill(yREC);
        h_x_REC->Fill(BjorkenX_REC);

        // Q2 reco
        double resQ2_beforeCuts = (Q2-Q2REC)/Q2;
        h_Q2_res_beforeCuts->Fill(Q2,resQ2_beforeCuts);
        h_Q2_res2_beforeCuts->Fill((Q2REC-Q2)/Q2);
        h_Q2_response_beforeCuts->Fill(Q2,Q2REC);
        h_Q2_migration_beforeCuts->Fill(Q2,Q2REC);
        h_Q2_vs_x_REC_beforeCuts->Fill(BjorkenX_REC,Q2REC);
                                
        // x reco
        double resx_beforeCuts = (BjorkenX-BjorkenX_REC)/BjorkenX;
        h_x_res_cut_beforeCuts->Fill(BjorkenX,resx_beforeCuts);
        h_x_res_cut2_beforeCuts->Fill((BjorkenX_REC-BjorkenX)/BjorkenX);
        h_x_response_cut_beforeCuts->Fill(BjorkenX,BjorkenX_REC);       
        h_x_migration_beforeCuts->Fill(BjorkenX,BjorkenX_REC);       
        
        // y reco
        double resy_beforeCuts = (y-yREC)/y;
        h_y_res_beforeCuts->Fill(y,resy_beforeCuts);
        h_y_res2_beforeCuts->Fill((yREC-y)/y);
        h_y_response_beforeCuts->Fill(y,yREC);
        h_y_migration_beforeCuts->Fill(y,yREC);

        // E over p
        double resEoP_beforeCuts = (EoverP_MC-EoverP_REC)/EoverP_MC;
        h_EoverP_res_beforeCuts->Fill(EoverP_MC,resEoP_beforeCuts);
        h_EoverP_response_beforeCuts->Fill(EoverP_REC,EoverP_MC);
        h_EoverP_migration_beforeCuts->Fill(EoverP_MC,EoverP_REC);
        h_EoverP_res2_beforeCuts->Fill((EoverP_REC-EoverP_MC)/EoverP_MC);

        // theta
        h_theta_response_EEMC_beforeCuts->Fill(theta_REC,theta_MC);
        double restheta_beforeCuts = (theta_MC-theta_REC)/theta_MC;
        h_theta_resolution_beforeCuts->Fill(theta_MC,restheta_beforeCuts);
        h_theta_migration_beforeCuts->Fill(theta_MC,theta_REC);
        h_theta_res2_beforeCuts->Fill((theta_REC-theta_MC)/theta_MC);

        // energy
        double resE_beforeCuts= (energy_MC-energy_REC)/energy_MC;
        h_energy_res_EEMC_beforeCuts->Fill(energy_MC, resE_beforeCuts);
        h_energy_response_EEMC_beforeCuts->Fill(energy_REC,energy_MC);
        h_energy_migration_beforeCuts->Fill(energy_MC,energy_REC);
        h_energy_res2_beforeCuts->Fill((energy_REC-energy_MC)/energy_MC);

        h_eta_REC_had_in_ZDC_e->Fill(eta_e); 
        h_eta_REC_had_in_ZDC->Fill(eta_REC); 
        h_eta_REC_EEMC->Fill(eta_REC);
        h_phi_REC_EEMC->Fill(phi_REC);
        h_theta_REC_EEMC->Fill(theta_REC);
        h_theta_response_EEMC->Fill(theta_REC,theta_MC);
        h_theta_diff->Fill(theta_REC-theta_MC);
        h_phi_diff->Fill(phi_REC-phi_MC);
        h_eta_diff->Fill(eta_REC-eta_MC);

        // default clustering position
    	h_emClus_position_x_REC->Fill(xpos);
    	h_emClus_position_y_REC->Fill(ypos);
        h_XvsY_clus->Fill(xpos,ypos);
        
    	// self clustering position
    	h_emHits_position_x_REC->Fill(xClus);
    	h_emHits_position_y_REC->Fill(yClus);
        h_XvsY_hits->Fill(xClus,yClus);

        // track positions
        h_XvsY_trk->Fill(xtrk,ytrk);
        h_trk_position_x_REC->Fill(xtrk);
        h_trk_position_y_REC->Fill(ytrk);

        h_Xclus_minus_Xtrk->Fill(x_clus_trk_diff);
        h_Yclus_minus_Ytrk->Fill(y_clus_trk_diff);
        h_Xhits_minus_Xtrk->Fill(x_hit_trk_diff);
        h_Yhits_minus_Ytrk->Fill(y_hit_trk_diff);
        
    	// track-base e' energy
        double energy_REC_trk = scatREC.E();
    	double res_trk= (energy_MC-energy_REC_trk)/energy_MC;
        h_energy_res_trk->Fill(energy_MC,res_trk);
    	h_energy_REC_trk->Fill(energy_REC_trk);

        int incoherent_rec = 0;
        int total_hfs_particles = 0;
        int kaons = 0;
        int extra_hfs = 0;
        int eta_out = 0;
        double total_kaons = 0;
        double two_hfs = 0;
        double hfs_kaons = 0;
        bool pass_etaHFS  = true; 
        bool pass_e = true;
        bool pass_HFS = true;
        bool pass_eta = true;
        bool pass_kaons = true;

        // Reset per event
        vector<KaonCand> kaonCandidates;
        kaonCandidates.clear();
        // Loop over tracks
        for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++)
        {
        	TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
        	particle.SetVectM(trk,MASS_PION); // assume pions;
            //------------------------------------------------------
            //        -exclude e' from hadronic final state and kaons
            //        -select phi->kk daughters;
            //-------------------------------------------------------
        	if(itrk!=rec_elect_index) 
            {
            	hfs += particle;
                total_hfs_particles++;
            	h_eta_REC_trk->Fill(trk.Eta());
                //if(reco_charge_array[itrk]<0 || reco_charge_array[itrk]>0)
                //{
                    if(usePID)
                    {
                        // Identify kaons
                        if(reco_pdg_array[itrk] == 321 || reco_pdg_array[itrk] == -321)
                        {
                            TLorentzVector k;
                            k.SetVectM(trk, MASS_KAON);
                            kaonCandidates.push_back({k, reco_pdg_array[itrk]});
                        }
                    }
                    if(!usePID)
                    {
                        if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
                    	if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
                    }
                //}
        	}
        }
        if(usePID)
        {
            // Check if found exactly K+K-
            if (kaonCandidates.size() == 2) 
            {
                for (auto &kc : kaonCandidates) 
                { 
                    if (kc.pdg == 321) kplusREC = kc.pKaon; 
                    else if (kc.pdg == -321) kminusREC = kc.pKaon; 
                }
                if (kplusREC.E() > 0 && kminusREC.E() > 0) 
                { 
                    vmREC = kplusREC + kminusREC; 
                    kaons = 2;
                } 
            }
        }
        if(!usePID)
        {
            if (kplusREC.E() > 0 && kminusREC.E() > 0) 
            {
                vmREC=kplusREC+kminusREC;
                kaons = 2;
            }
        }
        h_total_hfs_particles->Fill(total_hfs_particles);
        double daughterPlus_mass_REC  = kplusREC.M();
        double daughterMinus_mass_REC = kminusREC.M();
        double vm_mass_REC = vmREC.M();
        double vm_pT_REC = vmREC.Pt();
        double vm_pz_REC = vmREC.Pz();
        double vm_p_REC = vmREC.P();
        double vm_Epz_REC = vmREC.E()-vmREC.Pz();

        h_VM_daughterPlus_mass_REC_before->Fill(daughterPlus_mass_REC);
        h_VM_daughterMinus_mass_REC_before->Fill(daughterMinus_mass_REC);
        h_VM_mass_REC_before->Fill(vm_mass_REC);
        h_VM_pt_REC_before->Fill(vm_pT_REC);
        h_VM_pz_REC_before->Fill(vm_pz_REC);
        h_VM_p_REC_before->Fill(vm_p_REC);
        h_VM_Epz_REC_before->Fill(vm_Epz_REC);
        h_t_after_VMrec->Fill(t_MC);

        double resEpzVM_beforeCuts = (vm_Epz_MC-vm_Epz_REC)/vm_Epz_MC;
        h_VM_Epz_response_beforeCuts->Fill(vm_Epz_REC,vm_Epz_MC);
        h_VM_Epz_migration_beforeCuts->Fill(vm_Epz_MC,vm_Epz_REC);
        h_VM_Epz_res_beforeCuts->Fill(vm_Epz_MC,resEpzVM_beforeCuts);
        h_VM_Epz_res2_beforeCuts->Fill((vm_Epz_REC-vm_Epz_MC)/vm_Epz_MC);

        h_VM_pt_response_beforeCuts->Fill(vm_pT_REC,vm_pT_MC);    
        double respTVM_beforeCuts = (vm_pT_MC-vm_pT_REC)/vm_pT_MC;
        h_VM_pt_resolution_beforeCuts->Fill(vm_pT_MC,respTVM_beforeCuts);
        h_VM_pt_migration_beforeCuts->Fill(vm_pT_MC,vm_pT_REC);
        h_VM_pt_res2_beforeCuts->Fill((vm_pT_REC-vm_pT_MC)/vm_pT_MC);
       
        nAll_before++;
        if(applyCuts)
        {
            // HFS selection
            if(total_hfs_particles != 2) {h_didntPasstotalHFS->Fill(0.5); nAfterHFS2++; continue;}
            else
            {
                two_hfs++;
                h_totalHFS->Fill(0.5);
            }
        }
        h_t_pi12_afterHFS->Fill(t_MC);

        if(applyCuts)
        {
            if(kaons != 2) {h_didntPassnKaons->Fill(0.5); nAfterKaons2++; continue;}
            else
            {
                total_kaons = 2;
                h_nKaons->Fill(0.5);
            }
        }   
        h_t_pi12_afterKaons->Fill(t_MC);
        
        // E-pz scat' e
        double Epz_REC_plusHFS = (scatClusEREC+hfs).E() - (scatREC+hfs).Pz();
        double resEpz_beforeCuts = (Epz_MC-Epz_REC_plusHFS)/Epz_MC;
        h_Epz_res_beforeCuts->Fill(Epz_MC,resEpz_beforeCuts);
        h_Epz_response_beforeCuts->Fill(Epz_REC_plusHFS,Epz_MC);
        h_Epz_migration_beforeCuts->Fill(Epz_MC,Epz_REC_plusHFS);
        h_Epz_res2_beforeCuts->Fill((Epz_REC_plusHFS-Epz_MC)/Epz_MC);
        h_Epz_beforeCut->Fill(Epz_REC_plusHFS);

        if(applyCuts)
        {    
            // 4vector of VM;
            if(vmREC.E()==0) {nAfterRecVM++; continue;}
            if(fabs(vmREC.M()-1.02)>=0.02) {h_didntPassMass->Fill(0.5); nPassMassRange++; continue;}
            else
            {
                h_passMass->Fill(0.5);
            }
        }
        h_t_pi12_afterMass->Fill(t_MC);

        if(applyCuts)
        {
            if(Q2REC<1.|| Q2REC>10.) {nQ2++; continue;}
            if(yREC<0.01 || yREC>0.85) {ny++; continue;}
        }
        h_t_pi12_afterQ2y->Fill(t_MC);
    
        nPassEventCuts++; 
        h_passEventSelection->Fill(0.5);

        if (useZDC) 
        {
            // --- Fill ECal ZDC clusters --------------------------------------------
            for (int i = 0; i < Ecal_zdc_energy_array.GetSize(); i++) 
            {
                Cluster_ZDC cl;
                cl.x = Ecal_zdc_x_array[i];
                cl.y = Ecal_zdc_y_array[i];
                cl.z = Ecal_zdc_z_array[i];
                cl.energy = Ecal_zdc_energy_array[i];
                event.zdc_clusters.push_back(cl);
    
                TVector3 zdc_Ecal_position(Ecal_zdc_x_array[i], Ecal_zdc_y_array[i], Ecal_zdc_z_array[i]);
                h_zdc_Ecal_eta->Fill(zdc_Ecal_position.Eta());
                h_zdc_Ecal_theta->Fill(zdc_Ecal_position.Theta());
            }
            int n_clusters_ECal = Ecal_zdc_energy_array.GetSize();
            h_zdc_Ecal_multiplicity->Fill(n_clusters_ECal);
        
            // --- Fill HCal ZDC clusters --------------------------------------------
            for (int i = 0; i < Hcal_zdc_energy_array.GetSize(); i++) 
            {
                Cluster_ZDC cl;
                cl.x = Hcal_zdc_x_array[i];
                cl.y = Hcal_zdc_y_array[i];
                cl.z = Hcal_zdc_z_array[i];
                cl.energy = Hcal_zdc_energy_array[i];
                event.zdc_clusters.push_back(cl);
    
                TVector3 zdc_Hcal_position(Hcal_zdc_x_array[i], Hcal_zdc_y_array[i], Hcal_zdc_z_array[i]);
                h_zdc_Hcal_eta->Fill(zdc_Hcal_position.Eta());
                h_zdc_Hcal_theta->Fill(zdc_Hcal_position.Theta());
            }
            int n_clusters_HCal = Hcal_zdc_energy_array.GetSize();
            h_zdc_Hcal_multiplicity->Fill(n_clusters_HCal);
            if(applyVetoes)
            {    
                if(event.zdc_clusters.size()>0) {h_didntPassZDC->Fill(0.5); nPassZDC++; continue;}
                else
                {
                    h_passZDC->Fill(0.5);
                }
            }
            h_t_pi12_afterZDC->Fill(t_MC);        
        }

        int nNeutrons_RECO = 0;
        for (auto& cl : event.zdc_clusters)
        {
            if (cl.energy < 21.0) continue; // [GeV]
            TVector3 pos(cl.x, cl.y, cl.z);
            double eta = pos.Eta();
            if (eta < 4.0 || eta > 4.8) continue;
            nNeutrons_RECO++;
        }

        // --- Fill EEMC clusters --------------------------------------------
        for(int iclus=0;iclus<em_energy_array.GetSize();iclus++)
        {
      		Cluster_EEMC cluster;
			cluster.energy=em_energy_array[iclus];
			cluster.x=em_x_array[iclus];
			cluster.y=em_y_array[iclus];
    	    event.clusters_eemc.push_back(cluster);

            hEEMC_energy->Fill(cluster.energy);
            hEEMC_xy->Fill(cluster.x, cluster.y);
            
            TVector3 EEMC_position(cluster.x, cluster.y, cluster.z);
            h_EEMC_eta->Fill(EEMC_position.Eta());
            h_EEMC_theta->Fill(EEMC_position.Theta());
    	}
        hEEMC_hits->Fill(event.clusters_eemc.size());

        // --- Fill RP clusters --------------------------------------------
    	for(int ihit=0;ihit<rp_x_array.GetSize();ihit++)
        {
      		Hit_RP hit;
			hit.x=rp_x_array[ihit];
			hit.y=rp_y_array[ihit];
			hit.z=rp_z_array[ihit];
    	    event.hit_rp.push_back(hit);
            
            TVector3 rp_position(rp_x_array[ihit], rp_y_array[ihit], rp_z_array[ihit]);
            h_rp_eta->Fill(rp_position.Eta());
            h_rp_theta->Fill(rp_position.Theta());
    	}
        hRP_hits->Fill(event.hit_rp.size());
        if(applyVetoes)
        {
            if(event.hit_rp.size()>0) {h_didntPassRP->Fill(0.5); nPassRP++; continue;}
            else
            {
                h_passRP->Fill(0.5);
            }
        }
        h_t_pi12_afterRP->Fill(t_MC);

    	// --- Fill OMD clusters --------------------------------------------
    	for(int ihit=0;ihit<omd_x_array.GetSize();ihit++)
        {
      		Hit_OMD hit;
			hit.x=omd_x_array[ihit];
			hit.y=omd_y_array[ihit];
			hit.z=omd_z_array[ihit];
    	    event.hit_omd.push_back(hit);
            
            TVector3 omd_position(omd_x_array[ihit], omd_y_array[ihit], omd_z_array[ihit]);
            h_omd_eta->Fill(omd_position.Eta());
            h_omd_theta->Fill(omd_position.Theta());
    	}            
        hOMD_hits->Fill(event.hit_omd.size());
        if(applyVetoes)
        {
        	if(event.hit_omd.size()>0) {h_didntPassOMD->Fill(0.5); nPassOMD++; continue;}
            else
            {
                h_passOMD->Fill(0.5);
            }
        }
        h_t_pi12_afterOMD->Fill(t_MC);

        double Q2REC_ac = Q2REC; 
        double yREC_ac = yREC;
        cout << "Q2 REC after cuts: " << Q2REC_ac << " y REC after cuts " << yREC_ac << endl;

        // Q2 reco
        h_Q2_afterCut->Fill(Q2REC);
        double resQ2 = (Q2-Q2REC)/Q2;
        h_Q2_res->Fill(Q2,resQ2);
        h_Q2_res2->Fill((Q2REC-Q2)/Q2);
        h_Q2_response->Fill(Q2,Q2REC);
        h_dQ2overQ2_REC->Fill(resQ2);
        h_Q2_migration->Fill(Q2,Q2REC);
        h_Q2_vs_x_REC->Fill(BjorkenX_REC,Q2REC);
        h_Q2_res_vs_counts->Fill((Q2REC-Q2)/Q2);
                                
        // x reco
        h_x_afterCut->Fill(BjorkenX_REC);
        double resx = (BjorkenX-BjorkenX_REC)/BjorkenX;
        h_x_res_cut->Fill(BjorkenX,resx);
        h_x_res_cut2->Fill((BjorkenX_REC-BjorkenX)/BjorkenX);
        h_x_response_cut->Fill(BjorkenX,BjorkenX_REC);       
        h_x_migration->Fill(BjorkenX,BjorkenX_REC);       
        h_x_res_vs_counts->Fill((BjorkenX_REC-BjorkenX)/BjorkenX);
        
        // y reco
        double resy = (y-yREC)/y;
        h_y_res->Fill(y,resy);
        h_y_res2->Fill((yREC-y)/y);
        h_y_response->Fill(y,yREC);
        h_dyOvery_REC->Fill(resy);
        h_y_migration->Fill(y,yREC);
        h_y_afterCut->Fill(yREC);

        // VM reco
        double theta_VM_REC = vmREC.Theta();
        double vm_phi_REC = vmREC.Phi();
        h_VM_theta_REC->Fill(theta_VM_REC);
        h_VM_phi_REC->Fill(vm_phi_REC);
        h_VM_mass_REC->Fill(vm_mass_REC);
        h_VM_daughterPlus_mass_REC->Fill(daughterPlus_mass_REC);
        h_VM_daughterMinus_mass_REC->Fill(daughterMinus_mass_REC);
        h_vm_px->Fill(vmREC.Px());
        h_vm_py->Fill(vmREC.Py());
        h_vm_pz->Fill(vmREC.Pz());

        double resEpzVM = (vm_Epz_MC-vm_Epz_REC)/vm_Epz_MC;
        h_VM_Epz_REC->Fill(vm_Epz_REC);
        h_VM_Epz_response->Fill(vm_Epz_REC,vm_Epz_MC);
        h_VM_Epz_migration->Fill(vm_Epz_MC,vm_Epz_REC);
        h_VM_Epz_res->Fill(vm_Epz_MC,resEpzVM);
        h_VM_Epz_res_counts->Fill((vm_Epz_REC-vm_Epz_MC)/vm_Epz_MC);

        h_VM_pt_REC->Fill(vm_pT_REC);
        h_VM_pt_response->Fill(vm_pT_REC,vm_pT_MC);    
        double respTVM = (vm_pT_MC-vm_pT_REC)/vm_pT_MC;
        h_VM_pt_resolution->Fill(vm_pT_MC,respTVM);
        h_VM_pt_migration->Fill(vm_pT_MC,vm_pT_REC);
        h_VM_pt_res_counts->Fill((vm_pT_REC-vm_pT_MC)/vm_pT_MC);

        h_VM_pz_REC->Fill(vm_pz_REC);
        double respzVM = (vm_pz_MC-vm_pz_REC)/vm_pz_MC;
        h_VM_pz_resolution->Fill(vm_pz_MC,respzVM);

        h_VM_p_REC->Fill(vm_p_REC);
        double respVM = (vm_p_MC-vm_p_REC)/vm_p_MC;
        h_VM_p_resolution->Fill(vm_p_MC,respVM);    

        double respKplus = (kplusMC.P()-kplusREC.P())/kplusMC.P();
        h_kPlus_p_resolution->Fill(kplusMC.P(),respKplus);
        double respKminus = (kminusMC.P()-kminusREC.P())/kminusMC.P();
        h_kMinus_p_resolution->Fill(kminusMC.P(),respKminus);
        h_kPlus_p->Fill(kplusREC.P());
        h_kPlus_phi->Fill(kplusREC.Phi());
        h_kMinus_p->Fill(kminusREC.P());
        h_kMinus_phi->Fill(kminusREC.Phi());
        h_kPlus_theta->Fill(kplusREC.Theta());
        h_kMinus_theta->Fill(kminusREC.Theta());
                                
        // e' p from EEMC reco
        h_e_pt_REC_EEMC->Fill(e_pT_REC_cal);
        h_e_pz_REC_EEMC->Fill(e_pz_REC_cal);
        h_e_p_REC_EEMC->Fill(e_p_REC_cal);
        h_e_pt_res->Fill(e_pT_MC,(e_pT_REC_cal-e_pT_MC)/e_pT_MC);
        h_e_pz_res->Fill(e_pz_MC,(e_pz_REC_cal-e_pz_MC)/e_pz_MC);
        h_e_p_res->Fill(e_p_MC,(e_p_REC_cal-e_p_MC)/e_p_MC);
        h_e_pz_response->Fill(e_pz_REC_cal,e_pz_MC);

        // e' p from track reco
        h_e_pt_REC_trk->Fill(e_pT_REC_trk);
        h_e_pz_REC_trk->Fill(e_pz_REC_trk);
        h_e_p_REC_trk->Fill(e_p_REC_trk);
        h_e_pz_trk_res->Fill(e_pz_MC,(e_pz_REC_trk-e_pz_MC)/e_pz_MC);
        h_e_p_trk_res->Fill(e_p_MC,(e_p_REC_trk-e_p_MC)/e_p_MC);

        // E-pz scat' e
        double resEpz = (Epz_MC-Epz_REC_plusHFS)/Epz_MC;
        h_Epz_res->Fill(Epz_MC,resEpz);
        h_Epz_REC->Fill(Epz_REC_plusHFS);
        h_Epz_response->Fill(Epz_REC_plusHFS,Epz_MC);
        h_Epz_migration->Fill(Epz_MC,Epz_REC_plusHFS);
        h_Epz_res_counts->Fill((Epz_REC_plusHFS-Epz_MC)/Epz_MC);

        // E over p
        double resEoP = (EoverP_MC-EoverP_REC)/EoverP_MC;
        h_EoverP_res->Fill(EoverP_MC,resEoP);
        h_EoverP_afterCut->Fill(EoverP_REC);
        h_EoverP_response->Fill(EoverP_REC,EoverP_MC);
        h_EoverP_migration->Fill(EoverP_MC,EoverP_REC);
        h_EoverP_res_counts->Fill((EoverP_REC-EoverP_MC)/EoverP_MC);

        // eta
        h_eta_REC_EEMC_after->Fill(eta_REC);
        h_eta_diff_after->Fill(eta_REC-eta_MC);

        // theta
        h_theta_REC_EEMC_after->Fill(theta_REC);
        h_theta_response_EEMC_after->Fill(theta_REC,theta_MC);
        h_theta_diff_after->Fill(theta_REC-theta_MC);
        double restheta = (theta_MC-theta_REC)/theta_MC;
        h_theta_resolution->Fill(theta_MC,restheta);
        h_theta_migration->Fill(theta_MC,theta_REC);
        h_theta_res_counts->Fill((theta_REC-theta_MC)/theta_MC);

        // phi
        h_phi_REC_EEMC_after->Fill(phi_REC);
        h_phi_diff_after->Fill(phi_REC-phi_MC);
        double resphi = (phi_MC-phi_REC)/phi_MC;
        h_phi_resolution->Fill(phi_MC,resphi);
        h_phi_response->Fill(phi_REC,phi_MC);

        // energy
        double resE= (energy_MC-energy_REC)/energy_MC;
        h_energy_res_EEMC_after->Fill(energy_MC, resE);
        h_energy_response_EEMC_after->Fill(energy_REC,energy_MC);
        h_energy_REC_EEMC_after->Fill(energy_REC);
        h_energy_migration->Fill(energy_MC,energy_REC);
        h_emClus_position_REC_cut->Fill(xClus,yClus); 
        h_energy_res_counts->Fill((energy_REC-energy_MC)/energy_MC);
                            
        // 2 versions: track and energy cluster:
        double t_trk_REC = giveme_t_method_L(ebeam,scatREC,abeam,vmREC);  // method L (track e')
        double t_REC = giveme_t_method_L(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
        h_theta_L->Fill(theta);

        h_t_REC_trk_cut->Fill(t_trk_REC);
        h_t_REC_EEMC_cut->Fill(t_REC);
        h_t_response_EEMC_cut->Fill(t_REC,t_MC);
        h_t_response_trk_cut->Fill(t_trk_REC,t_MC);

        // t track resolution 
        double resttrk = (t_MC-t_trk_REC)/t_MC;
        h_t_res_trk_cut->Fill(t_MC, resttrk);
	
        // t EEMC resolution;
        double res_percent = (t_MC-t_REC)/t_MC;
        h_t_res_EEMC_cut_percent->Fill(t_MC,res_percent);
        double rest = (t_MC-t_REC);
        h_t_res_EEMC_cut->Fill(t_MC,rest);
                                
        // projection method
        TLorentzVector T_rec = giveme_t_new_method(ebeam,scatClusEREC,abeam,vmREC); // method L (EEMC e')
          
        // define nHat direction
        TVector3 eScattered_momentum_rec = scatClusEREC.Vect(); 
        TVector3 nHat_rec = e_momentum.Cross(eScattered_momentum_rec); // nHat = p_e x p_e'
        nHat_rec = nHat_rec.Unit();
        h_nHat_phi->Fill(nHat_rec.Phi());
        double nHat_reco = nHat_rec.Mag2();

        // define projected VM (y-component)
        double pv_rec = vmREC.Vect().Dot(nHat_rec);
        TVector3 py_vector_rec = nHat_rec*pv_rec;
        TLorentzVector py_rec(py_vector_rec.X(),py_vector_rec.Y(),0,0);
        double ty_rec = -(py_rec).Mag2();
        double qy_rec = sqrt(ty_rec);
                                
        // define x and z-components 
        TLorentzVector px_rec = T_rec - py_rec; 
        TLorentzVector pz_rec(0,0,px_rec.Z(),px_rec.E());
        TLorentzVector px_rec_true = T_rec - py_rec - pz_rec;
        double tz_rec = -(pz_rec).Mag2();
        double tx_rec = -(px_rec_true).Mag2();
        double qx_rec = sqrt(tx_rec);
        double qz_rec = sqrt(tz_rec);

        if (tx_rec < 0 || ty_rec < 0) continue;

        // plot t = tx+ty+t//
        t_total_rec = tx_rec+ty_rec+tz_rec;
        TLorentzVector q_rec_vect(qx_rec,qy_rec,sqrt(tz_rec),0);
        double q_rec = q_rec_vect.Mag2();
        have_reco = true;

        bool veto_clust = (event.zdc_clusters.size() > 0);
        bool veto_RP   = (event.hit_rp.size() > 0);
        bool veto_OMD  = (event.hit_omd.size() > 0);
        bool pass_pi12 = false;

        h_t_REC_2d_wRES_cut->Fill(qx_rec,qy_rec);
        h_t_REC_wRES_cut->Fill(t_total_rec);

        TVector3 tvector(tx_rec,ty_rec,tz_rec);
        TVector3 qvector(qx_rec,qy_rec,qz_rec);
        h_t_phi->Fill(tvector.Phi());
        h_q_phi->Fill(qvector.Phi());
        h_t_theta->Fill(tvector.Theta());
        h_q_theta->Fill(qvector.Theta());
        h_t_perp->Fill(tvector.Pt());
        h_q_perp->Fill(qvector.Pt());
        h_q_mag->Fill(qvector.Mag2());
        h_q_x->Fill(qx_rec);
        h_q_y->Fill(qy_rec);

        double t_total_rec_perp = tvector.Pt();
        double q_perp = qvector.Pt();
        
        // apply cut with // subtracted out of t
        if (qy_rec == 0) continue;
        double theta_rec = atan(fabs(qx_rec)/fabs(qy_rec));
        if(fabs(theta_rec)<PI/2)
        {
            h_t_REC_wRES_cut_pi2->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi2->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/3)
        {
            h_t_REC_wRES_cut_pi3->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi3->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/4)
        {
            h_t_REC_wRES_cut_pi4->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi4->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/6)
        {
            h_t_REC_wRES_cut_pi6->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi6->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/9)
        {
            h_t_REC_wRES_cut_pi9->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi9->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/12)
        {   
            pass_pi12 = true;
            nPassRec_pi12++;
            h_passPi12->Fill(0.5);
            
            h_theta_wedge->Fill(theta_rec);
            h_res_qy->Fill(qy-qy_rec);
            h_res_q->Fill(q_MC-q_rec);
            h_res_qperp->Fill(q_MC_perp-q_perp);
            h_res_eScat->Fill(scatMC.Pt()-scatClusEREC.Pt());
            h_res_vm->Fill(vmMC.Pt()-vmREC.Pt());
            h_res_omega->Fill(theta-theta_rec);
            h_res_nHat->Fill(nHat_MC-nHat_reco);
            h_res_vm_phi->Fill(vmMC.Phi()-vmREC.Phi());
            h_res_vm_theta->Fill(vmMC.Theta()-vmREC.Theta());
            h_res_vm_mag->Fill((vmMC-vmREC).Mag2());
                                    
            h_t_REC_wRES_cut_pi12->Fill(t_total_rec);                                   
            h_t_REC_2d_wRES_cut_pi12->Fill(qx_rec,qy_rec);
            h_t_response_cut->Fill(t_total_rec,t_MC);
            double restcut = (t_MC-t_total_rec)/t_MC;
            h_t_res->Fill(t_MC,restcut);
            h_t_migration->Fill(t_MC,t_total_rec);
            double res_percent_proj_pi12 = (t_MC-t_total_rec)/t_MC;
            h_t_res_proj_percent_pi12->Fill(t_MC,res_percent_proj_pi12);
            double res_t_proj_pi12 = (t_MC-t_total_rec);
            h_t_res_proj_12->Fill(t_MC,res_t_proj_pi12);
            h_t_res_counts->Fill((t_total_rec-t_MC)/t_MC);

            h_res_t->Fill(t_MC-t_total_rec);
            h_res_tperp->Fill(t_MC_perp-t_total_rec_perp);

            // Diagnostics
            if (!veto_clust) h_t_pi12_notVetoed_clust->Fill(t_total_rec); 
            else h_t_pi12_vetoed_clust->Fill(t_total_rec);
            if (!veto_RP) h_t_pi12_notVetoed_RP->Fill(t_total_rec); 
            else h_t_pi12_vetoed_RP->Fill(t_total_rec); 
            if (!veto_OMD) h_t_pi12_notVetoed_OMD->Fill(t_total_rec); 
            else h_t_pi12_vetoed_OMD->Fill(t_total_rec); 
            if (!pass_etaHFS) h_t_pi12_vetoed_etaHFS->Fill(t_total_rec); 
            else h_t_pi12_notVetoed_etaHFS->Fill(t_total_rec);  

            h_REC_neutrons_vs_t_2D->Fill(nNeutrons_RECO, nNeutrons);
            h_pass_pi12_flag->Fill(1); 
            h_mult_neut_pi12->Fill(nNeutrons_RECO);
            h_mult_zdc_hits_pi12->Fill(event.zdc_REC_hits.size());
            h_mult_zdc_clusters_pi12->Fill(event.zdc_clusters.size());
            h_mult_rp_hits_pi12->Fill(event.hit_rp.size());
            h_mult_omd_hits_pi12->Fill(event.hit_omd.size());
            h_mult_eemc_clusters_pi12->Fill(event.clusters_eemc.size());
            h_mult_hfs_pi12->Fill(total_hfs_particles);
            h_mult_electron_pi12->Fill(scat_electron);
            h_mult_kaon_pi12->Fill(total_kaons);
            double HFS_plus_kaon = total_hfs_particles+total_kaons;
            h_mult_HFSplusKaon->Fill(HFS_plus_kaon);
        }
        if (!pass_pi12) h_pass_pi12_flag->Fill(0);
        if(fabs(theta_rec)<PI/16)
        {
            h_t_REC_wRES_cut_pi16->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi16->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/20)
        {
            h_t_REC_wRES_cut_pi20->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi20->Fill(qx_rec,qy_rec);
        }
        if(fabs(theta_rec)<PI/24)
        {
            h_t_REC_wRES_cut_pi24->Fill(t_total_rec);
            h_t_REC_2d_wRES_cut_pi24->Fill(qx_rec,qy_rec);
        }                       
    }//end while event loop

    cout << "Total events: " << nAll << endl;
    cout << "" << endl;
    cout << " ------------- REC Level Selections ---------------" << endl;
    cout << "Events killed by scatClusEREC==0:                  " << nAfterScatClus << endl;
    cout << "Events killed by 15<E-pz<25:                       " << nEpz << endl;
    cout << "Events killed by 0.8<E/p<1.2:                      " << nEoP << endl;
    cout << "Events killed by electron requirement:             " << nAfterElectron << endl;
    cout << "Events killed by incoherent REC:                   " << nAfterIncoherentRec << endl;
    //cout << "Events killed by HFS==2:                           " << nAfterHFS2 << endl;
    //cout << "Events killed by kaons==2:                         " << nAfterKaons2 << endl;
    cout << "Events killed by rec VM==0:                        " << nAfterRecVM << endl;
    cout << "Events killed by rec VM mass:                      " << nPassMassRange << endl;
    cout << "Events killed by 1<Q2<10:                          " << nQ2 << endl;
    cout << "Events killed by 0.01<y<0.85:                      " << ny << endl;
    cout << "" << endl;
    cout << "Events remaining after Q2 and y selection:         " << nPassEventCuts << endl;
    cout << "Events killed by ZDC:                              " << nPassZDC << endl;
    cout << "Events killed by RP:                               " << nPassRP  << endl;
    cout << "Events killed by OMD:                              " << nPassOMD << endl;
    cout << "Events remaining in #theta_{max} |t| Distribution: " << nPassRec_pi12 << endl;
    cout << "" << endl;

    output->Write();
    output->Close();

    return 0;
}