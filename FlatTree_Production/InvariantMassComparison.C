//Electron has PDGID "11", Muon is "13".

#include "zcandidate.h"

int main()
{
    TFile* f = TFile::Open("tag_1_delphes_events.root");
    f->ls();

    double pt, eta, phi, mass, iso_var;
    int charge, pdgid;
    int charge1, charge2;
    double Z_mass = 91.1876; //(GeV) from PDG page: http://pdg.lbl.gov/2019/listings/rpp2019-list-z-boson.pdf

    TTreeReader myReader("Delphes", f);
    int numberOfEvents = myReader.GetEntries(1);

    // Define the branches to be read

    // for electrons
    TTreeReaderValue<int> eSize_R(myReader, "Electron_size");
    TTreeReaderArray<float> ePt_R(myReader, "Electron.PT");
    TTreeReaderArray<float> eEta_R(myReader, "Electron.Eta");
    TTreeReaderArray<float> ePhi_R(myReader, "Electron.Phi");
    TTreeReaderArray<int> eQ_R(myReader, "Electron.Charge");
    TTreeReaderArray<float> eIsoVar_R(myReader, "Electron.IsolationVar");

    // for muons
    TTreeReaderValue<int> muSize_R(myReader, "Muon_size");
    TTreeReaderArray<float> muPt_R(myReader, "Muon.PT");
    TTreeReaderArray<float> muEta_R(myReader, "Muon.Eta");
    TTreeReaderArray<float> muPhi_R(myReader, "Muon.Phi");
    TTreeReaderArray<int> muQ_R(myReader, "Muon.Charge");
    TTreeReaderArray<float> muIsoVar_R(myReader, "Muon.IsolationVar");

    //for generated particles
    TTreeReaderValue<int> Size_R(myReader, "Particle_size");
    TTreeReaderArray<int> PID_R(myReader, "Particle.PID");
    TTreeReaderArray<int> Status_R(myReader, "Particle.Status");
    TTreeReaderArray<int> Mother_R(myReader, "Particle.M1");
    TTreeReaderArray<int> Q_R(myReader, "Particle.Charge");
    TTreeReaderArray<float> M_R(myReader, "Particle.Mass");
    TTreeReaderArray<float> Pt_R(myReader, "Particle.PT");
    TTreeReaderArray<float> Eta_R(myReader, "Particle.Eta");
    TTreeReaderArray<float> Phi_R(myReader, "Particle.Phi");

    //Initialise containers
    leptons aux_lepton;
    zcandidate aux_zcandidate;

    vector<leptons> v_electrons;
    vector<leptons> v_muons;
    vector<leptons> v_gen;

    vector<zcandidate> event_zcandidates;
    vector<zcandidate> best_zcandidates;
    vector<zcandidate> gen_z_daughters;

    vector<int> indexes;

    TLorentzVector daughter1, daughter2;

    //Histogram to fill with invariant masses.
    TH1F* detected_histo = new TH1F("Detected Histogram", "Detected Z-->l+ l- Invariant Mass; Invariant Mass [GeV]; #Events", 100, 65, 115);
    TH1F* generated_histo = new TH1F("Generated Histogram", "Generated Z-->l+ l- Invariant Mass; Invariant Mass [GeV]; #Events", 100, 65, 115);

    int progress, particleSize, eventSelector;
    progress = 0;

    //*****************************************************************************************************************************************************************************
    //*****************************************************************************************************************************************************************************
    //*****************************************************************************************************************************************************************************
    Int_t nRecLep;
    Int_t nGenLep;

    vector<float> recLepPt;
    vector<float> recLepEta;
    vector<float> recLepPhi;
    vector<float> recLepIsoVar;
    vector<int> recLepCharge;
    vector<int> recLepPDGID;
    float recLepInvMass;

    vector<float> genLepPt;
    vector<float> genLepEta;
    vector<float> genLepPhi;
    vector<float> genLepIsoVar;
    vector<int> genLepCharge;
    vector<int> genLepPDGID;
    float genLepInvMass;

    TFile *file = TFile::Open("Z_Candidates_10k.root", "RECREATE");
    // Create the TTree and branches
    TTree *t = new TTree("candidateTree", "Tree with arrays");

    //Branches for the reconstruction data
    t->Branch("nRecLep", &nRecLep);
    t->Branch("recLepPt", &recLepPt);
    t->Branch("recLepEta", &recLepEta);
    t->Branch("recLepPhi", &recLepPhi);
    t->Branch("recLepIsoVar", &recLepIsoVar);
    t->Branch("recLepCharge", &recLepCharge);
    t->Branch("recLepPDGID", &recLepPDGID);
    t->Branch("recLepInvMass", &recLepInvMass);

    //Branches for the generated data
    t->Branch("nGenLep", &nGenLep);
    t->Branch("genLepPt", &genLepPt);
    t->Branch("genLepEta", &genLepEta);
    t->Branch("genLepPhi", &genLepPhi);
    t->Branch("genLepIsoVar", &genLepIsoVar);
    t->Branch("genLepCharge", &genLepCharge);
    t->Branch("genLepPDGID", &genLepPDGID);
    t->Branch("genLepInvMass", &genLepInvMass);
    //*****************************************************************************************************************************************************************************

    for (int i = 0; i < numberOfEvents; ++i)
        {
            eventSelector = 0;
            //std::cout << "Event No:" << i << std::endl;
            myReader.Next();

            v_electrons.clear();
            v_muons.clear();
            event_zcandidates.clear();
            v_gen.clear();
            indexes.clear();

            particleSize = *Size_R;
            // loop over the electrons (if there are at least 2 of them) and push them into v_electrons container
            if (*eSize_R > 1)
                {
                    for (int i_electron = 0; i_electron < *eSize_R; ++i_electron)
                        {
                            pt      = ePt_R.At(i_electron);
                            eta     = eEta_R.At(i_electron);
                            phi     = ePhi_R.At(i_electron);
                            mass    = 0.;
                            charge  = eQ_R.At(i_electron);
                            iso_var = eIsoVar_R.At(i_electron);
                            pdgid   = charge * (-11);
                            //if(progress<100){cout<<progress<<": "<<eta<<endl;}
                            aux_lepton.setAll(pt, eta, phi, mass, charge, iso_var, pdgid);
                            v_electrons.push_back(aux_lepton);
                        }
                }

            // loop over the muons (if there are at least 2 of them) and push them into v_muons container
            if (*muSize_R > 1)
                {
                    for (int i_muon = 0; i_muon < *muSize_R; ++i_muon)
                        {
                            pt      = muPt_R.At(i_muon);
                            eta     = muEta_R.At(i_muon);
                            phi     = muPhi_R.At(i_muon);
                            mass    = 0.;
                            charge  = muQ_R.At(i_muon);
                            iso_var = muIsoVar_R.At(i_muon);
                            pdgid   = charge * (-13);

                            aux_lepton.setAll(pt, eta, phi, mass, charge, iso_var, pdgid);
                            v_muons.push_back(aux_lepton);
                        }
                }

            // make pair of different electrons - calculate the invariant mass - construct the Z candidate - push Z candidate into event_zcandidates
            for (int i_pair = 0; i_pair < v_electrons.size(); ++i_pair)
                {
                    for (int j_pair = 0; j_pair < v_electrons.size(); ++j_pair)
                        {
                            if (i_pair != j_pair && i_pair <= j_pair) // do not pair lepton with itself and exclude the permutations i.e (1,2) and (2,1) are the same pairs
                                {
                                    charge1 = v_electrons.at(i_pair).getQ();
                                    charge2 = v_electrons.at(j_pair).getQ();

                                    if (charge1 * charge2 < 0)
                                        {
                                            daughter1 = v_electrons.at(i_pair).get4P();
                                            daughter2 = v_electrons.at(j_pair).get4P();
                                            mass      = (daughter1 + daughter2).M();
                                            //std::cout << mass << std::endl;
                                            aux_zcandidate.setAll(v_electrons.at(i_pair), v_electrons.at(j_pair), mass, i_pair, j_pair);
                                            event_zcandidates.push_back(aux_zcandidate);
                                        }

                                }
                        }
                }

            // make pair of different muons - calculate the invariant mass - construct the Z candidate - push Z candidate into event_zcandidates
            for (int i_pair = 0; i_pair < v_muons.size(); ++i_pair)
                {
                    for (int j_pair = 0; j_pair < v_muons.size(); ++j_pair)
                        {
                            if (i_pair != j_pair && i_pair <= j_pair) // do not pair lepton with itself and exclude the permutations i.e (1,2) and (2,1) are the same pairs
                                {
                                    charge1 = v_muons.at(i_pair).getQ();
                                    charge2 = v_muons.at(j_pair).getQ();

                                    if (charge1 * charge2 < 0)
                                        {
                                            daughter1 = v_muons.at(i_pair).get4P();
                                            daughter2 = v_muons.at(j_pair).get4P();
                                            mass      = (daughter1 + daughter2).M();
                                            //std::cout << mass << std::endl;
                                            aux_zcandidate.setAll(v_muons.at(i_pair), v_muons.at(j_pair), mass, i_pair, j_pair);
                                            event_zcandidates.push_back(aux_zcandidate);
                                        }
                                }
                        }
                }

            // sort zcandidate container wrt abs(candidate mass - Z mass) (ascending)
            std::sort(event_zcandidates.begin(), event_zcandidates.end(), [=](zcandidate first, zcandidate second) {return abs(first.getMass() - Z_mass) < abs(second.getMass() - Z_mass);});
            if (event_zcandidates.size() > 0)
                {
                    //*****************************************************************************************************************************************************************************
                    //*****************************************************************************************************************************************************************************
                    //*****************************************************************************************************************************************************************************
                    //Fill a dummy vector with all the leptons with the two daughters coming first
                    v_gen.push_back(event_zcandidates.at(0).getDaughter1());
                    v_gen.push_back(event_zcandidates.at(0).getDaughter2());
                    for (int i = 0; i < v_electrons.size(); ++i)
                        {
                            if (abs(v_electrons.at(i).getPDGID()) != abs(v_gen.at(0).getPDGID())) {v_gen.push_back(v_electrons.at(i));}
                            if (abs(v_electrons.at(i).getPDGID()) == abs(v_gen.at(0).getPDGID()) && i != event_zcandidates.at(0).getDaughter1_idx() && i != event_zcandidates.at(0).getDaughter2_idx()) {v_gen.push_back(v_electrons.at(i));}
                        }
                    for (int i = 0; i < v_muons.size(); ++i)
                        {
                            if (abs(v_muons.at(i).getPDGID()) != abs(v_gen.at(0).getPDGID())) {v_gen.push_back(v_muons.at(i));}
                            if (abs(v_muons.at(i).getPDGID()) == abs(v_gen.at(0).getPDGID()) && i != event_zcandidates.at(0).getDaughter1_idx() && i != event_zcandidates.at(0).getDaughter2_idx()) {v_gen.push_back(v_muons.at(i));}
                        }

                    //Clear the root input vectors
                    recLepPt.clear();
                    recLepEta.clear();
                    recLepPhi.clear();
                    recLepIsoVar.clear();
                    recLepCharge.clear();
                    recLepPDGID.clear();

                    //Push the information from the dummy leptons vector we just filled into the branch vectors for the root file.
                    nRecLep = v_gen.size();
                    for (int i = 0; i < v_gen.size(); ++i)
                        {
                            recLepPt.push_back(v_gen.at(i).getPt());
                            recLepEta.push_back(v_gen.at(i).getEta());
                            recLepPhi.push_back(v_gen.at(i).getPhi());
                            recLepIsoVar.push_back(v_gen.at(i).getIsoVar());
                            recLepCharge.push_back(v_gen.at(i).getQ());
                            recLepPDGID.push_back(v_gen.at(i).getPDGID());
                        }

                    mass = event_zcandidates.at(0).getMass();
                    best_zcandidates.push_back(event_zcandidates.at(0));
                    detected_histo->Fill(best_zcandidates.back().getMass());
                    eventSelector = 1;
                    recLepInvMass = mass;
                    //*****************************************************************************************************************************************************************************
                }

            v_gen.clear();
            //Find the pairs of z daughters in each event
            if (eventSelector == 1)
                {
                    for (int k = 0; k < particleSize; ++k)
                        {
                            //cout<<Status_R.At(k)<<" "<<PID_R.At(k)<<endl;
                            if (Mother_R.At(k) > -1)
                                {
                                    if (PID_R.At(Mother_R.At(k)) == 23 && Status_R.At(k) > 0)
                                        {
                                            if (abs(PID_R.At(k)) == 11 || abs(PID_R.At(k)) == 13)
                                                {
                                                    pt      = Pt_R.At(k);
                                                    eta     = Eta_R.At(k);
                                                    phi     = Phi_R.At(k);
                                                    mass    = 0.;
                                                    charge  = Q_R.At(k);
                                                    iso_var = 0.;
                                                    pdgid   = PID_R.At(k);
                                                    indexes.push_back(k);
                                                    aux_lepton.setAll(pt, eta, phi, mass, charge, iso_var, pdgid);
                                                    v_gen.push_back(aux_lepton);
                                                }
                                        }
                                }
                        }
                    if (indexes.size() > 1)
                        {
                            daughter1 = v_gen.at(0).get4P();
                            daughter2 = v_gen.at(1).get4P();
                            mass = (daughter1 + daughter2).M();
                            aux_zcandidate.setAll(v_gen.at(0), v_gen.at(1), mass, indexes.at(0), indexes.at(1));
                            gen_z_daughters.push_back(aux_zcandidate);
                            if (abs(mass - Z_mass) >= 0) {generated_histo->Fill(mass);}
                            //cout<<i<<" "<<l<<" "<<pdgid<<" "<<Status_R.At(l)<<" "<<Mother_R.At(l)<<" "<<PID_R.At(Mother_R.At(l))<<endl;
                            //*****************************************************************************************************************************************************************************
                            //*****************************************************************************************************************************************************************************
                            //*****************************************************************************************************************************************************************************
                            //Clear the root input vectors
                            genLepPt.clear();
                            genLepEta.clear();
                            genLepPhi.clear();
                            genLepIsoVar.clear();
                            genLepCharge.clear();
                            genLepPDGID.clear();

                            nGenLep = 2;
                            genLepPt.push_back(v_gen.at(0).getPt()); genLepPt.push_back(v_gen.at(1).getPt());
                            genLepEta.push_back(v_gen.at(0).getEta()); genLepEta.push_back(v_gen.at(1).getEta());
                            genLepPhi.push_back(v_gen.at(0).getPhi()); genLepPhi.push_back(v_gen.at(1).getPhi());
                            genLepIsoVar.push_back(v_gen.at(0).getIsoVar()); genLepIsoVar.push_back(v_gen.at(1).getIsoVar());
                            genLepCharge.push_back(v_gen.at(0).getQ()); genLepCharge.push_back(v_gen.at(1).getQ());
                            genLepPDGID.push_back(v_gen.at(0).getPDGID()); genLepPDGID.push_back(v_gen.at(1).getPDGID());
                            genLepInvMass = mass;

                            t->Fill();
                            //*****************************************************************************************************************************************************************************
                        }
                }
            ++progress;
            if (progress % (numberOfEvents / 10) == 0) {cout << (100 * progress) / numberOfEvents << "% of events have been checked." << endl;}

        }
    file->Write();

    delete file;

    auto canvas3 = new TCanvas("canvas3", "canvas title", 400, 400);
    generated_histo->SetLineColor(kBlue);
    detected_histo->SetLineColor(kRed);
    generated_histo->Draw();
    detected_histo->Draw("same");

    auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->AddEntry(generated_histo, "Generated Mass", "l");
    legend->AddEntry(detected_histo, "Detected Mass", "l");
    legend->Draw();
    canvas3->Draw("HIST");

    return 0;
}

void InvariantMassComparison()
{
    main();
}