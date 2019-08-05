class leptons
{
    double pt, eta, phi, m, q, isovar;
    int n, pdgid;
    
    public:
    leptons() = default;
    leptons(double, double, double, double, double, double, int);
    void setAll(double, double, double, double, double, double, int);
    void setPt(double);
    void setEta(double);
    void setPhi(double);
    void setM(double);
    void setIsoVar(double);
    void setQ(double);
    void setPDGID(int);
    double getPt();
    double getEta();
    double getPhi();
    double getM();
    double getIsoVar();
    double getQ();
    int getPDGID();
    TLorentzVector get4P();
    static bool ComparePt(leptons, leptons);
    static bool CompareEta(leptons, leptons);
    static bool ComparePhi(leptons, leptons);
    static bool CompareM(leptons, leptons);
    static bool CompareQ(leptons, leptons);
    static bool CompareIsoVar(leptons, leptons);
    static bool ComparePDGID(leptons, leptons);
};

inline leptons::leptons(double Pt, double Eta, double Phi, double M, double Q, double IsoVar, int PDGID)
{
    pt=Pt;eta=Eta;phi=Phi;m=M;q=Q;isovar=IsoVar;pdgid=PDGID;
};
void leptons::setAll(double Pt, double Eta, double Phi, double M, double Q, double IsoVar, int PDGID)
{
    pt=Pt;eta=Eta;phi=Phi;m=M;q=Q;isovar=IsoVar;pdgid=PDGID;
};
void leptons::setPt(double Pt){pt=Pt;}
void leptons::setEta(double Eta){eta=Eta;}
void leptons::setPhi(double Phi){phi=Phi;}
void leptons::setM(double M){m=M;}
void leptons::setQ(double Q){q=Q;}
void leptons::setIsoVar(double IsoVar){isovar=IsoVar;}
void leptons::setPDGID(int PDGID){pdgid=PDGID;}

double leptons::getPt(){return pt;}
double leptons::getEta(){return eta;}
double leptons::getPhi(){return phi;}
double leptons::getM(){return m;}
double leptons::getQ(){return q;}
double leptons::getIsoVar(){return isovar;}
int leptons::getPDGID(){return  pdgid;}

TLorentzVector leptons::get4P()
{
    TLorentzVector lorentz;
    lorentz.SetPtEtaPhiM(pt, eta, phi, m);
    return lorentz;
}

bool leptons::ComparePt(leptons a, leptons b){return a.getPt() < b.getPt();}
bool leptons::CompareEta(leptons a, leptons b){return a.getEta() < b.getEta();}
bool leptons::ComparePhi(leptons a, leptons b){return a.getPhi() < b.getPhi();}
bool leptons::CompareM(leptons a, leptons b){return a.getM() < b.getM();}
bool leptons::CompareQ(leptons a, leptons b){return a.getQ() < b.getQ();}
bool leptons::CompareIsoVar(leptons a, leptons b){return a.getIsoVar() < b.getIsoVar();}
bool leptons::ComparePDGID(leptons a, leptons b){return a.getPDGID() < b.getPDGID();}

/*TTreeReader electrons::openFile(std::string a);
{
    filename=a;
    f = TFile::Open(filename);
    TTreeReader myReader(a, f);
    return myReader;
}*/

/*void electrons::readEvent(std::string a, std::int i, std::int v=0)
{
    int j=0;
    myReader.Restart();
    if(a=="electron")
    {
        TTreeReaderValue<int> Size_R(myReader, "Electron_size");
        TTreeReaderArray<float> PT_R(myReader, "Electron.PT");
        TTreeReaderArray<float> Eta_R(myReader, "Electron.Eta");
        TTreeReaderArray<float> Phi_R(myReader, "Electron.Phi");
        TTreeReaderArray<int> q_R(myReader, "Electron.Charge");
        TTreeReaderArray<float> iso_R(myReader, "Electron.IsolationVar");
    }
    if(a=="muon")
    {
        TTreeReaderValue<int> Size_R(myReader, "Muon_size");
        TTreeReaderArray<float> PT_R(myReader, "Muon.PT");
        TTreeReaderArray<float> Eta_R(myReader, "Muon.Eta");
        TTreeReaderArray<float> Phi_R(myReader, "Muon.Phi");
        TTreeReaderArray<int> q_R(myReader, "Muon.Charge");
        TTreeReaderArray<float> iso_R(myReader, "Muon.IsolationVar");
    }
    else{cout<<"Please provide either electron or muon as the first argument."<<endl;}
    while(myReader.next() && j<i-1) //I need to find a way to just access the i-th entry directly.
    {
        ++j;
    }
    myReader.next();
    n=*Size_R;
    if(v==1 && n==0){cout<<"This event produced no "<<a<<"s."<<endl;}
    if(n>0)
    {
        TLorentzVector lorentz;
        for(int k=0;k<n;++k)
        {
            lorentz.SetPtEtaPhiM(TPT_R[k],Eta_R[k],Phi_R[k], m_e);
            pt.push_back(lorentz);
            vector<double> values(q_R[k],iso_R[k]);
            v.push_back(values);
            values.clear();
            //lorentz.clear(); This doesn't work with the lorentz vector class, I'm trying to figure out an alternative.
        }
    }
}*/

//int electrons::n(){return n;} 

//double electrons::getP(int i=0){return math::sqrt(pt[i].Px)*(pt[i].Px)+pt[i].Py)*(pt[i].Py)+pt[i].Pz)*(pt[i].Pz));}