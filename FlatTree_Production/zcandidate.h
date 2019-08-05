#include "leptons.h"

class zcandidate
{
    leptons daughter1;
    leptons daughter2;
    double mass;
    int daughter1_idx, daughter2_idx;
    
    public:
    zcandidate() = default;;
    zcandidate(leptons, leptons, double, int, int);
    void setAll(leptons, leptons, double, int, int);
    leptons getDaughter1();
    leptons getDaughter2();
    double getMass();
    int getDaughter1_idx();
    int getDaughter2_idx();
};

inline zcandidate::zcandidate(leptons Daughter1, leptons Daughter2, double Mass, int idx1, int idx2)
{
    daughter1.setAll(Daughter1.getPt(), Daughter1.getEta(), Daughter1.getPhi(), Daughter1.getM(), Daughter1.getQ(), Daughter1.getIsoVar(), Daughter1.getPDGID());
    
    daughter2.setAll(Daughter2.getPt(), Daughter2.getEta(), Daughter2.getPhi(), Daughter2.getM(), Daughter2.getQ(), Daughter2.getIsoVar(), Daughter2.getPDGID());
    
    mass = Mass;
    daughter1_idx = idx1;
    daughter2_idx = idx2;
};

void zcandidate::setAll(leptons Daughter1, leptons Daughter2, double Mass, int idx1, int idx2)
{
    daughter1.setAll(Daughter1.getPt(), Daughter1.getEta(), Daughter1.getPhi(), Daughter1.getM(), Daughter1.getQ(), Daughter1.getIsoVar(), Daughter1.getPDGID());
    
    daughter2.setAll(Daughter2.getPt(), Daughter2.getEta(), Daughter2.getPhi(), Daughter2.getM(), Daughter2.getQ(), Daughter2.getIsoVar(), Daughter2.getPDGID());
    
    mass = Mass;
    daughter1_idx = idx1;
    daughter2_idx = idx2;
};
    
leptons zcandidate::getDaughter1()    {return daughter1;};  
leptons zcandidate::getDaughter2()    {return daughter2;};
double  zcandidate::getMass()         {return mass;};
int     zcandidate::getDaughter1_idx(){return daughter1_idx;};
int     zcandidate::getDaughter2_idx(){return daughter2_idx;};