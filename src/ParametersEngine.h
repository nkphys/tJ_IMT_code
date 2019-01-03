#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int lx, ly, ns, IterMax, RandomSeed;
    double Convergence_Error;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus,Fill,pi;
    double J_Exchange;
    double U_onsite;
    double Disorder_Strength, RandomDisorderSeed;

    bool Read_OPs;
    string File_OPs_in, File_OPs_out;

    double t_hopping;

    bool Simple_Mixing;
    bool Broyden_Mixing;
    double alpha_n, alpha_Sz, alpha_Sx, alpha_Sy;

    double temp,beta,Eav,maxmoment;

    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){

    t_hopping=1.0;

    double Simple_Mixing_double, Broyden_Mixing_double;
    double Read_OPs_double;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;


    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));

    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;

    Fill = matchstring(inputfile_,"Fill");
    cout << "TotalNumberOfParticles = "<< ns*Fill*2.0 << endl;

    IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
    Convergence_Error=matchstring(inputfile_,"Convergence_Error");
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    J_Exchange = matchstring(inputfile_,"J_Exchange");
    U_onsite = matchstring(inputfile_,"U_Onsite");

    alpha_n = matchstring(inputfile_,"alpha_n");
    alpha_Sz = matchstring(inputfile_,"alpha_Sz");
    alpha_Sx = matchstring(inputfile_,"alpha_Sx");
    alpha_Sy = matchstring(inputfile_,"alpha_Sy");


    Dflag = 'N';

    Simple_Mixing_double=double(matchstring(inputfile_,"Simple_Mixing"));
    Broyden_Mixing_double=double(matchstring(inputfile_,"Broyden_Mixing"));

    if(Broyden_Mixing_double==1.0){
        Broyden_Mixing=true;
        Simple_Mixing=false;

    }
    else if(Broyden_Mixing_double==0.0){
        Broyden_Mixing=false;
        Simple_Mixing=true;
    }
    else{
        cout<<"Broyden_Mixing should be either 0(false) or 1(true)"<<endl;
        assert((Broyden_Mixing_double==0.0) || (Broyden_Mixing_double==1.0));
    }



    Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
    if(Read_OPs_double==1.0){
        Read_OPs=true;
    }
    else{
        Read_OPs=false;
    }


    File_OPs_in=matchstring2(inputfile_,"Read_initial_OPvalues_file");
    File_OPs_out=matchstring2(inputfile_,"Write_Final_OPvalues_file");



    pi=4.00*atan(double(1.0));
    Eav=0.0;

    mus=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



