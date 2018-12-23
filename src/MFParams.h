#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:

    // Define Fields
    Matrix<double> Sz, Sx, Sy;
    Matrix<double> Local_density;
    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator1__ , mt19937_64& Generator2__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }


    double random1();
    double random2();
    void initialize();


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_,ly_,ns_;

    uniform_real_distribution<double> dis1_;//for random fields
    uniform_real_distribution<double> dis2_;//for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};



double MFParams::random1(){

    return dis1_(Generator1_);

}

double MFParams::random2(){

    return dis2_(Generator2_);

}


void MFParams::initialize(){

    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);

    Disorder.resize(lx_,ly_);

    Local_density.resize(lx_,ly_);

    Sz.resize(lx_,ly_);
    Sx.resize(lx_,ly_);
    Sy.resize(lx_,ly_);

    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#ix   iy    Dis[ix,iy]"<<endl;

    ofstream Initial_OrderParams_file("Initial_OrderParams_values");

    if(!Parameters_.Read_OPs){
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                Sz(i,j)=random1();
                Sx(i,j)=random1();
                Sy(i,j)=random1();
                Local_density(i,j)=random1();
            }
        }
        Initial_OrderParams_file<<"#seed="<<Parameters_.RandomSeed<<
                                  " for mt19937_64 Generator is used"<<endl;
    }
    else{
        //output_Local_Sx.txt
        //output_Local_Sy.txt
        //output_Local_Sz.txt
        //output_Local_density_obs.txt
        vector<string> OPstring;
        OPstring.clear();
        OPstring.push_back("Sz");OPstring.push_back("Sx");
        OPstring.push_back("Sy");OPstring.push_back("density_obs");

        for(int op_no=0;op_no<OPstring.size();op_no++){
            string fl_initial_OP_in = Parameters_.File_OPs_in + OPstring[op_no] + ".txt";
            ifstream file_initial_OP_in(fl_initial_OP_in.c_str());
            string temp1,temp2,temp3;
            int x,y;
            double val;

            file_initial_OP_in>>temp1>>temp2>>temp3;
            //cout<<temp1<<"  "<<temp2<<" "<<temp3<<endl;
            for (int i=0;i<lx_;i++){
                for (int j=0;j<ly_;j++){
                    file_initial_OP_in>>x>>y>>val;
                    //cout<<ix<<"\t"<<iy<<"\t"<<x<<"\t"<<y<<"\t"<<val<<endl;
                    if(OPstring[op_no]=="Sz"){
                        Sz(i,j)=val;
                    }
                    else if(OPstring[op_no]=="Sx"){
                        Sx(i,j)=val;
                    }
                    else if(OPstring[op_no]=="Sy"){
                        Sy(i,j)=val;
                    }
                    else if(OPstring[op_no]=="density_obs"){
                        Local_density(i,j)=val;
                    }
                }
                //file_in_Sz>>temp;
                //cout<<endl;
            }

        }

        Initial_OrderParams_file<<"#OParams are read from "<<Parameters_.File_OPs_in<<"---.txt files"<<endl;

    }



    Initial_OrderParams_file<<"#ix   iy   Sz(x,y)    Sx(x,y)   Sy(x,y)   Local_density(x,y)"<<endl;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            Initial_OrderParams_file<<i<<setw(15)<<j<<setw(15)<<Sz(i,j)<<setw(15)<<Sx(i,j)
                                   <<setw(15)<<Sy(i,j)<<setw(15)<<Local_density(i,j)<<endl;
        }
        Initial_OrderParams_file<<endl;
    }



    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //RANDOM Disorder
            Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
            Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
        }
        Disorder_conf_file<<endl;
    }


} // ----------

#endif
