#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw();
    void Calculate_Nw();
    void Calculate_IPR();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();

    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void Calculate_Local_Density();
    void Calculate_Local_Szxy();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);

    double Omega(int i);

    double BandWidth;
    double nia_t,nib_t,nic_t,n_t;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;

    // Declare Fields
    Matrix<double> Sz_obs, Sx_obs, Sy_obs;
    Matrix<double> Local_density_obs;
    double Error_Sz_;
    double Error_Sx_;
    double Error_Sy_;
    double Error_Local_n_;

    double AVG_Total_Energy, AVG_Total_Energy_sqr;


    //Declare Broyden_Mixing vectors
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;


};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/




void Observables::Get_OrderParameters_diffs(){

    Error_Local_n_=0.0;
    Error_Sz_=0.0;
    Error_Sx_=0.0;
    Error_Sy_=0.0;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){

            Error_Local_n_+= abs(Local_density_obs(i,j) - MFParams_.Local_density(i,j));
            Error_Sz_+= abs(Sz_obs(i,j) - MFParams_.Sz(i,j));
            Error_Sx_+= abs(Sx_obs(i,j) - MFParams_.Sx(i,j));
            Error_Sy_+= abs(Sy_obs(i,j) - MFParams_.Sy(i,j));

        }}


}

void Observables::Update_OrderParameters(int iter){


    //Simple mixing
    double alpha_n=Parameters_.alpha_n;
    double alpha_Sz=Parameters_.alpha_Sz;
    double alpha_Sx=Parameters_.alpha_Sx;
    double alpha_Sy=Parameters_.alpha_Sy;


    if(Parameters_.Simple_Mixing==true){

        if(iter==0){
            cout<<"Using Simple Mixing to gain Self-Consistency"<<endl;
        }

        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){

                MFParams_.Local_density(i,j) = (1-alpha_n)*MFParams_.Local_density(i,j)  + alpha_n*Local_density_obs(i,j);

                MFParams_.Sz(i,j) = (1-alpha_Sz)*MFParams_.Sz(i,j)  + alpha_Sz*Sz_obs(i,j);
                MFParams_.Sx(i,j) = (1-alpha_Sx)*MFParams_.Sx(i,j)  + alpha_Sx*Sx_obs(i,j);
                MFParams_.Sy(i,j) = (1-alpha_Sy)*MFParams_.Sy(i,j)  + alpha_Sy*Sy_obs(i,j);

            }}

    }

    /*
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;
     */

    vector<double> vec_V, vec_U;
    double Denominator_;
    vector<double> vec_L;
    int site;


    if(Parameters_.Broyden_Mixing==true){
        if(iter==0){
            cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;

            //Get Jinv_np1
            for(int i=0;i<4*ns_;i++){
                for(int j=0;j<4*ns_;j++){
                    if(i==j){
                        Jinv_np1[i][j]=-1.0*alpha_n;
                    }
                    else{
                        Jinv_np1[i][j]==0.0;
                    }
                }}

            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<lx_;j++){
                    site=Coordinates_.Nc(i,j);

                    F_n[site]=Local_density_obs(i,j) - MFParams_.Local_density(i,j);
                    F_n[site + ns_]=Sz_obs(i,j) - MFParams_.Sz(i,j);
                    F_n[site + (2*ns_)]=Sx_obs(i,j) - MFParams_.Sx(i,j);
                    F_n[site + (3*ns_)]=Sy_obs(i,j) - MFParams_.Sy(i,j);

                }
            }

            for(int i=0;i<4*ns_;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<4*ns_;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    MFParams_.Local_density(i,j) +=  Delta_x_n[site];
                    MFParams_.Sz(i,j) += Delta_x_n[site+ns_];
                    MFParams_.Sx(i,j) += Delta_x_n[site+(2*ns_)];
                    MFParams_.Sy(i,j) += Delta_x_n[site+(3*ns_)];
                }
            }

            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;


        }

        else{
            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<lx_;j++){
                    site=Coordinates_.Nc(i,j);

                    F_n[site]=Local_density_obs(i,j) - MFParams_.Local_density(i,j);
                    F_n[site + ns_]=Sz_obs(i,j) - MFParams_.Sz(i,j);
                    F_n[site + (2*ns_)]=Sx_obs(i,j) - MFParams_.Sx(i,j);
                    F_n[site + (3*ns_)]=Sy_obs(i,j) - MFParams_.Sy(i,j);

                }
            }

            //Get DeltaF_n
            for (int i=0;i<4*ns_;i++){
                DeltaF_n[i] = F_n[i] - F_nm1[i];
            }

            //Get vec_V = Jinv_n*DeltaF_n
            vec_V.clear();
            vec_V.resize(4*ns_);

            for(int i=0;i<4*ns_;i++){
                vec_V[i] =0.0;
                for(int j=0;j<4*ns_;j++){
            vec_V[i] += Jinv_n[i][j]*DeltaF_n[j];
                }
            }

            //Get vec_U = Delta_x_n^dagg*Jinv_n
            vec_U.clear();
            vec_U.resize(4*ns_);

            for(int i=0;i<4*ns_;i++){
                vec_U[i] =0.0;
                for(int j=0;j<4*ns_;j++){
            vec_U[i] += Delta_x_n[j]*Jinv_n[j][i];
                }
            }

            // Get Denominator_=<Delta_x_n|vec_V>
            Denominator_=0.0;
            for(int i=0;i<4*ns_;i++){
            Denominator_ +=Delta_x_n[i]*vec_V[i];
            }


            //Get vec_L=  Delta_x_n - vec_V;
            vec_L.clear();
            vec_L.resize(4*ns_);
            for(int i=0;i<4*ns_;i++){
                vec_L[i] = Delta_x_n[i] - vec_V[i];
            }


            //Get Mat_Temp [Remember to clear later on];
            Mat_2_doub Mat_Temp;
            Mat_Temp.resize(4*ns_);
            for(int i=0;i<4*ns_;i++){
               Mat_Temp[i].resize(4*ns_);
               for(int j=0;j<4*ns_;j++){
                Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_);
               }
            }


            //Get Jinv_np1
            for(int i=0;i<4*ns_;i++){
               for(int j=0;j<4*ns_;j++){
               Jinv_np1[i][j]  = Jinv_n[i][j]  + Mat_Temp[i][j];
               }
            }



            for(int i=0;i<4*ns_;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<4*ns_;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    MFParams_.Local_density(i,j) +=  Delta_x_n[site];
                    MFParams_.Sz(i,j) += Delta_x_n[site+ns_];
                    MFParams_.Sx(i,j) += Delta_x_n[site+(2*ns_)];
                    MFParams_.Sy(i,j) += Delta_x_n[site+(3*ns_)];
                }
            }

            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;


            //Clear Mat_Temp
            for(int i=0;i<4*ns_;i++){
                Mat_Temp[i].clear();
            }
            Mat_Temp.clear();

        }


    }

}

void Observables::Calculate_Local_Density(){

    int c1, c2;
    int site;

    Local_density_obs.resize(lx_,ly_);


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            Local_density_obs(i,j)=0.0;


            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                c1 = site; c2 = site;
                Local_density_obs(i,j) +=( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                           (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                           ).real();

                c1 = site + ns_; c2 = site + ns_;
                Local_density_obs(i,j) += ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                            ).real();

            }
        }
    }

}



void Observables::Calculate_IPR(){

    double IPR;
    int c1, site;
    double eta = 0.001;
    int n_chosen=Parameters_.ns*Parameters_.Fill*2.0;
    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                //  IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                IPR += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen))*
                        abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));




            }
        }
    }

    cout<<"IPR for state no. "<<Parameters_.ns*Parameters_.Fill*2.0<<" = "<<IPR<<endl;



    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){


                    IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                }
            }
        }
    }

    cout<<"IPR for near mu(with eta = "<<eta<<") = "<<IPR<<endl;


}

void Observables::Calculate_Local_Szxy(){

    int c1, c2;
    int site;
    complex<double> splus_, sminus_;

    Sz_obs.resize(lx_,ly_);
    Sx_obs.resize(lx_,ly_);
    Sy_obs.resize(lx_,ly_);

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            Sz_obs(i,j)=0.0;
            Sx_obs(i,j)=0.0;
            Sy_obs(i,j)=0.0;


            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                c1 = site; c2 = site;
                Sz_obs(i,j) += 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                     (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                     ).real();

                c1 = site + ns_; c2 = site + ns_;
                Sz_obs(i,j) -= 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                     (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                     ).real();


                c1 = site; c2 = site + ns_;
                splus_ = conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                        ;

                c1 = site + ns_; c2 = site;
                sminus_ = conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                        ;

                Sx_obs(i,j) += (0.5*(splus_ + sminus_)).real();
                Sy_obs(i,j) += (-0.5*iota_complex*(splus_ - sminus_)).real();


            }
        }
    }

}

void Observables::Calculate_Akw(){


    //---------Read from input file-----------------------//
    string fileout="Akw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.001;
    omega_min=-1.6;omega_max=2.6;d_omega=0.03;
    //---------------------------------------------------//


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_up_00;
    Mat_3_Complex_doub A_dn_00;
    A_up_00.resize(Parameters_.ns);
    A_dn_00.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        A_up_00[i].resize(Parameters_.ns);
        A_dn_00[i].resize(Parameters_.ns);

        for(int j=0;j<Parameters_.ns;j++){
            A_up_00[i][j].resize(omega_index_max);
            A_dn_00[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_up_00[j][l][omega_ind]=zero_complex;
                A_dn_00[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;

                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                    c1 = l + ns_; c2 = j+ ns_;
                    A_dn_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l; c2 = j;
                    A_up_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                }

                if(j==l){
                    Nup_check += (A_up_00[j][l][omega_ind])*d_omega;
                    Ndn_check += (A_dn_00[j][l][omega_ind])*d_omega;
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_up_00;
    complex<double> temp_dn_00;
    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    double k22_offset=0;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_up_00=zero_complex;
            temp_dn_00=zero_complex;

            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_up_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_00[j][l][omega_ind];

                    temp_dn_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_00[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp_up_00.real()<<"    "<<temp_dn_00.real()<<"    "
                       <<temp_up_00.imag()<<"     "<<temp_dn_00.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



}


void Observables::Calculate_Nw(){

    //---------Read from input file-----------------------//
    string fileout="Nw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=-100;omega_max=100.0;d_omega=0.0001;
    //---------------------------------------------------//

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Nw_out(fileout.c_str());

    double temp_val ;


    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

        temp_val=0.0;

        for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

            temp_val +=  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

        }


        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"     "<<temp_val<<"     "<<endl;

    }

    file_Nw_out<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

}

void Observables::Get_Non_Interacting_dispersion(){

}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

void Observables::DensityOfStates(){
    //-----------Calculate Bandwidth------//
    BandWidth=2.0;
    //-----------------------------------//

} // ----------


void Observables::OccDensity(){

} // ----------


void Observables::TotalOccDensity(){

} // ----------



void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}

void Observables::OccDensity(int tlabel){

} // ----------



void Observables::Initialize(){

    complex<double> zero(0.0,0.0);
    int space=2*ns_;
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);


    F_n.resize(ns_*4); //F_n=x_n_out - x_n_in []
    F_nm1.resize(ns_*4);
    DeltaF_n.resize(ns_*4); //DeltaF_n=F_n - F-nm1;
    Delta_x_n.resize(ns_*4); //Delta_x_n= x_n_in - x_nm1_in;
    Jinv_n.resize(ns_*4);
    Jinv_np1.resize(ns_*4);

    for(int i=0;i<ns_*4;i++){
        Jinv_n[i]. resize(ns_*4);
        Jinv_np1[i]. resize(ns_*4);
    }


} // ----------


double Observables::Omega(int i){
    return -20.0+double(i)*dosincr_;
} // ----------









#endif // OBSERVABLES_H
