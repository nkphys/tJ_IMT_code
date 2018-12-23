#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"

#ifndef SelfConsistencyEngine_H
#define SelfConsistencyEngine_H


class SelfConsistencyEngine{
public:
    SelfConsistencyEngine(Parameters& Parameters__, Coordinates& Coordinates__,
                          MFParams& MFParams__, Hamiltonian& Hamiltonian__,
                          Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {

    }

    void RUN_SelfConsistencyEngine();
    double Prob (double muu, double mu_new);
    double ProbCluster (double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_;

};

/*
 * ***********
 *  Functions in Class SelfConsistencyEngine ------
 *  ***********
*/

void SelfConsistencyEngine::RUN_SelfConsistencyEngine(){

    complex<double> zero(0.0,0.0);

    double muu_prev;
    double Curr_QuantE;
    double Prev_QuantE;
    double Curr_ClassicalE;
    double Prev_ClassicalE;

    double temp_=0.0001;
    cout << "Temperature = " << temp_<<" is being done"<<endl;
    Parameters_.temp=temp_;
    Parameters_.beta=double(11604.0/temp_);

    int x,y,act;

    string File_Out_progress;

    double initial_mu_guess;
    int n_states_occupied_zeroT;

    //starting with a random guess
    File_Out_progress = "output.txt";
    ofstream file_out_progress(File_Out_progress.c_str());


    file_out_progress<< "Maximum no of self consistency iterations = "<<Parameters_.IterMax<<"."<<endl;
    file_out_progress<<"Convergence error targetted = "<<Parameters_.Convergence_Error<<endl;

    Parameters_.Dflag='V'; // flag to calculate only Eigenvalue


    file_out_progress<<"Iter"<<setw(15)<<
                       "Error_Sz"<<setw(15)<<"Error_Sx"<<setw(15)<<"Error_Sy"<<setw(17)<<"Error_Local_n"<<setw(17)<<
                       "mu"<<setw(17)<<
                       "E_CL"<<setw(17)<<"E_Quantum"<<endl;



/*
    Prev_ClassicalE = Hamiltonian_.GetCLEnergy();
    Hamiltonian_.InteractionsCreate();
    Hamiltonian_.Diagonalize(Parameters_.Dflag);
    n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*2.0;
    initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
    Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);
    Prev_QuantE = Hamiltonian_.E_QM();
    Hamiltonian_.copy_eigs(1);
    cout<<"Initial Classical Energy = "<<Prev_ClassicalE<<endl;
    cout<<"Initial Quantum Energy = "<<Prev_QuantE<<endl;
    cout<<"Initial mu="<<Parameters_.mus<<endl;
    */


    double Error_Sz=100;
    double Error_Sx=100;
    double Error_Sy=100;
    double Error_Local_n=100;

    for(int count=0;count<Parameters_.IterMax;count++){
        if(
                Error_Sz>Parameters_.Convergence_Error
                ||
                Error_Sx>Parameters_.Convergence_Error
                ||
                Error_Sy>Parameters_.Convergence_Error
                ||
                Error_Local_n>Parameters_.Convergence_Error
                ){


            Curr_ClassicalE = Hamiltonian_.GetCLEnergy();
            Hamiltonian_.InteractionsCreate();
  //         Hamiltonian_.Ham_.print();
            Hamiltonian_.Check_Hermiticity();
            Hamiltonian_.Diagonalize(Parameters_.Dflag);
            n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*2.0;
            initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
            Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);
            Curr_QuantE = Hamiltonian_.E_QM();


            Observables_.Calculate_Local_Density();
            Observables_.Calculate_Local_Szxy();

            Observables_.Get_OrderParameters_diffs();

            Error_Sz=Observables_.Error_Sz_;
            Error_Sx=Observables_.Error_Sx_;
            Error_Sy=Observables_.Error_Sy_;
            Error_Local_n=Observables_.Error_Local_n_;


            file_out_progress<<count<<setw(15)<<
                               Error_Sz<<setw(15)<<Error_Sx<<setw(15)<<Error_Sy<<setw(17)<<Error_Local_n<<setw(17)<<
                               Parameters_.mus<<setw(17)<<
                               Curr_ClassicalE<<setw(17)<<Curr_QuantE<<endl;

            Observables_.Update_OrderParameters();

        }

    }

    string File_Out_Local_Sz = Parameters_.File_OPs_out + "Sz.txt";
    ofstream file_out_Local_Sz(File_Out_Local_Sz.c_str());
    file_out_Local_Sz<<"#ix    iy    Sz"<<endl;

    string File_Out_Local_Sx = Parameters_.File_OPs_out + "Sx.txt";
    ofstream file_out_Local_Sx(File_Out_Local_Sx.c_str());
    file_out_Local_Sx<<"#ix    iy    Sx"<<endl;

    string File_Out_Local_Sy = Parameters_.File_OPs_out + "Sy.txt";
    ofstream file_out_Local_Sy(File_Out_Local_Sy.c_str());
    file_out_Local_Sy<<"#ix    iy    Sy"<<endl;

    string File_Out_Local_density_obs = Parameters_.File_OPs_out + "density_obs.txt";
    ofstream file_out_Local_density_obs(File_Out_Local_density_obs.c_str());
    file_out_Local_density_obs<<"#ix    iy    n"<<endl;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            file_out_Local_Sz<<i<<setw(15)<<j<<setw(15)<<Observables_.Sz_obs(i,j)<<endl;
            file_out_Local_Sx<<i<<setw(15)<<j<<setw(15)<<Observables_.Sx_obs(i,j)<<endl;
            file_out_Local_Sy<<i<<setw(15)<<j<<setw(15)<<Observables_.Sy_obs(i,j)<<endl;
            file_out_Local_density_obs<<i<<setw(15)<<j<<setw(15)<<Observables_.Local_density_obs(i,j)<<endl;
        }
        file_out_Local_Sz<<endl;
        file_out_Local_Sx<<endl;
        file_out_Local_Sy<<endl;
        file_out_Local_density_obs<<endl;
    }


} // ---------



#endif // SelfConsistencyEngine_H
