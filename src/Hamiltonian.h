#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);


class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__, MFParams& MFParams__ )
        :Parameters_(Parameters__),Coordinates_(Coordinates__),MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void Check_up_down_symmetry();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double filling);    //::DONE

    double TotalDensity();   //::DONE
    double E_QM();   //::DONE

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    int lx_, ly_, ns_, orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

    double HS_factor;

};



double Hamiltonian::chemicalpotential(double muin,double filling){
    double mu_out;
    double n1,N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.05*(eigs_[nstate-1] - eigs_[0])/nstate;
    N=filling*double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged=false;


    if(1==1){
        for(int i=0;i<50000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_out)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                mu_out += (N-n1)*dMubydN;
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;

            }
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<endl;
        }

    }


    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if(1==2){
        mu1=eigs_[0];
        mu2=eigs_[nstate-1];
        for(int i=0;i<40000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_temp)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                if(n1 >N){
                    mu2=mu_temp;
                    mu_temp=0.5*(mu1 + mu_temp);
                }
                else{
                    mu1=mu_temp;
                    mu_temp=0.5*(mu2 + mu_temp);
                }

            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if(!converged){
            //cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }

    return mu_out;
} // ----------


void Hamiltonian::Initialize(){


    int ns =(Parameters_.lx_cluster)*(Parameters_.ly_cluster);


    ly_=Parameters_.ly;
    lx_=Parameters_.lx;
    ns_=Parameters_.ns;


    int space=2*ns_;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);

} // ----------

double Hamiltonian::TotalDensity(){
    double n1=0.0;
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return n1;
} // ----------



double Hamiltonian::E_QM(){
    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;
} // ----------



double Hamiltonian::GetCLEnergy(){

    double EClassical;
    int neigh_site;
    double local_den_val_i;
    double sz_val_i, sx_val_i, sy_val_i;
    double local_den_val_j;
    double sz_val_j, sx_val_j, sy_val_j;


    EClassical=0.0;
    for(int i=0;i<ns_;i++) {  // For each site
        local_den_val_i=MFParams_.Local_density(Coordinates_.indx(i),Coordinates_.indy(i));
        sz_val_i = MFParams_.Sz(Coordinates_.indx(i),Coordinates_.indy(i));
        sx_val_i = MFParams_.Sx(Coordinates_.indx(i),Coordinates_.indy(i));
        sy_val_i = MFParams_.Sy(Coordinates_.indx(i),Coordinates_.indy(i));

        // +x neighbour, +y neighbour
        for(int neigh_no=0;neigh_no<4;neigh_no++){
            if(neigh_no==0 || neigh_no==2){

                neigh_site = Coordinates_.neigh(i,neigh_no);

                local_den_val_j=MFParams_.Local_density(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sz_val_j = MFParams_.Sz(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sx_val_j = MFParams_.Sx(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sy_val_j = MFParams_.Sy(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));



                EClassical += -1.0*Parameters_.J_Exchange*( sz_val_i*sz_val_j + sx_val_i*sx_val_j  + sy_val_i*sy_val_j );
                EClassical +=  0.25*Parameters_.J_Exchange*(local_den_val_i*local_den_val_j);

            }
        }


    }

    return EClassical;
} // ----------



void Hamiltonian::InteractionsCreate(){

    int space=2*ns_;
    int neigh_site;
    int x_ng, y_ng;
    int a;
    double local_den_val;
    double sz_val, sx_val, sy_val;
    complex<double> splus_val, sminus_val;

    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            Ham_(i,j)=HTB_(i,j);
        }
    }

    for(int i=0;i<ns_;i++) {  // For each site

        // +\-x neighbour, +\-y neighbour
        for(int neigh_no=0;neigh_no<4;neigh_no++){

            if( (neigh_no<2 && Coordinates_.lx_!=1)
                    ||
                    (neigh_no>=2 && Coordinates_.ly_!=1)
                    ){
                neigh_site = Coordinates_.neigh(i,neigh_no);

                local_den_val=MFParams_.Local_density(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sz_val = MFParams_.Sz(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sx_val = MFParams_.Sx(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                sy_val = MFParams_.Sy(Coordinates_.indx(neigh_site),Coordinates_.indy(neigh_site));
                splus_val=one_complex*(sx_val) + iota_complex*(sy_val);
                sminus_val=one_complex*(sx_val) - iota_complex*(sy_val);

                Ham_(i,i) += (0.5*Parameters_.J_Exchange*(sz_val)) + (-0.25*Parameters_.J_Exchange*(local_den_val));
                Ham_(i+ns_,i+ns_) += (-0.5*Parameters_.J_Exchange*(sz_val)) + (-0.25*Parameters_.J_Exchange*(local_den_val));

                Ham_(i,i+ns_) +=  0.5*Parameters_.J_Exchange*sminus_val;
                Ham_(i+ns_,i) +=  0.5*Parameters_.J_Exchange*splus_val;
            }


        }



    }



} // ----------


void Hamiltonian::Check_up_down_symmetry()

{
    complex<double> temp(0,0);
    complex<double> temp2;

    for(int i=0;i<ns_;i++) {
        for(int j=0;j<ns_;j++) {
            temp2 = Ham_(i,j) - Ham_(i+ns_,j+ns_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp +=temp2*conj(temp2);
        }
    }

    cout<<"Assymetry in up-down sector: "<<temp<<endl;
}




void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<2*ns_;i++) {
        for(int j=0;j<2*ns_;j++) {
            if(Ham_(i,j) != conj(Ham_(j,i))) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(Ham_(i,j)==conj(Ham_(j,i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}


void Hamiltonian::HTBCreate(){

    int mx=Parameters_.TBC_mx;
    int my=Parameters_.TBC_my;
    complex<double> phasex, phasey;
    int l,m,a,b;

    HTB_.fill(0.0);


    for(l=0;l<ns_;l++) {

        // * +x direction Neighbor
        if(Coordinates_.lx_ != 1){
            if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
                phasex=exp(iota_complex*2.0*(1.0*mx)*PI/(1.0*Parameters_.TBC_cellsX));
                phasey=one_complex;
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,0);
            for(int spin=0;spin<2;spin++) {

                a = l + ns_*spin;
                b = m + ns_*spin;
                assert (a!=b);
                if(a!=b){
                    HTB_(a,b)=complex<double>(1.0*Parameters_.t_hopping,0.0)*phasex;
                    HTB_(b,a)=conj(HTB_(a,b));
                }


            }
        }


        // * +y direction Neighbor
        if(Coordinates_.ly_ != 1){
            if(Coordinates_.indy(l)==(Coordinates_.ly_ -1)){
                phasex=one_complex;
                phasey=exp(iota_complex*2.0*(1.0*my)*PI/(1.0*Parameters_.TBC_cellsY));
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,2);
            for(int spin=0;spin<2;spin++) {

                a = l + ns_*spin;
                b = m + ns_*spin;
                assert (a!=b);
                if(a!=b){

                    HTB_(a,b)=complex<double>(1.0*Parameters_.t_hopping,0.0)*phasey;

                    HTB_(b,a)=conj(HTB_(a,b));
                }

            }
        }

    }

    //On-Site potential
    for(int i=0;i<ns_;i++) {
        for(int spin=0;spin<2;spin++) {
            a = i + ns_*spin;
            HTB_(a,a)=complex<double>(1.0,0.0)*MFParams_.Disorder(Coordinates_.indx(i), Coordinates_.indy(i));
        }
    }


    //HTB_.print();
} // ----------



void Hamiltonian::Hoppings(){
    //DOES SOMETHING EXACT i.e NOTHING :)

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
