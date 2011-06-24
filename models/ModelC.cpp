//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include<silo.h>
#include<iostream>
#include"MersenneTwister.h"
#include<string>
#include<cstdlib>
#include<cstring>
#include<blitz/array.h>
#include"densityOperator.h"
#include"utility.h"

#define PI 3.14159265

using namespace std;
using namespace blitz;

//Run parameters
int nsteps = 5;
int REPORT = 1;
int PICTURE = 1;
string output_prefix = ".";
double dt = 1;
int seed = time(NULL);

//Model parameters
double bB = 0.1;
double bA = 0.1;
double L = 10;
int Ns = 64;
int Nx = 32;
int Ny = 32;
int Nz = 32;
double chi = 1.0;
double rho_0 = 1.0;
double Na = 1;
double Nb = 1;

MTRand rng;
Array<complex<double>,3> wA(shape(Nz,Ny,Nx)), wB(shape(Nz,Ny,Nx));

void cl_parse(int argc, char * argv[])
{
    for(int i=1;i<argc;i+=2)
    {
        if(strcmp(argv[i],"-seed")==0)   seed = atoi(argv[i+1]);

        else if(strcmp(argv[i],"-bB")==0)     bB = atof(argv[i+1]);
        else if(strcmp(argv[i],"-bA")==0)     bA = atof(argv[i+1]);
        else if(strcmp(argv[i],"-L")==0)     L = atof(argv[i+1]);
        else if(strcmp(argv[i],"-Ns")==0)    Ns = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Nx")==0)    Nx = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Ny")==0)    Ny = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Nz")==0)    Nz = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-chi")==0)    chi = atof(argv[i+1]);
        else if(strcmp(argv[i],"-rho_0")==0)    rho_0 = atof(argv[i+1]);
        else if(strcmp(argv[i],"-nsteps")==0)    nsteps = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Na")==0)    Na = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Nb")==0)    Nb = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-REPORT")==0)    REPORT = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-PICTURE")==0)    PICTURE = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-output_prefix")==0)        output_prefix = argv[i+1];
        else if(strcmp(argv[i],"-dt")==0)    dt = atof(argv[i+1]);

        else cout<<"Warning: Bad input parameter: "<<argv[i]<<endl;
    }

}

void print_stats()
{
    cout<<"nsteps:\t"<<nsteps<<"\n";
    cout<<"PICTURE:\t"<<PICTURE<<"\n";
    cout<<"REPORT:\t"<<REPORT<<"\n";
    cout<<"output_prefix:\t"<<output_prefix<<"\n";

    cout<<"seed:\t"<<seed<<"\n";
    cout<<"bA:\t"<<bA<<"\n";
    cout<<"bB:\t"<<bB<<"\n";
    cout<<"L:\t"<<L<<"\n";
    cout<<"Ns:\t"<<Ns<<"\n";
    cout<<"Nx:\t"<<Nx<<"\n";
    cout<<"Ny:\t"<<Ny<<"\n";
    cout<<"Nz:\t"<<Nz<<"\n";
    cout<<"chi:\t"<<chi<<"\n";
    cout<<"rho_0:\t"<<rho_0<<"\n";
    cout<<"Na:\t"<<Na<<"\n";
    cout<<"Nb:\t"<<Nb<<"\n";
    cout<<"dt:\t"<<dt<<"\n";
    cout<<endl<<endl;
}

void initRandom(Array<complex<double>,3>& arr)
{
    complex<double> total = 0;

    for(blitz::Array<complex<double>,3>::iterator it = arr.begin(); it != arr.end(); ++it)
    {
        double r1 = 2.0*(rng.rand53() - 0.5);
        double r2 = 2.0*(rng.rand53() - 0.5);

        *it = complex<double>(r1,r2);
        total += *it;
    }

    arr -= total / double(Nx*Ny*Nz);
}

int saveSilo(DBfile* db)
{
    //save mu and densOp here
    int err = DBMkDir(db,"/ModelC");
    err += DBMkDir(db,"/ModelC/data");
    err += DBSetDir(db,"/ModelC/data");
    err += quadMesh2Silo(db, Nx, Ny, Nz, "Mesh",L/double(Nx));
    err += quadVar2Silo(db,real(wA).data(),Nx,Ny,Nz,"real_wA","Mesh");
    err += quadVar2Silo(db,real(wB).data(),Nx,Ny,Nz,"real_wB","Mesh");

    return err;  
}

int main(int argc, char* argv[])
{
    cl_parse(argc,argv);
    print_stats();

    rng.seed(seed);

    densityOperator::Params pA;
    pA.b = bA;
    pA.L = L;
    pA.Ns = Ns;
    pA.Nz = Nz;
    pA.Ny = Ny;
    pA.Nx = Nx;
    
    densityOperator::Params pB;
    pB.b = bB;
    pB.L = L;
    pB.Ns = Ns;
    pB.Nz = Nz;
    pB.Ny = Ny;
    pB.Nx = Nx;

    densityOperator densOpA(pA);
    densityOperator densOpB(pB);
    
    initRandom(wA);
    initRandom(wB);

    double chi_inv = 1.0 / chi;

    int padlength = tostring(nsteps).length() + 3;
    for(long long counter = 0; counter <= nsteps;counter++)
    {
        if(counter % REPORT == 0)
        {
            timer();
            cout<<"% Done: "<<fixed<<setprecision(5)<<counter<<"/"<<nsteps<<"\n";
        }

        if(counter % PICTURE == 0)
        {
            int err = writeSilo(output_prefix+"/"+tostring(counter,padlength,'0')+".silo","w",&saveSilo);
            if(err) std::cout<<"Warning: writeSilo returned errors!\n";
        }

        Array<complex<double>,3> factor1((0.5*rho_0/chi_inv)*(wB - wA));
        
        wA += dt*( Na*densOpA.solve(&wA) - factor1 );
        wB += dt*( Nb*densOpB.solve(&wB) + factor1 );
    }

    return 0;
}
