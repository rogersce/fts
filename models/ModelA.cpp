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
double b = 0.1;
double L = 10;
int Ns = 64;
int Nx = 32;
int Ny = 32;
int Nz = 32;
double Npol = 100;
double u0 = 1.0;

MTRand rng;
densityOperator* densOp;
Array<complex<double>,3> mu , dens;

void cl_parse(int argc, char * argv[])
{
    for(int i=1;i<argc;i+=2)
    {
        if(strcmp(argv[i],"-seed")==0)   seed = atoi(argv[i+1]);

        else if(strcmp(argv[i],"-b")==0)     b = atof(argv[i+1]);
        else if(strcmp(argv[i],"-L")==0)     L = atof(argv[i+1]);
        else if(strcmp(argv[i],"-Ns")==0)    Ns = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Nx")==0)    Nx = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Ny")==0)    Ny = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Nz")==0)    Nz = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-Npol")==0)  Npol = atoi(argv[i+1]);
        else if(strcmp(argv[i],"-u0")==0)    u0 = atof(argv[i+1]);
        else if(strcmp(argv[i],"-nsteps")==0)    nsteps = atoi(argv[i+1]);
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
    cout<<"b:\t"<<b<<"\n";
    cout<<"L:\t"<<L<<"\n";
    cout<<"Ns:\t"<<Ns<<"\n";
    cout<<"Nx:\t"<<Nx<<"\n";
    cout<<"Ny:\t"<<Ny<<"\n";
    cout<<"Nz:\t"<<Nz<<"\n";
    cout<<"Npol:\t"<<Npol<<"\n";
    cout<<"u0:\t"<<u0<<"\n";
    cout<<"dt:\t"<<dt<<"\n";
    cout<<endl<<endl;
}

void initRandom()
{
    complex<double> total = 0;

    for(blitz::Array<complex<double>,3>::iterator it = mu.begin(); it != mu.end(); ++it)
    {
        double r1 = 2.0*(rng.rand53() - 0.5);
        double r2 = 2.0*(rng.rand53() - 0.5);
        
        *it = complex<double>(r1,r2);
        total += *it;
    }

    mu -= total / double(Nx*Ny*Nz);
}

int saveSilo(DBfile* db)
{
    //save mu and densOp here
    int err = DBMkDir(db,"/ModelA");
    err += DBMkDir(db,"/ModelA/data");
    err += DBSetDir(db,"/ModelA/data");
    err += quadMesh2Silo(db, Nx, Ny, Nz, "Mesh",L/double(Nx));
    
    Array<double,3> data(dens.shape());

    data = real(dens);
    err += quadVar2Silo(db,data.data(),Nx,Ny,Nz,"real_dens","Mesh");
    
    data = imag(dens);
    err += quadVar2Silo(db,data.data(),Nx,Ny,Nz,"imag_dens","Mesh");

    data = real(mu);
    err += quadVar2Silo(db,data.data(),Nx,Ny,Nz,"real_mu","Mesh");
    
    data = imag(mu);
    err += quadVar2Silo(db,data.data(),Nx,Ny,Nz,"imag_mu","Mesh");
    
    return err;  
}

int main(int argc, char* argv[])
{
    cl_parse(argc,argv);
    print_stats();

    rng.seed(seed);

    densityOperator::Params p;
    p.b = b;
    p.L = L;
    p.Ns = Ns;
    p.Nz = Nz;
    p.Ny = Ny;
    p.Nx = Nx;

    densOp = new densityOperator(p);
    mu.resize(shape(Nz,Ny,Nx));
    dens.resize(shape(Nz,Ny,Nx));
    initRandom();

    int padlength = tostring(nsteps).length() + 3;
    for(long long counter = 0; counter <= nsteps;counter++)
    {
        if(counter % REPORT == 0)
        {
            timer();
            cout<<"% Done: "<<fixed<<setprecision(5)<<counter<<"/"<<nsteps<<"\n";
            cout<<"Total number real particles: "<<sum(real(dens))<<"\n";
            cout<<"Total number imag particles: "<<sum(imag(dens))<<"\n";
        }

        if(counter % PICTURE == 0)
        {
            int err = writeSilo(output_prefix+"/"+tostring(counter,padlength,'0')+".silo","w",&saveSilo);
            if(err) std::cout<<"Warning: writeSilo returned errors!\n";
        }

        dens = Npol*densOp->solve(&mu); 
        mu += dt*(dens - mu/u0);
    }

    delete densOp;
    return 0;
}
