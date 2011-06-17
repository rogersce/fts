//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include"psmFpeSolver.h"
#include<cmath>

#define PI 3.14159265

using namespace blitz;
using namespace std;

psmFpeSolver::psmFpeSolver(Params _p) : p(_p)
{
   initialize();
}
psmFpeSolver::~psmFpeSolver()
{
   fftw_destroy_plan(forward);
   fftw_destroy_plan(reverse);
   fftw_free(u.data()); 
   fftw_free(u_k.data());
}

void psmFpeSolver::initialize()
{
   int nx = p.Nx/2 + 1;
   TinyVector<double,3> realShape(p.Nz,p.Ny,p.Nx);
   TinyVector<double,3> fourierShape(p.Nz,p.Ny,p.Nx);
   
   u.resize(realShape);
   u_k.resize(fourierShape);
   K2.resize(fourierShape);
   
   int nelsReal = product(realShape);
   int nelsFourier = product(fourierShape);
   double* u_raw = (double*) fftw_malloc(nelsReal*sizeof(double)); 
   fftw_complex* uk_raw = (fftw_complex*) fftw_malloc(nelsFourier*sizeof(fftw_complex));
         
   int _shape[] = {p.Nz,p.Ny,p.Nx};
   forward = fftw_plan_dft_r2c(3,_shape, u_raw, uk_raw, FFTW_MEASURE);
   reverse = fftw_plan_dft_c2r(3,_shape, uk_raw, u_raw, FFTW_MEASURE);
   
   //tricky. wrap blitz Arrays around the raw pointers above.
   //Blitz doesnt like fftw_complex, but likes complex<double>, so reinterpret fftw_complex* as complex<double>*
   u.reference(Array<double,3>(u_raw,realShape,neverDeleteData));
   u_k.reference(Array<complex<double>,3>(reinterpret_cast<complex<double>* >(uk_raw),fourierShape,neverDeleteData));
   
   //fill out K2   
   for(int k = 0; k < p.Nz; k++)
   {
      for(int j = 0;j < p.Ny; j++)
      {
         for(int i = 0;i < nx; i++)
         {
            int jp = j > p.Ny/2 ? -j + p.Ny : j;
            int kp = k > p.Nz/2 ? -k + p.Nz : k;
            K2(k,j,i) = 4.0*PI*PI*(i*i + jp*jp + kp*kp)/(p.L*p.L);
         }
      }
   }
   
   norm = 1.0 / double(nelsReal);
}

void psmFpeSolver::update(const blitz::Array<double,3>* ext_field)
{
   //real space
   if(ext_field) u *= exp(-0.5*p.dt*(*ext_field));
   
   //fourier space
   fftw_execute(forward);
   u_k *= exp(-p.D*p.dt*K2);
   
   //real space
   fftw_execute(reverse);
   u *= norm;
   
   if(ext_field) u *= exp(-0.5*p.dt*(*ext_field));
}

blitz::Array<complex<double>,3>& psmFpeSolver::get_uk()
{
   fftw_execute(forward);
   return u_k;
}
