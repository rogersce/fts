//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#ifndef PSMFPESOLVER
#define PSMFPESOLVER

#include<blitz/array.h>
#include<complex>
#include<fftw3.h>

class psmFpeSolver
{

public:

   struct Params
   {
      Params() {}
      Params(double _D, double _dt, double _L, int _Nz, int _Ny, int _Nx) : D(_D), dt(_dt), L(_L), Nz(_Nz), Ny(_Ny), Nx(_Nx) {}
      double D, dt, L;
      int Nx , Ny , Nz;
   } p;

   psmFpeSolver(Params _p);
   ~psmFpeSolver();

   void update(const blitz::Array<double,3>* ext_field = 0);
   blitz::Array<std::complex<double>,3>& get_uk();
   blitz::Array<double,3> u;   

private:

   void initialize();

   blitz::Array<std::complex<double>,3> u_k;
   blitz::Array<double,3> K2;
   
   fftw_plan forward, reverse;
   double norm;
};

#endif
