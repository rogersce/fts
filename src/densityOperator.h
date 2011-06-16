//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#ifndef DENSITYOPERATOR_H_
#define DENSITYOPERATOR_H_

#include"psmFpeSolver.h"
#include<blitz/array.h>

class densityOperator
{
public:

   struct Params
   {
      double b , L;
      int Nx, Ny, Nz, Ns;
   } p;

   densityOperator(Params _p);
   ~densityOperator();

   blitz::Array<double,3> solve(const blitz::Array<double,3>* ext_field = 0);
   
private:

   void initialize();

   psmFpeSolver psm;

   blitz::Array<double,3>* q;
};

#endif
