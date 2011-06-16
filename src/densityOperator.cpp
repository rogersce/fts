//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include"densityOperator.h"
#include<iostream>
#include<cmath>

#define PI 3.14159265

using namespace std;
using namespace blitz;

densityOperator::densityOperator(Params _p) : p(_p), 
   psm(psmFpeSolver::Params(_p.b*_p.b/6.0,1.0/double(_p.Ns),_p.L,_p.Nz,_p.Ny,_p.Nx))
{
   initialize();
}

densityOperator::~densityOperator()
{
   delete[] q;
}

void densityOperator::initialize()
{   
   q = new Array<double,3>[p.Ns+1];
   for(int i = 0;i <= p.Ns; i++) q[i].resize(shape(p.Nz,p.Ny,p.Nx));
}

Array<double,3> densityOperator::solve(const Array<double,3>* ext_field)
{
   //initialize the fpeSolver's propagator
   psm.u = 1;
   
   //solve for the propagator
   q[0] = psm.u;
   for(int s = 0; s < p.Ns; s++)
   {
      psm.update(ext_field);
      q[s+1] = psm.u;
   }
   
   //do Simpson's rule on the propagator solution at each point in space
   Array<double,3> densOp(q[0].shape());
   densOp = (2.0*q[0]*q[p.Ns]);
   for(int s = 1;s < p.Ns; s++) densOp += (s & 1 ? 4.0 : 2.0 )*q[s]*q[p.Ns - s];
   
   double VQw = real(psm.get_uk()(0,0,0));
   densOp *= 1.0/(3.0*p.Ns*VQw);
   
   return densOp;
}
