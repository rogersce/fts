#ifndef UTILITY_H_
#define UTILITY_H_

#include<silo.h>
#include<blitz/array.h>
#include"MersenneTwister.h"
#include<cnpy.h>

void timer();
int type2Silo(const std::type_info& t);
int box2Silo(DBfile* db, const blitz::TinyVector<double,3>& box);
int silo2Box(DBfile* db, blitz::TinyVector<double,3>& box);
int quadMesh2Silo(DBfile* db, unsigned int Nx, unsigned int Ny, unsigned int Nz, const char* meshname, float scale = 1);
void MTRand2Npy(std::string filename, const MTRand& rng, std::string mode = "w");
void npy2MTRand(std::string filename, MTRand& rng);

template<typename T>
std::string tostring(T i, int pad = 0, char padval = ' ')
{
	std::stringstream s;
	s << std::setfill(padval) << std::setw(pad) << i;
	return s.str();
}

template<typename F>
int writeSilo(std::string fname, std::string mode, F f)
{
   DBfile* db = 0;
   if(mode == "w") db = DBCreate(fname.c_str(),DB_CLOBBER,DB_LOCAL,NULL,DB_PDB);
   else if(mode == "a")
   {
      db = DBOpen(fname.c_str(),DB_PDB,DB_APPEND);
      if(!db) db = DBCreate(fname.c_str(), DB_CLOBBER,DB_LOCAL,NULL,DB_PDB);
   }
   else
   {
      std::cout<<"error: unknown mode passed for saveSilo! \n";
      return -1;
   }

   int err = f(db);

   DBClose(db);
   return err;
}

template<typename dtype>
int quadVar2Silo(DBfile* db, const dtype* data, unsigned int Nx, unsigned int Ny, unsigned int Nz, const char* varname, const char* meshname)
{
   int SILO_TYPE = type2Silo(typeid(dtype));  
   int dims[3] = {Nx,Ny,Nz};
   int err = DBPutQuadvar1(db,varname,meshname,(void*)data,dims,3,NULL,0,SILO_TYPE,DB_NODECENT,NULL);

   return err;
}

template<typename dtype>
int silo2QuadVar(DBfile* db, const char* meshname, dtype** data, unsigned int* Nx, unsigned int* Ny, unsigned int* Nz)
{
   DBquadvar* mesh = DBGetQuadvar(db,meshname);
   
   if(!mesh)
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<mesh<<" cannot be found!\n";
      return -1;
   }

   if(mesh->ndims != 3)
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<meshname<<" has dimensionality "<<mesh->ndims<<" (3 expected)!\n";
      return -1;
   }
   
   if(mesh->datatype != type2Silo(typeid(dtype)))
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<meshname<<" does not have the same dtype as input array! \n";
      return -1;
   }

   *Nx = mesh->dims[0];
   *Ny = mesh->dims[1];
   *Nz = mesh->dims[2];

   *data = new dtype[mesh->nels];

   memcpy(*data,(dtype*) mesh->vals[0],sizeof(dtype)*mesh->nels);

   DBFreeQuadvar(mesh);
   
   return 0;
}

template<typename dtype, int Ndims>
int silo2PointMesh(DBfile* db, const char* meshname, blitz::TinyVector<dtype,Ndims>** data, int* N)
{
   DBpointmesh* ptmesh = DBGetPointmesh(db,meshname);

   if(!ptmesh)
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<meshname<<" of type Pointmesh cannot be found!\n";
      return -1;
   }

   if(ptmesh->ndims != Ndims)
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<meshname<<" has dimensionality "<<ptmesh->ndims<<" ("<<Ndims<<" expected)!\n";
      return -1;
   }
   
   if(ptmesh->datatype != type2Silo(typeid(dtype)))
   {
      std::cout<<"error: "<<db->pub.name<<"/"<<meshname<<" does not have the same dtype as input array! \n";
      return -1;
   }

   *N = ptmesh->nels;

   *data = new blitz::TinyVector<double,Ndims>[*N];

   for(int i = 0;i < Ndims;i++)
   {
      dtype* ptr = (dtype*) ptmesh->coords[i];
      for(int j = 0;j < *N;j++) (*data)[j][Ndims-1-i] = ptr[j];
   }

   DBFreePointmesh(ptmesh);

   return 0;
}

#endif
