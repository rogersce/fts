#include"utility.h"

using namespace blitz;

void timer()
{
    static time_t time_0 = time(NULL);
    time_t elapsed = time(NULL) - time_0;
    time_t days = elapsed / 86400;
    time_t hours = (elapsed % 86400) / 3600;
    time_t minutes = (elapsed % 3600) / 60;
    time_t seconds = (elapsed % 60);
    std::cout<<"Simulation run time: "<<days<<" days, "<<hours<<" hours, "<<minutes<<" minutes, "<<seconds<<" seconds.\n";
}

int type2Silo(const std::type_info& t)
{
    if(t == typeid(float) )     return DB_FLOAT;
    if(t == typeid(double) )    return DB_DOUBLE;
    if(t == typeid(int) )       return DB_INT;
    if(t == typeid(short) )     return DB_SHORT;
    if(t == typeid(long) )      return DB_LONG;
    if(t == typeid(char) )      return DB_CHAR;
    if(t == typeid(long long) ) return DB_LONG_LONG;

    return -1;
}

int quadMesh2Silo(DBfile* db, unsigned int Nx, unsigned int Ny, unsigned int Nz, const char* meshname, float scale)
{   
    int dims[3] = {Nx,Ny,Nz};
    float* arr[3] = {new float[Nx] , new float[Ny] , new float[Nz]};

    for(unsigned int i = 0;i < Nx;i++)  arr[0][i] = scale*i;
    for(unsigned int i = 0;i < Ny;i++)  arr[1][i] = scale*i;
    for(unsigned int i = 0;i < Nz;i++)  arr[2][i] = scale*i;

    int err = DBPutQuadmesh(db,meshname,NULL,arr,dims,3,DB_FLOAT,DB_COLLINEAR,NULL);

    delete[] arr[0];
    delete[] arr[1];
    delete[] arr[2];

    return err;
}

int box2Silo(DBfile* db, const TinyVector<double,3>& box)
{
    double boxcoords[] = { 0, 0,      box[0], box[0], 0,      0,         box[0],    box[0],
        0, box[1], box[1], 0,      0,      box[1],    box[1],    0,
        0, 0,      0,      0,      box[2], box[2],    box[2],    box[2]};

    double* coord_ptrs[] = {boxcoords, boxcoords + 8, boxcoords + 16};
    int err = DBPutPointmesh(db, "Box", 3, coord_ptrs, 8, DB_DOUBLE, NULL);

    return err;
}

int silo2Box(DBfile* db, TinyVector<double,3>& box)
{
    TinyVector<double,3>* b = 0;
    int N;

    int err = silo2PointMesh(db,"Box", &b,&N);

    box = b[3][2] , b[1][1] , b[4][0];
    delete[] b;

    return err;
}

void MTRand2Npy(std::string filename, const MTRand& rng, std::string mode)
{
    int N = rng.SAVE;
    MTRand::uint32* state = new MTRand::uint32[N];
    rng.save(state);

    unsigned int shape[] = {N};
    cnpy::npz_save(filename,"MTRand::state",state,shape,1,mode);

    delete[] state;
}

void npy2MTRand(std::string filename, MTRand& rng) 
{
    cnpy::NpyArray arr = cnpy::npz_load(filename,"MTRand::state");
    assert(arr.word_size == sizeof(MTRand::uint32));
    assert(arr.shape.at(0) == rng.SAVE);
    rng.load((MTRand::uint32*)arr.data);
    delete[] arr.data;
}
