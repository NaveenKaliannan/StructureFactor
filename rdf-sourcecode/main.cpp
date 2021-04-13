/**
 *
 * @author  Naveen Kumar Kaliannan
 * @ Reach me via naveenkumar5892@gmail.com
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <vector> 
#include <numeric>
#include <sstream>  
#include <cmath> 


#define PI 3.14159265
#define deltat 20

using namespace std;
struct Vector
{
  double x,y,z;
};
struct Atom
{
  string gro;
  string symbol;
  double x,y,z;
  double vx,vy,vz;
  uint index;
  string segname;
  uint resid;
  string resname;
  string atomname;
  string atomtype;
  float charge;
  float atomicmass;
};

void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L) ;
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L) ;
double norm(double x,double y,double z){  return sqrt(x*x+y*y+z*z);}
double min_distance(double r, float l) ;
double mindis(double dx,double dy,double dz, const vector<float> & L) ;
void Printrdf(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  float x = 0, y = 0, z = 0, rij = 0, V = L[0] * L[1] * L[2], len = max(L[0], L[1]), RDF_h = 0.1, count = 0 ;
  uint idi = 0, idj = 0, ncell = 0, count_atoms = 0;
  len = max(len, L[2]);
  len = len/2.0;
  uint RDF_size = len*100;
  vector<float> RDF(RDF_size,0.0), rad(RDF_size,0.0);
  for(uint i = 1;i < RDF_size ;i++)
    {
      rad[i]          = i * 0.01;
    }  
  float no_of_residue_first  = 0;
  float no_of_residue_second = 0; 
  cout << "g \u03B1 - \u03B2 (rij)" << endl; ;
  cout << "Enter the first character of \u03B1 species : " ;
  char alpha ;
  cin >>  alpha ; 
  cout << "Enter number of \u03B1 species : " ;
  cin >>  no_of_residue_first ;
  cout << "Enter the first character of \u03B2 species : " ;
  char beta ;
  cin >>  beta ; 
  cout << "Enter number of \u03B2 species : " ;
  cin >>  no_of_residue_second ;

  for(uint t = 0; t < nsteps; t += deltat )  
    {
      cout << t << endl;
      count += 1;
      for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
        {
          count_atoms = 0;
          for(uint i = 0;i < natoms;++i)
            {
              idi = natoms*t+i;          
              for(uint j = 0;j < natoms;++j)
                {
                  idj = natoms*t+j;
                  if (r[idi].symbol[0] == alpha && r[idj].symbol[0] == beta && i != j)
                    {                                       
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L);       
                      if (rij <= rad[RDF_i] + (RDF_h/2.0) && rij > rad[RDF_i] - (RDF_h/2.0))
                        {
                          count_atoms += 1;  
                        }
                    }
                }

            }
          RDF[RDF_i]  += count_atoms  / ( ( no_of_residue_second / V) * 4.0 * PI * rad[RDF_i] * rad[RDF_i] *  RDF_h);
        }
    }

  ofstream outfile(filename);
  for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
    {
      outfile << rad[RDF_i]  << "  " << RDF[RDF_i] / (count * no_of_residue_first )  << endl;
    }
  outfile.close();
  outfile.clear();
}




/************ main implementation ************/
int main ( int argc, char** argv )
{
  uint natoms = 0, nsteps = 0, id = 0, nmol = 0;
  string temp, xyzfilename, psffilename, fieldfilename;
  float dt = 0;
  vector<float> L (3,0.0);

  //input arguments
  xyzfilename = argv[1];
  nsteps = atoi(argv[2]);
  dt = atof(argv[3]);
  L[0] = atof(argv[4]);
  L[1] = atof(argv[5]);
  L[2] = atof(argv[6]);
  natoms = atoi(argv[7]);
  vector<Atom> r(natoms*nsteps);

  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  BringintoBox(r, nsteps, natoms, L);
  Printrdf(r, nsteps, natoms, L, dt, argv[8]);
  return 0;
}



void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L)
{
  string temp;
  ifstream xyzfile(xyzfilename);
  for(uint t = 0; t < nsteps; ++t )
    { 
      xyzfile >> natoms; 
      getline(xyzfile, temp); 
      getline(xyzfile, temp); 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          xyzfile >> r[id].symbol >> r[id].x  >> r[id].y  >> r[id].z;           
          //cout << r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << endl;
        }
    }
  xyzfile.close();
  xyzfile.clear();
}


/* Broken bonds due to PBC are solved here, Minimum image convention was applied to bring them together */
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L)
{
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          if(r[id].x > L[0] || r[id].x < 0 ) { r[id].x  = min_distance(r[id].x, L[0]); }
          if(r[id].y > L[1] || r[id].y < 0 ) { r[id].y  = min_distance(r[id].y, L[1]); }
          if(r[id].z > L[2] || r[id].z < 0 ) { r[id].z  = min_distance(r[id].z, L[2]); }

          if(r[id].x > L[0] ) { r[id].x  -= L[0] ; }
          if(r[id].y > L[1] ) { r[id].y  -= L[1] ; }
          if(r[id].z > L[2] ) { r[id].z  -= L[2] ; }

          if(r[id].x < 0 ) { r[id].x += L[0] ; }
          if(r[id].y < 0 ) { r[id].y += L[1] ; }
          if(r[id].z < 0 ) { r[id].z += L[2] ; }
        }
    }

  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;
 
          // Unites SO4
          if(r[id].symbol[0] == 'S')
            {
              for(uint j = 0;j < natoms;++j)
                {
                  uint id2 = natoms*t+j; 
                  if(i == j){}
                  else if(r[id2].symbol[0] == 'O')
                    {
                      float rij = mindis(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z, L);
                      if(rij < 2.25)
                        { 
                          float rij2 = norm(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z);
                          if(abs(rij2 -  rij) > 0.1)
                            {  
                             float rij3 = norm(r[id].x - r[id2].x, 0, 0);
                             if(abs(rij3) > 2.5)
                               {   
                                 if(r[id].x > r[id2].x ) { r[id2].x += L[0]; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L[0]; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y += L[1]; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L[1]; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L[2]; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L[2]; }                                                                                                        
                               }                                      
                            }
                        }
                    }
                } 
            }


          // Unites H2O
          if(r[id].symbol[0] == 'O')
            {
              for(uint j = 0;j < natoms;++j)
                {
                  uint id2 = natoms*t+j; 
                  if(i == j){}
                  else if(r[id2].symbol[0] == 'H')
                    {
                      float rij = mindis(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z, L);
                      if(rij < 1.25)
                        { 
                          float rij2 = norm(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z);
                          if(abs(rij2 -  rij) > 0.1)
                            {  
                             float rij3 = norm(r[id].x - r[id2].x, 0, 0);
                             if(abs(rij3) > 2.5)
                               {   
                                 if(r[id].x > r[id2].x ) { r[id2].x += L[0]; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L[0]; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y  =  r[id2].y += L[1]; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L[1]; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L[2]; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L[2]; }                                                                                                        
                               }                                      
                            }
                        }
                    }
                } 
            }
        }
    }
}



double min_distance(double r, float l){ return r - l * round(r/l);}
double mindis(double dx,double dy,double dz, const vector<float> & L)
{ 
  return norm(dx - L[0] * round(dx/L[0]),dy - L[1] * round(dy/L[1]), dz - L[2] * round(dz/L[2]));
}








