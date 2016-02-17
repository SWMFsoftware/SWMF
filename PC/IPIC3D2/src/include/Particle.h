#ifndef _Particle_
#define _Particle_
#include "ipicdefs.h" // for longid

// Depends on width of vector unit;
// need to be known at compile time.
//
#define AoS_PCLS_AT_A_TIME 2

namespace ParticleType
{
  enum Type
  {
    AoS = 0,
    SoA,
    synched
  };
}

template <class T>
class Larray;
// intended to occupy 64 bytes
//
// particle for a specific species
class SpeciesParticle
{
  double u[3];
  double q;
  double x[3];
  double t;
 public:
  SpeciesParticle(){}
  SpeciesParticle(
    double u_,
    double v_,
    double w_,
    double q_,
    double x_,
    double y_,
    double z_,
    double t_)
  {
    u[0]=u_;
    u[1]=v_;
    u[2]=w_;
    q=q_;
    x[0]=x_;
    x[1]=y_;
    x[2]=z_;
    t=t_;
  }
  // accessors
  // double component(int i){ return u[i]; } // a hack
  double get_u(int i)const{ return u[i]; }
  double get_q()const{ return q; }
  double get_x(int i)const{ return x[i]; }
  double get_t()const{ return t; }
  void set_u(double* in, int n=3) { for(int i=0;i<n;i++) u[i] = in[i]; }
  void set_u(int i, double in) { u[i] = in; }
  void set_q(double in) { q = in; }
  void set_x(int i, double in) { x[i] = in; }
  void set_t(double in){ t=in; }
  // tracking particles would actually use q for the ID
  longid get_ID()const{ return longid(t); }
  void set_ID(longid in){ t = double(in); }
  // alternative accessors
  double get_x()const{ return x[0]; }
  double get_y()const{ return x[1]; }
  double get_z()const{ return x[2]; }
  double get_u()const{ return u[0]; }
  double get_v()const{ return u[1]; }
  double get_w()const{ return u[2]; }
  double& fetch_x(){ return x[0]; }
  double& fetch_y(){ return x[1]; }
  double& fetch_z(){ return x[2]; }
  double& fetch_q(){ return q; }
  double& fetch_u(){ return u[0]; }
  double& fetch_v(){ return u[1]; }
  double& fetch_w(){ return u[2]; }
  double& fetch_t(){ return t; }
  void set_x(double in){ x[0]=in; }
  void set_y(double in){ x[1]=in; }
  void set_z(double in){ x[2]=in; }
  void set_u(double in){ u[0]=in; }
  void set_v(double in){ u[1]=in; }
  void set_w(double in){ u[2]=in; }
  void set_to_zero()
  {
    for(int i=0;i<8;i++) u[i]=0;
  }
  void set(
    double _u, double _v, double _w, double _q,
    double _x, double _y, double _z, double _t
    )
  {
    u[0] = _u; u[1] = _v; u[2] = _w; q = _q;
    x[0] = _x; x[1] = _y; x[2] = _z; t = _t;
  }
};

// to support SoA notation
//
// this class will simply be defined differently
// when underlying representation is SoA
//
//class FetchPclComponent
//{
//  int offset;
//  Larray<SpeciesParticle>& list;
// public:
//  FetchPclComponent( Larray<SpeciesParticle>& _list, int _offset)
//  : list(_list), offset(_offset)
//  { }
//  double operator[](int i)
//  {
//    return list[i].component(offset);
//    // return component(offset)[i];
//  }
//};

// intended to occupy 64 bytes
//
// dust particle for second-order-accuracy implicit advance
#if 0
class IDpcl
{
  int c[3]; // cell
  float q; // charge
  float x[3]; // position
  float t; // subcycle time
  float hdx[3]; // xavg = x + hdx
  float qom; // charge to mass ratio of particle
  float u[3];
  float m; // mass of particle
};
#endif

#endif
