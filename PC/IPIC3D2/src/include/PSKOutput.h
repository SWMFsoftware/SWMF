/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*******************************************************************************************
  PSKOutput.h  -  Framework classes for PARSEK output
  -------------------
developers: D. Burgess, June/July 2006
 ********************************************************************************************/
#ifndef _PSK_OUTPUT_H_
#define _PSK_OUTPUT_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

#include "errors.h"
#include "Grid.h"
#include "PSKException.h"
#include "Particles3Dcomm.h"
#include "Field.h"
#include "Collective.h"
#include "VCtopology3D.h"
#include "MPIdata.h"
#include "ipicdefs.h"

using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

/**
 * 
 * Framework classes for PARSEK output
 * @date June/July 2006
 * @author David Burgess
 * @version 2.0
 *
 */

namespace PSK {

  /** class for handling IO exception*/
  class OutputException:public Exception {
  public:
  OutputException(const std::string & err_str, const std::string fn_str = "", int sys_errno = 0):Exception(err_str, fn_str, sys_errno) {
      _type_str += "::OutputException";
  } OutputException(const char *err_str, const char *fn_str = "", int sys_errno = 0):Exception(err_str, fn_str, sys_errno) {
      _type_str += "::OutputException";
  }};

  /** 
   *
   * !\brief Class for dimensionality of rectangular arrays
   *
   */
  class Dimens {

    std::vector < int >_dimens;
    friend class OutputAdaptor;

  public:
    Dimens(void) {;
    } Dimens(const Dimens & dimens):_dimens(dimens._dimens) {;
    }
    Dimens(const std::vector < int >&dimens):_dimens(dimens) {;
    }
    Dimens(int d1):_dimens(1) {
      _dimens[0] = d1;
    }
    Dimens(int d1, int d2):_dimens(2) {
      _dimens[0] = d1;
      _dimens[1] = d2;
    }
    Dimens(int d1, int d2, int d3):_dimens(3) {
      _dimens[0] = d1;
      _dimens[1] = d2;
      _dimens[2] = d3;
    }
    Dimens(int d1, int d2, int d3, int d4):_dimens(4) {
      _dimens[0] = d1;
      _dimens[1] = d2;
      _dimens[2] = d3;
      _dimens[3] = d4;
    }

    int size(void) const {
      return _dimens.size();
    } int operator[] (const int i) const {
      return _dimens[i];
    } int nels(void) const {
      int n = 1;
      for (int i = 0; i < _dimens.size(); ++i)
        n *= _dimens[i];
      return n;
    }
  };

  /**
   *  !\brief Virtual base class  
   *
   */
  class OutputAdaptor {

  public:
    OutputAdaptor(void) {;
    } virtual void open(const std::string & outf) {
      eprintf("Function not implemented");
      eprintf("Function not implemented");
    }
    virtual void close(void) {
      eprintf("Function not implemented");
    }

    // write int functions
    virtual void write(const std::string & objname, int i) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const int *i_array) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const longid *i_array) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const std::vector < int >&i_array) {
      eprintf("Function not implemented");
    }

    // write float functions
    virtual void write(const std::string & objname, float f) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const float *f_array) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const std::vector < float >&f_array) {
      eprintf("Function not implemented");
    }

    // write double functions
    virtual void write(const std::string & objname, double d) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const double *d_array) {
      eprintf("Function not implemented");
    }
    virtual void write(const std::string & objname, const Dimens dimens, const std::vector < double >&d_array) {
      eprintf("Function not implemented");
    }

  };
  /** bse class for output agent */
  class OutputAgentBase {

  public:
    OutputAgentBase(void) {;
    } OutputAgentBase(const OutputAgentBase & a) {;
    }

    virtual void output(const std::string & tag, int cycle) = 0;
    virtual void output(const std::string & tag, int cycle, int sample) = 0;

    virtual void open(const std::string & outf) = 0;
    virtual void close(void) = 0;


  };

  /** \brief Base class for OutputAgents using template for output adaptor */
template < class Toa > class OutputAgent:public OutputAgentBase {

  protected:
    Toa output_adaptor;

  public:
    OutputAgent(void) {;
    }
    OutputAgent(const OutputAgent & a) {;
    }

    virtual void output(const std::string & tag, int cycle) = 0;
    virtual void output(const std::string & tag, int cycle, int sample) = 0;

    void open(const std::string & outf) {
      output_adaptor.open(outf);
    }
    void open_append(const std::string & outf, bool doEraseFile=false) {
      output_adaptor.open_append(outf,doEraseFile);
    }

    void close(void) {
      output_adaptor.close();
    }
    // grid
    Grid *mygrid;
  };

  /** \brief Container (list) class for OutputAgents */
  template < class Toa > class OutputManager {
    std::list < OutputAgentBase * >agents_list;

  public:
    OutputManager(void) {;
    }

    void push_back(OutputAgentBase * a_p) {
      agents_list.push_back(a_p);
    }

    void output(const std::string & tag, int cycle) {
      typename std::list < OutputAgentBase * >::iterator p = agents_list.begin();
      while (p != agents_list.end())
        (*p++)->output(tag, cycle);
    }

    void output(const std::string & tag, int cycle, int sample) {
      typename std::list < OutputAgentBase * >::iterator p = agents_list.begin();
      while (p != agents_list.end())
        (*p++)->output(tag, cycle, sample);
    }

  };


  // ======================================================================


  class coutOutputAdaptor:public OutputAdaptor {

  public:
    coutOutputAdaptor(void) {;
    } void open(const std::string & outf) {
      std::cout << "coutPSKOutputAdaptor open() file: " << outf << "\n";
    }
    void close(void) {
      std::cout << "coutPSKOutputAdaptor close()\n";
    }


    // write int functions
    void write(const std::string & objname, int i) {
      std::cout << "coutPSKOutputAdaptor write int: <" << objname << "> : " << i << "\n";
    }

    void write(const std::string & objname, const Dimens dimens, const int *i_array) {
      std::cout << "coutPSKOutputAdaptor write int* array: <" << objname << "> : " << "\n";
    }
    void write(const std::string & objname, const Dimens dimens, const std::vector < int >&i_array) {
      std::cout << "coutPSKOutputAdaptor write vector<int> array: <" << objname << "> : " << "\n";
    }


    // write float functions
    void write(const std::string & objname, float f) {
      std::cout << "coutPSKOutputAdaptor write float: <" << objname << "> : " << f << "\n";
    }
    void write(const std::string & objname, const Dimens dimens, const float *f_array) {
      std::cout << "coutPSKOutputAdaptor write float* array: <" << objname << "> : " << "\n";
    }
    void write(const std::string & objname, const Dimens dimens, const std::vector < float >&f_array) {
      std::cout << "coutPSKOutputAdaptor write vector<float> array: <" << objname << "> : " << "\n";
    }

    // write double functions
    void write(const std::string & objname, double d) {
      std::cout << "coutPSKOutputAdaptor write double: <" << objname << "> : " << d << "\n";
    }
    void write(const std::string & objname, const Dimens dimens, const double *d_array) {
      std::cout << "coutPSKOutputAdaptor write double* array: <" << objname << "> : " << "\n";
    }
    void write(const std::string & objname, const Dimens dimens, const std::vector < double >&d_array) {
      std::cout << "coutPSKOutputAdaptor write vector<double> array: <" << objname << "> : " << "\n";
    }

  };

}                               // end namespace PSK


template < class Toa > class myOutputAgent:public PSK::OutputAgent < Toa > {
  Field *_field;
  Grid *_grid;
  VCtopology3D *_vct;
  Collective *_col;
  int ns;
  std::vector < Particles * >_part;

public:
  myOutputAgent(void) {;
  }

  void set_simulation_pointers(Field * field, Grid * grid, VCtopology3D * vct, Collective * col) {
    _field = field;
    _grid = grid;
    _vct = vct;
    _col = col;
  }

  void set_simulation_pointers_part(Particles * part) {
    _part.push_back(part);
  }

//  void set_simulation_pointers_testpart(Particles * part) {
//    _testpart.push_back(part);
//  }

  /** method to write on disk. Acceptable tags are:

    collective
    total_topology 
    proc_topology
    Ball --> to write all B components
    Bx,By,Bz
    Eall --> to write all E components
    Ex,Ey,Ez
    phi --> scalar vector
    Jall --> to write all J (current density) components
    Jx,Jy,Jz
    Jsall --> to write all Js (current densities for each species) components
    Jxs,Jys,Jzs
    rho -> net charge density
    rhos -> charge densities for each species
    pressure -> pressure tensor for each species
    position -> particle position (x,y)
    velocity -> particle velocity (u,v,w)
    q -> particle charge
    ID -> particle ID (note: TrackParticleID has to be set true in Collective)
    k_energy -> kinetic energy for each species
    B_energy -> energy of magnetic field
    E_energy -> energy of electric field
*/

  void output(const string & tag, int cycle) {
    stringstream ss;
    stringstream cc;
    stringstream ii;
    ss << MPIdata::instance().get_rank();
    cc << cycle;
    const int ns = _col->getNs();
    if (tag.find("last_cycle", 0) != string::npos)
      this->output_adaptor.write("/last_cycle", cycle);
    if (tag.find("collective", 0) != string::npos) {

      this->output_adaptor.write("/collective/Lx", _col->getLx());
      this->output_adaptor.write("/collective/Ly", _col->getLy());
      this->output_adaptor.write("/collective/Lz", _col->getLz());
      this->output_adaptor.write("/collective/x_center", _col->getx_center());
      this->output_adaptor.write("/collective/y_center", _col->gety_center());
      this->output_adaptor.write("/collective/z_center", _col->getz_center());
      this->output_adaptor.write("/collective/L_square", _col->getL_square());
      this->output_adaptor.write("/collective/Bx0", _col->getB0x());
      this->output_adaptor.write("/collective/By0", _col->getB0y());
      this->output_adaptor.write("/collective/Bz0", _col->getB0z());
      this->output_adaptor.write("/collective/Nxc", _col->getNxc());
      this->output_adaptor.write("/collective/Nyc", _col->getNyc());
      this->output_adaptor.write("/collective/Nzc", _col->getNzc());
      this->output_adaptor.write("/collective/Dx", _col->getDx());
      this->output_adaptor.write("/collective/Dy", _col->getDy());
      this->output_adaptor.write("/collective/Dz", _col->getDz());
      this->output_adaptor.write("/collective/Dt", _col->getDt());
      this->output_adaptor.write("/collective/Th", _col->getTh());
      this->output_adaptor.write("/collective/Ncycles", _col->getNcycles());
      this->output_adaptor.write("/collective/Ns", _col->getNs());
      this->output_adaptor.write("/collective/NsTestPart", _col->getNsTestPart());
      this->output_adaptor.write("/collective/c", _col->getC());
      this->output_adaptor.write("/collective/Smooth", _col->getSmooth());

      this->output_adaptor.write("/collective/bc/PfaceXright", _col->getBcPfaceXright());
      this->output_adaptor.write("/collective/bc/PfaceXleft", _col->getBcPfaceXleft());
      this->output_adaptor.write("/collective/bc/PfaceYright", _col->getBcPfaceYright());
      this->output_adaptor.write("/collective/bc/PfaceYleft", _col->getBcPfaceYleft());
      this->output_adaptor.write("/collective/bc/PfaceZright", _col->getBcPfaceZright());
      this->output_adaptor.write("/collective/bc/PfaceZleft", _col->getBcPfaceZleft());

      this->output_adaptor.write("/collective/bc/PHIfaceXright", _col->getBcPHIfaceXright());
      this->output_adaptor.write("/collective/bc/PHIfaceXleft", _col->getBcPHIfaceXleft());
      this->output_adaptor.write("/collective/bc/PHIfaceYright", _col->getBcPHIfaceYright());
      this->output_adaptor.write("/collective/bc/PHIfaceYleft", _col->getBcPHIfaceYleft());
      this->output_adaptor.write("/collective/bc/PHIfaceZright", _col->getBcPHIfaceZright());
      this->output_adaptor.write("/collective/bc/PHIfaceZleft", _col->getBcPHIfaceZleft());


      this->output_adaptor.write("/collective/bc/EMfaceXright", _col->getBcEMfaceXright());
      this->output_adaptor.write("/collective/bc/EMfaceXleft", _col->getBcEMfaceXleft());
      this->output_adaptor.write("/collective/bc/EMfaceYright", _col->getBcEMfaceYright());
      this->output_adaptor.write("/collective/bc/EMfaceYleft", _col->getBcEMfaceYleft());
      this->output_adaptor.write("/collective/bc/EMfaceZright", _col->getBcEMfaceZright());
      this->output_adaptor.write("/collective/bc/EMfaceZleft", _col->getBcEMfaceZleft());


      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/collective/species_" + ii.str() + "/Npcelx", _col->getNpcelx(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/Npcely", _col->getNpcely(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/Npcelz", _col->getNpcelz(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/qom", _col->getQOM(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/uth", _col->getUth(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/vth", _col->getVth(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/wth", _col->getWth(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/u0", _col->getU0(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/v0", _col->getV0(i));
        this->output_adaptor.write("/collective/species_" + ii.str() + "/w0", _col->getW0(i));
      };

      const int nstestpart = _col->getNsTestPart();
      for (int i = 0; i < nstestpart; ++i) {
        stringstream ii;
        ii << (i+ns);
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/Npcelx", _col->getNpcelx(i+ns));
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/Npcely", _col->getNpcely(i+ns));
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/Npcelz", _col->getNpcelz(i+ns));
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/qom", 	  _col->getQOM(i+ns));
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/pitch_angle", _col->getPitchAngle(i));
        this->output_adaptor.write("/collective/testspecies_" + ii.str() + "/energy", 	   _col->getEnergy(i));
      };

    }

    if (tag.find("total_topology", 0) != string::npos) {
      this->output_adaptor.write("/topology/XLEN", _vct->getXLEN());
      this->output_adaptor.write("/topology/YLEN", _vct->getYLEN());
      this->output_adaptor.write("/topology/ZLEN", _vct->getZLEN());
      this->output_adaptor.write("/topology/Nprocs", _vct->getNprocs());
      this->output_adaptor.write("/topology/periodicX", _vct->getPERIODICX());
      this->output_adaptor.write("/topology/periodicY", _vct->getPERIODICY());
      this->output_adaptor.write("/topology/periodicZ", _vct->getPERIODICZ());

    }

    if (tag.find("proc_topology", 0) != string::npos) {
      int *coord = new int[3];
      coord[0] = _vct->getCoordinates(0);
      coord[1] = _vct->getCoordinates(1);
      coord[2] = _vct->getCoordinates(2);
      this->output_adaptor.write("/topology/cartesian_coord", PSK::Dimens(3), coord);
      this->output_adaptor.write("/topology/cartesian_rank", _vct->getCartesian_rank());
      this->output_adaptor.write("/topology/Xleft_neighbor", _vct->getXleft_neighbor());
      this->output_adaptor.write("/topology/Xright_neighbor", _vct->getXright_neighbor());
      this->output_adaptor.write("/topology/Yleft_neighbor", _vct->getYleft_neighbor());
      this->output_adaptor.write("/topology/Yright_neighbor", _vct->getYright_neighbor());
      this->output_adaptor.write("/topology/Zleft_neighbor", _vct->getZleft_neighbor());
      this->output_adaptor.write("/topology/Zright_neighbor", _vct->getZright_neighbor());
      delete[]coord;
    }

    if (tag.find("Bcall", 0) != string::npos) {
      this->output_adaptor.write("/fields/Bxc/cycle_" + cc.str(), PSK::Dimens(_grid->getNXC() - 2, _grid->getNYC() - 2, _grid->getNZC() - 2), _field->getBxc());
      this->output_adaptor.write("/fields/Byc/cycle_" + cc.str(), PSK::Dimens(_grid->getNXC() - 2, _grid->getNYC() - 2, _grid->getNZC() - 2), _field->getByc());
      this->output_adaptor.write("/fields/Bzc/cycle_" + cc.str(), PSK::Dimens(_grid->getNXC() - 2, _grid->getNYC() - 2, _grid->getNZC() - 2), _field->getBzc());
    }
    
    // Bfield is written without ghost cells and defined in nodes
    if (tag.find("Ball", 0) != string::npos) {
      this->output_adaptor.write("/fields/Bx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getBxTot());//_field->getBx()
      this->output_adaptor.write("/fields/By/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getByTot());//_field->getBy()
      this->output_adaptor.write("/fields/Bz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getBzTot());//_field->getBz()
    }
    else if (tag.find("Bx", 0) != string::npos) {
      this->output_adaptor.write("/fields/Bx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getBxTot());//_field->getBx()

    }
    else if (tag.find("By", 0) != string::npos) {
      this->output_adaptor.write("/fields/By/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getByTot());//_field->getBy()

    }
    else if (tag.find("Bz", 0) != string::npos) {
      this->output_adaptor.write("/fields/Bz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getBzTot());//_field->getBz()

    }


    // Efield is written without ghost cells and defined in nodes

    if (tag.find("Eall", 0) != string::npos) {
      this->output_adaptor.write("/fields/Ex/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEx());

      this->output_adaptor.write("/fields/Ey/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEy());
      this->output_adaptor.write("/fields/Ez/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEz());
    }
    else if (tag.find("Ex", 0) != string::npos) {
      this->output_adaptor.write("/fields/Ex/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEx());
    }

    else if (tag.find("Ey", 0) != string::npos) {
      this->output_adaptor.write("/fields/Ey/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEy());
    }
    else if (tag.find("Ez", 0) != string::npos) {
      this->output_adaptor.write("/fields/Ez/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getEz());
    }

    // PHI is written without ghost cells and defined in centers 
    if (tag.find("phi", 0) != string::npos) {
      this->output_adaptor.write("/potentials/phi/cycle_" + cc.str(), PSK::Dimens(_grid->getNXC() - 2, _grid->getNYC() - 2, _grid->getNZC() - 2), _field->getPHI());
    }

    // J (current density) is written without ghost cells and defined in nodes
    if (tag.find("Jall", 0) != string::npos) {
      _field->sumOverSpeciesJ();
      this->output_adaptor.write("/moments/Jx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJx());
      this->output_adaptor.write("/moments/Jy/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJy());
      this->output_adaptor.write("/moments/Jz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJz());
    }
    else if (tag.find("Jx", 0) != string::npos) {
      _field->sumOverSpeciesJ();
      this->output_adaptor.write("/moments/Jx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJx());
    }
    else if (tag.find("Jy", 0) != string::npos) {
      _field->sumOverSpeciesJ();
      this->output_adaptor.write("/moments/Jy/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJy());
    }
    else if (tag.find("Jz", 0) != string::npos) {
      _field->sumOverSpeciesJ();
      this->output_adaptor.write("/moments/Jz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getJz());
    }

    // Js (current density for species s) is written without ghost cells and defined in nodes
    if (tag.find("Jsall", 0) != string::npos) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJxs());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jy/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJys());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJzs());
      }
    }
    else if (tag.find("Jxs", 0) != string::npos) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jx/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJxs());
      }
    }
    else if (tag.find("Jys", 0) != string::npos) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jy/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJys());
      }
    }
    else if (tag.find("Jzs", 0) != string::npos) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/Jz/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getJzs());
      }
    }

    // rhos (number density for species s) is written without ghost cells and defined in nodes
    if (tag.find("rhos", 0) != string::npos) {

      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/rho/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getRHOns());
      }
    }

    // rho (number density ) is written without ghost cells and defined in nodes
//    if (tag.find("rho", 0) != string::npos) {
//      this->output_adaptor.write("/moments/rho/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), _field->getRHOn());
//    }
    // pressure for species s is written without ghost cells and defined in nodes
    if (tag.find("pressure", 0) != string::npos) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pXX/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpXXsn());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pXY/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpXYsn());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pXZ/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpXZsn());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pYY/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpYYsn());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pYZ/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpYZsn());
        this->output_adaptor.write("/moments/species_" + ii.str() + "/pZZ/cycle_" + cc.str(), PSK::Dimens(_grid->getNXN() - 2, _grid->getNYN() - 2, _grid->getNZN() - 2), i, _field->getpZZsn());
      }
    }

    // kinetic energy for species s (normalized on q)
    if (tag.find("k_energy", 0) != string::npos) {
      double K_en;
      for (int i = 0; i < ns; ++i) {
        K_en = 0;
        stringstream ii;
        ii << i;
        double qom = fabs(_col->getQOM(i));
        for (int n = 0; n < _part[i]->getNOP(); n++) {
          K_en += fabs(_part[i]->getQ(n)) * _part[i]->getU(n) * _part[i]->getU(n);
          K_en += fabs(_part[i]->getQ(n)) * _part[i]->getV(n) * _part[i]->getV(n);
          K_en += fabs(_part[i]->getQ(n)) * _part[i]->getW(n) * _part[i]->getW(n);
        }
        K_en *= 0.5 / qom;
        this->output_adaptor.write("/energy/kinetic/species_" + ii.str() + "/cycle_" + cc.str(), K_en);
      }
    }

    // magnetic energy (taken on nodes)
    if (tag.find("B_energy", 0) != string::npos) {
      double B_en;
      B_en = 0;
      for (int i = 1; i < _grid->getNXN(); i++) // get rid of ghost cells
        for (int j = 1; j < _grid->getNYN(); j++)
          for (int k = 1; k < _grid->getNYN(); k++)
            //B_en += _field->getBx(i, j, k) * _field->getBx(i, j, k) + _field->getBy(i, j, k) * _field->getBy(i, j, k) + _field->getBz(i, j, k) * _field->getBz(i, j, k);
        	  B_en += _field->getBxTot(i, j, k) * _field->getBxTot(i, j, k) + _field->getByTot(i, j, k) * _field->getByTot(i, j, k) + _field->getBzTot(i, j, k) * _field->getBzTot(i, j, k);
      B_en = B_en / 32 / atan(1.0); // here there should be a getfourPI for the right value of pi 
      this->output_adaptor.write("/energy/magnetic/cycle_" + cc.str(), B_en);
    }


    // electric energy (taken on nodes)
    if (tag.find("E_energy", 0) != string::npos) {
      double E_en;
      E_en = 0;
      for (int i = 1; i < _grid->getNXN(); i++) // get rid of ghost cells
        for (int j = 1; j < _grid->getNYN(); j++)
          for (int k = 1; k < _grid->getNYN(); k++)
            E_en += _field->getEx(i, j, k) * _field->getEx(i, j, k) + _field->getEy(i, j, k) * _field->getEy(i, j, k) + _field->getEz(i, j, k) * _field->getEz(i, j, k);
      E_en = E_en / 32 / atan(1.0); // here there should be a getfourPI for the right value of pi 
      this->output_adaptor.write("/energy/electric/cycle_" + cc.str(), E_en);
    }

  }

  void output(const string & tag, int cycle, int sample) {
    stringstream cc;
    cc << cycle;
    const int ns = _col->getNs();
    const int nstestpart = _col->getNsTestPart();

    /* Store the seed for the random number generator to be able to 
       reproduce the same distribution on ghost particles*/
    if (tag.find("pseudo_random_seed",0) != string::npos){
      for ( int i=0; i<ns; ++i ){
        stringstream ii;
        ii << i;
        this->output_adaptor.write ("/particles/species_"+ii.str()+"/pseudo_random_seed",_part[i]->getPseudoRandomSeed());
      }    
    }    
    
    // Particle position
    if (tag.find("position", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/x/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getXall());
        this->output_adaptor.write("/particles/species_" + ii.str() + "/y/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getYall());
        this->output_adaptor.write("/particles/species_" + ii.str() + "/z/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getZall());
      }
    }
    // Test Particle position
    else if (tag.find("testpartpos", 0) != string::npos & sample == 0) {
        for (int i = 0; i < nstestpart; ++i) {
            stringstream ii;
            ii << (_part[i+ns]->get_species_num());
            this->output_adaptor.write("/testparticles/species_" + ii.str() + "/x/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getXall());
            this->output_adaptor.write("/testparticles/species_" + ii.str() + "/y/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getYall());
            this->output_adaptor.write("/testparticles/species_" + ii.str() + "/z/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getZall());
        }
   }
    else if (tag.find("x", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/x/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getXall());
      }
    }
    else if (tag.find("x", 0) != string::npos & sample != 0) {
      std::vector < double >X;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const int num_samples = _part[i]->getNOP()/sample;
        X.reserve(num_samples);
        for (int n = 0; n < _part[i]->getNOP(); n += sample) {
          X.push_back(_part[i]->getX(n));
        }
        this->output_adaptor.write("/particles/species_" + ii.str() + "/x/cycle_" + cc.str(), PSK::Dimens(X.size()), X);
      }
    }
    else if (tag.find("position", 0) != string::npos & sample != 0) {
      std::vector < double >X, Y, Z;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const int num_samples = _part[i]->getNOP()/sample;
        X.reserve(num_samples);
        Y.reserve(num_samples);
        Z.reserve(num_samples);
        for (int n = 0; n < _part[i]->getNOP(); n += sample) {
          X.push_back(_part[i]->getX(n));
          Y.push_back(_part[i]->getY(n));
          Z.push_back(_part[i]->getZ(n));
        }
        this->output_adaptor.write("/particles/species_" + ii.str() + "/x/cycle_" + cc.str(), PSK::Dimens(X.size()), X);
        this->output_adaptor.write("/particles/species_" + ii.str() + "/y/cycle_" + cc.str(), PSK::Dimens(Y.size()), Y);
        this->output_adaptor.write("/particles/species_" + ii.str() + "/z/cycle_" + cc.str(), PSK::Dimens(Z.size()), Z);
      }
    }


    // Particle velocity
    if (tag.find("velocity", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/u/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getUall());
        this->output_adaptor.write("/particles/species_" + ii.str() + "/v/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getVall());
        this->output_adaptor.write("/particles/species_" + ii.str() + "/w/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getWall());
      }
    }else if (tag.find("testpartvel", 0) != string::npos & sample == 0) {
        for (int i = 0; i < nstestpart; ++i) {
          stringstream ii;
          ii << (_part[i+ns]->get_species_num());
          this->output_adaptor.write("/testparticles/species_" + ii.str() + "/u/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getUall());
          this->output_adaptor.write("/testparticles/species_" + ii.str() + "/v/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getVall());
          this->output_adaptor.write("/testparticles/species_" + ii.str() + "/w/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getWall());
        }
      }
    else if (tag.find("u", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/u/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getUall());
      }
    }
    else if (tag.find("u", 0) != string::npos & sample != 0) {
      std::vector < double >U;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const int num_samples = _part[i]->getNOP()/sample;
        U.reserve(num_samples);
        for (int n = 0; n < _part[i]->getNOP(); n += sample) {
          U.push_back(_part[i]->getU(n));
        }
        this->output_adaptor.write("/particles/species_" + ii.str() + "/u/cycle_" + cc.str(), PSK::Dimens(U.size()), U);
      }
    }
    else if (tag.find("velocity", 0) != string::npos & sample != 0) {
      std::vector < double >U, V, W;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const int num_samples = _part[i]->getNOP()/sample;
        U.reserve(num_samples);
        V.reserve(num_samples);
        W.reserve(num_samples);
        for (int n = 0; n < _part[i]->getNOP(); n += sample) {
          U.push_back(_part[i]->getU(n));
          V.push_back(_part[i]->getV(n));
          W.push_back(_part[i]->getW(n));
        }
        this->output_adaptor.write("/particles/species_" + ii.str() + "/u/cycle_" + cc.str(), PSK::Dimens(U.size()), U);
        this->output_adaptor.write("/particles/species_" + ii.str() + "/v/cycle_" + cc.str(), PSK::Dimens(V.size()), V);
        this->output_adaptor.write("/particles/species_" + ii.str() + "/w/cycle_" + cc.str(), PSK::Dimens(W.size()), W);
      }
    }


    // Particle charge
    if (tag.find("q", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/q/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getQall());
      }
    }
    // Test Particle charge
    else if (tag.find("testpartcharge", 0) != string::npos & sample == 0) {
      for (int i = 0; i < nstestpart; ++i) {
        stringstream ii;
        ii <<  (_part[i+ns]->get_species_num());
        this->output_adaptor.write("/testparticles/species_" + ii.str() + "/q/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getQall());
      }
    }
    else if (tag.find("q", 0) != string::npos & sample != 0) {
      std::vector < double >Q;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const int num_samples = _part[i]->getNOP()/sample;
        Q.reserve(num_samples);
        for (int n = 0; n < _part[i]->getNOP(); n += sample) {
          Q.push_back(_part[i]->getQ(n));
        }
        this->output_adaptor.write("/particles/species_" + ii.str() + "/q/cycle_" + cc.str(), PSK::Dimens(Q.size()), Q);
      }
    }


    // Particle ID
    //
    // (why was this using "long")?
    
    if (tag.find("ID", 0) != string::npos & sample == 0) {
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        this->output_adaptor.write("/particles/species_" + ii.str() + "/ID/cycle_" + cc.str(), PSK::Dimens(_part[i]->getNOP()), _part[i]->getParticleIDall());
      }
    }
    // Test Particle ID
    else if (tag.find("testparttag", 0) != string::npos & sample == 0) {
      for (int i = 0; i < nstestpart; ++i) {
        stringstream ii;
        ii <<  (_part[i+ns]->get_species_num());
        this->output_adaptor.write("/testparticles/species_" + ii.str() + "/ID/cycle_" + cc.str(), PSK::Dimens(_part[i+ns]->getNOP()), _part[i+ns]->getParticleIDall());
      }
    }
    else if (tag.find("ID", 0) != string::npos & sample != 0) {
      std::vector <double>ID;
      for (int i = 0; i < ns; ++i) {
        stringstream ii;
        ii << i;
        const double* pclID = _part[i]->getParticleIDall();
        const int num_samples = _part[i]->getNOP()/sample;
        ID.reserve(num_samples);

        for (int n = 0; n < _part[i]->getNOP(); n += sample)
          ID.push_back(pclID[n]);
        this->output_adaptor.write("/particles/species_" + ii.str() + "/ID/cycle_" + cc.str(), PSK::Dimens(ID.size()), &ID[0]);

      }
    }
  }

#ifdef BATSRUS
  double weightedValue(double ****V, const int ix, const int iy, const int iz, const int is,
		       const double w000, const double w001, const double w010, const double w011, 
		       const double w100, const double w101, const double w110, const double w111){
    double tmp;
    tmp =0.0;
    
    tmp += w000 * V[is][ix][iy][iz];
    tmp += w001 * V[is][ix][iy][iz - 1];
    tmp += w010 * V[is][ix][iy - 1][iz];
    tmp += w011 * V[is][ix][iy - 1][iz - 1];
    tmp += w100 * V[is][ix - 1][iy][iz];
    tmp += w101 * V[is][ix - 1][iy][iz - 1];
    tmp += w110 * V[is][ix - 1][iy - 1][iz];
    tmp += w111 * V[is][ix - 1][iy - 1][iz - 1];
    
    return(tmp);
  }
  
  
  double weightedValue(double ***V, const int ix, const int iy, const int iz, 
		       const double w000, const double w001, const double w010, const double w011, 
		       const double w100, const double w101, const double w110, const double w111){
    double tmp;
    tmp =0.0;
    
    tmp += w000 * V[ix][iy][iz];
    tmp += w001 * V[ix][iy][iz - 1];
    tmp += w010 * V[ix][iy - 1][iz];
    tmp += w011 * V[ix][iy - 1][iz - 1];
    tmp += w100 * V[ix - 1][iy][iz];
    tmp += w101 * V[ix - 1][iy][iz - 1];
    tmp += w110 * V[ix - 1][iy - 1][iz];
    tmp += w111 * V[ix - 1][iy - 1][iz - 1];
     
    return(tmp);
  }

  void getFluidSIatPoint(int nDim, int nPoint, double *Xyz_I, double *data_I, int nVarIn){

    const int nVar             = _col->getnVarFluid();
    const int nIon             = _col->getnIon();
    const bool useAnisoP       = _col->getUseAnisoP();
    const bool useMhdPe        = _col->getUseMhdPe();
    const bool useMultiFluid   = _col->getUseMultiFluid();
    const bool useMultiSpecies = _col->getUseMultiSpecies();
          
    assert_eq(nVar, nVarIn);

    unsigned long n,i,j;	
    double *localState_IV;
    
    double ****Pxx, ****Pyy,****Pzz, ****Pxy, ****Pxz, ****Pyz;
    
    double QoMi, Rhoi, Uxi, Uyi, Uzi, Pi;
    double QoMe, Rhoe, Uxe, Uye, Uze, Pe;
    double Rho, Trace;
    double PeXX, PeYY, PeZZ, PeXY, PeXZ, PeYZ, PiXX, PiYY, PiZZ, PiXY, PiXZ, PiYZ;
    double PtXX, PtYY, PtZZ, PtXY, PtXZ, PtYZ, BX, BY, BZ;
    double PitXX, PitYY, PitZZ, PitXY, PitXZ, PitYZ;
    double Mx, My, Mz; // Total momentum.
    double Mix, Miy, Miz; // Momentum of i species.

    double xmin, ymin, zmin;

    double xp,yp,zp; 
    
    Pxx = _field->getp0XXsn();
    Pyy = _field->getp0YYsn();
    Pzz = _field->getp0ZZsn();
    Pxy = _field->getp0XYsn();
    Pxz = _field->getp0XZsn();
    Pyz = _field->getp0YZsn();

    xmin = _col->getPhyXMin(); 
    ymin = _col->getPhyYMin();
    zmin = _col->getPhyZMin();

    for(int i = 0; i<nPoint*nDim; i++){
      Xyz_I[i] *= _col->getSi2NoL();
    }



    for(int iPoint = 0; iPoint < nPoint; iPoint++){
      if(_col->isThisRun(&Xyz_I[iPoint*nDim])){
	double pos_D[3];
	pos_D[0] = Xyz_I[iPoint*nDim ] - xmin;
	pos_D[1] = Xyz_I[iPoint*nDim+1] - ymin;
	pos_D[2] = Xyz_I[iPoint*nDim+2] - zmin;
	if(! _col->isThisProc(&pos_D[0],true)){
	  cout<<" proc= "<<_vct->getCartesian_rank()
	      <<" xp= "<<Xyz_I[iPoint*nDim]
	      <<" yp= "<<Xyz_I[iPoint*nDim+1]
	      <<" zp= "<<Xyz_I[iPoint*nDim+2]
	      <<" xmin= "<<xmin
	      <<" ymin= "<<ymin
	      <<" zmin= "<<zmin
	      <<endl;
	}
	xp = Xyz_I[iPoint*nDim] - xmin;
	yp = (nDim > 1) ? Xyz_I[iPoint*nDim + 1] - ymin : 0.0;
	zp = (nDim > 2) ? Xyz_I[iPoint*nDim + 2] - zmin : 0.0;

	
	
	double weight_I[8];
	int ix, iy, iz;
	_grid->getInterpolateWeight(xp,yp,zp,ix,iy,iz,weight_I);
	
	const double w000 = weight_I[0];
	const double w001 = weight_I[1];
	const double w010 = weight_I[2];
	const double w011 = weight_I[3];
	const double w100 = weight_I[4];
	const double w101 = weight_I[5];
	const double w110 = weight_I[6];
	const double w111 = weight_I[7];

	n = iPoint*nVar;	

	// Electron
	QoMe = _col->getQOM(0);
	Rhoe = weightedValue(_field->getRHOns(),ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111)/QoMe;
	//cout<<"ix = "<<ix<<" iy = "<<iy<<" iz = "<<iz<<" rhoe = "<<Rhoe<<endl;
	if(fabs(Rhoe)>0){
	  Uxe  = weightedValue(_field->getJxs(),ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMe*Rhoe);
	  Uye  = weightedValue(_field->getJys(),ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMe*Rhoe);
	  Uze  = weightedValue(_field->getJzs(),ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMe*Rhoe);
	}else{
	  Uxe  = 0;
	  Uye  = 0;
	  Uze  = 0;
	}
	PeXX  = weightedValue(Pxx,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	PeYY  = weightedValue(Pyy,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	PeZZ  = weightedValue(Pzz,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	Pe   = (PeXX + PeYY + PeZZ)/3.0;
	PeXY = weightedValue(Pxy,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	PeXZ = weightedValue(Pxz,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	PeYZ = weightedValue(Pyz,ix,iy,iz,0,w000,w001,w010,w011,w100,w101,w110,w111);
	
	BX  = weightedValue(_field->getBx(),ix,iy,iz,w000,w001,w010,w011,w100,w101,w110,w111);
	BY  = weightedValue(_field->getBy(),ix,iy,iz,w000,w001,w010,w011,w100,w101,w110,w111);
	BZ  = weightedValue(_field->getBz(),ix,iy,iz,w000,w001,w010,w011,w100,w101,w110,w111);
	
	for (int i=0; i<nVar; i++){
	  // The variable for hyperbolic clean, is not known by iPIC3D,
	  // but it is needed to be passed back. So data_I need to be
	  // initilized.
	  data_I[n + i] = 0;
	}
	
	data_I[n + _col->iBx] = BX; 
	data_I[n + _col->iBy] = BY;
	data_I[n + _col->iBz] = BZ;
	if(useMhdPe)data_I[n + _col->iPe] = Pe;
	
	Rho  = Rhoe;
	PtXX = PeXX;
	PtYY = PeYY;
	PtZZ = PeZZ;
	PtXY = PeXY;
	PtXZ = PeXZ;
	PtYZ = PeYZ;
	Mx   = Rhoe*Uxe;
	My   = Rhoe*Uye;
	Mz   = Rhoe*Uze;
	
	int iSpecies;
	for(int iIon=0; iIon<nIon; ++iIon){
	  iSpecies = iIon + 1; // iSpecies = 0 is electron.	  
	  QoMi = _col->getQOM(iSpecies);
	  Rhoi = weightedValue(_field->getRHOns(),ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111)/QoMi;
	  if(fabs(Rhoi)>0){
	    Uxi  = weightedValue(_field->getJxs(),ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMi*Rhoi);
	    Uyi  = weightedValue(_field->getJys(),ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMi*Rhoi);
	    Uzi  = weightedValue(_field->getJzs(),ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111)/(QoMi*Rhoi);
	  }else{
	    Uxi  = 0;
	    Uyi  = 0;
	    Uzi  = 0;
	  }
	  PiXX  = weightedValue(Pxx,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);
	  PiYY  = weightedValue(Pyy,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);
	  PiZZ  = weightedValue(Pzz,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);
	  Pi   = (PiXX + PiYY + PiZZ)/3.0;
	  PiXY = weightedValue(Pxy,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);
	  PiXZ = weightedValue(Pxz,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);
	  PiYZ = weightedValue(Pyz,ix,iy,iz,iSpecies,w000,w001,w010,w011,w100,w101,w110,w111);

	  Mix = Rhoi*Uxi;
	  Miy = Rhoi*Uyi;
	  Miz = Rhoi*Uzi;

	  // Sum to total density/pressure.
	  Rho  += Rhoi;
	  PtXX += PiXX;
	  PtYY += PiYY;
	  PtZZ += PiZZ;
	  PtXY += PiXY;
	  PtXZ += PiXZ;
	  PtYZ += PiYZ;
	  Mx   += Mix;
	  My   += Miy;
	  Mz   += Miz;

	  // Density
	  if(useMultiFluid || useMultiSpecies){
	    data_I[n + _col->iRho_I[iIon]] = Rhoi;	  
	  }else{
	    // Only one ion species.Rho = Rhoi + Rhoe
	    data_I[n + _col->iRho_I[iIon]] = Rho;	  
	  }

	  // Pressure.
	  if(useMultiFluid){
	    // ONLY works for iso pressure so far!!!!!
	    data_I[n + _col->iP_I[iIon]] = Pi;
	    if(useAnisoP) {
	      cout<<"Multi-fluid model can not work with aniso pressure now!!"
		  <<endl;
	      abort();
	    }
	  }

	  // Momentum.
	  if(useMultiFluid){
	    data_I[n + _col->iRhoUx_I[iIon]]  = Mix;
	    data_I[n + _col->iRhoUy_I[iIon]]  = Miy;
	    data_I[n + _col->iRhoUz_I[iIon]]  = Miz; 
	  }
	  
	}// iIon

	// Do not includes electron density. The total density passed
	// in MHD side is useless, so it doesnot matter whether Rho include
	// electron or not. -- Yuxi
	if(useMultiSpecies)
	  data_I[n + _col->iRhoTotal] = Rho - Rhoe; 
	
	// Momentum
	if(!useMultiFluid){
	  // Include electron momentum.
	  data_I[n + _col->iRhoUx_I[0]]  = Mx;
	  data_I[n + _col->iRhoUy_I[0]]  = My;
	  data_I[n + _col->iRhoUz_I[0]]  = Mz; 
	}


	//Sum of ion pressure.
	PitXX = PtXX - PeXX;
	PitYY = PtYY - PeYY;
	PitZZ = PtZZ - PeZZ;
	PitXY = PtXY - PeXY;
	PitXZ = PtXZ - PeXZ;
	PitYZ = PtYZ - PeYZ;

	// Pressure
	if(!useMultiFluid){
	  // AnisoP
	  if(useAnisoP){
	    if(useMhdPe)
	      data_I[n + _col->iPpar_I[0]]
		= (BX*PitXX*BX + BY*PitYY*BY + BZ*PitZZ*BZ +
		   2.0*BX*PitXY*BY + 2.0*BX*PitXZ*BZ +
		   2.0*BY*PitYZ*BZ)/(BX*BX+BY*BY+BZ*BZ);
	    else
	      data_I[n + _col->iPpar_I[0]]
		= (BX*PtXX*BX + BY*PtYY*BY + BZ*PtZZ*BZ +
		   2.0*BX*PtXY*BY + 2.0*BX*PtXZ*BZ +
		   2.0*BY*PtYZ*BZ)/(BX*BX+BY*BY+BZ*BZ);	    
	  }// useAnisoP

	  // Isotropic Pressure. 
	  if(useMhdPe) data_I[n + _col->iP_I[0]]  = (PitXX + PitYY + PitZZ)/3;
	  else data_I[n + _col->iP_I[0]]  = (PtXX + PtYY + PtZZ)/3;
	}
	
      }

    }
    
    // Convert to SI units
    for(int iPoint = 0; iPoint < nPoint; iPoint++){
      if(_col->isThisRun(&Xyz_I[iPoint*nDim])){	
	n = iPoint*nVar;
	for (int iVar=0; iVar<nVar; ++iVar){
	  data_I[n+iVar] *= _col->getNo2Si_V(iVar);
	}
      }
    }

    // if it used later (remove if it is useless)
    for(int i = 0; i<nPoint*nDim; i++){
      Xyz_I[i] *= _col->getNo2SiL();
    }
      
  }
#endif
  

};

#endif // _PSK_OUTPUT_H_
