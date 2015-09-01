/*!************************************************************************* EMfields3D.h - ElectroMagnetic fields definition ------------------- begin : May 2008 copyright : KU Leuven developers : Stefano Markidis, Giovanni Lapenta ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H

#include "asserts.h"
#include "ipicfwd.h"
#include "Alloc.h"
#include "Basic.h"


/*! Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver @date May 2008 @par Copyright: (C) 2008 KUL @author Stefano Markidis, Giovanni Lapenta. @version 3.0 */

// dimension of vectors used in fieldForPcls
const int DFIELD_3or4=4; // 4 pads with garbage but is needed for alignment

class Particles3Dcomm;
class Moments10;
class EMfields3D                // :public Field
{
  public:
    /*! constructor */
    EMfields3D(Collective * col, Grid * grid, VirtualTopology3D *vct);
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    void init();
    /*! init beam */
    void initBEAM(double x_center, double y_center, double z_center, double radius);
    /*! initialize GEM challenge */
    void initGEM();
    void initOriginalGEM();
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation();
    /*! initialize GEM challenge with dipole-like tail without perturbation */
    void initGEMDipoleLikeTailNoPert();
    /*! initialize GEM challenge with no Perturbation */
    void initGEMnoPert();
#ifdef BATSRUS
    /*! initialize from BATSRUS */
    void initBATSRUS();
#endif
    /*! Random initial field */
    void initRandomField();
    /*! Init Force Free (JxB=0) */
    void initForceFree();
    /*! initialized with rotated magnetic field */
    void initEM_rotate(double B, double theta);
    /*! add a perturbattion to charge density */
    void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid);
    /*! add a perturbattion to the EM field */
    void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid);
    /*! Initialise a combination of magnetic dipoles */
    void initDipole();
    void initDipole2D();
    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE();
    /*! Image of Poisson Solver (for SOLVER) */
    void PoissonImage(double *image, double *vector);
    /*! Image of Maxwell Solver (for Solver) */
    void MaxwellImage(double *im, double *vector);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov);
    /*! Impose a constant charge inside a spherical zone of the domain */
    void ConstantChargePlanet(double R, double x_center, double y_center, double z_center);
    void ConstantChargePlanet2DPlaneXZ(double R, double x_center, double z_center);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBC();
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBCv2();
    /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
    void calculateB();
    /*! fix B on the boundary for gem challange */
    void fixBgem();
    /*! fix B on the boundary for gem challange */
    void fixBforcefree();

    /*! Calculate the three components of Pi(implicit pressure) cross image vector */
    void PIdot(arr3_double PIdotX, arr3_double PIdotY, arr3_double PIdotZ,
      const_arr3_double vectX, const_arr3_double vectY, const_arr3_double vectZ, int ns);
    /*! Calculate the three components of mu (implicit permeattivity) cross image vector */
    void MUdot(arr3_double MUdotX, arr3_double MUdotY, arr3_double MUdotZ,
      const_arr3_double vectX, const_arr3_double vectY, const_arr3_double vectZ);
    /*! Calculate rho hat, Jx hat, Jy hat, Jz hat */
    void calculateHatFunctions();


    /*! communicate ghost for densities and interp rho from node to center */
    void interpDensitiesN2C();
    /*! set to 0 all the densities fields */
    void setZeroDensities();
    /*! set to 0 primary moments */
    void setZeroPrimaryMoments();
    /*! set to 0 all densities derived from primary moments */
    void setZeroDerivedMoments();
    /*! Sum rhon over species */
    void sumOverSpecies();
    /*! Sum current over different species */
    void sumOverSpeciesJ();
    /*! Smoothing after the interpolation* */
    void smooth(arr3_double vector, int type);
    /*! SPECIES: Smoothing after the interpolation for species fields* */
    void smooth(double value, arr4_double vector, int is, int type);
    /*! smooth the electric field */
    void smoothE();

    /*! copy the field data to the array used to move the particles */
    void set_fieldForPcls();
    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns);
    /*! sum moments (interp_P2G) versions */
    void sumMoments(const Particles3Dcomm* part);
    void sumMoments_AoS(const Particles3Dcomm* part);
    void sumMoments_AoS_intr(const Particles3Dcomm* part);
    void sumMoments_vectorized(const Particles3Dcomm* part);
    void sumMoments_vectorized_AoS(const Particles3Dcomm* part);
    void sumMomentsOld(const Particles3Dcomm& pcls);
    /*! add accumulated moments to the moments for a given species */
    //void addToSpeciesMoments(const TenMoments & in, int is);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction X to current density field at node X,Y,Z */
    void addJx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Y to current density field at node X,Y,Z */
    void addJy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Z to current density field at node X,Y,Z */
    void addJz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! add an amount of pressure density - direction XX to current density field at node X,Y,Z */
    void addPxx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XY to current density field at node X,Y,Z */
    void addPxy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
    void addPxz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YY to current density field at node X,Y,Z */
    void addPyy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
    void addPyz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
    void addPzz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! adjust densities on boundaries that are not periodic */
    void adjustNonPeriodicDensities(int is);


    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
      int dir);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(
      arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX,
      const_arr3_double vectorY,
      const_arr3_double vectorZ,
      int dir);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorRightX(double **susxx, double **susyx, double **suszx);
    void sustensorLeftX (double **susxx, double **susyx, double **suszx);
    void sustensorRightY(double **susxy, double **susyy, double **suszy);
    void sustensorLeftY (double **susxy, double **susyy, double **suszy);
    void sustensorRightZ(double **susxz, double **susyz, double **suszz);
    void sustensorLeftZ (double **susxz, double **susyz, double **suszz);

    /*** accessor methods ***/

    /*! get Potential array */
    arr3_double getPHI() {return PHI;}

    // field components defined on nodes
    //
    double getEx(int X, int Y, int Z) const { return Ex.get(X,Y,Z);}
    double getEy(int X, int Y, int Z) const { return Ey.get(X,Y,Z);}
    double getEz(int X, int Y, int Z) const { return Ez.get(X,Y,Z);}
    double getBx(int X, int Y, int Z) const { return Bxn.get(X,Y,Z);}
    double getBy(int X, int Y, int Z) const { return Byn.get(X,Y,Z);}
    double getBz(int X, int Y, int Z) const { return Bzn.get(X,Y,Z);}
    //
    const_arr4_pfloat get_fieldForPcls() { return fieldForPcls; }
    arr3_double getEx() { return Ex; }
    arr3_double getEy() { return Ey; }
    arr3_double getEz() { return Ez; }
    arr3_double getBx() { return Bxn; }
    arr3_double getBy() { return Byn; }
    arr3_double getBz() { return Bzn; }

    
    //for parallel vtk
    arr3_double getBxc(){return Bxc;};
    arr3_double getByc(){return Byc;};
    arr3_double getBzc(){return Bzc;};


    //arr3_double getRHOc() { return rhoc; }
    //arr3_double getRHOn() { return rhon; }
    //double getRHOc(int X, int Y, int Z) const { return rhoc.get(X,Y,Z);}
    //double getRHOn(int X, int Y, int Z) const { return rhon.get(X,Y,Z);}

    // densities per species:
    //
    double getRHOcs(int X,int Y,int Z,int is)const{return rhocs.get(is,X,Y,Z);}
    double getRHOns(int X,int Y,int Z,int is)const{return rhons.get(is,X,Y,Z);}
    arr4_double getRHOns(){return rhons;}
    arr4_double getRHOcs(){return rhocs;}


    double getBx_ext(int X, int Y, int Z) const{return Bx_ext.get(X,Y,Z);}
    double getBy_ext(int X, int Y, int Z) const{return By_ext.get(X,Y,Z);}
    double getBz_ext(int X, int Y, int Z) const{return Bz_ext.get(X,Y,Z);}
    
    arr3_double getBx_ext() { return Bx_ext; }
    arr3_double getBy_ext() { return By_ext; }
    arr3_double getBz_ext() { return Bz_ext; }


    //B_tot = B + B_ext
    arr3_double getBxTot() { addscale(1.0,Bxn,Bx_ext,Bx_tot,nxn,nyn,nzn); return Bx_tot; }
    arr3_double getByTot() { addscale(1.0,Byn,By_ext,By_tot,nxn,nyn,nzn); return By_tot; }
    arr3_double getBzTot() { addscale(1.0,Bzn,Bz_ext,Bz_tot,nxn,nyn,nzn); return Bz_tot; }
    double getBxTot(int X, int Y, int Z) const{return Bxn.get(X,Y,Z)+Bx_ext.get(X,Y,Z);;}
    double getByTot(int X, int Y, int Z) const{return Byn.get(X,Y,Z)+By_ext.get(X,Y,Z);}
    double getBzTot(int X, int Y, int Z) const{return Bzn.get(X,Y,Z)+Bz_ext.get(X,Y,Z);}

    arr4_double getpXXsn() { return pXXsn; }
    double getpXXsn(int X,int Y,int Z,int is)const{return pXXsn.get(is,X,Y,Z);}

    arr4_double getpXYsn() { return pXYsn; }
    double getpXYsn(int X,int Y,int Z,int is)const{return pXYsn.get(is,X,Y,Z);}

    arr4_double getpXZsn() { return pXZsn; }
    double getpXZsn(int X,int Y,int Z,int is)const{return pXZsn.get(is,X,Y,Z);}

    arr4_double getpYYsn() { return pYYsn; }
    double getpYYsn(int X,int Y,int Z,int is)const{return pYYsn.get(is,X,Y,Z);}

    arr4_double getpYZsn() { return pYZsn; }
    double getpYZsn(int X,int Y,int Z,int is)const{return pYZsn.get(is,X,Y,Z);}

    arr4_double getpZZsn() { return pZZsn; }
    double getpZZsn(int X,int Y,int Z,int is)const{return pZZsn.get(is,X,Y,Z);}


    double getJx(int X, int Y, int Z) const { return Jx.get(X,Y,Z);}
    double getJy(int X, int Y, int Z) const { return Jy.get(X,Y,Z);}
    double getJz(int X, int Y, int Z) const { return Jz.get(X,Y,Z);}
    arr3_double getJx() { return Jx; }
    arr3_double getJy() { return Jy; }
    arr3_double getJz() { return Jz; }
    arr4_double getJxs() { return Jxs; }
    arr4_double getJys() { return Jys; }
    arr4_double getJzs() { return Jzs; }

    double getJxs(int X,int Y,int Z,int is)const{return Jxs.get(is,X,Y,Z);}
    double getJys(int X,int Y,int Z,int is)const{return Jys.get(is,X,Y,Z);}
    double getJzs(int X,int Y,int Z,int is)const{return Jzs.get(is,X,Y,Z);}

    /*! get the electric field energy */
    double getEenergy();
    /*! get the magnetic field energy */
    double getBenergy();

    /*! fetch array for summing moments of thread i */
    Moments10& fetch_moments10Array(int i){
      assert_le(0,i);
      assert_lt(i,sizeMomentsArray);
      return *(moments10Array[i]);
    }
    int get_sizeMomentsArray() { return sizeMomentsArray; }

    /*! print electromagnetic fields info */
    void print(void) const;
    
    
    //get MPI Derived Datatype
    MPI_Datatype getYZFacetype(bool isCenterFlag){return isCenterFlag ?yzFacetypeC : yzFacetypeN;}
    MPI_Datatype getXZFacetype(bool isCenterFlag){return isCenterFlag ?xzFacetypeC : xzFacetypeN;}
    MPI_Datatype getXYFacetype(bool isCenterFlag){return isCenterFlag ?xyFacetypeC : xyFacetypeN;}
    MPI_Datatype getXEdgetype(bool isCenterFlag){return  isCenterFlag ?xEdgetypeC : xEdgetypeN;}
    MPI_Datatype getYEdgetype(bool isCenterFlag){return  isCenterFlag ?yEdgetypeC : yEdgetypeN;}
    MPI_Datatype getZEdgetype(bool isCenterFlag){return  isCenterFlag ?zEdgetypeC : zEdgetypeN;}
    MPI_Datatype getXEdgetype2(bool isCenterFlag){return  isCenterFlag ?xEdgetypeC2 : xEdgetypeN2;}
    MPI_Datatype getYEdgetype2(bool isCenterFlag){return  isCenterFlag ?yEdgetypeC2 : yEdgetypeN2;}
    MPI_Datatype getZEdgetype2(bool isCenterFlag){return  isCenterFlag ?zEdgetypeC2 : zEdgetypeN2;}
    MPI_Datatype getCornertype(bool isCenterFlag){return  isCenterFlag ?cornertypeC : cornertypeN;}



    MPI_Datatype getProcview(){return  procview;}
    MPI_Datatype getXYZeType(){return xyzcomp;}
    MPI_Datatype getProcviewXYZ(){return  procviewXYZ;}
    MPI_Datatype getGhostType(){return  ghosttype;}

    void freeDataType();
    bool isLittleEndian(){return lEndFlag;};

  public: // accessors
    const Collective& get_col()const{return _col;}
    const Grid& get_grid()const{return _grid;};
    const VirtualTopology3D& get_vct()const{return _vct;}
    /* ********************************* // VARIABLES ********************************* */
    
  private:
    // access to global data
    const Collective& _col;
    const Grid& _grid;
    const VirtualTopology3D&_vct;
    /*! light speed */
    double c;
    /* 4*PI for normalization */
    double FourPI;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    int SmoothNiter;
    /*! delt = c*th*dt */
    double delt;
    /*! number of particles species */
    int ns;
    /*! GEM challenge parameters */
    double B0x, B0y, B0z, delta;
    /** Earth Model parameters */
    double B1x, B1y, B1z;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! Boundary electron speed */
    double ue0, ve0, we0;


    // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
    /*! number of cells - X direction, including + 2 (guard cells) */
    int nxc;
    /*! number of nodes - X direction, including + 2 extra nodes for guard cells */
    int nxn;
    /*! number of cell - Y direction, including + 2 (guard cells) */
    int nyc;
    /*! number of nodes - Y direction, including + 2 extra nodes for guard cells */
    int nyn;
    /*! number of cell - Z direction, including + 2 (guard cells) */
    int nzc;
    /*! number of nodes - Z direction, including + 2 extra nodes for guard cells */
    int nzn;
    /*! local grid boundaries coordinate */
    double xStart, xEnd, yStart, yEnd, zStart, zEnd;
    /*! grid spacing */
    double dx, dy, dz, invVOL;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;
    /** source center - X direction   */
    double x_center;
    /** source center - Y direction   */
    double y_center;
    /** source center - Z direction   */
    double z_center;
    /** Characteristic length */
    double L_square;

    /*! PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes */
    array3_double PHI;

    // Electric field component used to move particles
    // organized for rapid access in mover_PC()
    // [This is the information transferred from cluster to booster].
    array4_pfloat fieldForPcls;

    // Electric field components defined on nodes
    //
    array3_double Ex;
    array3_double Ey;
    array3_double Ez;

    // implicit electric field components defined on nodes
    //
    array3_double Exth;
    array3_double Eyth;
    array3_double Ezth;

    // magnetic field components defined on central points between nodes
    //
    array3_double Bxc;
    array3_double Byc;
    array3_double Bzc;

    // magnetic field components defined on nodes
    //
    array3_double Bxn;
    array3_double Byn;
    array3_double Bzn;

    // *************************************
    // TEMPORARY ARRAY
    // ************************************
    /*!some temporary arrays (for calculate hat functions) */
    array3_double tempXC;
    array3_double tempYC;
    array3_double tempZC;
    array3_double tempXN;
    array3_double tempYN;
    array3_double tempZN;
    /*! other temporary arrays (in MaxwellSource) */
    array3_double tempC;
    array3_double tempX;
    array3_double tempY;
    array3_double tempZ;
    array3_double temp2X;
    array3_double temp2Y;
    array3_double temp2Z;
    /*! and some for MaxwellImage */
    array3_double imageX;
    array3_double imageY;
    array3_double imageZ;
    array3_double Dx;
    array3_double Dy;
    array3_double Dz;
    array3_double vectX;
    array3_double vectY;
    array3_double vectZ;
    array3_double divC;
    //array3_double arr;
    /* temporary arrays for summing moments */
    int sizeMomentsArray;
    Moments10 **moments10Array;

    // *******************************************************************************
    // *********** SOURCES **
    // *******************************************************************************

    /*! Charge density, defined on central points of the cell */
    array3_double rhoc;
    /*! Charge density, defined on nodes */
    array3_double rhon;
    /*! Implicit charge density, defined on central points of the cell */
    array3_double rhoh;
    /*! SPECIES: charge density for each species, defined on nodes */
    array4_double rhons;
    /*! SPECIES: charge density for each species, defined on central points of the cell */
    array4_double rhocs;

    // current density defined on nodes
    //
    array3_double Jx;
    array3_double Jy;
    array3_double Jz;

    // implicit current density defined on nodes
    //
    array3_double Jxh;
    array3_double Jyh;
    array3_double Jzh;

    // species-specific current densities defined on nodes
    //
    array4_double Jxs;
    array4_double Jys;
    array4_double Jzs;

    // magnetic field components defined on nodes
    //
    array3_double   Bx_ext;
    array3_double   By_ext;
    array3_double   Bz_ext;

    array3_double   Bx_tot;
    array3_double   By_tot;
    array3_double   Bz_tot;

    // external current, defined on nodes
    array3_double   Jx_ext;
    array3_double   Jy_ext;
    array3_double   Jz_ext;

    // pressure tensor components, defined on nodes
    array4_double pXXsn;
    array4_double pXYsn;
    array4_double pXZsn;
    array4_double pYYsn;
    array4_double pYZsn;
    array4_double pZZsn;

    /*! Field Boundary Condition
      0 = Dirichlet Boundary Condition: specifies the
          value on the boundary of the domain
      1 = Neumann Boundary Condition: specifies the value of
          derivative on the boundary of the domain
      2 = Periodic boundary condition */

    // boundary conditions for electrostatic potential
    //
    int bcPHIfaceXright;
    int bcPHIfaceXleft;
    int bcPHIfaceYright;
    int bcPHIfaceYleft;
    int bcPHIfaceZright;
    int bcPHIfaceZleft;

    /*! Boundary condition for electric field 0 = perfect conductor 1 = magnetic mirror */
    //
    // boundary conditions for EM field
    //
    int bcEMfaceXright;
    int bcEMfaceXleft;
    int bcEMfaceYright;
    int bcEMfaceYleft;
    int bcEMfaceZright;
    int bcEMfaceZleft;


    /*! GEM Challenge background ion */
    double *rhoINIT;
    /*! Drift of the species */
    bool *DriftSpecies;

    /*! boolean for divergence cleaning */
    bool PoissonCorrection;
    /*! RESTART BOOLEAN */
    int restart1;

    /*! CG tolerance criterium for stopping iterations */
    double CGtol;
    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;


    //MPI Derived Datatype for Center Halo Exchange
    MPI_Datatype yzFacetypeC;
    MPI_Datatype xzFacetypeC;
    MPI_Datatype xyFacetypeC;
    MPI_Datatype xEdgetypeC;
    MPI_Datatype yEdgetypeC;
    MPI_Datatype zEdgetypeC;
    MPI_Datatype xEdgetypeC2;
    MPI_Datatype yEdgetypeC2;
    MPI_Datatype zEdgetypeC2;
    MPI_Datatype cornertypeC;

    //MPI Derived Datatype for Node Halo Exchange
    MPI_Datatype yzFacetypeN;
    MPI_Datatype xzFacetypeN;
    MPI_Datatype xyFacetypeN;
    MPI_Datatype xEdgetypeN;
    MPI_Datatype yEdgetypeN;
    MPI_Datatype zEdgetypeN;
    MPI_Datatype xEdgetypeN2;
    MPI_Datatype yEdgetypeN2;
    MPI_Datatype zEdgetypeN2;
    MPI_Datatype cornertypeN;
    
    //for VTK output
    MPI_Datatype  procviewXYZ,xyzcomp,procview,ghosttype;
    bool lEndFlag;
    
    void OpenBoundaryInflowB(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,

      int nx, int ny, int nz);
    void OpenBoundaryInflowE(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
      int nx, int ny, int nz);
    void OpenBoundaryInflowEImage(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
      const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
      int nx, int ny, int nz);
};

inline void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        rhons[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of charge density to current density - direction X to current density field on the node */
inline void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jxs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Y to current density field on the node */
inline void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jys[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Z to current density field on the node */
inline void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        Jzs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XX to current density field on the node */
inline void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXXsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XY to current density field on the node */
inline void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XZ to current density field on the node */
inline void EMfields3D::addPxz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YY to current density field on the node */
inline void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YZ to current density field on the node */
inline void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction ZZ to current density field on the node */
inline void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pZZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}

inline void get_field_components_for_cell(
  const double* field_components[8],
  const_arr4_double fieldForPcls,
  int cx,int cy,int cz)
{
  // interface to the right of cell
  const int ix = cx+1;
  const int iy = cy+1;
  const int iz = cz+1;

  // is this faster?
  //
  //field_components[0] = fieldForPcls[ix][iy][iz]; // field000
  //field_components[1] = fieldForPcls[ix][iy][cz]; // field001
  //field_components[2] = fieldForPcls[ix][cy][iz]; // field010
  //field_components[3] = fieldForPcls[ix][cy][cz]; // field011
  //field_components[4] = fieldForPcls[cx][iy][iz]; // field100
  //field_components[5] = fieldForPcls[cx][iy][cz]; // field101
  //field_components[6] = fieldForPcls[cx][cy][iz]; // field110
  //field_components[7] = fieldForPcls[cx][cy][cz]; // field111
  //
  // or is this?
  //
  // creating these aliases seems to accelerate this method (by about 30%?)
  // on the Xeon host processor, suggesting deficiency in the optimizer.
  //
  arr3_double_get field0 = fieldForPcls[ix];
  arr3_double_get field1 = fieldForPcls[cx];
  arr2_double_get field00 = field0[iy];
  arr2_double_get field01 = field0[cy];
  arr2_double_get field10 = field1[iy];
  arr2_double_get field11 = field1[cy];
  field_components[0] = field00[iz]; // field000 
  field_components[1] = field00[cz]; // field001 
  field_components[2] = field01[iz]; // field010 
  field_components[3] = field01[cz]; // field011 
  field_components[4] = field10[iz]; // field100 
  field_components[5] = field10[cz]; // field101 
  field_components[6] = field11[iz]; // field110 
  field_components[7] = field11[cz]; // field111 
}

typedef EMfields3D Field;

#endif
