/*******************************************************************************************
  debug.h  -  definitions for debug and diagnostics (written by Alec Johnson)
 ********************************************************************************************/
#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <cstdarg>
#include <cstdio>

#include "errors.h"

#ifdef BATSRUS
#include <string>
#include <iostream>
#include <fstream>
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid3DCU.h"
#endif

void fprintf_fileLine(FILE * fptr, const char *type, const char *func,
  const char *file, int line_number, const char *format, ...);

#define dprintf(args...) fprintf_fileLine(stdout, "DEBUG", __func__, __FILE__, __LINE__,## args)
#define dprint(var) printvar_fileLine(__func__, __FILE__,__LINE__,#var,var);
#define declare_dprintvar_fileLine(type) \
void printvar_fileLine(const char*,const char*,int,const char*,type);

declare_dprintvar_fileLine(int);
declare_dprintvar_fileLine(long long);
declare_dprintvar_fileLine(double);
declare_dprintvar_fileLine(const char *);
declare_dprintvar_fileLine(const void *);

#ifdef BATSRUS
extern Collective * _col0;
extern Grid3DCU * _grid0;
extern VCtopology3D * _vct0;
extern string testFuncs;
extern int iTest,jTest,kTest;
extern int iProc;
extern ofstream outfile;
extern string filename;

void init_debug_SWMF(Collective *col, Grid3DCU *grid,
			    VCtopology3D *vct,
			    string testFuncsIn,
			    int iTestIn, int jTestIn, int kTestIn);
void finalize_debug_SWMF();
bool do_test_func(ofstream *&outfileOut, string func);
bool do_test_cell(ofstream *&outfileOut,
			 int i=-1, int j=-1, int k=-1);
#endif

#endif
