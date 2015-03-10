// see naming correspondence 
// for description see the reference provided below

// the parameters below are set at initialization
double  EpsMin, EpsMax, SqrtEpsMin_x_5;
double  LocalTimeStepMin;
int     nIterCorrector;
double* NDensityMin;


// auxilary function, see meaning in the reference below
inline double _alpha(double ratio){
  return (180. + ratio * (60. + ratio * (11. + ratio))) / 
         (360. + ratio * (60. + ratio * (12. + ratio)));
}


void chemeq2(double GlobalTimeStep, 
	     double (*get_reaction_rates)(double*, double*, double*), 
	     int nSpec, 
	     int *NDensity){
  /*
    code taken from
    David R. Mott, Elaine S. Oran,
    CHEMEQ2: A Solver for the Stiff Ordinary Differential Equations 
    of Chemical Kinetics,
    NRL Memorandum Report (Naval Research Laboratory, Washington, DC)

    Translated from FORTRAN77,  naming correspondence
      ORIGINAL---HERE---------------PAPER      
      **** INPUT ARGUMENTS **************************
      dtg--------GlobalTimeStep-----\Delta t_g
      gsub-------get_reaction_rates-
      n----------nSpec--------------n
      y----------NDensity-----------y
      **** MAJOR PARAMETERS *************************
      q----------ProdRate-----------q,   Eq.(1)
      d----------LossRate-----------p y, Eq.(1)
      scrarray---FullRate-----------
      rtau-------TimeScaleRatio-----\Delta t / \tau
      **** CONTROL PARAMETERS ***********************
      tn---------LocalTime----------t - t^0
      dt---------LocalTimeStep------\Delta t
      dtmin------LocalTimeStepMin---
      iter-------iIterCorrector-----
      itermax----nIterCorrector-----
      eps--------Eps----------------\sigma,   Eq.(48)
      epsmax-----EpsMax-------------c,        Eq.(47)
      epsmin-----EpsMin-------------\epsilon, Eq.(46)
      
    integration method: alpha-QSS as described in
    D.R.Mott, E.S.Oran, B van Leer,
    A Quasi-Steady-State Solver for the Stiff Ordinary Differential Equations
    of Reaction Kinetics,
    Journal of Computational Physics 164, 407-428, 2000
   */

  // local (in-solver) time & time step
  double LocalTime = 0.0;
  double LocalTimeStep;
  
  int iIterCorrector;

  // lists of reaction rates 
  static double* ProdRate;
  static double* LossRate;
  static double* FullRate;
  // starting values and first prediction of concentrations
  static double* NDensityStart, NDensityPred;
  // ratio of the local time step to reaction time scale
  static double* TimeScaleRatio, TimeScaleRatioStart;

  // init == 0 => the 1st call of the function 
  static int init = 0;

  // allocate arrays at the first call
  if(init == 0){
    ProdRate           = new double [nSpec];
    LossRate           = new double [nSpec];
    FullRate           = new double [nSpec];
    NDensityStart      = new double [nSpec];
    NDensityPred       = new double [nSpec];
    ProdRateStart      = new double [nSpec];
    TimeScaleRatio     = new double [nSpec];
    TimeScaleRatioStart= new double [nSpec];
    init++;
  }


  // get production and loss rate at initial values
  (*get_reaction_rates)(NDensity, ProdRate, LossRate);

  {// set initial value for time step
    double scrtch, ascr, scr0, scr1, scr2;
    for(int iSpec = 0; iSpec < nSpec; iSpec++){
      ascr   = abs(ProdRate[iSpec]);
      scr2   = (0.1*EpsMin*ascr - LossRate[iSpec] >= 0) ? 1.0 : -1.0;
      scr2  /= NDensity[iSpec];
      scr1   = scr2 * Lossrate[iSpec];
      scr0   =-scr2 * abs(ascr - LossRate[iSpec]);
      scrtch = (scrtch > scr1) ? scrtch : src1;
      scrtch = (scrtch > scr0) ? scrtch : src0;
    }
    LocalTimeStep = (SqrtEpsMin_x_5 < scrtch) ? SqrtEpsMin_x_5 : scrtch;
  }

  // find time scale ratios and save starting values
  for(int iSpec = 0; iSpec < nSpec; iSpec++){
    TimeScaleRatio[     iSpec]= LocalTimeStep * LossRate / NDensity[iSpec];
    TimeScaleRatioStart[iSpec]= TimeScaleRatio[iSpec];
    NDensityStart[      iSpec]= NDensity[iSpec];
    ProdRateStart[      iSpec]= ProdRate[iSpec]; 
  }

  // integrate up to GlobalTimeStep
  while(GlobalTimeStep > LocalTime){
    // prediction
    for(int iSpec = 0; iSpec < nSpec; iSpec++){
      double rtaui    = TimeScaleRatio[iSpec]; 
      double alpha    = _alpha(rtaui);
      FullRate[iSpec] = (ProdRate[iSpec]-LossRate[iSpec]) / (1. + alpha*rtaui);
    }
    
    // correction, may be several iterations
    for(iIterCorrector = 0; iIterCorrector < nIterCorrector; iIterCorrector++){
      // compute predictor values
      for(int iSpec = 0; iSpec < nSpec; iSpec++){
	NDensity[iSpec] = NDensityStart[iSpec]+LocalTimeStep * FullRate[iSpec];
	if(NDensity[iSpec] < NDensityMin[iSpec])
	  NDensity[iSpec] = NDensityMin[iSpec];
      }
      
      // save predictor values of 1st iteration
      if(iIterCorrector == 0){
	LocalTime += LocalTimeStep;
	for(int iSpec = 0; iSpec < nSpec; iSpec++){
	  NDensityPred[iSpec] = NDensity[iSpec];
	}
      }
      
      // get rates at predictor value
      (*get_reaction_rates)(NDensity, ProdRate, LossRate);
      
      // compute corrector value
      for(int iSpec = 0; iSpec < nSpec; iSpec++){
	double rtaub = 
	  0.5*(TimeScaleRatioStart[iSpec] + 
	       LocalTimeStep * LossRate[iSpec] / NDensity[iSpec]);
	double alpha = _alpha(rtaub);
	double qt = 
	  ProdRateStart[iSpec] * (1. - alpha) +
	  ProdRate[     iSpec] * alpha;
	FullRate[iSpec] = 
	  (qt - rtaub / LocalTimeStep * NDensityStart[iSpec]) / 
	  (1. + alpha * rtaub);
      }
    }
 
    double Eps = 1E-10;   
    {// compute final value of NDensity, compute correction term Eps
      double scr0, scr1, scr2;
      for(int iSpec = 0; iSpec < nSpec; iSpec++){
	scr2 = NDensityStart[iSpec] + LocalTimeStep * FullRate[iSpec];
	scr2 = (scr2 < 0.0) ? 0.0 : scr2;
	scr1 = abs(scr2 - NDensityPred[iSpec]);
	NDensity[iSpec] = (scr2>NDensityMin[iSpec])? scr2 : NDensityMin[iSpec];
	if(0.25 * (NDensityStart[iSpec] + NDensity[iSpec]) > NDensity[iSpec]){
	  scr1/= NDensity[iSpec];
	  scr0 = 
	    abs(ProdRate[iSpec] - LossRate[iSpec]) /
	    (ProdRate[iSpec] + LossRate[iSpec] + 1E-30);
	  scr0 = 0.5 * ( ((scr0 < scr1) ? scr0 : scr1) + scr1);
	  Eps = (scr0 > Eps) ? scr0 : Eps;
	}
      }
    }
    Eps /= EpsMin;
    
    // correction Eps exceeds EpsMax => reset time and repeat
    if(Eps <= EpsMax){
      if(LocalTime >= GlobalTimeStep) return;
    }
    else LocalTime = 0.0;
    
    // time step modification
    double dto = LocalTimeStep;
    LocalTimeStep *= 0.005 + pow(Eps,-0.5);
    if(LocalTimeStep < LocalTimeStepMin)
      LocalTimeStep = LocalTimeStepMin;
    if(LocalTimeStep > GlobalTimeStep - LocalTime)
      LocalTimeStep = GlobalTimeStep - LocalTime;
    
    if(Eps > EpsMax){
      dto = LocalTimeStep / dto;
      
      for(int iSpec = 0; iSpec < nSpec; iSpec++){
	TimeScaleRatioStart[iSpec] *= dto;
      } 
      else {
	// get production and loss rate at initial values
	(*get_reaction_rates)(NDensity, ProdRate, LossRate);
	
	// find time scale ratios and save starting values
	for(int iSpec = 0; iSpec < nSpec; iSpec++){
	  TimeScaleRatio[iSpec]     = LocalTimeStep *LossRate /NDensity[iSpec];
	  TimeScaleRatioStart[iSpec]= TimeScaleRatio[iSpec];
	  NDensityStart[iSpec]      = NDensity[iSpec];
	  ProdRateStart[iSpec]      = ProdRate[iSpec]; 
	}
      }
    }
  }
}
