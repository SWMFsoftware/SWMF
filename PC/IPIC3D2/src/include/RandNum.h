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

#ifndef RANDNUM_H
#define RANDNUM_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double PseudoRand(long *idum){

  // "minimal" random number generator of Park and Miller. Taken from
  // "numerical recipes in C"
  // Return a uniform random deviate between 0.0 and 1.0. Set or reset idum 
  // to any integer value (except the unlikely value MASK) it initialize 
  // the sequence; idum must not be alterd between calls or successive 
  // deviates in a sequence.
  
  long k;
  double ans;
  
  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum <0) *idum +=IM;
  ans=AM*(*idum);
  *idum ^=MASK;
  return ans;
  
};

// cleaning
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

#endif
