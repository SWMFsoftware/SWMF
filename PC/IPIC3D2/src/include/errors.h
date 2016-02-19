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

#ifndef ipic_errors_H
#define ipic_errors_H

//void errmsg_printf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
//void eprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
//void Wprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void eprintf_fileLine(FILE * fptr, const char *type,
  const char *func, const char *file, int line_number,
  const char *format, ...);
void fprintf_fileLine(FILE * fptr,
  const char *type, const char *func, const char *file, int line_number,
  const char *format, ...);

#define eprintf(args...) \
  eprintf_fileLine(stdout,"ERROR",__func__, __FILE__, __LINE__, ## args);
#define error_printf(args...) \
  eprintf_fileLine(stdout,"ERROR",__func__, __FILE__, __LINE__, ## args);
//#define eprintf(args...) \
//  eprintf_fileLine("ERROR",__func__, __FILE__, __LINE__, ## args);
#define warning_printf(args...) \
  fprintf_fileLine(stdout,"WARNING",__func__, __FILE__, __LINE__, ## args);
#define declare_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val);
declare_invalid_value_error(double);
declare_invalid_value_error(int);
declare_invalid_value_error(const char*);
#define unsupported_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "unsupported", #val, val);
#define invalid_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "invalid", #val, val);

#endif
