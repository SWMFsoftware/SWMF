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

#ifndef _PSK_HDF5_ADAPTOR_H_
#ifndef NO_HDF5

#define _PSK_HDF5_ADAPTOR_H_

#include "PSKOutput.h"
#include <algorithm>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "arraysfwd.h"

namespace PSK {

  class HDF5OutputAdaptor:public OutputAdaptor {

    std::string _hdf5_file_name;
    hid_t _hdf5_file_id;

    static std::string purify_object_name(const std::string & objname);
    static void split_name(const std::string & name, std::vector < std::string > &elements);

    void get_dataset_context(const std::string & name, std::vector < hid_t > &hid_array, std::string & dataset_name);

  public:
      HDF5OutputAdaptor(void) {}
      void open(const std::string & name);
    void open_append(const std::string & name);
    void close(void);

    void write(const std::string & tag, int i_value);
    void write(const std::string & tag, long i_value);
    void write(const std::string & tag, const Dimens dimens, const int *i_array);
    void write(const std::string & tag, const Dimens dimens, const long *i_array);
    void write(const std::string & tag, const Dimens dimens, const longid *i_array);

    void write(const std::string & tag, const Dimens dimens, const std::vector < int >&i_array);

    void write(const std::string & objname, const Dimens dimens, const int ***i_array);


    // write float functions
    void write(const std::string & objname, float f);
    void write(const std::string & objname, const Dimens dimens, const float *f_array);
    void write(const std::string & objname, const Dimens dimens, const std::vector < float >&f_array);
    void write(const std::string & objname, const Dimens dimens, const float ***f_array);

    // write double functions
    void write(const std::string & objname, double d);
    void write(const std::string & objname, const Dimens dimens, const double *d_array);
    void write(const std::string & objname, const Dimens dimens, const std::vector < double >&d_array);
    void write(const std::string & objname, const Dimens dimens, const_arr3_double d_array);
    void write(const std::string & objname, const Dimens dimens, const int i, const_arr4_double d_array);

    void write(const std::string & objname, const Dimens dimens, double **d_array);

    void write(const std::string & objname, const Dimens dimens, const int i, const_arr3_double d_array);

  };

}                               // namespace
#endif
#endif
