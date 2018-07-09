/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/ts_serial.hpp
 *
 *  \brief Header file for ts_serial.cpp
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ems.c 5832 2018-06-26 23:49:51Z riz008 $
 */

#ifndef _TS_SERIAL_HPP
#define _TS_SERIAL_HPP

extern "C" {
#include "ems.h"
}

using namespace std;

/*
 * Main class
 */
class tsSerial
{
public:
  /* Default constructor */
  tsSerial(char *fname);

  /* Destructor */
  ~tsSerial(void);

  /* Sparsify function */
  void serialise(char *vname);

  inline timeseries_t *get_ts(void) {
    return(ts);
  }

  /* Getters */
  inline int getSize(void) {
    return(s2x.size());
  }

  inline double getX(int s) {
    return(s2x[s]);
  }

  inline double getY(int s) {
    return(s2y[s]);
  }

  inline int getC0(int s) {
    return(s2c0[s]);
  }

  inline int getC1(int s) {
    return(s2c1[s]);
  }

private:
  timeseries_t *ts;

  /* Geographical data */
  vector<double> s2x;
  vector<double> s2y;

  /* Index space */
  vector<int> s2c0;
  vector<int> s2c1;
};


#endif /* _TS_SERIAL_HPP */
