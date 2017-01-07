#pragma once
/*********************************************************\
 * ignore specific warnings                              *
\*********************************************************/
#pragma warning(disable : 4996)
#pragma warning(disable : 4244)

/*********************************************************\
 * define header file                                    *
\*********************************************************/
#ifndef Correlation_h
#define Correlation_h

/*********************************************************\
 * define variable                                       *
\*********************************************************/
#define MAX_ARRAY_SIZE  1000

/*********************************************************\
 * include libraries                                     *
\*********************************************************/
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include <iostream>
#include <sstream>
#include "cpl_conv.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"
#include <Eigen/Dense>

/*********************************************************\
 * use namespaces                                        *
\*********************************************************/
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace tbb;
using namespace std;

/*********************************************************\
 * declare functions                                     *
\*********************************************************/
int parse_command_file_char(FILE*, char*, char*);
int main(int, char*[]);
long ld2l(long double);
FILE * openTextfile_read(char*);
int parse_command_file_char(FILE*, char*, char*);
GDALDataset* openImagefile(char*);
OGRGeometry* buildGeometry (double, double, double, double, OGRSpatialReference *);
OGREnvelope* getOverlap(double, double, double, double, double, double, double, double);
float*** arrayOfImage3D(GDALDataset*, OGREnvelope*, double[]);

/*********************************************************\
 * end of source code                                    *
\*********************************************************/
#endif