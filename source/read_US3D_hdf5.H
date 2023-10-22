#pragma once
void readGridArraySizes(const char* gridfile, int* nn, int* nf, int* nc, int* nz);
void readGridIntegerArray(const char* gridfile, const char* dataName, std::vector<int>& data);
void readGridDoubleArray(const char* gridfile, const char* dataName, std::vector<double>& data);
std::vector<double> readUS3DSolutionFile(const char* datafile, const char* dataName);
void readUS3DSolutionNSV(const char* datafile, int* nsv);
int getUS3DSolutionNSV(const char* datafile);
int getUS3DSolutionNBV(const char* datafile);
std::vector<std::string> readUS3DDataBVnames(const char* datafile);
std::vector<std::string> readUS3DDataSVnames(const char* datafile);