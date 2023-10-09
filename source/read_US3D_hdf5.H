#pragma once
std::tuple<int, int, int, int>  readGridArraySizes(const char* gridfile);
std::tuple< std::vector<double>, std::vector<int>, std::vector<int> > readGridArrays(const char* gridfile);
std::vector<int> readGridIntegerArray(const char* gridfile, const char* dataName);
std::vector<double> readGridDoubleArray(const char* gridfile, const char* dataName);