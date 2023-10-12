#pragma once
#include <list>
#include <vector>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>

std::list<int> appendAndUniqify(std::list<int>& inputList, const std::list<int>& valuesToAppend);


void write_VTUsurfaceFile(const char* filename, bool writeAscii, vtkNew<vtkUnstructuredGrid>& ugrid);
void add_polygon(vtkNew<vtkUnstructuredGrid>& ugrid, const int& faceID, const std::vector<int>& ifn);
void add_wall_faces(const char* gridfile, const char* datafile, const int& izn, vtkNew<vtkUnstructuredGrid>& ugrid);
void add_points_poly(vtkNew<vtkUnstructuredGrid>& ugrid, const std::vector<double>& xcn, const int& nn);

void write_VTUFile(const char* filename, bool writeAscii, const vtkNew<vtkUnstructuredGrid>& ugrid);
void add_cell_chunk(vtkNew<vtkUnstructuredGrid>& ugrid);
void add_cell(vtkNew<vtkUnstructuredGrid>& ugrid, const int& cellID,  const std::vector<int>& ief, const std::vector<int>& ifn);
void add_points(vtkNew<vtkUnstructuredGrid>& ugrid, const std::vector<double>& xcn, const int& nn );
void add_cell_data(vtkNew<vtkUnstructuredGrid>& ugrid, const std::vector<double>& solution, const char* variables[], const int& nv, const int& nel );