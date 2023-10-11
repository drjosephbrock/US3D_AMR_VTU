// #include <vtkDataSetMapper.h>
// #include <vtkIdList.h>
// #include <vtkProperty.h>
// #include <vtkUnstructuredGrid.h>
// #include <string>
// #include <list>
#include "write_us3d_vtu.H"
#include <algorithm>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <iostream>
#include <thread>



std::list<int> appendAndUniqify(std::list<int>& inputList, const std::list<int>& valuesToAppend) {
    // Append values to the input list
    inputList.insert(inputList.end(), valuesToAppend.begin(), valuesToAppend.end());

    // Sort the list to bring duplicates together
    inputList.sort();

    // Remove duplicates
    inputList.unique();

    return inputList;
}

void write_VTUFile(const char* filename,
                  bool writeAscii,
                  const vtkNew<vtkUnstructuredGrid>& ugrid){

  // Here we write out the cube.
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetInputData(ugrid);
  writer->SetFileName(filename);
  if (writeAscii){
    writer->SetDataModeToAscii();
  }
  writer->Update();
}

void add_cell(vtkNew<vtkUnstructuredGrid>& ugrid, const int& cellID, 
          const std::vector<int>& ief,
          const std::vector<int>& ifn){
  

  int node, face;
  int nFaceNodes;
  int ncellFaces = ief[25*cellID];
  std::list<int> CellNodes;
  vtkNew<vtkPoints> points;
  vtkNew<vtkIdList> faces;
  vtkIdType face3[3] = {-1,  -1,  -1};
  vtkIdType face4[4] = {-1,  -1,  -1,  -1};
  vtkIdType face5[5] = {-1,  -1,  -1,  -1, -1};
  vtkIdType face6[6] = {-1,  -1,  -1,  -1, -1, -1};
  vtkIdType face7[7] = {-1,  -1,  -1,  -1, -1, -1, -1};
  vtkIdType face8[8] = {-1,  -1,  -1,  -1, -1, -1, -1, -1};

  auto addFace = [&](const vtkIdType* face, int size) {
      faces->InsertNextId(size);
      for (int i = 0; i < size; ++i) {
          faces->InsertNextId(face[i]);
      }
  };


  for (int cellface = 0; cellface < ncellFaces; cellface++){
    face = ief[25*cellID + cellface + 1]-1;
    nFaceNodes = ifn[9*face];
    // faces->InsertNextId(nFaceNodes);
    for (int faceNode = 0; faceNode < nFaceNodes; faceNode++){
      node = ifn[9*face + faceNode + 1]-1;
      CellNodes.push_back(node);
      // faces->InsertNextId(face[faceNode])
      if (nFaceNodes == 4){
        face4[faceNode] = node;
      } else if (nFaceNodes == 5){
        face5[faceNode] = node;
      } else if (nFaceNodes == 6){
        face6[faceNode] = node;
      } else if (nFaceNodes == 7){
        face7[faceNode] = node;
      } else if (nFaceNodes == 8){
        face8[faceNode] = node;
      }
    }
    if (nFaceNodes == 4){
      addFace(face4, 4);
    } else if (nFaceNodes == 5){
      addFace(face5, 5);
    } else if (nFaceNodes == 6){
      addFace(face6, 6);
    } else if (nFaceNodes == 7){
      addFace(face7, 7);
    } else if (nFaceNodes == 8){
      addFace(face8, 8);
    }
  }

  std::list<int> UniqueCellNodes;
  appendAndUniqify(UniqueCellNodes, CellNodes);
  vtkIdType pointIds[UniqueCellNodes.size()];
  int index = 0;
  for (auto it = UniqueCellNodes.begin(); it != UniqueCellNodes.end(); ++it) {
      pointIds[index++] = *it;
  }
  ugrid->InsertNextCell(VTK_POLYHEDRON, UniqueCellNodes.size(), pointIds, ncellFaces, faces->GetPointer(0));
}

void add_points_parallel(vtkPoints* points, const std::vector<double>& xcn, int start, int end){
    for (int i = start; i < end; i++){
        // Calculate the indices once to avoid repeated calculations
        int index = 3 * i;
        double x = xcn[index];
        double y = xcn[index + 1];
        double z = xcn[index + 2];
      points->InsertNextPoint(x, y, z);
    }
}
void add_points(vtkNew<vtkUnstructuredGrid>& ugrid, const std::vector<double>& xcn, const int& nn ){
    std::cout <<"add_points "<< nn << std::endl;

    vtkNew<vtkPoints> points;
    points->Allocate(nn);

    // number of threads to use
    int num_threads = std::thread::hardware_concurrency();

    // Create and launch threads
    std::vector<std::thread> threads;
    int points_per_thread = nn / num_threads;
    int start = 0;
    
    for (int i = 0; i < num_threads; i++){
      int end = (i == num_threads - 1)? nn : start + points_per_thread;
      threads.push_back(std::thread(add_points_parallel, points.GetPointer(), std::ref(xcn), start, end));
      start = end;
    }

    // Join the threads
    for (auto& thread : threads){
      thread.join();
    }

    ugrid->SetPoints(points);
}

void add_cell_data(vtkNew<vtkUnstructuredGrid>& ugrid, const std::vector<double>& solution, const char* variables[], const int& nv, const int& nel) {
    for (int var = 0; var < nv; var++) {
        vtkNew<vtkDoubleArray> scalars = vtkNew<vtkDoubleArray>();
        scalars->SetName(variables[var]);
        std::cout << "Writing '" << variables[var] << "' variable" << std::endl;
        for (int cell = 0; cell < nel; cell++) {
            double value = solution[nv * cell + var];
            scalars->InsertNextValue(value);
        }
        ugrid->GetCellData()->AddArray(scalars);
    }
}
