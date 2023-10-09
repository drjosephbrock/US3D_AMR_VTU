#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <string>
#include "read_US3D_hdf5.H"
#include <list>
#include <algorithm>

// Function to display a progress bar
void displayProgressBar(int progress, int total) {
    const int barWidth = 50;
    float fraction = static_cast<float>(progress) / total;
    int numBars = static_cast<int>(fraction * barWidth);

    std::cout << "[";

    for (int i = 0; i < barWidth; ++i) {
        if (i < numBars) {
            std::cout << "=";
        } else {
            std::cout << " ";
        }
    }

    std::cout << "] " << int(fraction * 100.0) << "%\r";
    std::cout.flush();
}

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

int main(int, char*[])
{

  // const char* gridfile = "/h/group/FILESHARE/amr_grid.h5";
  // const char* gridfile = "../../pre-amr/test/amr_grid.h5";
  // const char* gridfile = "../../2cell_amr_grid.h5";
  // const char* gridfile = "../../2d_msl_amr_grid.h5";
  const char* gridfile = "amr_grid.h5";
  vtkNew<vtkPoints> points;
  vtkNew<vtkUnstructuredGrid> ugrid;

  
  auto [nn, nf, nc, nz] = readGridArraySizes(gridfile);
  printf("  --  Nodes: %d\n", nn);
  printf("  --  Faces: %d\n", nf);
  printf("  --  Cells: %d\n", nc);
  printf("  --  Zones: %d\n", nz);
  
  // Calculate memory requirements for each vector
  // Assuming sizeof(int) = 4 bytes and sizeof(double) = 8 bytes
  double sizeof_int = 4.0;  // bytes
  double sizeof_double = 8.0;  // bytes
  double ief_memory = (double)nc * 25 * sizeof_int;
  double ifn_memory = (double)nf * 8 * sizeof_int;
  double xcn_memory = nn * 3 *sizeof_double;
  printf(" -- ief mem=%e GB\n", ief_memory/ (1024*1024*1024));
  printf(" -- ifn mem=%e GB\n", ifn_memory/ (1024*1024*1024));
  printf(" -- xcn mem=%e GB\n", xcn_memory/ (1024*1024*1024));

  // Append points list
  {
    auto xcn =readGridDoubleArray(gridfile, "/xcn");
    int nn = xcn.size()/3;
    for (int i = 0; i < nn; i++){
      points->InsertNextPoint(xcn[3*i], xcn[3*i + 1], xcn[3*i + 2]);
    }
    ugrid->SetPoints(points);
  }

  // Add face/cell connectivity data
  {
    auto ief = readGridIntegerArray(gridfile, "/iefpoly");
    auto ifn = readGridIntegerArray(gridfile, "/ifnpoly");
    int ncellFaces, nFaceNodes, face, node;
    std::list<int> CellNodes;

    // Loop over cells
    for (int cell = 0; cell < nc; cell ++){
      displayProgressBar(cell, nc);
      add_cell(ugrid, cell, ief, ifn);
    }
  }

  std::cout << "writing file" << std::endl;
  write_VTUFile("amr_mesh.vtu", 0, ugrid);
  std::cout << "Complete" << std::endl;

  return EXIT_SUCCESS;
}
