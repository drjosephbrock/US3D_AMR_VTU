#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <string>
#include "read_US3D_hdf5.H"
#include "write_us3d_vtu.H"

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


int main(int, char*[])
{

  // const char* gridfile = "/h/group/FILESHARE/amr_grid.h5";
  // const char* gridfile = "../../pre-amr/test/amr_grid.h5";
  // const char* gridfile = "../../2cell_amr_grid.h5";
  // const char* gridfile = "../../2d_msl_amr_grid.h5";
  const char* gridfile = "amr_grid.h5";
  const char* datafile = "data.h5";
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
    add_points(ugrid, xcn, nn);
    // for (int i = 0; i < nn; i++){
    //   points->InsertNextPoint(xcn[3*i], xcn[3*i + 1], xcn[3*i + 2]);
    // }
    // ugrid->SetPoints(points);
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
    std::cout << std::endl;
  }

  // Read Data and append to file
  {
    int nv, nel;
    const char* varnames[6] = {"rho","u","v","w","T","res"};
    auto solution = readUS3DSolutionFile(datafile, "/solution/run_1/interior", nv, nel);
    // std::cout << solution.size()/6 << std::endl;
    // std::cout << std::endl<< nv << ", " << nel << std::endl;
    add_cell_data(ugrid, solution, varnames, 6, nc);
  }

  std::cout << "writing file" << std::endl;
  write_VTUFile("amr_mesh.vtu", 0, ugrid);
  std::cout << "Complete" << std::endl;

  return EXIT_SUCCESS;
}
