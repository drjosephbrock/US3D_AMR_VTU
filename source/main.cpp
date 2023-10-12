#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <string>
#include <cstring>
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


int main(int argc, char* argv[])
{
  bool plotSolution = true;
  for (int i = 1; i < argc; i++){
     if (std::strcmp(argv[i], "-plot_solution=false") == 0) {
        plotSolution = false;
        break;
     }
  }
  if (plotSolution){
     std::cout << "Plotting mesh and solution" << std::endl;
  } else{
     std::cout << "Plotting mesh only" << std::endl;
  }

  const char* gridfile = "grid.h5";
  const char* datafile = "data.h5";
  vtkNew<vtkPoints> points;
  vtkNew<vtkUnstructuredGrid> ugrid;

  int nn, nf, nc, nz;
  readGridArraySizes(gridfile, &nn, &nf, &nc, &nz);

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
    printf("  --  reading points... ");
    std::vector<double> xcn(3*nn);
    readGridDoubleArray(gridfile, "/xcn", xcn);
    add_points(ugrid, xcn, nn);
    printf(" done\n");
  }
  add_cell_chunk(ugrid);

  // // Add face/cell connectivity data
  // {
  //   std::vector<int> ief(25*nc);
  //   std::vector<int> ifn(9*nf);
  //   readGridIntegerArray(gridfile, "/iefpoly", ief);
  //   readGridIntegerArray(gridfile, "/ifnpoly", ifn);

  //   // Loop over cells
  //   for (int cell = 0; cell < nc; cell ++){
  //     displayProgressBar(cell, nc);
  //     add_cell(ugrid, cell, ief, ifn);
  //   }
  //   std::cout << std::endl;
  // }

  // // Read Data and append to file
  // if (plotSolution){
  //   int nv = 6;
  //   const char* varnames[6] = {"rho","u","v","w","T","res"};
  //   // std::vector<double> solution(nv*nc);
  //   auto solution = readUS3DSolutionFile(datafile, "/solution/run_1/interior", nv, nc);
  //   // std::cout << solution.size()/6 << std::endl;
  //   // std::cout << std::endl<< nv << ", " << nel << std::endl;
  //     add_cell_data(ugrid, solution, varnames, 6, nc);
  // }

  std::cout << "writing file" << std::endl;
  write_VTUFile("amr_mesh.vtu", 0, ugrid);
  std::cout << "Complete" << std::endl;

  return EXIT_SUCCESS;
}
