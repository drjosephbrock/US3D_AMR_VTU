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

  // Append points list
  {
    printf("  --  reading points... ");
    std::vector<double> xcn(3*nn);
    readGridDoubleArray(gridfile, "/xcn", xcn);
    add_points(ugrid, xcn, nn);
    printf(" done\n");
  }
  add_cell_chunk(ugrid);

  std::cout << "writing file" << std::endl;
  write_VTUFile("amr_mesh.vtu", 0, ugrid);
  std::cout << "Complete" << std::endl;

  return EXIT_SUCCESS;
}
