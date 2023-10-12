#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
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

  int izn = 3;
  const char* gridfile = "grid.h5";
  const char* datafile = "data.h5";
  vtkNew<vtkUnstructuredGrid> ugrid;
  vtkNew<vtkPolyData> polygonPolyData;


  int nn, nf, nc, nz;
  readGridArraySizes(gridfile, &nn, &nf, &nc, &nz);

  printf("  --  Nodes: %d\n", nn);
  printf("  --  Faces: %d\n", nf);
  printf("  --  Cells: %d\n", nc);
  printf("  --  Zones: %d\n", nz);

  // Volume mesh
  // // Append points list
  // {
  //   printf("  --  reading points... ");
  //   std::vector<double> xcn(3*nn);
  //   readGridDoubleArray(gridfile, "/xcn", xcn);
  //   add_points(ugrid, xcn, nn);
  //   printf(" done\n");
  // }
  // add_cell_chunk(ugrid);
  // std::cout << "writing file" << std::endl;
  // write_VTUFile("amr_mesh.vtu", 0, ugrid);
  // std::cout << "Complete" << std::endl;

  // Surface mesh
  add_wall_faces(gridfile, datafile, izn, ugrid);

  write_VTUsurfaceFile("amr_surface.vtu", 1, ugrid);

  return EXIT_SUCCESS;
}
