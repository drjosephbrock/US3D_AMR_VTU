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
#include <unistd.h>

int main(int argc, char *argv[])
{
   int zone = 0;
   const char *filename = "amr_surface.vtu"; // Default string value
   int opt;
   bool plotSolution = false;
   bool plotSurface = false;
   bool plotVolume = false;
   while ((opt = getopt(argc, argv, "z:f:s:v")) != -1)
   {
      switch (opt)
      {
      case 'z':
         zone = std::atoi(optarg); // Convert the argument to an integer
         plotSurface = true;
         break;
      case 'f':
         filename = optarg; // Store the argument as a string
         break;
      case 's':
         plotSolution = true;
         break;
      case 'v':
         plotVolume = true;
         break;
      default:
         std::cout << opt << std::endl;
         // Handle invalid options or display usage information
         std::cerr << "Usage: " << argv[0] << " -z <zone_value> -f <filename> -s true (default to faults)" << std::endl;
         return 1;
      }
   }
   if (plotSolution)
   {
      std::cout << "Plotting mesh and solution" << std::endl;
   }
   else
   {
      std::cout << "Plotting mesh only" << std::endl;
   }

   const char *gridfile = "grid.h5";
   const char *datafile = "data.h5";
   vtkNew<vtkUnstructuredGrid> ugrid;
   vtkNew<vtkPolyData> polygonPolyData;

   int nn, nf, nc, nz;
   readGridArraySizes(gridfile, &nn, &nf, &nc, &nz);

   printf("  --  Nodes: %d\n", nn);
   printf("  --  Faces: %d\n", nf);
   printf("  --  Cells: %d\n", nc);
   printf("  --  Zones: %d\n", nz);

   // Volume mesh

   // Read Data and append to file
   if (plotVolume)
   {
      {
         printf("  --  reading points... ");
         std::vector<double> xcn(3*nn);
         readGridDoubleArray(gridfile, "/xcn", xcn);
         add_points(ugrid, xcn, nn);
         printf(" done\n");
      }
      
      add_cell_chunk(ugrid);
      
      std::vector<std::string> bvnames = readUS3DDataBVnames(datafile);
      std::vector<std::string> svnames = readUS3DDataSVnames(datafile);

      int nsv = getUS3DSolutionNSV(datafile);
      std::cout << "solution size: " << nsv << " x " << nc << std::endl;
      std::vector<double> solution = readUS3DSolutionFile(datafile, "/solution/run_1/interior");
      
      
      const char* varnames[nsv];
      
      // Initialize the const char* array with the C-style strings from the vector
      for (size_t i = 0; i < svnames.size(); ++i) {
         varnames[i] = svnames[i].c_str();
      }
      add_cell_data(ugrid, solution, varnames, nsv, nc);

      std::cout << "writing file" << std::endl;
      write_VTUFile("amr_mesh.vtu", 0, ugrid);
      std::cout << "Complete" << std::endl;
   }

   // Read Data and append to file
   if (plotSurface){
      // Surface mesh
      add_wall_faces(gridfile, datafile, zone, ugrid, plotSolution);
      write_VTUsurfaceFile(filename, 1, ugrid);
   }

   return EXIT_SUCCESS;
}
