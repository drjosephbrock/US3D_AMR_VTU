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

#include <vtkCellArray.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <iostream>
#include <hdf5.h>
#include "read_US3D_hdf5.H"

// Function to display a progress bar
void displayProgressBar2(int progress, int total)
{
  const int barWidth = 50;
  float fraction = static_cast<float>(progress) / total;
  int numBars = static_cast<int>(fraction * barWidth);

  std::cout << "[";

  for (int i = 0; i < barWidth; ++i)
  {
    if (i < numBars)
    {
      std::cout << "=";
    }
    else
    {
      std::cout << " ";
    }
  }

  std::cout << "] " << int(fraction * 100.0) << "%\r";
  std::cout.flush();
}

std::list<int> appendAndUniqify(std::list<int> &inputList, const std::list<int> &valuesToAppend)
{
  // Append values to the input list
  inputList.insert(inputList.end(), valuesToAppend.begin(), valuesToAppend.end());

  // Sort the list to bring duplicates together
  inputList.sort();

  // Remove duplicates
  inputList.unique();

  return inputList;
}

void write_VTUFile(const char *filename, bool writeAscii, const vtkNew<vtkUnstructuredGrid> &ugrid)
{

  // Here we write out the cube.
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetInputData(ugrid);
  writer->SetFileName(filename);
  if (writeAscii)
  {
    writer->SetDataModeToAscii();
  }
  writer->Update();
}

void write_VTUsurfaceFile(const char *filename, bool writeAscii, vtkNew<vtkUnstructuredGrid> &ugrid)
{

  // Here we write out the cube.
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetInputData(ugrid);
  writer->SetFileName(filename);
  if (writeAscii)
  {
    writer->SetDataModeToAscii();
  }
  writer->Update();
}

void add_wall_faces(const char *gridfile, const char *datafile, const int &izn, vtkNew<vtkUnstructuredGrid> &ugrid, const bool &plotSolution)
{

  vtkNew<vtkCellArray> polygons;
  int chunk_size = 100000;

  // Open gridfile
  hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
  // Open xcn (nodes), ifn (for interior cell data) and ifnpoly (for surface nodes)
  hid_t zdef_dataset = H5Dopen2(file_id, "/zones/zdefs", H5P_DEFAULT);
  hid_t xcn_dataset = H5Dopen2(file_id, "xcn", H5P_DEFAULT);
  hid_t ifn_dataset = H5Dopen2(file_id, "ifn", H5P_DEFAULT);
  hid_t ifnpoly_dataset = H5Dopen2(file_id, "ifnpoly", H5P_DEFAULT);
  // Get dataspaces for hyperslab
  hid_t zdef_dataspace = H5Dget_space(zdef_dataset);
  hid_t xcn_dataspace = H5Dget_space(xcn_dataset);
  hid_t ifn_dataspace = H5Dget_space(ifn_dataset);
  hid_t ifnpoly_dataspace = H5Dget_space(ifnpoly_dataset);

  hsize_t xcn_dims[2], ifn_dims[2], ifnpoly_dims[2], dataset_dims[2];

  H5Sget_simple_extent_dims(zdef_dataspace, dataset_dims, NULL);
  std::cout << "zdefs dimensions: " << dataset_dims[0] << " x " << dataset_dims[1] << std::endl;
  int zdef[dataset_dims[0]][dataset_dims[1]];

  H5Dread(zdef_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, zdef);
  std::cout << "Selected zone: " << izn << " has range: " << zdef[izn - 1][3] << " - " << zdef[izn - 1][4] << std::endl;

  int faceIDX;
  int nZoneFaces = zdef[izn - 1][4] - zdef[izn - 1][3] + 1;

  std::vector<int> selected_ifn_data(8);
  std::vector<int> selected_ifnpoly_data(9);
  std::vector<int> ifn_zone(9 * nZoneFaces);
  std::vector<int> FaceNodes;

  // Get surface nodes using ifnpoly (in case surface is refined)
  // Also add nodes from each ifn array into FaceNode vector....
  H5Sget_simple_extent_dims(ifnpoly_dataspace, ifn_dims, NULL);
  for (int i = 0; i < nZoneFaces; i++)
  {
    faceIDX = i + zdef[izn - 1][3] - 1;
    hsize_t start[2] = {(hsize_t)faceIDX, (hsize_t)0};
    hsize_t count[2] = {(hsize_t)1, (hsize_t)ifn_dims[1]};
    H5Sselect_hyperslab(ifnpoly_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t memspace = H5Screate_simple(2, count, NULL);
    hsize_t memstart[2] = {(hsize_t)0, (hsize_t)0};
    hsize_t memcount[2] = {(hsize_t)1, (hsize_t)ifn_dims[1]};
    hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, NULL, memcount, NULL);
    herr_t read_status = H5Dread(ifnpoly_dataset, H5T_NATIVE_INT, memspace, ifnpoly_dataspace, H5P_DEFAULT, selected_ifnpoly_data.data());
    for (int node = 0; node < 9; node++)
    {
      ifn_zone[9 * i + node] = selected_ifnpoly_data[node];
      if (selected_ifnpoly_data[node] > 0 && node > 0)
      {
        FaceNodes.push_back(selected_ifnpoly_data[node]);
      }
    }
  }

  // Sort and uniqify FaceNode vector.
  std::sort(FaceNodes.begin(), FaceNodes.end());
  auto last = std::unique(FaceNodes.begin(), FaceNodes.end());
  FaceNodes.erase(last, FaceNodes.end());
  std::cout << FaceNodes.size() << std::endl;

  H5Sget_simple_extent_dims(xcn_dataspace, xcn_dims, NULL);
  std::cout << "xcn dimensions: " << xcn_dims[0] << " x " << xcn_dims[1] << std::endl;
  std::vector<double> xcn_reduced(3 * FaceNodes.size());
  std::vector<double> selected_xcn_data(3);

  // Using reduced FaceNodes vector, build mapping from global to surface local numbering
  std::vector<int> nodeMap(xcn_dims[0]); // global to local mapping
  for (int i = 0; i < FaceNodes.size(); i++)
  {
    nodeMap[FaceNodes[i]] = i;

    hsize_t start[2] = {(hsize_t)(FaceNodes[i] - 1), (hsize_t)0};
    hsize_t count[2] = {(hsize_t)1, (hsize_t)3};
    H5Sselect_hyperslab(xcn_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t memspace = H5Screate_simple(2, count, NULL);
    hsize_t memstart[2] = {(hsize_t)0, (hsize_t)0};
    hsize_t memcount[2] = {(hsize_t)1, (hsize_t)3};
    hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, NULL, memcount, NULL);
    herr_t read_status = H5Dread(xcn_dataset, H5T_NATIVE_DOUBLE, memspace, xcn_dataspace, H5P_DEFAULT, selected_xcn_data.data());
    for (int direction = 0; direction < 3; direction++)
    {
      xcn_reduced[3 * i + direction] = selected_xcn_data[direction];
    }
  }

  add_points_poly(ugrid, xcn_reduced, FaceNodes.size());

  // Loop over faces of izn and save polygon data using nodeMap to convert to local numbering
  std::vector<int> localifn(9);
  for (int face = 0; face < nZoneFaces; face++)
  {
    localifn[0] = ifn_zone[9 * face];
    for (int j = 1; j < ifn_zone[9 * face] + 1; j++)
    { // ifn_zone[9*i] + s1 for first index being number of nodes for face
      localifn[j] = nodeMap[ifn_zone[9 * face + j]];
    }
    add_polygon(ugrid, face, localifn);
  }

  if (plotSolution)
  {
    hid_t datafile_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t bvars_dataset = H5Dopen2(datafile_id, "/solution/run_1/boundaries", H5P_DEFAULT);
    hid_t bvars_dataspace = H5Dget_space(bvars_dataset);
    H5Sget_simple_extent_dims(bvars_dataspace, dataset_dims, NULL);
    std::cout << "bvars dimensions: " << dataset_dims[0] << " x " << dataset_dims[1] << std::endl;

    int nv = dataset_dims[1];
    std::vector<double> faceSolution(nv * nZoneFaces);
    std::vector<int> selected_faceSolution(nv);
    for (int i = 0; i < nZoneFaces; i++)
    {
      faceIDX = i + zdef[izn - 1][3] - 1;
      hsize_t start[2] = {(hsize_t)faceIDX, (hsize_t)0};
      hsize_t count[2] = {(hsize_t)1, (hsize_t)nv};
      H5Sselect_hyperslab(bvars_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
      hid_t memspace = H5Screate_simple(2, count, NULL);
      hsize_t memstart[2] = {(hsize_t)0, (hsize_t)0};
      hsize_t memcount[2] = {(hsize_t)1, (hsize_t)nv};
      hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, NULL, memcount, NULL);
      herr_t read_status = H5Dread(bvars_dataset, H5T_NATIVE_INT, memspace, bvars_dataspace, H5P_DEFAULT, selected_faceSolution.data());
      for (int var = 0; var < nv; var++)
      {
        faceSolution[nv * i + var] = selected_faceSolution[var];
      }
    }
    
    int nbv = getUS3DSolutionNBV(datafile);
    std::vector<std::string> bvnames = readUS3DDataBVnames(datafile);
    const char* varnames[nbv];
    
    // Initialize the const char* array with the C-style strings from the vector
    for (size_t i = 0; i < bvnames.size(); ++i) {
        varnames[i] = bvnames[i].c_str();
    }
    add_cell_data(ugrid, faceSolution, varnames, nv, nZoneFaces);
    H5Fclose(datafile_id);
  }

  H5Sclose(xcn_dataspace);
  H5Sclose(ifn_dataspace);
  H5Sclose(ifnpoly_dataspace);
  H5Dclose(xcn_dataset);
  H5Dclose(ifn_dataset);
  H5Dclose(ifnpoly_dataset);
  H5Fclose(file_id);
}

void add_polygon(vtkNew<vtkUnstructuredGrid> &ugrid, const int &faceID, const std::vector<int> &ifn){

  // Create the polygon
  vtkNew<vtkPolygon> polygon;
  polygon->GetPointIds()->SetNumberOfIds(ifn[0]);
  for (int i = 0; i < ifn[0]; i++)
  {
    polygon->GetPointIds()->SetId(i, ifn[i + 1]);
  }
  ugrid->InsertNextCell(polygon->GetCellType(), polygon->GetPointIds());
}

void add_points_poly(vtkNew<vtkUnstructuredGrid> &ugrid, const std::vector<double> &xcn, const int &nn){

  vtkNew<vtkPoints> points;
  points->Allocate(nn);

  for (int i = 0; i < nn; i++)
  {
    // Calculate the indices once to avoid repeated calculations
    int index = 3 * i;
    double x = xcn[index];
    double y = xcn[index + 1];
    double z = xcn[index + 2];
    points->InsertNextPoint(x, y, z);
  }
  ugrid->SetPoints(points);
}

void add_cell_chunk(vtkNew<vtkUnstructuredGrid> &ugrid){

  hid_t file_id = H5Fopen("grid.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
  {
    fprintf(stderr, "Error opening HDF5 file\n");
    return;
  }

  hid_t ief_dataset = H5Dopen2(file_id, "iefpoly", H5P_DEFAULT);
  if (ief_dataset < 0)
  {
    fprintf(stderr, "Error opening ifn dataset\n");
    H5Fclose(file_id);
    return;
  }

  hid_t ifn_dataset = H5Dopen2(file_id, "ifnpoly", H5P_DEFAULT);
  if (ifn_dataset < 0)
  {
    fprintf(stderr, "Error opening xcn dataset\n");
    H5Dclose(ief_dataset);
    H5Fclose(file_id);
    return;
  }

  hid_t ief_dataspace = H5Dget_space(ief_dataset);
  if (ief_dataspace < 0)
  {
    fprintf(stderr, "Error getting ifn dataspace\n");
    H5Dclose(ief_dataset);
    H5Dclose(ifn_dataset);
    H5Fclose(file_id);
    return;
  }

  hsize_t ief_dims[2], ifn_dims[2];
  H5Sget_simple_extent_dims(ief_dataspace, ief_dims, NULL);
  int num_elements = (int)ief_dims[0];
  int num_faces = (int)ief_dims[1];

  int chunk_size = 100000;
  int faceIDX;

  std::vector<int> ifn_chunk(25 * 9); // max nFaces/cell
  std::vector<int> selected_ifn_data(9);

  hsize_t ifn_start[2] = {0, 0}; // Initialize the hyperslab start
  hsize_t ifn_count[2] = {1, 9}; // Adjust for your dataset dimensions (assuming 2D)

  // Loop to read and process data in chunks
  for (int start_idx = 0; start_idx < num_elements; start_idx += chunk_size) {
    // Calculate the end index for the current chunk
    int end_idx = start_idx + chunk_size;
    if (end_idx > num_elements) {
      end_idx = num_elements;
    }
    int slab = (end_idx - start_idx);

    std::vector<int> ief_chunk(slab * ief_dims[1]);

    hsize_t start[2] = {(hsize_t)start_idx, (hsize_t)0};
    hsize_t count[2] = {(hsize_t)slab, (hsize_t)ief_dims[1]};

    H5Sselect_hyperslab(ief_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t memspace = H5Screate_simple(2, count, NULL);

    if (memspace < 0) {
      std::cout << "Error creating memspace." << std::endl;
      H5Sclose(ief_dataspace);
      H5Dclose(ief_dataset);
      H5Fclose(file_id);
      return;
    }

    hsize_t memstart[2] = {(hsize_t)0, (hsize_t)0};
    hsize_t memcount[2] = {(hsize_t)slab, (hsize_t)ief_dims[1]};

    // hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, NULL, memcount, NULL);

    if (status < 0) {
      std::cout << "Error selecting hyperslab." << std::endl;
      H5Sclose(ief_dataspace);
      H5Dclose(ief_dataset);
      H5Fclose(file_id);
      return;
    }
    hssize_t numel = H5Sget_select_npoints(ief_dataspace);

    herr_t read_status = H5Dread(ief_dataset, H5T_NATIVE_INT, memspace, ief_dataspace, H5P_DEFAULT, ief_chunk.data());

    if (read_status < 0) {
      std::cout << "Error reading data. status=" << read_status << std::endl;

      H5Sclose(ief_dataspace);
      H5Dclose(ief_dataset);
      H5Fclose(file_id);
      return;
    }

    for (int cell = 0; cell < slab; cell++) {
      // Process the data chunk here
      int ncellFaces = ief_chunk[25 * cell];
      ifn_count[0] = 1;
      ifn_count[1] = 9;

      hid_t ifn_dataspace = H5Dget_space(ifn_dataset);
      H5Sget_simple_extent_dims(ifn_dataspace, ifn_dims, NULL);
      for (int cellface = 0; cellface < ncellFaces; cellface++) {
        faceIDX = (ief_chunk[25 * cell + cellface + 1] - 1);

        hsize_t start[2] = {(hsize_t)faceIDX, (hsize_t)0};
        hsize_t count[2] = {(hsize_t)1, (hsize_t)ifn_dims[1]};
        H5Sselect_hyperslab(ifn_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
        hid_t memspace = H5Screate_simple(2, count, NULL);
        hsize_t memstart[2] = {(hsize_t)0, (hsize_t)0};
        hsize_t memcount[2] = {(hsize_t)1, (hsize_t)ifn_dims[1]};
        hid_t status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, NULL, memcount, NULL);
        herr_t read_status = H5Dread(ifn_dataset, H5T_NATIVE_INT, memspace, ifn_dataspace, H5P_DEFAULT, selected_ifn_data.data());
        for (int node = 0; node < 9; node++)
        {
          ifn_chunk[9 * cellface + node] = selected_ifn_data[node];
        }
      }
      add_cell(ugrid, cell, ief_chunk, ifn_chunk);
      // H5Sclose(ifn_dataspace); // Close the dataspace
      // selected_ifn_data.clear(); // Clear and reuse the vector
    }
    displayProgressBar2(start_idx, num_elements);
  }

  // Close the HDF5 objects
  // delete[] ief_chunk;
  H5Dclose(ifn_dataset);
  H5Sclose(ief_dataspace);
  H5Dclose(ief_dataset);
  H5Fclose(file_id);
}

void add_cell(vtkNew<vtkUnstructuredGrid> &ugrid, const int &cellID, const std::vector<int> &ief, const std::vector<int> &ifn)
{

  int node, face;
  int nFaceNodes;
  int ncellFaces = ief[25 * cellID];
  std::list<int> CellNodes;
  vtkNew<vtkPoints> points;
  vtkNew<vtkIdList> faces;
  vtkIdType face8[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

  auto addFace = [&](const vtkIdType *face, int size)
  {
    faces->InsertNextId(size);
    for (int i = 0; i < size; ++i)
    {
      faces->InsertNextId(face[i]);
    }
  };

  for (int cellface = 0; cellface < ncellFaces; cellface++)
  {
    nFaceNodes = ifn[9 * cellface];
    for (int faceNode = 0; faceNode < nFaceNodes; faceNode++)
    {
      node = ifn[9 * cellface + faceNode + 1] - 1;
      CellNodes.push_back(node);
      face8[faceNode] = node;
    }
    addFace(face8, nFaceNodes);
  }

  std::list<int> UniqueCellNodes;
  appendAndUniqify(UniqueCellNodes, CellNodes);
  vtkIdType pointIds[UniqueCellNodes.size()];
  int index = 0;
  for (auto it = UniqueCellNodes.begin(); it != UniqueCellNodes.end(); ++it)
  {
    pointIds[index++] = *it;
  }
  ugrid->InsertNextCell(VTK_POLYHEDRON, UniqueCellNodes.size(), pointIds, ncellFaces, faces->GetPointer(0));
}

void add_points(vtkNew<vtkUnstructuredGrid> &ugrid, const std::vector<double> &xcn, const int &nn)
{
  vtkNew<vtkPoints> points;
  points->Allocate(nn);

  for (int i = 0; i < nn; i++)
  {
    // Calculate the indices once to avoid repeated calculations
    int index = 3 * i;
    double x = xcn[index];
    double y = xcn[index + 1];
    double z = xcn[index + 2];
    points->InsertNextPoint(x, y, z);
  }
  ugrid->SetPoints(points);
}

void add_cell_data(vtkNew<vtkUnstructuredGrid> &ugrid,
                   const std::vector<double> &solution,
                   const char *variables[], const int &nv, const int &nel)
{
  for (int var = 0; var < nv; var++)
  {
    vtkNew<vtkDoubleArray> scalars = vtkNew<vtkDoubleArray>();
    scalars->SetName(variables[var]);
    std::cout << "Writing '" << variables[var] << "' variable" << std::endl;
    for (int cell = 0; cell < nel; cell++)
    {
      double value = solution[nv * cell + var];
      scalars->InsertNextValue(value);
    }
    ugrid->GetCellData()->AddArray(scalars);
  }
}
