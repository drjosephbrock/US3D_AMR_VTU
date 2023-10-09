#include <iostream>
#include <cstdio>
#include <fstream>
#include <hdf5.h>
#include <vector>
#include <math.h>
#include <tuple>
#include "read_US3D_hdf5.H"

void HDF5Log(const std::string& message){
    std::cout << message << std::endl;
}

std::tuple<int, int, int, int>  readGridArraySizes(const char* gridfile){
    int nNodes, nFaces, nCells, nZones;
    hid_t grid_id, group_id;
    hid_t nn_id, nc_id, nf_id, nz_id;
    herr_t status;

    grid_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(grid_id, "/info/grid", H5P_DEFAULT);
    
    nn_id = H5Aopen(group_id, "nn", H5P_DEFAULT);
    nc_id = H5Aopen(group_id, "nc", H5P_DEFAULT);
    nf_id = H5Aopen(group_id, "nf", H5P_DEFAULT);

    status = H5Gclose(group_id);
    group_id = H5Gopen(grid_id, "/zones", H5P_DEFAULT);
    nz_id = H5Aopen(group_id, "nz", H5P_DEFAULT);
    
    // Check for errors in opening datasets
    if (nn_id < 0 || nc_id < 0 || nf_id < 0 || nz_id < 0) {
        HDF5Log("Failed to open one or more datasets.");
        H5Fclose(grid_id); // Close file before returning
        return std::make_tuple(-1, -1, -1, -1);
    }

    // Read attributes
    // HDF5Log("Reading nn");
    status = H5Aread(nn_id, H5T_NATIVE_INT, &nNodes);
    if (status < 0) {
        HDF5Log("Error reading nn attribute.");
    }

    // HDF5Log("Reading nc");
    status = H5Aread(nc_id, H5T_NATIVE_INT, &nCells);
    if (status < 0) {
        HDF5Log("Error reading nc attribute.");
    }

    // HDF5Log("Reading nf");
    status = H5Aread(nf_id, H5T_NATIVE_INT, &nFaces);
    if (status < 0) {
        HDF5Log("Error readingnf attribute.");
    }

    // HDF5Log("Reading nz");
    status = H5Aread(nz_id, H5T_NATIVE_INT, &nZones);
    if (status < 0) {
        HDF5Log("Error reading nz attribute.");
    }

    status = H5Aclose(nn_id);
    status = H5Aclose(nc_id);
    status = H5Aclose(nf_id);
    status = H5Aclose(nz_id);

    status = H5Gclose(group_id);
    status = H5Fclose(grid_id);

    return std::make_tuple(nNodes, nFaces, nCells, nZones);
}

std::vector<double> readGridDoubleArray(const char* gridfile, const char* dataName){
    
    herr_t status;

    hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    int numRows = dims[0];
    int numCols = dims[1];
    std::vector<double> data(numRows * numCols);

    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return data;

}

std::vector<int> readGridIntegerArray(const char* gridfile, const char* dataName){
    
    herr_t status;

    hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    int numRows = dims[0];
    int numCols = dims[1];
    std::vector<int> data(numRows * numCols);

    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return data;

}

std::tuple< std::vector<double>, std::vector<int>, std::vector<int> > readGridArrays(const char* gridfile) {

    // int nn, nf, nc, nz;

    auto [nn, nf, nc, nz] = readGridArraySizes(gridfile);
    printf("  --  Nodes: %d\n", nn);
    printf("  --  Faces: %d\n", nf);
    printf("  --  Cells: %d\n", nc);
    printf("  --  Zones: %d\n", nz);
    std::vector<double> xcn(nn * 3);
    std::vector<int> ifn(nf * 8);
    std::vector<int> ief(nc*25);
    // hid_t grid_id, dataset_id, group_id, attr_id;
    hid_t grid_id, xcn_id, ief_id, ifn_id;
    herr_t status;

    grid_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (grid_id < 0) {
        HDF5Log("Failed to open grid file.");
        return std::make_tuple(xcn, ifn, ief);
    }

    // Open datasets
    xcn_id = H5Dopen(grid_id, "/xcn", H5P_DEFAULT);
    ief_id = H5Dopen(grid_id, "/iefpoly", H5P_DEFAULT);
    ifn_id = H5Dopen(grid_id, "/ifnpoly", H5P_DEFAULT);
    
    // Check for errors in opening datasets
    if (xcn_id < 0 || ief_id < 0 || ifn_id < 0) {
        HDF5Log("Failed to open one or more datasets.");
        H5Fclose(grid_id); // Close file before returning
        return std::make_tuple(xcn, ifn, ief);
    }

    // Read datasets
    // HDF5Log("Reading xcn");
    status = H5Dread(xcn_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xcn.data());
    if (status < 0) {
        HDF5Log("Error reading xcn dataset.");
    }

    // HDF5Log("Reading ief");
    status = H5Dread(ief_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ief.data());
    if (status < 0) {
        HDF5Log("Error reading ief dataset.");
    }

    // HDF5Log("Reading ifn");
    status = H5Dread(ifn_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ifn.data());
    if (status < 0) {
        HDF5Log("Error reading ifn dataset.");
    }
    HDF5Log("Complete");

    HDF5Log("Closing datasets");
    // Close datasets
    H5Dclose(xcn_id);
    H5Dclose(ief_id);
    H5Dclose(ifn_id);

    HDF5Log("Closing file");
    // Close file
    H5Fclose(grid_id);

    HDF5Log("returning tuple....");
    return std::make_tuple(xcn, ifn, ief);
}