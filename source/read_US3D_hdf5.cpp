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

void readGridArraySizes(const char* gridfile, int* nn, int* nf, int* nc, int* nz){
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
        return;
    }

    // Read attributes
    // HDF5Log("Reading nn");
    status = H5Aread(nn_id, H5T_NATIVE_INT, nn);
    if (status < 0) {
        HDF5Log("Error reading nn attribute.");
    }

    // HDF5Log("Reading nc");
    status = H5Aread(nc_id, H5T_NATIVE_INT, nc);
    if (status < 0) {
        HDF5Log("Error reading nc attribute.");
    }

    // HDF5Log("Reading nf");
    status = H5Aread(nf_id, H5T_NATIVE_INT, nf);
    if (status < 0) {
        HDF5Log("Error readingnf attribute.");
    }

    // HDF5Log("Reading nz");
    status = H5Aread(nz_id, H5T_NATIVE_INT, nz);
    if (status < 0) {
        HDF5Log("Error reading nz attribute.");
    }

    status = H5Aclose(nn_id);
    status = H5Aclose(nc_id);
    status = H5Aclose(nf_id);
    status = H5Aclose(nz_id);

    status = H5Gclose(group_id);
    status = H5Fclose(grid_id);

    return;
}

void readGridDoubleArray(const char* gridfile, const char* dataName, std::vector<double>& data){

    herr_t status;

    hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    data.resize(dims[0] * dims[1]);
    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return;

}

void readGridIntegerArray(const char* gridfile, const char* dataName, std::vector<int>& data){
    
    herr_t status;

    hid_t file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    data.resize(dims[0] * dims[1]);
    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return;

}

std::vector<double> readUS3DSolutionFile(const char* datafile, const char* dataName, int numRows, int numCols){

    herr_t status;

    hid_t file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    numCols = (int)dims[0];
    numRows = (int)dims[1];
    std::vector<double> data(numRows * numCols);

    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return data;

}
