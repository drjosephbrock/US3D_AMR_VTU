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

// Function to trim trailing white space from a string
std::string trim(const std::string& str) {
    // size_t end = str.find_last_not_of(" \t\n\r\f\v");
    size_t end = str.find_last_not_of("\n\r\f\v");
    if (end == std::string::npos) {
        // String contains only white space
        return "";
    }
    return str.substr(0, end + 1);
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

std::vector<double> readUS3DSolutionFile(const char* datafile, const char* dataName){

    herr_t status;

    hid_t file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, dataName, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);
    
    // Get the dimensions of the data
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    std::vector<double> data(dims[0] * dims[1]);

    // Read the data from the dataset
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());


    // Close the dataset, dataspace, and file
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return data;

}

void readUS3DSolutionNSV(const char* datafile, int* nsv){

    hid_t file_id, group_id;
    hid_t nsv_id;
    herr_t status;

    file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/info/solver", H5P_DEFAULT);
    nsv_id = H5Aopen(group_id, "nsv", H5P_DEFAULT);
    status = H5Aread(nsv_id, H5T_NATIVE_INT, nsv);
    status = H5Aclose(nsv_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);

}

int getUS3DSolutionNSV(const char* datafile){

    hid_t file_id, group_id;
    hid_t nsv_id;
    herr_t status;
    int nsv;

    file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/info/solver", H5P_DEFAULT);
    nsv_id = H5Aopen(group_id, "nsv", H5P_DEFAULT);
    status = H5Aread(nsv_id, H5T_NATIVE_INT, &nsv);
    status = H5Aclose(nsv_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
    return nsv;
}

int getUS3DSolutionNBV(const char* datafile){

    hid_t file_id, group_id;
    hid_t nsv_id;
    herr_t status;
    int nsv;

    file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id = H5Gopen(file_id, "/info/solver", H5P_DEFAULT);
    nsv_id = H5Aopen(group_id, "nbv", H5P_DEFAULT);
    status = H5Aread(nsv_id, H5T_NATIVE_INT, &nsv);
    status = H5Aclose(nsv_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
    return nsv;
}

std::vector<std::string> readUS3DDataSVnames(const char* datafile){

    std::vector<std::string> svnames;

    // Open the HDF5 file
    hid_t file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Unable to open file 'data.h5'\n");
        return svnames;
    }

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, "/info/solver/svnames", H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Unable to open dataset '/info/solver/svnames'\n");
        H5Fclose(file_id);
        return svnames;
    }

    // Get the datatype and dataspace of the dataset
    hid_t datatype = H5Dget_type(dataset_id);
    hid_t dataspace = H5Dget_space(dataset_id);

    // Get the number of elements in the dataset
    hsize_t num_elements;
    H5Sget_simple_extent_dims(dataspace, &num_elements, NULL);

    // Read the data from the dataset
    char data[num_elements][20];
    char dest[256];     // Destination string
    svnames.resize(num_elements);

    H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    // Print the values
    for (int i = 0; i < num_elements; i++) {
        std::strncpy(dest, data[i], 20);
        dest[20] = 0; // null terminate destination
        svnames[i] = dest;
    }

    H5Tclose(datatype);
    H5Sclose(dataspace);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return svnames;
}

std::vector<std::string> readUS3DDataBVnames(const char* datafile){

    std::vector<std::string> bvnames;

    // Open the HDF5 file
    hid_t file_id = H5Fopen(datafile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Unable to open file 'data.h5'\n");
        return bvnames;
    }

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, "/info/solver/bvnames", H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Unable to open dataset '/info/solver/bvnames'\n");
        H5Fclose(file_id);
        return bvnames;
    }

    // Get the datatype and dataspace of the dataset
    hid_t datatype = H5Dget_type(dataset_id);
    hid_t dataspace = H5Dget_space(dataset_id);

    // Get the number of elements in the dataset
    hsize_t num_elements;
    H5Sget_simple_extent_dims(dataspace, &num_elements, NULL);
    std::cout << num_elements << std::endl;

    // Read the data from the dataset
    char data[num_elements][20];
    char dest[256];     // Destination string
    bvnames.resize(num_elements);

    H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    // Print the values
    for (int i = 0; i < num_elements; i++) {
        std::strncpy(dest, data[i], 20);
        dest[20] = 0; // null terminate destination
        bvnames[i] = dest;
    }

    for (int i = 0; i < num_elements; i++){
        std::cout << "bvnames[" << i << "]:" << trim(bvnames[i]) << std::endl;
    }

    H5Tclose(datatype);
    H5Sclose(dataspace);
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return bvnames;
}