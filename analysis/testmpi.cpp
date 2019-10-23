/*
 * Analysis code for the Gray-Scott simulation.
 * Reads variable U and and extracts the iso-surface using VTK.
 * Writes the extracted iso-surface using ADIOS.
 *
 * Keichi Takahashi <keichi@is.naist.jp>
 *
 */
#include <putgetMeta/metaclient.h>
#include <iostream>
#include <sstream>

#include <adios2.h>

#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include "../common/timer.hpp"
#include <thread>
#include "../simulation/settings.h"

int main(int argc, char *argv[])
{
    std::cout << "debug 0.01" << std::endl;
    MPI_Init(&argc, &argv);

    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    std::cout << "debug 0.02" << std::endl;

    //const unsigned int color = 5;
    //MPI_Comm comm;
    //MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);
    std::cout << "debug 0.03" << std::endl;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    std::cout << "debug 0.04" << "proc " << procs << "rank " << rank <<  std::endl;
    
    /*
    int dims[3] = {0};
    MPI_Dims_create(procs, 3, dims);
    size_t npx = dims[0];
    size_t npy = dims[1];
    size_t npz = dims[2];

    std::cout << "debug 0.05" << std::endl;
    int coords[3] = {0};
    int periods[3] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    size_t px = coords[0];
    size_t py = coords[1];
    size_t pz = coords[2];

    std::cout << "debug 0.06" << std::endl;

    if (argc < 4)
    {
        if (rank == 0)
        {
            std::cerr << "Too few arguments" << std::endl;
            std::cout << "Usage: isosurface input output isovalue" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    std::cout << "debug 0" << std::endl;

    const std::string input_fname(argv[1]);
    const std::string output_fname(argv[2]);
    const double isovalue = std::stod(argv[3]);

    std::cout << "debug 0.1" << std::endl;

    */

   sleep(5);

   MPI_Finalize();

    return 0;
}
