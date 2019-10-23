/*
 * Analysis code for the Gray-Scott application.
 * Reads variable U and V, and computes the PDF for each 2D slices of U and V.
 * Writes the computed PDFs using ADIOS.
 *
 * Norbert Podhorszki, pnorbert@ornl.gov
 *
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

#include "adios2.h"
#include "../common/timer.hpp"
#include "../simulation/settings.h"
#include <unistd.h>

#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <putgetMeta/metaclient.h>

bool epsilon(double d) { return (d < 1.0e-20); }
bool epsilon(float d) { return (d < 1.0e-20); }

void write_vtk(const std::string &fname,
               const vtkSmartPointer<vtkPolyData> polyData)
{
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fname.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

vtkSmartPointer<vtkPolyData>
compute_isosurface(const adios2::Variable<double> &varField,
                   const std::vector<double> &field, double isovalue)
{
    // Convert field values to vtkImageData
    auto importer = vtkSmartPointer<vtkImageImport>::New();
    importer->SetDataSpacing(1, 1, 1);
    importer->SetDataOrigin(varField.Start()[2], varField.Start()[1],
                            varField.Start()[0]);
    importer->SetWholeExtent(0, varField.Count()[2] - 1, 0,
                             varField.Count()[1] - 1, 0,
                             varField.Count()[0] - 1);
    importer->SetDataExtentToWholeExtent();
    importer->SetDataScalarTypeToDouble();
    importer->SetNumberOfScalarComponents(1);
    importer->SetImportVoidPointer(const_cast<double *>(field.data()));

    // Run the marching cubes algorithm
    auto mcubes = vtkSmartPointer<vtkMarchingCubes>::New();
    mcubes->SetInputConnection(importer->GetOutputPort());
    mcubes->ComputeNormalsOn();
    mcubes->SetValue(0, isovalue);
    mcubes->Update();

    // Return the isosurface as vtkPolyData
    return mcubes->GetOutput();
}

/*
 * Function to compute the PDF of a 2D slice
 */
template <class T>
void compute_pdf(const std::vector<T> &data,
                 const std::vector<std::size_t> &shape, const size_t start,
                 const size_t count, const size_t nbins, const T min,
                 const T max, std::vector<T> &pdf, std::vector<T> &bins)
{
    if (shape.size() != 3)
        throw std::invalid_argument("ERROR: shape is expected to be 3D\n");

    size_t slice_size = shape[1] * shape[2];
    pdf.resize(count * nbins);
    bins.resize(nbins);

    size_t start_data = 0;
    size_t start_pdf = 0;

    T binWidth = (max - min) / nbins;
    for (auto i = 0; i < nbins; ++i)
    {
        bins[i] = min + (i * binWidth);
    }

    if (nbins == 1)
    {
        // special case: only one bin
        for (auto i = 0; i < count; ++i)
        {
            pdf[i] = slice_size;
        }
        return;
    }

    if (epsilon(max - min) || epsilon(binWidth))
    {
        // special case: constant array
        for (auto i = 0; i < count; ++i)
        {
            pdf[i * nbins + (nbins / 2)] = slice_size;
        }
        return;
    }

    for (auto i = 0; i < count; ++i)
    {
        // Calculate a PDF for 'nbins' bins for values between 'min' and 'max'
        // from data[ start_data .. start_data+slice_size-1 ]
        // into pdf[ start_pdf .. start_pdf+nbins-1 ]
        for (auto j = 0; j < slice_size; ++j)
        {
            if (data[start_data + j] > max || data[start_data + j] < min)
            {
                std::cout << " data[" << start * slice_size + start_data + j
                          << "] = " << data[start_data + j]
                          << " is out of [min,max] = [" << min << "," << max
                          << "]" << std::endl;
            }
            size_t bin = static_cast<size_t>(
                std::floor((data[start_data + j] - min) / binWidth));
            if (bin == nbins)
            {
                bin = nbins - 1;
            }
            ++pdf[start_pdf + bin];
        }
        start_pdf += nbins;
        start_data += slice_size;
    }
    return;
}

/*
 * Print info to the user on how to invoke the application
 */
void printUsage()
{
    std::cout
        << "Usage: pdf_calc input output [N] [output_inputdata]\n"
        << "  input:   Name of the input file handle for reading data\n"
        << "  output:  Name of the output file to which data must be written\n"
        << "  N:       Number of bins for the PDF calculation, default = 1000\n"
        << "  output_inputdata: YES will write the original variables besides "
           "the analysis results\n\n";
}

/*
 * MAIN
 */
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, comm_size, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 2;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (argc != 3)
    {
        std::cout << " <infileName> <isovalue>\n";
        if (rank == 0)
            printUsage();
        MPI_Finalize();
        return 0;
    }

    std::string in_filename;
    size_t nbins = 100;
    in_filename = argv[1];
    int isovalue = static_cast<size_t>(std::stoi(argv[2]));

    //get grid size from setting
    Settings settings = Settings::from_json("settings.json");
    size_t grid_size = settings.L;
    std::cout << "------use grid size " << grid_size << std::endl;

    std::size_t u_global_size,
        v_global_size;
    std::size_t u_local_size, v_local_size;

    std::vector<std::size_t> shape;

    std::vector<double> u;
    std::vector<double> v;
    int simStep = -5;

    std::vector<double> pdf_u;
    std::vector<double> pdf_v;
    std::vector<double> bins_u;
    std::vector<double> bins_v;

    // adios2 variable declarations
    adios2::Variable<double> var_u_in, var_v_in;
    adios2::Variable<int> var_step_in;
    adios2::Variable<double> var_u_pdf, var_v_pdf;
    adios2::Variable<double> var_u_bins, var_v_bins;
    adios2::Variable<int> var_step_out;
    adios2::Variable<double> var_u_out, var_v_out;

    // adios2 io object and engine init
    adios2::ADIOS ad("adios2.xml", comm, adios2::DebugON);

    // IO objects for reading and writing
    adios2::IO reader_io = ad.DeclareIO("SimulationOutput");
    if (!rank)
    {
        std::cout << "PDF analysis reads from Simulation using engine type:  "
                  << reader_io.EngineType() << std::endl;
    }

    // Engines for reading and writing
    adios2::Engine reader =
        reader_io.Open(in_filename, adios2::Mode::Read, comm);

    bool shouldIWrite = (!rank || reader_io.EngineType() == "HDF5");

    // read data per timestep
    int stepAnalysis = 0;

#ifdef ENABLE_TIMERS
    Timer timer_compute;
    Timer timer_read;

    std::ostringstream log_fname;
    log_fname << "pdf_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    log << "step\tread_pdf" << std::endl;
#endif

    while (true)
    {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_read.start();
#endif

        // Begin step
        adios2::StepStatus read_status =
            reader.BeginStep(adios2::StepMode::Read, 10.0f);

        if (read_status == adios2::StepStatus::OtherError)
        {
            std::cout << "---in pdfwithana adios2 status is unknown---" << read_status << std::endl;
            break;
        }
        if (read_status != adios2::StepStatus::OK)
        {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout << "get status for pdf " << read_status << std::endl;
            if (read_status == adios2::StepStatus::EndOfStream)
            {
                break;
            }
            continue;
        }

        int stepSimOut = reader.CurrentStep();

        // Inquire variable and set the selection at the first step only
        // This assumes that the variable dimensions do not change across
        // timesteps

        // Inquire variable
        var_u_in = reader_io.InquireVariable<double>("U");
        var_step_in = reader_io.InquireVariable<int>("step");

        std::pair<double, double> minmax_u = var_u_in.MinMax();

        shape = var_u_in.Shape();

        // Calculate global and local sizes of U and V
        u_global_size = shape[0] * shape[1] * shape[2];
        u_local_size = u_global_size / comm_size;

        size_t count1 = shape[0] / comm_size;
        size_t start1 = count1 * rank;
        if (rank == comm_size - 1)
        {
            // last process need to read all the rest of slices
            count1 = shape[0] - count1 * (comm_size - 1);
        }

        /*std::cout << "  rank " << rank << " slice start={" <<  start1
            << ",0,0} count={" << count1  << "," << shape[1] << "," << shape[2]
            << "}" << std::endl;*/

        // Set selection
        var_u_in.SetSelection(adios2::Box<adios2::Dims>(
            {start1, 0, 0}, {count1, shape[1], shape[2]}));

        // Read adios2 data
        reader.Get<double>(var_u_in, u);
        reader.Get<int>(var_step_in, simStep);

        // End adios2 step
        reader.EndStep();
#ifdef ENABLE_TIMERS
        double timer_read_t = timer_read.stop();
        log << simStep << "\t" << timer_read_t << std::endl;
        MPI_Barrier(comm);
#endif

        std::cout << "read_time for " << simStep << " " << timer_read_t << std::endl;

        if (!rank)
        {
            std::cout << "PDF Analysis step " << stepAnalysis
                      << " processing sim output step " << stepSimOut
                      << " sim compute step " << simStep << std::endl;
        }

        // HDF5 engine does not provide min/max. Let's calculate it
        //        if (reader_io.EngineType() == "HDF5")
        {
            auto mmu = std::minmax_element(u.begin(), u.end());
            minmax_u = std::make_pair(*mmu.first, *mmu.second);
        }

        // Compute PDF
        std::vector<double> pdf_u;
        std::vector<double> bins_u;

        compute_pdf(u, shape, start1, count1, nbins, minmax_u.first,
                    minmax_u.second, pdf_u, bins_u);

        //start analytics if the checking indicator is ok
        //in real case, this results should be aggregated from different rank, than start checking at one place
        //this value shoule be decided based on the output of the compute pdf in real case
        bool indicator = false;
        //
        //if (simStep == 1)
        //if(simStep >= 0)
        //if (simStep % 2 == 0)
        //if (simStep == 1 || simStep == 10)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        //if(simStep%10==1)
        //if(simStep%2==1)
        if(simStep%5==1)
        {
            indicator = true;
        }

        if (indicator)
        {
            //caculate pdf
            //todo add the checking operation for the pdf
            auto polyData = compute_isosurface(var_u_in, u, isovalue);

            //write_adios(writer, polyData, varPoint, varCell, varNormal, varOutStep,
            //            step, comm);
            //replase by sleep for experiment to avoid the impact of IO
            Settings settings = Settings::from_json("./settings.json");
            //sim<ana+c
            std::this_thread::sleep_for(std::chrono::milliseconds(10 * settings.L));

            //sim>ana+c
            //std::this_thread::sleep_for(std::chrono::milliseconds(3*settings.L));
            /*
            std::string dir = "./vtkdata";

            char countstr[50];
            sprintf(countstr, "%04d", simStep);

            std::string fname = dir + "/vtkiso_" + std::string(countstr) + ".vtk";
            //the format here is the vtk
            write_vtk(fname, polyData);
            */
            std::cout << "ok for ts " << simStep << std::endl;

            //sleep adjusted time
            /*
            usleep(1000 * 5 * grid_size);
            std::cout << "adjusted time " << 1000 * 5 * grid_size << std::endl;
            */
        }

        ++stepAnalysis;
    }

    // cleanup
    reader.Close();
    if (rank == 0)
    {
        //tick finish
        MetaClient metaclient = getMetaClient();
        string reply = metaclient.Recordtimetick("WFTIMER");
    }
    MPI_Finalize();

    return 0;
}
