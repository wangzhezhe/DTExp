

#include <putEvent/eventclient.h>
#include <putgetMeta/metaclient.h>
#include <utils/ThreadPool.h>
#include <unistd.h>
#include <queue>

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

bool epsilon(double d) { return (d < 1.0e-20); }
bool epsilon(float d) { return (d < 1.0e-20); }

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

    if (argc != 2)
    {
        std::cout << " <infileName>\n";
        MPI_Finalize();
        return 0;
    }

    //init the metaclient
    MetaClient metaclient = getMetaClient();
    std::string keyDataOk = "INDICATOR1";
    std::string keySimFinish = "SIMFINISH";
    EventClient eventclient = getEventClient();

    std::string in_filename;
    size_t nbins = 100;
    in_filename = argv[1];

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

    std::ostringstream log_fname;
    log_fname << "pdf_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    log << "step\tcompute_pdf" << std::endl;
#endif

    while (true)
    {

        // Begin step
        adios2::StepStatus read_status =
            reader.BeginStep(adios2::StepMode::Read, 10.0f);
        if (read_status == adios2::StepStatus::NotReady)
        {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK)
        {
            break;
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
        if (shouldIWrite)
        {
            reader.Get<int>(var_step_in, &simStep);
        }

        // End adios2 step
        reader.EndStep();

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
        //if (simStep  == 1)
        //if(simStep>=0)
        //if (simStep % 2 == 0)
        //if (simStep >= 0)
        //if (simStep == 1 || simStep == 10)
        //if (simStep >= 0)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        //if(simStep%10==1)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        //if(simStep%5==1)
        //if(simStep%2==1)
        //if(simStep%5==1)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        //if(true)
        //if(simStep%5==1 || simStep%5==2 || simStep%5==3 || simStep%5==4)
        if(true)
        {
            indicator = true;
        }
        if (rank == 0)
        {
            if (indicator)
            {
                //set the metadata to the metadata server
                //assume the consumer know how to generate the variable

                std::string metainfots = std::to_string(simStep - 1);
                //vector<string> eventList, string source, string metadata, string matchType
                std::vector<std::string> pushList;
                pushList.push_back("INTERESTINGTOPIC1");
                std::string reply = eventclient.Publish(pushList, "testid", metainfots, "NAME");
                std::cout << "event publish recieve: " << reply << " for ts " << simStep << std::endl;
            }
        }

        ++stepAnalysis;
    }

    // cleanup
    reader.Close();

    if (rank == 0)
    {
        //program finish, putmeta
        std::string reply = metaclient.Putmetaspace(keySimFinish, "OK");
        std::cout << "Program finish request recieve: " << reply << std::endl;
    }
    MPI_Finalize();

    return 0;
}