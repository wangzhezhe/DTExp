#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <adios2.h>
#include <mpi.h>

#include "../common/timer.hpp"
#include "gray-scott.h"
#include "writer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <thread>
#include <putgetMeta/metaclient.h>

bool epsilon(double d) { return (d < 1.0e-20); }
bool epsilon(float d) { return (d < 1.0e-20); }

adios2::Variable<double> var_u_pdf, var_v_pdf;
adios2::Variable<double> var_u_bins, var_v_bins;
adios2::Variable<int> var_step_out;
adios2::Variable<double> var_u_out, var_v_out;

/*
 * Function to compute the PDF of a 2D slice
 */
template <class T>
void compute_pdf(const std::vector<T> data,
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

bool DoCheck(int rank, int procs, unsigned int step, GrayScott &sim)
{
    unsigned int stepAnalysis = step;
    unsigned int simStep = stepAnalysis;

    std::vector<double> pdf_u;
    std::vector<double> bins_u;

    size_t nbins = 100;

    std::vector<double> u = sim.u_noghost();

    std::vector<std::size_t> shape(3, 0);
    shape[0] = sim.size_x;
    shape[1] = sim.size_y;
    shape[2] = sim.size_z;

    if (rank == 0)
    {
        std::cout << "PDF Analysis step " << stepAnalysis
                  << " processing sim output step " << simStep
                  << " sim compute step " << simStep << std::endl;
    }

    // HDF5 engine does not provide min/max. Let's calculate it
    //        if (reader_io.EngineType() == "HDF5")

    auto mmu = std::minmax_element(u.begin(), u.end());
    std::pair<double, double> minmax_u = std::make_pair(*mmu.first, *mmu.second);

    int comm_size = procs;
    size_t count1 = shape[0] / comm_size;
    size_t start1 = count1 * rank;

    // Compute PDF
    compute_pdf(u, shape, start1, count1, nbins, minmax_u.first,
                minmax_u.second, pdf_u, bins_u);

    //in real case, this results should be aggregated from different rank, than start checking at one place
    //if (step >= 0)
    //if (step % 2 == 0)
    //if (step == 1)
    //if (simStep == 1 || simStep == 10)
    //
    //if(step%5==1 || step%5==2 || step%5==3 || simStep%5==4)
    //if (step >= 0)
    //if(step%2==1)
    if(simStep%5==1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void print_io_settings(const adios2::IO &io)
{
    std::cout << "Simulation writes data using engine type:              "
              << io.EngineType() << std::endl;
}

void print_settings(const Settings &s)
{
    std::cout << "grid:             " << s.L << "x" << s.L << "x" << s.L
              << std::endl;
    std::cout << "steps:            " << s.steps << std::endl;
    std::cout << "plotgap:          " << s.plotgap << std::endl;
    std::cout << "F:                " << s.F << std::endl;
    std::cout << "k:                " << s.k << std::endl;
    std::cout << "dt:               " << s.dt << std::endl;
    std::cout << "Du:               " << s.Du << std::endl;
    std::cout << "Dv:               " << s.Dv << std::endl;
    std::cout << "noise:            " << s.noise << std::endl;
    std::cout << "output:           " << s.output << std::endl;
    std::cout << "adios_config:     " << s.adios_config << std::endl;
}

void print_simulator_settings(const GrayScott &s)
{
    std::cout << "process layout:   " << s.npx << "x" << s.npy << "x" << s.npz
              << std::endl;
    std::cout << "local grid size:  " << s.size_x << "x" << s.size_y << "x"
              << s.size_z << std::endl;
}

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 1;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    if (argc < 2)
    {
        if (rank == 0)
        {
            std::cerr << "Too few arguments" << std::endl;
            std::cerr << "Usage: gray-scott settings.json" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (rank == 0)
    {
        //send request to timer that wf start
        MetaClient metaclient = getMetaClient();
        string reply = metaclient.Recordtimestart("WFTIMER");
        std::cout << "Timer received: " << reply << std::endl;
    }

    Settings settings = Settings::from_json(argv[1]);

    GrayScott sim(settings, comm);
    sim.init();

    adios2::ADIOS adios(settings.adios_config, comm, adios2::DebugON);

    adios2::IO io_main = adios.DeclareIO("SimulationOutput");

    Writer writer_main(settings, sim, io_main);

    writer_main.open(settings.output);

#ifdef ENABLE_TIMERS
    Timer timer_compute;
    Timer timer_check;
    Timer timer_writer;

    std::ostringstream log_fname;
    log_fname << "gray_scott_pdf_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
#endif

    if (rank == 0)
    {
        print_settings(settings);
        print_simulator_settings(sim);
        std::cout << "========================================" << std::endl;

        std::cout << "PDF analysis reads from Simulation use function call:  " << std::endl;
    }

    for (int i = 0; i < settings.steps;)
    {

#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_compute.start();
#endif
        for (int j = 0; j < settings.plotgap; j++)
        {
            sim.iterate();
            //this time is used to construct different simulating-analysis patterns
            std::this_thread::sleep_for(std::chrono::milliseconds(400));
            i++;
        }

        if (i == settings.plotgap)
        {
            std::cout << "engine type for sim output is " << io_main.EngineType() << std::endl;
        }

#ifdef ENABLE_TIMERS
        double timer_compute_value = timer_compute.stop();
        MPI_Barrier(comm);
        log << "timer compute sim: " << timer_compute_value << std::endl;
        timer_check.start();
#endif
        bool indicator = DoCheck(rank, procs, i, sim);

#ifdef ENABLE_TIMERS
        double timer_check_value = timer_check.stop();
        MPI_Barrier(comm);
        log << "timer check pdf: " << timer_check_value << std::endl;

        timer_check.start();
#endif
        if (rank == 0)
        {
            std::cout << "indicator is " << indicator << " for do check at ts: " << i << std::endl;
        }

        if (indicator)
        {
#ifdef ENABLE_TIMERS

            timer_writer.start();
#endif
            MPI_Barrier(comm);
            writer_main.write(i, sim);

#ifdef ENABLE_TIMERS
            double timer_writer_value = timer_writer.stop();
            MPI_Barrier(comm);
            log << "timer writer: " << timer_writer_value << std::endl;

#endif
        }
    }

    writer_main.close();


    if (rank == 0)
    {
        MetaClient metaclient = getMetaClient();
        string reply = metaclient.Recordtimetick("WFTIMER");
        std::cout << "Timer received for simwithpdf finish: " << reply << std::endl;
    }
    MPI_Finalize();
}
