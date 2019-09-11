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

bool epsilon(double d) { return (d < 1.0e-20); }
bool epsilon(float d) { return (d < 1.0e-20); }

// MPI info
int rank, procs, wrank;

adios2::Variable<double> var_u_pdf, var_v_pdf;
adios2::Variable<double> var_u_bins, var_v_bins;
adios2::Variable<int> var_step_out;
adios2::Variable<double> var_u_out, var_v_out;

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

void DoCheck(adios2::IO &pdfwriter_io, adios2::Engine &pdfwriter, int rank, unsigned int step, bool iffirstStep, GrayScott &sim)
{
    unsigned int stepAnalysis = step;
    unsigned int simStep = stepAnalysis;

    std::vector<double> pdf_u;
    std::vector<double> pdf_v;
    std::vector<double> bins_u;
    std::vector<double> bins_v;

    size_t nbins = 100;

    std::vector<double> u = sim.u_noghost();
    std::vector<double> v = sim.v_noghost();

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
    auto mmv = std::minmax_element(v.begin(), v.end());
    std::pair<double, double> minmax_v = std::make_pair(*mmv.first, *mmv.second);

    int comm_size = procs;
    size_t count1 = shape[0] / comm_size;
    size_t start1 = count1 * rank;

    // Compute PDF
    compute_pdf(u, shape, start1, count1, nbins, minmax_u.first,
                minmax_u.second, pdf_u, bins_u);

    compute_pdf(v, shape, start1, count1, nbins, minmax_v.first,
                minmax_v.second, pdf_v, bins_v);

    // write U, V, and their norms out
    if (iffirstStep)
    {
        var_u_pdf = pdfwriter_io.DefineVariable<double>(
            "U/pdf", {shape[0], nbins}, {start1, 0}, {count1, nbins});
        var_v_pdf = pdfwriter_io.DefineVariable<double>(
            "V/pdf", {shape[0], nbins}, {start1, 0}, {count1, nbins});

        var_u_bins = pdfwriter_io.DefineVariable<double>("U/bins", {nbins},
                                                         {0}, {nbins});
        var_v_bins = pdfwriter_io.DefineVariable<double>("V/bins", {nbins},
                                                         {0}, {nbins});
        var_step_out = pdfwriter_io.DefineVariable<int>("step");
    }

    pdfwriter.BeginStep();
    //check data
    std::cout << "u size " << pdf_u.size() << std::endl;
    pdfwriter.Put<double>(var_u_pdf, pdf_u.data());
    pdfwriter.Put<double>(var_v_pdf, pdf_v.data());

    pdfwriter.Put<double>(var_u_bins, bins_u.data());
    pdfwriter.Put<double>(var_v_bins, bins_v.data());
    pdfwriter.Put<int>(var_step_out, simStep);

    //if (write_inputvars)
    //{
    //    writer.Put<double>(var_u_out, u.data());
    //    writer.Put<double>(var_v_out, v.data());
    //}
    pdfwriter.EndStep();

    return;
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

    Settings settings = Settings::from_json(argv[1]);

    GrayScott sim(settings, comm);
    sim.init();

    adios2::ADIOS adios(settings.adios_config, comm, adios2::DebugON);

    adios2::IO io_ckpt = adios.DeclareIO("SimulationCheckpoint");

    //the wrapper for the writer enginge
    Writer writer_ckpt(settings, sim, io_ckpt);

    adios2::IO pdfwriter_io = adios.DeclareIO("PDFAnalysisOutput");

    std::string out_filename = "sim-pdf.bp";

    adios2::Engine pdfwriter =
        pdfwriter_io.Open(out_filename, adios2::Mode::Write, comm);

    if (rank == 0)
    {
        print_settings(settings);
        print_simulator_settings(sim);
        std::cout << "========================================" << std::endl;

        std::cout << "PDF analysis reads from Simulation use function call:  " << std::endl;
        std::cout << "PDF analysis writes using engine type:                 "
                  << pdfwriter_io.EngineType() << std::endl;
    }

#ifdef ENABLE_TIMERS
    Timer timer_total;
    Timer timer_compute;
    Timer timer_write;

    std::ostringstream log_fname;
    log_fname << "gray_scott_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    log << "step\ttotal_gs\tcompute_gs\twrite_gs" << std::endl;
#endif

    bool iffirstStep = true;

    for (int i = 0; i < settings.steps;)
    {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_total.start();
        timer_compute.start();
#endif

        for (int j = 0; j < settings.plotgap; j++)
        {
            sim.iterate();
            i++;
        }

#ifdef ENABLE_TIMERS
        double time_compute = timer_compute.stop();
        MPI_Barrier(comm);
        timer_write.start();
#endif

        if (rank == 0)
        {
            std::cout << "Simulation at step " << i
                      << " writing output step     " << i / settings.plotgap
                      << std::endl;
        }

        if (settings.checkpoint &&
            i % (settings.plotgap * settings.checkpoint_freq) == 0)
        {
            writer_ckpt.open(settings.checkpoint_output);
            writer_ckpt.write(i, sim);
            writer_ckpt.close();
        }

#ifdef ENABLE_TIMERS
        double time_write = timer_write.stop();
        double time_step = timer_total.stop();
        MPI_Barrier(comm);

        log << i << "\t" << time_step << "\t" << time_compute << "\t"
            << time_write << std::endl;
#endif

        //if the inline engine is used, read data and generate the vtkm data here
        //the adis needed to be installed before using

        if (i == settings.plotgap)
        {
            std::cout << "---using the inline engine, caculate the pdf---" << std::endl;
        }
        else
        {
            iffirstStep = false;
        }

        DoCheck(pdfwriter_io, pdfwriter, rank, i, iffirstStep, sim);
        std::cout << "ok for do check at ts: " << i << std::endl;
    }

    pdfwriter.Close();

#ifdef ENABLE_TIMERS
    log << "total\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed()
        << "\t" << timer_write.elapsed() << std::endl;

    log.close();
#endif

    MPI_Finalize();
}
