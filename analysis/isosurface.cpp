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

void write_vtk(const std::string &fname,
               const vtkSmartPointer<vtkPolyData> polyData)
{
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fname.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

void write_adios(adios2::Engine &writer,
                 const vtkSmartPointer<vtkPolyData> polyData,
                 adios2::Variable<double> &varPoint,
                 adios2::Variable<int> &varCell,
                 adios2::Variable<double> &varNormal,
                 adios2::Variable<int> &varOutStep, int step, MPI_Comm comm)
{
    int numCells = polyData->GetNumberOfPolys();
    int numPoints = polyData->GetNumberOfPoints();
    int rank;

    MPI_Comm_rank(comm, &rank);

    std::vector<double> points(numPoints * 3);
    std::vector<double> normals(numPoints * 3);
    std::vector<int> cells(numCells * 3); // Assumes that cells are triangles

    double coords[3];

    auto cellArray = polyData->GetPolys();

    cellArray->InitTraversal();

    // Iterate through cells
    for (int i = 0; i < polyData->GetNumberOfPolys(); i++)
    {
        auto idList = vtkSmartPointer<vtkIdList>::New();

        cellArray->GetNextCell(idList);

        // Iterate through points of a cell
        for (int j = 0; j < idList->GetNumberOfIds(); j++)
        {
            auto id = idList->GetId(j);

            cells[i * 3 + j] = id;

            polyData->GetPoint(id, coords);

            points[id * 3 + 0] = coords[0];
            points[id * 3 + 1] = coords[1];
            points[id * 3 + 2] = coords[2];
        }
    }

    auto normalArray = polyData->GetPointData()->GetNormals();

    // Extract normals
    for (int i = 0; i < normalArray->GetNumberOfTuples(); i++)
    {
        normalArray->GetTuple(i, coords);

        normals[i * 3 + 0] = coords[0];
        normals[i * 3 + 1] = coords[1];
        normals[i * 3 + 2] = coords[2];
    }

    int totalPoints, offsetPoints;
    MPI_Allreduce(&numPoints, &totalPoints, 1, MPI_INT, MPI_SUM, comm);
    MPI_Scan(&numPoints, &offsetPoints, 1, MPI_INT, MPI_SUM, comm);

    writer.BeginStep();

    varPoint.SetShape({static_cast<size_t>(totalPoints),
                       static_cast<size_t>(totalPoints > 0 ? 3 : 0)});
    varPoint.SetSelection({{static_cast<size_t>(offsetPoints - numPoints), 0},
                           {static_cast<size_t>(numPoints),
                            static_cast<size_t>(numPoints > 0 ? 3 : 0)}});

    varNormal.SetShape(varPoint.Shape());
    varNormal.SetSelection({varPoint.Start(), varPoint.Count()});

    if (numPoints)
    {
        writer.Put(varPoint, points.data());
        writer.Put(varNormal, normals.data());
    }

    int totalCells, offsetCells;
    MPI_Allreduce(&numCells, &totalCells, 1, MPI_INT, MPI_SUM, comm);
    MPI_Scan(&numCells, &offsetCells, 1, MPI_INT, MPI_SUM, comm);

    for (int i = 0; i < cells.size(); i++)
    {
        cells[i] += (offsetPoints - numPoints);
    }

    varCell.SetShape({static_cast<size_t>(totalCells),
                      static_cast<size_t>(totalCells > 0 ? 3 : 0)});
    varCell.SetSelection({{static_cast<size_t>(offsetCells - numCells), 0},
                          {static_cast<size_t>(numCells),
                           static_cast<size_t>(numCells > 0 ? 3 : 0)}});

    if (numCells)
    {
        writer.Put(varCell, cells.data());
    }

    if (!rank)
    {
        std::cout << "isosurface at step " << step << " writing out "
                  << totalCells << " cells and " << totalPoints << " points"
                  << std::endl;
    }

    writer.Put(varOutStep, step);

    writer.EndStep();
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 5;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);


    int dims[3] = {0};
    MPI_Dims_create(procs, 3, dims);
    size_t npx = dims[0];
    size_t npy = dims[1];
    size_t npz = dims[2];

    int coords[3] = {0};
    int periods[3] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    size_t px = coords[0];
    size_t py = coords[1];
    size_t pz = coords[2];


    if (argc < 4)
    {
        if (rank == 0)
        {
            std::cerr << "Too few arguments" << std::endl;
            std::cout << "Usage: isosurface input output isovalue" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    const std::string input_fname(argv[1]);
    const std::string output_fname(argv[2]);
    const double isovalue = std::stod(argv[3]);


    adios2::ADIOS adios("adios2.xml", comm, adios2::DebugON);

    adios2::IO inIO = adios.DeclareIO("SimulationOutput");
    
    adios2::Engine reader = inIO.Open(input_fname, adios2::Mode::Read);
    
    
    std::vector<double> u;
    int step;

#ifdef ENABLE_TIMERS
    Timer timer_total;
    Timer timer_read;
    Timer timer_compute;
    Timer timer_write;

    std::ostringstream log_fname;
    log_fname << "isosurface_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    //log << "step\ttotal_iso\tread_iso\tcompute_write_iso" << std::endl;
#endif

    if (rank == 0)
    {
        MetaClient metaclient = getMetaClient();
        std::string reply = metaclient.Recordtimestart("PROCESS");
        std::cout << "start PROCESS reply " << reply <<std::endl;
    }

    while (true)
    {
#ifdef ENABLE_TIMERS

        MPI_Barrier(comm);
        timer_total.start();
        timer_read.start();
#endif
        adios2::StepStatus read_status = reader.BeginStep();

        if (read_status == adios2::StepStatus::OtherError)
        {
            std::cout << "---in iso surface adios status is unknown---" << read_status << std::endl;
            break;
        }

        if (read_status != adios2::StepStatus::OK)
        {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout << "get status for iso " << read_status << std::endl;
            if (read_status == adios2::StepStatus::EndOfStream)
            {
                break;
            }
            continue;
        }


        adios2::Variable<double> varU = inIO.InquireVariable<double>("U");
        const adios2::Variable<int> varStep = inIO.InquireVariable<int>("step");

        adios2::Dims shape = varU.Shape();

        size_t size_x = (shape[0] + npx - 1) / npx;
        size_t size_y = (shape[1] + npy - 1) / npy;
        size_t size_z = (shape[2] + npz - 1) / npz;

        size_t offset_x = size_x * px;
        size_t offset_y = size_y * py;
        size_t offset_z = size_z * pz;

        if (px == npx - 1)
        {
            size_x -= size_x * npx - shape[0];
        }
        if (py == npy - 1)
        {
            size_y -= size_y * npy - shape[1];
        }
        if (pz == npz - 1)
        {
            size_z -= size_z * npz - shape[2];
        }

        varU.SetSelection({{offset_x, offset_y, offset_z},
                           {size_x + (px != npx - 1 ? 1 : 0),
                            size_y + (py != npy - 1 ? 1 : 0),
                            size_z + (pz != npz - 1 ? 1 : 0)}});

        reader.Get<double>(varU, u);
        reader.Get<int>(varStep, step);
        reader.EndStep();
        if (rank==0){
        MetaClient metaclient = getMetaClient();
        std::string reply = metaclient.Recordtimetick("PROCESS");
        std::cout << "ok for getting data " << step << " tick reply: " << reply << std::endl;
        }


#ifdef ENABLE_TIMERS
        double time_read = timer_read.stop();
        MPI_Barrier(comm);
        timer_compute.start();
#endif

        //todo add the checking operation for the pdf
        auto polyData = compute_isosurface(varU, u, isovalue);

        //write_adios(writer, polyData, varPoint, varCell, varNormal, varOutStep,
        //            step, comm);
        //replase by sleep for experiment to avoid the impact of IO
        Settings settings = Settings::from_json("./settings.json");

        //sim<ana+c
        std::this_thread::sleep_for(std::chrono::milliseconds(10*settings.L));

        //sim>ana+c
        //std::this_thread::sleep_for(std::chrono::milliseconds(3*settings.L));

        /*
        std::string dir = "./vtkdata";

        char countstr[50];
        sprintf(countstr, "%04d", step);

        std::string fname = dir + "/vtkiso_" + std::string(countstr) + ".vtk";
        //the format here is the vtk
        write_vtk(fname, polyData);
        */
        std::cout << "ok for ts " << step << std::endl;

#ifdef ENABLE_TIMERS
        double time_compute = timer_compute.stop();
        MPI_Barrier(comm);
        //timer_write.start();
#endif

#ifdef ENABLE_TIMERS
        //double time_write = timer_write.stop();
        double time_total_t = timer_total.stop();
        MPI_Barrier(comm);

        std::cout << "timestep: " << step << " total: " << time_total_t << " "
                  << " read: " << time_read << " "
                  << "compute: " << time_compute << std::endl;
#endif
    }

#ifdef ENABLE_TIMERS
    //log << "total\t" << timer_total.elapsed() << "\t" << timer_read.elapsed()
    //   << "\t" << timer_compute.elapsed() << std::endl;

    log.close();
#endif

    reader.Close();

    if (rank == 0)
    {
        //add adjusted time
        /*
        Settings settings = Settings::from_json("./settings.json");
        usleep(1000 * 5 * settings.L);
        std::cout << "adjusted time: " << 1000 * 5 * settings.L << std::endl;
        */
        //tick finish
        MetaClient metaclient = getMetaClient();
        std::string reply = metaclient.Recordtimetick("WFTIMER");
        std::cout << "tick reply: " << reply << std::endl;
    }

    MPI_Finalize();

    return 0;
}
