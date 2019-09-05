#include <putgetMeta/metaclient.h>
#include <utils/ThreadPool.h>
#include <unistd.h>
#include <stdlib.h>
#include <queue>
#include <thread>
#include <mutex>

#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include "adios2.h"
#include "../simulation/settings.h"

std::string keyDataBase = "INDICATOR";
std::string keySimFinish = "SIMFINISH";
std::queue<std::future<void>> results;

size_t npx;
size_t npy;
size_t npz;

size_t px;
size_t py;
size_t pz;

std::mutex taskNeedToFinishMutex;
size_t taskNeedToFinish = 0;

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

void startLocalTask(int ts)
{
    char command[200];
    sprintf(command, "%s %d", "mpirun -n 1 ./anastartbyevent", ts);
    printf("execute command by local way:(%s)\n", command);
    //test using
    system(command);
    return;
}

void pullandStartAna(int ts, adios2::IO inIO, adios2::Engine reader)
{

    std::vector<double> u;

    adios2::Variable<double> varU = inIO.InquireVariable<double>("U");

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
    //TODO, there are some bugs if set this value as ts
    //use 0 to replace here, this will not influence the testing
    varU.SetStepSelection({0, 1});

    reader.Get<double>(varU, u, adios2::Mode::Sync);

    int isovalue = 0.5;

    auto polyData = compute_isosurface(varU, u, isovalue);

    //write_adios(writer, polyData, varPoint, varCell, varNormal, varOutStep,
    //            step, comm);
    std::string dir = "./vtkdata";

    char countstr[50];
    sprintf(countstr, "%04d", ts);

    std::string fname = dir + "/vtkiso_" + std::string(countstr) + ".vtk";
    //the format here is the vtk
    write_vtk(fname, polyData);
    std::cout << "pullandStartAna, ok for ts: " << ts << std::endl;

    //sleep adjusted time
    Settings settings = Settings::from_json("./settings.json");
    usleep(1000 * 5 * settings.L);
    std::cout << "adjusted time " << 5 * settings.L << std::endl;

    taskNeedToFinishMutex.lock();
    taskNeedToFinish--;
    taskNeedToFinishMutex.unlock();
}

void checkResults(ThreadPool &pool, MetaClient &metaclient)
{

    while (true)
    {
        if (pool.getTaskSize() > 0)
        {
            results.front().get();
            results.pop();
        }
        else
        {
            break;
        }

        sleep(0.1);
    }
}

void checkMetaAndEnqueue(ThreadPool &pool, int taskNum, MetaClient &metaclient, const adios2::IO &inIO, adios2::Engine &reader)
{
    //check the meta when there is info, then trigure specific tasks
    while (true)
    {
        //check finish is empty
        //return

        std::string simFinishReply = metaclient.Getmetaspace(keySimFinish);

        //Attention!!! remeber to check the finished task !!! instead of only the pool size
        if (simFinishReply.compare("OK") == 0 && pool.getTaskSize() == 0 && taskNeedToFinish == 0)
        {

            std::cout << "simulation finish and all tasks finish" << std::endl;
            break;
        }
        else
        {
            for (int i = 0; i < taskNum; i++)

            {
                std::string keyDataOk = keyDataBase + std::to_string(i);
                //check the ts
                std::string indicatorReply = metaclient.Getmeta(keyDataOk);

                if (indicatorReply.compare("NULL") != 0)
                {
                    
                    std::cout << "process data with key " << keyDataOk << "for ts: " << indicatorReply << std::endl;
                    int ts = std::stoi(indicatorReply);
                    taskNeedToFinishMutex.lock();
                    taskNeedToFinish++;
                    taskNeedToFinishMutex.unlock();
                    results.push(
                        pool.enqueue([ts, inIO, reader] {
                            //return pullandStartAna(ts, inIO, reader);
                            //ana is started by separate program
                            return startLocalTask(ts);
                        }));
                }
                else
                {
                    // std::cout << "simFinishReply: " << simFinishReply <<std::endl;
                    // std::cout << "pool.getTaskSize()" << pool.getTaskSize() <<std::endl;
                    //std::cout << "taskNeedToFinish" << taskNeedToFinish <<std::endl;
                    //std::cout << "simFinishReply: " << simFinishReply <<std::endl;
                    //std::cout << "pool.getTaskSize()" << pool.getTaskSize() <<std::endl;
                    //std::cout << "while loop" <<std::endl;
                    //this time should be smaller than the frequency of event updating
                    sleep(0.01);
                }
            }
        }
    }
    return;
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

    std::cout << "ana info: rank " << rank << "proces " << procs << std::endl;

    int dims[3] = {0};
    MPI_Dims_create(procs, 3, dims);
    npx = dims[0];
    npy = dims[1];
    npz = dims[2];

    int coords[3] = {0};
    int periods[3] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    px = coords[0];
    py = coords[1];
    pz = coords[2];

    if (argc != 3)
    {
        std::cout << "<binary> <poolsize> <processId>" << std::endl;
        return 0;
    }

    int poolSize = std::stoi(argv[1]);
    int taskNum = std::stoi(argv[2]);

    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO inIO = adios.DeclareIO("SimulationOutput");
    std::string input_fname = "gs.bp";

    adios2::Engine reader = inIO.Open(input_fname, adios2::Mode::Read);

    adios2::StepStatus read_status =
        reader.BeginStep(adios2::StepMode::Read, 10.0f);

    while (true)
    {
        if (read_status == adios2::StepStatus::NotReady)
        {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK)
        {
            std::cout << "unknown read status: " << read_status << std::endl;
            return 0;
        }
        else
        {
            break;
        }
    }

    //init the thread pool

    ThreadPool pool(poolSize);
    MetaClient metaclient = getMetaClient();

    std::thread checkMetaThread(checkMetaAndEnqueue, std::ref(pool), taskNum, std::ref(metaclient), std::ref(inIO), std::ref(reader));
    std::thread checkResultsThread(checkResults, std::ref(pool), std::ref(metaclient));

    checkMetaThread.join();
    checkResultsThread.join();

    //not sure why there is bug when add this
    //reader.EndStep();
    //reader.Close();

    //tick finish
    string reply = metaclient.Recordtimetick("WFTIMER");

    return 0;
}
