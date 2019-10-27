#include <putgetMeta/metaclient.h>
#include <utils/ThreadPool.h>
#include <unistd.h>
#include <queue>
#include <thread>

#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include "../simulation/settings.h"
#include "../common/timer.hpp"
#include "adios2.h"

std::string keyDataOk = "INDICATOR1";
std::string keySimFinish = "SIMFINISH";
std::queue<std::future<void>> results;

typedef struct partitionInfo
{
    size_t npx;
    size_t npy;
    size_t npz;

    size_t px;
    size_t py;
    size_t pz;

} partitionInfo;

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

void pullandStartAna(int ts, int rank, partitionInfo &pinfo, adios2::IO inIO, adios2::Engine reader)
{

    //todo use the setselection to get the meta and start the analytics
    //indicate anyway
    //caculate pdf
    //todo add the checking operation for the pdf
    //todo get var_u

    std::vector<double> u;

    adios2::Variable<double> varU = inIO.InquireVariable<double>("U");

    adios2::Dims shape = varU.Shape();

    size_t npx = pinfo.npx;
    size_t npy = pinfo.npy;
    size_t npz = pinfo.npz;

    size_t px = pinfo.px;
    size_t py = pinfo.py;
    size_t pz = pinfo.pz;

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

    varU.SetStepSelection({ts, 1});

    reader.Get<double>(varU, u, adios2::Mode::Sync);

    int isovalue = 0.5;

    auto polyData = compute_isosurface(varU, u, isovalue);

    //write_adios(writer, polyData, varPoint, varCell, varNormal, varOutStep,
    //            step, comm);
    //replase by sleep for experiment to avoid the impact of IO
    Settings settings = Settings::from_json("./settings.json");
    //sim<ana+c
    //std::this_thread::sleep_for(std::chrono::milliseconds(10 * settings.L));

    //sim>ana+c
    std::this_thread::sleep_for(std::chrono::milliseconds(3*settings.L));

    /*
    std::string dir = "./vtkdata";

    char countstr[50];
    sprintf(countstr, "%02d_%04d", rank, ts);

    std::string fname = dir + "/vtkiso_" + std::string(countstr) + ".vtk";
    //the format here is the vtk
    write_vtk(fname, polyData);
    */
    std::cout << "pullandStartAna,ok for ts: " << ts << std::endl;

    //sleep adjusted time
    /*
    Settings settings = Settings::from_json("./settings.json");
    usleep(1000 * 5 * settings.L);
    */
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

    std::cout << "ana info: rank " << rank << " proces " << procs << std::endl;

    int dims[3] = {0};
    MPI_Dims_create(procs, 3, dims);
    partitionInfo pinfo;
    pinfo.npx = dims[0];
    pinfo.npy = dims[1];
    pinfo.npz = dims[2];

    int coords[3] = {0};
    int periods[3] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    pinfo.px = coords[0];
    pinfo.py = coords[1];
    pinfo.pz = coords[2];

    if (argc != 2)
    {
        std::cout << "<binary> eventStr" << std::endl;
        return 0;
    }

    MetaClient metaclient = getMetaClient();
    std::string recordkey = "process" + std::string(argv[1]);
    if (rank == 0)
    {
        metaclient.Recordtime(recordkey);
    }

    //parse the event notification
    int ts = std::stoi(argv[1]);

    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO inIO = adios.DeclareIO("SimulationOutput");
    std::string input_fname = "gs.bp";

    adios2::Engine reader = inIO.Open(input_fname, adios2::Mode::Read);

    reader.BeginStep();
    pullandStartAna(ts, rank, pinfo, inIO, reader);
    reader.EndStep();
    reader.Close();

    //tick finish
    if (rank == 0)
    {

        string reply = metaclient.Recordtimetick("WFTIMER");
        std::cout << "Timer received for anastartbyevent finish: " << reply << std::endl;

        metaclient.Recordtime(recordkey);
    }

    MPI_Finalize();

    return 0;
}