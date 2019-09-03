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

std::string keyDataOk = "INDICATOR1";
std::string keySimFinish = "SIMFINISH";
std::queue<std::future<void>> results;


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


void pullandStartAna(int ts)
{

    //todo use the setselection to get the meta and start the analytics
    //indicate anyway
    //caculate pdf
    //todo add the checking operation for the pdf
    //todo get var_u
    int isovalue = 0.5;

    auto polyData = compute_isosurface(var_u_in, u, isovalue);

    //write_adios(writer, polyData, varPoint, varCell, varNormal, varOutStep,
    //            step, comm);
    std::string dir = "./vtkdata";

    char countstr[50];
    sprintf(countstr, "%04d", simStep);

    std::string fname = dir + "/vtkiso_" + std::string(countstr) + ".vtk";
    //the format here is the vtk
    write_vtk(fname, polyData);
    std::cout << "ok for ts " << simStep << std::endl;

    //sleep adjusted time
    sleep(1);
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
    }
}

void checkMetaAndEnqueue(ThreadPool &pool, MetaClient &metaclient)
{
    //check the meta when there is info, then trigure specific tasks
    while (true)
    {
        //check finish is empty
        //return

        std::string simFinishReply = metaclient.Getmetaspace(keySimFinish);
        if (simFinishReply.compare("OK") == 0 && pool.getTaskSize() == 0)
        {

            std::cout << "simulation finish and all tasks finish" << std::endl;
            break;
        }
        else
        {

            //check the ts
            std::string indicatorReply = metaclient.Getmeta(keyDataOk);

            if (indicatorReply.compare("NULL") != 0)
            {
                std::cout << "process data with ts: " << indicatorReply << std::endl;
                int ts = std::stoi(indicatorReply);

                results.push(
                    pool.enqueue([ts] {
                        return pullandStartAna(ts);
                    }));
            }
            else
            {
                sleep(0.05);
            }
        }
    }
    return;
}

int main(int argc, char *argv[])
{

    std::cout << "Program finish request recieve: " << reply << std::endl;

    if (argc != 2)
    {
        std::cout << "<binary> <poolsize>" << std::endl;
        return 0;
    }

    int poolSize = std::stoi(argv[1]);

    //init the thread pool

    ThreadPool pool(poolSize);

    std::thread checkMetaThread(checkMetaAndEnqueue, std::ref(pool), std::ref(metaclient));
    std::thread checkResultsThread(checkResults, std::ref(pool), std::ref(metaclient));

    checkMetaThread.join();
    checkResultsThread.join();

    return 0;
}