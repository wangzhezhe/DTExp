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

#include "../simulation/settings.h"

std::string keyDataBase = "INDICATOR";
std::string keySimFinish = "SIMFINISH";
std::queue<std::future<void>> results;

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


void startLocalTask(int ts)
{
    char command[200];
    sprintf(command, "%s %d", "srun --mpi=pmix_v2 --mem-per-cpu=1000 -n 8 ./anastartbyevent", ts);
    printf("execute command by local way:(%s)\n", command);
    //test using
    system(command);
    return;
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

void checkMetaAndEnqueue(ThreadPool &pool, int taskNum, MetaClient &metaclient)
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
                        pool.enqueue([ts] {
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

    if (argc != 3)
    {
        std::cout << "<binary> <poolsize> <processId>" << std::endl;
        return 0;
    }

    int poolSize = std::stoi(argv[1]);
    int taskNum = std::stoi(argv[2]);

    ThreadPool pool(poolSize);
    MetaClient metaclient = getMetaClient();

    std::thread checkMetaThread(checkMetaAndEnqueue, std::ref(pool), taskNum, std::ref(metaclient));
    std::thread checkResultsThread(checkResults, std::ref(pool), std::ref(metaclient));

    checkMetaThread.join();
    checkResultsThread.join();

    //not sure why there is bug when add this
    //reader.EndStep();
    //reader.Close();

    //tick finish
    
    string reply = metaclient.Recordtimetick("WFTIMER");
    std::cout << "anapullmeta finish reply: " << reply << std::endl;
   


    return 0;
}
