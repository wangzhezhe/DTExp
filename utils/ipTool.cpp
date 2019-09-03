
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <net/if.h>
#include <unistd.h>
#include <arpa/inet.h>

#include <fstream>
#include <iostream>

#include <dirent.h>
#include <system_error>

#include <sys/stat.h>
#include <vector>
#include <map>
#include <dirent.h>
#include <system_error>
#include <cstring>
#include <set>

#include "ipTool.h"


using namespace std;

vector<string> loadAddrInDir(string Dir)
{
    // only init the client in clusterDir
    string dir = Dir;
    DIR *dp;
    vector<string> AddrList;
    struct dirent *entry;
    struct stat statbuf;
    try
    {
        if ((dp = opendir(dir.data())) == NULL)
        {
            printf("Can`t open directory %s\n", dir.data());
            return AddrList;
        }

        while ((entry = readdir(dp)) != NULL)
        {

            if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
            {
                continue;
            }

            //printf("load addr name %s\n", entry->d_name);
            AddrList.push_back(string(entry->d_name));
        }
        //change to the upper dir
        //printf("debug id %d load finish\n", GETIPCOMPONENTID);
        closedir(dp);
        return AddrList;
    }
    catch (const std::system_error &e)
    {
        std::cout << "getip Caught system_error with code " << e.code()
                  << " meaning " << e.what() << '\n';
        exit(1);
    }
}

void recordIPPortWithoutFile(string &ipstr, string INTERFACE)
{
    int n;
    struct ifreq ifr;
    //assume the network interface exist

    n = socket(AF_INET, SOCK_DGRAM, 0);
    //Type of address to retrieve - IPv4 IP address
    ifr.ifr_addr.sa_family = AF_INET;
    //Copy the interface name in the ifreq structure
    strncpy(ifr.ifr_name, INTERFACE.data(), IFNAMSIZ - 1);
    ioctl(n, SIOCGIFADDR, &ifr);
    close(n);
    //display result

    char *ip = inet_ntoa(((struct sockaddr_in *)&ifr.ifr_addr)->sin_addr);

    ipstr = string(ip);
    return;
}

void deleteDir(string DirPath)
{

    if (rmdir(DirPath.data()) != 0)
    {
        throw std::runtime_error("failed to remove dir");
    }

    exit(0);
}

void createDir(string DirPath)
{

    //if there is no Dir in current Dir Path, create
    fstream dirfile;
    dirfile.open(DirPath.data(), ios::in);
    if (!dirfile)
    {
        while (1)
        {
            const int dir_err = mkdir(DirPath.data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (dir_err < 0)
            {
                printf("Error creating directory %s, retry...", DirPath.data());
                usleep(10000);
                dirfile.open(DirPath.data(), ios::in);
                if (dirfile)
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }
    dirfile.close();
    return;
}

int getFreePortNum()
{

    int socket_desc, new_socket, c;
    struct sockaddr_in server;
    char message[100];

    //Create socket
    socket_desc = socket(AF_INET, SOCK_STREAM, 0);
    if (socket_desc == -1)
    {
        printf("Could not create socket");
    }

    //Prepare the sockaddr_in structure
    server.sin_family = AF_INET;
    server.sin_addr.s_addr = INADDR_ANY;
    server.sin_port = htons(0);

    //server function: bind, listen , accept

    //Bind

    if (bind(socket_desc, (struct sockaddr *)&server, sizeof(server)) < 0)
    {
        puts("bind failed");
        return 1;
    }

    socklen_t len = sizeof(server);

    if (getsockname(socket_desc, (struct sockaddr *)&server, &len) != -1)
    {

        //printf("port number %d\n", ntohs(server.sin_port));
        close(socket_desc);
        return int(ntohs(server.sin_port));
    }
    else
    {

        printf("failed to get port");
        return -1;
    }
}

string parsePort(string peerURL)
{

    int startPosition = peerURL.find("ipv4:");

    //delete ipv4:
    peerURL.erase(startPosition, 5);

    //delete the port suffix
    int len = peerURL.length();

    int endPosition = peerURL.find(":");

    //printf("start position of : (%d) original str %s\n", startPosition, peerURL.data());

    peerURL.erase(0, endPosition + 1);

    return peerURL;
}

string parseIP(string peerURL)
{

    int startPosition = peerURL.find("ipv4:");

    //delete ipv4:
    peerURL.erase(startPosition, 5);

    //delete the port suffix
    int len = peerURL.length();

    startPosition = peerURL.find(":");

    //printf("start position of : (%d) original str %s\n", startPosition, peerURL.data());

    peerURL.erase(startPosition, len - startPosition);

    return peerURL;
}