
/*
 *
 * Copyright 2015 gRPC authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <memory>
#include <string>
#include <grpcpp/grpcpp.h>
#include "unistd.h"
#include "eventclient.h"
#include "utils/ipTool.h"

const std::string metaserverDir = "./Metaserver";

//using matchType equals NAME here
std::string EventClient::Publish(std::vector<std::string> eventList, std::string source, std::string metadata, std::string matchType)
{
    // Container for the data we expect from the server.
    PubSubRequest request;
    PubSubReply reply;
    int size = eventList.size();
    int i = 0;
    for (i = 0; i < size; i++)
    {
        //attention the use here, the request could be transfered into a specific type with specific function
        request.add_pubsubmessage(eventList[i]);
    }

    request.set_source(source);

    request.set_metadata(metadata);

    request.set_matchtype(matchType);

    ClientContext context;

    Status status = stub_->Publish(&context, request, &reply);
    
    if (status.ok())
    {
        return reply.returnmessage();
    }
    else
    {
        std::cout << "rpc fail " << status.error_code() << ": " << status.error_message()
             << std::endl;

        return "RPC failed";
    }
    return reply.returnmessage();
}

std::string getserverAddr(){
    //defualt addr is 
    const std::string defaultWorkflowServerAddr = "./multinodeip/cluster0/coordinator/";
    std::vector<std::string> addrList = loadAddrInDir(defaultWorkflowServerAddr);
    return addrList[0];
}

EventClient getEventClient()
{
    std::string serverAddr = getserverAddr();

    if (serverAddr == "")
    {

        throw std::runtime_error("the addr for MMserver should not be empty");
    }

    EventClient eventclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));

    return eventclient;
}