
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

#include "metaserver.grpc.pb.h"
#include "metaclient.h"

const string metaserverDir = "./Metaserver";

std::string MetaClient::Recordtime(const std::string &key)
{
    // Data we are sending to the server.
    TimeRequest request;
    request.set_key(key);

    // Container for the data we expect from the server.
    TimeReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Recordtime(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Recordtimestart(const std::string &key)
{
    // Data we are sending to the server.
    TimeRequest request;
    request.set_key(key);

    // Container for the data we expect from the server.
    TimeReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Recordtimestart(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Recordtimetick(const std::string &key)
{
    // Data we are sending to the server.
    TimeRequest request;
    request.set_key(key);

    // Container for the data we expect from the server.
    TimeReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Recordtimetick(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Getmeta(const std::string &key)
{
    // Data we are sending to the server.
    GetRequest request;
    request.set_key(key);

    // Container for the data we expect from the server.
    GetReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Getmeta(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Putmeta(const std::string &key, const std::string &value)
{
    // Data we are sending to the server.
    PutRequest request;
    request.set_key(key);
    request.set_value(value);

    // Container for the data we expect from the server.
    PutReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Putmeta(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Getmetaspace(const std::string &key)
{
    // Data we are sending to the server.
    GetRequest request;
    request.set_key(key);

    // Container for the data we expect from the server.
    GetReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Getmetaspace(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

std::string MetaClient::Putmetaspace(const std::string &key, const std::string &value)
{
    // Data we are sending to the server.
    PutRequest request;
    request.set_key(key);
    request.set_value(value);

    // Container for the data we expect from the server.
    PutReply reply;

    // Context for the client. It could be used to convey extra information to
    // the server and/or tweak certain RPC behaviors.
    ClientContext context;

    // The actual RPC.
    Status status = stub_->Putmetaspace(&context, request, &reply);

    // Act upon its status.
    if (status.ok())
    {
        return reply.message();
    }
    else
    {
        std::cout << status.error_code() << ": " << status.error_message()
                  << std::endl;
        return "RPC failed";
    }
}

string getAddr()
{

    vector<string> addrList = loadAddrInDir(metaserverDir);
    if (addrList.size() > 0)
    {
        return addrList[0];
    }

    return "";
}

void recordKey(string key)
{

    string serverAddr = getAddr();

    if (serverAddr == "")
    {
        printf("failed to get server addr\n");
        return;
    }

    MetaClient metaclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));

    string reply = metaclient.Recordtime(key);
    std::cout << "Timer received: " << reply << std::endl;

    return;
}

void recordKeyStart(string key)
{
    //if use parallel file system
    //load the addr from the file
    string serverAddr = getAddr();

    if (serverAddr == "")
    {
        printf("failed to get server addr\n");
        return;
    }

    MetaClient metaclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));

    string reply = metaclient.Recordtimestart(key);
    std::cout << "Timer received: " << reply << std::endl;

    return;
}

void recordKeyTick(string key)
{

    string serverAddr = getAddr();

    if (serverAddr == "")
    {
        printf("failed to get server addr\n");
        return;
    }

    MetaClient metaclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));

    string reply = metaclient.Recordtimetick(key);
    std::cout << "Timer received: " << reply << std::endl;

    return;
}

MetaClient getMetaClient()
{
    string serverAddr = getAddr();

    if (serverAddr == "")
    {

        throw std::runtime_error("the addr for MMserver should not be empty");
    }

    MetaClient metaclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));

    return metaclient;
}
/*
int main(int argc, char **argv)
{
    // Instantiate the client. It requires a channel, out of which the actual RPCs
    // are created. This channel models a connection to an endpoint (in this case,
    // localhost at port 50051). We indicate that the channel isn't authenticated
    // (use of InsecureChannelCredentials()).

    string serverAddr = getAddr();

    if (serverAddr == "")
    {

        printf("server addr is 0\n");
        return 0;
    }

    MetaClient metaclient(grpc::CreateChannel(
        serverAddr, grpc::InsecureChannelCredentials()));
    //std::string user("world");
    //std::string reply = metaclient.SayHello(user);
    //std::cout << "Greeter received: " << reply << std::endl;

    std::string key("world");
    string reply = metaclient.Recordtime(key);
    std::cout << "Timer received: " << reply << std::endl;

    sleep(1);

    reply = metaclient.Recordtime(key);
    std::cout << "Timer received: " << reply << std::endl;

    key = "pattern1";
    string meta1("meta1");
    string meta2("meta2");
    string meta3("meta3");

    reply = metaclient.Getmeta(key);
    std::cout << "Get pattern1 recieve: " << reply << std::endl;

    reply = metaclient.Putmeta(key, meta1);
    std::cout << "Put meta1 recieve: " << reply << std::endl;
    reply = metaclient.Putmeta(key, meta2);
    std::cout << "Put meta2 recieve: " << reply << std::endl;
    reply = metaclient.Putmeta(key, meta3);
    std::cout << "Put meta3 recieve: " << reply << std::endl;

    reply = metaclient.Getmeta(key);
    std::cout << "Get pattern1 recieve: " << reply << std::endl;
    reply = metaclient.Getmeta(key);
    std::cout << "Get pattern1 recieve: " << reply << std::endl;
    reply = metaclient.Getmeta(key);
    std::cout << "Get pattern1 recieve: " << reply << std::endl;
    reply = metaclient.Getmeta(key);
    std::cout << "Get pattern1 recieve: " << reply << std::endl;

    printf("------test tick------\n");
    key = "pattern_tick";
    reply = metaclient.Recordtimestart(key);
    std::cout << "Put pattern_tick 1st recieve: " << reply << std::endl;
    reply = metaclient.Recordtimestart(key);
    std::cout << "Put pattern_tick 2st recieve: " << reply << std::endl;
    reply = metaclient.Recordtimetick(key);
    std::cout << "Put pattern_tick 3st recieve: " << reply << std::endl;

    reply = metaclient.Recordtimetick(key);
    std::cout << "Put pattern_tick 4st recieve: " << reply << std::endl;

    return 0;
}
*/