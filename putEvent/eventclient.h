

#ifndef PUBSUBCLIENT_H
#define PUBSUBCLIENT_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <unistd.h>
#include <pthread.h>

#include <grpcpp/grpcpp.h>

#include "workflowserver.grpc.pb.h"

using grpc::Channel;
using grpc::ClientContext;
using grpc::Status;
using workflowserver::Greeter;
using workflowserver::HelloReply;
using workflowserver::HelloRequest;
using workflowserver::PubSubReply;
using workflowserver::PubSubRequest;


class EventClient
{
public:
    EventClient(std::shared_ptr<Channel> channel)
        : stub_(Greeter::NewStub(channel)) {}

    std::string Publish(std::vector<std::string> eventList, std::string source, std::string metadata, std::string matchType);

private:
    std::unique_ptr<Greeter::Stub> stub_;
};


EventClient getEventClient();


#endif