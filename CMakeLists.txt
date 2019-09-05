cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
# ENABLE ADIOS2_ROOT
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project(gray-scott C CXX)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)


option(USE_TIMERS "Use profiling timers")
#if(USE_TIMERS)
  message(STATUS "Enabling profiling timers")
  add_definitions(-DENABLE_TIMERS)
#endif()

# We are not using the C++ API of MPI, this will stop the compiler look for it
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)


add_executable(pdf_calc analysis/pdf_calc.cpp)
target_link_libraries(pdf_calc adios2::adios2 MPI::MPI_C)



option(VTK "Build VTK apps")
if (VTK_ROOT)
  set(VTK ON)
endif(VTK_ROOT)

if (VTK)
  message(STATUS "Configuring VTK apps")

  find_package(VTK COMPONENTS
    FiltersCore
    IOImage
    IOXML
  )

  if(VTK_FOUND)

    add_executable(pdfwithana analysis/pdfwithana.cpp simulation/settings.cpp)
    target_link_libraries(pdfwithana adios2::adios2 ${VTK_LIBRARIES}
      MPI::MPI_C)

  endif(VTK_FOUND)

  find_package(VTK COMPONENTS
    FiltersCore
    FiltersGeometry
  )

  if(VTK_FOUND)
    add_executable(find_blobs analysis/find_blobs.cpp)
    target_link_libraries(find_blobs adios2::adios2 ${VTK_LIBRARIES}
      MPI::MPI_C)
  endif(VTK_FOUND)

  find_package(VTK COMPONENTS
    FiltersGeneral
  )

  if(VTK_FOUND)
    add_executable(compute_curvature analysis/curvature.cpp)
    target_link_libraries(compute_curvature adios2::adios2 ${VTK_LIBRARIES}
      MPI::MPI_C)
  endif(VTK_FOUND)


  find_package(VTK COMPONENTS
    RenderingOpenGL2
    ViewsInfovis
  )

  if(VTK_FOUND)
    add_executable(render_isosurface plot/render_isosurface.cpp)
    target_link_libraries(render_isosurface adios2::adios2 ${VTK_LIBRARIES}
      MPI::MPI_C)
  endif(VTK_FOUND)
endif(VTK)

# for meta client


# GRPC and Protocol Buffers libraries location
list(APPEND CMAKE_PREFIX_PATH "~/cworkspace/build/build_grpc" "~/cworkspace/build/build_protobuf")

# Cmake find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Protobuf REQUIRED)
find_package(GRPC REQUIRED)

set(PROTOS
    ${CMAKE_CURRENT_SOURCE_DIR}/protos/metaserver.proto
    ${CMAKE_CURRENT_SOURCE_DIR}/protos/workflowserver.proto
)

set(PROTO_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/proto-src)
file(MAKE_DIRECTORY ${PROTO_SRC_DIR})
include_directories(${PROTO_SRC_DIR})

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${PROTO_SRC_DIR} ${PROTOS})
grpc_generate_cpp(GRPC_SRCS GRPC_HDRS ${PROTO_SRC_DIR} ${PROTOS})

add_executable(pdfpushmeta 
analysis/pdfpushmeta.cpp 
putgetMeta/metaclient.cpp 
utils/ipTool.cpp
simulation/settings.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(pdfpushmeta adios2::adios2 MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf)
target_include_directories(pdfpushmeta PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

add_executable(pdfpushevent
analysis/pdfpushevent.cpp
putgetMeta/metaclient.cpp  
putEvent/eventclient.cpp 
utils/ipTool.cpp
simulation/settings.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(pdfpushevent adios2::adios2 MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf)
target_include_directories(pdfpushevent PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_executable(gray-scott
  simulation/main.cpp
  simulation/gray-scott.cpp
  simulation/settings.cpp
  simulation/writer.cpp
  putgetMeta/metaclient.cpp 
utils/ipTool.cpp
simulation/settings.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(gray-scott adios2::adios2 MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf)
target_include_directories(gray-scott PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


if (VTK)
find_package(VTK COMPONENTS
    FiltersCore
    IOImage
    IOXML
)

add_executable(anapullmeta 
analysis/anapullmeta.cpp 
putgetMeta/metaclient.cpp 
utils/ipTool.cpp
simulation/settings.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(anapullmeta adios2::adios2 MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf ${VTK_LIBRARIES})
target_include_directories(anapullmeta PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

add_executable(anastartbyevent 
analysis/anastartbyevent.cpp 
putgetMeta/metaclient.cpp 
utils/ipTool.cpp
simulation/settings.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(anastartbyevent adios2::adios2 MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf ${VTK_LIBRARIES})
target_include_directories(anastartbyevent PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

add_executable(simwithpdf
simulation/simwithpdf.cpp
simulation/gray-scott.cpp
simulation/settings.cpp
simulation/writer.cpp
putgetMeta/metaclient.cpp
utils/ipTool.cpp
${PROTO_SRCS} 
${GRPC_SRCS}) 

target_link_libraries(simwithpdf 
adios2::adios2 
MPI::MPI_C 
gRPC::grpc++_reflection 
protobuf::libprotobuf)
target_include_directories(simwithpdf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_executable(isosurface 
analysis/isosurface.cpp 
simulation/settings.cpp
putgetMeta/metaclient.cpp
utils/ipTool.cpp
${PROTO_SRCS} 
${GRPC_SRCS})

target_link_libraries(isosurface adios2::adios2 ${VTK_LIBRARIES} MPI::MPI_C gRPC::grpc++_reflection protobuf::libprotobuf)
target_include_directories(isosurface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

endif(VTK)
