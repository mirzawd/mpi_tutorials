#------------------------------------------------------------
# CMakeLists for CortexNemac
#------------------------------------------------------------
file(GLOB SOURCE_CPP *.cpp)

get_filename_component(exName ${CMAKE_CURRENT_LIST_DIR} NAME)

set(exName "${exName}")

message(STATUS "exName:  ${exName} ")

add_executable(${exName} ${SOURCE_CPP})

set_target_properties( ${exName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${APPS_INSTALL_DIR})

install(TARGETS ${exName} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)                                        
