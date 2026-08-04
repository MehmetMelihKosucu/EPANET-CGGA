# Unit tests for EPANET-CGGA.
#
# Add to the root CMakeLists.txt with:
#     enable_testing()
#     add_subdirectory(tests)

add_executable(test_flowhistory test_flowhistory.cpp)

target_include_directories(test_flowhistory PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Core
    ${CMAKE_SOURCE_DIR}/src/Elements
)

# Link against the same library target the CLI uses.
target_link_libraries(test_flowhistory epanet3)

add_test(NAME flow_history COMMAND test_flowhistory)
