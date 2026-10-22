# Tests for EPANET-CGGA.
#
# Add to the root CMakeLists.txt with:
#     enable_testing()
#     add_subdirectory(tests)
#
# Labels
#   unit        unit tests of the flow history subsystem (fast)
#   regression  complete simulations compared with stored reference results
#   long        the 24 h Hadimkoy simulation; excluded from the default CI run
#
# Examples
#   ctest --test-dir build -C Release -LE long        (everything except long)
#   ctest --test-dir build -C Release -L long         (only the long test)

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
add_executable(test_flowhistory test_flowhistory.cpp)

target_include_directories(test_flowhistory PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Core
    ${CMAKE_SOURCE_DIR}/src/Elements
)

# Link against the same library target the CLI uses.
target_link_libraries(test_flowhistory epanet3)

add_test(NAME flow_history COMMAND test_flowhistory)
set_tests_properties(flow_history PROPERTIES LABELS "unit")

# ---------------------------------------------------------------------------
# Regression tests
# ---------------------------------------------------------------------------
add_executable(regression_test regression/regression_test.cpp)
target_include_directories(regression_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(regression_test epanet3)

set(REG_REF_DIR ${CMAKE_SOURCE_DIR}/tests/regression/reference)
set(REG_OUT_DIR ${CMAKE_BINARY_DIR}/regression_results)
file(MAKE_DIRECTORY ${REG_OUT_DIR})
file(MAKE_DIRECTORY ${REG_REF_DIR})

# add_regression_case(<name> <input file> <comma-separated node ids> <timeout s> <labels>)
set(REG_UPDATE_COMMANDS "")
function(add_regression_case name inp nodes timeout labels)
    add_test(NAME regression_${name}
             COMMAND regression_test
                     --inp   ${inp}
                     --ref   ${REG_REF_DIR}/${name}.csv
                     --out   ${REG_OUT_DIR}/${name}.csv
                     --nodes ${nodes}
             WORKING_DIRECTORY ${REG_OUT_DIR})
    set_tests_properties(regression_${name} PROPERTIES
                         TIMEOUT ${timeout}
                         LABELS "${labels}")
    set(REG_UPDATE_COMMANDS ${REG_UPDATE_COMMANDS}
        COMMAND $<TARGET_FILE:regression_test>
                --inp   ${inp}
                --ref   ${REG_REF_DIR}/${name}.csv
                --out   ${REG_OUT_DIR}/${name}.csv
                --nodes ${nodes}
                --update
        PARENT_SCOPE)
endfunction()

add_regression_case(onizuka_d10 ${CMAKE_SOURCE_DIR}/Networks/Onizuka1986-EPA3-D10.inp "J1,J2,T1" 1800 "regression")
add_regression_case(onizuka_d1  ${CMAKE_SOURCE_DIR}/Networks/Onizuka1986-EPA3.inp     "J1,J2,T1" 1800 "regression")
add_regression_case(nault2016   ${CMAKE_SOURCE_DIR}/Networks/EPA3-Nault2016.inp        "J16,J24"  1800 "regression")
add_regression_case(hadimkoy    ${CMAKE_SOURCE_DIR}/Networks/EPA3-hk-burst.inp         "612,3234" 14400 "regression;long")

# Regenerates every reference file from the current build:
#     cmake --build build --config Release --target regression_update
add_custom_target(regression_update
    ${REG_UPDATE_COMMANDS}
    WORKING_DIRECTORY ${REG_OUT_DIR}
    COMMENT "Regenerating regression reference files"
    VERBATIM)
add_dependencies(regression_update regression_test)
