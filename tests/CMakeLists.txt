# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

set(TEST_SOURCES_UTILS
    test_utils.c
    unity/unity.c
    ../src/utils.c
    ../src/pvalue.c
)

# set(TEST_SOURCES_SPLIT
#     test_split.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_SPLIT_WIDE
#     test_split_wide.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_QUASI_SPLIT
#     test_quasi_split.c
#     ../src/split.c
#     ../src/utils.c
#     ../src/gee.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_SPLIT_ASYM
#     test_split_asym.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_SPLIT_GEE
#     test_split_gee.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_QUASI_PHI
#     test_phi_est.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_TREE
#     test_tree.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_TREE_WIDE
#     test_tree_wide.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_FOREST
#     test_forest.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/forest.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_FOREST_WIDE
#     test_forest_wide.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/forest.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_FOREST_ASYM
#     test_forest_asym.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/forest.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_FOREST_GEE
#     test_forest_gee.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/utils.c
#     ../src/tree.c
#     ../src/forest.c
#     ../src/surv.c
#     ../src/pvalue.c
# )

# set(TEST_SOURCES_FUNCTIONPOINTER
#     test_function_pointer.c
# )

# set(TEST_SURV

#     test_surv.c
#     ../src/surv.c
#     ../src/utils.c
# )

# set(TEST_DYNAMIC

#     test_dynamic.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/pvalue.c
#     ../src/gee.c)

# set(TEST_TREE_DYNAMIC

#     test_tree_dynamic.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/pvalue.c
#     ../src/tree.c
#     ../src/forest.c)

# set(TEST_TREE_ASYM

#     test_tree_asym.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/pvalue.c
#     ../src/gee.c
#     ../src/tree.c
#     ../src/forest.c)

# set(TEST_TREE_GEE

#     test_tree_gee.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/pvalue.c
#     ../src/gee.c
#     ../src/tree.c
#     ../src/forest.c)

# set(TEST_PATH

#     ./test_paths_print.c)

# set(TEST_TREE_PATH

#     ./test_tree_path.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/pvalue.c
#     ../src/tree.c)

# set(TEST_GLOBAL_VARIABLES

#     ./test_global_variable/helper
#     ./test_global_variable/main.c)

# set(TEST_TREE_BOOTSTRAP

#     test_tree_wt_bootstrap.c
#     ../src/utils.c
#     ../src/surv.c
#     ../src/split.c
#     ../src/gee.c
#     ../src/pvalue.c
#     ../src/tree.c
#     ../src/forest.c)

# set(TEST_GEE

#     test_gee_pure.c
#     ../src/utils.c
#     ../src/gee.c
#     ../src/pvalue.c
#     ../src/surv.c
#     ../src/split.c)

add_executable(test_utils ${TEST_SOURCES_UTILS})
# add_executable(test_split ${TEST_SOURCES_SPLIT})
# add_executable(test_split_wide ${TEST_SOURCES_SPLIT_WIDE})
# add_executable(test_quasi_split ${TEST_SOURCES_QUASI_SPLIT})
# add_executable(test_quasi_phi ${TEST_SOURCES_QUASI_PHI})
# add_executable(test_split_asym ${TEST_SOURCES_SPLIT_ASYM})
# add_executable(test_split_gee ${TEST_SOURCES_SPLIT_GEE})
# add_executable(test_tree ${TEST_SOURCES_TREE})
# add_executable(test_tree_wide ${TEST_SOURCES_TREE_WIDE})
# add_executable(test_forest ${TEST_SOURCES_FOREST})
# add_executable(test_forest_wide ${TEST_SOURCES_FOREST_WIDE})
# add_executable(test_forest_asym ${TEST_SOURCES_FOREST_ASYM})
# add_executable(test_forest_gee ${TEST_SOURCES_FOREST_GEE})
# add_executable(test_function_pointer ${TEST_SOURCES_FUNCTIONPOINTER})
# add_executable(test_surv ${TEST_SURV})
# add_executable(test_dynamic ${TEST_DYNAMIC})
# add_executable(test_paths_print ${TEST_PATH})
# add_executable(test_tree_dynamic ${TEST_TREE_DYNAMIC})
# add_executable(test_tree_asym ${TEST_TREE_ASYM})
# add_executable(test_tree_gee ${TEST_TREE_GEE})
# add_executable(test_tree_path ${TEST_TREE_PATH})
# add_executable(test_global_variables ${TEST_GLOBAL_VARIABLES})
# add_executable(test_tree_bootstrap ${TEST_TREE_BOOTSTRAP})
# add_executable(test_gee ${TEST_GEE})

target_link_libraries(test_utils m gsl gslcblas)
# target_link_libraries(test_split m gsl gslcblas)
# target_link_libraries(test_split_wide m gsl gslcblas)
# target_link_libraries(test_quasi_split m gsl gslcblas)
# target_link_libraries(test_quasi_phi m gsl gslcblas)
# target_link_libraries(test_split_asym m gsl gslcblas)
# target_link_libraries(test_split_gee m gsl gslcblas)
# target_link_libraries(test_tree m gsl gslcblas)
# target_link_libraries(test_tree_wide m gsl gslcblas)
# target_link_libraries(test_forest PRIVATE OpenMP::OpenMP_C m gsl gslcblas)
# target_link_libraries(test_forest_wide PRIVATE OpenMP::OpenMP_C m gsl gslcblas)
# target_link_libraries(test_forest_asym PRIVATE OpenMP::OpenMP_C m gsl gslcblas)
# target_link_libraries(test_forest_gee PRIVATE OpenMP::OpenMP_C m gsl gslcblas)
# target_link_libraries(test_function_pointer)
# target_link_libraries(test_surv m gsl gslcblas)
# target_link_libraries(test_dynamic m gsl gslcblas)
# target_link_libraries(test_tree_dynamic m gsl gslcblas)
# target_link_libraries(test_tree_asym m gsl gslcblas)
# target_link_libraries(test_tree_gee m gsl gslcblas)
# target_link_libraries(test_tree_path m gsl gslcblas)
# target_link_libraries(test_tree_bootstrap m gsl gslcblas)
# target_link_libraries(test_gee m gsl gslcblas)
# target_link_libraries(test_global_variables PRIVATE OpenMP::OpenMP_C)

# add_test(NAME MyTests COMMAND test_utils)
add_custom_target(TestUtils
    COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_utils
    DEPENDS test_utils
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

# add_custom_target(TestSplit
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_split
#     DEPENDS test_split
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestSplitWide
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_split_wide
#     DEPENDS test_split_wide
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestQuasiSplit
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_quasi_split
#     DEPENDS test_quasi_split
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestQuasiPhiEst
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_quasi_phi
#     DEPENDS test_quasi_phi
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestSplitAsym
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_split_asym
#     DEPENDS test_split_asym
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestSplitGee
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_split_gee
#     DEPENDS test_split_gee
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestTree
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree
#     DEPENDS test_tree
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestTreeWide
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_wide
#     DEPENDS test_tree_wide
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestForest
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_forest
#     DEPENDS test_forest
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestForestWide
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_forest_wide
#     DEPENDS test_forest_wide
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestForestAsym
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_forest_asym
#     DEPENDS test_forest_asym
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestForestGee
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_forest_gee
#     DEPENDS test_forest_gee
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestFunctionPointer
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_function_pointer
#     DEPENDS test_function_pointer
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestSurv
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_surv
#     DEPENDS test_surv
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestDynamic
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_dynamic
#     DEPENDS test_dynamic
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestTreeDynamic
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_dynamic
#     DEPENDS test_tree_dynamic
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
# )

# add_custom_target(TestTreeAsym
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_asym
#     DEPENDS test_tree_asym
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestTreeGee
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_gee
#     DEPENDS test_tree_gee
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestTreePath
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_path
#     DEPENDS test_tree_path
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestGlobalVariables
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_global_variables
#     DEPENDS test_global_variables
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestTreeBootstrap
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_tree_bootstrap
#     DEPENDS test_tree_bootstrap
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# add_custom_target(TestGee
#     COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes --log-file=vg_logfile_%p.out ./test_gee
#     DEPENDS test_gee
#     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})

# ## add the test for all
# add_custom_target(TestAll
#     DEPENDS TestUtils TestSplit TestSplitWide TestQuasiSplit TestQuasiPhiEst TestSplitAsym TestSplitGee TestTree TestTreeWide TestForestWide TestForestAsym TestForestGee TestFunctionPointer TestSurv TestDynamic TestTreeDynamic TestTreeAsym TestTreeGee TestTreePath TestGlobalVariables TestTreeBootstrap
# )