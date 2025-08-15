#--------------------------------------------------------------------------------#
# MASA analyzer                                                                  #
# ----------------------------------                                             #
# (C) 2024 Nokia Technologies Ltd. See LICENSE.md for license.                  #
#                                                                                #
#--------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.8.2)
project(MasaAnalyzerProject C)
set(ANALYZER_NAME masaAnalyzer)

set(CMAKE_C_STANDARD 11)

if(WIN32)
    add_definitions(
        -D_CRT_SECURE_NO_WARNINGS
        -Drestrict=__restrict
        /wd4028
        /wd4100
        /wd4244
        /wd4310
        /wd4996
        /MP
    )
endif()

set(ANALYZER_SOURCES
    masaAnalyzer.c
    masaCommonTables.c
    sphericalIndexing.c
    cldfb.c
    cldfb_prototype_filters.c
    pocketfft/pocketfft.c)

set(ANALYZER_HEADERS
    masaCommonHeader.h
    sphericalIndexing.h
    cldfb.h
    cldfb_prototype_filters.h
    pocketfft/pocketfft.h)


add_executable(${ANALYZER_NAME} ${ANALYZER_SOURCES} ${ANALYZER_HEADERS})

if(NOT WIN32)
    target_link_libraries(${ANALYZER_NAME} m)
endif()


add_custom_command(
    TARGET ${ANALYZER_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/eigen_to_foa_cldfb_domain_filters.bin
    $<TARGET_FILE_DIR:${ANALYZER_NAME}>/eigen_to_foa_cldfb_domain_filters.bin)

add_custom_command(
    TARGET ${ANALYZER_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/eigen_to_hoa2_cldfb_domain_filters.bin
    $<TARGET_FILE_DIR:${ANALYZER_NAME}>/eigen_to_hoa2_cldfb_domain_filters.bin)

add_custom_command(
    TARGET ${ANALYZER_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/sector_filters.bin
    $<TARGET_FILE_DIR:${ANALYZER_NAME}>/sector_filters.bin)

