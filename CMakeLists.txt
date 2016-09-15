cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project (MeshTallyReader)

INCLUDE_DIRECTORIES ( inc )


# looking for boost
FIND_PACKAGE(Boost)
IF(Boost_FOUND)
    INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIR})
ELSE()
    MESSAGE( FATAL ERROR " --  Boost not found. Run cmake again with parameter -DBOOST_INCLUDE_DIR:FILEPATH=/path/to/your/includes ...(/boost) <- stop before boost ^_^")
ENDIF()

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin )

ADD_SUBDIRECTORY(src)


