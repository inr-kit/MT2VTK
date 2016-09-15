# SET( TESTING_ACTIVATED "true" )

SET( PROGRAM_SOURCES main.cxx  MeshTallyReader.cxx  McnpTally.cxx )
SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -g )

ADD_EXECUTABLE( mt2vtk ${PROGRAM_SOURCES} )

############################
# compile and link testing #
############################
if(TESTING_ACTIVATED)
	SET( TEST_SOURCES test.cpp  MeshTallyReader.cxx  McnpTally.cxx )
	ADD_EXECUTABLE( test ${TEST_SOURCES} )
	TARGET_LINK_LIBRARIES ( test UnitTest++ )
endif()