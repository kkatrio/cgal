# This is the CMake script for compiling a CGAL application.

project (Generator_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

  add_definitions(-DQT_NO_KEYWORDS)
  #--------------------------------
  # Demo: Generator_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui( DT_UI_FILES Generator_2.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Generator_2.qrc )

  # use the Qt MOC preprocessor on classes that derives from QObject
  qt5_generate_moc( "Generator_2.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Generator_2.moc" )

  # The executable itself.
  add_executable  ( Generator_2 Generator_2.cpp Generator_2.moc ${DT_UI_FILES}  ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Generator_2 )

  target_link_libraries( Generator_2 PRIVATE
    CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Generator_2)
else()

  message(STATUS "NOTICE: This demo requires CGAL, and Qt5, and will not be compiled.")

endif()
