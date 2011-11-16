project (VisItReaderDohp)

cmake_minimum_required (VERSION 2.8)

find_package (ParaView REQUIRED)
include (${PARAVIEW_USE_FILE})

set (AVTDOHP_SOURCES
  avtDohpFileFormat.C
  )
include_directories(/home/jed/dohp/bompi/include /home/jed/dohp/include /home/jed/petsc/ompi/include /home/jed/petsc/include /home/jed/usr/itaps-ompi/include)

add_visit_plugin_reader (VisItReaderDohp "1.0"
  VISIT_READER_NAME "avtDohpFileFormat"
  VISIT_READER_TYPE "MTMD"
  VISIT_READER_FILE_PATTERN "dhm"
  SERVER_SOURCES ${AVTDOHP_SOURCES}
  )

target_link_libraries (VisItReaderDohp vtkVisItAVTAlgorithms /home/jed/dohp/bompi/lib/libdohp.so)
