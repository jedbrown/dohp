# Should only be included when Dohp_BUILD_SANDBOX

set (Dohp_SANDBOX_EXECS
  dohpblock)

foreach (EXE ${Dohp_SANDBOX_EXECS})
  add_executable (${EXE} "${EXE}.c")
  #target_link_libraries (${EXE} ${Dohp_LIBRARIES})
  target_link_libraries (${EXE} dohp)
endforeach (EXE)

install (TARGETS ${Dohp_SANDBOX_EXECS}
  RUNTIME DESTINATION bin)
