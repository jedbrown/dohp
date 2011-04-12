# Should only be included when Dohp_BUILD_TOOLS

set (Dohp_TOOLS_EXECS
  dohpblock
  dhmdump
  )

foreach (EXE ${Dohp_TOOLS_EXECS})
  add_executable (${EXE} "${EXE}.c")
  target_link_libraries (${EXE} dohp)
endforeach (EXE)

install (TARGETS ${Dohp_TOOLS_EXECS}
  RUNTIME DESTINATION bin)
