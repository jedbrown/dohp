install_files (/include FILES
  dohpconfig.h
  # public headers
  dohptype.h                    # mandatory and included by all components below
  dohpfs.h
  dohpjacobi.h
  dohpmesh.h
  dohpvec.h
  dohpviewer.h
  dohpviewerdhm.h
  dohpstring.h                  # safer string handling, entirely optional
  dohpsys.h                     # needed for user 'main', not usually otherwise
  dohp.h                        # macros, optional for users, used heavily internally
  dohpgeom.h                    # inlined geometry primitives
  # private headers, we install them only so that users can extend Dohp
  dohpfsimpl.h
  dohpjacimpl.h
  dohpmeshimpl.h
  dohpunits.h
  dohputil.h
  )
