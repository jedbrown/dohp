add_custom_command (OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/elastexact.h" "${CMAKE_CURRENT_BINARY_DIR}/elastexact.c"
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/elastexact.py" DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/elastexact.py" "${CMAKE_CURRENT_SOURCE_DIR}/dohpexact.py")
set_source_files_properties (elast.c PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_BINARY_DIR}")

add_custom_command (OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/stokesexact.c"
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/stokesexact.py" DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/stokesexact.py" "${CMAKE_CURRENT_SOURCE_DIR}/dohpexact.py")
set_source_files_properties (stokesexact.c PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}") # for stokesimpl.h

add_custom_command (OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/vhtexact.c"
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/vhtexact.py" DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/vhtexact.py" "${CMAKE_CURRENT_SOURCE_DIR}/dohpexact.py")
set_source_files_properties (vhtexact.c PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}") # for vhtimpl.h

dohp_link_executable (fs-ex1 ex1.c)
dohp_link_executable (fs-ex4 fsex4.c)
dohp_link_executable (ellip ellip.c)
dohp_link_executable (plapeig plapeig.c)
# dohp_link_executable (cellip cellip.c)
dohp_link_executable (bu bu.c)
dohp_link_executable (cunit cunit.c)
dohp_link_executable (elast elast.c elastexact.c)

set (STOKES_SOURCES stokes.c stokesexact.c)
set (VHT_SOURCES vht.c vhtexact.c)
if (dHAVE_GDAL) # GDAL is required in these cases
  set_source_files_properties (stokesjako.c PROPERTIES COMPILE_FLAGS "-I${GDAL_INCLUDE_DIR}")
  list (APPEND STOKES_SOURCES stokesjako.c)
  set_source_files_properties (vhtjako.c PROPERTIES COMPILE_FLAGS "-I${GDAL_INCLUDE_DIR}")
  list (APPEND VHT_SOURCES vhtjako.c)
endif ()

# dohp_link_executable does not currently work to link GDAL
add_executable (stokes ${STOKES_SOURCES})
target_link_libraries (stokes dohp ${GDAL_LIBRARY})
add_executable (vht ${VHT_SOURCES})
target_link_libraries (vht dohp ${GDAL_LIBRARY})

dohp_add_test (fs-mf-projection-v1 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 1)
#dohp_add_test (fs-mf-projection-v2 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 2)
dohp_add_test (fs-mf-projection-v3 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 3)
dohp_add_test (fs-mf-projection-v3-sparse 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 3 -proj_qmethod sparse)
dohp_add_test (fs-mf-op-proj-0  1 fs-ex1 -snes_mf_operator -ksp_type minres -pc_type jacobi -const_BDeg 9 -snes_monitor_short -ksp_converged_reason -ksp_rtol 1e-10 -frequency 5,4,3 -require_ptwise 2e-6,3e-6,2.3e-5 -require_grad 7e-7,3e-6,2e-4 -dfs_ordering_type natural)
dohp_add_test (fs-mf-op-proj-1  1 fs-ex1 -const_BDeg 6 -snes_mf_operator -ksp_converged_reason -pc_type jacobi -ksp_type minres -ksp_rtol 1e-10 -snes_monitor_short -frequency 2,2,2 -require_ptwise 2e-6,2e-6,9e-6 -require_grad 3e-6,7e-6,1.1e-4 -dfs_ordering_type natural)

dohp_add_test (cunit-identity-0 1 cunit -dfs_ordering_type natural -deform_type identity -expanded_view -cexp_view -nce_view -nce_compare)
dohp_add_test (cunit-affine-0 1 cunit -dfs_ordering_type natural -deform_type affine -nce_view -nce_compare)

dohp_add_test (bu-0 1 bu -dmesh_gen block -dmeshgen_block_mnp 3,2,2 -dfs_ordering_type natural)

dohp_add_test (fsex4-submesh-p1 1 fs-ex4 -bdeg 1 -dfs_view -dfs_ordering_type natural)
dohp_add_test (fsex4-submesh-p1-read 1 fs-ex4 -bdeg 1 -dfs_view -dfs_ordering_type natural -read_back)
dohp_add_test (fsex4-submesh-p2-read 1 fs-ex4 -bdeg 2 -dfs_view -dfs_ordering_type natural -read_back -read_back_vec)

dohp_add_test (ellip-e0-b1-p16 1 ellip -dmesh_gen block -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 1 -ellip_p 1.6 -dfs_ordering_type natural)
dohp_add_test (ellip-e0-b3-p16 1 ellip -dmesh_gen block -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 3 -ellip_p 1.6 -djac_tensor_no_unroll -dfs_ordering_type natural)
dohp_add_test (ellip-e0-b3-p14-sse 1 ellip -dmesh_gen block -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 3 -ellip_p 1.4)
dohp_add_test (ellip-e0-b3-p17-lam1-sse 1 ellip -dmesh_gen block -snes_monitor_short -ksp_monitor_short -pc_type lu -exact 0 -const_bdeg 3 -ellip_p 1.7 -ellip_lam 1)
dohp_add_test (ellip-morph 1 ellip -dmesh_gen block -morph 1 -twist 0.4 -ksp_monitor_short -snes_max_it 1 -ksp_rtol 1e-8 -ksp_type cg -pc_type lu)
dohp_add_test (ellip-morph-qfast 1 ellip -dmesh_gen block -morph 0.5 -twist 0.1 -ellip_p 1.9 -exact 3 -ksp_monitor_short -snes_max_it 1 -ksp_type richardson -pc_type lu -ellip_jac_qmethod fast)

dohp_add_test (elast-e0-b3-g01 1 elast -dmesh_gen block -snes_monitor_short -ksp_monitor_short -pc_type ilu -exact_a 0.98 -exact_b 1.01 -exact_c 1.01 -exact 0 -const_bdeg 3 -dfs_ordering_type natural)

dohp_add_test (stokes-q3q2-scr 1 stokes -dmesh_gen block -snes_max_it 1 -snes_monitor_short -ksp_type preonly
  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_real_diagonal -fieldsplit_p_ksp_monitor_short -fieldsplit_p_ksp_type cg
  -fieldsplit_u_ksp_converged_reason -stokes_Ap_mat_type sbaij -stokes_case Exact1 -const_bdeg 3 -pressure_codim 1)
dohp_add_test (stokes-q3q2-fgmres 1 stokes -dmesh_gen block -snes_monitor_short -snes_rtol 1e-7
  -ksp_monitor_short -ksp_type fgmres -ksp_rtol 1e-2
  -pc_type fieldsplit -pc_fieldsplit_real_diagonal -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type upper
  -fieldsplit_u_ksp_type cg -fieldsplit_u_ksp_converged_reason -fieldsplit_u_ksp_rtol 1e-2
  -fieldsplit_p_ksp_type preonly
  -stokes_case Exact1 -const_bdeg 3 -pressure_codim 1 -rheo_p 1.2 -rheo_eps 0.1)

dohp_add_test (vht-e1-q3q2q3-fgmres 1 vht -dmesh_gen block -snes_monitor_short -snes_converged_reason -ksp_converged_reason
  -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type additive -fieldsplit_e_pc_type lu
  -vht_case Exact1 -rheo_kappa0 1 -rheo_kappa1 0.1 -rheo_p 1.5 -rheo_eps 0.1 -rheo_T0 6 -rheo_Q 5 -rheo_R 0.33)
