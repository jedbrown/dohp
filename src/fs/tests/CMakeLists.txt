add_custom_command (OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/elastexact.h" "${CMAKE_CURRENT_BINARY_DIR}/elastexact.c"
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/elastexact.py"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/elastexact.py"
  )
set_source_files_properties (elast.c PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_BINARY_DIR}")

dohp_link_executable (fs-ex1 ex1.c)
dohp_link_executable (fs-ex4 fsex4.c)
dohp_link_executable (ellip ellip.c)
# dohp_link_executable (cellip cellip.c)
dohp_link_executable (bu bu.c)
dohp_link_executable (cunit cunit.c)
dohp_link_executable (elast elast.c elastexact.c)
dohp_link_executable (stokes stokes.c)

dohp_add_test (fs-mf-projection-v1 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 1)
#dohp_add_test (fs-mf-projection-v2 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 2)
dohp_add_test (fs-mf-projection-v3 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 3)
dohp_add_test (fs-mf-projection-v3-sparse 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 5 -snes_monitor_short -ksp_converged_reason -dfs_ordering_type natural -proj_version 3 -proj_qmethod sparse)
dohp_add_test (fs-mf-op-proj-0  1 fs-ex1 -snes_mf_operator -ksp_type minres -pc_type jacobi -const_BDeg 9 -snes_monitor_short -ksp_converged_reason -ksp_rtol 1e-10 -frequency 5,4,3 -require_ptwise 2e-6,3e-6,2.3e-5 -require_grad 7e-7,3e-6,2e-4 -dfs_ordering_type natural)
dohp_add_test (fs-mf-op-proj-1  1 fs-ex1 -const_BDeg 6 -snes_mf_operator -ksp_converged_reason -pc_type jacobi -ksp_type minres -ksp_rtol 1e-10 -snes_monitor_short -frequency 2,2,2 -require_ptwise 2e-6,2e-6,9e-6 -require_grad 3e-6,7e-6,1.1e-4 -dfs_ordering_type natural)

dohp_add_test (cunit-identity-0 1 cunit -dfs_ordering_type natural -deform_type identity -expanded_view -cexp_view -nce_view -nce_compare)
dohp_add_test (cunit-affine-0 1 cunit -dfs_ordering_type natural -deform_type affine -nce_view -nce_compare)

dohp_add_test (bu-0 1 bu -dmesh_in "${Dohp_DATA_DIR}/dblock211.h5m" -dfs_ordering_type natural)

dohp_add_test (fsex4-submesh-p1 1 fs-ex4 -bdeg 1 -dfs_view -dfs_ordering_type natural)
dohp_add_test (fsex4-submesh-p1-read 1 fs-ex4 -bdeg 1 -dfs_view -dfs_ordering_type natural -read_back)
dohp_add_test (fsex4-submesh-p2-read 1 fs-ex4 -bdeg 2 -dfs_view -dfs_ordering_type natural -read_back -read_back_vec)

dohp_add_test (ellip-e0-b1-p16 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 1 -ellip_p 1.6 -dfs_ordering_type natural)
dohp_add_test (ellip-e0-b3-p16 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 3 -ellip_p 1.6 -djac_tensor_no_unroll -dfs_ordering_type natural)
dohp_add_test (ellip-e0-b3-p14-sse 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 3 -ellip_p 1.4)
dohp_add_test (ellip-e0-b3-p17-lam1-sse 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact 0 -const_bdeg 3 -ellip_p 1.7 -ellip_lam 1)
dohp_add_test (ellip-morph 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -morph 1 -twist 0.4 -ksp_monitor_short -snes_max_it 1 -ksp_rtol 1e-8 -ksp_type cg -pc_type lu)
dohp_add_test (elast-e0-b3-g01 1 elast -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type ilu -exact_a 0.98 -exact_b 1.01 -exact_c 1.01 -exact 0 -const_bdeg 3 -elast_scale 0.1 -dfs_ordering_type natural)
