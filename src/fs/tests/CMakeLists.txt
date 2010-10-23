dohp_link_executable (fs-ex1 ex1.c)
dohp_link_executable (fs-ex4 fsex4.c)
dohp_link_executable (ellip ellip.c)
# dohp_link_executable (cellip cellip.c)
dohp_link_executable (bu bu.c)
dohp_link_executable (cunit cunit.c)
dohp_link_executable (elast elast.c)
dohp_link_executable (stokes stokes.c)

dohp_add_test (fs-mf-projection 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 6 -snes_monitor -ksp_converged_reason -dfs_ordering_type natural -malloc_dump)
dohp_add_test (fs-mf-op-proj-0  1 fs-ex1 -snes_mf_operator -ksp_type minres -pc_type jacobi -const_BDeg 10 -snes_monitor -ksp_converged_reason -ksp_rtol 1e-10 -frequency 5,4,3 -require_ptwise 2e-6,3e-6,2.3e-5 -require_grad 7e-7,3e-6,2e-4 -dfs_ordering_type natural -malloc_dump)
dohp_add_test (fs-mf-op-proj-1  1 fs-ex1 -const_BDeg 7 -snes_mf_operator -ksp_converged_reason -pc_type jacobi -ksp_type minres -ksp_rtol 1e-10 -snes_monitor -frequency 2,2,2 -require_ptwise 2e-6,2e-6,9e-6 -require_grad 3e-6,7e-6,1.1e-4 -dfs_ordering_type natural -malloc_dump)

dohp_add_test (cunit-identity-0 1 cunit -dfs_ordering_type natural -deform_type identity -expanded_view -cexp_view -nce_view -nce_compare)
dohp_add_test (cunit-affine-0 1 cunit -dfs_ordering_type natural -deform_type affine -nce_view -nce_compare)

dohp_add_test (bu-0 1 bu -dmesh_in "${Dohp_DATA_DIR}/dblock211.h5m" -dfs_ordering_type natural -malloc_dump)

dohp_add_test (fsex4-submesh-p2 1 fs-ex4 -bdeg 2 -dfs_view -dfs_ordering_type natural -malloc_dump)
dohp_add_test (fsex4-submesh-p2-read 1 fs-ex4 -bdeg 2 -dfs_view -dfs_ordering_type natural -read_back -malloc_dump)
dohp_add_test (fsex4-submesh-p3-read 1 fs-ex4 -bdeg 3 -dfs_view -dfs_ordering_type natural -read_back -read_back_vec -malloc_dump)

dohp_add_test (ellip-e0-b2-p16 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 2 -ellip_p 1.6 -dfs_ordering_type natural -malloc_dump)
dohp_add_test (ellip-e0-b4-p16 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 4 -ellip_p 1.6 -djac_tensor_no_unroll -dfs_ordering_type natural -malloc_dump)
dohp_add_test (ellip-e0-b4-p14-sse 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 4 -ellip_p 1.4 -malloc_dump)
dohp_add_test (ellip-e0-b4-p17-lam1-sse 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type lu -exact 0 -const_bdeg 4 -ellip_p 1.7 -ellip_lam 1 -malloc_dump)
dohp_add_test (ellip-morph 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -morph 1 -twist 0.4 -ksp_monitor_short -snes_max_it 1 -ksp_rtol 1e-8 -ksp_type cg -pc_type lu -malloc_dump)
dohp_add_test (elast-e0-b4-g01 1 elast -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor_short -ksp_monitor_short -pc_type ilu -exact_a 0.98 -exact_b 1.01 -exact_c 1.01 -exact 0 -const_bdeg 4 -elast_gamma 0.1 -malloc_dump)
