dohp_link_executable (fs-ex1 ex1.c)
dohp_link_executable (ellip ellip.c)
# dohp_link_executable (cellip cellip.c)
dohp_link_executable (bu bu.c)
dohp_link_executable (elast elast.c)

dohp_add_test (fs-mf-projection 1 fs-ex1 -snes_mf -ksp_type minres -const_BDeg 6 -snes_monitor -ksp_converged_reason)
dohp_add_test (fs-mf-op-proj-0  1 fs-ex1 -snes_mf_operator -ksp_type minres -pc_type jacobi -const_BDeg 10 -snes_monitor -ksp_converged_reason -ksp_rtol 1e-10 -frequency 5,4,3 -require_ptwise 2e-6,3e-6,2.3e-5 -require_grad 7e-7,3e-6,2e-4)
dohp_add_test (fs-mf-op-proj-1  1 fs-ex1 -const_BDeg 7 -snes_mf_operator -ksp_converged_reason -pc_type jacobi -ksp_type minres -ksp_rtol 1e-10 -snes_monitor -frequency 2,2,2 -require_ptwise 2e-6,2e-6,9e-6 -require_grad 3e-6,7e-6,1.1e-4)

dohp_add_test (bu-0 1 bu -dmesh_in "${Dohp_DATA_DIR}/dblock211.h5m")

dohp_add_test (ellip-e0-b2-p18 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor -ksp_monitor -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 2 -ellip_p 1.8)
dohp_add_test (ellip-e0-b4-p18 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor -ksp_monitor -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 4 -ellip_p 1.8 -djac_tensor_no_unroll)
dohp_add_test (ellip-e0-b4-p17-sse 1 ellip -dmesh_in "${Dohp_DATA_DIR}/dblock4.h5m" -snes_monitor -ksp_monitor -pc_type lu -exact_a 1.6 -exact_b 1.7 -exact_c 1.8 -exact 0 -const_bdeg 4 -ellip_p 1.7)
