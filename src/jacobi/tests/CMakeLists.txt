dohp_link_executable(jacobi-ex1 ex1.c)
dohp_add_test (jacobi-1 1 jacobi-ex1 -min_rdeg 10 -max_rdeg 11 -max_bdeg 10 -malloc_dump)
dohp_add_test (jacobi-1-lobatto 1 jacobi-ex1 -min_rdeg 9 -max_rdeg 11 -max_bdeg 10 -dquad_tensor_gauss_family lobatto -malloc_dump)

dohp_link_executable(jacobi-ex2 ex2.c)
dohp_add_test (jacobi-2-modal-fast-gauss 1 jacobi-ex2 -djac_type modal -dquad_tensor_gauss_family gauss -dquad_tensor_method fast -rule_degree 6 -malloc_dump)
dohp_add_test (jacobi-2-modal-fast-lobatto 1 jacobi-ex2 -djac_type modal -dquad_tensor_gauss_family lobatto -dquad_tensor_method fast -rule_degree 4 -malloc_dump)
dohp_add_test (jacobi-2-modal-sparse-gauss 1 jacobi-ex2 -djac_type modal -dquad_tensor_gauss_family gauss -dquad_tensor_method sparse -rule_degree 6 -malloc_dump)
