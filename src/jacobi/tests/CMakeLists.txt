dohp_link_executable(jacobi-ex1 ex1.c)
dohp_add_test (jacobi-1 1 jacobi-ex1 -min_rdeg 11 -max_rdeg 12 -max_bdeg 11 -malloc_dump)
