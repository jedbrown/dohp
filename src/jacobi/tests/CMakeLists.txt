dohp_link_executable(jacobi-ex1 ex1.c)
add_test (jacobi1 jacobi-ex1 -min_rdeg 11 -max_rdeg 12 -max_bdeg 11)
