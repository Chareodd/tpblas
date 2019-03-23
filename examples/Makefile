all: test_dot test_dot_dyn test_complexe test_complexe_dyn test_copy test_copy_dyn test_asum test_asum_dyn test_gemm test_gemm_dyn test_gemv test_gemv_dyn

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-O2 -Wall -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  = -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN  = -L$(LIB_DIR) $(LIBDYN)

#PROBLEME DE LINK ET AUTRE !
#test_iamin:  test_iamin.o flop.o vector.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_iamin test_iamin.o flop.o vector.o $(OPTIONS_LINK_STATIC)
#
#test_iamin_dyn: test_iamin.o flop.o matrix.o
#	gcc -o test_iamin_dyn test_iamin.o flop.o matrix.o $(OPTIONS_LINK_DYN)
#
#test_iamin.o: test_iamin.c $(INC_DIR)/mnblas.h
#	gcc $(OPTIONS_COMPIL) -c test_iamin.c
#
#test_iamax:  test_iamax.o flop.o vector.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_iamax test_iamax.o flop.o vector.o $(OPTIONS_LINK_STATIC)
#
#test_iamax_dyn: test_iamax.o flop.o matrix.o
#	gcc -o test_iamax_dyn test_iamax.o flop.o matrix.o $(OPTIONS_LINK_DYN)
#
#test_iamax.o: test_iamax.c $(INC_DIR)/mnblas.h
#	gcc $(OPTIONS_COMPIL) -c test_iamax.c
#
#test_nrm2:  test_nrm2.o flop.o vector.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_nrm2 test_nrm2.o flop.o vector.o $(OPTIONS_LINK_STATIC)
#
#test_nrm2_dyn: test_nrm2.o flop.o matrix.o
#	gcc -o test_nrm2_dyn test_nrm2.o flop.o matrix.o $(OPTIONS_LINK_DYN)
#
#test_nrm2.o: test_nrm2.c $(INC_DIR)/mnblas.h
#	gcc $(OPTIONS_COMPIL) -c test_nrm2.c

test_gemv:  test_gemv.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemv test_gemv.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_gemv_dyn: test_gemv.o flop.o vec_mat.o
	gcc -o test_gemv_dyn test_gemv.o flop.o vec_mat.o $(OPTIONS_LINK_DYN)

test_gemv.o: test_gemv.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemv.c


test_gemm:  test_gemm.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemm test_gemm.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_gemm_dyn: test_gemm.o flop.o vec_mat.o
	gcc -o test_gemm_dyn test_gemm.o flop.o vec_mat.o $(OPTIONS_LINK_DYN)

test_gemm.o: test_gemm.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemm.c

test_asum: test_asum.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_asum test_asum.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_asum_dyn: test_asum.o flop.o vec_mat.o
	gcc -o test_asum_dyn test_asum.o flop.o vec_mat.o $(OPTIONS_LINK_DYN)

test_asum.o: test_asum.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_asum.c

test_complexe: test_complexe.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_complexe test_complexe.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_complexe_dyn: test_complexe.o flop.o vec_mat.o
	gcc -o test_complexe_dyn test_complexe.o flop.o vec_mat.o $(OPTIONS_LINK_DYN)

test_complexe.o: test_complexe.c
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_copy: test_copy.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_copy test_copy.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_copy_dyn: test_copy.o flop.o vec_mat.o
	gcc -o test_copy_dyn test_copy.o flop.o vec_mat.o $(OPTIONS_LINK_DYN)

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_dot: test_dot.o flop.o vec_mat.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o vec_mat.o $(OPTIONS_LINK_STATIC)

test_dot_dyn: test_dot.o flop.o vec_mat.o
	gcc -o test_dot_dyn flop.o test_dot.o vec_mat.o $(OPTIONS_LINK_DYN)

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_dot.c

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

vec_mat.o: vec_mat.c $(INC_DIR)/vec_mat.h
	gcc $(OPTIONS_COMPIL) -c vec_mat.c	

clean:
	rm -f *.o test_dot test_dot_dyn test_complexe test_complexe_dyn test_copy test_copy_dyn test_asum test_asum_dyn test_gemm test_gemm_dyn test_gemv test_gemv_dyn *~
