# une bibliotheque statique (libmnblas.a)
# une bibliotheque dynamique (libmnblasdyn.so)

all: libmnblas.a libmnblasdyn.so

LIB_DIR=../lib
INC_DIR=../include

OPTIONS=-Wall -fPIC $(OPTIM) -I$(INC_DIR)

# Bibliothèque statique

libmnblas.a: copy.o swap.o dot.o iamin.o iamax.o axpy.o nrm2.o asum.o complexe.o gemm.o aux_blas3.o absolute_value.o gemv.o
	rm -f libmnblas.a $(LIB_DIR)/libmnblas.a
	ar -r libmnblas.a copy.o swap.o dot.o iamin.o iamax.o axpy.o nrm2.o asum.o complexe.o gemm.o aux_blas3.o absolute_value.o gemv.o
	cp libmnblas.a $(LIB_DIR)

# Bibliothèque dynamique

libmnblasdyn.so:copy.o swap.o dot.o iamin.o iamax.o axpy.o nrm2.o asum.o complexe.o gemm.o aux_blas3.o absolute_value.o gemv.o
	rm -f libmnblasdyn.so $(LIB_DIR)/libmnblasdyn.so
	gcc -shared -o libmnblasdyn.so copy.o swap.o dot.o iamin.o iamax.o axpy.o nrm2.o asum.o complexe.o gemm.o aux_blas3.o absolute_value.o gemv.o -lm
	cp libmnblasdyn.so $(LIB_DIR)

# Fichiers auxiliaires
complexe.o: $(INC_DIR)/complexe.h complexe.c
	gcc $(OPTIONS) -c complexe.c

absolute_value.o: $(INC_DIR)/absolute_value.h absolute_value.c
	gcc $(OPTIONS) -c absolute_value.c

aux_blas3.o: $(INC_DIR)/aux_blas3.h aux_blas3.c
	gcc $(OPTIONS) -c aux_blas3.c

# BLAS 1

copy.o : copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c copy.c

swap.o : swap.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c swap.c

dot.o : dot.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c dot.c

iamin.o : iamin.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c iamin.c

iamax.o : iamax.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c iamax.c

axpy.o : axpy.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c axpy.c

nrm2.o : nrm2.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c nrm2.c -lm

asum.o : asum.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c asum.c

# BLAS 2

gemv.o: gemv.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemv.c

# BLAS 3

gemm.o: gemm.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemm.c

clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a
