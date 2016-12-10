include ../Make.include

SOLVERS: 
	cd LU; make LU
	cd Invert; make Inv
	cd Cholesky; make Cholesky
	cd QR; make QR
	cd Reduction; make Red
	cd Eigensolvers; make Eig
	cd SVD; make SVD

clean:
	cd LU; make clean
	cd Invert; make clean
	cd Cholesky; make clean
	cd QR; make clean
	cd Reduction; make clean
	cd Eigensolvers; make clean
	cd SVD; make clean
	cd LINPACK; make clean
	rm -f $(solvers) SOLVERS *.o *~

