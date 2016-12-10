include Make.include


all:
	cd TEMPLATE;      ${MAKE} TEMP
	cd API;           ${MAKE} API
	cd BLAS1;         ${MAKE} BLAS1
	cd BLAS2;         ${MAKE} BLAS2
	cd BLAS3;         ${MAKE} BLAS3
	cd COPY;          ${MAKE} COPY
	cd Local_BLAS;    ${MAKE} Local_BLAS
	cd OBJ;           ${MAKE} OBJ
	cd REDUCE;        ${MAKE} REDUCE
	cd UTIL;          ${MAKE} UTIL
	cd IO;            ${MAKE} IO
	cd SOLVERS;       ${MAKE} SOLVERS
	cd POOCLAPACK;     ${MAKE} POOCLAPACK
	cd FORTRAN_interface;     ${MAKE} FORTRAN_interface

removeall: clean
	rm -f *.a

clean: 
	cd TEMPLATE;     ${MAKE} clean
	cd API;          ${MAKE} clean
	cd BLAS1;        ${MAKE} clean
	cd BLAS2;        ${MAKE} clean
	cd BLAS3;        ${MAKE} clean
	cd Local_BLAS;   ${MAKE} clean
	cd COPY;         ${MAKE} clean
	cd OBJ;          ${MAKE} clean
	cd REDUCE;       ${MAKE} clean
	cd UTIL;         ${MAKE} clean
	cd SOLVERS;      ${MAKE} clean
	cd EXAMPLES;     ${MAKE} clean
	cd IO;           ${MAKE} clean
	cd TEST;         ${MAKE} clean
	cd FORTRAN_interface  ;  ${MAKE} clean
	cd POOCLAPACK;    ${MAKE} clean
	rm -f *~ core




