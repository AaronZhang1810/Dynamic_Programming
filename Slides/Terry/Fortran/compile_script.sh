#######
#compile_script.sh
#
#This script compiles the base library and the main code
#######

FCFLAGS="-O3 -fopenmp -m64"
PROJNAME=firm_price_wrapper
COMP=gfortran

echo "Removing old files."

rm base_lib.o ${PROJNAME}.o

echo "Compiling the base module."

${COMP} ${FCFLAGS} -c -o base_lib.o base_lib.f90
${COMP} ${FCFLAGS} -c -o firm_price_solve.o firm_price_solve.f90

echo "Compiling the main program."

${COMP} ${FCFLAGS} -c -o ${PROJNAME}.o -I . ${PROJNAME}.f90

echo "Creating the executable."

${COMP} ${FCFLAGS} base_lib.o firm_price_solve.o ${PROJNAME}.o -o ${PROJNAME}.exe
