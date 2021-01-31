gfortran -c main_1dheatconduction.f90 -llapack -lblas
gfortran -o m_1dhc main_1dheatconduction.o -llapack -lblas
./m_1dhc
