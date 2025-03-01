COMPILFORTAN  = gfortran
COMPILCPP     = gcc

NETCDFPATH    = /usr
HDF5PATH      = /usr
LAPACKPATH    = /usr

OPTIFLAGSFORTRAN   =  -O2
OPTIFLAGSCPP  =  -O2

FORTANFLAGS   =  -cpp -ffree-line-length-none -fconvert=big-endian $(OPTIFLAGSFORTRAN)
CPPFLAGS =  $(OPTIFLAGSCPP)

INCS          = -Isrc -I. -I$(NETCDFPATH)/include
LIBSFORTRAN   = -L$(LAPACKPATH)/lib -llapack -lblas -L$(NETCDFPATH)/lib -lnetcdf -lnetcdff
LIBSCPP       = -L$(NETCDFPATH)/lib -lnetcdf_c++ -lnetcdf -lnetcdff

ifeq ($(HOSTNAME),lemaitre2.cism.ucl.ac.be)
        COMPILFORTAN  = ifort
        COMPILCPP     = icc

        NETCDFPATH    = /usr/local/netcdf/4.1.3/intel
        LAPACKPATH    = /usr/local/lapack/3.2.1/gcc
        BLASPATH      = /usr/local/blas/3.2.1/gcc
#        CURLPATH      = /home/astr/sallaz/opt/curl/last

        OPTIFLAGSFORTRAN   = -O2 -W1 -warn all
#        OPTIFLAGSFORTRAN   =  -O0 -g -W1 -warn all -static -CB
        OPTIFLAGSCPP  = -O2
#        OPTIFLAGSCPP  =  -O0 -g -static

        FORTANFLAGS   =  -cpp -convert big_endian -assume byterecl -align dcommon -extend_source -Dlinux86  $(OPTIFLAGSFORTRAN)
        CPPFLAGS =  $(OPTIFLAGSCPP)

        INCS          = -Isrc -I. -I$(NETCDFPATH)/include
        LIBSFORTRAN   = -L$(LAPACKPATH)/lib64 -llapack -L$(BLASPATH)/lib64 -lblas -L$(NETCDFPATH)/lib -lnetcdf -lnetcdff -lstdc++
        LIBSCPP       = -L$(NETCDFPATH)/lib -lnetcdf_c++
endif

ifeq ($(HOSTNAME),green)
        COMPILFORTAN  = ifort
        COMPILCPP     = icc

        NETCDFPATH    = /home/astr/sallaz/opt/netcdf/last
        HDF5PATH      = /home/astr/sallaz/opt/hdf5/last
        LAPACKPATH    = /home/astr/sallaz/opt/lapack/3.3.0/intel/10.1/64
        CURLPATH      = /home/astr/sallaz/opt/curl/last

        OPTIFLAGSFORTRAN   =  -O0 -g -W1 -warn all -static
        OPTIFLAGSCPP  =  -O0 -g -static

        FORTANFLAGS   =  -cpp -convert big_endian -assume byterecl -align dcommon -extend_source -Dlinux86  $(OPTIFLAGSFORTRAN)
        CPPFLAGS =  $(OPTIFLAGSCPP)

        INCS          = -Isrc -I. -I$(NETCDFPATH)/include
        LIBSFORTRAN   = -L$(LAPACKPATH)/lib -llapack -lblas -L$(NETCDFPATH)/lib -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -lhdf5_hl -lhdf5 -lz -lm -lstdc++ -L$(CURLPATH)/lib -lcurl -lrt
        LIBSCPP       = -L$(NETCDFPATH)/lib -lnetcdf_c++4 -lnetcdf_c++ -lnetcdf -L$(HDF5PATH)/lib -lhdf5_hl -lhdf5 -lz -lm -lstdc++  -L$(CURLPATH)/lib -lcurl -lrt
endif

EXE =  moddata_cov

all:$(EXE)

moddata_cov: moddata_cov.f90 readvarcdfCov.f90 month2SaisonAvgCov.f90 covariance.f90
	$(COMPILFORTAN) $(FORTANFLAGS) $(INCS) -o $(<D)/$@ $^ $(LIBSFORTRAN)

%.o: %.cpp
	$(COMPILCPP) -o $@ -c $<

%.o: %.f90
	$(COMPILFORTAN) $(FORTANFLAGS) $(INCS) -module $(<D) -o $@ -c $<

clean2:
	resu=`find . -name "*__genmod.*"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\

clean:
	@for i in ${EXE}; do \
	  resu=`find . -name "$${i}"`; \
	  if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
        done; \
	resu=`find . -name "*.o"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
	resu=`find . -name "*.mod"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
	
	
