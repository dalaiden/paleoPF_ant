#
include make.macros

VPATH = .:./filter:./covarianceMatrix:./fcost

EXE =  moddata_co PartFilter

all:$(EXE)

moddata_co: moddata_co.f90 filter_scale.f90 readvarcdf.f90 month2SaisonAvg.f90 cost_euclide.f90 covariance.f90 makename.o
	$(COMPILFORTAN) $(FORTANFLAGS) $(INCS) -o $(<D)/$@ $^ $(LIBSFORTRAN)

PartFilter: PartFilter.f90 ./ran_mod.o ./makename.o
	$(COMPILFORTAN) $(FORTANFLAGS) $(INCS) -o $(<D)/$@ $^ $(LIBSFORTRAN)

%.o: %.cpp
	$(COMPILCPP) $(CPPFLAGS) -o $@ -c $<

%.o: %.f90
	$(COMPILFORTAN) $(FORTANFLAGS) $(INCS) -module $(<D) -o $@ -c $<

clean2:
	resu=`find . -name "*__genmod.*"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\

clean:
	@for i in ${EXE}; do \
	  resu=`find . -name "$${i}"`; \
	  if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
        done; \
	resu=`find . -name "*.mod"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
	resu=`find . -name "*.o"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;\
	resu=`find . -name "*__genmod.*"`; if [ "$${resu}" != "" ]; then echo "rm -f $${resu}"; rm -f $${resu}; fi;
