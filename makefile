#**********************************************************
# Author:
#
#  Jinsong Chen
#  Lawrence Berkeley National Lab
#  Earth Sciences Division
#  Berkeley, CA 94720
#  jchen@lbl.gov
#
#Notices: 

#  (1) The software was developed at Lawrence Berkeley 
#      Laboratory by Jinsong Chen of Earth Sciences Division.
#      The work is supported by the U. S. Department of 
#      Energy. Any commercial use of the software requires
#      PRIOR AGREEMENT with Lawrence Berkeley Laboratory. For
#      further information, please contact Jinsong Chen at
#      jchen@lbl.gov.
#
#  (2) Any use of the software in source and binary forms
#     (with or without modification) is permitted with
#      agreement of citing the following paper:
#
#  "Chen, J., A. Kemna, and S. Hubbard (2008), A comparison
#  between Gauss-Newton and Markov Chain Monte Carlo based
#  methods for inverting spectral induced polarization data,
#  for Cole-Cole parameters, Geophysics, Vol. 73, No. 6."
#***********************************************************
#
#Directory
WORK         =.
SUBFUN       =$(WORK)/SubFun

#Compiler
COMPILE.C    =  g++ -c
COMPILE.PLUS =  g++ -c
LINK.PLUS    =  g++ -v  

FFLAGS       = 
LIBS         = 
##########################################################
#Main-functions
MYMAIN     = $(WORK)/sisip_main.o
ESTMAIN    = $(WORK)/sisip_est.o

#Sub-functions
SAMP    = $(WORK)/sisip_samp.o
POST    = $(WORK)/sisip_post.o
COLE    = $(WORK)/sisip_cole.o
MISC    = $(WORK)/sisip_misc.o
BASE    = $(SUBFUN)/basefuns.o
RAND    = $(SUBFUN)/randfuns.o

########################################################
mainOBJ = $(MYMAIN)     $(SAMP)    $(POST)  \
          $(COLE)       $(MISC)             \
          $(BASE)       $(RAND) 

estOBJ =  $(ESTMAIN)    $(COLE)    $(MISC)    \
          $(BASE)       $(RAND)
#########################################################
INCLUDE   = $(SUBFUN)/basefuns.h \
            $(SUBFUN)/randfuns.h \
            $(SUBFUN)/timer.h    \
            $(WORK)/sisip.h 

##########################################################
# Link objective functions to  obtain executive files

mymain: $(INCLUDE)  $(mainOBJ) 
	$(LINK.PLUS) $(FFLAGS) $(mainOBJ) -o $(WORK)/mymain $(LIBS)
myest: $(INCLUDE)  $(estOBJ) 
	$(LINK.PLUS) $(FFLAGS) $(estOBJ) -o $(WORK)/myest $(LIBS)

##########################################################
# Compile all other sub-functions

$(MYMAIN): $(INCLUDE) $(WORK)/sisip_main.cpp
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_main.cpp -o $@
$(ESTMAIN): $(INCLUDE) $(WORK)/sisip_est.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_est.cpp -o $@
##########################################################
allOBJ= $(mainOBJ)  $(estOBJ)

clean:
	rm -f $(allOBJ) $(allRUN)

#############################################################
# Compile sub-functions

$(COLE): $(INCLUDE)   $(WORK)/sisip_cole.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_cole.cpp -o $@
$(POST): $(INCLUDE)   $(WORK)/sisip_post.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_post.cpp -o $@
$(SAMP): $(INCLUDE)   $(WORK)/sisip_samp.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_samp.cpp -o $@
$(MISC): $(INCLUDE) $(WORK)/sisip_misc.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/sisip_misc.cpp -o $@
$(BASE): $(INCLUDE) $(SUBFUN)/basefuns.c  
	$(COMPILE.C) $(FFLAGS) $(SUBFUN)/basefuns.c -o $@
$(RAND): $(INCLUDE) $(SUBFUN)/randfuns.c 
	$(COMPILE.C) $(FFLAGS) $(SUBFUN)/randfuns.c -o $@
