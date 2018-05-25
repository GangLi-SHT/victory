FC          = mpiifort
FC_FLAGS    = -traceback -check bounds -warn all 
DIR         = SET YOUR VICTORY ROOT DIRECTORY HERE
BIN         = $(DIR)/bin
FLIB        = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -liomp5 -lpthread 
MODPATH     = -I$(DIR)/src/SRC_MOD

LLIB = $(DIR)/lib/lib.a $(DIR)/lib/libdfftpack.a

LIB = $(LLIB) $(FLIB)

#
#  Archiver and its flags
#

ARCH   = ar
ARFLAG = cr
RANLIB = ranlib

# set flags 
FLAGS = $(FC_FLAGS)

