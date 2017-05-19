
# REQUIRED

# BASEDIR is the location of the package. 
# if not defined, assumed to be PWD
export BASEDIR=$PWD
#export BASEDIR=/home/twongjirad/working/uboone/qosc

# -----------------------------------
# Prob3++
export PROBTPP_DIR=/home/twongjirad/working/t2k/Prob3
export PROBTPP_LIBNAME=ThreeProb_2.10

# -----------------------------------
# T2K REWIGHT
# SUPPORTED VERSIONS: 0, v1r17
# 0 means that T2KREWEIGHT NOT used
export T2KRW_VER=0
# Location of T2KRewight
export T2KREWEIGHT_DIR=
# NEUT_ROOT

# -----------------------------------
# Google Performance Tools
# Used to trace the program
# If installed and want to be used
# Of course, as with all profiling tools, there will be a cost in performance if activated
export USE_GPERFTOOLS=0

# -----------------------------------
# Python Tools
# Cython support available
export BUILD_CYTHON_LIBRARIES=0
export PYTHONINC=/usr/local/python-2.7/include/python2.7
export PYTHONLIB=" -L/usr/local/python-2.7/lib -lpython2.7"

# ---------------------------------------
# Set Paths
export LD_LIBRARY_PATH=$PWD/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${PROBTPP_DIR}:${LD_LIBRARY_PATH}
