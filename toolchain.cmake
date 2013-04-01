# Specify the target system (this allows cross-compiling)
SET(CMAKE_SYSTEM_NAME CRAYXT_COMPUTE_LINUX)  

# specify the cross compiler
set(CMAKE_C_COMPILER /opt/cray/xt-asyncpe/5.17/bin/cc)
set(CMAKE_CXX_COMPILER /opt/cray/xt-asyncpe/5.17/bin/CC)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search 
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
