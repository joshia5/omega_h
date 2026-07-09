rm -rf *.txt *.tcl *.cmake Makefile install src Testing CMakeFiles
cmake /users/joshia5/lore.scorec.rpi.edu/develop/mfem_omega/omega_h \
-DCMAKE_INSTALL_PREFIX=$PWD/install \
-DBUILD_SHARED_LIBS=ON \
-DOmega_h_USE_CUDA=OFF \
-DOmega_h_USE_Kokkos=OFF \
-DCMAKE_BUILD_TYPE=Debug \
-DCMAKE_CXX_COMPILER=g++ \
-DOmega_h_CHECK_BOUNDS=ON \
-DOmega_h_USE_MPI=ON \
-DOmega_h_USE_SimModSuite=ON \
-DSIM_MPI="mpich4.1.1" \
-DBUILD_TESTING=ON 
make curve_test
