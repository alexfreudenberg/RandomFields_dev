# Cracked version of RandomFields 3.6
Enhancements: 
* Add slot for extra optimizer
* Printing tryCatch erros
* Printing name of optimiser

How to get fast R:
* Install Intel MKL: `sudo apt-get install intel-mkl-full`
* Follow https://cran.r-project.org/doc/manuals/r-release/R-admin.html#MKL
* Compile R from source with `./configure --with-blas='-mkl=parallel' CXXFLAGS='-march=native -g -O2 -fopenmp' CFLAGS='-march=native -g -O2 -fopenmp' LDFLAGS='-fopenmp'` and `make`
* Enter `export MKL_INTERFACE_LAYER=GNU,LP64` and `export MKL_THREADING_LAYER=GNU` again
