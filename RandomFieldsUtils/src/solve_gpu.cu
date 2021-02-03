#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cuda_runtime_api.h>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <chrono>
#include <vector>

#include "errors_messages.h"
#include "RandomFieldsUtils.h"
#include "solve_gpu.h"
#include "options.h"




__global__ void logdet_kernel(double *d_matrix, Uint *d_size, double *d_logdet){
    __shared__ double logdet_loc;
    __shared__ double submatrix[THREADS_PER_BLOCK];
    logdet_loc = 0.0;
    *d_logdet = 0.0;
    int idx = blockDim.x * blockIdx.x + threadIdx.x,
        thread = threadIdx.x;
    if(idx < *d_size){
        submatrix[thread] = d_matrix[idx * (*d_size +1)];
    }

    __syncthreads();
    atomicAdd(&logdet_loc, idx >= *d_size ? 0 : (log(submatrix[thread])));

    __syncthreads();
    if(threadIdx.x ==0){atomicAdd(d_logdet, logdet_loc);
    };
};

int cholGPU(bool copy, double *matrix, Uint size, double *B, Uint rhs_cols,
     double *LogDet, double *RESULT){
    /*
        This function solves the problem
            A x = b
        on   an available GPU and writes the solution to the original memory
        Input: 
            matrix: pointer to rowwise allocated matrix A
            individuals: number of individuals in matrix, i.e. dimension
            vector: pointer to vector b
        Ouput:
            vector: contains solution x after the function has been called
    */

    //declare/define process variables
    int bufferSize = 0;
    int *info = NULL;
    int h_info = 0;
    double *buffer = NULL;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;

    //declare device variables
    double *d_matrix = NULL;
    double *d_B = NULL;
    double *d_logdet = NULL;
    Uint *d_size = NULL;

    //initialize handle and stream, calculate buffer size needed for cholesky
    cusolverDnCreate(&handle);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cusolverDnSetStream(handle, stream);

    cusolverDnDpotrf_bufferSize(handle, uplo, size, matrix,
        size, &bufferSize);
    //PRINTF("Buffersize: %f", ((float) bufferSize)/1073741824.0);
    cudaMalloc(&info, sizeof(int));
    cudaMalloc(&buffer, sizeof(double) * bufferSize);
    //allocate memory on device  
    cudaMalloc((void**)&d_matrix, sizeof(double) * size * size);
    cudaMalloc((void **)&d_B, sizeof(double) * size * rhs_cols);
    cudaMemset(info, 0, sizeof(int));

    //copy data to device
    cudaMemcpy(d_matrix, matrix, sizeof(double) * size * size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, sizeof(double) * size * rhs_cols, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    //write cholesky factorization to device copy of A
    cusolverDnDpotrf(handle, uplo, size,
            d_matrix, size, buffer, bufferSize, info);
            
    //Synchronize is necessary, otherwise error code "info" returns nonsense 
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("%s\n", cudaGetErrorString(err));

    //check for errors
    cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    if (0 != h_info) {
        PRINTF("Error: Cholesky factorization failed\n");
    }
    //calculate x = A\b
    cusolverDnDpotrs(handle, uplo, size, rhs_cols, 
            d_matrix, size, d_B,
             size, info);

    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Potrs: %s\n", cudaGetErrorString(err));
    
    if(LogDet != NULL){
        cudaMalloc((void**)&d_logdet, sizeof(double));
        cudaMalloc((void**)&d_size, sizeof(Uint));
        cudaMemcpy(d_size, &size, sizeof(Uint), cudaMemcpyHostToDevice);
        logdet_kernel <<< (size - 1)/THREADS_PER_BLOCK +1 ,THREADS_PER_BLOCK>>> (d_matrix, d_size, d_logdet);
        cudaDeviceSynchronize();
        cudaMemcpy(LogDet, d_logdet, sizeof(double), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaFree(d_size);
        cudaFree(d_logdet);
    }
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Err at Logdet: %s\n", cudaGetErrorString(err));

     //*LogDet = 1.0;
    //copy  solution from device to vector on host
    cudaMemcpy(RESULT, d_B, sizeof(double) * size * rhs_cols, cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Memcpy: %s\n", cudaGetErrorString(err));
    //free allocated memory
    cudaFree(info);
    cudaFree(buffer);
    cudaFree(d_matrix);
    cudaFree(d_B);
    cusolverDnDestroy(handle);
    cudaStreamDestroy(stream);
    return 0;
};



void mgpuSolve(double *matrix, Uint individuals, double *vector){
    /*
        This function solves the problem
            A x = b
        on an MULTIPLE GPUs and writes the solution to the original memory of b
        Input: 
            matrix: pointer to rowwise allocated matrix A
            individuals: number of individuals in matrix, i.e. dimension
            vector: pointer to vector b
        Ouput:
            vector: contains solution x after the function has been called
    */

    // Define auxiliary variables
    cusolverMgHandle_t handle = NULL;
    const int max_devices = 8; // Maximum number of devices to be used
    int nbGpus = 0;
    std::vector<int> deviceList;
    const int N = individuals, lda = N; // Dimension of matrix
    const int IA  = 1;
    const int JA  = 1;
    const int T_A = 256; //Tile size
    const int IB  = 1;
    const int JB  = 1;
    const int T_B = 1000, ldb = N; 
    int info = 0;
    int64_t lwork_potrf = 0, lwork_potrs = 0, lwork = 0 ;

    cudaLibMgMatrixDesc_t descrA, descrB;
    cudaLibMgGrid_t grid;
    double **array_d_A = NULL;
    double **array_d_B = NULL;
    double **array_d_work = NULL;

    // Create handles and select devices
    cusolverStatus_t status = cusolverMgCreate(&handle);
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Handle couldn't be created");
    
    cudaError_t cudaStat = cudaGetDeviceCount( &nbGpus );
    nbGpus = (nbGpus < max_devices)? nbGpus : max_devices;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    int cc_major = prop.major, cc_minor = prop.minor;
    for(int i = 0; i< nbGpus; i++){
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        if(prop.major == cc_major & prop.minor == cc_minor)
                deviceList.push_back(i);
    }
    nbGpus = deviceList.size();
    status = cusolverMgDeviceSelect(
        handle,
        nbGpus,
        &deviceList[0]);
    if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Devices couldn't be selected.");

    // Enable peer access for selected devices
    for(int i = 0; i < nbGpus; i++){
        cudaSetDevice(deviceList[i]);
        for(int j = 0; j< nbGpus; j++){
            if(i==j)continue;
            cudaStat = cudaDeviceEnablePeerAccess(deviceList[j],0);
            if(cudaStat != cudaSuccess)PRINTF("Device %d can't access device %d.",deviceList[i],deviceList[j]);
            PRINTF("Access enabled for devices (%d,%d)",deviceList[i],deviceList[j]);
        }
    }
    // Create device grid for vectors A, B
    status = cusolverMgCreateDeviceGrid(&grid, 1, nbGpus, &deviceList[0], CUDALIBMG_GRID_MAPPING_COL_MAJOR );
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Grid couldn't be created.");

    // Creeate matrix descriptions
    status = cusolverMgCreateMatrixDesc(
        &descrA,
        N,   /* nubmer of rows of (global) A */
        N,   /* number of columns of (global) A */
        N,   /* number or rows in a tile */
        T_A, /* number of columns in a tile */
        CUDA_R_64F,
        grid );
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Matrix descriptions couldn't be created.");
    status = cusolverMgCreateMatrixDesc(
        &descrB,
        N,    /* nubmer of rows of (global) B */
        1, /* number of columns of (global) B */
        N,    /* number or rows in a tile */
        T_B,  /* number of columns in a tile */
        CUDA_R_64F,
        grid );
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Matrix description B couldn't be created.");


    // Allocate arrays of device pointers which point at the memory allocated on each device
    array_d_A = (double**) malloc (sizeof(double*) * nbGpus );
    array_d_B = (double**)malloc(sizeof(double*)*nbGpus);
    array_d_work = (double**)malloc(sizeof(double*)*nbGpus);
    memset(array_d_work, 0, sizeof(void*)*nbGpus);

    // Calculate block size on device
    const int A_num_blks = ( N + T_A - 1) / T_A;
    const int B_num_blks = ( N + T_B - 1) / T_B;
    const int A_blks_per_device = (A_num_blks + nbGpus-1)/nbGpus;
    const int B_blks_per_device = (B_num_blks + nbGpus-1)/nbGpus;

    // Allocate memory on each device
    for( int p = 0 ; p < nbGpus ; p++){
        cudaSetDevice(deviceList[p]);
        cudaStat = cudaMalloc( &(array_d_A[p]), sizeof(double)*lda*T_A*A_blks_per_device );
        if(cudaSuccess != cudaStat)PRINTF("Memory for matrix A couldn't be allocated on device %d.",deviceList[p]);
        cudaStat = cudaMalloc( &(array_d_B[p]), sizeof(double)*ldb*T_B*B_blks_per_device );
        if(cudaSuccess != cudaStat)PRINTF("Memory for matrix B couldn't be allocated on device %d.",deviceList[p]);
    }

    // Copy arrays A and B to device
    for( int k = 0 ; k < A_num_blks ; k++){
    /* k = ibx * nbGpus + p */
        const int p   = (k % nbGpus);
        const int ibx = (k / nbGpus);
        double *h_Ak = matrix + (size_t)lda*T_A*k;
        double *d_Ak = array_d_A[p] + (size_t)lda*T_A*ibx;
        const int width = MIN( T_A, (N - T_A*k) );
        cudaStat = cudaMemcpy(d_Ak, h_Ak, sizeof(double)*lda*width, cudaMemcpyHostToDevice);
        if(cudaSuccess != cudaStat)PRINTF("Matrix A couldn't be copied at block (%d, %d).", p,ibx);
    }
    for( int k = 0 ; k < B_num_blks ; k++){
    /* k = ibx * nbGpus + p */
        const int p   = (k % nbGpus);
        const int ibx = (k / nbGpus);
        double *h_Bk = vector + (size_t) T_B*k;
        double *d_Bk = array_d_B[p] + (size_t) T_B*ibx;
        cudaStat = cudaMemcpy(d_Bk, h_Bk, sizeof(double)*T_B, cudaMemcpyHostToDevice);
        if(cudaSuccess != cudaStat)PRINTF("Matrix B couldn't be copied at block (%d, %d).", p,ibx);
    }

    // Calculate buffersizes necessary for potrf and potrs
    cudaDeviceSynchronize();
    status = cusolverMgPotrf_bufferSize(
        handle,
		CUBLAS_FILL_MODE_LOWER,
        N,
        (void**)array_d_A,
        IA, /* base-1 */
        JA, /* base-1 */
        descrA,
        CUDA_R_64F,
        &lwork_potrf);
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Buffer size potrf couldn't  be calculated");    
    cudaDeviceSynchronize();
    status = cusolverMgPotrs_bufferSize(
        handle,
		CUBLAS_FILL_MODE_LOWER,
        N,
        1, /* number of columns of B */
        (void**)array_d_A,
        IA,
        JA,
        descrA,
        (void**)array_d_B,
        IB,
        JB,
        descrB,
        CUDA_R_64F,
        &lwork_potrs);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Buffersize calculation: %s\n", cudaGetErrorString(err));
    if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Buffer size potrs couldn't  be calculated");    

    lwork = (lwork_potrf > lwork_potrs)? lwork_potrf : lwork_potrs;

    // Allocate workspace size
    for(int idx = 0 ; idx < nbGpus ; idx++){
        int deviceId = deviceList[idx];
        cudaSetDevice( deviceId );
        void *d_workspace = NULL;
        cudaStat = cudaMalloc(&d_workspace, sizeof(double)*lwork);
        if( cudaSuccess != cudaStat )PRINTF("Workspace couldn't be allocated.");
        ((void**)array_d_work )[idx] = d_workspace;
    }

    // Calculate potrf to workspace
    status = cusolverMgPotrf(
        handle,
		CUBLAS_FILL_MODE_LOWER,
        N,   
        (void**)array_d_A,
        IA,
        JA,
        descrA,
        CUDA_R_64F,
        (void**)array_d_work,
        lwork,
        &info  /* host */
    );
    cudaDeviceSynchronize();
    if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Potrf couldn't be calculated");
    if(info != 0)PRINTF("Info code %d", info);
    // Calculate potrs to B
    status = cusolverMgPotrs(
        handle,
		CUBLAS_FILL_MODE_LOWER,
        N,
        1, /* number of columns of B */
        (void**)array_d_A,
        IA,
        JA,
        descrA,
        (void**)array_d_B,
        IB,
        JB,
        descrB,
        CUDA_R_64F,
        (void**)array_d_work,
        lwork,
        &info  /* host */
    );
    cudaDeviceSynchronize();
    if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Potrs couldn't be calculated");
    if(info != 0)PRINTF("Info code %d", info);

    // Copy solution B back to host
    for( int k = 0 ; k < B_num_blks ; k++){
    /* k = ibx * nbGpus + p */
        const int p   = (k % nbGpus);
        const int ibx = (k / nbGpus);
        double *h_Bk = vector + (size_t) T_B*k;
        double *d_Bk = array_d_B[p] + (size_t) T_B*ibx;
        cudaStat = cudaMemcpy(h_Bk, d_Bk, sizeof(double)*T_B, cudaMemcpyDeviceToHost);
        if(cudaSuccess != cudaStat)PRINTF("Matrix B couldn't be copied at block (%d, %d).", p,ibx);
    }

    // Free memory on device and host
    for(int i = 0; i< nbGpus; i++){
        cudaSetDevice(deviceList[i]);
        cudaDeviceReset();
    }
    free(array_d_A); free(array_d_B); free(array_d_work);
}

void gpu_relmat_cublas(Uint* M, double* A, Uint snps, Uint individuals){
    /*
        Calculates the crossproduct of M with cublas and stores the result in A.
        Input:
            M: non-encoded matrix of dimension snps x indiv (k x n) storing genomic information
            A: pointer to result matrix
            snps: Number of snps
            individuals: number of individuals
        Output:
            A: matrix containing the type-casted result of M^T * M
        
        Note: cublas is fortran based and therefore assumes M is column-major. Therefore to calculate
            A we instruct cublasgemmex to calculate M * M^T and adjust its parameters.
            Furthermore, gemmex requires the matrix M to have a row number that is a multiple of four
            Therefore this function implements a zero-padding to add extra rows
    */
    
    //Define auxiliary variables as needed for gemmex
        Uint n = individuals;
        Uint m = individuals;
        Uint k = snps;
    
    //Auxiliary padding variables for padding
        Uint k_pad_diff = (PADDIM - k % PADDIM) % PADDIM;
        Uint k_pad = k + k_pad_diff;
        Uint dim = m * k_pad;
    
    //Start timing copy and calculation time
    #ifdef DEBUG
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_start;
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_stop;
        timer_start = std::chrono::high_resolution_clock::now();
    #endif
    //Declare cublas variables and allocate memory
        cublasHandle_t handle;
        cublasCreate(&handle);
        int8_t *d_M, *h_M;
        int32_t *d_C, *h_C;
        int32_t alpha = 1.f;
        int32_t beta = 0.f;
        cudaMalloc(&d_M, sizeof(int8_t) * dim);
        cudaMalloc(&d_C, sizeof(int32_t) * n * m );
        cudaMallocHost((void **)&h_M, sizeof(int8_t) * dim);
        cudaMallocHost((void **)&h_C, sizeof(int32_t) * n * m);
    
    
    
    //Type-cast matrix M to int8 and store the result in page-locked memory
    //Zero-pad matrix to get a row number that is a multiple of four
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(CORES)   
    #endif
        for(int i = 0; i < n; i++){
            for(int j = 0; j < k_pad; j++){
            h_M[j + i * k_pad] = (int8_t) (j< k ?  M[j + i * k] : 0 );
            }
        }
    
    
    //Copy int8 matrix to device
    cudaMemcpy(d_M, h_M, sizeof(int8_t) * dim, cudaMemcpyHostToDevice);  

    //Calculate the crossproduct and check for errros
        cublasStatus_t stat = cublasGemmEx(handle,
            CUBLAS_OP_T,
            CUBLAS_OP_N,
            n,
            m,
            k_pad,
            &alpha,
            d_M,
            CUDA_R_8I,
            k_pad, // I have no idea why this doesnt need to be individuals, same below
            d_M,
            CUDA_R_8I,
            k_pad,
            &beta,
            d_C,
            CUDA_R_32I,
            n,
            CUDA_R_32I, //CUBLAS_COMPUTE_32I,
            CUBLAS_GEMM_DEFAULT
            );
        
        if(stat) PRINTF("GemmEx failed.");
        cudaDeviceSynchronize();
    
    
    //copy result back to host
        cudaMemcpy(h_C, d_C, sizeof(int32_t) * n * m, cudaMemcpyDeviceToHost);
    
    //Convert result to double and store it in output matrix A
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(CORES)   
    #endif
        for (int i = 0; i < n * m; i++) A[i] = (double) h_C[i];
    
    //Free memory 
        cublasDestroy(handle);
        cudaFree(d_M);
        cudaFree(d_C);
        cudaFreeHost(h_C);
        cudaFreeHost(h_M);
    
    //Stop timer
    #ifdef DEBUG
        timer_stop = std::chrono::high_resolution_clock::now();
        PRINTF("Time: %.3f s\n", ((float) std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count())/1000000.0 );
    #endif
    } 
