//////////////////////
// copyright Filipe Oliveira and Sérgio Caldas
// Universidade do Minho
// Parallel Computing Paradigms
// 2015, December
// ///////////////////

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include "mpi.h"
#include <omp.h>

#define HIST_SIZE 256
#define TIME_RESOLUTION 1000000 // time measuring resolution (us)

//OMP THREADS
#define MAX_THREADS 48

//MPI definitions
#define MASTER 0       /* id of the first process */ 

#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    cerr << "MPI error calling \""#call"\"\n"; \
    exit(-1); };

#define TIMER_HIST 1;
#define TIMER_ACCU 2;
#define TIMER_TRAN 3;

using namespace std;

// time control
double initial_time, final_time, temporary_time_start, temporary_time_stop, total_duration;
double hist_call_time, hist_exit_time, hist_duration, hist_compute_duration, hist_transmit_duration, hist_iddle_duration; 
double accum_call_time, accum_exit_time, accum_duration, accum_compute_duration, accum_transmit_duration, accum_iddle_duration;
double transform_call_time, transform_exit_time, transform_duration, transform_compute_duration, tranform_transmit_duration, transform_iddle_duration;

// case study info
char mapped_by[40];
char comm_type[40];
char property [40];

// image matrixes
int * worker_initial_image; 
int * worker_final_image; 
int * initial_image, * final_image;
int * worker_local_histogram;

// image histograms
int histogram[HIST_SIZE];
float histogram_accumulated[HIST_SIZE];
int* master_histograms;
int* master_partial_images;
// mpi
MPI_Status status; 
int process_id;
int number_processes; 
int elements_per_worker;

// open mp
int thread_count;

/////////////////////////////////////////
/////////////////////////////////////////

void init_memory ( long long int total_pixels  ){

  worker_local_histogram = (int*) malloc ( HIST_SIZE * sizeof (int));
  initial_image = (int*) malloc( total_pixels * sizeof ( int ) );
  worker_final_image = ( int* ) malloc( elements_per_worker * sizeof ( int ) );

  memset ( histogram , 0 , sizeof(int) * HIST_SIZE);
  memset ( histogram_accumulated , 0 , sizeof(float) * HIST_SIZE);
  memset ( worker_local_histogram , 0 , sizeof(int) * HIST_SIZE );

}

void free_memory (){
  free( worker_local_histogram ); 
  free( worker_initial_image );
  free( worker_final_image );
  if ( process_id == MASTER ){
    free( initial_image );
    free( final_image );
  }
}

void fillMatrices ( long long int total_pixels  ) {

  final_image = (int*) malloc(total_pixels * sizeof ( int ) );

  for (long long int pixel_number = 0; pixel_number < total_pixels; ++pixel_number) {
    initial_image[pixel_number] = (((int) rand()) % ((int) 255));
  }
}

void clearCache(){
  double clearcache [30000000];
  for (unsigned i = 0; i < 30000000; ++i)
    clearcache[i] = i;
}

void writeResults (int matrix_side , char* comm_type, char* mapped_by ,  char* property , int number_nodes , int number_omp_threads ) {
  std::ostringstream filenameSS;
  filenameSS << "timing/" << number_nodes << "nodes_EXTRA_OMP_"<< matrix_side << "_" << comm_type << "_" <<  property <<".csv";
        std::string filename = filenameSS.str(); // get string out of stream  
ofstream file (filename.c_str() , ios::out | ios::app );
  file << number_processes << " , " << number_omp_threads <<" , "<< mapped_by  << " , " << matrix_side << " x " << matrix_side << " , " << hist_duration << " , " << accum_duration << " , "<< transform_duration << " , " << total_duration << endl;
  file.close();
}

void start_time ( void ){
  initial_time = MPI_Wtime();
}

void stop_time ( void ) {
  final_time = MPI_Wtime();
  hist_duration = hist_exit_time - hist_call_time;
  accum_duration = accum_exit_time - accum_call_time;
  transform_duration = transform_exit_time - transform_call_time;
  total_duration =  final_time - initial_time;
}

///////////////////////////////////////////
/************ master process *************/
///////////////////////////////////////////

void calculate_histogram ( long long int total_pixels , int thread_count ) {
  hist_call_time = MPI_Wtime();

  MPI_Bcast( initial_image , total_pixels , MPI_INT , MASTER , MPI_COMM_WORLD);

#pragma omp parallel num_threads( thread_count ) 
  {
    int thread_id = omp_get_thread_num();
    int thread_histogram[MAX_THREADS][HIST_SIZE];
#pragma omp for nowait schedule (static)
  for (long long int pixel_number = 0; pixel_number < total_pixels ; ++pixel_number) {
        thread_histogram[thread_id][initial_image[pixel_number]]++;
}
    for ( unsigned pos_hist_local = 0; pos_hist_local < HIST_SIZE; ++pos_hist_local ){
#pragma omp atomic
 worker_local_histogram[pos_hist_local] += thread_histogram[thread_id][pos_hist_local];
  }
}
  hist_exit_time = MPI_Wtime();
}

void calculate_accum ( long long int total_pixels  ){
  accum_call_time = MPI_Wtime();
    int accumulated_value = 0;
    for ( unsigned i = 0 ; i < HIST_SIZE ; i++ ){
      accumulated_value += histogram[i];
      histogram_accumulated[i] = accumulated_value * 255.0f / total_pixels ;
    }
  accum_exit_time = MPI_Wtime();
}

void transform_image( int thread_count ){
  transform_call_time = MPI_Wtime();
int range_top, range_min;
range_min = elements_per_worker*process_id;
#pragma omp parallel num_threads( thread_count ) 
  {
#pragma omp for nowait schedule (static)

for (long long int pixel_number = 0; pixel_number < elements_per_worker; ++pixel_number ) {
    worker_final_image[pixel_number] = ( int )( histogram_accumulated [ initial_image[pixel_number+range_min] ] );
  }
}  
  // Gather all partial images down to the root process
  MPI_Gather( worker_final_image , elements_per_worker , MPI_INT , final_image , elements_per_worker , MPI_INT , MASTER , MPI_COMM_WORLD);
  transform_exit_time = MPI_Wtime();
}

int main (int argc, char *argv[]) {
  int matrix_side = atoi(argv[1]);
  int total_pixels = matrix_side * matrix_side;
int number_nodes = 0;
int number_threads = 1;
if ( argc >=2 ){
    if (argc >= 3 ){
      strcpy (mapped_by,argv[2]);
    }
if (argc >= 4 ){
      strcpy (comm_type,argv[3]);
    }
if (argc >= 5 ){
      strcpy (property,argv[4]);
    }
if (argc >= 6 ){
        number_nodes = atoi(argv[5]);
}

 if (argc >= 7 ){
    number_threads = atoi(argv[6]);
}

    /**** MPI ****/  
    MPI_CHECK(  MPI_Init(&argc, &argv) );
    MPI_CHECK(  MPI_Comm_rank(MPI_COMM_WORLD, &process_id) );
    MPI_CHECK(  MPI_Comm_size(MPI_COMM_WORLD, &number_processes) );
    elements_per_worker = total_pixels / number_processes;
    //initialize memory
    init_memory( total_pixels );

    if( process_id == MASTER ){
      fillMatrices(total_pixels);
      clearCache();
    }
    start_time();
    /**** FIRST METHOD ****/
    calculate_histogram ( total_pixels , number_threads );
    /**** SECOND METHOD ****/
    calculate_accum ( total_pixels );
    /**** THIRD METHOD ****/
    transform_image( number_threads );
    stop_time();
    if( process_id == MASTER ){
      writeResults( matrix_side , comm_type, mapped_by ,  property , number_nodes , number_threads );

}
    //free memory
    free_memory();
    MPI_CHECK( MPI_Finalize() );
    return 0;
  }
  else {
    return 1;
  }
}

