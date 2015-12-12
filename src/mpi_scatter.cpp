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

// node info
char node_name[40];
char mapped_by[40];

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

void init_memory (){

  worker_local_histogram = (int*) malloc ( HIST_SIZE * sizeof (int));
  worker_initial_image = (int*) malloc(elements_per_worker * sizeof ( int ) );
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

  initial_image = (int*) malloc(total_pixels * sizeof ( int ) );
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

void writeResults (int matrix_side , int number_nodes,  char* mapped_by ,  char * node_name ) {
  ofstream file ("timing/timings.dat" , ios::out | ios::app );
  file << number_nodes << " , " << number_processes <<" , "<< mapped_by  << " , " << matrix_side << " x " << matrix_side << " , " << hist_duration << " , " << accum_duration << " , "<< transform_duration << " , " << total_duration <<" , " << node_name << endl;
  file.close();
}
void start_time ( void ){
  initial_time = MPI_Wtime();
}

void stop_time ( void ) {
  final_time = MPI_Wtime();
  hist_duration = hist_exit_time - hist_call_time;
  accum_duration = accum_exit_time - hist_call_time;
  transform_duration = transform_exit_time - accum_call_time;
  total_duration =  final_time - initial_time;
}

///////////////////////////////////////////
/************ master process *************/
///////////////////////////////////////////

void calculate_histogram ( ) {
  hist_call_time = MPI_Wtime();
  MPI_Scatter ( initial_image, elements_per_worker , MPI_INT , worker_initial_image, elements_per_worker, MPI_INT, MASTER, MPI_COMM_WORLD );

  for (long long int pixel_number = 0; pixel_number < elements_per_worker ; ++pixel_number) { 
    worker_local_histogram[ worker_initial_image[pixel_number] ]++;
  }

  if ( process_id == MASTER ){
    master_histograms = (int*) malloc(number_processes * HIST_SIZE * sizeof ( int ) );
  }
  // Gather all partial histograms down to the root process
  MPI_Gather( worker_local_histogram , HIST_SIZE , MPI_INT , master_histograms , HIST_SIZE , MPI_INT , MASTER , MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( process_id == MASTER ){
    for (int pos_hist = 0; pos_hist < HIST_SIZE; ++pos_hist ){
      for ( int p_num = 0; p_num < number_processes; ++p_num ){
        histogram[pos_hist] += master_histograms[p_num*pos_hist];
      }
    }
  }
  hist_exit_time = MPI_Wtime();
}

void calculate_accum ( long long int total_pixels  ){
  accum_call_time = MPI_Wtime();
  if ( process_id == MASTER ){
    int accumulated_value = 0;
    for ( unsigned i = 0 ; i < HIST_SIZE ; i++ ){
      accumulated_value += histogram[i];
      histogram_accumulated[i] = accumulated_value * 255.0f / total_pixels ;
    }
  }
  MPI_Bcast( histogram_accumulated , HIST_SIZE, MPI_FLOAT, MASTER , MPI_COMM_WORLD);
  accum_exit_time = MPI_Wtime();
}

void transform_image( ){
  transform_call_time = MPI_Wtime();
  for (long long int pixel_number = 0; pixel_number < elements_per_worker; ++pixel_number ) {
    worker_final_image[pixel_number] = ( int )( histogram_accumulated [ worker_initial_image[pixel_number] ] );
  }

  // Gather all partial images down to the root process
  MPI_Gather( worker_final_image , elements_per_worker , MPI_INT , final_image , elements_per_worker , MPI_INT , MASTER , MPI_COMM_WORLD);
  transform_exit_time = MPI_Wtime();
}

int main (int argc, char *argv[]) {

  int matrix_side = atoi(argv[1]);
	int number_nodes = atoi(argv[3]);
  int total_pixels = matrix_side * matrix_side;

  if ( argc >=3 ){
    if (argc >= 5 ){
      strcpy (mapped_by,argv[4]);
      strcpy (node_name,argv[5]);
    }


    /**** MPI ****/  
    MPI_CHECK(  MPI_Init(&argc, &argv) );
    MPI_CHECK(  MPI_Comm_rank(MPI_COMM_WORLD, &process_id) );
    MPI_CHECK(  MPI_Comm_size(MPI_COMM_WORLD, &number_processes) );
    elements_per_worker = total_pixels / number_processes;

    //initialize memory
    init_memory();

    if( process_id == MASTER ){
      fillMatrices(total_pixels);
      clearCache();
    }
    start_time();
    /**** FIRST METHOD ****/
    calculate_histogram ( );
    /**** SECOND METHOD ****/
    calculate_accum ( total_pixels );
    /**** THIRD METHOD ****/
    transform_image( );
    stop_time();
    if( process_id == MASTER ){
	writeResults( matrix_side, number_nodes,  mapped_by, node_name );
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

