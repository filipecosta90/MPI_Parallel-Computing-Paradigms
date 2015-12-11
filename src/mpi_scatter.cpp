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

#include <sys/time.h>
#include "mpi.h"
#include <omp.h>

#define HIST_SIZE 256
#define TIME_RESOLUTION 1000000 // time measuring resolution (us)

//OMP THREADS
#define MAX_THREADS 48

//MPI definitions
#define MASTER 0       /* id of the first process */ 
#define FROM_MASTER 1  /* setting a message type */ 
#define FROM_WORKER 2  /* setting a message type */ 

#define MPI_CHECK(call) \
  if((call) != MPI_SUCCESS) { \
    cerr << "MPI error calling \""#call"\"\n"; \
    exit(-1); };

using namespace std;

// time control
long long unsigned initial_time, final_time, hist_time, accum_time, transform_time, temporary_time;
long long unsigned hist_duration, accum_duration, transform_duration, total_duration;
timeval t;

// node info
char node_name[40];

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
int number_workers;
int number_processes; 
int elements_per_worker;
int offset;
int message_type;

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

void writeResults (int number_threads , int rows, int columns ,  char * node_name ) {
  ofstream file ("timing/timings.dat" , ios::out | ios::app );
  file << number_threads << " , " << hist_duration << " , " << accum_duration << " , "<< transform_duration << " , " << total_duration << " , " << rows <<" x "<<  columns  <<" , " << node_name << endl;
  file.close();
}

void start (void) {
  gettimeofday(&t, NULL);
  initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

void mark_time ( int break_num ) {
  gettimeofday(&t, NULL);
  temporary_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
  if ( break_num == 1 ){
    hist_time = temporary_time;
  } else if ( break_num == 2 ){
    accum_time = temporary_time;
  }
  else if ( break_num == 3 ){
    transform_time = temporary_time;
  }
}

void stop ( void ) {
  gettimeofday(&t, NULL);
  final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
  hist_duration = hist_time - initial_time;
  accum_duration = accum_time - hist_time;
  transform_duration = transform_time - accum_time;
  total_duration =  final_time - initial_time;
}

///////////////////////////////////////////
/************ master process *************/
///////////////////////////////////////////

void calculate_histogram ( ) {

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
}

void calculate_accum ( long long int total_pixels  ){
  if ( process_id == MASTER ){
    int accumulated_value = 0;
    for ( unsigned i = 0 ; i < HIST_SIZE ; i++ ){
      accumulated_value += histogram[i];
      histogram_accumulated[i] = accumulated_value * 255.0f / total_pixels ;
    }
  }
  message_type = FROM_MASTER;
  MPI_Bcast( histogram_accumulated , HIST_SIZE, MPI_FLOAT, message_type, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}

void transform_image( ){

  for (long long int pixel_number = 0; pixel_number < elements_per_worker; ++pixel_number ) {
    worker_final_image[pixel_number] = ( int )( histogram_accumulated [ worker_initial_image[pixel_number] ] );
  }

  // Gather all partial images down to the root process
  MPI_Gather( worker_final_image , elements_per_worker , MPI_INT , final_image , elements_per_worker , MPI_INT , MASTER , MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}

int main (int argc, char *argv[]) {

  int rows = atoi(argv[2]);
  int columns = atoi(argv[3]);
  int total_pixels = rows * columns;
  int number_threads = atoi(argv[1]);

  if ( argc > 3 ){
    if (argc > 4 ){
      strcpy (node_name,argv[4]);
    }

    if (number_threads > MAX_THREADS ){
      number_threads = MAX_THREADS;
    }

    /**** MPI ****/  
    MPI_CHECK(  MPI_Init(&argc, &argv) );
    MPI_CHECK(  MPI_Comm_rank(MPI_COMM_WORLD, &process_id) );
    MPI_CHECK(  MPI_Comm_size(MPI_COMM_WORLD, &number_processes) );
    elements_per_worker = total_pixels / number_processes;

    printf("\tProcess: %d initialized!\n", process_id);

    //initiaze accum and hist 
    init_memory ();

    if( process_id == MASTER ){
      fillMatrices(total_pixels);
      clearCache();
      start();
    }
    /**** FIRST METHOD ****/
    calculate_histogram ( );
    if ( process_id == MASTER ){ 
      mark_time(1);
    }
    /**** SECOND METHOD ****/
    calculate_accum ( total_pixels );
    if ( process_id == MASTER ){
      mark_time(2);
    }
    /**** THIRD METHOD ****/
    transform_image( );
    if ( process_id == MASTER ){
      mark_time(3);
      stop();
      writeResults(number_threads , rows, columns, node_name );
    }
    MPI_CHECK( MPI_Finalize() );
    printf("\tProcess: %d finalized!\n", process_id);
    return 0;
  }
  else {
    return 1;
  }
}

