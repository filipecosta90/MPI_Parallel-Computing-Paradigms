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

//MPI
#define MASTER 0       /* id of the first process */ 
#define FROM_MASTER 1  /* setting a message type */ 
#define FROM_WORKER 2  /* setting a message type */ 

using namespace std;

// globals
int  histogram[HIST_SIZE];
float acumulado[HIST_SIZE];
timeval t;
long long int * initial_image, * final_image;
long long unsigned initial_time, final_time, hist_time, accum_time, transform_time, temporary_time;
long long unsigned hist_duration, accum_duration, transform_duration, total_duration;

MPI_Status status; 
int process_id;
int number_workers;
int number_processes; 
int elements_per_send_rcv;
int extra_send_rcv;
int offset;
int message_type;
int num_elements_to_send_rcv;

void fillMatrices ( long long int total_pixels  ) {

  initial_image = (long long int*) malloc(total_pixels * sizeof ( long long int ) );
  final_image = (long long int*) malloc(total_pixels * sizeof ( long long int ) );

  for (long long int pixel_number = 0; pixel_number < total_pixels; ++pixel_number) {
    initial_image[pixel_number] = (((int) rand()) % ((int) 255));
  }

  for ( unsigned  pos = 0; pos < HIST_SIZE ; pos++ ){
    histogram[pos]=0;
    acumulado[pos] = 0.0f;
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

void calcula_histograma ( long long int total_pixels , int thread_count ){

  if (process_id == MASTER){

    /************ master process *************/

    elements_per_send_rcv = total_pixels / number_workers;
    //the extra will go to first worker
    extra_send_rcv = total_pixels % number_workers;
    offset = 0;
    message_type = FROM_MASTER;
    num_elements_to_send_rcv = elements_per_send_rcv + extra_send_rcv;

    for (int dest_worker = 1 ; dest_worker <= number_workers; dest_worker++) {
      printf("sending to process %d histogram from %d to %d\n", dest_worker, offset, offset+num_elements_to_send_rcv);
      //send the number of elements to read
      MPI_Send(&num_elements_to_send_rcv , 1 ,MPI_INT, dest_worker , message_type , MPI_COMM_WORLD);
      //send the real histogram data
      MPI_Send( &initial_image[offset] , num_elements_to_send_rcv , MPI_INT , dest_worker , message_type , MPI_COMM_WORLD);
      
      offset = offset + num_elements_to_send_rcv;
      
      //after first worker set the number of elements to send to the normal
      num_elements_to_send_rcv = elements_per_send_rcv;
    }

    message_type = FROM_WORKER;

    int local_histogram[number_workers][HIST_SIZE];
    for (int source_worker = 1; source_worker<=number_workers; source_worker++) {
      //    MPI_Recv(&num_elements_to_send_rcv , 1 , MPI_INT , source_worker , message_type , MPI_COMM_WORLD , &status);
      MPI_Recv(&local_histogram[source_worker-1][0] , HIST_SIZE , MPI_INT , source_worker , message_type , MPI_COMM_WORLD, &status);
    }
  }

  else{
    /************ worker processes *************/

    message_type = FROM_MASTER;
    int source = MASTER;
    long long int * worker_image; 
    int mpi_histogram[HIST_SIZE];

    // reeive number of elements
    MPI_Recv(&num_elements_to_send_rcv,1,MPI_INT, source , message_type , MPI_COMM_WORLD,&status);
    printf("\t\tworker %d going to handle  %d elements \n", process_id,  num_elements_to_send_rcv );

    worker_image = (long long int*) malloc(num_elements_to_send_rcv * sizeof ( long long int ) );
    MPI_Recv(&worker_image[0], num_elements_to_send_rcv , MPI_INT , source , message_type , MPI_COMM_WORLD , &status);
#pragma omp parallel num_threads( thread_count ) 
    {

      int thread_id = omp_get_thread_num();
      int local_histogram[MAX_THREADS][HIST_SIZE];
#pragma omp for nowait schedule (static)
      for (long long int pixel_number = 0; pixel_number < num_elements_to_send_rcv ; ++pixel_number) { 
        local_histogram[thread_id][ worker_image[pixel_number] ]++;
      }
      
      for ( unsigned pos_hist_local = 0; pos_hist_local < HIST_SIZE; ++pos_hist_local ){
#pragma omp atomic 
        mpi_histogram[pos_hist_local] += local_histogram[thread_id][pos_hist_local];
      }
    }
    printf("\t\t\t###%dcaculated local histogram\t ### going to send o papa!!\n", process_id);
    message_type = FROM_WORKER;
    //send the real histogram data
    MPI_Send( &mpi_histogram , HIST_SIZE  , MPI_INT , MASTER , message_type , MPI_COMM_WORLD);
  }
}

void calcula_acumulado ( long long int total_pixels  ){
  int valor_acumulado = 0;
  for ( unsigned i = 0 ; i < HIST_SIZE ; i++ ){
    valor_acumulado += histogram[i];
    acumulado[i] = valor_acumulado * 255.0f / total_pixels ;
  }
}

void transforma_imagem( long long int total_pixels , int thread_count  ){
  if(process_id == MASTER){
#pragma omp parallel num_threads( thread_count ) 
    {
#pragma omp for nowait schedule (static)
      for (long long int pixel_number = 0; pixel_number < total_pixels; ++pixel_number ) {
        final_image[pixel_number] = (int )( acumulado[ initial_image[pixel_number]] );
      }
    } 
  }
  else {
    printf("worker %d in TRANSFORMA_IMAGEM \n", process_id );
  }
}

int main (int argc, char *argv[]) {
  /**** MPI ****/  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
  MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
  number_workers = number_processes-1; 

  if ( argc > 3 ){

    printf("\tprocess: %d initialized \n", process_id);
    int number_threads = atoi(argv[1]);
    int rows = atoi(argv[2]);
    int columns = atoi(argv[3]);

    long long int total_pixels = rows * columns;
    char node_name[40];
    if (argc > 4 ){
      strcpy (node_name,argv[4]);
    }
    if (number_threads > MAX_THREADS ){
      number_threads = MAX_THREADS;
    }
    if( process_id == MASTER ){
      fillMatrices(total_pixels);
      clearCache();
      printf("\tmatrix: %d*%d \n", rows, columns);
      start();

    }
    calcula_histograma( total_pixels , number_threads );
    if ( process_id == MASTER ){  
      mark_time(1);
      calcula_acumulado( total_pixels );
      mark_time(2);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    transforma_imagem( total_pixels , number_threads );
    MPI_Barrier(MPI_COMM_WORLD);
    if ( process_id == MASTER){
      mark_time(3);
      stop();
      writeResults(number_threads , rows, columns, node_name );
    }
    printf("process: %d finalized \n", process_id);
    MPI_Finalize();
    return 0;
  }
  else {
    return 1;
  }
}

