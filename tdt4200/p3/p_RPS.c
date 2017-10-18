#include "RPS.h"
#include <time.h>
#include <pthread.h>

void swap_petris();

cell* petri_A;
cell* petri_B;
int thread_count;

int main(int argc, char** argv){
  long thread;
  pthread_t* thread_handles;

  thread_count = strtol(argv[1], NULL, 10);

  thread_handles = malloc(thread_count * sizeof(pthread_t));

  printf("running %d iterations\n",ITERATIONS);

  srand(time(NULL));
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

  int seed = rand();

  // Seed some CAs
  for(int ii = 0; ii < 100; ii++){
    int rx = rand() % (IMG_X - 1);
    int ry = rand() % (IMG_Y - 1);
    int rt = rand() % 4;

    petri_A[TRANS(rx,ry)].color = rt;
    petri_A[TRANS(rx,ry)].strength = 1;
  }

  for(int ii = 0; ii < ITERATIONS; ii++){

    // This should be parallelized somehow
    iterate_image(petri_A, petri_B);

    swap_petris();
  }
  make_bmp(petri_A, 0);

}


void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}
