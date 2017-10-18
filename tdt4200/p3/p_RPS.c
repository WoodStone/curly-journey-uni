#include "RPS.h"
#include <time.h>
#include <pthread.h>

void swap_petris();
void iterate(void* rank);

cell* petri_A;
cell* petri_B;
int thread_count;
int seed;

int main(int argc, char** argv){
  int thread;
  pthread_t* thread_handles;

  thread_count = strtol(argv[1], NULL, 10);

  thread_handles = malloc(thread_count * sizeof(pthread_t));

  printf("running %d iterations\n",ITERATIONS);

  srand(time(NULL) + thread_count*2);
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

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
    seed = rand();
    seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;

    for (thread = 0; thread < thread_count; thread++) {
      pthread_create(&thread_handles[thread], NULL, iterate, (void*) thread);
    }

    for (thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
    }

    swap_petris();
  }

  make_bmp(petri_A, "RPS_pthread");

}

void iterate(void* rank) {
  int local_rank = (int) rank;

  for(int xx = local_rank + 1; xx < IMG_X - 2; xx = xx + thread_count){
    for(int yy = 1; yy < IMG_Y - 2; yy++){
      petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
    }
  }

  return NULL;
}

void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}
