#include "RPS.h"
#include <time.h>
#include <omp.h>

void swap_petris();

cell* petri_A;
cell* petri_B;

int main(int argc, char** argv){

  int thread_count = strtol(argv[1], NULL, 10);
  int rank = omp_get_thread_num();

  printf("running %d iterations\n",ITERATIONS);

  srand(time(NULL) + rank * thread_count);

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

    int seed = rand();
    seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;

    # pragma omp parallel for num_threads(thread_count)
    for(int xx = 1; xx < IMG_X - 2; xx++){
      for(int yy = 1; yy < IMG_Y - 2; yy++){
        petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
      }
    }

    swap_petris();
  }
  
  if (rank == 0) {
    make_bmp(petri_A, 0);
  }

}

void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}
