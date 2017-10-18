#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "RPS_MPI.h"

void initialize();
void initialize_petri();
void exchange_borders();
void iterate_CA();
void gather_petri();
void create_types();
void free_all();
void debug();
void add_petri_to_buff(cell** buff, cell* petri, int rec_dim_x, int rec_dim_y, int coords[2]);
void alloc_petri_buff();
void free_petri_buff();
cell pick_nabo();
bool beats(cell me, cell other);
cell update_cell(int x, int y, int dim_x, cell* petri);

int rank;
int size;

// I denote mpi process specific values with hungarian notation, adding a p

// The dimensions of the processor grid. Same for every process
int p_x_dims;
int p_y_dims;

// The location of a process in the process grid. Unique for every process
int p_my_x_dim;
int p_my_y_dim;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int p_local_petri_x_dim;
int p_local_petri_y_dim;
int p_local_petri_size;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t;  // TODO: Implement this
MPI_Datatype border_col_t;  // TODO: Implement this
MPI_Datatype local_petri_t; // Already implemented
MPI_Datatype mpi_cell_t;    // Already implemented

// Each process is responsible for one part of the petri dish.
// Since we can't update the petri-dish in place each process actually
// gets two petri-dishes which they update in a lockstep fashion.
// dish A is updated by writing to dish B, then next step dish B updates dish A.
// (or you can just swap them inbetween iterations)
cell* local_petri_A;
//cell* local_petri_B;

cell** petri_buff;
cell* north;
cell* south;
cell* east;
cell* west;


int main(int argc, char** argv){

  // Ask MPI what size (number of processors) and rank (which process we are)
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  srand(time(NULL) + rank * size);
  rand();

  ////////////////////////////////
  // Create cartesian communicator
  int dims[2];
  dims[1] = p_x_dims;
  dims[0] = p_y_dims;

  int periods[2]; // we set these to 0 because we are not interested in wrap-around
  periods[0] = 0;
  periods[1] = 0;

  int coords[2];
  coords[1] = p_my_x_dim;
  coords[0] = p_my_y_dim;

  //(int nnodes, int ndims, int *dims)
  MPI_Dims_create(size, 2, dims);
  //(MPI_Comm comm_old, int ndims, int *dims, int *periods, int reorder, MPI_Comm *comm_cart)
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  //(MPI_Comm comm, int rank, int maxdims, int *coords)
  MPI_Cart_coords(cart_comm, rank, 2, coords);

  //(MPI_Comm comm, int direction, int displ, int *source, int *dest)
  MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
  MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

  p_x_dims = dims[1];
  p_y_dims = dims[0];

  p_my_x_dim = coords[1];
  p_my_y_dim = coords[0];
  ////////////////////////////////
  ////////////////////////////////

  initialize();


  create_types();

  // A super basic example sending some data:

  // cell* my_test_cell = malloc(10*sizeof(cell));
  // for(int ii = 0; ii < 10; ii++){
  //   my_test_cell[ii].strength = ii;
  //   my_test_cell[ii].color = rank;
  // }

  // if(rank == 0){
  //   cell* rec_buf = malloc(sizeof(cell)*10);
  //   for(int ii = 0; ii < size - 1; ii++){
  //     MPI_Recv(rec_buf, 10, mpi_cell_t, ii+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     printf("receiving from rank %d: \n", ii+1);
  //     for(int jj = 0; jj < 10; jj++){
  //       printf("[%d, %d]  ", rec_buf[jj].color, rec_buf[jj].strength);
  //     }
  //     printf("\n");
  //   }
  // }
  // else{
  //   MPI_Send(my_test_cell, 10, mpi_cell_t, 0, 0, MPI_COMM_WORLD);
  // }

  for (int i = 0; i < ITERATIONS; i++) {
    exchange_borders();
    iterate_CA();
  }

  gather_petri();

  MPI_Finalize();

  if(rank==0){
    // TODO: Write the petri to an image
    make_bmp(petri_buff, 0);
    
    free_petri_buff();
  }

  // You should probably make sure to free your memory here
  // We will dock points for memory leaks, don't let your hard work go to waste!
  // free_stuff()
  free_all();

  exit(0);
}


void create_types(){

  ////////////////////////////////
  ////////////////////////////////
  // cell type
  const int    nitems=2;
  int          blocklengths[2] = {1,1};
  MPI_Datatype types[2] = {MPI_INT, MPI_INT};
  MPI_Aint     offsets[2];

  offsets[0] = offsetof(cell, color);
  offsets[1] = offsetof(cell, strength);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
  MPI_Type_commit(&mpi_cell_t);
  ////////////////////////////////
  ////////////////////////////////



  ////////////////////////////////
  ////////////////////////////////
  // A message for a local petri-dish
  MPI_Type_contiguous(p_local_petri_x_dim * p_local_petri_y_dim,
                      mpi_cell_t,
                      &local_petri_t);
  MPI_Type_commit(&local_petri_t);
  ////////////////////////////////
  ////////////////////////////////


  //TODO: Create MPI types for border exchange

  ////////////////////////////////
  ////////////////////////////////
  // Row data
  //MPI_Type_vector(p_local_petri_x_dim, 1, 1, mpi_cell_t, &border_row_t);
  MPI_Type_contiguous(p_local_petri_x_dim, mpi_cell_t, &border_row_t);
  MPI_Type_commit(&border_row_t);
  ////////////////////////////////
  ////////////////////////////////

  ////////////////////////////////
  ////////////////////////////////
  // Column data
  MPI_Type_vector(p_local_petri_y_dim, 1, p_local_petri_x_dim, mpi_cell_t, &border_col_t);
  MPI_Type_commit(&border_col_t);
  ////////////////////////////////
  ////////////////////////////////

}


void initialize(){
  //TODO: assign the following to something more useful than 0
  p_local_petri_x_dim = p_my_x_dim < p_x_dims - 1 ? (IMG_X / p_x_dims) : (IMG_X / p_x_dims) + (IMG_X % p_x_dims);
  p_local_petri_y_dim = p_my_y_dim < p_y_dims - 1 ? (IMG_Y / p_y_dims) : (IMG_Y / p_y_dims) + (IMG_Y % p_y_dims);
  p_local_petri_size = p_local_petri_x_dim * p_local_petri_y_dim;

  //printf("rank%d; dim: %d:%d; ldim: %d:%d\n", rank, p_my_x_dim, p_my_y_dim, p_local_petri_x_dim, p_local_petri_y_dim);

  // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
  // than just your piece of the petri.
  local_petri_A = malloc(p_local_petri_size * sizeof(cell));
  //local_petri_B = malloc(p_local_petri_size * sizeof(cell));

  north = calloc(p_local_petri_x_dim, sizeof(cell));
  south = calloc(p_local_petri_x_dim, sizeof(cell));
  east = calloc(p_local_petri_y_dim, sizeof(cell));
  west = calloc(p_local_petri_y_dim, sizeof(cell));

  // TODO: Randomly perturb the local dish. Only perturb cells that belong to your process,
  // leave border pixels white.
  for (int i = 0; i < p_local_petri_size; i++) {
    local_petri_A[i].color = rand() % 4;
    local_petri_A[i].strength = 1;

    //////////////
    // If border, set to WHITE
    if (i % p_local_petri_x_dim == 0 || (i + 1) % p_local_petri_x_dim == 0) {
      local_petri_A[i].color = WHITE;
    }
    if (i < p_local_petri_x_dim || i > p_local_petri_x_dim * (p_local_petri_y_dim - 1)) {
      local_petri_A[i].color = WHITE;     
    }
    //////////////
  }

}


void exchange_borders(){
  //TODO: Exchange borders inbetween each step

  int total_reqs = 8;
  int offset_row = p_local_petri_x_dim * (p_local_petri_y_dim - 1);
  int offset_col = p_local_petri_x_dim - 1;

  MPI_Request reqs[total_reqs];
  MPI_Status stats[total_reqs];

  ////////////
  //North-South
  MPI_Irecv(north, 1, border_row_t, p_north, 0, MPI_COMM_WORLD, &reqs[0]);
  MPI_Irecv(south, 1, border_row_t, p_south, 0, MPI_COMM_WORLD, &reqs[1]);

  MPI_Isend(&local_petri_A[0], 1, border_row_t, p_north, 0, MPI_COMM_WORLD, &reqs[2]);
  MPI_Isend(&local_petri_A[offset_row], 1, border_row_t, p_south, 0, MPI_COMM_WORLD, &reqs[3]);
  ////////////

  ////////////
  //East-West
  MPI_Irecv(east, p_local_petri_y_dim, mpi_cell_t, p_east, 0, MPI_COMM_WORLD, &reqs[4]);
  MPI_Irecv(west, p_local_petri_y_dim, mpi_cell_t, p_west, 0, MPI_COMM_WORLD, &reqs[5]);

  MPI_Isend(&local_petri_A[0], 1, border_col_t, p_west, 0, MPI_COMM_WORLD, &reqs[6]);
  MPI_Isend(&local_petri_A[offset_col], 1, border_col_t, p_east, 0, MPI_COMM_WORLD, &reqs[7]);
  ////////////
  
  MPI_Waitall(total_reqs, reqs, stats);

}

void iterate_CA(){
  //TODO: Iterate the cellular automata one step
  int temp_x_dim = p_local_petri_x_dim + 2;
  int temp_y_dim = p_local_petri_y_dim + 2;

  cell* temp_a = calloc(temp_x_dim * temp_y_dim, sizeof(cell));
  cell* temp_b = calloc(temp_x_dim * temp_y_dim, sizeof(cell));

  //////////////
  //Move local petri to temporaly buffer for calculation
  for (int i = 0; i < p_local_petri_size; i++) {
    temp_a[temp_x_dim + 1 + ((i / p_local_petri_x_dim) * 2) + i] = local_petri_A[i];
  }
  //////////////

  //////////////
  //Add north and south border data to buffer
  for (int i = 0; i < p_local_petri_x_dim; i++) {
    temp_a[1 + i] = north[i];
    temp_a[temp_x_dim * (temp_y_dim - 1) + 1 + i] = south[i];
  }
  //////////////

  //////////////
  //Add east and west border data to buffer
  for (int i = 0; i < p_local_petri_y_dim; i++) {
    temp_a[temp_x_dim + (temp_x_dim * i)] = west[i];
    temp_a[temp_x_dim + (temp_x_dim * i) + temp_x_dim - 1] = east[i];
  }
  //////////////

  //////////////
  //Update every cell except the border
  for (int y = 1; y < p_local_petri_y_dim + 1; y++) {
    for (int x = 1; x < p_local_petri_x_dim + 1; x++) {
      temp_b[y * temp_x_dim + x] = update_cell(x, y, temp_x_dim, temp_a);
    }
  }
  //////////////

  //////////////
  //Move buffer to local petri
  for (int i = 0; i < p_local_petri_size; i++) {
    local_petri_A[i] = temp_b[temp_x_dim + 1 + ((i / p_local_petri_x_dim) * 2) + i];
  }
  //////////////

  free(temp_a);
  free(temp_b);

}

cell pick_nabo(int x, int y, cell* petri) {
  int chosen = rand() % 8;
  if (chosen == 4) chosen++;
  int c_x = chosen % 3;
  int c_y = chosen / 3;

  return petri[(y + c_y - 1) * (p_local_petri_x_dim + 2) + (x + c_x - 1)];
}

cell update_cell(int x, int y, int dim_x, cell* petri) {
  cell neighbor_cell = pick_nabo(x, y, petri);
  cell my_cell = petri[y * dim_x + x];
  if(neighbor_cell.color == WHITE){
    return my_cell;
  }
  cell next_cell = my_cell;

  if(my_cell.color == WHITE){
    next_cell.strength = 1;
    next_cell.color = neighbor_cell.color;
    return next_cell;
  }
  else {
    if(beats(my_cell, neighbor_cell)){
      next_cell.strength++;
    }
    else{
      next_cell.strength--;
    }
  }

  if(next_cell.strength == 0){
    next_cell.color = neighbor_cell.color;
    next_cell.strength = 1;
  }

  if(next_cell.strength > 4){
    next_cell.strength = 4;
  }

  return next_cell;
}

void gather_petri(){
  //TODO: Gather the final petri for process rank 0
  if (rank == 0) {
    int tmp_coords[] = {0, 0};

    alloc_petri_buff();

    add_petri_to_buff(petri_buff, local_petri_A, p_local_petri_x_dim, p_local_petri_y_dim, tmp_coords);

    for (int s_rank = 1; s_rank < size; s_rank++) {
      //////////////
      int recv_dim_x;
      int recv_dim_y;
      MPI_Status status;
      MPI_Recv(&recv_dim_x, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&recv_dim_y, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, &status);
      //////////////

      //////////////
      int recv_size = recv_dim_x * recv_dim_y;
      cell* recv_buff = malloc(recv_size * sizeof(cell));
      MPI_Recv(recv_buff, recv_size, local_petri_t, s_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //////////////

      //////////////
      int recv_coords[2];
      MPI_Cart_coords(cart_comm, s_rank, 2, recv_coords);
      //////////////

      add_petri_to_buff(petri_buff, recv_buff, recv_dim_x, recv_dim_y, recv_coords);

      free(recv_buff);
    }
  } else {
    MPI_Send(&p_local_petri_x_dim, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&p_local_petri_y_dim, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(local_petri_A, 1, local_petri_t, 0, 0, MPI_COMM_WORLD);
  }
}

void add_petri_to_buff(cell** buff, cell* petri, int dim_x, int dim_y, int petri_coords[2]) {
  int normal_dim_x = IMG_X / p_x_dims;
  int normal_dim_y = IMG_Y / p_y_dims;
  int start_dim_x = normal_dim_x * petri_coords[1];
  int start_dim_y = normal_dim_y * petri_coords[0];

  for (int x = start_dim_x; x < start_dim_x + dim_x; x++) {
    for (int y = start_dim_y; y < start_dim_y + dim_y; y++) {
      buff[x][y] = petri[x - start_dim_x + (y - start_dim_y) * dim_x];
    }
  }
}

void alloc_petri_buff() {
  petri_buff = malloc(IMG_X * sizeof(cell*));
    for (int i = 0; i < IMG_X; i++) {
      petri_buff[i] = malloc(IMG_Y * sizeof(cell));
  }
}

void free_petri_buff() {
  for (int i = 0; i < IMG_X; i++) {
    free(petri_buff[i]);
  }
  free(petri_buff);
}

void free_all() {
  free(local_petri_A);

  free(north);
  free(south);
  free(east);
  free(west);
}


/*void debug() {
  printf("Dims: %d:%d\n", p_x_dims, p_y_dims);
  printf("Tdims: %d:%d\n", p_local_petri_x_dim, p_local_petri_y_dim);

  printf("North:%d\n", p_north);
  printf("South:%d\n", p_south);
  printf("East:%d\n", p_east);
  printf("West:%d\n", p_west);
}*/
