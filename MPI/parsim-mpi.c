#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <mpi.h>
#define G 6.67408e-11
#define EPSILON2 (0.005*0.005)
#define DELTAT 0.1

unsigned int seed;
long long *work_size; // number of cells that the process computes
int rank; // id of the process
int psize; // number of processes
int start_point; // global id of first cell in the process


MPI_Datatype MPI_CENTER_MASS;
MPI_Datatype MPI_PARTICLE_T;

typedef struct particle_{
    long long id; // particle id
    double x; //x coordinate
    double y; //y coordinate
    double vx; //velocity along x axis
    double vy; //velocity along y axis
    double Fx; //resultant of the forces along x axis
    double Fy; //resultant of the forces along y axis
    double ax; //particle acceleration along x axis
    double ay; //particle acceleration along y axis
    double m; //mass of the particle
    int alive; //1 if alive, 0 if evaporated/collided
}particle_t;

typedef struct parcell_t{
    int n_particles;
    particle_t *par;
    int size;
}parcell;

typedef struct cm_{
    double X; //X center of mass
    double Y; //Y of center of mass
    double M; //sum of the mass of each particle
}center_mass;

typedef struct {
    int a;
    int b;
    int cell;
} ParColision;

void create_mpi_center_mass_type() {
    int block_lengths[3] = {1, 1, 1}; 
    MPI_Aint offsets[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; 

    offsets[0] = offsetof(center_mass, X);
    offsets[1] = offsetof(center_mass, Y);
    offsets[2] = offsetof(center_mass, M);

    MPI_Type_create_struct(3, block_lengths, offsets, types, &MPI_CENTER_MASS);
    MPI_Type_commit(&MPI_CENTER_MASS);
}

void create_mpi_particle_t(){
    MPI_Aint displacements[11];
    int blocklengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[11] = {
        MPI_LONG_LONG,  // id
        MPI_DOUBLE,     // x
        MPI_DOUBLE,     // y
        MPI_DOUBLE,     // vx
        MPI_DOUBLE,     // vy
        MPI_DOUBLE,     // Fx
        MPI_DOUBLE,     // Fy
        MPI_DOUBLE,     // ax
        MPI_DOUBLE,     // ay
        MPI_DOUBLE,     // m
        MPI_INT         // alive
    };

    displacements[0] = offsetof(particle_t, id);
    displacements[1] = offsetof(particle_t, x);
    displacements[2] = offsetof(particle_t, y);
    displacements[3] = offsetof(particle_t, vx);
    displacements[4] = offsetof(particle_t, vy);
    displacements[5] = offsetof(particle_t, Fx);
    displacements[6] = offsetof(particle_t, Fy);
    displacements[7] = offsetof(particle_t, ax);
    displacements[8] = offsetof(particle_t, ay);
    displacements[9] = offsetof(particle_t, m);
    displacements[10] = offsetof(particle_t, alive);

    MPI_Type_create_struct(11, blocklengths, displacements, types, &MPI_PARTICLE_T);
    MPI_Type_commit(&MPI_PARTICLE_T);
}

void init_r4uni(int input_seed){
    seed = input_seed + 987654321;
}

double rnd_uniform01(){
    int seed_in = seed;
    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);
    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

double rnd_normal01(){
    double u1, u2, z, result;
    do {
        u1 = rnd_uniform01();
        u2 = rnd_uniform01();
        z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        result = 0.5 + 0.15 * z;      // Shift mean to 0.5 and scale
    } while (result < 0 || result >= 1);
    return result;
}

void init_particles(long seed, double side, long ncside, long long n_part, parcell *par){   
    double (*rnd01)() = rnd_uniform01;
    long long i;
    particle_t aux;
    double cell_size = (double) side / ncside;

    if(seed < 0) {
        rnd01 = rnd_normal01;
        seed = -seed;
    }

    for(int i=0; i < work_size[rank]; i++){
        par[i].n_particles = 0;
        par[i].size = n_part/(ncside*ncside) +1;
        par[i].par = malloc(par[i].size * sizeof(particle_t));
    }

    init_r4uni(seed);

    for(i = 0; i < n_part; i++) {

        aux.x = rnd01() * side;
        aux.y = rnd01() * side;
        aux.vx = (rnd01() - 0.5) * side / ncside / 5.0;
        aux.vy = (rnd01() - 0.5) * side / ncside / 5.0;
        aux.m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;

        double grid_x_aux =  aux.x / cell_size;
        int grid_x = (int)grid_x_aux;

        double grid_y_aux = aux.y / cell_size;
        int grid_y = (int)grid_y_aux;
        
        int id_aux = grid_x * ncside + grid_y;
        
        if (id_aux >= start_point && id_aux < start_point+work_size[rank]){
            int local_cell = id_aux - start_point;
            par[local_cell].par[par[local_cell].n_particles].id = i;
            par[local_cell].par[par[local_cell].n_particles].x = aux.x; 
            par[local_cell].par[par[local_cell].n_particles].y = aux.y; 
            par[local_cell].par[par[local_cell].n_particles].vx = aux.vx; 
            par[local_cell].par[par[local_cell].n_particles].vy = aux.vy;
            par[local_cell].par[par[local_cell].n_particles].m = aux.m;
            par[local_cell].par[par[local_cell].n_particles].alive = 1;

            par[local_cell].n_particles++;

            if (par[local_cell].n_particles == par[local_cell].size){
                par[local_cell].par = realloc(par[local_cell].par, par[local_cell].size * 2 * sizeof(particle_t));
                par[local_cell].size *= 2;
            }
        }
    }

    
}


void calc_center_mass(center_mass * cm, long long num_particles, parcell* par, long grid_size){
    int grid_x;
    int grid_y;

    for(int i=start_point; i < start_point+work_size[rank] ;i++){
        
        cm[i].M = 0;
        cm[i].X = 0;
        cm[i].Y = 0;

        for(int j=0; j < par[i-start_point].n_particles; j++){
            par[i-start_point].par[j].Fx = 0;
            par[i-start_point].par[j].Fy = 0;

            if ( par[i-start_point].par[j].alive == 0)
                continue;

            cm[i].M +=  par[i-start_point].par[j].m;
            cm[i].X += ( par[i-start_point].par[j].m *  par[i-start_point].par[j].x);
            cm[i].Y += ( par[i-start_point].par[j].m *  par[i-start_point].par[j].y);
        
        }

        if(cm[i].M == 0){
            cm[i].X = 0;
            cm[i].Y = 0;
            continue;
        }

        cm[i].X /= cm[i].M;
        cm[i].Y /= cm[i].M;

    }
    center_mass* send;
    send = &cm[0];
    int aux_start=0;
    
    for (int i=0; i<psize ; i++){
        MPI_Bcast(&cm[aux_start], work_size[i], MPI_CENTER_MASS, i, MPI_COMM_WORLD);
        aux_start += work_size[i];

    }
    
}

void cell_calculation(parcell* st_par, long grid_size, double space_size, int t){
    double cell_size = (double)space_size / grid_size;
    double x;
    double y;
    int n_prev = 0, n_next = 0;
    parcell to_send_prev;
    parcell to_send_next;
    parcell to_send_prev_prev;
    parcell to_send_next_next;
    

    to_send_next.par = malloc(work_size[rank] * sizeof(particle_t));
    to_send_next.size = work_size[rank];
    to_send_next.n_particles = 0;
    to_send_prev.par = malloc(work_size[rank] * sizeof(particle_t));
    to_send_prev.size = work_size[rank];
    to_send_prev.n_particles = 0;
    to_send_next_next.par = malloc(2 * sizeof(particle_t));
    to_send_next_next.size = 2;
    to_send_next_next.n_particles = 0;
    to_send_prev_prev.par = malloc(2 * sizeof(particle_t));
    to_send_prev_prev.size = 2;
    to_send_prev_prev.n_particles = 0;

    int recv_count = 0;
    int send_count = 0;

    int prev_rank = 0;
    int next_rank = 0;

    int recv_count_count = 0;
    int send_count_count = 0;

    int prev_prev_rank = 0;
    int next_next_rank = 0;

    if(rank==0){
        prev_rank = psize -1;
        next_rank = rank + 1;
        if (psize > 2){
            prev_prev_rank = psize - 2;
            next_next_rank = rank + 2; 
        }else{
            prev_prev_rank = rank;
            next_next_rank = rank;
        }
        
    }
    else if(rank == psize-1){
        next_rank = 0;
        prev_rank = rank -1;
        if (psize > 2){
            prev_prev_rank = rank - 2;
            next_next_rank = 1;
        }else{
            prev_prev_rank = rank;
            next_next_rank = rank;
        }
    }
    else if(rank==1){
        prev_rank = 0;
        next_rank = rank + 1;
        if (psize > 2){
            prev_prev_rank = psize - 1;
            next_next_rank = rank + 2; 
        }else{
            prev_prev_rank = rank;
            next_next_rank = rank;
        }
    }
    else if(rank == psize-2){
        next_rank = psize-1;
        prev_rank = rank -1;
        if (psize > 2){
            prev_prev_rank = rank - 2;
            next_next_rank = 0;
        }else{
            prev_prev_rank = rank;
            next_next_rank = rank;
        }
    }
    else{
        next_rank = rank + 1;
        prev_rank = rank -1;
        if (psize > 2){
            prev_prev_rank = rank - 2;
            next_next_rank = rank + 2;
        }else{
            prev_prev_rank = rank;
            next_next_rank = rank;
        }
    }

    for(int cell = 0; cell < work_size[rank]; cell++){
        for (int id_par=0; id_par < st_par[cell].n_particles; id_par++){

            if(st_par[cell].par[id_par].alive != 1){
                continue;
            }

            x = st_par[cell].par[id_par].x;
            y = st_par[cell].par[id_par].y;

            if(x<0){
                x += space_size;
                st_par[cell].par[id_par].x += space_size;
            }
                
            if(y<0){
                y += space_size;
                st_par[cell].par[id_par].y += space_size;
            }
        
            if(x > space_size){
                x -= space_size;
                st_par[cell].par[id_par].x -= space_size;
            }
            
            if(y > space_size){
                y -= space_size;
                st_par[cell].par[id_par].y -= space_size;
            }
        
            double grid_x_aux = x / cell_size;
            int grid_x = (int)grid_x_aux;
        
            double grid_y_aux = y / cell_size;
            int grid_y = (int)grid_y_aux;
        
            int new_cell = (grid_x * grid_size + grid_y) - start_point;
            
            if (new_cell != cell){
            
                if(rank == 0 && new_cell + start_point >= grid_size*grid_size - work_size[psize-1]){

                    to_send_prev.par[to_send_prev.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_prev.par[to_send_prev.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_prev.par[to_send_prev.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_prev.par[to_send_prev.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_prev.par[to_send_prev.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_prev.par[to_send_prev.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_prev.par[to_send_prev.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_prev.n_particles ++;
                    
                    if(to_send_prev.n_particles == to_send_prev.size){
                        to_send_prev.par = realloc(to_send_prev.par, to_send_prev.size * 2 * sizeof(particle_t));
                        to_send_prev.size *= 2;
                    }

                }else if(rank == 0 && new_cell + start_point >= grid_size*grid_size-work_size[psize-1]-work_size[psize-2]){

                    to_send_prev_prev.par[to_send_prev_prev.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_prev_prev.n_particles ++;
                    
                    if(to_send_prev_prev.n_particles == to_send_prev_prev.size){
                        to_send_prev_prev.par = realloc(to_send_prev_prev.par, to_send_prev_prev.size * 2 * sizeof(particle_t));
                        to_send_prev_prev.size *= 2;
                    }
                    

////////////////////////////////////////////////////////////////////////////////////////////////////////
                }else if(rank == psize-1 && new_cell + start_point < work_size[0]){

                    to_send_next.par[to_send_next.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_next.par[to_send_next.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_next.par[to_send_next.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_next.par[to_send_next.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_next.par[to_send_next.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_next.par[to_send_next.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_next.par[to_send_next.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_next.n_particles ++;

                    
                    if(to_send_next.n_particles == to_send_next.size){
                        to_send_next.par = realloc(to_send_next.par, to_send_next.size * 2 * sizeof(particle_t));
                        to_send_next.size *= 2;
                    }
                }else if(rank == psize - 1 && new_cell + start_point < work_size[0] + work_size[1]){

                    to_send_next_next.par[to_send_next_next.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_next_next.par[to_send_next_next.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_next_next.par[to_send_next_next.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_next_next.par[to_send_next_next.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_next_next.par[to_send_next_next.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_next_next.par[to_send_next_next.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_next_next.par[to_send_next_next.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_next_next.n_particles ++;

                    
                    if(to_send_next_next.n_particles == to_send_next_next.size){
                        to_send_next_next.par = realloc(to_send_next_next.par, to_send_next_next.size * 2 * sizeof(particle_t));
                        to_send_next_next.size *= 2;
                    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////// 

                }else if(new_cell < -work_size[rank - 1]-1){

                    to_send_prev_prev.par[to_send_prev_prev.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_prev_prev.par[to_send_prev_prev.n_particles].alive = st_par[cell].par[id_par].alive;


                    to_send_prev_prev.n_particles ++;
                    
                    if(to_send_prev_prev.n_particles == to_send_prev_prev.size){
                        to_send_prev_prev.par = realloc(to_send_prev_prev.par, to_send_prev_prev.size * 2 * sizeof(particle_t));
                        to_send_prev_prev.size *= 2;
                    }

                }else if (new_cell < 0){

                    to_send_prev.par[to_send_prev.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_prev.par[to_send_prev.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_prev.par[to_send_prev.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_prev.par[to_send_prev.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_prev.par[to_send_prev.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_prev.par[to_send_prev.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_prev.par[to_send_prev.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_prev.n_particles ++;
                    
                    if(to_send_prev.n_particles == to_send_prev.size){
                        to_send_prev.par = realloc(to_send_prev.par, to_send_prev.size * 2 * sizeof(particle_t));
                        to_send_prev.size *= 2;
                    }
////////////////////////////////////////////////////////////////////////////////////////////////

                }else if(new_cell > work_size[rank] + work_size[rank + 1]-1){

                    to_send_next_next.par[to_send_next_next.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_next_next.par[to_send_next_next.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_next_next.par[to_send_next_next.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_next_next.par[to_send_next_next.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_next_next.par[to_send_next_next.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_next_next.par[to_send_next_next.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_next_next.par[to_send_next_next.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_next_next.n_particles ++;
                    
                    if(to_send_next_next.n_particles == to_send_next_next.size){
                        to_send_next_next.par = realloc(to_send_next_next.par, to_send_next_next.size * 2 * sizeof(particle_t));
                        to_send_next_next.size *= 2;
                    }


                }else if (new_cell >= work_size[rank]){

                    to_send_next.par[to_send_next.n_particles].id = st_par[cell].par[id_par].id;
                    to_send_next.par[to_send_next.n_particles].x = st_par[cell].par[id_par].x;
                    to_send_next.par[to_send_next.n_particles].y = st_par[cell].par[id_par].y;
                    to_send_next.par[to_send_next.n_particles].vx = st_par[cell].par[id_par].vx;
                    to_send_next.par[to_send_next.n_particles].vy = st_par[cell].par[id_par].vy;
                    to_send_next.par[to_send_next.n_particles].m = st_par[cell].par[id_par].m;
                    to_send_next.par[to_send_next.n_particles].alive = st_par[cell].par[id_par].alive;

                    to_send_next.n_particles ++;
                    
                    if(to_send_next.n_particles == to_send_next.size){
                        to_send_next.par = realloc(to_send_next.par, to_send_next.size * 2 * sizeof(particle_t));
                        to_send_next.size *= 2;
                    }

                }else{

                    st_par[new_cell].par[st_par[new_cell].n_particles].id = st_par[cell].par[id_par].id;
                    st_par[new_cell].par[st_par[new_cell].n_particles].x = st_par[cell].par[id_par].x;
                    st_par[new_cell].par[st_par[new_cell].n_particles].y = st_par[cell].par[id_par].y;
                    st_par[new_cell].par[st_par[new_cell].n_particles].vx = st_par[cell].par[id_par].vx;
                    st_par[new_cell].par[st_par[new_cell].n_particles].vy = st_par[cell].par[id_par].vy;
                    st_par[new_cell].par[st_par[new_cell].n_particles].m = st_par[cell].par[id_par].m;
                    st_par[new_cell].par[st_par[new_cell].n_particles].alive = st_par[cell].par[id_par].alive;

                    st_par[new_cell].n_particles++;
                    

                    if(st_par[new_cell].n_particles == st_par[new_cell].size){
                        st_par[new_cell].par = realloc(st_par[new_cell].par, st_par[new_cell].size * 2 * sizeof(particle_t));
                        st_par[new_cell].size *= 2;
                    }
                }
        
                if (id_par != st_par[cell].n_particles-1){
                    st_par[cell].par[id_par] = st_par[cell].par[st_par[cell].n_particles-1]; 
                    id_par--;
                }
                st_par[cell].n_particles--;
            }
        }
        
    }

    particle_t *rcv_prev_par = NULL;
    particle_t *rcv_next_par = NULL;
    particle_t *rcv_prev_prev_par = NULL;
    particle_t *rcv_next_next_par = NULL;


    MPI_Request recv_requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    MPI_Request num_send_requests[2], num_recv_requests[2];
    MPI_Status num_statuses[2];
    MPI_Request send_requests[2];

    MPI_Request recv_requests_requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    MPI_Request num_send_requests_requests[2], num_recv_requests_requests[2];
    MPI_Status num_statuses_statuses[2];
    MPI_Request send_requests_requests[2];

    int incoming_prev_prev_count = 0, incoming_next_next_count = 0;
    int incoming_prev_count = 0, incoming_next_count = 0;
    MPI_Irecv(&incoming_prev_count, 1, MPI_INT, prev_rank, 1, MPI_COMM_WORLD, &num_recv_requests[0]);
    MPI_Irecv(&incoming_next_count, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD, &num_recv_requests[1]);
    
    MPI_Irecv(&incoming_prev_prev_count, 1, MPI_INT, prev_prev_rank, 4, MPI_COMM_WORLD, &num_recv_requests_requests[0]);
    MPI_Irecv(&incoming_next_next_count, 1, MPI_INT, next_next_rank, 5, MPI_COMM_WORLD, &num_recv_requests_requests[1]);
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    int prev_count = to_send_prev.n_particles;
    int next_count = to_send_next.n_particles;

    int prev_prev_count = to_send_prev_prev.n_particles;
    int next_next_count = to_send_next_next.n_particles;    
        
    MPI_Isend(&prev_count, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &num_send_requests[0]);
    MPI_Isend(&next_count, 1, MPI_INT, next_rank, 1, MPI_COMM_WORLD, &num_send_requests[1]);

    MPI_Isend(&prev_prev_count, 1, MPI_INT, prev_prev_rank, 5, MPI_COMM_WORLD, &num_send_requests_requests[0]);
    MPI_Isend(&next_next_count, 1, MPI_INT, next_next_rank, 4, MPI_COMM_WORLD, &num_send_requests_requests[1]);

    /////////////////////////////////////////////////////////////////////////////////////////
    MPI_Waitall(2, num_recv_requests, num_statuses);
    MPI_Waitall(2, num_send_requests, MPI_STATUSES_IGNORE);

    MPI_Waitall(2, num_recv_requests_requests, num_statuses_statuses);
    MPI_Waitall(2, num_send_requests_requests, MPI_STATUSES_IGNORE);
    ///////////////////////////////////////////////////////////////////////////////////////
    // Allocate receive buffers
    if (incoming_prev_count > 0) {
        recv_count++;
        rcv_prev_par = malloc(incoming_prev_count * sizeof(particle_t));
        MPI_Irecv(rcv_prev_par, incoming_prev_count, MPI_PARTICLE_T, prev_rank, 2, MPI_COMM_WORLD, &recv_requests[0]);  
    }    
    if (incoming_next_count > 0) {
        recv_count++;
        rcv_next_par = malloc(incoming_next_count * sizeof(particle_t));
        MPI_Irecv(rcv_next_par, incoming_next_count, MPI_PARTICLE_T,next_rank, 3, MPI_COMM_WORLD, &recv_requests[1]);
    }

    if (incoming_prev_prev_count > 0) {
        recv_count_count++;
        rcv_prev_prev_par = malloc(incoming_prev_prev_count * sizeof(particle_t));
        MPI_Irecv(rcv_prev_prev_par, incoming_prev_prev_count, MPI_PARTICLE_T, prev_prev_rank, 6, MPI_COMM_WORLD, &recv_requests_requests[0]);  
    }    
    if (incoming_next_next_count > 0) {
        recv_count_count++;
        rcv_next_next_par = malloc(incoming_next_next_count * sizeof(particle_t));
        MPI_Irecv(rcv_next_next_par, incoming_next_next_count, MPI_PARTICLE_T,next_next_rank, 7, MPI_COMM_WORLD, &recv_requests_requests[1]);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (prev_count > 0) {
        send_count++;
        
        MPI_Isend(to_send_prev.par, prev_count, MPI_PARTICLE_T, prev_rank, 3, MPI_COMM_WORLD, &send_requests[0]);
    }
    
    if (next_count > 0) {
        send_count++;
        
        MPI_Isend(to_send_next.par, next_count, MPI_PARTICLE_T, next_rank, 2, MPI_COMM_WORLD, &send_requests[1]);
    }

    if (prev_prev_count > 0) {
        send_count_count++;
        
        MPI_Isend(to_send_prev_prev.par, prev_prev_count, MPI_PARTICLE_T, prev_prev_rank, 7, MPI_COMM_WORLD, &send_requests_requests[0]);
    }
    
    if (next_next_count > 0) {
        send_count_count++;
        
        MPI_Isend(to_send_next_next.par, next_next_count, MPI_PARTICLE_T, next_next_rank, 6, MPI_COMM_WORLD, &send_requests_requests[1]);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (recv_count > 0) {

        if (incoming_prev_count > 0) {
            MPI_Wait(&recv_requests[0], MPI_STATUS_IGNORE);
        }
        if (incoming_next_count > 0) {
            MPI_Wait(&recv_requests[1], MPI_STATUS_IGNORE);
        }

    }

    if (recv_count_count > 0) {

        if (incoming_prev_prev_count > 0) {
            MPI_Wait(&recv_requests_requests[0], MPI_STATUS_IGNORE);
        }
        if (incoming_next_next_count > 0) {
            MPI_Wait(&recv_requests_requests[1], MPI_STATUS_IGNORE);
        }

    }

    ////////////////////////////////////////////////////////

    if (incoming_prev_count > 0){
        for (int i = 0; i < incoming_prev_count; i++) {
            x = rcv_prev_par[i].x;
            y = rcv_prev_par[i].y;
            
            double grid_x_aux = x / cell_size;
            int grid_x = (int)grid_x_aux;
        
            double grid_y_aux = y / cell_size;
            int grid_y = (int)grid_y_aux;
        
            int new_cell = (grid_x * grid_size + grid_y) - start_point;

            st_par[new_cell].par[st_par[new_cell].n_particles].id = rcv_prev_par[i].id;
            st_par[new_cell].par[st_par[new_cell].n_particles].x = rcv_prev_par[i].x;
            st_par[new_cell].par[st_par[new_cell].n_particles].y = rcv_prev_par[i].y;
            st_par[new_cell].par[st_par[new_cell].n_particles].vx = rcv_prev_par[i].vx;
            st_par[new_cell].par[st_par[new_cell].n_particles].vy = rcv_prev_par[i].vy;
            st_par[new_cell].par[st_par[new_cell].n_particles].m = rcv_prev_par[i].m;
            st_par[new_cell].par[st_par[new_cell].n_particles].alive = rcv_prev_par[i].alive;
            st_par[new_cell].n_particles++;
            
            if (st_par[new_cell].n_particles >= st_par[new_cell].size) {
                particle_t *new_par = realloc(st_par[new_cell].par, st_par[new_cell].size * 2 * sizeof(particle_t));
                st_par[new_cell].par = new_par;
                st_par[new_cell].size *= 2;
            }
        }
    }


    if (incoming_prev_prev_count > 0){
        for (int i = 0; i < incoming_prev_prev_count; i++) {
            x = rcv_prev_prev_par[i].x;
            y = rcv_prev_prev_par[i].y;
            
            double grid_x_aux = x / cell_size;
            int grid_x = (int)grid_x_aux;
        
            double grid_y_aux = y / cell_size;
            int grid_y = (int)grid_y_aux;
        
            int new_cell = (grid_x * grid_size + grid_y) - start_point;

            st_par[new_cell].par[st_par[new_cell].n_particles].id = rcv_prev_prev_par[i].id;
            st_par[new_cell].par[st_par[new_cell].n_particles].x = rcv_prev_prev_par[i].x;
            st_par[new_cell].par[st_par[new_cell].n_particles].y = rcv_prev_prev_par[i].y;
            st_par[new_cell].par[st_par[new_cell].n_particles].vx = rcv_prev_prev_par[i].vx;
            st_par[new_cell].par[st_par[new_cell].n_particles].vy = rcv_prev_prev_par[i].vy;
            st_par[new_cell].par[st_par[new_cell].n_particles].m = rcv_prev_prev_par[i].m;
            st_par[new_cell].par[st_par[new_cell].n_particles].alive = rcv_prev_prev_par[i].alive;
            
            st_par[new_cell].n_particles++;
            
            if (st_par[new_cell].n_particles >= st_par[new_cell].size) {
                particle_t *new_par = realloc(st_par[new_cell].par, st_par[new_cell].size * 2 * sizeof(particle_t));
                st_par[new_cell].par = new_par;
                st_par[new_cell].size *= 2;
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (incoming_next_count > 0) {
        for (int i = 0; i < incoming_next_count; i++) {
            x = rcv_next_par[i].x;
            y = rcv_next_par[i].y;
            
            double grid_x_aux = x / cell_size;
            int grid_x = (int)grid_x_aux;
        
            double grid_y_aux = y / cell_size;
            int grid_y = (int)grid_y_aux;
        
            int new_cell = (grid_x * grid_size + grid_y) - start_point;

            st_par[new_cell].par[st_par[new_cell].n_particles].id = rcv_next_par[i].id;
            st_par[new_cell].par[st_par[new_cell].n_particles].x = rcv_next_par[i].x;
            st_par[new_cell].par[st_par[new_cell].n_particles].y = rcv_next_par[i].y;
            st_par[new_cell].par[st_par[new_cell].n_particles].vx = rcv_next_par[i].vx;
            st_par[new_cell].par[st_par[new_cell].n_particles].vy = rcv_next_par[i].vy;
            st_par[new_cell].par[st_par[new_cell].n_particles].m = rcv_next_par[i].m;
            st_par[new_cell].par[st_par[new_cell].n_particles].alive = rcv_next_par[i].alive;
            
            st_par[new_cell].n_particles++;

            if (st_par[new_cell].n_particles >= st_par[new_cell].size) {
                particle_t *new_par = realloc(st_par[new_cell].par, st_par[new_cell].size * 2 * sizeof(particle_t));
                st_par[new_cell].par = new_par;
                st_par[new_cell].size *= 2;
            }
        }
    }

    if (incoming_next_next_count > 0) {
        for (int i = 0; i < incoming_next_next_count; i++) {
            x = rcv_next_next_par[i].x;
            y = rcv_next_next_par[i].y;
            
            double grid_x_aux = x / cell_size;
            int grid_x = (int)grid_x_aux;
        
            double grid_y_aux = y / cell_size;
            int grid_y = (int)grid_y_aux;
        
            int new_cell = (grid_x * grid_size + grid_y) - start_point;
                        
            st_par[new_cell].par[st_par[new_cell].n_particles].id = rcv_next_next_par[i].id;
            st_par[new_cell].par[st_par[new_cell].n_particles].x = rcv_next_next_par[i].x;
            st_par[new_cell].par[st_par[new_cell].n_particles].y = rcv_next_next_par[i].y;
            st_par[new_cell].par[st_par[new_cell].n_particles].vx = rcv_next_next_par[i].vx;
            st_par[new_cell].par[st_par[new_cell].n_particles].vy = rcv_next_next_par[i].vy;
            st_par[new_cell].par[st_par[new_cell].n_particles].m = rcv_next_next_par[i].m;
            st_par[new_cell].par[st_par[new_cell].n_particles].alive = rcv_next_next_par[i].alive;
            st_par[new_cell].n_particles++;

            if (st_par[new_cell].n_particles >= st_par[new_cell].size) {
                particle_t *new_par = realloc(st_par[new_cell].par, st_par[new_cell].size * 2 * sizeof(particle_t));
                st_par[new_cell].par = new_par;
                st_par[new_cell].size *= 2;
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Wait for sends to complete (if any)
    if (prev_count > 0 && send_requests[0] != MPI_REQUEST_NULL) {
        MPI_Wait(&send_requests[0], MPI_STATUS_IGNORE);
    }
    
    if (next_count > 0 && send_requests[1] != MPI_REQUEST_NULL) {
        MPI_Wait(&send_requests[1], MPI_STATUS_IGNORE);
    }

    // Wait for sends to complete (if any)
    if (prev_prev_count > 0 && send_requests_requests[0] != MPI_REQUEST_NULL) {
        MPI_Wait(&send_requests_requests[0], MPI_STATUS_IGNORE);
    }
    
    if (next_next_count > 0 && send_requests_requests[1] != MPI_REQUEST_NULL) {
        MPI_Wait(&send_requests_requests[1], MPI_STATUS_IGNORE);
    }
    
    free(rcv_next_par);
    free(to_send_next.par);
    free(to_send_prev.par);
    free(rcv_prev_par);

    free(rcv_next_next_par);
    free(to_send_next_next.par);
    free(to_send_prev_prev.par);
    free(rcv_prev_prev_par);
}


int simulation(center_mass *cells, double space_size, long grid_size, long long num_particles, long long num_timesteps, parcell *st_par){

    double delta_x = 0, delta_y = 0; //displacement of the particle in x and y
    int collision_count = 0; //count collisions
    double cell_size = (double)space_size / grid_size;
    ParColision *colision ;
    //ParColision colision[num_particles];
    particle_t *px,*py;

    for(int t = 0; t < num_timesteps; t++){
        
        calc_center_mass(cells, num_particles, st_par, grid_size);

        for(int j=0; j< work_size[rank] ;j++){

            for(int k = 0; k<st_par[j].n_particles; k++){  

                px = &st_par[j].par[k];

                if(px->alive == 1){

                    for(int l=k+1; l<st_par[j].n_particles; l++){

                        py = &st_par[j].par[l];

                        if (py->alive == 1){

                            delta_x = py->x - px->x;
                            delta_y = py->y - px->y;
                            double distance2 = delta_x * delta_x + delta_y * delta_y;
                            double distance = sqrt(distance2); 
                            double distance3 = distance2 * distance;
                            
                            double force = G * (py->m * px->m) / distance3;

                            double fx = force * delta_x;
                            double fy = force * delta_y;

                            px->Fx += fx;
                            px->Fy += fy;

                            py->Fx -= fx;
                            py->Fy -= fy;

                        }
                    }

                    for(int w =-1; w<2; w++){
                        for(int y =-1; y<2; y++){
                            if(y==0 && w== 0){  // Skip the cell of the particle
                                continue;
                            }
                            int cell_x = (j + start_point) / grid_size;
                            int cell_y = (j + start_point) % grid_size;

                            int cellad_x = cell_x + w;
                            int cellad_y = cell_y + y;

                            double x_aux, y_aux;

                            if(cellad_x<0){
                                cellad_x = (int) grid_size-1;
                            }else if(cellad_x >= (int) grid_size){
                                cellad_x = 0;
                            }

                            if(cellad_y<0){
                                cellad_y = (int) grid_size-1;
                            }else if(cellad_y >= (int)grid_size){
                                cellad_y = 0;
                            }

                            if(cell_x - cellad_x < -1){
                                x_aux = cells[cellad_x*grid_size + cellad_y].X - space_size;
                            }else if(cell_x - cellad_x > 1){
                                x_aux = cells[cellad_x*grid_size + cellad_y].X + space_size;
                            }else{
                                x_aux = cells[cellad_x*grid_size + cellad_y].X;
                            }

                            if(cell_y - cellad_y < -1){
                                y_aux = cells[cellad_x*grid_size + cellad_y].Y - space_size;
                            }else if(cell_y - cellad_y > 1){
                                y_aux = cells[cellad_x*grid_size + cellad_y].Y + space_size;
                            }else{
                                y_aux = cells[cellad_x*grid_size + cellad_y].Y;
                            }

                            delta_x = x_aux - px->x;
                            delta_y = y_aux - px->y;

                            double distance2_cm = delta_x * delta_x + delta_y * delta_y;  // r^2
                            double distance_cm = sqrt(distance2_cm);
                            double distance3_cm = distance2_cm * distance_cm;

                            double force = G * (cells[cellad_x*grid_size + cellad_y].M * px->m) / distance3_cm;  // G * M * m / r^3

                            px->Fx += force * delta_x;
                            px->Fy += force * delta_y;
                        } 
                    } 

                    px->ax = px->Fx/px->m; 
                    px->ay = px->Fy/px->m;
                    px->x = px->x + px->vx*DELTAT + 0.5*px->ax*DELTAT*DELTAT; //calculate new position x
                    px->y = px->y + px->vy*DELTAT + 0.5*px->ay*DELTAT*DELTAT; //calculate new position y
                    px->vx = px->vx + px->ax*DELTAT; //calculate new velocity along x
                    px->vy = px->vy + px->ay*DELTAT; //calculate new velocity along y
                }

            }
        }

        cell_calculation(st_par, grid_size, space_size, t);
        int aux_size = (num_particles / grid_size/ grid_size) * 0.1;
        colision = malloc(aux_size * sizeof(ParColision));
        int aux = 0;
        for(int j=0; j<work_size[rank]; j++){
            for (int idx_a=0; idx_a < st_par[j].n_particles; idx_a++){
                for (int idx_b=idx_a+1; idx_b < st_par[j].n_particles; idx_b++){
                    if ( st_par[j].par[idx_a].alive == 0 || st_par[j].par[idx_b].alive == 0)
                        continue;
                        
                    delta_x = st_par[j].par[idx_a].x - st_par[j].par[idx_b].x;
                    delta_y = st_par[j].par[idx_a].y - st_par[j].par[idx_b].y;
                    double distance2 = delta_x * delta_x + delta_y * delta_y;

                    if(distance2 <= EPSILON2){
                        colision[aux].a = idx_a;
                        colision[aux].b = idx_b;
                        colision[aux].cell = j;
                        aux++;
                        collision_count++;

                        if (aux == aux_size){
                            colision = realloc (colision, aux_size*2*sizeof(ParColision));
                            aux_size *= 2; 
                        }
                        
                    }
                }
            }
        }
        
        for(int n=0; n<aux; n++){
            st_par[colision[n].cell].par[colision[n].a].alive = 0;
            st_par[colision[n].cell].par[colision[n].b].alive = 0;
            for(int m=n+1; m<aux; m++){
                if((colision[n].cell == colision[m].cell) && (colision[n].b == colision[m].a || colision[n].b == colision[m].b || colision[n].a == colision[m].a)){
                    collision_count--;
                    continue;
                }
            }
        }
        free(colision);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    
    return collision_count;
}

void print_result(parcell* st_par, int local_collisions,double exec_time){

    int total_collisions; 
    double positions[2];

    for (int i=0; i<work_size[rank]; i++){
        for(int j=0; j<st_par[i].n_particles; j++){
            if (st_par[i].par[j].id == 0){
                positions[0] = st_par[i].par[j].x;
                positions[1] = st_par[i].par[j].y;
                // Envia para o rank 0
                MPI_Send(positions, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
                //fprintf(stdout, "%.3f %.3f\n", st_par[i].par[j].x, st_par[i].par[j].y);
        }
    }
    MPI_Reduce(&local_collisions, &total_collisions, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        MPI_Status status;
        MPI_Recv(positions, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        fprintf(stdout, "%.3f %.3f\n", positions[0], positions[1]);
        fprintf(stdout, "%d\n", total_collisions);
        fprintf(stderr, "%.1fs\n", exec_time); 
        
    }
}

int main(int argc, char *argv[]){
    //read args values
    double exec_time;
    int colisions;
    long sseed;
    
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <seed> <space_size> <grid_size> <num_particles> <num_timesteps>\n", argv[0]);
        return EXIT_FAILURE;
    }

    sseed = atoi(argv[1]);
    double space_size = atof(argv[2]);
    long grid_size = atol(argv[3]);
    long long num_particles = atoll(argv[4]);
    long long num_timesteps = atoll(argv[5]);

    parcell* particles;

    // Inicialization of MPI and types to send
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize); 
    
    create_mpi_center_mass_type();
    create_mpi_particle_t();

    int aux_size = (grid_size*grid_size)/psize;
    int remain = (grid_size*grid_size)%psize;

    // Compute the number of cells in the process
    work_size = malloc(psize*sizeof(long long));
    for (int i=0; i < psize; i++){
        work_size[i] = aux_size;
        if (i < remain){
            work_size[i] ++;
        }
    }

    start_point=0;
    for (int i=0; i<rank; i++){
        start_point += work_size[i];
    }
    
    particles = malloc(work_size[rank] * sizeof(parcell));
    center_mass* cells = malloc((grid_size*grid_size) * sizeof(center_mass)); 
    
    init_particles(sseed, space_size, grid_size, num_particles, particles);

    exec_time = -omp_get_wtime();
    int local_colisions = simulation(cells, space_size, grid_size, num_particles, num_timesteps, particles);
    exec_time += omp_get_wtime();
    print_result(particles, local_colisions,exec_time);
    

    MPI_Finalize();

    free(particles);
    free(cells);
    free(work_size);
}