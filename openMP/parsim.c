#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#define G 6.67408e-11
#define EPSILON2 (0.005*0.005)
#define DELTAT 0.1

unsigned int seed;

typedef struct particle_{

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

typedef struct cm_{
        double X; //X center of mass
        double Y; //Y of center of mass
        double M; //sum of the mass of each particle
        int * par_index; //Index of particles in cell
        int n_par; //Number of particles in cell
        int cap;
}center_mass;

typedef struct {
    int a;
    int b;
} ParColision;

void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

double rnd_uniform01()
{
    int seed_in = seed;
    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);
    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

double rnd_normal01()
{
    double u1, u2, z, result;
    do {
        u1 = rnd_uniform01();
        u2 = rnd_uniform01();
        z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        result = 0.5 + 0.15 * z;      // Shift mean to 0.5 and scale
    } while (result < 0 || result >= 1);
    return result;
}

void init_particles(long seed, double side, long ncside, long long n_part, particle_t *par)
{   
    double (*rnd01)() = rnd_uniform01;
    long long i;

    if(seed < 0) {
        rnd01 = rnd_normal01;
        seed = -seed;
    }
    
    init_r4uni(seed);

    for(i = 0; i < n_part; i++) {
        par[i].x = rnd01() * side;
        par[i].y = rnd01() * side;
        par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
        par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;

        par[i].m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
        par[i].alive = 1; //all particles begin alive
    }
}

void calc_center_mass(center_mass ** cm, long long num_particles, particle_t* par, double cell_size, long grid_size){
    int grid_x;
    int grid_y;
    
    // Set M of each cell to 0
    #pragma omp parallel
    #pragma omp for schedule(dynamic, 10) collapse(2) nowait
        for(int i=0; i<grid_size;i++){
            for(int j=0; j<grid_size; j++){
                cm[i][j].M = 0;
                cm[i][j].X = 0;
                cm[i][j].Y = 0;
            }
        }

    // Compute the M of each cell
    #pragma omp parallel
    #pragma omp for schedule(dynamic, 10) private(grid_x, grid_y) nowait
        for(int i = 0; i< num_particles; i++){

            par[i].Fx = 0;
            par[i].Fy = 0;

            if (par[i].alive == 0)
                continue;

            double grid_x_aux = par[i].x / cell_size;
            grid_x = (int) grid_x_aux; 
        
            double grid_y_aux = par[i].y / cell_size;
            grid_y = (int) grid_y_aux;

            #pragma omp critical
            {
                cm[grid_x][grid_y].M += par[i].m;
                cm[grid_x][grid_y].X += (par[i].m * par[i].x);
                cm[grid_x][grid_y].Y += (par[i].m * par[i].y);
            }
            
        }

        #pragma omp parallel
        #pragma omp for collapse(2) nowait
        for(int i = 0; i< grid_size; i++){
            for(int j= 0; j< grid_size; j++){
                if(cm[j][i].M == 0){
                    cm[j][i].X = 0;
                    cm[j][i].Y = 0;
                    continue;
                }

                #pragma omp atomic
                cm[j][i].X /= cm[j][i].M;
                #pragma omp atomic
                cm[j][i].Y /= cm[j][i].M;
            }
        }
}

void grid_calculation(center_mass ** cm, long long num_particles, double cell_size, particle_t *par, long grid_size){

        // Allocate and initialize locks for each grid cell
        omp_lock_t **grid_locks = (omp_lock_t **)malloc(grid_size * sizeof(omp_lock_t *));
        for (int i = 0; i < grid_size; i++) {
            grid_locks[i] = (omp_lock_t *)malloc(grid_size * sizeof(omp_lock_t));
            for (int j = 0; j < grid_size; j++) {
                omp_init_lock(&grid_locks[i][j]);
            }
        }
    
        // Initialize grid cells
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                cm[i][j].n_par = 0;
            }
        }
    
        // Assign particles to grid cells
        #pragma omp parallel for
        for (int p = 0; p < num_particles; p++) {
            if (par[p].alive) {
                const int grid_x = (int)(par[p].x / cell_size);
                const int grid_y = (int)(par[p].y / cell_size);
    
                // Handle periodic boundary conditions
                const int wrapped_x = (grid_x + grid_size) % grid_size;
                const int wrapped_y = (grid_y + grid_size) % grid_size;
    
                omp_set_lock(&grid_locks[wrapped_x][wrapped_y]);
    
                // Check if resize needed
                if (cm[wrapped_x][wrapped_y].n_par >= cm[wrapped_x][wrapped_y].cap - 1) {
                    const int new_cap = cm[wrapped_x][wrapped_y].cap * 2;
                    int *temp = realloc(cm[wrapped_x][wrapped_y].par_index, new_cap * sizeof(int));
                    if (!temp) {
                        fprintf(stderr, "Realloc failed for grid cell [%d][%d]\n", wrapped_x, wrapped_y);
                        omp_unset_lock(&grid_locks[wrapped_x][wrapped_y]);
                        continue;
                    }
                    cm[wrapped_x][wrapped_y].par_index = temp;
                    cm[wrapped_x][wrapped_y].cap = new_cap;
                }
    
                // Store particle index
                const int idx = cm[wrapped_x][wrapped_y].n_par++;
                cm[wrapped_x][wrapped_y].par_index[idx] = p;
                cm[wrapped_x][wrapped_y].par_index[idx + 1] = -1;
    
                omp_unset_lock(&grid_locks[wrapped_x][wrapped_y]);
            }
        }
    
        // Cleanup locks
        for (int i = 0; i < grid_size; i++) {
            for (int j = 0; j < grid_size; j++) {
                omp_destroy_lock(&grid_locks[i][j]);
            }
            free(grid_locks[i]);
        }
        free(grid_locks);
    
}

int simulation(double space_size, long grid_size, long long num_particles, long long num_timesteps, particle_t *par){
    
    center_mass ** cm;
    double delta_x = 0, delta_y = 0; //displacement of the particle in x and y
    int collision_count = 0; //count collisions
    double cell_size = (double)space_size / grid_size;
    ParColision colision[num_particles];

    //Allocation of memory to save the center of mass of each cell
    cm = malloc(grid_size * sizeof(center_mass*));
        for (int i=0; i < grid_size; i++){
            cm[i] = malloc(grid_size * sizeof(center_mass));
            for (int j= 0; j<grid_size; j++){
                cm[i][j].par_index = malloc((int)(num_particles/(grid_size*grid_size)) * sizeof(int));
                cm[i][j].cap = (int)(num_particles/(grid_size*grid_size)) ;
            }
        }

    grid_calculation(cm,num_particles, cell_size, par, grid_size);

    //Simulation loop
    for(int i = 0; i < num_timesteps; i++){

        calc_center_mass(cm, num_particles, par, cell_size, grid_size); // Compute the center of mass at the current instant for every cell
        
        #pragma omp parallel 
        #pragma omp for private(delta_x, delta_y) schedule(dynamic, 10) collapse(2) nowait
        for(int idx_x = 0; idx_x < grid_size; idx_x++){
            for(int idx_y = 0; idx_y < grid_size; idx_y++){               

                for(int j = 0; j<cm[idx_x][idx_y].n_par; j++){

                    int px = cm[idx_x][idx_y].par_index[j];

                    if(par[px].alive == 1){
                        for(int k = j+1; k<cm[idx_x][idx_y].n_par; k++){
                            int py = cm[idx_x][idx_y].par_index[k];

                            if (par[py].alive == 1){

                                delta_x = par[py].x - par[px].x;
                                delta_y = par[py].y - par[px].y;
                                double distance2 = delta_x * delta_x + delta_y * delta_y;
                                double distance = sqrt(distance2); 
                                double distance3 = distance2 * distance;
                                
                                double force = G * (par[py].m * par[px].m) / distance3;

                                double fx = force * delta_x;
                                double fy = force * delta_y;

                                par[px].Fx += fx;
                                par[px].Fy += fy;

                                par[py].Fx -= fx;
                                par[py].Fy -= fy;
        
                            }
                        }
                        
                        // Force component from the centers of mass
                        for(int w =-1; w<2; w++){
                            for(int y =-1; y<2; y++){
                                if(y==0 && w== 0){  // Skip the cell of the particle
                                    continue;
                                }

                                int cell_x = idx_x + w;
                                int cell_y = idx_y + y;

                                double x_aux, y_aux;

                                if(cell_x<0){
                                    cell_x = (int) grid_size-1;
                                }else if(cell_x >= (int) grid_size){
                                    cell_x = 0;
                                }

                                if(cell_y<0){
                                    cell_y = (int) grid_size-1;
                                }else if(cell_y >= (int)grid_size){
                                    cell_y = 0;
                                }

                                if(idx_x - cell_x < -1){
                                    x_aux = cm[cell_x][cell_y].X - space_size;
                                }else if(idx_x - cell_x > 1){
                                    x_aux = cm[cell_x][cell_y].X + space_size;
                                }else{
                                    x_aux = cm[cell_x][cell_y].X;
                                }
                                

                                if(idx_y - cell_y < -1){
                                    y_aux = cm[cell_x][cell_y].Y - space_size;
                                }else if(idx_y - cell_y > 1){
                                    y_aux = cm[cell_x][cell_y].Y + space_size;
                                }else{
                                    y_aux = cm[cell_x][cell_y].Y;
                                }

                                delta_x = x_aux - par[px].x;
                                delta_y = y_aux - par[px].y;


                                double distance2_cm = delta_x * delta_x + delta_y * delta_y;  // r^2
                                double distance_cm = sqrt(distance2_cm);
                                double distance3_cm = distance2_cm * distance_cm;

                                double force = G * (cm[cell_x][cell_y].M * par[px].m) / distance3_cm;  // G * M * m / r^3

                                par[px].Fx += force * delta_x;
                                par[px].Fy += force * delta_y;
                                
                            } 
                        } 

                        par[px].ax = par[px].Fx/par[px].m; 
                        par[px].ay = par[px].Fy/par[px].m;
                        par[px].x = par[px].x + par[px].vx*DELTAT + 0.5*par[px].ax*DELTAT*DELTAT; //calculate new position x
                        par[px].y = par[px].y + par[px].vy*DELTAT + 0.5*par[px].ay*DELTAT*DELTAT; //calculate new position y
                        par[px].vx = par[px].vx + par[px].ax*DELTAT; //calculate new velocity along x
                        par[px].vy = par[px].vy + par[px].ay*DELTAT; //calculate new velocity along y  
                        
                        if(par[px].x<0)
                            par[px].x = space_size + par[px].x;
                
                        if(par[px].y<0)
                            par[px].y = space_size + par[px].y;

                        if(par[px].x > space_size)
                            par[px].x = par[px].x - space_size;
                        
                        if(par[px].y > space_size)
                            par[px].y = par[px].y - space_size;

                    }
                
                }
                
            }
        }

        grid_calculation(cm,num_particles, cell_size, par, grid_size);

        #pragma omp parallel 
        #pragma omp for private(delta_x, delta_y) schedule(dynamic, 10) collapse(2) nowait
        for(int k=0; k<grid_size; k++){
            for(int w = 0; w<grid_size; w++){

                    for (int idx_a=0; idx_a < cm[k][w].n_par; idx_a++){
                        for (int idx_b=idx_a+1; idx_b < cm[k][w].n_par; idx_b++){
                            if ( par[cm[k][w].par_index[idx_a]].alive == 0 || par[cm[k][w].par_index[idx_b]].alive == 0)
                                continue;
                                
                            delta_x = par[cm[k][w].par_index[idx_a]].x - par[cm[k][w].par_index[idx_b]].x;
                            delta_y = par[cm[k][w].par_index[idx_a]].y - par[cm[k][w].par_index[idx_b]].y;
                            double distance2 = delta_x * delta_x + delta_y * delta_y;

                            if(distance2 <= EPSILON2){
                                if (collision_count < num_particles) {
                                    colision[collision_count].a = cm[k][w].par_index[idx_a];
                                    colision[collision_count].b = cm[k][w].par_index[idx_b];
                                    #pragma omp atomic
                                    collision_count++;
                                }
                                continue;
                            } 
                        } 
                    }                
            }
        }    

        #pragma omp parallel 
        #pragma omp for schedule(dynamic, 10) nowait
        for(int n=0; n<collision_count; n++){
            par[colision[n].a].alive = 0;
            par[colision[n].b].alive = 0;
            for(int m=n+1; m<collision_count; m++){
                if(colision[n].b == colision[m].a || colision[n].b == colision[m].b || colision[n].a == colision[m].a){
                    #pragma omp atomic
                    collision_count--;
                }
            }
        }
    }
    
    #pragma omp parallel 
    #pragma omp for schedule(dynamic, 10) nowait
    for(int i=0; i < grid_size; i++) {
        free(cm[i]);
        for(int j = 0; j < grid_size; j++) {
            free(cm[i][j].par_index);
        }
    }
    free(cm);

    return collision_count;
}

 void print_result(double x, double y, int collisions){
    fprintf(stdout, "%.3f %.3f\n", x, y);
    fprintf(stdout, "%d\n", collisions);
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

    particle_t* particles = malloc(num_particles * sizeof(particle_t));
    
    init_particles(sseed, space_size, grid_size, num_particles, particles);
    exec_time = -omp_get_wtime();
    
    colisions = simulation(space_size, grid_size, num_particles, num_timesteps, particles);
    exec_time += omp_get_wtime();
    
    print_result(particles[0].x, particles[0].y, colisions);
    fprintf(stderr, "%.1fs\n", exec_time);
    
    free(particles);
}