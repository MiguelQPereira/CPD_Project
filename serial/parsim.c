#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define G 6.67408e-11
#define EPSILON2 (0.005*0.005)
#define DELTAT 0.1

unsigned int seed;

typedef struct particle_{

    double x; //x coordinate
    double y; //y coordinate
    double vx; //velocity along x axis
    double vy; //velocity along y axis
    double F; //resultant of the forces
    double a; //particle acceleration
    double m; //mass of the particle
    
}particle_t;

typedef struct cm_{
        double X; //X center of mass
        double Y; //Y of center of mass
        double M; //sum of the mass of each particle
}center_mass;

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
    }
}

void calc_center_mass(center_mass** cm, long long num_particles,particle_t* par, double space_size, long grid_size){
    int grid_x;
    int grid_y;
    
    // Set M of each cell to 0
    for(int i=0; i<grid_size;i++){
        for(int j=0; j<grid_size; j++){
            cm[i][j].M = 0;
        }
    }

    // Compute the M of each cell
    for(int i; i< num_particles; i++){
        grid_x = par[i].x / (space_size/grid_size);
        grid_y = par[i].y / (space_size/grid_size);

        cm[grid_x][grid_y].M += par[i].m;
    }

    for(int i; i< num_particles; i++){

        
        //Depois acabo
    }
  
}

void simulation(double space_size, long grid_size, long long num_particles, double num_timesteps, particle_t *par){
    
    center_mass cm[grid_size][grid_size];
    int grid_x;
    int grid_y;
    double distance;
    int position_x;
    int position_y;

    
    for(int i = 0; i < num_timesteps; i++){

        calc_center_mass(cm, num_particles, par, space_size, grid_size); // Compute the center  mass at the current instant for every cell
        
        for(int j = 0; j < num_particles; j++){
            grid_x = par[j].x / (space_size/grid_size);
            grid_y = par[j].y / (space_size/grid_size);
            for(int k; k < num_particles; k++){
                if(k == j){
                    continue;
                }
                position_x = par[k].x / (space_size/grid_size);
                position_y = par[k].y / (space_size/grid_size);

                if(position_x == grid_x && position_y == grid_y){
                   distance = sqrt((par[j].x - par[k].x))  
                }

            }
            
        }
    }

}

int main(int argc, char *argv[])
{
    //read args values
    double exec_time;
    
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <seed> <space_size> <grid_size> <num_particles> <num_timesteps>\n", argv[0]);
        return EXIT_FAILURE;
    }

    seed = atoi(argv[1]);
    double space_size = atof(argv[2]);
    long grid_size = atol(argv[3]);
    long long num_particles = atoll(argv[4]);
    double num_timesteps = atof(argv[5]);
    
    particle_t* particles = malloc(num_particles * sizeof(particle_t));
    
    init_particles(seed, space_size, grid_size, num_particles, particles);
    exec_time = -omp_get_wtime();
    simulation(space_size, grid_size, num_particles, num_timesteps, particles);
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result();
    // to the stdout!
}