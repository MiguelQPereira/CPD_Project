// Compilar:
// gcc -o parsim parsim.c -lm -fopenmp
// Ex para executar:
// ./parsim 1 2 3 10 1
// Resultado suposto dar:
// 1.570 0.056
// 0


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
}center_mass;

typedef struct {
    int a;
    int b;
} ParColision;

int colisao_count = 0;

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

void calc_center_mass(center_mass ** cm, long long num_particles,particle_t* par, double space_size, long grid_size){
    int grid_x;
    int grid_y;
    
    double cell_size = (double)space_size / grid_size;

    // Set M of each cell to 0
    for(int i=0; i<grid_size;i++){
        for(int j=0; j<grid_size; j++){
            cm[i][j].M = 0;
        }
    }

    // Compute the M of each cell
    for(int i = 0; i< num_particles; i++){
        printf("p: %f\n", par[i].x);

        double grid_x_aux = par[i].x / cell_size;
        grid_x = (int) grid_x_aux; 
    
        double grid_y_aux = par[i].y/ cell_size;
        grid_y = (int) grid_y_aux;

        printf("%d %d\n", grid_x, grid_y);
        cm[grid_x][grid_y].M += par[i].m;
        cm[grid_x][grid_y].X += (par[i].m * par[i].x);
        cm[grid_x][grid_y].Y += (par[i].m * par[i].y);
    }

    for(int i = 0; i< grid_size; i++){
        for(int j= 0; j< grid_size; j++){
            if(cm[j][i].M == 0){
                cm[j][i].X = 0;
                cm[j][i].Y = 0;
                printf("CM - %f %f %f\n", cm[j][i].X, cm[j][i].Y, cm[j][i].M);
                continue;
            }

            cm[j][i].X /= cm[j][i].M;
            cm[j][i].Y /= cm[j][i].M;

            printf("CM -%d %d %f %f %f\n",i,j, cm[j][i].X, cm[j][i].Y, cm[j][i].M);
        }
    }
  
}

int simulation(double space_size, long grid_size, long long num_particles, long long num_timesteps, particle_t *par){
    
    center_mass ** cm; //(line vs columns)
    int grid_x, grid_y; //position x and y of the particle that we want to calculate the forces exerted
    double distance_cell = 0; //distance to the particles in the same cell
    double distance_cm = 0; //distance to the center of mass of the other cells
    int position_x, position_y; //position of the other particles in the grid
    double delta_x = 0, delta_y = 0; //displacement of the particle in x and y
    double rx = 0, ry = 0; //direction along x and y
    int count = 0;
    int colisao_count = 0;
    double cell_size = (double)space_size / grid_size;
    ParColision colision[num_particles];

    // Allocation of memory to save the center of mass of each cell
    cm = malloc(grid_size * sizeof(center_mass*));
    for (int i=0; i < grid_size; i++){
        cm[i] = malloc(grid_size * sizeof(center_mass));
    }
    
    // Simulation loop
    for(int i = 0; i < num_timesteps; i++){
        calc_center_mass(cm, num_particles, par, space_size, grid_size); // Compute the center of mass at the current instant for every cell
        for(int j = 0; j < num_particles; j++){
            printf("Part: %d", j);
            par[j].Fx = 0;
            par[j].Fy = 0;

            double grid_x_aux = par[j].x / cell_size;
            grid_x = (int) grid_x_aux; 
        
            double grid_y_aux = par[j].y/ cell_size;
            grid_y = (int) grid_y_aux;
            
            for(int k = 0; k < num_particles; k++){
                if(j>=k){
                    continue;
                }

                double pos_x_aux = par[k].x / cell_size;
                position_x = (int) pos_x_aux; 
            
                double pos_y_aux = par[k].y/ cell_size;
                position_y = (int) pos_y_aux;

                if(position_x == grid_x && position_y == grid_y){   // Tenho quase a certeza que isto pode ser simplificado para:
                    distance_cell = sqrt(pow((par[j].x - par[k].x),2)+pow((par[j].y-par[k].y),2));  // Esta linha desaparece
                    
                    if(distance_cell <= 0.005){
                        if (colisao_count < num_particles) {
                            colision[colisao_count].a = j;
                            colision[colisao_count].b = k;
                            colisao_count++;
                        }
                        continue;
                    }
                    delta_x = position_x - grid_x;
                    delta_y = position_y - grid_y;
                    rx = delta_x/distance_cell;     // Esta linha desaparece
                    ry = delta_y/distance_cell;     // Esta linha desaparece
                    par[j].Fx += G*((par[k].m*par[j].m)/(pow(distance_cell,2)))*rx;     // par[j].Fx += G*((par[k].m*par[j].m)/(pow(delta_x,2)));
                    par[j].Fy += G*((par[k].m*par[j].m)/(pow(distance_cell,2)))*ry;     // par[j].Fx += G*((par[k].m*par[j].m)/(pow(delta_y,2)));
                }
            }

            delta_x = 0;
            delta_y = 0;
            rx = 0;
            ry = 0;
            

            for(int w =-1; w<2; w++){
                for(int y =-1; y<2; y++){

                    if(w==y && w== 0){
                        //printf("aqui");
                        continue;
                    }
                    if((grid_x + w)>=0 && (grid_y +y)>=0 && (grid_x+w)<grid_size && (grid_y+y)<grid_size){
                        distance_cm = sqrt(pow((par[j].x - cm[grid_x + w][grid_y + y].X),2)+pow((par[j].y-cm[grid_x + w][grid_y + y].Y),2));
                        delta_x = cm[grid_x + w][grid_y + y].X - grid_x;
                        delta_y = cm[grid_x  +w][grid_y + y].Y - grid_y;
                        rx = delta_x/distance_cm;
                        ry = delta_y/distance_cm;
                        par[j].Fx += G*((cm[w+ grid_x][y + grid_y].M*par[j].m)/(pow(distance_cm,2)))*rx;
                        par[j].Fy += G*((cm[w + grid_x][y+ grid_y].M*par[j].m)/(pow(distance_cm,2)))*ry;
                        printf("    %d %d ola\n", (grid_x + w), (grid_y +y));
                        continue;
                        
                    }
                    
                    printf("    %d %d adeus\n", (grid_x + w), (grid_y +y));
                }
            }
            
            par[j].ax = par[j].Fx/par[j].m; 
            par[j].ay = par[j].Fy/par[j].m;
            par[j].vx = par[j].vx + par[j].ax*DELTAT; //calculate new velocity along x
            par[j].vy = par[j].vy + par[j].ay*DELTAT; //calculate new velocity along y
            par[j].x = par[j].x + par[j].vx*DELTAT + 0.5*par[j].ax*pow(DELTAT,2); //calculate new position x
            par[j].y = par[j].y + par[j].vy*DELTAT + 0.5*par[j].ay*pow(DELTAT,2); //calculate new position y

        }      
    }

    for(int n=0; n<colisao_count; n++){
        for(int m=0; m<colisao_count; m++){
            if(n>=m){
                continue;
            }
            if(colision[n].b == colision[m].a){
                colisao_count--;
            }
        }
    }
    printf("x: %.3f ; y: %.3f \n",par[0].x, par[0].y);  
    
    printf("Colisions: %d \n", colisao_count);
    

    for(int i=0; i < grid_size; i++)
        free(cm[i]);
    free(cm);

    return colisao_count;

}

 /* void print_result(particle_t *par, int collisions){

    printf("x: %.4f ; y: %.4f \n",par[0].x, par[0].y);
    printf("Colisions: %d \n", collisions);
} */

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
    long long num_timesteps = atoll(argv[5]);
    
    printf("%lu", seed);

    particle_t* particles = malloc(num_particles * sizeof(particle_t));
    
    fprintf(stderr, "Generating initial state...\n");
    init_particles(seed, space_size, grid_size, num_particles, particles);

    printf("MAssa 0 -- %f", particles[0].m);
    
    fprintf(stderr, "Starting Simulation...\n");
    exec_time = -omp_get_wtime();
    
    simulation(space_size, grid_size, num_particles, num_timesteps, particles);
    
    exec_time += omp_get_wtime();
    fprintf(stderr, "Simulation ended.\n");

    fprintf(stderr, "%.1fs\n", exec_time);
    
    //print_result();
    
    free(particles);
}
