/*
 * Lattice-Boltzmann simulation definitions
 */
#ifndef __SIMULATION_H__
#define __SIMULATION_H__

/*
 * Simulation constants
 *
 * (some are declared using #define so they can be used to specify array sizes)
 */
							
// Number of time steps to take
extern unsigned int n_t_steps;

// Simulation configuration
#define CONFIG_PIPE            1
#define CONFIG_PERIODIC_BLOCK  2
#define CONFIG_USER1           3
#define CONFIG_USER2           4
#define CONFIG_MULTIPHASE	   5
#define CONFIG_CIRCLE          6
#define CONFIG_WETTINGCHANNEL  7 
extern unsigned int config;

// Size of flow domain
#define n_x_points   64
#define n_z_points   32

// Number of points for approximating the probability distribution
// (see comments before function propagation() for details)
#define n_value     9

/*
 *	Simulation parameters defined external (in main-cvi.c) as they 
 *	are made accessible to be changed via the user interface
 */

// Properties of the fluid
extern double global_density, viscosity;

// Body force in x- and z-direction
extern double f_x, f_z;

// Shan Chen force coefficient
extern double g_sc;

// Initial Perturbation
extern double init_perturb;

// Pseudo wall density used for wetting simulations
extern double obstacle_density;

// Flag to (de-) activate Shan Chen forces
extern unsigned int shan_chen_model;

/*
 * Function prototypes
 *
 * See simulation.c for comments.
 */

void   	simulation_init(void);
void   	simulation_step(void);
void   	simulation_clear_block(void);
void   	simulation_writeoutput(unsigned int field_counter);

void   	init_field(void);
void   	collision(void);
void   	boundary_walls_pipe(void);
void   	boundary_walls_periodic(void);
void   	clear_block(unsigned int bx1, unsigned int bx2, unsigned int bz1, unsigned int bz2);
void   	boundary_block(unsigned int bx1, unsigned int bx2, unsigned int bz1, unsigned int bz2);
int    	in_circle(double cx, double cz, double r, int x, int z);
void   	init_circle(double cx, double cz, double r);
void   	clear_obstacle_sites(void);
void   	boundary_surface_links(void);
void   	propagation(void);

double 	getdensity(unsigned int i, unsigned int k);
void   	calc_velocity(unsigned int i, unsigned int k, double *u_x, double *u_z);
void   	calc_shan_chen_force(unsigned int i, unsigned int k, double *f_sc_x,double *f_sc_z);

void   	write_velocities(unsigned int field_counter);
void   	write_pressures(unsigned int field_counter);
void   	write_dump(unsigned int field_counter);
void   	read_dump(const char *filename);

/* This function displays a message on screen and/or the log file. It is
 * defined in the main source file. */
void   message(const char *msg, ...);

#endif /* __SIMULATION_H__ */
