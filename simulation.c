/*
 * Lattice-Boltzmann simulation
 *
 * Based on fortran code by J.J Derksen
 * Adapted by B.J. Geurts for the course computational physics
 * Converted to C, added LabWindows/CVI GUI by L.M. Vergaij-Huizer and W. van Engen [02-2007]
 * Cleanup and comments by Z. Pouransari and W. van Engen [02-2008]
 * Extra functionalities added by A.W. van Cuijk [03-2010]
 * Extra functionalities and comments added by F. Janoschek and S.Schmieschek [02-2012]  
 * Comments added by S.Schmieschek [01-2013] 
 *
 */
#include "toolbox.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "simulation.h"

/*
 * Simulation data
 */

// Relaxation parameter 
double omega;

// Density distribution on grid points + boundaries
double mx[n_value][n_x_points+2][n_z_points+2];

// Structure to hold all the surface links of an arbitrarily shaped obstacle.
struct surface_link_t {
	int fluid_x,fluid_z;
	int wall_x,wall_z;
	int f2w;
};							
struct surface_link_t surface_links[(n_value-1)*n_x_points*n_z_points];

// Total number of bounce back rules aka surface links
int n_surface_links = 0;

// Inverse lattice directions
const int neg[n_value] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

// Structure to hold all the sites of an arbitrarily shaped obstacle.
struct obstacle_site_t {
	int x,z;
};
struct obstacle_site_t obstacle_sites[n_x_points*n_z_points];

// Total number of obstacle nodes defined in the surface link boundary scheme 
int n_obstacle_sites = 0;

// Array to hold weights of lattice directions for equilibrium calculation
double tp_i[n_value];	

// Array to hold weights of lattice directions for force calculation
double tf_i[n_value];

// Direction/probability matrix, x- and z-component.
int c_i[n_value][2];


/*
 * Setup data structures for simulation and print parameters
 */

void simulation_init()
{
	message("-----------------------------------------------------------------------\n");
	message("Initializing simulation\n");
	
	// Display simulation parameters
	message("          n_t_steps = %5d\n", n_t_steps);
	message("          viscosity = %10.8f\n", viscosity);								
	message("     global_density = %10.8f\n", global_density);
	message("        flow domain = %4d x %4d\n", n_x_points, n_z_points);
	message("            n_value = %4d\n", n_value);
	
	message("               g_sc = %10.8f\n", g_sc);
	message("       init_perturb = %10.8f\n", init_perturb);	
	// Calculate interal parameters
	omega = 2.0/(6*viscosity+1);
	//omega = 1.0;
	init_field();
//	clear_block(28, 36, 12, 20); 
}

/*
 * Execute a simulation step
 */
void simulation_step()
{
	unsigned int i, k;
	boundary_walls_periodic();
	collision();
	simulation_clear_block();
	
	switch(config) {
	case CONFIG_PIPE:
		boundary_walls_pipe();
		break;
		
	case CONFIG_PERIODIC_BLOCK:
		boundary_walls_periodic();
		boundary_block(20, 22,  12, 17);
		break;
		
	case CONFIG_MULTIPHASE:
		boundary_walls_periodic();
		break;

	case CONFIG_CIRCLE:
		boundary_walls_periodic();
		boundary_surface_links();
		boundary_walls_periodic();
		break;		
		
	case CONFIG_WETTINGCHANNEL:
		boundary_walls_periodic();
		boundary_block(1, n_x_points, 1, 1);
		boundary_block(1, n_x_points, n_z_points, n_z_points); 
		boundary_walls_periodic();
		break;	

	case CONFIG_USER1:
		// Put your own configuration here!
		break;
		
	case CONFIG_USER2:
		// And another configuration here.
		break;

	default:
		message("Error: unrecognized configuration!!!\n");
		break;
	}
	
	propagation();

	simulation_clear_block();
}

/*
 * Clear the area of the block
 *
 * This function can be called at any point during the simulation
 * without affecting the result.
 * This is used in the simulation and just before visualisation.
 */
void simulation_clear_block()
{
	unsigned int i, k;
	
	switch(config) {
	case CONFIG_PIPE:
		// The pipe configuration has no block present.
		break;
		
	case CONFIG_PERIODIC_BLOCK:
		clear_block(20, 22,  12, 17);
		break;
		
	case CONFIG_MULTIPHASE:
		break;

	case CONFIG_CIRCLE:
		clear_obstacle_sites();
		break;

	case CONFIG_WETTINGCHANNEL:
		clear_block(1, n_x_points, 1, 1);
		clear_block(1, n_x_points, n_z_points, n_z_points);
	break;

	case CONFIG_USER1:
		// Clear your own block here!
		break;
		
	case CONFIG_USER2:
		// Clear another block configuration here
		break;
	}
}

/*
 * Write simulation results to data files.
 *
 * The field_counter is a file number, so multiple files can be written
 * during one program run. The output files are named with the three-digit
 * number indicated by field_counter.
 * The program as supplied only uses a field_counter of zero.
 */
void simulation_writeoutput(unsigned int field_counter)
{

	write_velocities(field_counter);
	write_pressures(field_counter);
	write_dump(field_counter);

}

/*
 * Initialize the field
 * - Calculate the values of c_ix, c_iy, c_iz, tp_i, tf_i
 * - Calcluate the initial density distribution
 */
void init_field()
{
	unsigned int i, j, k, q;

	// Initialize c_ix, c_iz, tp_i and tf_i
	tp_i[0] = 4.0/9.0;
	tf_i[0] = 0.0;
	for (k=1; k<5; k++) {
		tp_i[k] = 1.0/9.0;
		tf_i[k] = 1.0/3.0;
	}
	for (k=5; k<n_value; k++) {
		tp_i[k] = 1.0/36.0;
		tf_i[k] = 1.0/12.0;
	}

	// Direction/probability matrix for q/n_value in x-direction
	c_i[8][0] =  1;  c_i[4][0] =  0;  c_i[7][0] = -1;
	c_i[1][0] =  1;  c_i[0][0] =  0;  c_i[3][0] = -1;
	c_i[5][0] =  1;  c_i[2][0] =  0;  c_i[6][0] = -1;
	// Direction/probability matrix for q/n_value in z-direction
	c_i[8][1] = -1;  c_i[4][1] = -1;  c_i[7][1] = -1;
	c_i[1][1] =  0;  c_i[0][1] =  0;  c_i[3][1] =  0;
	c_i[5][1] =  1;  c_i[2][1] =  1;  c_i[6][1] =  1;

	// Calculate initial density distribution
	for (i=0; i<n_x_points+2; i++) {
		for (k=0; k<n_z_points+2; k++) {
			for (q=0; q<n_value; q++) {
				mx[q][i][k] = tp_i[q] * global_density + (double) Random(-1,1) * init_perturb;
			}
		}
	}

	switch(config) {
	case CONFIG_PIPE:
		break;
	case CONFIG_PERIODIC_BLOCK:
		break;
	case CONFIG_CIRCLE:
		n_surface_links = 0;
		n_obstacle_sites = 0;
		init_circle(15.5, 15.5, 5.0);
		break;
	case CONFIG_USER1:
		break;
	case CONFIG_USER2:
		break;
	}
	
}

/*
 * Calculate the boundary conditions at the walls
 * for a completely periodic configuration.
 */
void boundary_walls_periodic()
{
	unsigned int i, k, q;

	for (k=0; k<=n_z_points+1; k++) {
		for (q=0; q<n_value; q++) {
			mx[q][0][k] = mx[q][n_x_points][k];
			mx[q][n_x_points+1][k] = mx[q][1][k];
		}
	}
	for (i=0; i<=n_x_points+1; i++) {
		for (q=0; q<n_value; q++) {
			mx[q][i][0] = mx[q][i][n_z_points];
			mx[q][i][n_z_points+1] = mx[q][i][1];
		}
	}
}

/*
 * Calculate the boundary conditions at the walls
 * for a pipe configuration.
 */
void boundary_walls_pipe()
{
	unsigned int i, k, q;

	for (k=0; k<=n_z_points+1; k++) {
		for (q=0; q<n_value; q++) {
			mx[q][0][k] = mx[q][n_x_points][k];
			mx[q][n_x_points+1][k] = mx[q][1][k];
		}
	}
	for (i=1; i<=n_x_points; i++) {
		// No slip at the bottom of the channel
		mx[2][i  ][0] = mx[4][i][1];
		mx[5][i-1][0] = mx[7][i][1];
		mx[6][i+1][0] = mx[8][i][1];
		// And neither at the top of the channel		
		mx[4][i  ][n_z_points+1] = mx[2][i][n_z_points];
		mx[7][i+1][n_z_points+1] = mx[5][i][n_z_points];
		mx[8][i-1][n_z_points+1] = mx[6][i][n_z_points];
	}
}

/*
 * Calculate the boundary conditions at the specified block
 */
void boundary_block(unsigned int bx1, unsigned int bx2,
	unsigned int bz1, unsigned int bz2)
{
	unsigned int i, k;

	for (i=bx1; i<=bx2; i++) {
		mx[2][i][bz2] = mx[4][i  ][bz2+1];
		mx[5][i][bz2] = mx[7][i+1][bz2+1];
		mx[6][i][bz2] = mx[8][i-1][bz2+1];

		mx[4][i][bz1] = mx[2][i  ][bz1-1];
		mx[7][i][bz1] = mx[5][i-1][bz1-1];
		mx[8][i][bz1] = mx[6][i+1][bz1-1];
	}

	for (k=bz1; k<=bz2; k++) {
		mx[1][bx2][k] = mx[3][bx2+1][k];
		mx[5][bx2][k] = mx[7][bx2+1][k+1];
		mx[8][bx2][k] = mx[6][bx2+1][k-1];

		mx[3][bx1][k] = mx[1][bx1-1][k];
		mx[6][bx1][k] = mx[8][bx1-1][k+1];
		mx[7][bx1][k] = mx[5][bx1-1][k-1];
	}
}

/*
 * Set the density distribution of a block to zero.
 */
void clear_block(unsigned int bx1, unsigned int bx2,
	unsigned int bz1, unsigned int bz2)
{
	unsigned int i, k, q;
	
	for (i=bx1; i<=bx2; i++) {
		for (k=bz1; k<=bz2; k++) {
			if (obstacle_density == 0.0) { 
				mx[0][i][k] = global_density;
			} else {
				mx[0][i][k] = obstacle_density;
			}
			for (q=1; q<n_value; q++) {
				mx[q][i][k] = 0.0;
			}
		}
	}
}

/*
 * Returns 1 if (x,z) is within a circle around (cx,cz) with radius r or 0 otherwise.
 * Periodic boundaries are taken account of.
 */
int in_circle(double cx, double cz, double r, int x, int z)
{
	// minimum image criterion
	double dx = cx - (double) x;
	double dz = cz - (double) z;	
	
	if (dx >= 0.5 * (double) n_x_points) dx -= (double) n_x_points;
	else if (dx < -0.5 * (double) n_x_points) dx += (double) n_x_points;

	if (dz >= 0.5 * (double) n_z_points) dz -= (double) n_z_points;
	else if (dz < -0.5 * (double) n_z_points) dz += (double) n_z_points;
	
	return (dx * dx + dz * dz <= r * r) ? 1 : 0;
}

/*
 * Sets up obstacle_sites and surface_links for a circular obstacle with
 * radius r at position (cx,cz).
 */
void init_circle(double cx, double cz, double r)
{
	int x, z, q;
	for (x = 1; x <= n_x_points; ++x) {
		for (z = 1; z <= n_z_points; ++z) {
			if (in_circle(cx, cz, r, x, z)) {
				struct obstacle_site_t * o = &(obstacle_sites[n_obstacle_sites]);
				o->x = x;
				o->z = z;
				++n_obstacle_sites;
			} else {
				for (q = 1; q < n_value; ++q) {
					int tx = x + c_i[q][0];
					int tz = z + c_i[q][1];					
					if (tx == 0) tx = n_x_points;
					else if (tx == n_x_points + 1) tx = 1;
					if (tz == 0) tz = n_z_points;
					else if (tz == n_z_points + 1) tz = 1;
					if (in_circle(cx, cz, r, tx, tz)) {
						struct surface_link_t * l = &(surface_links[n_surface_links]);
						l->fluid_x = x;
						l->fluid_z = z;
						l->wall_x = tx;
						l->wall_z = tz;
						l->f2w = q;
						++n_surface_links;
					}
				}
			}
		}
	}
}

/*
 * Implements the bounce back boundary condition for all arbitrarily shaped
 * obstacles.
 */
void boundary_surface_links(void)
{
	int i;
	for (i = 0; i < n_surface_links; ++i) {
		struct surface_link_t * l = &(surface_links[i]);
		mx[neg[l->f2w]][l->wall_x][l->wall_z] = mx[l->f2w][l->fluid_x][l->fluid_z];
	}
}

/*
 * Removes the velocity from the interior of all aritrarily shaped obstacles.
 */
void clear_obstacle_sites(void)
{
	int i, q;
	for (i = 0; i < n_obstacle_sites; ++i) {
		struct obstacle_site_t * o = &(obstacle_sites[i]);
		mx[0][o->x][o->z] = global_density;
		for (q = 1; q < n_value; ++q) mx[q][o->x][o->z] = 0.0;
	}
}

/*
 * Propagate probability to the next cells, depending on direction.
 * n_value specifies the number of directions and their meanings are
 * as follows (with the default n_value of 9):
 *
 *       6    2    5
 *            ^
 *          \ | /
 *       3 <--0--> 1
 *   z      / | \
 *   ^        v  
 *   |   7    4    8
 *   |
 *   +---> x
 *
 * The direction corresponding to the numbers are specified in c_i[],
 * which is defined in init_field().
 * Note that this function and init_field() need to be modified when
 * this number is changed.
 */
void propagation()
{
	unsigned int i, j, k;

	// loop backward in z-direction for positive z-direction
	for (k=n_z_points+1; k>0; k--) {
		// loop backward in x-direction too for positive x-direction
		for (i=n_x_points+1; i>0; i--) {
			mx[1][i][k] = mx[1][i-1][k];
			mx[5][i][k] = mx[5][i-1][k-1];
			mx[2][i][k] = mx[2][i  ][k-1];
		}
		// loop forward in x-direction for negative x-direction
		for (i=0; i<n_x_points+1; i++) {
			mx[6][i][k] = mx[6][i+1][k-1];
		}
	}

	// loop forward in z-direction for negative z-direction
	for (k=0; k<n_z_points+1; k++) {
		// loop forward in x-direction too for negative x-direction
		for (i=0; i<n_x_points+1; i++) {
			mx[3][i][k] = mx[3][i+1][k];
			mx[4][i][k] = mx[4][i  ][k+1];
			mx[7][i][k] = mx[7][i+1][k+1];
		}
		// loop backward in x-direction for positive x-direction
		for (i=n_x_points+1; i>0; i--) {
			mx[8][i][k] = mx[8][i-1][k+1];
		}
	}
}

/*
 * The collision step (LBGK 9-speed)
 */
void collision()
{
	unsigned int i, j, k, q;
	double f_c[n_value];
	double f_sc_x,f_sc_z;
	double density, u_x, u_z, c_u, u_2, ni_eq;

	for (q=0; q<n_value; q++) {
		f_c[q] = f_x * c_i[q][0] + f_z * c_i[q][1];
		f_c[q] *= tf_i[q];
	}

	for (k=1; k<n_z_points+1; k++) {
		for (i=1; i<n_x_points+1; i++) {
			density = getdensity(i, k);
			calc_velocity(i, k, &u_x, &u_z);

			if (shan_chen_model == 1) {
				// Calculate the Shan Chen force
				calc_shan_chen_force(i, k, &f_sc_x, &f_sc_z);
				// Add a Shan Chen force contribution to the equilibrium velocity
				u_x += (1.0 / omega) * f_sc_x / density;
				u_z += (1.0 / omega) * f_sc_z / density;
			}
			
			for (q=0; q<n_value; q++) {
				// Calculate equilibrium distribution
				c_u = c_i[q][0] * u_x + c_i[q][1] * u_z;
				u_2 = u_x*u_x + u_z*u_z;
				// ----- make changes here -----
				ni_eq = 1.0 + 3.0 * c_u + 4.5 * c_u*c_u - 1.5 * u_2;
				// -----------------------------
				ni_eq *= tp_i[q] * density;
				mx[q][i][k] += -omega * ( mx[q][i][k] - ni_eq ) + f_c[q];
			}
		}
	}
}


/*
 * Get the density at a gridpoint
 *
 * This is the sum of the densities in all directions.
 */
double getdensity(unsigned int i, unsigned int k)
{
	unsigned int q;
	double density=0.0;

	for (q=0; q<n_value; q++)
		density += mx[q][i][k];

	return density;
}

/*
 * Return the velocity at a gridpoint
 */
void calc_velocity(unsigned int i, unsigned int k, double *u_x, double *u_z)
{
	unsigned int q;
	double loc_density = 0.0;
	*u_x = 0.0;
	*u_z = 0.0;
	
	for (q=0; q<n_value; q++) {
		*u_x += mx[q][i][k] * c_i[q][0];
		*u_z += mx[q][i][k] * c_i[q][1];
		loc_density += mx[q][i][k];
	}
	*u_x /= loc_density;
	*u_z /= loc_density;
}

/*
 *	Calculated the Shan Chen force at a grid point
 */

void calc_shan_chen_force(unsigned int i, unsigned int k, double *f_sc_x, double *f_sc_z)
{
	unsigned int 	q;
	unsigned int 	ii,kk;
	double 			f_sc[n_value];
	double 			rho_loc = 0.0;
	double 			rho_neighbour = 0.0;
	double 			psi_loc = 0.0;
	double 			psi_neighbour = 0.0;

	// Initialize referenced variables
	*f_sc_x = 0.0;
	*f_sc_z = 0.0;
	
	// Calculate the local density functional
	rho_loc = getdensity(i,k);
	psi_loc = (1-exp(-rho_loc));
	
	// Loop over all neighbouring lattice sites
	for (q=0; q<n_value; q++) {
		// Get coordinate of neighbour in direction q
		ii = i + c_i[q][0];
		kk = k + c_i[q][1];

		// Calculate the neighbour density functional
		rho_neighbour = getdensity(ii,kk);
		psi_neighbour = (1-exp(-rho_neighbour));
		
		// Calculate the Shan Chen force
		f_sc[q] = -g_sc * psi_loc * tp_i[q] * psi_neighbour;
	}
	
	// Project the force onto Cartesian coordinates
	for (q=0; q<n_value; q++) {
		*f_sc_x += f_sc[q] * c_i[q][0];
		*f_sc_z += f_sc[q] * c_i[q][1];
	}
}

/*
 * Write the velocities to a text file.
 *
 * The field_counter is a file number, so multiple velocity files can be
 * written during one program run. The output files are velocityXXX.txt, with
 * XXX the number indicated by field_counter.
 * The program as supplied only uses a field_counter of zero.
 */
void write_velocities(unsigned int field_counter)
//void write_velocities(char* suffix)
{
	double u_x, u_z;
	unsigned int i, k;
	char filename[255];
	FILE *f;

//	sprintf(filename, "velocity_%s.txt", suffix);
	sprintf(filename, "velocity_%03d.txt", field_counter); 
	f = fopen(filename, "w");
	if (!f) {
		message("Error: could not open velocity output file '%s'\n", filename);
		return;
	}
	// write header, start with '%' to make matlab happy
	fprintf(f, "%%  x               z               u_x             u_z\n");
	
	// ... and data
	for (i=1; i<n_x_points+1; i++) {
		for (k=1; k<n_z_points+1; k++) {
			calc_velocity(i, k, &u_x, &u_z);
			fprintf(f, "%15e %15e %15e %15e\n", i*1.0-0.5, k*1.0-0.5, u_x, u_z);
		
		}
	}
	
	
	fclose(f);
}


/*
 * Write the pressures (=densities) to a text file.
 *
 * The field_counter is a file number, so multiple pressure files can be
 * written during one program run. The output files are pressureXXX.txt, with
 * XXX the number indicated by field_counter.
 * The program as supplied only uses a field_counter of zero.
 */ 
//void write_pressures(char* suffix)
void write_pressures(unsigned int field_counter)
{
	unsigned int i, k;
	char filename[255];
	FILE *f;

	//sprintf(filename, "pressure_%s.txt", suffix);
	sprintf(filename, "density_%03d.txt", field_counter);     
	f = fopen(filename, "w");
	if (!f) {
		message("Error: could not open pressure output file '%s'\n", filename);
		return;
	}
	
	// write header ...
	fprintf(f, "%%  x               z               p\n");
	
	// ... and data
	for (i=1; i<n_x_points+1; i++) {
		for (k=1; k<n_z_points+1; k++) {
			fprintf(f, "%15e %15e %15e\n", i*1.0-0.5, k*1.0-0.5, getdensity(i,k));
		}
	}
	
	fclose(f);
}

/*
 * Write simulation data to a binary data file.
 *
 * The field_counter is a file number, so multiple dump files can be
 * written during one program run. The output files are dumpXXX.dat, with
 * XXX the number indicated by field_counter.
 * The program as supplied only uses a field_counter of zero.
 */
//void write_dump(char* suffix)
void write_dump(unsigned int field_counter)
{
	unsigned int nx = n_x_points, nz = n_z_points;
	char filename[255];
	FILE *f;

//	sprintf(filename, "dump_%s.dat", suffix);
	sprintf(filename, "dump_%03d.dat", field_counter);
	f = fopen(filename, "wb");
	if (!f) {
		message("Error: could not open dump output file '%s'\n", filename);
		return;
	}

	if ( fwrite((void*)&nx, sizeof(unsigned int), 1, f) != 1 ||
		 fwrite((void*)&nz, sizeof(unsigned int), 1, f) != 1 ||
		 fwrite((void*)mx, sizeof(mx), 1, f) != 1 ||
		 fwrite((void*)tp_i, sizeof(tp_i), 1, f)!=1 ||
		 fwrite((void*)tf_i, sizeof(tf_i), 1, f)!=1
		 	) {
		message("Error: dump file '%s' write error\n", filename);
		return;
	}
	
	message("Simulation data written to file '%s'\n", filename);

	fclose(f);
}

/*
 * Read simulation data from a binary data file.
 *
 * All relevant data is read from the binary file, so a previous simulation
 * can be processed or continued. Note that the program has to be compiled
 * with the same number of x and z gridpoints.
 */
void read_dump(const char *filename)
{
	FILE *f;
	char c;
	unsigned int nx, nz;

	f = fopen(filename, "rb");
	if (!f) {
		message("Error: could not open dump input file '%s'\n", filename);
		exit(1);
	}

	// Read the file and make sure filesize is correct
	if (fread((void*)&nx, sizeof(unsigned int), 1, f) != 1 ||
		fread((void*)&nz, sizeof(unsigned int), 1, f) != 1 ) {
		message("Error: dump file '%s' header read error\n", filename);
		return;
	}
	if (nx != n_x_points || nz != n_z_points) {
		message("Error: grid mismatch; recompile with n_x_points=%u and n_z_points=%u\n", nx, nz);
		return;
	}
	
	if (fread((void*)mx, sizeof(mx), 1, f) != 1 ||
		fread((void*)tp_i, sizeof(tp_i), 1, f) != 1 ||
		fread((void*)tf_i, sizeof(tf_i), 1, f) != 1
			) {
		message("Error: dump file '%s' read error\n", filename);
		return;
	}
	
	if (fread((void*)&c, 1, 1, f)!=0) {
		message("Error: dump file '%s' is too large\n", filename);
		return;
	}
	
	message("Simulation data read from file '%s'\n", filename);
	
	fclose(f);
}

