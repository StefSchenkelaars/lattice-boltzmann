/*
 * Lattice-Boltzmann LabWindows/CVI GUI program
 *
 * GUI front-end for Lattice-Boltzmann simulation
 * by L.M. Vergaij-Huizer and W. van Engen
 */
#include <ansi_c.h>
#include <toolbox.h>
#include <userint.h>
#include <stdarg.h>
#include <float.h>

#include "simulation.h"											
#include "LB1_UI.h"


// Local function prototypes
void gui_simulation_run(void);
void gui_simulation_display(void);
void gui_set_graph_attributes(void);

#define TYPE_VELOCITY	1	// calculate and plot velocity graphs
#define TYPE_DENSITY	2	// calculate and plot density graph
#define TYPE_VELX		3
#define TYPE_VELZ		4

void plotdata(unsigned int type);

// File global variables
static int hPanel;
FILE *logfile = NULL;

// Definition of starting parameters
unsigned int n_t_steps = 10000;
unsigned int config = CONFIG_PIPE;
int           dovisual = 1;
int           dovectors = 0;
unsigned int  n_drawing_steps = 50;
unsigned char init_needed = 0;
unsigned char running = 0;
double 		  global_density = 1.0;
double 		  viscosity = 1.6667E-1;
double		  f_x = 1E-4;
double		  f_z = 0.0;
unsigned int  shan_chen_model = 0;
double		  g_sc = -4.1;
double		  obstacle_density = 0.0;    
double		  init_perturb = 0.0;

/*
 * Program entry point: execution starts here
 */
int main(int argc, char *argv[])
{
	int hTab;
	
	if ((hPanel = LoadPanel (0, "LB1_UI.uir", PANEL)) < 0)
        return -1;
	
	logfile = fopen("log.txt", "w");
	if (!logfile)
		message("Warning: could not open 'log.txt' for writing\n");

	// Set controls according to values
	SetCtrlVal(hPanel, PANEL_CONFIG, config);
	SetCtrlVal(hPanel, PANEL_NUM_DRAWPERIOD, n_drawing_steps);
	SetCtrlVal(hPanel, PANEL_TOG_DRAW, dovisual);
	SetCtrlVal(hPanel, PANEL_TOG_VEC, dovectors);
	SetCtrlVal(hPanel, PANEL_NUM_TIMESTEPS, n_t_steps);
	SetCtrlVal(hPanel, PANEL_TOG_SC, shan_chen_model);

	gui_set_graph_attributes();
	
	DisplayPanel (hPanel);
	
	simulation_init();

	RunUserInterface();
	
	if (logfile) fclose(logfile);
}


// We use one callback function for all buttons and use the control argument
// to distinguish between different buttons.
int CVICALLBACK button_pressed(int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	char textbuffer[MAX_PATHNAME_LEN];
	
	if (event!= EVENT_COMMIT)
		return 0;
	
	switch(control) {
	case PANEL_BTN_INIT:
		simulation_init();
		init_needed = 0;
		break;
		
	case PANEL_BTN_RUN:
		if (!running)
			gui_simulation_run();
		else
			running = 0;
		break;

	case PANEL_BTN_DISPLAY_FILE:
		if(FileSelectPopup ("","*.dat","*.dat","Select data dump input file",VAL_SELECT_BUTTON,0,1,1,0,textbuffer)) {
			read_dump(textbuffer);
			gui_simulation_display();
		}

		break;
		
	case PANEL_BTN_DISPLAY:
		gui_simulation_display();
		break;
		
	case PANEL_BTN_RESET:
		gui_set_graph_attributes();
		gui_simulation_display();

	default:
		break;
	}
	
	return 0;
}

// We use one callback function for all value controls and use the control
// argument to distinguish between different controls.
int CVICALLBACK control_changed(int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	if (event!=EVENT_COMMIT)
		return 0;

	switch (control) {
	case PANEL_CONFIG:
		GetCtrlVal(panel, control, &config);
		init_needed = 1;
		break;
		
	case PANEL_NUM_TIMESTEPS:
		GetCtrlVal(panel, control, &n_t_steps);
		break;
		
	case PANEL_TOG_DRAW:
		GetCtrlVal(panel, control, &dovisual);
		break;
		
	case PANEL_TOG_VEC:
		GetCtrlVal(panel, control, &dovectors);
		break;
		
	case PANEL_NUM_DRAWPERIOD:
		GetCtrlVal(panel, control, &n_drawing_steps);
		break;
		
	case PANEL_GLOBAL_DENSITY:
		GetCtrlVal(panel, control, &global_density);
		break;
		
	case PANEL_VISCOSITY:
		GetCtrlVal(panel, control, &viscosity);
		break;
		
	case PANEL_F_X:
		GetCtrlVal(panel, control, &f_x);
		break;
		
	case PANEL_F_Z:
		GetCtrlVal(panel, control, &f_z);
		break;
		
	case PANEL_TOG_SC:
		GetCtrlVal(panel, control, &shan_chen_model);
		break;
		
	case PANEL_G_SC:
		GetCtrlVal(panel, control, &g_sc);
		break;
		
	case PANEL_INIT_PERTURBATION:
		GetCtrlVal(panel, control, &init_perturb);
		break;
		
	case PANEL_WALL_DENS:
		GetCtrlVal(panel, control, &obstacle_density);
		break;
	
	}
	return 0;
}

// We use one callback function for all controls related to the cursor and
// use the control argument to distinguish between different controls.
int CVICALLBACK cursor_changed (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	double cursor_x,cursor_z;
	int cursor_x_int,cursor_z_int;
	int hTab;

	// Cursor on 2D velocity plot was changed
	if (event==EVENT_VAL_CHANGED && control==TABTOP0_GRAPH_VELO_COL) {
		GetGraphCursor(panel,TABTOP0_GRAPH_VELO_COL,1,&cursor_x,&cursor_z);
		SetCtrlVal(hPanel,PANEL_CURSOR_X,(int)(cursor_x+0.5));
		SetCtrlVal(hPanel,PANEL_CURSOR_Z,(int)(cursor_z+0.5));
		gui_simulation_display();
		
	// Numeric controls change
	} else if (event==EVENT_COMMIT) {
		switch (control) {
		case PANEL_CURSOR_X:
		case PANEL_CURSOR_Z:
			GetCtrlVal(panel,PANEL_CURSOR_X,&cursor_x_int);
			GetCtrlVal(panel,PANEL_CURSOR_Z,&cursor_z_int);
			cursor_x=cursor_x_int-0.5;
			cursor_z=cursor_z_int-0.5;
			GetPanelHandleFromTabPage (hPanel, PANEL_TAB_TOP, 0, &hTab);
			SetGraphCursor(hTab,TABTOP0_GRAPH_VELO_COL,1,cursor_x,cursor_z);
			gui_simulation_display();
			break;
		}
	}
	
	return 0;
}

int CVICALLBACK end (int panel, int event, void *callbackData,
		int eventData1, int eventData2)
{
	switch (event)
	{
	case EVENT_CLOSE:
		QuitUserInterface(0);
		break;
	}
	return 0;
}

// Set active axes and scaling modes and other attributes
void gui_set_graph_attributes()
{
	int hTab;
	
	GetPanelHandleFromTabPage (hPanel, PANEL_TAB_TOP, 0, &hTab);
	SetCtrlAttribute(hTab,TABTOP0_LEGEND_VELO,ATTR_ACTIVE_YAXIS,VAL_RIGHT_YAXIS);
	SetCtrlAttribute(hTab, TABTOP0_LEGEND_VELO, ATTR_REFRESH_GRAPH, 0);
	SetAxisScalingMode(hTab,TABTOP0_GRAPH_VELO_COL,VAL_BOTTOM_XAXIS,VAL_MANUAL,0,n_x_points);
	SetAxisScalingMode(hTab,TABTOP0_GRAPH_VELO_COL,VAL_LEFT_YAXIS,VAL_MANUAL,0,n_z_points);
	SetCtrlAttribute(hTab, TABTOP0_GRAPH_VELO_COL, ATTR_REFRESH_GRAPH, 0);
	
	GetPanelHandleFromTabPage (hPanel, PANEL_TAB_BOTTOM, 0, &hTab);
	SetCtrlAttribute(hTab,TABBOTTOM0_LEGEND_DENS,ATTR_ACTIVE_YAXIS,VAL_RIGHT_YAXIS);
	SetCtrlAttribute(hTab, TABBOTTOM0_LEGEND_DENS, ATTR_REFRESH_GRAPH, 0);
	SetAxisScalingMode(hTab,TABBOTTOM0_GRAPH_DENS_COL,VAL_BOTTOM_XAXIS,VAL_MANUAL,0,n_x_points);
	SetAxisScalingMode(hTab,TABBOTTOM0_GRAPH_DENS_COL,VAL_LEFT_YAXIS,VAL_MANUAL,0,n_z_points);
	SetCtrlAttribute(hTab, TABBOTTOM0_GRAPH_DENS_COL, ATTR_REFRESH_GRAPH, 0);

	
	GetPanelHandleFromTabPage (hPanel, PANEL_TAB_BOTTOM, 1, &hTab);
	SetCtrlAttribute(hTab, TABBOTTOM1_GRAPH_VELO_PROF_X, ATTR_REFRESH_GRAPH, 0);
	SetAxisScalingMode(hTab,TABBOTTOM1_GRAPH_VELO_PROF_X,VAL_BOTTOM_XAXIS,VAL_MANUAL,0,n_x_points);
	SetCtrlAttribute(hPanel,PANEL_CURSOR_X,ATTR_MAX_VALUE,n_x_points);

	GetPanelHandleFromTabPage (hPanel, PANEL_TAB_BOTTOM, 2, &hTab);
	SetCtrlAttribute(hPanel,PANEL_CURSOR_Z,ATTR_MAX_VALUE,n_z_points);
	SetAxisScalingMode(hTab,TABBOTTOM2_GRAPH_VELO_PROF_Z,VAL_BOTTOM_XAXIS,VAL_MANUAL,0,n_z_points);
	SetCtrlAttribute(hTab, TABBOTTOM2_GRAPH_VELO_PROF_Z, ATTR_REFRESH_GRAPH, 0);
}

// Create a colormap for the intensityplot
void gui_colormap_create(ColorMapEntry colormap[64], double min_data, double max_data)
{
	unsigned int i;
	
	for (i=0;i<64;i++)
	{
		colormap[i].dataValue.valDouble=min_data+i*(max_data-min_data)/64.0;
		
		if (i<=8)
			colormap[i].color=MakeColor(0,0,255*(0.5+i/16.0));
		else if (i<=24)
			colormap[i].color=MakeColor(0,(i-8)*255/16,255);
		else if (i<=40)
			colormap[i].color=MakeColor((i-24)*255/16,255,(41-i)*255/16);
		else if (i<=57)
			colormap[i].color=MakeColor(255,(57-i)*255/16,0);
		else 
			colormap[i].color=MakeColor((73-i)*255/16,0,0);
	}			 

	colormap[0].color=MakeColor(0,0,10);
}

// Update the legend range for an intensity plot
void gui_intensity_legend_update(int hPanel, unsigned int control, double min_data, double max_data, ColorMapEntry colormap[64])
{
	static double range[65][2];
	unsigned int i, k;

	// fill the range array
	for (k=0; k<65; k++)
		for (i=0; i<2; i++)
			range[k][i]=min_data+k*(max_data-min_data)/64.0;

	if (min_data==max_data) max_data=min_data+1e-8; // to avoid error in SetAxisScalingMode
	if (max_data>=min_data)
		SetAxisScalingMode(hPanel,control,VAL_RIGHT_YAXIS,VAL_MANUAL,min_data,max_data);
	else
		SetAxisScalingMode(hPanel,control,VAL_RIGHT_YAXIS,VAL_MANUAL,min_data-DBL_MIN,max_data+DBL_MIN);
	PlotScaledIntensity(hPanel,control,range,2,65,VAL_DOUBLE,(max_data-min_data)/64.0,min_data,1,0,colormap,MakeColor(10,0,0),64,0,0);
}

// Put vectors on a plot
void gui_plot_vectors(int hPanel, unsigned int idCtrl, 
	double mag[n_z_points][n_x_points], double x[n_z_points][n_x_points], double z[n_z_points][n_x_points])
{
	unsigned int i, k;
	double u_x, u_z, u_mag;
	double x1, x2, z1, z2, alpha;
	
	const double arrow_length = 0.6;
	const double  head_length = 0.3;
	const double  head_angle  = PI/6.0;

	for (i=0; i<n_x_points; i++) {
		for (k=0; k<n_z_points; k++) {
			u_mag = mag[k][i];
			u_x =   x[k][i];
			u_z =   z[k][i];
			
			// only plot points if velocity isn't close to zero
			if ( u_mag == 0.0 ) // warning: floating pt comparisons are a bad idea
				continue;
			
			// begin position of vector
			x1 = i + 0.5;
			z1 = k + 0.5;
			// end position of vector
			x2 = x1 + u_x/u_mag*arrow_length;
			z2 = z1 + u_z/u_mag*arrow_length;
			// angle for arrow head
			if (z2==z1)
				alpha = (x2>x1) ? PI/2.0 : -PI/2.0;
			else
				alpha = atan( (x2-x1)/(z2-z1) );
			if (z2<z1) alpha = alpha+PI;
			
			PlotLine (hPanel, idCtrl, x1, z1, x2, z2, VAL_WHITE);
			PlotLine (hPanel, idCtrl, x2, z2, x2-head_length*sin(alpha+head_angle), z2-head_length*cos(alpha+head_angle), VAL_WHITE);
			PlotLine (hPanel, idCtrl, x2, z2, x2-head_length*sin(alpha-head_angle), z2-head_length*cos(alpha-head_angle), VAL_WHITE);
			
			//PlotPoint(hPanel, idCtrl, i+0.5, k+0.5, VAL_SOLID_DIAMOND, VAL_RED);
		}
	}
	
}

// Plot an intensity plot with legend
void gui_plot_intensity(int hPanel, int idGraph, int idLegend, double mag[n_z_points][n_x_points])
{
	double min_data=DBL_MAX, max_data=0.0;
	int i, k;
	ColorMapEntry colormap[256];
	
	for (i=0; i<n_x_points; i++) {
		for (k=0; k<n_z_points; k++) {
			if ( mag[k][i] > max_data ) max_data = mag[k][i];
			if ( mag[k][i] < min_data ) min_data = mag[k][i];
		}
	}
	
	gui_colormap_create(colormap, min_data-(min_data/10.0), max_data+(max_data/10.0));
	PlotScaledIntensity(hPanel, idGraph, mag, n_x_points, n_z_points, VAL_DOUBLE, 1,0.5,1,0.5,colormap,MakeColor(10,0,0),64,0,0);
	gui_intensity_legend_update(hPanel, idLegend, min_data, max_data, colormap);
}

// Plot X- and Z-profiles. This is done in one function to avoid dynamic
// memory allocation.
void gui_plot_profiles(int hPanelX, int idGraphX, int hPanelZ, int idGraphZ,
	double mag[n_z_points][n_x_points], double x[n_z_points][n_x_points], double z[n_z_points][n_x_points],
	int cursor_x, int cursor_z)
{
	double profx_val[n_x_points]; // x-axis in x profile
	double profx_x[n_x_points], profx_z[n_x_points], profx_mag[n_x_points]; // y-axes in x profile
	double profz_val[n_z_points]; // x-axis in z profile
	double profz_x[n_z_points], profz_z[n_z_points], profz_mag[n_z_points]; // y-axes in z profile
	double parabola[n_z_points];
	double parabola_val[n_z_points]; 
	double delta[n_z_points];
	double max_z;
	int line_width = VAL_FAT_LINE;
	int plotid;
	unsigned int i,k;
	
	max_z = 0.0;
	// extract profile data
	for (i=0; i<n_x_points; i++) {
		profx_val[i] = i+0.5; // center data point horizontally
		
		profx_x[i]   =   x[cursor_z][i];
		profx_z[i]   =   z[cursor_z][i];
		profx_mag[i] = mag[cursor_z][i];
		
	}
	for (k=0; k<n_z_points; k++) {
		parabola_val[k] = k+1.5;
		profz_val[k] = k+0.5; // center data point horizontally
		profz_x[k]   =   x[k][cursor_x];
		profz_z[k]   =   z[k][cursor_x];
		profz_mag[k] = mag[k][cursor_x];
	}

	// plot x profiles
	plotid = PlotXY(hPanelX, idGraphX, profx_val, profx_x,   n_x_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
	SetPlotAttribute(hPanelX, idGraphX, plotid, ATTR_PLOT_LG_TEXT, "x");
	plotid = PlotXY(hPanelX, idGraphX, profx_val, profx_z,   n_x_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_GREEN);
	SetPlotAttribute(hPanelX, idGraphX, plotid, ATTR_PLOT_LG_TEXT, "z");
	plotid = PlotXY(hPanelX, idGraphX, profx_val, profx_mag, n_x_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
	SetPlotAttribute(hPanelX, idGraphX, plotid, ATTR_PLOT_LG_TEXT, "length");

	// plot z profiles
	plotid = PlotXY(hPanelZ, idGraphZ, profz_val, profz_x,   n_z_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
	SetPlotAttribute(hPanelZ, idGraphZ, plotid, ATTR_PLOT_LG_TEXT, "x");
	plotid = PlotXY(hPanelZ, idGraphZ, profz_val, profz_z,   n_z_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_GREEN);
	SetPlotAttribute(hPanelZ, idGraphZ, plotid, ATTR_PLOT_LG_TEXT, "z");
	plotid = PlotXY(hPanelZ, idGraphZ, profz_val, profz_mag, n_z_points, VAL_DOUBLE, VAL_DOUBLE, line_width, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
	SetPlotAttribute(hPanelZ, idGraphZ, plotid, ATTR_PLOT_LG_TEXT, "length");
	
	
}

// Plot velocity plots
void gui_plot_velocity()
{
	double velo_mag[n_z_points][n_x_points];
	double velo_x[n_z_points][n_x_points];
	double velo_z[n_z_points][n_x_points];
	double u_x, u_z;
	double cursor_x, cursor_z;
	int hTabIntensity, hTabProfX, hTabProfZ;
	int i, k;
	
	// get data
	for (i=0; i<n_x_points; i++) {
		for (k=0; k<n_z_points; k++) {
			// store velocity magnitude and components
			calc_velocity(i+1, k+1, &u_x, &u_z);
			velo_mag[k][i] = sqrt(u_x*u_x + u_z*u_z);
			velo_x[k][i]   = u_x;
			velo_z[k][i]   = u_z;
		}
	}
	
	// plot intensity & profiles
	GetPanelHandleFromTabPage(hPanel, PANEL_TAB_TOP,    0, &hTabIntensity);
	GetPanelHandleFromTabPage(hPanel, PANEL_TAB_BOTTOM, 1, &hTabProfX);
	GetPanelHandleFromTabPage(hPanel, PANEL_TAB_BOTTOM, 2, &hTabProfZ);
	GetGraphCursor(hTabIntensity, TABTOP0_GRAPH_VELO_COL, 1, &cursor_x, &cursor_z);

	DeleteGraphPlot(hTabIntensity, TABTOP0_GRAPH_VELO_COL, -1, VAL_DELAYED_DRAW);
	gui_plot_intensity(hTabIntensity, TABTOP0_GRAPH_VELO_COL, TABTOP0_LEGEND_VELO, velo_mag);
	if (dovectors)
		gui_plot_vectors(hTabIntensity, TABTOP0_GRAPH_VELO_COL, velo_mag, velo_x, velo_z);
	RefreshGraph(hTabIntensity, TABTOP0_GRAPH_VELO_COL);
	RefreshGraph(hTabIntensity, TABTOP0_LEGEND_VELO);
	
	DeleteGraphPlot(hTabProfX, TABBOTTOM1_GRAPH_VELO_PROF_X, -1, VAL_DELAYED_DRAW);
	DeleteGraphPlot(hTabProfZ, TABBOTTOM2_GRAPH_VELO_PROF_Z, -1, VAL_DELAYED_DRAW);
	gui_plot_profiles(hTabProfX, TABBOTTOM1_GRAPH_VELO_PROF_X, hTabProfZ, TABBOTTOM2_GRAPH_VELO_PROF_Z, velo_mag, velo_x, velo_z, cursor_x-0.5, cursor_z-0.5);
	RefreshGraph(hTabProfX, TABBOTTOM1_GRAPH_VELO_PROF_X);
	RefreshGraph(hTabProfZ, TABBOTTOM2_GRAPH_VELO_PROF_Z);
}

// Plot density plot
void gui_plot_density()
{
	double density[n_z_points][n_x_points];
	int hTab;
	unsigned int i, k;
	
	// get data
	for (i=0; i<n_x_points; i++) {
		for (k=0; k<n_z_points; k++) {
			density[k][i] = getdensity(i+1, k+1);
		}
	}

	// plot intensity
	GetPanelHandleFromTabPage (hPanel, PANEL_TAB_BOTTOM, 0, &hTab);
	DeleteGraphPlot(hTab, TABBOTTOM0_GRAPH_DENS_COL, -1, VAL_DELAYED_DRAW);
	gui_plot_intensity(hTab, TABBOTTOM0_GRAPH_DENS_COL, TABBOTTOM0_LEGEND_DENS, density);
	RefreshGraph(hTab, TABBOTTOM0_GRAPH_DENS_COL);
	RefreshGraph(hTab, TABBOTTOM0_LEGEND_DENS);
}

// Update all graphs
void gui_simulation_display()
{
	simulation_clear_block();
	gui_plot_velocity();
	gui_plot_density();
}

// Run the simulation, with or without periodic visualisation
void gui_simulation_run()
{
	unsigned int t, k;
	double u_avg, u_x, u_z;
	
	// Lock controls
	SetCtrlAttribute(hPanel, PANEL_NUM_TIMESTEPS, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_CONFIG, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_GLOBAL_DENSITY, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_VISCOSITY, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_F_X, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_F_Z, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_TOG_SC, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_G_SC, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_INIT_PERTURBATION, ATTR_DIMMED, 1);
	SetCtrlAttribute(hPanel, PANEL_WALL_DENS, ATTR_DIMMED, 1);

	// Get values from controls
	GetCtrlVal(hPanel, PANEL_NUM_TIMESTEPS,&n_t_steps);
	SetCtrlVal(hPanel, PANEL_BTN_RUN, 1);
	GetCtrlVal(hPanel, PANEL_GLOBAL_DENSITY, &global_density);
	GetCtrlVal(hPanel, PANEL_VISCOSITY, &viscosity);
	GetCtrlVal(hPanel, PANEL_F_X, &f_x);
	GetCtrlVal(hPanel, PANEL_F_Z, &f_z);
	GetCtrlVal(hPanel, PANEL_TOG_SC, &shan_chen_model);
	GetCtrlVal(hPanel, PANEL_G_SC, &g_sc);
	GetCtrlVal(hPanel, PANEL_INIT_PERTURBATION, &init_perturb);
	GetCtrlVal(hPanel, PANEL_WALL_DENS, &obstacle_density);
	ProcessDrawEvents();
	
	// Initialise simulation if needed
	if (init_needed) {
		simulation_init();
		init_needed=0;
	}
	
	message("Starting simulation, %d steps\n", n_t_steps);
	ProcessDrawEvents();
	
	
	// Run simulation
	running = 1;
	for (t=1, u_avg=0; t<=n_t_steps && running; t++) {
		
		// do one simulation step
		simulation_step();

		// Print something once in every 'n_drawing_steps' steps
		if ( t%n_drawing_steps == 0 ) {
			// Calculate the average velocity at x=n_x_points
			for (k=1; k<n_z_points+1; k++) {
				calc_velocity(n_x_points, k, &u_x, &u_z);
				u_avg += u_x;
			}
			u_avg /= n_z_points;
			// Write a message in the lower output box of the window.
			message("%10d / %-10d   average u_x(x=%d)=%e\n", t, n_t_steps, n_x_points, u_avg);
			
			// If simulation is exploding, we better stop to avoid math errors
			if (IsNotANumber(u_avg)) {
				message("Warning: simulation is exploding\n");
				running = 0;
				continue;
			}
			
			if (dovisual) {
				gui_simulation_display();
			}
			
			/* Process system events once in a while.
			 * If your GUI is unresponsive, you may want to lower this number.
			 * Large grids or other computationally intensive settings will
			 * probabely need this. */
			if (t%50==0) {
				ProcessSystemEvents();
			}
		}
	}
	
	// Write output files
	simulation_writeoutput(0);
	
	// Now enable the controls again
	SetCtrlAttribute(hPanel, PANEL_NUM_TIMESTEPS, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_CONFIG, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_GLOBAL_DENSITY, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_VISCOSITY, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_F_X, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_F_Z, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_TOG_SC, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_G_SC, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_INIT_PERTURBATION, ATTR_DIMMED, 0);
	SetCtrlAttribute(hPanel, PANEL_WALL_DENS, ATTR_DIMMED, 0);
	SetCtrlVal(hPanel, PANEL_BTN_RUN, 0);
	if (!running)
		message("Simulation aborted\n");
	
	running = 0;
}


/* A general function to print messages, can also used in simulation.c!
 * You don't need to understand how this works, you only need to know that you
 * can use something like
 * 	message("the value of counter i is now %d\n", i);
 * to put something on screen.
 * This function is using the printf formatting syntax, for a reference see e.g.
 * http://pubs.opengroup.org/onlinepubs/009695399/functions/fprintf.html
 */
void message(const char *msg, ...)
{
	char buf[4096];
	va_list ap;

	// Append message to text control and log file
	va_start(ap, msg);
	vsprintf(buf, msg, ap);
	SetCtrlVal(hPanel, PANEL_MESSAGES, buf);
	if (logfile) {
		vfprintf(logfile, msg, ap);
		fflush(logfile);
	}
	va_end(ap);
}
