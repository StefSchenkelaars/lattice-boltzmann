/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/* Copyright (c) National Instruments 2013. All Rights Reserved.          */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                            1       /* callback function: end */
#define  PANEL_NUM_TIMESTEPS              2       /* control type: numeric, callback function: control_changed */
#define  PANEL_CONFIG                     3       /* control type: ring, callback function: control_changed */
#define  PANEL_TOG_DRAW                   4       /* control type: radioButton, callback function: control_changed */
#define  PANEL_NUM_DRAWPERIOD             5       /* control type: numeric, callback function: control_changed */
#define  PANEL_BTN_RESET                  6       /* control type: command, callback function: button_pressed */
#define  PANEL_BTN_DISPLAY                7       /* control type: command, callback function: button_pressed */
#define  PANEL_BTN_INIT                   8       /* control type: command, callback function: button_pressed */
#define  PANEL_BTN_RUN                    9       /* control type: textButton, callback function: button_pressed */
#define  PANEL_BTN_DISPLAY_FILE           10      /* control type: command, callback function: button_pressed */
#define  PANEL_TAB_TOP                    11      /* control type: tab, callback function: (none) */
#define  PANEL_CURSOR_X                   12      /* control type: numeric, callback function: cursor_changed */
#define  PANEL_CURSOR_Z                   13      /* control type: numeric, callback function: cursor_changed */
#define  PANEL_TAB_BOTTOM                 14      /* control type: tab, callback function: (none) */
#define  PANEL_MESSAGES                   15      /* control type: textBox, callback function: (none) */
#define  PANEL_TEXTMSG                    16      /* control type: textMsg, callback function: (none) */
#define  PANEL_DECORATION                 17      /* control type: deco, callback function: (none) */
#define  PANEL_TOG_VEC                    18      /* control type: radioButton, callback function: control_changed */
#define  PANEL_G_SC                       19      /* control type: numeric, callback function: control_changed */
#define  PANEL_F_X                        20      /* control type: numeric, callback function: control_changed */
#define  PANEL_F_Z                        21      /* control type: numeric, callback function: control_changed */
#define  PANEL_GLOBAL_DENSITY             22      /* control type: numeric, callback function: control_changed */
#define  PANEL_VISCOSITY                  23      /* control type: numeric, callback function: control_changed */
#define  PANEL_INIT_PERTURBATION          24      /* control type: numeric, callback function: control_changed */
#define  PANEL_TOG_SC                     25      /* control type: radioButton, callback function: control_changed */
#define  PANEL_WALL_DENS                  26      /* control type: numeric, callback function: control_changed */

     /* tab page panel controls */
#define  TABBOTTOM0_GRAPH_DENS_COL        2       /* control type: graph, callback function: (none) */
#define  TABBOTTOM0_LEGEND_DENS           3       /* control type: graph, callback function: (none) */

     /* tab page panel controls */
#define  TABBOTTOM1_GRAPH_VELO_PROF_X     2       /* control type: graph, callback function: (none) */

     /* tab page panel controls */
#define  TABBOTTOM2_GRAPH_VELO_PROF_Z     2       /* control type: graph, callback function: (none) */

     /* tab page panel controls */
#define  TABTOP0_GRAPH_VELO_COL           2       /* control type: graph, callback function: cursor_changed */
#define  TABTOP0_LEGEND_VELO              3       /* control type: graph, callback function: (none) */


     /* Control Arrays: */

          /* (no control arrays in the resource file) */


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK button_pressed(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK control_changed(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK cursor_changed(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK end(int panel, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
