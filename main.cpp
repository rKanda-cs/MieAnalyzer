//=============================================================================                              
// Copyright (C) 2013 Nyam-Erdene Odontsengel.
// Last updated by on Nov, 2013.
//=============================================================================

// Define macro for circle ratio "M_PI"
#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif

// Include header files
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <iostream>
#include <glut.h>
#include <algorithm>
#include "MIE.h"
using namespace std;

// Link library
#pragma comment( lib, "DLL_MIE.lib" )

// Load libraries of Mie solutions 
// (r,th: polar coordinate, k: wavenumber, radius: radius, index: refractive index) 
__declspec(dllimport) 
complex<double> TM_scattered_internal( double r, double th, double k, double radius, double index );
__declspec(dllimport) 
complex<double> TM_scattered_external( double r, double th, double k, double radius, double index );

// Simulation space size
const int Nx = 64;//64;//80;//128;
const int Ny = 64;//64;//80;//128;

// wavelength
double lam;
// scatterer radius
double radius;
// scatterer refractive index
double index;
// wavenumber
double k;
// scatterer center
double Cx, Cy;

// color coefficient
double colorCoef = 1;

// intensity field
double* int_f;
//double (*intensity)[Ny];						// intensity array
double (*amplitude)[Ny];						// amplitude array
// output file
//FILE* out;

#define FILENAME_BUF 50
char filename[FILENAME_BUF];
int filename_len;
double err;

#define SHOW_SCT 0
#define HIDE_SCT 1

#define SHOW_INTENSITY 0
#define SHOW_AMPLITUDE 1

// define show scatterer or not
int show_scatterer;
// define show intensity or amplitude
int show_intensity_amplitude;

enum modeEnum{
	TM_WGM,
	TE_WGM,
	TM_Mie,
	TE_Mie,
};

int CMode;			// WGM mode number for circumference direction
int RMode;			// WGM mode number for radius direction

// TM or TE mode
modeEnum mode;

// Get index
int idx( int i, int j )
{
    return ( i*Ny + j );
}

// Cylinder
bool IsInScatterer( double x, double y, double r )
{
    if ( x*x + y*y <= r*r ) return true;
    else return false;
}

// Create Color from Intensity
void ColorIntensity( double v, double* r, double* g, double* b )
{
	v = v * colorCoef;
    // Limit value (0 <= nv <= 1)
    double nv = max( 0.0, min( 1.0, abs(v) ) );
    // Get color
    if ( v >= 0.75 ) { *r = 1.0; *g = 4.0*(1.0-nv); *b = 0.0; }
    else if ( v >= 0.5 ) { *r = 4.0*(nv-0.5); *g = 1.0; *b = 0.0; }
    else if ( v >= 0.25 ) { *r = 0.0; *g = 1.0; *b = 4.0*(0.5-nv); }
    else { *r = 0.0; *g = 4.0*nv; *b = 1.0; }
}
// Create Color from Amplitude
void ColorAmplitude( double v, double* r, double* g, double* b )
{
	// Limit value (0 <= lv <= 1)
	//v = v;												// in Mie
	if(mode == TM_WGM){
		v = v / 10.0;										// in TM WGM
	}
	else if(mode == TE_WGM){
		v = v / 5.0;										// in TE WGM
	}
	double lv = max( -0.5, min( 0.5, v ) ) + 0.5;

	// Get color
	if ( v >= 0.25 ) {
		*r = 2.0*v;
		*g = 0.0;//4.0*(1.0-lv);
		*b = 0.0;
	}
	else if ( v >= 0.0 ) {
		*r = 2.0*(lv-0.5);
		*g = 0.0;//1.0;
		*b = 0.0;
	}
	else if ( v >= -0.25 ) {
		*r = 0.0;
		*g = 0.0;//1.0;
		*b = 2.0*(0.5-lv);
	}
	else {
		*r = 0.0;
		*g = 0.0;//4.0*lv;
		*b = -2.0*v;//1.0;
	}
}

// Keyboard event
void keyboard( unsigned char key, int x, int y )
{
	FILE *fp;
	// angle
	double q;
	double interpolation;

    switch ( key ) {
		// press s to change status of show_scatterer ( want to see scatterer or not )
		case 's':
			if(show_scatterer == SHOW_SCT){
				show_scatterer = HIDE_SCT;
			}
			else{
				show_scatterer = SHOW_SCT;
			}

			// Call display event
			glutPostRedisplay();
			break;

		// press a to see intensity or amplitude
		case 'a':
			if(show_intensity_amplitude == SHOW_AMPLITUDE){
				show_intensity_amplitude = SHOW_INTENSITY;
			}
			else{
				show_intensity_amplitude = SHOW_AMPLITUDE;
			}

			// Call display event
			glutPostRedisplay();
			break;

		// print to file
		case 'p':
			char filename[100];
			int filename_len;
			double coef;	// intensity radius coefficient
			// Write out Intensity around scatterer a2 = coef*radius
			double a2;

			coef = 1.2;//0.9;//1.2;
			printf ( "Insert coef = " );
			scanf_s( "%lf", &coef );
			printf ( "coef = %f\n", coef );
			a2 = coef * radius;

			if(mode == TM_WGM){
				filename_len = sprintf_s(filename, 100, "../outputFile/TM_WGM_MIE_SOLUTION_Ez_Intensity_rad%3.1f_lam%3.1f_coef%2.1f.txt", radius, lam, coef);
			}
			else if( mode == TE_WGM ){
				filename_len = sprintf_s(filename, 100, "../outputFile/TE_WGM_MIE_SOLUTION_Ey_Intensity_rad%3.1f_lam%3.1f_coef%2.1f.txt", radius, lam, coef);
			}
			else if( mode == TM_Mie ){
				filename_len = sprintf_s(filename, 100, "../outputFile/TM_MIE_SOLUTION_Ez_Intensity_rad%3.1f_lam%3.1f_coef%2.1f.txt", radius, lam, coef);
			}
			else if( mode == TE_Mie ){
				filename_len = sprintf_s(filename, 100, "../outputFile/TE_MIE_SOLUTION_Ey_Intensity_rad%3.1f_lam%3.1f_coef%2.1f.txt", radius, lam, coef);
			}
			fopen_s( &fp, filename, "w");

			q = 0.0;
			while ( q <= 180 ) {
				// distance for x and y direction from center
				double ci = a2 * cos( q * M_PI/180 );
				double cj = a2 * sin( q * M_PI/180 );
				// coordinate
				double ii = Cx + ci;
				double jj = Cy + cj;
				int i = (int)ii;
				int j = (int)jj;
				
				printf("q = %f\t", q);
				double alfa = ii - i; printf("alfa = %f\t",alfa);
				double beta = jj - j; printf("beta = %f\n",beta);

				// calculation like numerical sumilation using grid points
				interpolation = (1-alfa) * (1-beta) * int_f[idx(i,j)] + alfa * (1-beta) * int_f[idx(i+1,j)]
								+ (1-alfa) * beta * int_f[idx(i,j+1)] + alfa * beta * int_f[idx(i+1,j+1)];

			    // calculation by real distance and angle
				//if(mode == TM_WGM || mode == TM_Mie){
				//	interpolation = norm(TM_scattered_external( a2, q * M_PI/180, k, radius, index ));
				//}
				//else if( mode == TE_WGM || mode == TE_Mie){
				//	interpolation = norm(TEy_scattered_external( a2, q * M_PI/180, k, radius, index ));
				//}

				fprintf( fp, "%f\n", interpolation );
				q++;
			}
			fclose(fp);
			printf( "saved to file %s\n", filename );
			break;

		// quit
		case 'q':
			// Release memory
			delete [] int_f;
			exit( 0 );
        default:
            break;
    }
}


// Display event
void display( void )
{
	double x, y;
    double r, g, b;

	// Clear background
    glClear( GL_COLOR_BUFFER_BIT );

	// Begin to draw points
    glBegin( GL_POINTS );

	// Draw distribution
    for ( int i=0; i < Nx; i++ ) {
        for ( int j=0; j < Ny; j++ ) {
            // Get index
            int id = idx( i, j );
            // Set coordinate of grid points (windows range: x=-1...+1, y=-1...+1)
            x = -1.0 + 2*i/(double)(Nx-1);
            y = -1.0 + 2*j/(double)(Ny-1);

			// show amplitude
			if(show_intensity_amplitude == SHOW_AMPLITUDE){
				ColorAmplitude( amplitude[i][j], &r, &g, &b );
			}
			// show intensity
			else{
				// Get color of intensity (you may multiplay intensity[id] by scaling value e.g. 0.1)
				ColorIntensity( int_f[id], &r, &g, &b );
			}

			double ii = i - Cx;
            double jj = j - Cy;
            // Get distance from center
            double dis_to_center = sqrt( ii*ii + jj*jj );
			if ( dis_to_center < radius ) {
				r = r*5.0/8.0; g = g*5.0/8.0; b = b*5.0/8.0;
			}

			// show scatterer by line
			if ( (radius - 0.5 <= dis_to_center  && dis_to_center <= radius + 0.5) && show_scatterer == SHOW_SCT ){
				r = 1.0;//r*0.0/8.0;
				g = 1.0;//g*0.0/8.0;
				b = 1.0;//b*0.0/8.0;
			}

			// Set color
            glColor4d( r, g, b, 1 );
            // Draw points
            glVertex2d( x, y );
        }
    }

	// End to draw points
    glEnd();

	// Update view
    glFlush();
}

// Initialize OpenGL
void InitGL( int argc, char** argv )
{
    // window position
    glutInitWindowPosition( 10, 25 );
    // window size
    glutInitWindowSize( 600, 600 );
    // display mode
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
    // Initialize GLUT
    glutInit( &argc, argv );

    // Create window
    glutCreateWindow( "Mie Scattering analitic solution" );

    // Register events
    glutDisplayFunc( display );
    glutKeyboardFunc( keyboard );
	glPointSize( 8 );
} 

// main function
int main(int argc, char* argv[])
{
	// ##### Choose simulation mode
	//mode = TM_WGM;		// TM WGM mode
	//mode = TE_WGM;		// TE WGM mode
	mode = TM_Mie;
	//mode = TE_Mie;
	
	// ##### choose CMode and RMode
	CMode = 6; // 10, 12, ...
	RMode = 1; // 2, 3, ...

	FILE *fp;
    //-----------------------------------------------------
    // Initialize
    //-----------------------------------------------------

	// in default show scatterer, show amplitude and show scattered filed
	show_scatterer = SHOW_SCT;
	show_intensity_amplitude = SHOW_INTENSITY;

	// Center of Cylinder
	div_t divx = div( Nx-1, 2 );
	div_t divy = div( Ny-1, 2 );
	//div_t divx = div( Nx, 2 );
	//div_t divy = div( Ny, 2 );
	Cx = divx.quot;
	Cy = divy.quot;

	// refractive index
	if( mode == TM_WGM ){
	    radius = 32.0;

		Cx = divx.quot + 0.5;
		Cy = divy.quot + 0.5;

		colorCoef = 1.0 / 10.0;

		switch (CMode) {
		case 6:
			lam = 2.0*radius;

			switch (RMode) {
			case 1: index = 2.745; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;

		case 10:
			lam = radius;

			switch (RMode) {
			case 1:	index = 2.088; break;
			case 2:	index = 2.717; break;
			case 3:	index = 3.290; break;
			case 4:	index = 3.839; break;
			case 5:	index = 4.374; break;
			case 6:	index = 4.901; break;
			case 7:	index = 5.423; break;
			case 8:	index = 5.941; break;
			case 9:	index = 6.456; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;

		case 12:
			lam = 0.75*radius;

			switch (RMode) {
			case 1:	index = 1.821; break;
			case 2:	index = 2.313; break;
			case 3:	index = 2.755; break;
			case 4:	index = 3.176; break;
			case 5:	index = 3.585; break;
			case 6:	index = 3.987; break;
			case 7:	index = 4.383; break;
			case 8:	index = 4.775; break;
			case 9:	index = 5.165; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;

		default: cout << "wrong value for CMode" << endl;
		}
	}else if(mode == TE_WGM){
	    radius = 32.0;

		colorCoef = 1/5.0;

		switch (CMode) {
		case 6:
			lam = 2.0*radius;

			switch (RMode) {
			case 1: index = 2.683; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;
		/*
		case 10: // need to correct values of index
			lam = radius;

			switch (RMode) {
			case 1:	index = 2.088; break;
			case 2:	index = 2.717; break;
			case 3:	index = 3.290; break;
			case 4:	index = 3.839; break;
			case 5:	index = 4.374; break;
			case 6:	index = 4.901; break;
			case 7:	index = 5.423; break;
			case 8:	index = 5.941; break;
			case 9:	index = 6.456; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;

		case 12: // need to correct values of index
			lam = 0.75*radius;

			switch (RMode) {
			case 1:	index = 1.821; break;
			case 2:	index = 2.313; break;
			case 3:	index = 2.755; break;
			case 4:	index = 3.176; break;
			case 5:	index = 3.585; break;
			case 6:	index = 3.987; break;
			case 7:	index = 4.383; break;
			case 8:	index = 4.775; break;
			case 9:	index = 5.165; break;
			default: cout << "wrong value for RMode" << endl;
			}

			break;
		//*/
		default: cout << "wrong value for CMode" << endl;
		}
	}
	else if(mode == TM_Mie || mode == TE_Mie){			//ƒpƒ‰ƒ[ƒ^Ý’è
		radius = 8;
		lam = radius;

		index = 1.6;
		colorCoef = 1.0;

		if(mode == TM_Mie){
			Cx = divx.quot + 0.5;
			Cy = divy.quot + 0.5;
		}
	}
	else{
		printf("error: choose TM or TE mode\n");
		exit(EXIT_FAILURE);
	}

    cout << "Scatterer center-x = " << Cx << endl; 
    cout << "Scatterer center-y = " << Cy << endl; 
/*
    cout << "Grid width                 = "; cin >> Nx;
    cout << "Grid height                = "; cin >> Ny;
    cout << "Grid interval              = "; cin >>  h;
    cout << "Wavelength                 = "; cin >> lam;
    cout << "Scatterer radius           = "; cin >>  radius;
    cout << "Scatterer Refractive index = "; cin >>  index;
    cout << "Scatterer center-x         = "; cin >> Cx; 
    cout << "Scatterer center-y         = "; cin >> Cy; 
*/

    // Calculate wavenumber
    k = 2.0*M_PI / lam;

    // Get memory
    if ( (int_f = new double [Nx*Ny]) == NULL ) return 1;
	if ( (amplitude = new double[Nx][Ny]) == NULL) return 1;

    // Open output file
	if(mode == TM_WGM){
		fopen_s( &fp, "../outputFile/TM_WGM_MIE_SOLUTION_Ez_Field.txt", "w" );
	}
	else if(mode == TE_WGM){
		fopen_s( &fp, "../outputFile/TE_WGM_MIE_SOLUTION_Ey_Field.txt", "w" );
	}
	else if( mode == TM_Mie ){
		fopen_s( &fp, "../outputFile/TM_MIE_SOLUTION_Ez_Field.txt", "w");
	}
	else if( mode == TE_Mie ){
		fopen_s( &fp, "../outputFile/TE_MIE_SOLUTION_Ey_Field.txt", "w");
	}
	printf( "iteration = %d\n", ITERATION );

    //-----------------------------------------------------
    // Calculate analytic solution
    //-----------------------------------------------------
    for ( int i=0; i<Nx; i++ ) {
        for ( int j=0; j<Ny; j++ ) {
            // Get memory index
            int idx = i*Ny + j;
            // Get position from center
            double x = i - Cx;
            double y = j - Cy;
            // Get distance from center
            double r = sqrt( x*x + y*y );
            // Get angle by radian
            double th = atan2( y, x );
            // Compute intensity
            complex<double> scat;

			if(mode == TM_WGM || mode == TM_Mie){
				if ( r < radius ) scat = TM_scattered_internal( r, th, k, radius, index );
				else         scat = TM_scattered_external( r, th, k, radius, index );
			}
			else if(mode == TE_WGM || mode == TE_Mie){
				if ( r < radius ) scat = TEy_scattered_internal( r, th, k, radius, index );
				else         scat = TEy_scattered_external( r, th, k, radius, index );
			}

            // Set data
			amplitude[i][j] = scat.real();
			//if(scat.real() >=0)
			//	amplitude[i][j] = abs(scat);
			//else
			//	amplitude[i][j] = -abs(scat);
            
			int_f[idx] = norm( scat );
        }
    }

	//-----------------------------------------------------
    // Output result
    //-----------------------------------------------------

    // Output data    
    for ( int i=0; i<Ny; i++ ) {
		for ( int j=0; j<Nx; j++ ) {
			fprintf( fp, " %10.5lf", int_f[i] );
		}
		fprintf( fp, "\n", int_f[i] );
    }

    //-----------------------------------------------------
    // Finalize
    //-----------------------------------------------------

    // Close output file
    fclose( fp );

	// Initialize GL
    InitGL( argc, argv );
    // GLUT loop
    glutMainLoop();
    // Release memory
    delete [] int_f;
	delete [] amplitude;

    return 0;
}
