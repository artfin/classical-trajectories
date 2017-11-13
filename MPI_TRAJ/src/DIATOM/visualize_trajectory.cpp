#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using std::cout;
using std::endl;

using std::vector;
using std::tuple;
using std::string;

// angle of rotation for the camera direction
float angle = 0.0;

// actual vector representing the camera's direction
float lx = 0.0f, lz = -1.0f;

// XZ position of the camera
float x = 0.0f, z = 5.0f;

void read_coordinates( vector<tuple<double, double, double>> &he_coords,
					   vector<tuple<double, double, double>> &ar_coords,
					   string filename )
{
	std::ifstream input( filename );
	
	string line;

	tuple<double, double, double> he_temp;
	tuple<double, double, double> ar_temp;

	if ( input.is_open() )
	{
		while ( std::getline( input, line ) )
		{
			std::stringstream iss( line );

			iss >> std::get<0>( he_temp );
			iss >> std::get<1>( he_temp );
			iss >> std::get<2>( he_temp );

			iss >> std::get<0>( ar_temp );
			iss >> std::get<1>( ar_temp );
			iss >> std::get<2>( ar_temp );

			he_coords.push_back( he_temp );
			ar_coords.push_back( ar_temp );	
		}
	
		input.close();
	} 
	else 
	{
		cout << "Unable to open the file." << endl;
	}
}

void changeSize( int w, int h )
{
	float ratio = w / h;

	// use the projection 
	glMatrixMode( GL_PROJECTION );

	// reset matrix
	glLoadIdentity();

	// set the viewport to be the the entire window
	glViewport( 0, 0, w, h );

	// set the correct perspective
	gluPerspective( 45.0f, ratio, 0.1f, 100.0f );

	// get back to the modelview
	glMatrixMode( GL_MODELVIEW );
}

void renderScene( void )
{
	// clear color and depth buffers
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// reset transformations
	glLoadIdentity();
	// set the camera
	gluLookAt( x, 1.0f, z,
			   x + lx, 1.0f, z + lz,
			   0.0f, 1.0f, 0.0f 
			 );
	glTranslatef( 1.0, 0.0, 1.0 );
	glutSolidSphere( 0.05f, 10, 10 );
	
	glTranslatef( 3.0, 0.0, 3.0 );
	glutSolidSphere( 0.05f, 10, 10 );

	glutSwapBuffers();
}

void processSpecialKeys( int key, int xx, int yy )
{
	float fraction = 1.0f;

	switch( key )
	{
		case GLUT_KEY_LEFT:
			angle -= 0.01f;
			lx = sin( angle );
			lz = -cos( angle );
			break;
		case GLUT_KEY_RIGHT:
			angle += 0.01f;
			lx = sin( angle );
			lz = -cos( angle );
			break;
		case GLUT_KEY_UP:
			x += lx * fraction;
			z += lz * fraction;
			break;
		case GLUT_KEY_DOWN:
			x -= lx * fraction;
			z -= lz * fraction;
			break;
	}
}

int main( int argc, char* argv[] )
{
	vector<tuple<double, double, double>> he_coords;
	vector<tuple<double, double, double>> ar_coords;

	read_coordinates( he_coords, ar_coords, "coords.txt" );

	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize( 400, 400 );
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( "Trajectory" );

	// register callbacks
	glutDisplayFunc( renderScene );
	glutReshapeFunc( changeSize );
	glutIdleFunc( renderScene );

	// processing keys callback
	glutSpecialFunc( processSpecialKeys );

	// opengl init
	glEnable( GL_DEPTH_TEST );

	glutMainLoop();

	return 0;
}
