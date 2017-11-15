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
using std::get;

// vectors of coordinates
vector<tuple<double, double, double>> he_coordinates;
vector<tuple<double, double, double>> ar_coordinates;

static int cur_coord = 0;

// angle of rotation for the camera direction
float angle = 0.0;

// actual vector representing the camera's direction
float lx = 0.0f, lz = 1.0f;

// XZ position of the camera
float x = 0.0f, z = 0.0f;

// the key states
// these variables will be zero the no key is pressed
float deltaAngle = 0.0f;
float deltaMove = 0;

GLfloat whiteSpecularLight[] = { 1.0, 1.0, 1.0 };
GLfloat whiteDiffuseLight[] = { 1.0, 1.0, 1.0 };

void light( void )
{
	glLightfv( GL_LIGHT0, GL_SPECULAR, whiteSpecularLight );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, whiteDiffuseLight );
}

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

void computePos( float deltaMove )
{
	x += deltaMove * lx * 0.1f;
	z += deltaMove * lz * 0.1f;
}

void computeDir( float deltaAngle )
{
	angle += deltaAngle;
	lx = sin( angle );
	lz = - cos( angle );
}

void print_text( const int &x, const int &y, const int &z, const string msg )
{
	glColor3f( 1.0f, 1.0f, 1.0f );
	glLoadIdentity();
	glRasterPos3f( x, y, z );

	for( int i = 0; i < msg.length(); i++ )
	{
			//cout << "msg[" << i << "] = " << msg[i] << endl;
		glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, msg[i] );
	}	
}

void to_string( const double &x, string &msg )
{
	std::ostringstream convert;
	convert << x;
	msg = convert.str();
}

void renderScene( void )
{
	if ( deltaMove )
	{
		computePos( deltaMove );
	}
	if ( deltaAngle )
	{
		computeDir( deltaAngle );
	}

	// clear color and depth buffers
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// reset transformations
	glLoadIdentity();

	light();

	// set the camera
	gluLookAt( x, 0.0f, z,
			   x + lx, 1.0f, z + lz,
			   0.0f, 1.0f, 0.0f 
			 );

	glTranslatef( get<0>( he_coordinates[cur_coord] ), get<0>( he_coordinates[cur_coord]), get<2>( he_coordinates[cur_coord] ) );
	glColor3f( 1.0, 0.0, 0.0 );
	glutSolidSphere( 0.05f, 10, 10 );

	glLoadIdentity();

	glTranslatef( get<0>( ar_coordinates[cur_coord]), get<1>( ar_coordinates[cur_coord]), get<2>( ar_coordinates[cur_coord]) );
	glColor3f( 0.0, 1.0, 0.0 );
	glutSolidSphere( 0.05f, 10, 10 );

	glColor3f( 1.0f, 1.0f, 1.0f );
	GLfloat axe_width = 3.0;

	string x_camera_coord, z_camera_coord;
	to_string( x, x_camera_coord );
	to_string( z, z_camera_coord );

	string camera_coords = "X: " + x_camera_coord + "; Z: " + z_camera_coord;
	//cout << "camera_coords: " << camera_coords << endl;

	print_text( x + 1.0, 0.0, z + 1.0, camera_coords );

	//glLoadIdentity();
	glLineWidth( (GLfloat) axe_width );
	glBegin( GL_LINES );
		glVertex3f( 0.0f, 0.0f, -50.0f );
		glVertex3f( 0.0f, 0.0f, 50.0f );
	glEnd();

	glLoadIdentity();
	glLineWidth( (GLfloat) axe_width );
	glBegin( GL_LINES );
		glVertex3f( 50.0f, 0.0f, 0.0f );
		glVertex3f( -50.0f, 0.0f, 0.0f );
	glEnd();

	glutSwapBuffers();
}

void pressKey( int key, int xx, int yy )
{
	switch( key )
	{
		case GLUT_KEY_LEFT: 
				deltaAngle = -0.01f; 
				break;
		case GLUT_KEY_RIGHT: 
				deltaAngle = 0.01f; 
				break;
		case GLUT_KEY_UP: 
				deltaMove = 0.5f; 
				break;
		case GLUT_KEY_DOWN:
				deltaMove = -0.5f;
				break;
	}
}

void releaseKey( int key, int x, int y )
{
	switch( key )
	{
		case GLUT_KEY_LEFT:
		case GLUT_KEY_RIGHT: deltaAngle = 0.0f; break;
		case GLUT_KEY_UP:
		case GLUT_KEY_DOWN: deltaMove = 0; break;
	}
}

void keyboard( unsigned char key, int x, int y )
{
	switch( key )
	{
		case 'n':
			cur_coord ++;
			break;
	}
}

int main( int argc, char* argv[] )
{
	read_coordinates( he_coordinates, ar_coordinates, "coords.txt" );

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
	glutKeyboardFunc( keyboard );
	glutSpecialFunc( pressKey );
	
	// advanced keyboard functions
	glutIgnoreKeyRepeat( 1 );
	glutSpecialUpFunc( releaseKey );


	// opengl init
	glEnable( GL_DEPTH_TEST );
	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );
	glEnable(GL_COLOR_MATERIAL);

	glutMainLoop();

	return 0;
}
