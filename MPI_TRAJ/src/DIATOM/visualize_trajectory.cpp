#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;

using std::vector;
using std::tuple;
using std::string;

static double camera_x = 0.0;
static double camera_y = 0.0;
static double camera_z = 5.0;

static double orientation_x = 0.0;
static double orientation_y = 0.0;

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

void init( void )
{
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 100.0 };
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_SMOOTH );

	glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular );
	glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess );
	glLightfv( GL_LIGHT0, GL_POSITION, light_position );

	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );
	glEnable( GL_DEPTH_TEST );
}

void display( void )
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glColor3f( 1.0, 1.0, 1.0 );
	
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
	gluLookAt( camera_x, camera_y, camera_z, orientation_x, orientation_y, 0.0, 0.0, 1.0, 0.0 );
	
	cout << "camera_z: " << camera_z << endl;

	//glScalef( 1.0, 1.0, 1.0 );
	//glutSolidSphere( 1.0, 20, 16 );
   	glutWireCube( 1.0 );
	glFlush();	
}

void reshape( int w, int h )
{
	glViewport( 0, 0, (GLsizei) w, (GLsizei) h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	if ( w <= h )
	{
		glOrtho( -1.5, 1.5, -1.5 * (GLfloat) h / (GLfloat) w, 
			1.5 * (GLfloat) h / (GLfloat) w, -10.0, 10.0 );
	}
	else
	{
		glOrtho( -1.5 * (GLfloat) w / (GLfloat) h,
			1.5 * (GLfloat) w / (GLfloat) h, -1.5, 1.5, -10.0, 10.0 );
	}

	//glMatrixMode( GL_MODELVIEW );
	//glLoadIdentity();
}

void keyboard( unsigned char key, int x, int y )
{
	switch( key )
	{
		case 'w':
			cout << "Button 'w' has been detected" << endl;
			camera_z -= 0.1;
			glutPostRedisplay();
			break;
		case 's':
			cout << "Button 's' has been detected" << endl;
			camera_z += 0.1;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}


int main( int argc, char* argv[] )
{
	vector<tuple<double, double, double>> he_coords;
	vector<tuple<double, double, double>> ar_coords;

	read_coordinates( he_coords, ar_coords, "coords.txt" );

	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize( 800, 800 );
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( argv[0] );

	init();

	glutDisplayFunc( display );
	glutReshapeFunc( reshape );

	glutKeyboardFunc( keyboard );

	glutMainLoop();

	return 0;
}
