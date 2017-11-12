#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>

using std::cout;
using std::endl;

static int year = 0, day = 0;
static int moon_year = 0, moon_day = 0;
static double porbit_radius = 5.0;
static double morbit_radius = 0.5;

static double pradius = 0.15;
static double mradius = 0.05;

static int porbit_speed = 1;
static int morbit_speed = 3;

static int protate_speed = 5;
static int mrotate_speed = 5;

static double z_position = 10.0;

void init( void )
{
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_FLAT );
}

void display( void )
{
	glClear( GL_COLOR_BUFFER_BIT );
	glColor3f( 1.0, 1.0, 1.0 );

	glPushMatrix();
	glutWireSphere( 1.0, 20, 16 ); // sun
	glRotatef( (GLfloat) year, 0.0, 1.0, 0.0 );
	glTranslatef( porbit_radius, 0.0, 0.0 );
	glRotatef( (GLfloat) day, 0.0, 1.0, 0.0 );
	glutWireSphere( 0.2, 10, 8 ); // smaller planet
	//glPopMatrix();

	//glPushMatrix();
	glRotatef( (GLfloat) moon_year, 0.0, 1.0, 0.0 );
	glTranslatef( morbit_radius, 0.0, 0.0 );
	glRotatef( (GLfloat) moon_day, 0.0, 1.0, 0.0 );
	glutWireSphere( mradius, 10, 8 ); // moon
	glPopMatrix();

	glutSwapBuffers();
}

void reshape( int w, int h )
{
	glViewport( 0, 0, (GLsizei) w, (GLsizei) h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity ();
	gluPerspective( 60.0, (GLfloat) w / (GLfloat) h, 1.0, 20.0 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	gluLookAt( 0.0, 0.0, z_position, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );
}

void keyboard( unsigned char key, int x, int y )
{
	switch( key )
	{
		case 'd':
			day = ( day + 10 ) % 360;
			moon_day = ( moon_day + 10 ) % 360;
			glutPostRedisplay();
			break;
		case 'D':
			day = ( day - 10 ) % 360;
			moon_day = ( moon_day - 10 ) % 360;
			glutPostRedisplay();
			break;
		case 'y':
			year = ( year + 5 ) % 360;
			moon_year = ( moon_year + 5 ) % 360;
			glutPostRedisplay();
			break;
		case 'Y':
			year = ( year - 5 ) % 360;
			moon_year = ( moon_year - 5 ) % 360;
			glutPostRedisplay();
			break;
		case 't':
			// planet
			year = ( year + porbit_speed ) % 360;
			day = ( day + protate_speed ) % 360;

			// moon
			moon_year = ( moon_year + morbit_speed) % 360;
			moon_day = ( day + mrotate_speed ) % 360;
			glutPostRedisplay();
			break;
		case 'a':
			// move frame back
			z_position += 0.5;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}

int main( int argc, char** argv )
{
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
	glutInitWindowSize( 500, 500 );
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( argv[0] );
	init();
	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutKeyboardFunc( keyboard );
	glutMainLoop();

	return 0;
}
