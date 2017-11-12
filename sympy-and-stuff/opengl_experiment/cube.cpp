#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>

using std::cout;
using std::endl;

static double camera_x = 0.0;
static double camera_y = 0.0;
static double camera_z = 5.0;

static double orientation_x = 0.0;
static double orientation_y = 0.0;

void init( void )
{
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_FLAT );
}

void display( void )
{
	glClear( GL_COLOR_BUFFER_BIT );
	//glColor3f( 1.0, 1.0, 1.0 );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	gluLookAt( camera_x, camera_y, camera_z, orientation_x, orientation_y, 0.0, 0.0, 1.0, 0.0 );
	//glScalef( 1.0, 2.0, 1.0 );
	glutWireCube( 1.0 );
	glFlush();
}

void reshape( int w, int h )
{
	glViewport( 0, 0, (GLsizei) w, (GLsizei) h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glFrustum( -1.0, 1.0, -1.0, 1.0, 1.5, 20.0 );
	//gluPerspective( 60.0, 1.0, 1.5, 20.0 );
	glMatrixMode( GL_MODELVIEW );
}

void keyboard( unsigned char key, int x, int y )
{
	switch( key )
	{
		case 'w':
			camera_z -= 0.1;
			glutPostRedisplay();
			break;
		case 's':
			camera_z += 0.1;
			glutPostRedisplay();
			break;		
		case 'a':
			camera_x -= 0.1;
			glutPostRedisplay();
			break;
		case 'd':
			camera_x += 0.1;
			glutPostRedisplay();
			break;
		case 'z':
			camera_y += 0.1;
			glutPostRedisplay();
			break;
		case 'x':
			camera_y -= 0.1;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}

void process_arrows( int key, int x, int y )
{
	switch( key )
	{
		case 27: 
			break;
		// GLUT_KEY_LEFT == left arrow
		case 100:
			orientation_x -= 0.1;
			glutPostRedisplay();
			break;
		// GLUT_KEY_RIGHT == right arrow
		case 102:
			orientation_x += 0.1;
			glutPostRedisplay();
			break;
		// GLUT_KEY_UP == up arrow
		case 101:
			orientation_y += 0.1;
			glutPostRedisplay();
			break;
		// GLUT_KEY_DOWN == down arrow
		case 103:
			orientation_y -= 0.1;
			glutPostRedisplay();
			break;
	}
}

int main( int argc, char** argv )
{
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize( 800, 800 );
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( argv[0] );
	init();

	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	
	// process keys
	glutKeyboardFunc( keyboard );
	glutSpecialFunc( process_arrows );
	
	glutMainLoop();
	 
	return 0;
}
