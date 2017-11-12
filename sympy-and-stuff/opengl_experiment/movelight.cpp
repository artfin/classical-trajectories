#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>

using std::endl;
using std::cout;

static int phi = 0;
static int theta = 0;

void init( void )
{
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_SMOOTH );
	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );
	glEnable( GL_DEPTH_TEST );
}

// the light position is reset after the modeling transformation
void display( void )
{
	GLfloat position[] = { 0.0, 0.0, 1.5, 1.0 };

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glPushMatrix();
	glTranslatef( 0.0, 0.0, -5.0 );

	glPushMatrix();
	glRotated( (GLdouble) phi, 1.0, 0.0, 0.0 );
	glRotated( (GLdouble) theta, 0.0, 1.0, 0.0 );
	glLightfv( GL_LIGHT0, GL_POSITION, position );

	glTranslated( 0.0, 0.0, 1.5 );
	glDisable( GL_LIGHTING );
	glColor3f( 0.0, 1.0, 1.0 );
	glutWireCube( 0.1 );
	glEnable( GL_LIGHTING );
	glPopMatrix();

	glutSolidTorus( 0.275, 0.85, 8, 15 );
	//glutSolidSphere( 1.0, 20, 16 );
	glPopMatrix();
	glFlush();
}

void reshape( int w, int h )
{
		//GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glViewport( 0, 0, (GLsizei) w, (GLsizei) h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective( 40.0, (GLfloat) w / (GLfloat) h, 1.0, 20.0 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	//glLightfv( GL_LIGHT0, GL_POSITION, light_position );
}

void keyboard( unsigned char key, int x, int y )
{
	int angle_step = 2;

	switch( key )
	{
		case 's': 
			phi = ( phi + angle_step ) % 360;
			glutPostRedisplay();
			break;
		case 'w':
			phi = ( phi - angle_step ) % 360;
			glutPostRedisplay();
			break;
		case 'd':
			theta = ( theta + angle_step ) % 360;
			glutPostRedisplay();
			break;
		case 'a':
			theta = ( theta - angle_step ) % 360;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}

int main( int argc, char** argv )
{
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH );
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
