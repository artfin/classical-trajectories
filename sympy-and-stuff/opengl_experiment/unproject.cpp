#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>

using std::cout;
using std::endl;

void display( void )
{
	glClear( GL_COLOR_BUFFER_BIT );
	glFlush();
}

void reshape( int w, int h )
{
	glViewport( 0, 0, (GLsizei) w, (GLsizei) h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective( 45.0, (GLfloat) w / (GLfloat) h, 1.0, 100.0 );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
}

void mouse( int button, int state, int x, int y )
{
	GLint viewport[4];
	GLdouble mvmatrix[16], projmatrix[16];
	GLint realy; // Opengl y coordinate position
	GLdouble wx, wy, wz; // returned world x,y,z coords

	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			if ( state == GLUT_DOWN )
			{
				glGetIntegerv( GL_VIEWPORT, viewport );
				glGetDoublev( GL_MODELVIEW_MATRIX, mvmatrix );
				glGetDoublev( GL_PROJECTION_MATRIX, projmatrix );
				realy = viewport[3] - (GLint) y - 1;
				cout << "Coordinates at cursor are (" << x << "," << realy << ")" << endl;
				gluUnProject( (GLdouble) x, (GLdouble) realy, 0.0, mvmatrix, projmatrix, viewport, &wx, &wy, &wx );
				cout << "World coords at z = 0.0 are (" << wx << "," << wy << "," << wz << ")" << endl;
				gluUnProject( (GLdouble) x, (GLdouble) realy, 1.0, mvmatrix, projmatrix, viewport, &wx, &wy, &wz );
				cout << "World coords at z = 1.0 are (" << wx << "," << wy << "," << wz << ")" << endl;
			}
			
			break;
		case GLUT_RIGHT_BUTTON:
			if ( state == GLUT_DOWN )
			{
				exit( 0 );
			}
			break;
		default:
			break;
	}
}

int main( int argc, char** argv )
{
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize( 500, 500 );
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( argv[0] );
	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutMouseFunc( mouse );
	glutMainLoop();

	return 0;
}

