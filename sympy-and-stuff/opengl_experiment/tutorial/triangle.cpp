#include <GL/glut.h>
#include <stdlib.h>

float angle = 0.0f;

float red = 1.0f, blue = 1.0f, green = 1.0f;

void changeSize( int w, int h )
{
	float ratio = w / h;

	// use the projection matrix
	glMatrixMode( GL_PROJECTION );

	// reset matrix
	glLoadIdentity();

	// set the viewport to be the entire window
	glViewport( 0, 0, w, h );

	// set the correct perspective
	gluPerspective( 45.0f, ratio, 0.1f, 100.0f );

	// get back to the modelview 
	glMatrixMode( GL_MODELVIEW );
}

// process every key that has ASCII number
// key: the key that has been pressed (in ASCII)
// x, y -- position of mouse when that event happened
void processNormalKeys( unsigned char key, int x, int y )
{
	// exit when Escape is pressed
	if ( key == 27 )
	{
		exit( 0 );
	}
}

void processSpecialKeys( int key, int x, int y )
{
	switch( key )
	{
		case GLUT_KEY_F1:
			red = 1.0;
			green = 0.0;
			blue = 0.0;
			break;
		case GLUT_KEY_F2:
			red = 0.0;
			green = 1.0;
			blue = 0.0;
			break;
		case GLUT_KEY_F3:
			red = 0.0;
			green = 0.0;
			blue = 1.0;
			break;
	}
}


void renderScene( void )
{
	// clear color and depth buffers
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// reset transformations
	glLoadIdentity();
	// set the camera
	gluLookAt( 0.0f, 0.0f, 10.0f,
  			   0.0f, 0.0f, 0.0f,
			   0.0f, 1.0f, 0.0f );
	glRotatef( angle, 0.0f, 1.0f, 0.0f );

	// the function responsible for setting the color
	glColor3f( red, green, blue );

	glBegin( GL_TRIANGLES );
		glVertex3f( -2.0f, -2.0f, 0.0 );
		glVertex3f(  2.0f, 0.0f, 0.0 );
		glVertex3f(  0.0f, 2.0f, 0.0 );
	glEnd();

	angle += 0.1f;

	glutSwapBuffers();
}

int main( int argc, char** argv )
{
	//init GLUT and create window
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA );
	glutInitWindowPosition( 100, 100 );
	glutInitWindowSize( 320, 320 );
	glutCreateWindow( "Lighthouse3D" );

	// register callbacks
	glutDisplayFunc( renderScene );
	glutReshapeFunc( changeSize );

	// processing keyboard events
	glutKeyboardFunc( processNormalKeys );
	glutSpecialFunc( processSpecialKeys );

	// the idle func registration
	glutIdleFunc( renderScene );

	// enter Glut event processing cycle
	glutMainLoop();

	return 0;
}
