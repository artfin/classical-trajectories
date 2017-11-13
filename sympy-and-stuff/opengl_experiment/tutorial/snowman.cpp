#include <GL/gl.h>
#include <GL/glut.h>
#include <cmath>

// angle of rotation for the camera direction
float angle = 0.0;

// actual vector representing the camera's direction
float lx = 0.0f, lz = -1.0f;

// XZ position of the camera
float x = 0.0f, z = 5.0f;

// the key states
// these variables will be zero when no is pressed
float deltaAngle = 0.0f;
float deltaMove = 0;

GLfloat whiteSpecularLight[] = { 0.0, 1.0, 1.0 };
//GLfloat blackAmbientLight[] = { 1.0, 0.78, 1.0, 10.0 };
GLfloat whiteDiffuseLight[] = { 1.0, 1.0, 1.0 };

void light( void )
{
	glLightfv( GL_LIGHT0, GL_SPECULAR, whiteSpecularLight );
	//glLightfv( GL_LIGHT0, GL_AMBIENT, blackAmbientLight );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, whiteDiffuseLight );
}

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

void drawSnowMan( void )
{
	glColor3f( 1.0f, 1.0f, 1.0f );

	// draw body
	glTranslatef( 0.0f, 0.75f, 0.0f );
	glutSolidSphere( 0.75f, 20, 20 );

	// draw head
	glTranslatef( 0.0f, 1.0f, 0.0f );
	glutSolidSphere( 0.25f, 20, 20 );

	// draw eyes
	glPushMatrix();
	glColor3f( 0.0f, 0.0f, 0.0f );
	glTranslatef( 0.05f, 0.10f, 0.18f );
	glutSolidSphere( 0.05f, 10, 10 );
	glTranslatef( -0.1f, 0.0f, 0.0f );
	glutSolidSphere( 0.05f, 10, 10 );
	glPopMatrix();

	// draw nose
	glColor3f( 1.0f, 0.5f, 0.5f );
	glutSolidCone( 0.08f, 0.5f, 10, 2 );
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
	lz = -cos( angle );
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
	// first: position of camera
	// second: look at point = line of sight + camera position
	// third: up vector (the camera is not tilted)
	gluLookAt( x, 1.0f, z,
			   x + lx, 1.0f, z + lz,
			   0.0f, 1.0f, 0.0f 
			 );

	// draw ground
	glColor3f( 0.9f, 0.9f, 0.9f );
	glBegin( GL_QUADS );
		glVertex3f( -100.0f, 0.0f, -100.0f );
		glVertex3f( -100.0f, 0.0f,  100.0f );
		glVertex3f(  100.0f, 0.0f,  100.0f );
		glVertex3f(  100.0f, 0.0f, -100.0f );
	glEnd();

	// draw 36 snowmen
	for ( int i = -3; i < 3; i++ )
	{
		for ( int j = -3; j < 3; j++ )
		{
			glPushMatrix();
			glTranslatef( i * 10.0, 0.0, j * 10.0 );
			drawSnowMan();
			glPopMatrix();
		}
	}

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

void processSpecialKeys( int key, int xx, int yy )
{
	float fraction = 0.1f;

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

int main( int argc, char** argv )
{
	// init GLUT and crate window
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA );
	glutInitWindowPosition( 100, 100 );
	glutInitWindowSize( 320, 320 );
	glutCreateWindow( "Lighthouse3D" );

	// register callbacks
	glutDisplayFunc( renderScene );
	glutReshapeFunc( changeSize );
	glutIdleFunc( renderScene );
	//glutKeyboardFunc( processNormalKeys );
	glutSpecialFunc( pressKey );

	// advanced keyboard functions
	glutIgnoreKeyRepeat( 1 );
	glutSpecialUpFunc( releaseKey );

	// openg init
	glEnable( GL_DEPTH_TEST );
	
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 100.0 };
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 1.0 };
	
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_SMOOTH );

	glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular );
	glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess );
	glLightfv( GL_LIGHT0, GL_POSITION, light_position );

	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );
	glEnable(GL_COLOR_MATERIAL);
	// enter glut event processing cycle
	glutMainLoop();

	return 0;
}
