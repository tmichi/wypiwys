/**
* main.cpp
* (c)2008 Takashi Michikawa
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#include <OpenGL/gl.h>

#include "BufferManager.h"
#include "CameraAnimation.h"
#include "Camera.h"
#include "Mesh.h"

#define PAINT_IDLE 0
#define PAINT_READY 1
#define PAINT_ING 2
//toggles
static int togglePaint = PAINT_IDLE;
static int toggleROI = 0;
static int toggleObject = 1;
static int toggleDump = 0;


//viewpoint
CameraAnimation anim;
Camera *camera;
BufferManager* bufmng;
// regions
std::vector<Vector3d> pnt;
std::vector<double> weight;

//mouse
static int left_mouse;
static int middle_mouse;
static int oldx;
static int oldy;
static std::deque<int> paintx;
static std::deque<int> painty;
//mesh
static GLuint meshid;

void initLight()
{
        // init light
        static float light0_ambient[] =  {0.1f, 0.1f, 0.1f, 1.0f};
        static float light0_diffuse[] =  {1.0f, 1.0f, 1.0f, 0.0f};
        static float light0_position[] = {0.0f, 0.0f,10.0f, 0.0f};
        static float light0_specular[] = {0.4f, 0.4f, 0.4f, 1.0f};
        glLightfv( GL_LIGHT0, GL_AMBIENT, light0_ambient );
        glLightfv( GL_LIGHT0, GL_DIFFUSE, light0_diffuse );
        glLightfv( GL_LIGHT0, GL_SPECULAR, light0_specular );
        glLightfv( GL_LIGHT0, GL_POSITION, light0_position );
        glEnable( GL_LIGHT0 );
        glEnable( GL_LIGHTING );
}



void drawPivot()
{
        GLfloat mat_ambient[]    = { 0.2f, 0.8f, 0.2f, 1.0f };
        GLfloat mat_diffuse[]    = { 0.2f, 0.8f, 0.2f, 1.0f };
        GLfloat mat_specular[]   = { 0.8f, 0.8f, 0.8f, 1.0f };
        GLfloat mat_shininess[]  = { 100.0f };
        glMaterialfv( GL_FRONT, GL_AMBIENT,  mat_ambient );
        glMaterialfv( GL_FRONT, GL_DIFFUSE,  mat_diffuse );
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular );
        glMaterialfv( GL_FRONT, GL_SHININESS,mat_shininess );
        double rad = camera->getRadius();
        glutSolidSphere( rad * 0.01 , 8, 8 );
}

void display()
{
        if( togglePaint != PAINT_IDLE ) {
                glClear( GL_DEPTH_BUFFER_BIT );
                glAccum( GL_RETURN, 1.0 );
                glBegin( GL_POINTS );
                glPointSize( 9.0 );
                glColor3f( 1,0,0 );
                for( size_t i = 0 ; i < paintx.size() ; i++ ) {
                        glVertex2f( paintx[i],painty[i] );
                }
                glEnd();
        } else {
                glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
                glMatrixMode( GL_MODELVIEW );
                glLoadIdentity();
                camera->applyTransform();

                if( toggleObject ) {
                        glCallList( meshid );
                }
                if( toggleROI ) {
                        for( size_t i = 0 ; i < pnt.size() ; i++ ) {
                                glPushMatrix();
                                glTranslated( pnt.at( i ).x, pnt.at( i ).y, pnt.at( i ).z );
                                drawPivot();
                                glPopMatrix();
                        }
                }
                glAccum( GL_LOAD,1.0 );
        }
        glutSwapBuffers();
}

static void
reshape( int width, int height )
{
        glViewport( 0, 0, width, height );
        camera->updatePerspective();
        glutPostRedisplay();
}

void reset( void )
{

        left_mouse = 0;
        middle_mouse = 0;

        camera->init();
        camera->updatePerspective();
}

void
keyboard( unsigned char c, int x, int y )
{
        x; y; // for C4100

        if( c == 'b' ) {
                camera->zoom( 0.1 );
        }
        if( c == 'f' ) {
                camera->zoom( -0.1 );
        }
        if( c == 0x1b ) {

                bufmng->startViewMode();
                camera->updatePerspective();
                togglePaint = PAINT_IDLE;
        }
}

void
motion( int x, int y )
{
        int viewport[4];
        glGetIntegerv( GL_VIEWPORT,viewport );
        const int H = viewport[3];

        if ( left_mouse ) {
                if( togglePaint == PAINT_ING ) {
                        paintx.push_back( x );
                        painty.push_back( H-y );
                } else if ( oldx != x || oldy != y ) {
                        camera->rotate( oldx, oldy, x, y );
                }
        } else if ( middle_mouse ) {
                camera->zoom( oldx - x );
                camera->updatePerspective();
        }
        oldx = x;
        oldy = y;
        glutPostRedisplay();
}

void init_mouse()
{
        oldx = -100;
        oldy = -100;
        paintx.clear();
        painty.clear();
}
void set3DViewMode()
{
        glEnable( GL_LIGHTING );
        glDepthFunc( GL_LESS );
        camera->updatePerspective();
}
void
mouse( int b, int s, int x, int y )
{

        oldx = x;
        oldy = y;

        //ボタンを押している
        if ( s == GLUT_DOWN ) {
                switch ( b ) {
                case GLUT_LEFT_BUTTON:
                        left_mouse = GL_TRUE;

                        if( togglePaint == PAINT_READY ) {
                                togglePaint = PAINT_ING;
                        }
                        break;
                case GLUT_MIDDLE_BUTTON:
                        middle_mouse = GL_TRUE;
                        break;
                }
        } else { //ボタンを離した
                if ( left_mouse ) {
                        left_mouse = GL_FALSE;
                        if( togglePaint == PAINT_ING ) {
                                bufmng->readMask();
                                if( toggleDump == 1 ) dumpDepth();
                                bufmng->getPoints( pnt,weight );
                                int framenum = ( int ) ( 10*3 ); // 5sec. animation
                                anim.init( camera, pnt, weight, framenum );

                                bufmng->startViewMode();
                                camera->updatePerspective();
                                togglePaint = PAINT_IDLE;
                                init_mouse();
                        }
                } else if ( middle_mouse ) {
                        middle_mouse = GL_FALSE;
                }
        }
        glutPostRedisplay();
}

void
animate( void )
{
        if( anim.animate() )glutPostRedisplay();
}


void
init( void )
{
        initLight();
        ::glClearAccum( 0,0,0,0 );
        glEnable( GL_CULL_FACE );
        glCullFace( GL_BACK );
        glEnable( GL_DEPTH_TEST );
        glDepthFunc( GL_LESS );
        glClearColor( 0,0,0,1 );
        glShadeModel( GL_FLAT );
        glPointSize( 8.0 );

        togglePaint = PAINT_IDLE;
        reset();
        init_mouse();
}



void
menu( int choice )
{
        switch( choice ) {
        case 0:
                togglePaint = PAINT_READY;
                pnt.clear();
                weight.clear();
                if( toggleDump == 1 ) dumpDepth();

                bufmng->readDepth();
                bufmng->startPaintMode();
                break;
        case 1:
                toggleObject = 1- toggleObject;
                cerr<<"Toggle Object ";
                if( toggleObject ) cerr<<"on";
                else cerr<<"off";
                cerr<<endl;
                break;
        case 2:

                toggleROI = 1- toggleROI;
                cerr<<"Toggle ROI ";
                if( toggleROI ) cerr<<"on";
                else cerr<<"off";
                cerr<<endl;
                break;
        case 3:
                cerr<<"Reset"<<endl;
                reset();
                break;
        case 4:
                toggleDump = 1- toggleDump;
                cerr<<"Toggle Dump ";
                if( toggleDump ) cerr<<"on";
                else cerr<<"off";
                break;
        case 5:
                dumpImage();
                cerr<<"image dumpped"<<endl;
                display();
                break;
        case 6:
                exit( 0 );
                break;
        }

        glutPostRedisplay();
}

int
main( int argc, char **argv )
{
        int W = 640;
        int H = 480;

        glutInit( &argc, argv );
        glutInitWindowSize( W,H );
        glutInitDisplayMode( GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE |GLUT_ACCUM );
        glutCreateWindow( "wypiwys (c) 2008 Takashi Michikawa" );
        glutReshapeFunc( reshape );
        glutDisplayFunc( display );
        glutKeyboardFunc( keyboard );
        glutMotionFunc( motion );
        glutMouseFunc( mouse );
        glutIdleFunc( animate );

        //メニュー作成
        glutCreateMenu( menu );
        glutAddMenuEntry( "Start Painting", 0 );
        glutAddMenuEntry( "Toggle Object", 1 );
        glutAddMenuEntry( "Toggle ROI", 2 );
        glutAddMenuEntry( "Reset View", 3 );
        glutAddMenuEntry( "toggle Dump", 4 );
        glutAddMenuEntry( "Dump Image", 5 );
        glutAddMenuEntry( "Quit", 6 );
        glutAttachMenu( GLUT_RIGHT_BUTTON );


        if( argc == 2 ) {
                double init_rad = 0;
                Vector3d init_center;

                if( checkFormat( argv[1], "obj" ) ) {
                        meshid = read_mesh( argv[1], init_rad, init_center );
                } else if( checkFormat( argv[1], "xyz" ) ) {
                        meshid = read_pointset( argv[1], init_rad, init_center );
                } else {
                        return -1;
                }
                camera = new Camera( init_rad, init_center );
                camera->setInit(); // set current state as initaial state
        } else {
                fprintf( stderr, "Usage : %s input.obj\n",argv[0] );
                return -1;
        }

        init();

        bufmng = new BufferManager();
        glutMainLoop();
        delete bufmng;
        delete camera;
        return 0;
}
