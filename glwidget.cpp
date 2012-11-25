#include <QtGui>
 #include <QtOpenGL>

 #include <math.h>

 #include "glwidget.h"
 #include "qtlogo.h"
#include "patch.h"


 #ifndef GL_MULTISAMPLE
 #define GL_MULTISAMPLE  0x809D
 #endif

 GLWidget::GLWidget(QWidget *parent)
     : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
 {
     logo = 0;
     xRot = 0;
     yRot = 0;
     zRot = 0;

     qtGreen = QColor::fromCmykF(0.40, 0.0, 1.0, 0.0);
     qtPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
     g=new Geometry();
     glp=new Patch(g);
     qSetColor(glp->faceColor, QColor(128,255,0));
 }

 GLWidget::~GLWidget()
 {
 }

 QSize GLWidget::minimumSizeHint() const
 {
     return QSize(50, 50);
 }

 QSize GLWidget::sizeHint() const
 {
     return QSize(400, 400);
 }

 static void qNormalizeAngle(int &angle)
 {
     while (angle < 0)
         angle += 360 * 16;
     while (angle > 360 * 16)
         angle -= 360 * 16;
 }

 void GLWidget::setXRotation(int angle)
 {
     qNormalizeAngle(angle);
     if (angle != xRot) {
         xRot = angle;
         emit xRotationChanged(angle);
         updateGL();
     }
 }

 void GLWidget::setYRotation(int angle)
 {
     qNormalizeAngle(angle);
     if (angle != yRot) {
         yRot = angle;
         emit yRotationChanged(angle);
         updateGL();
     }
 }

 void GLWidget::setZRotation(int angle)
 {
     qNormalizeAngle(angle);
     if (angle != zRot) {
         zRot = angle;
         emit zRotationChanged(angle);
         updateGL();
     }
 }

 void GLWidget::initializeGL()
 {
     qglClearColor(qtPurple.dark());

     logo = new QtLogo(this, 64);
     logo->setColor(qtGreen.dark());

     glEnable(GL_DEPTH_TEST);
     glEnable(GL_CULL_FACE);
     glShadeModel(GL_SMOOTH);
     glEnable(GL_LIGHTING);
     glEnable(GL_LIGHT0);
     glEnable(GL_MULTISAMPLE);
     //glEnable (GL_BLEND);
     //glBlendFunc (GL_ONE, GL_ONE);
     //getCamera()->setProjectionMatrixAsPerspective(
       //          45.0f, static_cast<double>(width())/static_cast<double>(height()), -5.,5.);
     glFrustum(-1,1,1,-1,-5,5);
     static GLfloat lightPosition[4] = { 1, 1.0, 1.0, -10.0 };
     const float amb = 0.5;
         const float LightAmbient[][4]  = {  { amb, amb, amb, 1.0f },
                                             { amb, amb, amb, 1.0f }
                                         };
         const float LightDiffuse[] [4] = {  { 1.0f, 1.0f, 1.0f, 1.0f },
                                             { 1.0f, 1.0f, 1.0f, 1.0f }
                                         };
         const float LightPosition[][4] = {  { 1.0f,  4.0f, 2.0f, 0.0f },
                                             { 0.0f, 10.0f, 0.0f, 1.0f }
                                         };

         const float LightPositionb[][4] = {  { 1.0f,  1.0f, -5.0f, 0.0f },
                                             { 5.0f, 10.0f, 0.0f, 1.0f }
                                         };
     glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
     //glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient[0]);
         //glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse[0]);
         glLightfv(GL_LIGHT0, GL_POSITION, LightPosition[0]);
         glLightfv(GL_LIGHT0, GL_POSITION, LightPositionb[0]);

 }

 void GLWidget::paintGL()
 {
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     glLoadIdentity();
     glTranslatef(0, 0, -15);
     glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
     glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
     glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
     //logo->draw();

     g->loadArrays();

     glEnableClientState(GL_VERTEX_ARRAY);
     glEnableClientState(GL_NORMAL_ARRAY);


     glp->draw();

     glDisableClientState(GL_VERTEX_ARRAY);
     glDisableClientState(GL_NORMAL_ARRAY);

 }

 void GLWidget::resizeGL(int width, int height)
 {
     int side = qMin(width, height);
     glViewport((width - side) / 2, (height - side) / 2, side, side);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
 #ifdef QT_OPENGL_ES_1
     glOrthof(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
 #else
     glOrtho(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
 #endif
     glMatrixMode(GL_MODELVIEW);
 }

 void GLWidget::mousePressEvent(QMouseEvent *event)
 {
     lastPos = event->pos();
 }

 void GLWidget::mouseMoveEvent(QMouseEvent *event)
 {
     int dx = event->x() - lastPos.x();
     int dy = event->y() - lastPos.y();

     if (event->buttons() & Qt::LeftButton) {
         setXRotation(xRot + 8 * dy);
         setYRotation(yRot + 8 * dx);
     } else if (event->buttons() & Qt::RightButton) {
         setXRotation(xRot + 8 * dy);
         setZRotation(zRot + 8 * dx);
     }
     lastPos = event->pos();
 }
