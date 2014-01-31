#include "interface.h"
#include <stdio.h>
#include <math.h>
#include <QTime>

#ifdef WIN32
#include "glee.h"
#define glActiveTexture glActiveTextureARB
#define glMultiTexCoord2d glMultiTexCoord2dARB
#endif

/*#include <iostream>
using namespace std;*/

vector4 GraphicsWidget::floattovec4(double x)
{
	double a,b,c,d;
	double logx = floor(log2(x));
	double xm = x/pow(2.0,logx) - 1.0;
	
	d = logx + 128.0;
	a = floor(xm*256.0);
	xm = xm*256.0 - a;
	b = floor(xm*256.0);
	xm = xm*256.0 - b;
	c = floor(xm*256.0);
	
	return vector4(a,b,c,d);
}

#define NRS_2 512
#define NBS_2 128
#define NRL_2 512
#define NBL_2 1024

GraphicsWidget::GraphicsWidget(QWidget* parent, Simulation* _sim, bool _doppler, double _fov)
	: QGLWidget(parent)
{
	sim = _sim;
	doppler = _doppler;
	fov = _fov;
	paused = true;
	
	FILE* fp;
	double* temp;
	fp = fopen("data/large_b.dat", "rb");
	fread(&nb, sizeof(int), 1, fp);
	fread(&nr, sizeof(int), 1, fp);
	large_b = new unsigned char[NRL_2*NBL_2*4];
	int i,j;
	for(i=0; i<nb; i++)
	{
		temp = new double[nr];
		fread(temp, sizeof(double), nr, fp);
		for(j=0; j<nr; j++)
		{
			vector4 a = floattovec4(temp[j]);
			large_b[i*NRL_2*4+j*4] = (unsigned char)a[0];
			large_b[i*NRL_2*4+j*4+1] = (unsigned char)a[1];
			large_b[i*NRL_2*4+j*4+2] = (unsigned char)a[2];
			large_b[i*NRL_2*4+j*4+3] = (unsigned char)a[3];
		}
		for(j=nr; j<NRL_2; j++)
		{
			large_b[i*NRL_2*4+j*4] = large_b[i*NRL_2*4+j*4+1] = large_b[i*NRL_2*4+j*4+2] = large_b[i*NRL_2*4+j*4+3] = 0;
		}
		delete[] temp;
	}
	for(i=nb; i<NBL_2; i++)
		for(j=0; j<NRL_2; j++)
			large_b[i*NRL_2*4+j*4] = large_b[i*NRL_2*4+j*4+1] = large_b[i*NRL_2*4+j*4+2] = large_b[i*NRL_2*4+j*4+3] = 0;
	fclose(fp);
	
	fp = fopen("data/small_b.dat", "rb");
	fread(&nbs, sizeof(int), 1, fp);
	fread(&nrs, sizeof(int), 1, fp);
	small_b = new unsigned char[NRS_2*NBS_2*4];
	
	for(i=0; i<nbs; i++)
	{
		temp = new double[nrs];
		fread(temp, sizeof(double), nrs, fp);
		for(j=0; j<nrs; j++)
		{
			vector4 a = floattovec4(temp[j]);
			small_b[i*NRS_2*4+j*4] = (unsigned char)a[0];
			small_b[i*NRS_2*4+j*4+1] = (unsigned char)a[1];
			small_b[i*NRS_2*4+j*4+2] = (unsigned char)a[2];
			small_b[i*NRS_2*4+j*4+3] = (unsigned char)a[3];
		}
		for(j = nrs; j < NRS_2; j++)
		{
			small_b[i*NRS_2*4+j*4] = small_b[i*NRS_2*4+j*4+1] = small_b[i*NRS_2*4+j*4+2] = small_b[i*NRS_2*4+j*4+3] = 0;
		}
		delete[] temp;
	}
	for(i=nbs; i<NBS_2; i++)
		for(j=0; j<NRS_2; j++)
			small_b[i*NRS_2*4+j*4] = small_b[i*NRS_2*4+j*4+1] = small_b[i*NRS_2*4+j*4+2] = small_b[i*NRS_2*4+j*4+3] = 0;
	fclose(fp);
}

GraphicsWidget::~GraphicsWidget()
{
	delete[] small_b;
	delete[] large_b;
}

void GraphicsWidget::shaderLog(QString log)
{
	QFile file("shader.log");
	file.open(QIODevice::Append);
	
	QTextStream out(&file);
	out << log << endl;
	
	file.close();
}

void GraphicsWidget::outText(double x, double y, QString txt)
{
	int i;
	
	glDisable(GL_LIGHTING);
	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tex_font);
	glColor4d(1.0, 1.0, 1.0, 1.0);
	
	for(i=0; i<txt.length(); i++)
	{
		unsigned char c = txt.toAscii().data()[i];
		double ty = 15.0/16.0 - floor((double)c/16.0)/16.0;
		double tx = ((double)c - floor((double)c/16.0)*16.0)/16.0;
		
		glBegin(GL_QUADS);
		glMultiTexCoord2d(GL_TEXTURE3, tx, ty);
		glVertex2d(x,y);
		glMultiTexCoord2d(GL_TEXTURE3, tx + 0.8/16.0, ty);
		glVertex2d(x+0.08, y);
		glMultiTexCoord2d(GL_TEXTURE3, tx + 0.8/16.0, ty + 1.0/16.0);
		glVertex2d(x+0.08, y+0.1);
		glMultiTexCoord2d(GL_TEXTURE3, tx, ty + 1.0/16.0);
		glVertex2d(x, y+0.1);
		glEnd();
		
		x += 0.08;
	}
}

void GraphicsWidget::renderGPU(State s)
{
	shader -> bind();
	
	QVector4D basis_f[4];
	int i;
	for(i=0; i<4; i++)
	{
		basis_f[i] = QVector4D(s.basis[i][0], s.basis[i][1], s.basis[i][2], s.basis[i][3]);
	}
	shader -> setUniformValueArray("basis",basis_f,4);
	vector4 pos = s.pos.toVector4();
	shader -> setUniformValue("ship_pos", (GLfloat)pos[0], (GLfloat)pos[1], (GLfloat)pos[2], (GLfloat)pos[3]);
	shader -> setUniformValue("mass", (GLfloat)sim->getMass());
	shader -> setUniformValue("fov", (float)fov);
	shader -> setUniformValue("doppler", doppler ? 1 : 0);
	
	shader -> setUniformValue("background", 0);
	shader -> setUniformValue("small_b", 1);
	shader -> setUniformValue("large_b", 2);
	
	glBegin(GL_QUADS);
	glVertex2d(-0.99*width()/height(),0.99);
	glVertex2d(0.99*width()/height(),0.99);
	glVertex2d(0.99*width()/height(),-0.99);
	glVertex2d(-0.99*width()/height(),-0.99);
	glEnd();
	
	shader -> release();
}

//------------------------------------------------------------------

#define N_CIRC 16

void GraphicsWidget::renderGauges(State s)
{
	//calculate velocity 
	Metric* g = s.m -> getMetric(s.pos.getCoordSystem());

	vector4 dt = vector4(1.0, 0.0, 0.0, 0.0);
	dt /= sqrt(g -> g(dt, dt, s.pos));
	
	double r = s.pos[1];
	double M = sim -> getMass();
	vector4 dr = vector4(r/(r-2*M), 1.0, 0.0, 0.0);
	dr /= sqrt(-g -> g(dr, dr, s.pos));
	
	double gamma = g -> g(s.basis[0], dt, s.pos);
	double v = sqrt(1.0 - 1.0/(gamma*gamma));
	
	double vx, vy, vz;
	vx = g -> g(s.basis[1], dt, s.pos)/gamma;
	vy = g -> g(s.basis[2], dt, s.pos)/gamma;
	vz = g -> g(s.basis[3], dt, s.pos)/gamma;
	
	double rx, ry, rz;
	rx = g -> g(s.basis[1], dr, s.pos);
	ry = g -> g(s.basis[2], dr, s.pos);
	rz = g -> g(s.basis[3], dr, s.pos);
	
	glPushMatrix();
	glTranslatef(0.65, -0.85, 0.0);
	glScalef(0.5, 0.5, 0.5);
	outText(0.0, 0.0, "v = " + QString::number(v) + " c");
	glPopMatrix();
		
	glPushMatrix();
	glTranslatef(0.65, -0.95, 0.0);
	glScalef(0.5, 0.5, 0.5);
	outText(0.0, 0.0, "R = " + QString::number(r/M/2.0) + " Rs");
	glPopMatrix();
	
	glPushMatrix();
	glTranslatef(-1.0, 0.87, 0.0);
	glScalef(0.5, 0.5, 0.5);
	outText(0.0, 0.0, "Time warp: " + QString::number(s.time_warp) + "x");
	glPopMatrix();
	
	if(paused)
	{	
		glPushMatrix();
		glTranslatef(-0.18, -0.04, 0.0);
		glScalef(0.8, 0.8, 0.8);
		outText(0.0, 0.0, "PAUSED");
		glPopMatrix();
	}
	
	//render HUD
	
	glColor4d(0.0, 1.0, 0.5, 1.0);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(2.0);
	
	//"crosshair"
	
	glBegin(GL_LINES);
		glVertex2d(-0.1, 0.0);
		glVertex2d(-0.04, 0.0);
		
		glVertex2d(0.1, 0.0);
		glVertex2d(0.04, 0.0);
		
		glVertex2d(-0.04, -0.04);
		glVertex2d(0.0, 0.0);
		
		glVertex2d(0.0, 0.0);
		glVertex2d(0.04, -0.04);
	glEnd();
	
	//velocity
	
	int i;
	
	double vzs = vz/vx/tan(fov*M_PI/360.0);
	double vys = -vy/vx/tan(fov*M_PI/360.0);
	
	if(r > 2*M)
	{	
		if(vx > 0)
		{
			glBegin(GL_LINE_STRIP);
				for(i = 0; i < N_CIRC+1; i++)
					glVertex2d(vys + 0.1*cos(i*2*M_PI/N_CIRC), vzs + 0.1*sin(i*2*M_PI/N_CIRC));
			glEnd();
		}
		glBegin(GL_LINES);
			glVertex2d(vys - 0.1, vzs);
			glVertex2d(vys + 0.1, vzs);
			glVertex2d(vys, vzs - 0.1);
			glVertex2d(vys, vzs + 0.1);
		glEnd();
	}
	
	//direction to black hole
	glColor4d(0.0, 0.5, 1.0, 1.0);
	glLineWidth(3.0);
	
	double rzs = rz/rx/tan(fov*M_PI/360.0);
	double rys = -ry/rx/tan(fov*M_PI/360.0);
	
	if(r > 2*M && rx > 0)
	{
		glBegin(GL_LINES);
		glVertex2d(rys - 0.08, rzs - 0.08);
		glVertex2d(rys + 0.08, rzs + 0.08);
		glVertex2d(rys + 0.08, rzs - 0.08);
		glVertex2d(rys - 0.08, rzs + 0.08);
		glEnd();
	}
}

//------------------------------------------------------------------

void GraphicsWidget::initializeGL()
{
	makeCurrent();
	glDisable(GL_LIGHTING);
	glClearColor(1.0,1.0,1.0,1.0);
	
	shader = new QGLShaderProgram;
	
	if(!shader -> addShaderFromSourceFile(QGLShader::Vertex, "data/vshader.glsl"))
		shaderLog(shader->log());
	
	if(!shader -> addShaderFromSourceFile(QGLShader::Fragment, "data/pshader.glsl"))
		shaderLog(shader->log());
	
	if(!shader -> link())
		shaderLog(shader->log());
	
	glEnable(GL_ALPHA);
	
	glActiveTexture(GL_TEXTURE0);
	tex_sky = bindTexture(QImage("data/sky.png"));
	
	glActiveTexture(GL_TEXTURE1);
	glGenTextures(1, &tex_smallb);
	glBindTexture(GL_TEXTURE_2D, tex_smallb);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NRS_2, NBS_2, 0, GL_RGBA, GL_UNSIGNED_BYTE, small_b);
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	
	glActiveTexture(GL_TEXTURE2);
	glGenTextures(1, &tex_largeb);
	glBindTexture(GL_TEXTURE_2D, tex_largeb);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NRL_2, NBL_2, 0, GL_RGBA, GL_UNSIGNED_BYTE, large_b);
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	
	glActiveTexture(GL_TEXTURE3);
	tex_font = bindTexture(QImage("data/font.png"));
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
}

void GraphicsWidget::resizeGL(int w, int h)
{
	if(w==0) w=1;
	if(h==0) h=1;
	glViewport(0,0,w,h);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0*width()/height(), 1.0*width()/height(), -1.0, 1.0, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GraphicsWidget::paintGL()
{	
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glClear(GL_COLOR_BUFFER_BIT);
	
	State s;
	sim -> getState(s);
	
	//draw frame
	glBegin(GL_LINE_LOOP);
	glColor3d(0.0,0.0,0.0);
	glVertex2d(-0.99*width()/height(),0.99);
	glVertex2d(0.99*width()/height(),0.99);
	glVertex2d(0.99*width()/height(),-0.99);
	glVertex2d(-0.99*width()/height(),-0.99);
	glEnd();
	
	renderGPU(s);
	renderGauges(s);
	
	glFlush();
}

void GraphicsWidget::toggleDoppler()
{
	doppler = !doppler;
}

void GraphicsWidget::increaseFov()
{
	if(fov < 120.0) fov += 1.0;
}

void GraphicsWidget::decreaseFov()
{
	if(fov > 10.0) fov -= 1.0;
}

void GraphicsWidget::setPaused(bool p)
{
	paused = p;
}

//************************************************************************

MainWindow::MainWindow()
{
	sim = new Simulation();
	gl = new GraphicsWidget(this, sim);
	gl->move(0,0);
	
	resize(800,600);
	move(100,100);
	show();
	gl->show();
	
	gl->updateGL();
	
	time = getTime();
	QTimer::singleShot(0, this, SLOT(tick()));
	running = false;
}

MainWindow::~MainWindow()
{
}

long MainWindow::getTime()
{
	QTime t = QTime::currentTime();
	return t.msec() + 1000*t.second() + 1000*60*t.minute() + 1000*3600*t.hour();
}

long MainWindow::deltaT(long t1, long t2)
{
	return (t2-t1)%(24*3600*1000);
}

void MainWindow::keyPressEvent(QKeyEvent* ev)
{
	keys[ev->key()] = 1;
	keys_changed[ev->key()] = 1;
}

void MainWindow::keyReleaseEvent(QKeyEvent* ev)
{
	keys[ev->key()] = 0;
	keys_changed[ev->key()] = 1;
}

void MainWindow::resizeEvent(QResizeEvent*)
{
	gl->resize(width(), height());
}

void MainWindow::tick()
{
	double dt = (double)deltaT(time, getTime())/1000.0;
	time = getTime();
	
	gl -> updateGL();
	if(running)
	{
		if(keys[Qt::Key_R] && keys_changed[Qt::Key_R]) sim -> decreaseTWarp();
		if(keys[Qt::Key_T] && keys_changed[Qt::Key_T]) sim -> increaseTWarp();
		
		double coeff = 1.0;
		State s;
		sim -> getState(s);
		
		if(keys[Qt::Key_Shift]) coeff *= 5.0;
		if(keys[Qt::Key_W]) sim -> applyAngVel(0.0, coeff*0.6/s.time_warp, 0.0);
		if(keys[Qt::Key_S]) sim -> applyAngVel(0.0, coeff*-0.6/s.time_warp, 0.0);
		if(keys[Qt::Key_A]) sim -> applyAngVel(0.0, 0.0, coeff*0.6/s.time_warp);
		if(keys[Qt::Key_D]) sim -> applyAngVel(0.0, 0.0, coeff*-0.6/s.time_warp);
		if(keys[Qt::Key_Q]) sim -> applyAngVel(coeff*-0.6/s.time_warp, 0.0, 0.0);
		if(keys[Qt::Key_E]) sim -> applyAngVel(coeff*0.6/s.time_warp, 0.0, 0.0);
		
		if(keys[Qt::Key_Y]) sim -> applyForce(coeff*0.3, 0.0, 0.0);
		if(keys[Qt::Key_H]) sim -> applyForce(coeff*-0.3, 0.0, 0.0);
		if(keys[Qt::Key_G]) sim -> applyForce(0.0, coeff*0.3, 0.0);
		if(keys[Qt::Key_J]) sim -> applyForce(0.0, coeff*-0.3, 0.0);
		
		if(keys[Qt::Key_Z]) gl -> decreaseFov();
		if(keys[Qt::Key_X]) gl -> increaseFov();
		if(keys[Qt::Key_P] && keys_changed[Qt::Key_P]) gl -> toggleDoppler();
	}
	else
	{
		double coeff = 1.0;
		if(keys[Qt::Key_Shift]) coeff *= 5.0;
		if(keys[Qt::Key_W]) sim -> rotate(coeff*-0.6*dt, 0.0, 0.0);
		if(keys[Qt::Key_S]) sim -> rotate(coeff*0.6*dt, 0.0, 0.0);
		if(keys[Qt::Key_A]) sim -> rotate(0.0, coeff*-0.6*dt, 0.0);
		if(keys[Qt::Key_D]) sim -> rotate(0.0, coeff*0.6*dt, 0.0);
		if(keys[Qt::Key_Q]) sim -> rotate(0.0, 0.0, coeff*-0.6*dt);
		if(keys[Qt::Key_E]) sim -> rotate(0.0, 0.0, coeff*0.6*dt);
	}
	
	if(keys[Qt::Key_Space] && keys_changed[Qt::Key_Space])
	{
		running = !running;
		if(running)
		{
			sim -> launch();
			gl -> setPaused(false);
		}
		else
		{
			sim -> stop();
			gl -> setPaused(true);
		}
	}
	
	//set all "changed" keys to false
	QMap<int, int>::iterator it;
	for(it = keys_changed.begin(); it != keys_changed.end(); it++)
		keys_changed[it.key()] = 0;
	
	QTimer::singleShot(0, this, SLOT(tick()));
}
