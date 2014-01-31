#ifndef __INTERFACE__
#define __INTERFACE__

#include "simulation.h"
#include <QtGui>
#include <QMouseEvent>
#include <QtOpenGL>

class GraphicsWidget : public QGLWidget
{
Q_OBJECT
	GLuint tex_sky, tex_largeb, tex_smallb;
	GLuint tex_font;

	int nb, nr, nbs, nrs;
	unsigned char *small_b, *large_b;
	
	bool doppler;
	double fov;
	bool paused;
	
	QPoint lastPos;
	
	Simulation* sim;
	QGLShaderProgram* shader;
	
	vector4 floattovec4(double x);
	void shaderLog(QString);
	void outText(double x, double y, QString txt);

protected:
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	
	void renderGPU(State);
	void renderGauges(State);
	
public:
	GraphicsWidget(QWidget*, Simulation*, bool _doppler = true, double _fov = 60.0);
	~GraphicsWidget();
	
	void toggleDoppler();
	void increaseFov();
	void decreaseFov();
	void setPaused(bool);
};

class MainWindow : public QMainWindow
{
Q_OBJECT
	Simulation* sim;
	GraphicsWidget* gl;
	
	QMap<int,int> keys;
	QMap<int,int> keys_changed;
	
	long time;
	bool running;
	
	long getTime();
	long deltaT(long, long);
protected:
	void keyPressEvent(QKeyEvent* event);
	void keyReleaseEvent(QKeyEvent* event);
	void resizeEvent(QResizeEvent*);
public:
	MainWindow();
	~MainWindow();
public slots:
	void tick();
};

#endif
