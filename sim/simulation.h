#ifndef __SIMULATION__
#define __SIMULATION__

#include <entity.h>
#include <schw.h>
#include <rk4integrator.h>
#include <QtGui>

struct State
{
	SchwManifold* m;
	Point pos;
	vector4 basis[4];
	double proper_time;
	double time_warp;
};

class Simulation : public QThread
{
Q_OBJECT
	SchwManifold* m;
	Entity* ship;
	RK4Integrator integrator;
	
	//simulation parameters
	double time_warp;
	long time;
	bool running;
	
	long getTime();
	long deltaT(long, long);
public:
	Simulation();
	~Simulation();
	
	void run();
	
	double getMass() { return m->getMass(); }
	
public slots:
	void increaseTWarp();
	void decreaseTWarp();
	
	void applyForce(double, double, double);
	void applyAngVel(double, double, double);
	void rotate(double, double, double);
	
	void launch();
	void stop();
	void reset();
	void performTimeStep();
	void getState(State&);
	
signals:
	void autoStop();
};

#endif
