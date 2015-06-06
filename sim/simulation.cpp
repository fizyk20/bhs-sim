#include "simulation.h"
#include <QTime>
#include <QTimer>
#include <math.h>

Simulation::Simulation()
{
	m = new SchwManifold(1.0);
	double r0 = 30.0;
	Point p0(EF, r0 + 2*log(0.5*(r0-2.0)), r0, 1.5707, 0.0);	// so that t = 0 at the beginning
	ship = new Entity(m, p0, 
						 vector4(1.0, 0.0, 0.0, 0.0),
						 vector4(-r0/(r0-2.0), -1.0, 0.0, 0.0),
						 vector4(0.0, 0.0, 0.0, -1.0),
						 vector4(0.0, 0.0, -1.0, 0.0));
	ship -> setIntegrator(&integrator);
	running = false;
	time_warp = 1.0;
}

Simulation::~Simulation()
{
	running = false;
	if(ship != NULL) delete ship;
	delete m;
}

long Simulation::getTime()
{
	QTime t = QTime::currentTime();
	return t.msec() + 1000*t.second() + 1000*60*t.minute() + 1000*3600*t.hour();
}

long Simulation::deltaT(long t1, long t2)
{
	return (t2-t1)%(24*3600*1000);
}

void Simulation::run()
{
	exec();
}

void Simulation::launch()
{
	time = getTime();
	running = true;
	QTimer::singleShot(0, this, SLOT(performTimeStep()));
}

void Simulation::stop()
{
	running = false;
}

void Simulation::reset()
{
	running = false;
	if(ship != NULL) delete ship;
	double r0 = 30.0;
	Point p0(EF, 0.0, r0, 1.5707, 0.0);
	ship = new Entity(m, p0, 
						 vector4(1.0, 0.0, 0.0, 0.0),
						 vector4(-r0/(r0-2.0), -1.0, 0.0, 0.0),
						 vector4(0.0, 0.0, 0.0, -1.0),
						 vector4(0.0, 0.0, -1.0, 0.0));
	ship -> setIntegrator(&integrator);
}

void Simulation::increaseTWarp()
{
	time_warp *= 2.0;
}

void Simulation::decreaseTWarp()
{
	time_warp /= 2.0;
}

void Simulation::applyForce(double x, double y, double z)
{
	ship -> applyForce(x,y,z);
}

void Simulation::applyAngVel(double x, double y, double z)
{
	ship -> applyAngVel(x,y,z);
}

void Simulation::rotate(double pitch, double yaw, double roll)
{
	ship -> rotate(pitch, yaw, roll);
}

void Simulation::performTimeStep()
{
	if(!running) return;
	
	double dt = (double)deltaT(time,getTime())/1000.0*time_warp;
	if(dt > 1.0) dt = 1.0;
	time = getTime();
	
	ship -> propagate(dt);
	
	if(ship -> getPos()[1] < 0.01*m->getMass())
	{
		stop();
		emit autoStop();
		return;
	}
	
	QTimer::singleShot(0, this, SLOT(performTimeStep()));
}

void Simulation::getState(State& s)
{
	s.m = m;
	Point p = ship -> getPos();
	s.basis[0] = m->convertVectorTo(ship->getVel(), p, EF);
	s.basis[1] = m->convertVectorTo(ship->getX(), p, EF);
	s.basis[2] = m->convertVectorTo(ship->getY(), p, EF);
	s.basis[3] = m->convertVectorTo(ship->getZ(), p, EF);
	s.pos = m->convertPointTo(p, EF);
	s.time_warp = time_warp;
	s.proper_time = 0.0; //ship -> getProperTime();
}

