#include <particle.h>
#include <schw.h>
#include <dpintegrator.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>

using namespace std;

/*
 * This generator works as follows:
 * - set a value for the r coordinate of the observer
 * - generate light rays directed to the past for the circle in the horizontal plane
 * - propagate them and find the resulting deflection
 */
#define N_R 100			// number of points in the R coordinate
#define N_POINTS 180	// number of points on each circle
#define HORIZON 2.6*M

const double M = 0.5;
SchwManifold manifold(M);

vector4 ray(double M, double r, double phi)
{
	vector4 T(-1.0/sqrt(2), (M/r + 0.5)/sqrt(2), 0., 0.);
	vector4 A(-1.0/sqrt(2), (M/r - 1.5)/sqrt(2), 0., 0.);

	return T + A*cos(phi) + vector4(0., 0., 0., -1./r)*sin(phi);
}

double deflected_final_phi(double r, double phi)
{
	Particle photon(&manifold, Point(EF, 0.0, r, M_PI/2, M_PI), ray(M, r, phi));
	DPIntegrator integrator;
	integrator.setStepSize(1e-4);
	integrator.setMaxErr(1e-12);
	photon.setIntegrator(&integrator);

	Metric* g = manifold.getMetric(EF);
	while(photon.getPos()[1] < 1000.0*M && photon.getPos()[1] > HORIZON)
	{
		//cout << "\tr = " << photon.getPos()[1] << "\tphi = " << photon.getPos()[3] << "\t" << g->g(photon.getVel(), photon.getVel(), photon.getPos()) << endl;
		photon.propagate();
	}

	if(photon.getPos()[1] <= HORIZON)
		return -1000.0;

	double final_phi = photon.getPos()[3];	// phi coordinate

	double factor = sqrt(1.0 - 2.0*M/photon.getPos()[1]);

	// approximate the rest of the path by a straight line - calculate the additional angle
	double vr = factor*photon.getVel()[0] + photon.getVel()[1]/factor;
	double vphi = photon.getPos()[1] * photon.getVel()[3];

	double add_phi = atan2(vphi, vr);

	return final_phi + add_phi;
}

double flat_final_phi(double phi)
{
	double cosphi = cos(phi);
	return acos((3*cosphi-1)/(3-cosphi));
}

int main()
{
	double r = 20.0;

	for(int i = 0; i < N_POINTS + 1; i++)
	{
		double phi = M_PI * i / N_POINTS;
		double result = deflected_final_phi(r, phi);
		cout << "Phi = " << phi << "\tFinal phi = " << result << "\tFlat final phi = " << flat_final_phi(phi) << endl;
	}
	return 0;
}