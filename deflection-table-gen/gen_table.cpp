#include <particle.h>
#include <schw.h>
#include <rk4integrator.h>
#include <iostream>
#include <iomanip>
#include <fstream>
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

double b(double r, vector4 v)
{
    return r*r*r*v[3]/(r-2*M)/v[0]/M;
}

double deflected_final_phi(double r, double phi)
{
    vector4 v0 = ray(M, r, phi);
    cout << setprecision(10);
    cout << "v0 = [" << v0[0] << ", " << v0[1] << ", " << v0[2] << ", " << v0[3] << "] " << b(r, v0) << endl;
	Particle photon(&manifold, Point(EF, 0.0, r, M_PI/2, M_PI), v0);
	RK4Integrator integrator(1e-2);
	//integrator.setStepSize(1e-4);
	//integrator.setMaxErr(1e-6);
	photon.setIntegrator(&integrator);

	Metric* g = manifold.getMetric(EF);
	while(photon.getPos()[1] < 1000.0*M && photon.getPos()[1] > HORIZON)
	{
        double r = photon.getPos()[1];
        double E = (1-2*M/r)*photon.getVel()[0];
        double L = r*r*photon.getVel()[3];

		cout << "\tr = " << r << "\tphi = " << photon.getPos()[3] << "\tlen = " << g->g(photon.getVel(), photon.getVel(), photon.getPos());
        cout << "\tb = " << b(r, photon.getVel()) << "\tE = " << E;
        cout << "\tL = " << L;
		photon.propagate();

        // E and L correction
        r = photon.getPos()[1];
        vector4 v = photon.getVel();
        cout << "\tv_b = [" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]";
        v[0] = E/(1-2*M/r);
        v[3] = L/r/r;
        v[1] = E/2 - L*L*(1-2*M/r)/2/E/r/r;
        v[2] = 0.0;
        cout << "\tv = [" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]" << endl;
        photon.setVel(v);
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
	/*ofstream fout("deflection.csv");
	cout << "Start" << endl;
	fout << "r,phi,final,flat,deflection" << endl;
	for(double r = 1.6; r <= 100.0; r += 0.4)
	{
		for(int i = 0; i < N_POINTS + 1; i++)
		{
			double phi = M_PI * i / N_POINTS;
			cout << "R = " << r << "\tPhi = " << phi << " (" << i+1 << "/" << N_POINTS << ")" << endl;
			double result = deflected_final_phi(r, phi);
			double flat = flat_final_phi(phi);
			fout << r << "," << phi << "," << result << "," << flat << "," << flat - result << endl;
		}
	}*/
    cout << deflected_final_phi(1.6, 0.82) << endl;
	return 0;
}