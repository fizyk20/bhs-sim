#include <particle.h>
#include <schw.h>
#include <dpintegrator.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>
#include <stdlib.h>

using namespace std;

/*
 * This generator works as follows:
 * - set a value for the r coordinate of the observer
 * - generate light rays directed to the past for the circle in the horizontal plane
 * - propagate them and find the resulting deflection
 */
#define PSPHERE 3.0*M

const double M = 1.0;
SchwManifold manifold(M);

// function generating a parametrized four-velocity of a ray of light
vector4 ray(double M, double r, double phi)
{
    vector4 T(-1.0/sqrt(2), (M/r + 0.5)/sqrt(2), 0., 0.);
    vector4 A(-1.0/sqrt(2), (M/r - 1.5)/sqrt(2), 0., 0.);

    return T + A*cos(phi) + vector4(0., 0., 0., -1./r)*sin(phi);
}

// this returns the minimal phi for which the ray doesn't fall into the black hole
double threshold_phi(double r)
{
    double a = 27.0/8*M*M;
    double b = r + 2*M;
    double c = 3*r - 2*M;
    double rm = r/M;
    double delta = 4*M*M*M*r*r*r*r*r*(rm*rm*rm - 27*rm + 54);
    if(r < 3.0*M)
    {
        return acos((a*b*c - 0.5*sqrt(delta))/(a*b*b + r*r*r*r));
    }
    else
    {
        return acos((a*b*c + 0.5*sqrt(delta))/(a*b*b + r*r*r*r));
    }
}

// function mapping 0-1 interval onto the set of possible angles
double map_phi(double r, double x)
{
    if(x <= 0.0) x = 0.001;
    if(x > 1.0) x = 1.0;

    double min_phi = threshold_phi(r);
    double max_phi = M_PI;
    return (max_phi-min_phi)*x + min_phi;
}

// function calculating the final phi coordinate of a deflected ray
// the input parameters are the starting r coordinate and the angle from the direction to the black hole
double deflected_final_phi(double r, double phi)
{
    bool start_below_psphere = r <= PSPHERE;

    Particle photon(&manifold, Point(EF, 0.0, r, M_PI/2, M_PI), ray(M, r, phi));
    DPIntegrator integrator;
    integrator.setStepSize(1e-4);
    integrator.setMaxErr(1e-12);
    integrator.setMaxStep(3.0);
    photon.setIntegrator(&integrator);

    Metric* g = manifold.getMetric(EF);
    while((photon.getPos()[1] < 1000.0*M || photon.getPos()[1] < 30*r) &&   // generate until R = 1000 M or R = 30 R0, whichever is greater
         // if starting from below photon sphere, treat passing the starting radius as falling in; else, it is passing the photon sphere
          ((!start_below_psphere && photon.getPos()[1] > PSPHERE) || (start_below_psphere && photon.getPos()[1] >= r)))
    {
        //cout << "    r = " << photon.getPos()[1] << " phi = " << photon.getPos()[3] << endl;
        photon.propagate();
    }

    if(photon.getPos()[1] <= PSPHERE)
        return -1000.0;

    double final_phi = photon.getPos()[3];    // phi coordinate

    double factor = sqrt(1.0 - 2.0*M/photon.getPos()[1]);

    // approximate the rest of the path by a straight line - calculate the additional angle
    double vr = photon.getVel()[1]/factor;
    double vphi = photon.getPos()[1] * photon.getVel()[3];

    double add_phi = atan2(vphi, vr);

    return final_phi + add_phi;
}

// function returning the final phi coordinate of a ray if it was travelling in a flat space
double flat_final_phi(double phi)
{
    double cosphi = cos(phi);
    return acos((3*cosphi-1)/(3-cosphi));
}

double calc_b(double r, double x)
{
    double phi = map_phi(r, x);
    vector4 v = ray(M, r, phi);
    return r*r*v[3]/((1-2*M/r)*v[0] - v[1]);
}

int main(int argc, char** argv)
{
    for(double x = 0.0; x < 0.99; x += 0.05)
    {
        char filename[30];
        sprintf(filename, "def-x%3.2f.dat", x);

        ofstream fout(filename);
        for(double rx = -50.0; rx <= 70.0; rx += 1.0)
        {
            // generate the values of the radius to focus on the area below the photon sphere
            double r = 3*M*exp(rx/10);
            double phi = map_phi(r, x);
            double result = deflected_final_phi(r, phi);
            double flat = flat_final_phi(phi);
            cout << "R = " << r << "\tx = " << x << "\tPhi = " << phi << "\tResult = " << result << "\tFlat = " << flat << "\tDefl = " << flat - result << endl;
            fout << r << "\t" << flat - result << "\t" << calc_b(r,x) << endl;
        }
        fout.close();
    }
    return 0;
}
