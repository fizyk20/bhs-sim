#include <engine/particle.h>
#include <engine/schw.h>
#include <engine/dpintegrator.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>

using namespace std;

#define N_LB 1000
#define N_SB 103

vector< vector<double> > pre_tab;	//phi = pretab[n_b][rho];
double* tab_lb[N_LB]; 	//final table for large b
double* tab_sb[N_SB];
SchwManifold* m;

double b_n(int n)
{
	return 0.1*sqrt(2701.0 + (double)n*n);	//no sense for b=sqrt(27)
}

int n_b(double b)
{
	return (int)sqrt(b*b*100.0-2700.0);
}

double r_n(int n)
{
	return 0.05*sqrt(400.0 + (double)n*n);
}

int n_r(double r)
{
	return (int)sqrt(r*r*400.0-400.0);
}

//large b

#define DPHIMAX 0.0001
#define N 4

double calculate_phi_lb(double b, double r, map<double, double> data)
{
	double x[N], y[N];
	double phi = acos(sqrt(27.0)/b);
	double r_min = b/sqrt(3.0)*(cos(phi/3) + sqrt(3.0)*sin(phi/3));
	
	map<double, double>::iterator it = data.begin();
	
	while(it->first < r && it != data.end()) 
		it++;
	
	if(it == data.end())
		return asin(b/r/r_min);
	
	it--;
	if(it != data.begin()) it--;
	if(it != data.begin()) it--;
	
	map<double, double>::iterator it2 = it;
	int n_to_end = 0;
	while(it2 != data.end() && N-n_to_end > 0)
	{
		it2++;
		n_to_end++;
	}
	
	for(;N - n_to_end > 0; n_to_end++)	//go back some steps
		it--;
	
	int i,j;
	double result = 0.0;
	
	for(i=0; i<N; i++)
	{
		x[i] = it->first;
		y[i] = it->second;
		it++;
		if(it == data.end()) break;
	}
	
	for(i=0; i<N; i++)
	{
		double res = y[i];
		for(j=0; j<N; j++)
			if(j != i) res *= (r-x[j])/(x[i]-x[j]);
		result += res;
	}
	
	return result;
}

vector<double> generate_single_b_lb(int n)
{
	vector<double> result;
	map<double, double> path;
	
	double b = b_n(n);
	
	double phi = acos(sqrt(27.0)/b);
	double r_min = b/sqrt(3.0)*(cos(phi/3) + sqrt(3.0)*sin(phi/3));
	
	double prev_r = r_min;
	double r = r_min;
	
	Particle p(m, Point(EF, 0.0, r_min, M_PI/2, 0.0), vector4(-1.0/(1.0-2.0/r_min), 0.0, 0.0, b/r_min/r_min));
	DPIntegrator integrator;
	p.setIntegrator(&integrator);

	do
	{
		prev_r = r;
		
		p.propagate();
		Point pos = p.getPos();
		r = pos[1];
		phi = pos[3];
		path[prev_r/r_min] = pos[3];
		
		pos[3] = 0.0;
		p.setPosVel(pos, p.getVel());
	}
	while(phi > DPHIMAX);
	
	//sumowanie
	map<double,double>::iterator it, it2;
	
	for(it = path.begin(); it != path.end(); it++)
	{
		double sum = 0.0;
		for(it2 = --path.end(); it2 != it; it2--)
			sum += it2->second;
		
		it->second += sum + asin(b/r_n(n_r(prev_r)));
	}
	
	//interpolacja
	
	int i;
	result.push_back(path.begin()->second); //push r=r_min without recalculating
	
	for(i = 1; r_n(i)*r_min < prev_r; i++)
	{
		double a = calculate_phi_lb(b, r_n(i), path);
		//cout << a << endl; //debug
		result.push_back(a);
	}
	
	return result;
}

//small b

double b_n_s(int n)
{
	return n*0.05;
}

int n_b_s(double b)
{
	return (int)(b/0.05);
}

double r_n_s(int n)
{
	return 0.1*(n+1);
}

int n_r_s(double r)
{
	return (int)(r/0.1)-1;
}

double calculate_phi_sb(double b, double r, map<double, double> data)
{
	double x[4], y[4];
	
	map<double, double>::iterator it = data.begin();
	
	while(it->first < r && it != data.end()) 
		it++;
	
	if(it == data.end())
		return asin(b/r);
	
	it--;
	if(it != data.begin()) it--;
	
	int i,j,k;
	double result = 0.0;
	
	for(i=0; i<4; i++)
	{
		x[i] = it->first;
		y[i] = it->second;
		it++;
		if(it == data.end()) break;
	}
	k = i;
	
	for(i=0; i<k; i++)
	{
		double res = y[i];
		for(j=0; j<k; j++)
			if(j != i) res *= (r-x[j])/(x[i]-x[j]);
		result += res;
	}
	
	return result;
}

vector<double> generate_single_b_sb(int n)
{
	vector<double> result;
	map<double, double> path;
	double phi;
	
	double b = b_n_s(n);
	
	double prev_r = r_n_s(0);
	double r = prev_r;
	double szatan = sqrt(r*r-b*b*(1.0-2.0/r));
	
	Particle p(m, Point(EF, 0.0, r, M_PI/2, 0.0), vector4((-r+szatan)/(1.0-2.0/r), szatan, 0.0, b/r));
	DPIntegrator integrator;
	p.setIntegrator(&integrator);
	
	do
	{
		prev_r = r;
		
		p.propagate();
		Point pos = p.getPos();
		r = pos[1];
		phi = pos[3];
		path[prev_r] = pos[3];
		
		pos[3] = 0.0;
		p.setPosVel(pos, p.getVel());
		
	}
	while(phi > DPHIMAX || r < 10.0);
	
	Metric* g = m->getMetric(EF);
	vector4 u = p.getVel();
	cout << g->g(u,u,p.getPos()) << endl;
	
	//sumowanie
	map<double,double>::iterator it, it2;
	
	for(it = path.begin(); it != path.end(); it++)
	{
		double sum = 0.0;
		for(it2 = --path.end(); it2 != it; it2--)
			sum += it2->second;
		
		it->second += sum + asin(b/r_n_s(n_r_s(prev_r)));
	}
	
	//interpolacja
	
	int i;
	result.push_back(path.begin()->second); //push r=r_min without recalculating
	
	for(i = 1; r_n_s(i) < prev_r; i++)
	{
		double a = calculate_phi_sb(b, r_n_s(i), path);
		//cout << a << endl; //debug
		result.push_back(a);
	}
	
	return result;
}

int main()
{
	m = new SchwManifold(1.0);
	
	//large b
	//----------------------------------
	
	cout << "Starting - large b" << endl;
	
	int i,j,k;
	for(i=0; i<N_LB; i++)
	{
		cout << "b_" << i << " = " << b_n(i) << endl;
		vector<double> data = generate_single_b_lb(i);
		pre_tab.push_back(data);
		//debug
		cout << "Total deflection: " << data[0]-M_PI/2 << endl;
	}
	
	int i_max, max_size = 0;
	for(i=0; i<N_LB; i++)
	{
		if(pre_tab[i].size() > max_size)
		{
			i_max = i;
			max_size = pre_tab[i].size();
		}
	}
	cout << "i_max = " << i_max << ";  max_size = " << max_size << endl;
	max_size++;
	
	for(i=0; i<N_LB; i++)
	{
		tab_lb[i] = new double[max_size];
		for(j=0; j<max_size; j++)
		{
			if(j<pre_tab[i].size())
				tab_lb[i][j] = pre_tab[i][j];
			else
			{
				double b = b_n(i);
				double phi = acos(sqrt(27.0)/b);
				double r_min = b/sqrt(3.0)*(cos(phi/3) + sqrt(3.0)*sin(phi/3));
				
				tab_lb[i][j] = asin(b/r_min/r_n(j));
			}
		}
	}
	
	//zapis do pliku
	FILE* fp = fopen("large_b.dat","wb");
	int a;
	a = N_LB;
	fwrite(&a, sizeof(int), 1, fp);
	fwrite(&max_size, sizeof(int), 1, fp);
	for(i=0; i<N_LB; i++)
		fwrite(tab_lb[i], sizeof(double), max_size, fp);
	fclose(fp);
	
	//---------------------------
	//koniec large_b
	
	for(i=0; i<N_SB; i++)
	{
		cout << "b_" << i << " = " << b_n_s(i) << endl;
		vector<double> data = generate_single_b_sb(i);
		pre_tab.push_back(data);
		//debug
		cout << "Total deflection: " << data[0] << endl;
	}
	
	max_size = 0;
	for(i=0; i<N_SB; i++)
	{
		if(pre_tab[i].size() > max_size)
		{
			i_max = i;
			max_size = pre_tab[i].size();
		}
	}
	cout << "i_max = " << i_max << ";  max_size = " << max_size << endl;
	max_size++;
	
	for(i=0; i<N_SB; i++)
	{
		tab_sb[i] = new double[max_size];
		for(j=0; j<max_size; j++)
		{
			if(j<pre_tab[i].size())
				tab_sb[i][j] = pre_tab[i][j];
			else
			{
				double b = b_n_s(i);
				tab_sb[i][j] = asin(b/r_n_s(j));
			}
		}
	}
	
	//zapis do pliku
	fp = fopen("small_b.dat","wb");

	a = N_SB;
	fwrite(&a, sizeof(int), 1, fp);
	fwrite(&max_size, sizeof(int), 1, fp);
	for(i=0; i<N_SB; i++)
		fwrite(tab_sb[i], sizeof(double), max_size, fp);
	fclose(fp);
	
	return 0;
}
