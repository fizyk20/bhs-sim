varying vec4 position;

uniform sampler2D large_b;
uniform sampler2D small_b;
uniform sampler2D background;

uniform vec4 basis[4];
uniform vec4 ship_pos;
uniform float mass;
uniform float fov;
uniform int doppler;

#define N 2
#define N1 0
#define N2 2

#define NRL 476
#define NBL 1000
#define NRS 505
#define NBS 103

#define NRL_2 512.0
#define NBL_2 1024.0
#define NRS_2 512.0
#define NBS_2 128.0

float vec4tofloat(in vec4 b)
{
	vec4 a = b*255.0;
	float x = a[2]/256.0;
	x += a[1];
	x /= 256.0;
	x += a[0];
	x /= 256.0;
	x += 1.0;
	x *= pow(2.0, a[3]-128.0);
	
	return x;
}

float getSmallB(in float b, in float r)
{
	int nrs, nbs;
	nrs = NRS; //textureSize(small_b, 0)[0];
	nbs = NBS;  //textureSize(small_b, 0)[1];

	float result = 0.0;
	int bi = int(b/0.05);
	int ri = int(r/0.1)-1;
	
	float d[N];
	int i, j, k;
	
	if(ri > nrs-1) return asin(b/r);
	
	if(ri > nrs-N2) ri = nrs-N2;	//ekstrapolacja
	if(bi > nbs-N2) bi = nbs-N2;
	if(ri < N1) ri = N1;
	if(bi < N1) bi = N1;
	for(i = 0; i < N; i++)
	{
		d[i] = 0.0;
		for(j = ri-N1; j < ri+N2; j++)
		{
			float a = vec4tofloat(texture2D(small_b, vec2( (float(j)+0.5)/NRS_2, (float(i+bi-N1)+0.5)/NBS_2 ) ));

			for(k = ri-N1; k < ri+N2; k++)
				if(k != j) a *= (r - 0.1*float(k+1))/(0.1*(float(j-k)));

			d[i] += a;
		}
		
		for(k = 0; k < N; k++)
			if(k != i) d[i] *= (b - 0.05*(float(k+bi-N1)))/(0.05*(float(i-k)));
			
		result += d[i];
	}
	
	return(result);
}

float b_n(in int n) 
{ 
	float n_f = float(n);
	return(0.1*sqrt(2701.0 + n_f*n_f));
}

float r_n(in int n)
{
	float n_f = float(n);
	return(0.05*sqrt(400.0 + n_f*n_f));
}

float getLargeB(in float b, in float r2, in bool ur)
{
	int nrl, nbl;
	nrl = NRL;  //textureSize(large_b, 0)[0];
	nbl = NBL; //textureSize(large_b, 0)[1];

	float result = 0.0;
	float r_min;
	float phi;
	phi = acos(sqrt(27.0)/b);
	r_min = b/sqrt(3.0)*(cos(phi/3.0) + sqrt(3.0)*sin(phi/3.0));
	float r = r2/r_min;

	int bi = int(sqrt(b*b*100.0-2701.0));
	int ri = int(sqrt(r*r*400.0-400.0));
	
	float d[N];
	int i, j, k;
	
	if(ri > nrl-1) return(asin(b/r/r_min));
	
	if(ri > nrl-N2) ri = nrl-N2;	//ekstrapolacja
	if(bi > nbl-N2) bi = nbl-N2;
	if(ri < N1) ri = N1;
	if(bi < N1) bi = N1;

	for(i = 0; i < N; i++)
	{
		d[i] = 0.0;
		for(j = ri-N1; j < ri+N2; j++)
		{
			float a = vec4tofloat(texture2D(large_b, vec2( (float(j)+0.5)/NRL_2, (float(i+bi-N1)+0.5)/NBL_2 ) ));
			
			for(k = ri-N1; k < ri+N2; k++)
				if(k != j) a *= (r - r_n(k))/(r_n(j) - r_n(k));

			d[i] += a;
		}
		
		for(k = 0; k < N; k++)
			if(k != i) d[i] *= (b - b_n(k+bi-N1))/(b_n(i+bi-N1) - b_n(k+bi-N1));
			
		result += d[i];
	}
	
	if(ur)
	{
		float a = 0.0;
		for(i = 0; i < N; i++)
		{
			d[i] = vec4tofloat(texture2D(large_b, vec2( 0.5/NRL_2, (float(i+bi-N1)+0.5)/NBL_2 ) ));
			
			for(j = 0; j < N; j++)
				if(j != i) d[i] *= (b - b_n(j+bi-N1))/(b_n(i+bi-N1) - b_n(j+bi-N1));
			a += d[i];
		}
		result = 2.0*a - result;
	}
	
	return(result);
}

float deflection(in float b, in float r, in bool ur)
{
	if(b >= sqrt(27.0))
		return(getLargeB(b,r,ur));
	else
		return(getSmallB(b,r));
}

float g(in vec4 u, in vec4 v)
{
	float r = ship_pos[1];
	float th = ship_pos[2];
	return((1.0-2.0*mass/r)*u[0]*v[0] - u[0]*v[1] - u[1]*v[0] - r*r*u[2]*v[2] - r*r*sin(th)*sin(th)*u[3]*v[3]);
}

float tanh(in float x)
{
	return((exp(x)-exp(-x))/(exp(x)+exp(-x)));
}

float r_f(in float x)
{
	return exp(-(x-650.0)*(x-650.0)/28800.0) + 0.5*(tanh((400.0-x)/50.0)+1.0) + 0.2/(1.0+exp((660.0-x)/100.0));
}

float g_f(in float x)
{
	
	return exp(-(x-510.0)*(x-510.0)/9800.0) + 0.5*(tanh((340.0-x)/30.0)+1.0);
}

float b_f(in float x)
{
	return 0.5*(tanh((440.0-x)/60.0)+1.0);
}

vec4 doppler_shift(in vec4 color, in float factor)
{
	float R = color[0];
	float G = color[1];
	float B = color[2];
	
	float wavelengths[3];
	wavelengths[0] = 430.0;
	wavelengths[1] = 520.0;
	wavelengths[2] = 650.0;
	
	float r[3];
	r[0] = r_f(wavelengths[0]);
	r[1] = r_f(wavelengths[1]); 
	r[2] = r_f(wavelengths[2]);

	float g[3];
	g[0] = g_f(wavelengths[0]);
	g[1] = g_f(wavelengths[1]);
	g[2] = g_f(wavelengths[2]);
	
	float b[3];
	b[0] = b_f(wavelengths[0]);
	b[1] = b_f(wavelengths[1]);
	b[2] = b_f(wavelengths[2]);
	
	float det = r[0]*g[1]*b[2] + r[1]*g[2]*b[0] + r[2]*g[0]*b[1] - r[0]*g[2]*b[1] - r[1]*g[0]*b[2] - r[2]*g[1]*b[0];
	float det0 = R*g[1]*b[2] + r[1]*g[2]*B + r[2]*G*b[1] - R*g[2]*b[1] - r[1]*G*b[2] - r[2]*g[1]*B;
	float det1 = r[0]*G*b[2] + R*g[2]*b[0] + r[2]*g[0]*B - r[0]*g[2]*B - R*g[0]*b[2] - r[2]*G*b[0];
	float det2 = r[0]*g[1]*B + r[1]*G*b[0] + R*g[0]*b[1] - r[0]*G*b[1] - r[1]*g[0]*B - R*g[1]*b[0];
	
	R = det0/det*r_f(wavelengths[0]/factor) + det1/det*r_f(wavelengths[1]/factor) + det2/det*r_f(wavelengths[2]/factor);
	G = det0/det*g_f(wavelengths[0]/factor) + det1/det*g_f(wavelengths[1]/factor) + det2/det*g_f(wavelengths[2]/factor);
	B = det0/det*b_f(wavelengths[0]/factor) + det1/det*b_f(wavelengths[1]/factor) + det2/det*b_f(wavelengths[2]/factor);
	
	vec3 res = vec3(R,G,B);
	if(res[0] > 1.0) res /= res[0];
	if(res[1] > 1.0) res /= res[1];
	if(res[2] > 1.0) res /= res[2];
	
	return vec4(res[0], res[1], res[2], 1.0);
}

void main()
{
	float dmax = tan(fov*3.14159265/360.0);
	float dx = position[0]*dmax;
	float dy = position[1]*dmax;

	vec4 dir = (basis[1] - dx*basis[2] + dy*basis[3])/sqrt(dx*dx+dy*dy+1.0);
	dir -= basis[0]; //vector directed to the past
			
	float r = ship_pos[1];
	float th = ship_pos[2];
	float ph = ship_pos[3];
	float b = r*r*sqrt(sin(th)*sin(th)*dir[3]*dir[3] + dir[2]*dir[2])/((1.0-2.0*mass/r)*dir[0] - dir[1])/mass;
	b = abs(b);
	
	r /= mass;
			
	bool ur = (r > 2.0 && dir[1] < 0.0) || (r < 2.0 && (1.0-2.0/r)*dir[0] - dir[1] > 0.0);
			
	if((r <= 3.0 && b >= sqrt(27.0)) || (b < sqrt(27.0) && ur))
	{
		gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
	}
	else
	{
		float def_ang = deflection(b, r, ur);
				
		vec3 vr = vec3(sin(th)*cos(ph), sin(th)*sin(ph), cos(th));
		vec3 dth = vec3(cos(th)*cos(ph), cos(th)*sin(ph), -sin(th));
		vec3 dph = vec3(-sin(th)*sin(ph), sin(th)*cos(ph), 0.0);
		vec3 vdir = dir[2]*dth + dir[3]*dph;
		
		vec3 axis = normalize(cross(vr,vdir));
		vec3 vr2 = normalize(cross(axis,vr));
		
		vec3 result = vr*cos(def_ang) + vr2*sin(def_ang);
		float th2 = acos(result[2]);
		float ph2 = atan(result[1], result[0]);
		
		while(ph2 < 0.0) ph2 += 2.0*3.14159265;

		vec4 pixel = texture2D(background, vec2(ph2/2.0/3.14159265, th2/3.14159265));
				
		if(doppler != 0)
		{
			float doppler_factor = g(dir, basis[0])/((1.0-2.0/r)*dir[0] - dir[1]);		
			gl_FragColor = doppler_shift(pixel, doppler_factor);
		}
		else
		{
			gl_FragColor = pixel;
		}	
	}
} 
