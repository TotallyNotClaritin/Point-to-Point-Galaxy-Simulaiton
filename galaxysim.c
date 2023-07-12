#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"
#include <string.h>
#include <math.h>
#include <omp.h>

#ifdef _WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif


/*Simulation is in units of parsecs, solar masses, and years*/

int num_teams = 32;
int num_blocks = 64;
int Norder = 6;
int threads = 8;
int softener = 100;
double accelcap = 4 * pow(10,-6.0170); /*pow(10,-6.0125);*/
double bigG = 4.3009 * pow(10,-3);

void leapfrog(double y[], double dydx[], int n, double x, double h, double yout[], double r[], double m[], int part, void (*derivs)(double, double [], double [], double[], double[], int));
double randn(double  mu, double sigma);

/*

void collide(double x[], double m[], int j, int k)
{
double ang[2];
double normal[3];
double vimpact;
double vone[3];
double vtwo[3];
double mass;
double sqdist;
double impulse;

int count = 0;

for(int i = 2; i <= 6; i = i + 2)
{
vone[count] = x[i + 6*j];
vtwo[count] = x[i + 6*k];
count = count + 1;
}

count = 0;

*calculate normal
for(int i = 1; i <= 5; i = i + 2)
{
        normal[count] = x[i + 6*k] - x[i + 6*j];
        sqdist += pow(normal[count], 2);
        count = count + 1;
}

normalize
for(int i = 0; i <= 2; i++)
{
        normal[i] /= sqrt(sqdist);
}

calculate reduced mass of the system

mass = 1 / (1/m[j+1] + 1/m[k+1]);

calculate the velocity of impact

for(int i = 0; i <= 2; i++)
{
       vimpact += normal[i] * (vone[i] - vtwo[i]);
}
printf("%f" , vimpact);
impulse = 2 * mass * vimpact;

double dvone[3];
double dvtwo[3];

for(int i = 0; i <= 2; i++)
{
dvone[i] = -impulse/m[j+1] * normal[i];
dvtwo[i] = impulse/m[k+1] * normal[i];
}

add the delta v
count = 0;

for(int i = 1; i <= 5; i = i +2)
{
x[i + 6*j] += dvone[count];
x[i + 6*k] += dvtwo[count];

count = count + 1;
}
printf("\n");


}
void collision(double x[], double r[], double m[], int part)
{

collision check if the sq distance between both is less than or qeual to the
 * sum of the radii squared they will be in collision

double sqdist = 0.0;
double radiisq = 0.0;

        for(int j = 0; j < part; j++)
        {
                for(int k = 0; k < part; k++)
                {
                        if(j != k)
                        {
                                for(int i = 1; i <= 5; i = i + 2)
                                {
                                       sqdist += pow(x[i + 6*j] - x[i + 6*k],2);
                                }
                                radiisq = pow(r[j+1]+r[k+1],2);

                                if(sqdist <= radiisq)
                                        {

                                                collide(x, m, j, k);

                                        }
                        }
                }
        }
}

*/
void derivs(double t, double x[], double dxdt[], double r[], double m[], int part)
{

double accel[3];

#pragma omp parallel for  num_threads(threads)
for(int i = 0; i < part; i++)
{
        double denom = 0.0;

        dxdt[1 + 6*i] = x[2 + 6*i];
        dxdt[3 + 6*i] = x[4 + 6*i];
        dxdt[5 + 6*i] = x[6 + 6*i];

        dxdt[2 + 6*i] = 0;
        dxdt[4 + 6*i] = 0;
        dxdt[6 + 6*i] = 0;

	accel[0] = 0;
	accel[1] = 0;
	accel[2] = 0;

       for(int j = 0; j < part; j++)
        {
        if (j != i)
                {
                denom = pow(sqrt(pow(x[1+6*i]-x[1+6*j],2)+pow(x[3+6*i]-x[3+6*j],2)+pow(x[5+6*i]-x[5+6*j],2)),3) + softener;
		
                accel[0] = (bigG * m[i+1]*m[j+1]*(x[1 + 6*j]-x[1 + 6*i]))/denom;
                accel[1] = (bigG * m[i+1]*m[j+1]*(x[3 + 6*j]-x[3 + 6*i]))/denom;
                accel[2] = (bigG * m[i+1]*m[j+1]*(x[5 + 6*j]-x[5 + 6*i]))/denom;


	       for(int z = 1; z <= 3; z++)
		        {
	      		 if(accel[z-1] > accelcap)
	       		 {
		        accel[z-1] = accelcap;
		        }
		       else if(accel[z-1] < -accelcap)
		        {
		        accel[z-1] = -accelcap;
		        }
		       dxdt[z*2+6*i] += accel[z-1];
			}

                }
        }

}
}

double calc(double x[], double m[], int part)
{
double vsq;
double separation;
double vr;
double energy;
double kinetic;
double gravitational;
double d;

#pragma omp parallel for reduction(+:kinetic) num_threads(threads)
for(int i =0; i < part; i++)
{
	kinetic += 0.5*m[i+1]*(x[2+6*i]*x[2+6*i]+x[4+6*i]*x[4+6*i]+x[6+6*i]+x[6+6*i]);
}


#pragma omp parallel for reduction(-:gravitational) collapse(2) num_threads(threads)
for(int i = 0; i < part; i++)
	{
		for(int j = 0; j < part; j++)
		{
			if(i != j)
			{
				d = pow(x[1+6*i]-x[1+6*j],2)+pow(x[3+6*i]-x[3+6*j],2)+pow(x[5+6*i]-x[5+6*j],2);
				if(d > pow(10,-12))
				{
				gravitational -= (bigG * m[i+1]*m[j+1])/sqrt(d);
				}
			}		
		}
	}

printf("%lf, %lf\n", kinetic, gravitational);
energy = kinetic + gravitational;
return energy;
}


int main(int argc, char *argv[])
{


int num = 1000;
int buf[num];    // initialize the buffer to 0s
#pragma acc parallel create(buf[num])
    for (int i = 0; i < num; ++i) {
        buf[i] = 1;
    }
    if (buf[0] == 1) {
        printf("Running on CPU\n");
    } else {
        printf("Running on GPU\n");
    }

int devices = omp_get_num_devices();
printf("%d\n", devices);






if(argc != 4)
{
printf("Usage: %s step size, number of particles, simulation time length \n", argv[0]);
return 1;
}

printf("%s step size = %s\n, number of particles = %s, simualtion time length = %s \n", argv[0],argv[1],argv[2],argv[3]);

/*frog goes {number of equations], {x initial conditions}, {v initial}, [t bottom], [t step], [t max], derivs(n,
 * x[], dvdt[]) */

struct timespec start, finish;
double elapsed;

clock_gettime(CLOCK_MONOTONIC, &start);

srand48(time(NULL));

char *input = argv[1];
char *inputtwo = argv[2];
char *inputthree = argv[3];
char *p;

double t;
double h = strtod(input, &p);
int part = strtol(inputtwo, &p, 10);
double t_max = strtod(inputthree, &p);
int t_steps = (int) t_max/h + 1;


int n = part*Norder;
double *m = dvector(1,part);
double *r = dvector(1,part);
double *x = dvector(1,n);
double *dxdt = dvector(1,n);

double min = -3066;
double max = -min;
double numin = -6132;
double numax = -numin;
double outmin = -32408;
double outmax = -outmin;
double nuzmin = -306;
double nuzmax = -nuzmin;
double outzmin = -270;
double outzmax = -outzmin;
double vmin = -0.00000000001;
double vmax = 0.000000000001;
double radialv;
double bulgemass;
double massinside;
double radius;
double secondradius;
double angle;
double ener;

/*random particle generation code for primordial composition (75% hydrogen, 25% helium)*/
for (int i = 0; i < part; i++)
{
x[1 + 6*i] = 100000;
x[3 + 6*i] = 100000;
x[5 + 6*i] = 100000;
}
/*place holder values*/
/*gal one*/

#pragma omp parallel for  num_threads(threads)
for(int z = 0; z <= 1; z++)
{

for(int i = 0 + z*(part/10)*5; i < part/10 + z*(part/10)*5; i++)
{
m[i+1] = 1;
r[i+1] = 7.35* pow(10,-8);
        while(sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2) + pow(x[5+6*i],2)) > max)
        {
        x[1 + 6*i] = min + ((double)drand48()) * (max - min);
        x[3 + 6*i] = min + ((double)drand48()) * (max - min);
        x[5 + 6*i] = min + ((double)drand48()) * (max - min);
        }
}


for(int i = part/10 + z*(part/10)*5; i < (part/10)*5 + z*(part/10)*5; i++)
{
m[i+1] = 10;
r[i+1] = 7.35*pow(10,-8);
        while(sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) < max || sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) > outmax)
        {
        x[1 + 6*i] = outmin + ((double)drand48()) * (outmax - outmin);
        x[3 + 6*i] = outmin + ((double)drand48()) * (outmax - outmin);
        x[5 + 6*i] = nuzmin + ((double)drand48()) * (nuzmax - nuzmin);
        }
}

/*
for(int i = part/10 + z*(part/10)*5; i < (part/10)*3 + z*(part/10)*5; i++)
{
m[i+1] = 1;
r[i+1] = 7.35*pow(10,-8);
        while(sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) < max || sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) > numax)
        {
        x[1 + 6*i] = numin + ((double)drand48()) * (numax - numin);
        x[3 + 6*i] = numin + ((double)drand48()) * (numax - numin);
        x[5 + 6*i] = nuzmin + ((double)drand48()) * (nuzmax - nuzmin);
        }
}

for(int i = (part/10)*3 + z*(part/10)*5; i < (part/10)*5 + z*(part/10)*5; i++)
{
m[i+1] = 1;
r[i+1] = 7.35*pow(10,-8);
        while(sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) < numax || sqrt(pow(x[1+6*i],2) + pow(x[3+6*i],2)) > outmax)
        {
        x[1 + 6*i] = outmin + ((double)drand48()) * (outmax - outmin);
        x[3 + 6*i] = outmin + ((double)drand48()) * (outmax - outmin);
        x[5 + 6*i] = outzmin + ((double)drand48()) * (outzmax - outzmin);
        }
}
*/
}
double sigma;

for(int z = 0; z <= 1; z++)
{

bulgemass = 0;


for(int i = 0 + z*(part/10)*5; i < part/10 + z*(part/10)*5; i++)
{
bulgemass += m[i+1];
}


for(int i = 0 + z*(part/10)*5; i < part/10 + z*(part/10)*5; i++)
{
radius = sqrt(pow(x[1 + 6*i],2) + pow(x[3 + 6*i],2) + pow(x[5 + 6*i],2));
massinside = 0;

for(int j = 0 + z*(part/10)*5; j < part/10 + z*(part/10)*5; j++)
{
    if(i != j)
    {
    secondradius = sqrt(pow(x[1 + 6*j],2) + pow(x[3 + 6*j],2) + pow(x[5 + 6*j],2));
        if(secondradius < radius)
        {
        massinside += m[j+1];
        }
    }
}
sigma = sqrt((bigG*massinside)/(3*radius));

x[2 + 6*i] = randn(0.0, sigma);
x[4 + 6*i] = randn(0.0, sigma);
x[6 + 6*i] = randn(0.0, sigma);
}

for(int i = part/10 + z*(part/10)*5;  i < (part/10)*5 + z*(part/10)*5; i++)
{
radius = sqrt(pow(x[1 + 6*i],2) + pow(x[3 + 6*i],2));
massinside = 0;

for(int j = part/10 + z*(part/10)*5;  j < (part/10)*5 + z*(part/10)*5; j++)
{
    if(i != j)
    {
    secondradius = sqrt(pow(x[1 + 6*j],2) + pow(x[3 + 6*j],2));
        if(secondradius < radius)
        {
        massinside += m[j+1];
        }
    }
}


radialv = sqrt((bigG*(bulgemass+massinside))/(radius));
angle = atan2(x[3+6*i],x[1+6*i]);

x[2 + 6*i] = -radialv * sin(angle);
x[4 + 6*i] = radialv * cos(angle);
x[6 + 6*i] = 0;
}
}


/*gal two*/


for(int i = part/2; i < part; i++)
{
x[1 + 6*i] += 90000;
x[2 + 6*i] -= 0.005;
}



/*galaxy gen complete*/
char identifier[100] = "-";
strcat(identifier, argv[1]);
strcat(identifier, argv[2]);
strcat(identifier, argv[3]);
strcat(identifier, ".dat");
char leapname[50] = "nbodyleap";

strcat(leapname,identifier);

remove(leapname);
FILE *frogdat = fopen(leapname, "a");


fprintf(frogdat, "%lf, ", t);
for(int j = 1; j<=n; j++)
{
        if(j == n)
        {
		ener = calc(x,m,part);
                fprintf(frogdat, "%lf, %lf", x[j], ener);
        }
        else
        {
                fprintf(frogdat, "%lf, ", x[j]);
        }
}
fprintf(frogdat, "\n");

for(int j=1; j<=t_steps; j++)

{
derivs(t,x,dxdt,r,m, part);
leapfrog(x, dxdt, n, t, h, x, r, m, part, *derivs);
//collision(x,r,m,part);
t = j *h;
fprintf(frogdat, "%lf, ", t);
for(int j = 1; j<=n; j++)
{
        if(j == n)
        {
		ener = calc(x,m,part);
                fprintf(frogdat, "%lf, %lf", x[j], ener);
        }
        else
        {
                fprintf(frogdat, "%lf, ", x[j]);
        }
}
fprintf(frogdat, "\n");
}

clock_gettime(CLOCK_MONOTONIC, &finish);
elapsed = (finish.tv_sec - start.tv_sec);
elapsed += (finish.tv_nsec - start.tv_nsec)/1000000000.0;

printf("The simulation was successful, it took %f seconds\n", elapsed);


return 0;
}

