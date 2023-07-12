#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


double randn(double mu, double sigma)
{
	double U1, U2, W, multiply;
	static double X1, X2;
	static int hunt = 0;

	if (hunt == 1)
		{
			hunt = !hunt;
			return (mu + sigma * (double) X2);
		}
	do
		{
			U1 = -1 + ((double) rand()/RAND_MAX) * 2;
			U2 = -1 + ((double) rand()/RAND_MAX) *2;
			W = pow(U1,2)+pow(U2,2);
		}
	while (W >= 1 || W == 0);

	multiply = sqrt((-2 *log(W))/W);
	X1 = U1 * multiply;
	X2 = U2 * multiply;

	hunt = !hunt;

	return(mu + sigma * (double) X1);
}

/*
int main(int argc, char *argv[])
{
srand(time(NULL));
char *p;
double input = strtod(argv[1], &p);
double number = randn(0, input);
printf("%f \n", number);

return 0;
}
*/

