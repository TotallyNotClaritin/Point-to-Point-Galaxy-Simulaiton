#define NRANSI

void leapfrog(double y[], double dydx[], int n, double x, double h, double yout[], double r[], double m[], int part,
        void (*derivs)(double, double [], double [], double[], double[], int))
{
        int i;
        double xh,hh;

        hh=h*0.5;
        xh=x+hh;
        for (i=1;i<=n;i=i+2) yout[i+1]=y[i+1]+hh*dydx[i+1];
        for (i=1;i<=n;i=i+2) yout[i]=y[i]+h*yout[i+1];
        (*derivs)(xh,yout,dydx,r,m,part);
        for (i=1;i<=n;i=i+2) yout[i+1]=yout[i+1]+hh*dydx[i+1];
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
