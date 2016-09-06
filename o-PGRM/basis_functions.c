double polynomial (double x, double y, int n)
{
    return pow(x,n%N)*pow(y,n/N);
}

double chebyshev_1d(double x, int n)
{
	//n++;
    if(n == 0)
        return 1.;
    if(n == 1)
        return x;
    if(n == 2)
        return 2.*x*x-1.;
    if(n == 3)
        return 4.*x*x*x-3.*x;
    if(n == 4)
        return 8.*x*x*(x*x-1.)+1.;
    if(n == 5)
        return 4.*x*x*x*(4.*x*x-5.)+5.*x;
    if(n == 6)
		return -1. + x*x*(18. + x*x*(-48. + 32.*x*x));
	if(n == 7)
        return x*(-7. + x*x*(56. + x*x*(-112. + 64.*x*x)));
    if(n == 8)
		return 1. + x*x*(-32. + x*x*(160. + x*x*(-256. + 128.*x*x)));
    if(n == 9)
		return x*(9. + x*x*(-120. + x*x*(432. + x*x*(-576. + 256.*x*x))));
    if(n == 10)
		return -1. + x*x*(50. + x*x*(-400. + x*x*(1120. + x*x*(-1280. + 512.*x*x))));
    if(n > 10)
        return 2.*x*chebyshev_1d(x,n-1)-chebyshev_1d(x,n-2);    
    return 0.;
}

double chebyshev_2d(double x, double y, int n)
{
    return chebyshev_1d(x/(X1-X0),n%N)*chebyshev_1d(y/(Y1-Y0),n/N);
}

double chebyshev_1dU(double x, int n)
{
	double xx = x*x;
	if(n == 0)
		return 1.;
	if(n == 1)
		return 2.*x;
	if(n == 2)
		return -1. + 4.*xx;
	if(n == 3)
		return x*(-4. + 8.*xx);
	if(n == 4)
		return 1. + xx*(-12. + 16.*xx);
	if(n == 5)
		return x*(6. + xx*(-32. + 32.*xx));
	if(n == 6)
		return -1. + xx*(24. + xx*(-80. + 64.*xx));
	if(n == 7)
		return x*(-8. + xx*(80. + xx*(-192. + 128.*xx)));
	if(n == 8)
		return 1. + xx*(-40. + xx*(240. + xx*(-448. + 256.*xx)));
	if(n == 9) 
		return x*(10. + xx*(-160. + xx*(672. + xx*(-1024. + 512.*xx))));
	if(n == 10)
		return -1. + xx*(60. + xx*(-560. + xx*(1792. + xx*(-2304. + 1024.*xx))));
    if(n > 10)
        return 2.*x*chebyshev_1dU(x,n-1)-chebyshev_1dU(x,n-2);
    return 0.;	
}
double chebyshev_2dU(double x, double y, int n)
{
    return chebyshev_1dU(x/(X1-X0),n%N)*chebyshev_1dU(y/(Y1-Y0),n/N);
}

double cubic_stepx, cubic_stepy;

double cubic_b_splines (double x, double y, int n)
{
    cubic_stepx = (X1-X0)/(double)(N-1);
    cubic_stepy = (Y1-Y0)/(double)(N-1);
    return f_B_3(0.5*cubic_stepx*(x-cubic_stepx*(double)(n%(N))))*
           f_B_3(0.5*cubic_stepy*(y-cubic_stepy*(double)(n/(N))));
}
double cubic_b_splines_2 (double x, double y, int n)
{
    //cubic_stepx = (X1-X0)/(double)(N-1);
    //cubic_stepy = (Y1-Y0)/(double)(N-1);
    return f_B_3((N-1)/(X1-X0)*(x-X0-cubic_stepx*(double)(n%(N))))*
           f_B_3((N-1)/(Y1-Y0)*(y-Y0-cubic_stepy*(double)(n/(N))));
}
double cubic_b_splines_3 (double x, double y, int n)
{
	cubic_stepx = (X1-X0)/(double)(N+1);
	cubic_stepy = (Y1-Y0)/(double)(N+1);
	return  f_B_3((x-X0-cubic_stepx*(double)(1+n/N))/cubic_stepx)*
			f_B_3((y-Y0-cubic_stepy*(double)(1+n%N))/cubic_stepy);
}
double fup3_poly (double x, double y, int n)
{
    cubic_stepx = (X1-X0)/(double)(N-1);
    cubic_stepy = (Y1-Y0)/(double)(N-1);
	return  f_fup3_poly((N-1)/(X1-X0)*(x-X0-cubic_stepx*(double)(n%(N))))*
			f_fup3_poly((N-1)/(Y1-Y0)*(y-Y0-cubic_stepy*(double)(n/(N))));
}



double fup_basis (double x, double y, int n)
{
    cubic_stepx = (X1-X0)/(double)(N-1);
    cubic_stepy = (Y1-Y0)/(double)(N-1);
    return f_fup3_poly(0.3333333333333*cubic_stepx*(x-cubic_stepx*(double)(n%(N))))*
           f_fup3_poly(0.3333333333333*cubic_stepy*(y-cubic_stepy*(double)(n/(N))));
           //0.6666666666*
}

void init_basis(int id)
{
    //	if(id == 1)
    phi = &cubic_b_splines;
    if(id == 2)
    {
        phi = &polynomial;
    }
    if(id == 3)
    {
        phi = &chebyshev_2d;
    }
    if(id == 4)
    {
        phi = &chebyshev_2dU;
    }
    if(id == 5) //the only normal B-s. basis
    {
        cubic_stepx = (X1-X0)/(double)(N-1);
        cubic_stepy = (Y1-Y0)/(double)(N-1);
        phi 		= &cubic_b_splines_2;
    }
    if(id == 6)
    {
	    phi 		= &cubic_b_splines_3;
	}
    if(id == 7)
    {
	    phi 		= &fup3_poly;
	}
}
