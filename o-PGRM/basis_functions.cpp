//double polynomial (double x, double y, int n)
//{
//    return pow(x,n%N)*pow(y,n/N);
//}

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
