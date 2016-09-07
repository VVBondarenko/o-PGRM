#include "af_poly.c"
#include "B-splines.c"
#include "basis_functions.cpp"

typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;

typedef struct basis_args {
    double	x, y;
    int		m, n;
} basis_args;


class basis{
public:
    basis(int basis_type, int nn, double x0, double x1, double y0, double y1);
    basis(int basis_type, int nn, rect_area Area);
    //    ~basis();
    int N;
    rect_area area;
    double value_temp   (double x, double y, int n);
    double value_temp   (basis_args args);
    double (basis::*value)     (double x, double y, int n);
private:
    int Basis_Type;
    double cubic_stepx, cubic_stepy;
    double simplePoly   (double x, double y, int n);
    double chebyshevT   (double x, double y, int n);
    double chebyshevU   (double x, double y, int n);
    double Bsplines     (double x, double y, int n);
    double fup3_poly    (double x, double y, int n);
};

double basis::simplePoly (double x, double y, int n)
{
    return pow(x,n%this->N)*pow(y,n/this->N);
}

double basis::chebyshevT(double x, double y, int n)
{
    return chebyshev_1d(x/(area.x1-area.x0),n%N)*chebyshev_1d(y/(area.y1-area.y0),n/N);
}

double basis::chebyshevU(double x, double y, int n)
{
    return chebyshev_1dU(x/(area.x1-area.x0),n%N)*chebyshev_1dU(y/(area.y1-area.y0),n/N);
}

double basis::Bsplines (double x, double y, int n)
{
    return f_B_3((N-1)/(area.x1-area.x0)*(x-area.x0-cubic_stepx*(double)(n%(N))))*
            f_B_3((N-1)/(area.y1-area.y0)*(y-area.y0-cubic_stepy*(double)(n/(N))));
}

double basis::fup3_poly (double x, double y, int n)
{
    return f_fup3_poly((N-1)/(area.x1-area.x0)*(x-area.x0-cubic_stepx*(double)(n%(N))))*
            f_fup3_poly((N-1)/(area.y1-area.y0)*(y-area.y0-cubic_stepy*(double)(n/(N))));
}

double basis::value_temp(double x, double y, int n)
{
    if(Basis_Type == 1)
        return Bsplines(x,y,n);
    if(Basis_Type == 2)
        return simplePoly(x,y,n);
    if(Basis_Type == 3)
        return chebyshevT(x,y,n);
    if(Basis_Type == 4)
        return chebyshevU(x,y,n);
    if(Basis_Type == 5)
        return fup3_poly(x,y,n);
    return 0.;
}

double basis::value_temp(basis_args args)
{
    double x = args.x, y = args.y;
    int n = args.n;
    if(Basis_Type == 1)
        return Bsplines(x,y,n);
    if(Basis_Type == 2)
        return simplePoly(x,y,n);
    if(Basis_Type == 3)
        return chebyshevT(x,y,n);
    if(Basis_Type == 4)
        return chebyshevU(x,y,n);
    if(Basis_Type == 5)
        return fup3_poly(x,y,n);
    return 0.;
}

basis::basis(int basis_type, int nn, double x0, double x1, double y0, double y1)
{
    this->N = nn;
    this->Basis_Type = basis_type;
    this->area.x0 = x0;
    this->area.x1 = x1;
    this->area.y0 = y0;
    this->area.y1 = y1;
    cubic_stepx = (area.x1-area.x0)/(double)(N-1);
    cubic_stepy = (area.y1-area.y0)/(double)(N-1);
    basis::value = &basis::Bsplines;
    if(basis_type == 2)
        basis::value = &basis::simplePoly;
    if(basis_type == 3)
        basis::value = &basis::chebyshevT;
    if(basis_type == 4)
        basis::value = &basis::chebyshevU;
    if(basis_type == 5)
        basis::value = &basis::fup3_poly;
}
basis::basis(int basis_type, int nn, rect_area Area)
{
    this->N = nn;
    this->Basis_Type = basis_type;
    this->area = Area;
    //    this->area.x0 = x0;
    //    this->area.x1 = x1;
    //    this->area.y0 = y0;
    //    this->area.y1 = y1;
    cubic_stepx = (area.x1-area.x0)/(double)(N-1);
    cubic_stepy = (area.y1-area.y0)/(double)(N-1);
    basis::value = &basis::Bsplines;
    if(basis_type == 2)
        basis::value = &basis::simplePoly;
    if(basis_type == 3)
        basis::value = &basis::chebyshevT;
    if(basis_type == 4)
        basis::value = &basis::chebyshevU;
    if(basis_type == 5)
        basis::value = &basis::fup3_poly;
}

