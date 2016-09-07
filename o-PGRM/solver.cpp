#include "basis.cpp"

class solver
{
private:
    gsl_matrix *sys;//inited
    gsl_vector *rightpart, *solution;//inited
    double nodes[4], weights[4];//inited
    double diff_step, glob_delta, intStep;//inited


    int type_of_problem;    //nearly not used, should be reduced
    basis *basis_of_system; //inited

//    virtual double rightpart_f  (double x, double y) = 0;
//    virtual double boundary_phi (double x, double y) = 0;
//    virtual double omega        (double x, double y) = 0;
//    virtual double omega2       (double x, double y) = 0;

    //this methods is used for testing
    double rightpart_f  (double x, double y);
    double boundary_phi (double x, double y);
    double omega        (double x, double y);
    double omega2       (double x, double y);

    double left_under_int(basis_args args);//inited
    double left_under_int(double, double, int, int);//inited

    double right_under_int(basis_args args);//inited
    double gauss_integral_right(int dimension, basis_args Args);//inited
    double structure            (double x, double y, int n);//inited but should be corrected
    double structure            (basis_args args);//inited but should be corrected

    double SubIntegralLeft(basis_args args);
    double IntegralLeft(basis_args args);
    double SubIntegralRight(basis_args args);
    double integralRight(basis_args args);



public:
    solver(int basisType, int basisN, rect_area area, int Boundary_problem);//inited

    double gauss_integral_left(int dimension, basis_args Args);//inited

    void form_matrix();//inited
    void form_rightpart();//inited
    void form_system();//inited
    void solve();//inited
    double value_at         (double x, double y);//tbd
    void plot(/*int format*/);//tbd
    void print_solution();

};


solver::solver(int basisType, int basisN, rect_area area, int Boundary_problem)
{
    this->type_of_problem = Boundary_problem;
    this->basis_of_system = new basis(basisType,basisN,area);

    //initializing node & weights
    nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
    nodes[0] = -nodes[3];
    nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
    nodes[1] = -nodes[2];

    weights[2] = 0.5+sqrt(30)/36;
    weights[3] = 0.5-sqrt(30)/36;
    weights[0] = weights[3];
    weights[1] = weights[2];

    sys 		= gsl_matrix_alloc(basisN*basisN,basisN*basisN);
    rightpart	= gsl_vector_alloc(basisN*basisN);
    solution	= gsl_vector_alloc(basisN*basisN);

    //initializing step parameters
    diff_step 	= pow(2.,-9);
    glob_delta 	= 1./diff_step;

}
double solver::structure(double x, double y, int n)
{
    return basis_of_system->value_temp(x,y,n)*omega(x,y);
}
double solver::structure(basis_args args)
{
    return basis_of_system->value_temp(args)*omega(args.x,args.y);
}

double solver::left_under_int (basis_args args)
{
    double 	x = args.x;
    double 	y = args.y;
    int 	m = args.m;
    int 	n = args.n;

    return  structure(x,y,m)*(
                structure(x+diff_step,y,n)+structure(x-diff_step,y,n)+
                structure(x,y+diff_step,n)+structure(x,y-diff_step,n)
                -4.*structure(x,y,n))*glob_delta*glob_delta;
}
double solver::left_under_int (double x, double y, int m, int n)
{
//    double 	x = args.x;
//    double 	y = args.y;
//    int 	m = args.m;
//    int 	n = args.n;

    return  structure(x,y,m)*(
                structure(x+diff_step,y,n)+structure(x-diff_step,y,n)+
                structure(x,y+diff_step,n)+structure(x,y-diff_step,n)
                -4.*structure(x,y,n))*glob_delta*glob_delta;
}

double solver::gauss_integral_left (int dimension, basis_args Args)
{
    double x0 = basis_of_system->area.x0;
    double x1 = basis_of_system->area.x1;
    double y0 = basis_of_system->area.y0;
    double y1 = basis_of_system->area.y1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep, stepy = (y1-y0)/intStep;

    basis_args temp_args = Args;
    //argument processing

    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
                res += weights[j]*gauss_integral_left(1,temp_args);
            }
        }

        return 0.5*res*stepy;
    }
    if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*left_under_int(temp_args);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}

double solver::right_under_int(basis_args args)
{
    double 	x = args.x;
    double 	y = args.y;
    int 	m = args.m;

    return rightpart_f(x,y)*structure(x,y,m);
}
double solver::gauss_integral_right(int dimension, basis_args Args)
{
    double x0 = basis_of_system->area.x0;
    double x1 = basis_of_system->area.x1;
    double y0 = basis_of_system->area.y0;
    double y1 = basis_of_system->area.y1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep, stepy = (y1-y0)/intStep;

    basis_args temp_args = Args;
    //argument processing

    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                res += weights[j]*gauss_integral_right(1, temp_args);
            }
        }

        return 0.5*res*stepy;
    }
    if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*right_under_int(temp_args);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}


void solver::form_matrix()
{
    int i, j, NN = basis_of_system->N * basis_of_system->N;

    basis_args args;// = {0,0,0,0};
    args.x = 0;
    args.y = 0;
    args.m = 0;
    args.n = 0;

    for(i = 0; i < NN; i++)
    {
        args.m = i;
        for(j = 0; j < NN; j++)
        {
            args.n = j;
            gsl_matrix_set(sys, i,j, gauss_integral_right(2,args));
        }
    }
}
void solver::form_rightpart()
{
    int i, NN = basis_of_system->N * basis_of_system->N;

    basis_args args;// = {0,0,0,0};
    args.x = 0;
    args.y = 0;
    args.m = 0;
    args.n = 0;


    for(i = 0; i < NN; i++)
    {
        args.m = i;
        gsl_vector_set(rightpart, i, gauss_integral_right(2,args));
    }
}
void solver::form_system()
{
    int i, j, NN = basis_of_system->N * basis_of_system->N;

    basis_args args;// = {0,0,0,0};
    args.x = 0;
    args.y = 0;
    args.m = 0;
    args.n = 0;


    for(i = 0; i < NN; i++)
    {
        args.m = i;
        //gsl_vector_set(rightpart, i, gauss_integral_right(2,args));
        gsl_vector_set(rightpart, i, right_under_int(args));
        for(j = 0; j < NN; j++)
        {
            args.n = j;
//            gsl_matrix_set(sys, i,j, gauss_integral_right(2,args));
            gsl_matrix_set(sys, i,j, left_under_int(args));
        }
    }
}

void solver::solve()
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (basis_of_system->N*basis_of_system->N);
    gsl_linalg_LU_decomp (sys, p, &i);
    gsl_linalg_LU_solve (sys, p, rightpart, solution);
}
double solver::value_at(double x, double y)
{
    int i;
    double result = 0.;

    const double X0 = basis_of_system->area.x0;
    const double X1 = basis_of_system->area.x1;
    const double Y0 = basis_of_system->area.y0;
    const double Y1 = basis_of_system->area.y1;

    if(x==X0 && y==Y0)
        return value_at(x+5.*diff_step,y+5.*diff_step);
    else if(x==X0 && y==Y1)
        return value_at(x+5.*diff_step,y-5.*diff_step);
    else if(x==X1 && y==Y0)
        return value_at(x-5.*diff_step,y+5.*diff_step);
    else if(x==X1 && y==Y1)
        return value_at(x-5.*diff_step,y-5.*diff_step);
    for(i=0; i<basis_of_system->N*basis_of_system->N; i++)
    {
        result += gsl_vector_get(solution, i)*structure(x,y,i);
    }

    result+=boundary_phi(x,y);
    return result;
}


double solver::boundary_phi(double x, double y)
{
    return x*y*0.;
}

double solver::rightpart_f(double x, double y)
{
    return -2.*sin(x)*sin(y);
}

double solver::omega(double x, double y)
{
    return (x-M_PI)*(x+M_PI)*(y-M_PI)*(y+M_PI);
}

double solver::omega2(double x, double y)
{
    return x*y*0.;
}


void solver::plot()
{
    double hx = (basis_of_system->area.x1-basis_of_system->area.x0)/64.,
           hy = (basis_of_system->area.y1-basis_of_system->area.y0)/64.,
           i,j;

    FILE * op;
    op = fopen("./plot_region", "w");
    for(i=basis_of_system->area.x0; i<=basis_of_system->area.x1; i+=hx)
        for(j=basis_of_system->area.y0; j<=basis_of_system->area.y1; j+=hy)
        {
//                fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, value_at(i,j));
                fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, structure(i,j,63));
        }
    fclose(op);
    i = system("../bin/plotter.py ./plot_region Numerical &");
}

void solver::print_solution()
{
    gsl_matrix_fprintf(stdout, sys, "%f");
}



double solver::SubIntegralLeft(basis_args args)
{
    int i,j;
    double x0 = basis_of_system->area.x0;
    double x1 = basis_of_system->area.x1;


    double res = 0., step = (x1-x0)/intStep;
    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 4; j++)
        {
            args.x = (double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step;
            res += weights[j]*left_under_int(args);
        }
    }

    return 0.5*res*step;
}

double solver::IntegralLeft(basis_args args)
{
    int i,j;
    double y0 = basis_of_system->area.y0;
    double y1 = basis_of_system->area.y1;



    double res = 0., step = (y1-y0)/intStep;
    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 4; j++)
        {
            args.y = (double)(i-1)*step + y0 + 0.5*(nodes[j]+1.)*step;
            res += weights[j]*SubIntegralLeft(args);
        }
    }

    return 0.5*res*step;
}

double solver::SubIntegralRight(basis_args args)
{
    int i,j;
    double x0 = basis_of_system->area.x0;
    double x1 = basis_of_system->area.x1;


    double res = 0., step = (x1-x0)/intStep;
    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 4; j++)
        {
            args.x = (double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step;
            res += weights[j]*right_under_int(args);
        }
    }

    return 0.5*res*step;
}

double solver::integralRight(basis_args args)
{
    int i,j;
    double y0 = basis_of_system->area.y0;
    double y1 = basis_of_system->area.y1;


    double res = 0., step = (y1-y0)/intStep;
    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 4; j++)
        {
            args.y = (double)(i-1)*step + y0 + 0.5*(nodes[j]+1.)*step;
            res += weights[j]*SubIntegralRight(args);
        }
    }

    return 0.5*res*step;
}
