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

    virtual double rightpart_f  (double x, double y) = 0;
    virtual double boundary_phi (double x, double y) = 0;
    virtual double omega        (double x, double y) = 0;
    virtual double omega2       (double x, double y) = 0;

    double left_under_int(basis_args args);//inited
    double right_under_int(basis_args args);//inited
    double gauss_integral_left(int dimension, basis_args Args);//inited
    double gauss_integral_rigth(int dimension, basis_args Args);//inited

    double structure            (double x, double y, int n);//inited but should be corrected

public:
    solver(int basisType, int basisN, rect_area area, int Boundary_problem);//inited

    void form_matrix();//inited
    void form_rightpart();//inited
    void form_system();//inited
    void solve();//inited
    double value_at         (double x, double y);//tbd
    void plot(int format);//tbd
};


solver::solver(int basisType, int basisN, rect_area area, int Boundary_problem)
{
    this->type_of_problem = Boundary_problem;
    this->basis_of_system = new basis(basisType,basisN,area);


    sys 		= gsl_matrix_alloc(basisN*basisN,basisN*basisN);
    rightpart	= gsl_vector_alloc(basisN*basisN);
    solution	= gsl_vector_alloc(basisN*basisN);

    //initializing step parameters
    diff_step 	= pow(2.,-9);
    glob_delta 	= 1./diff_step;

    //initializing node & weights
    nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
    nodes[0] = -nodes[3];
    nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
    nodes[1] = -nodes[2];

    weights[2] = 0.5+sqrt(30)/36;
    weights[3] = 0.5-sqrt(30)/36;
    weights[0] = weights[3];
    weights[1] = weights[2];


}
double solver::structure(double x, double y, int n)
{
    return basis_of_system->value_temp(x,y,n)*omega(x,y);
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
double solver::gauss_integral_rigth(int dimension, basis_args Args)
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
                res += weights[j]*gauss_integral_rigth(1, temp_args);
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

    basis_args args = {0,0,0,0};
    for(i = 0; i < NN; i++)
    {
        args.m = i;
        for(j = 0; j < NN; j++)
        {
            args.n = j;
            gsl_matrix_set(sys, i,j, gauss_integral_rigth(2,args));
        }
    }
}
void solver::form_rightpart()
{
    int i, NN = basis_of_system->N * basis_of_system->N;

    basis_args args = {0,0,0,0};
    for(i = 0; i < NN; i++)
    {
        args.m = i;
        gsl_vector_set(rightpart, i, gauss_integral_rigth(2,args));
    }
}
void solver::form_system()
{
    int i, j, NN = basis_of_system->N * basis_of_system->N;

    basis_args args = {0,0,0,0};
    for(i = 0; i < NN; i++)
    {
        args.m = i;
        gsl_vector_set(rightpart, i, gauss_integral_rigth(2,args));
        for(j = 0; j < NN; j++)
        {
            args.n = j;
            gsl_matrix_set(sys, i,j, gauss_integral_rigth(2,args));
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

