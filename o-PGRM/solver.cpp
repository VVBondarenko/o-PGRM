#include "basis.cpp"

class solver
{
private:
    gsl_matrix *sys;
    gsl_vector *rightpart, *solution;

    int type_of_problem;
    basis *basis_of_system;

    virtual double rightpart_f  (double x, double y) = 0;
    virtual double boundary_phi (double x, double y) = 0;
    virtual double omega        (double x, double y) = 0;
    virtual double omega2       (double x, double y) = 0;

    double left_under_int(double x, double y, int n, int m);
    double right_under_int(double x, double y, int n);
    //possibly, should use earlier defined stucture instead

    double structure    (double x, double y, int n);

public:
    solver(int basisType, int basisN, rect_area area, int Boundary_problem);

/*    solver(int basisN,
           //double (*rightpart_f)(double, double),
           double (*boundary_phi)(double, double),
           double (*omega)(double, double),
           double (*omega2)(double, double),
           int boundary_problem);
    ~solver();*/

    void form_matrix();
    void form_rightpart();
    void form_system();
    void solve();
    void plot(int format);

    double value_at         (double x, double y);
};
solver::solver(int basisType, int basisN, rect_area area, int Boundary_problem)
{
    this->type_of_problem = Boundary_problem;
    this->basis_of_system = new basis(basisType,basisN,area);


}
double solver::structure(double x, double y, int n)
{
        return basis_of_system->value_temp(x,y,n)*omega(x,y);
}


void solver::solve()
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (basis_of_system->N*basis_of_system->N);
    gsl_linalg_LU_decomp (sys, p, &i);
    gsl_linalg_LU_solve (sys, p, rightpart, solution);
}





