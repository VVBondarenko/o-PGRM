
class solver
{
private:
    gsl_matrix *sys;
    gsl_vector *rightpart, *solution;

    basis basis_of_system;

    double (*rightpart_f)  (double x, double y);
    double (*boundary_f)   (double x, double y);
    double (*omega)        (double x, double y);
    double (*omega2)       (double x, double y);
    double (*structure)    (double x, double y, int n);

public:
    solver();
    ~solver();

    void form_matrix();
    void form_rightpart();
    void form_system();
    void solve();
    double value_at         (double x, double y);
};

