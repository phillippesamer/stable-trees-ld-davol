#include "io.h"
#include "model.h"

#include <cstdlib>

using namespace std;

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cout << "usage: \t" << argv[0] << " [input graph file] fixed-cardinality" << endl << endl;
        return 0;
    }

    IO* data = new IO();
    
    // parse given input file and look for errors in it
    if (data->parse_gcclib(argv[1]) == false)
    {
        cout << "unable to parse input file" << endl;
        delete data;
        return 0;
    }

    Model *model = new Model(data, atoi(argv[2]));
    int solutions_count = model->solve_lp_relax();
    if (solutions_count > 0)
    {
        if (model->check_half_integer_solution() == false)
            cout << "FALSE" << endl;
        else
            cout << "OK" << endl;
    }

    delete data;
    delete model;
    return 0;
}
