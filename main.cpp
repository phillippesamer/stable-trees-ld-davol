#include "io.h"
#include "model.h"

#include <lemon/list_graph.h>

#include <cstdlib>

using namespace std;
using namespace lemon;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "usage: \t" << argv[0] << " [input file]" << endl << endl;
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

    Model *model = new Model(data);
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
