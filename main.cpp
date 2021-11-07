#include "io.h"
#include "model.h"
#include "ldda.h"

#include <cstdlib>

using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "usage: \t" << argv[0] << " [input file]" << endl << endl;
        return 0;
    }

    IO* instance = new IO();
    
    // parse given input file and look for errors in it
    if (instance->parse_gcclib(argv[1]) == false)
    {
        cout << "unable to parse input file" << endl;
        delete instance;
        return 0;
    }

    Model *model = new Model(instance);
    model->solve(true);

    /*
    cout << endl;
    for (long i=0; i < instance->num_edges(); i++)
    {
        pair<ModelStatus,double> probing = model->probe_var(i,1);
        if (probing.first == AT_OPTIMUM)
        {
            cout << "probing var x[" << i << "] = 1 gives ObjVal=" << probing.second 
            << " (runtime: " << model->runtime() << " s)" << endl;
        }
        else if (probing.first == IS_INFEASIBLE)
        {
            cout << "probing var x[" << i << "] = 1 gives an infeasible model" 
            << " (runtime: " << model->runtime() << " s)" << endl;
        }
    }
    cout << endl;
    */

    LDDA *lagrangean = new LDDA(instance, model);
    lagrangean->dual_ascent();

    // clean up
    delete lagrangean;
    delete instance;
    delete model;

    return 0;
}
