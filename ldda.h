#ifndef _LDDA_H_
#define _LDDA_H_

#include "model.h"

/***
 * \file ldda.h
 * 
 * Module for the operations within Lagrangean Decomposition based dual ascent.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 02.11.2021
 */
class LDDA
{
public:
    double dual_ascent();
private:
    double edge_deletion_bound();
    double edge_contraction_bound();
    double vertex_deletion_bound();
    double vertex_fix_bound();
};

#endif
