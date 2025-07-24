#pragma once

#include "numeric_wrapper.h"
#include "vector_3d.h"
#include "vector_3d_big_float.h"

double radius_edge_ratio(const pointType *p1, const pointType *p2,
                         const pointType *p3, const pointType *p4);

double energy_dirichlet(const pointType *p1, const pointType *p2,
                        const pointType *p3, const pointType *p4);

double energy_dirichlet_big_float(const pointType *p1, const pointType *p2,
                                  const pointType *p3, const pointType *p4);
