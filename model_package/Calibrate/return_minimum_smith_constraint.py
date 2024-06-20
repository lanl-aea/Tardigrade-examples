
import numpy

import linear_elastic_parameter_constraint_equations as constraints


def return_minimum_smith_constraint(parameters, svals=None):
    '''Evaluate the 13 Smith constraints for micromorphic linear elasticity and return the minimum value

    :params list parameters: The micromorphic linear elasticity parameters
    :params list svals: TODO! Figure out what this is

    :returns: Minimum value for the 13 Smith constraints
    '''

    elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                                  'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                                  'tau8', 'tau9', 'tau10', 'tau11']
    parameter_dictionary = dict(zip(elastic_parameter_ordering, parameters[:18]))
    consts = [constraints.evaluate_g1, 
              constraints.evaluate_g2, 
              constraints.evaluate_g3,
              constraints.evaluate_g4,
              constraints.evaluate_g5,
              constraints.evaluate_g6, 
              constraints.evaluate_g7, 
              constraints.evaluate_g8,
              constraints.evaluate_g9,
              constraints.evaluate_g10,
              constraints.evaluate_g11,
              constraints.evaluate_g12,
              constraints.evaluate_g13]
    if (svals is None):
        svals = dict([(f's{i+1}', 0) for i in range(len(consts))])
    parameter_dictionary.update(svals)
    consts = [const(**parameter_dictionary) for const in consts]
    cvals = numpy.array([c[0] for c in consts])

    return cvals.min()