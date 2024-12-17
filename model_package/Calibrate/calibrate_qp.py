#!python
import argparse
import pathlib
import sys
import yaml

from itertools import compress
import numpy
import matplotlib.pyplot
import pandas
import scipy

sys.path.append(r'/projects/tea/tardigrade_GED/tardigrade_micromorphic_element/src/python')
sys.path.append(r'/projects/tea/tardigrade-examples/model_package')
sys.path.append(r'/projects/tea/tardigrade-examples/model_package/Calibrate')
sys.path.append(f'/projects/tea/tardigrade_GED/tardigrade_micromorphic_linear_elasticity/src/python')

import calibrate_element
import calibration_tools
import micromorphic
import xdmf_reader_tools as XRT
import linear_elastic_parameter_constraint_equations as constraints


elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                              'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                              'tau8', 'tau9', 'tau10', 'tau11']


def calibrate_qp(input_file, output_file, case, Emod, nu, L, element=0, qp=0, increment=None, plot_file=None,
              average=False, UQ_file=None, cal_norm='L1', bound_half_width=1.e5, dev_norm_errors=False, input_elastic_parameters=None):
    ''' Unpack DNS data and run calibration routine

    :param str input_file: The homogenized XDMF file output by the Micromorphic Filter
    :param str output_file:  The resulting list of parameters stored in a yaml file
    :param int case: The calibration "case".
                     1: two parameter,
                     2: 7 parameter,
                     3: 7 parameter plus tau7 without error for M,
                     4: all 18 parameters,
                     5: 7 parameter plus tau7 with error for M,
                     6: 11 higher order parameters,
                     7: 7 parameters using fixed higher order parameters determined from case 6,
                     8: 7 parameters using initial guess and tighter bounds for higher order parameters determined from case 6
    :param float Emod: Estimate of a homogenized elastic modulus, used for initial parameter estimation
    :param float nu: Estimate of a homogenized Poisson ratio, used for initial parameter estimation
    :param float L: DNS max dimension (width, height, depth, etc.), used for initial parameter estimation
    :param int element: The macro (filter) element to calibration, default is zero
    :param int increment: An optional list of one or more increments to perform calibration
    :param str plot_file: Optional root filename to for plotting results
    :param bool average: Boolean whether or not homogenized DNS results will be averaged
    :param str UQ_file: Optional csv filename to store function evaluations and parameter sets for UQ
    :param str cal_norm: The type of norm to use for calibration ("L1", "L2", or "L1-L2")
    :param float bound_half_width: The uniform parameter bound "half-width" to apply for all parameters to be calibrated. Bounds for lambda will be [0., bound_half_width]. All other parameter bounds will be [-1*bound_half_width, bound_half_width]
    :param bool dev_norm_errors: Boolean whether to inclue deviatoric stress norms during calibration
    :param str input_elastic_parameters: Yaml file containing previously calibrated elastic parameters

    :returns: calibrated parameters by minimizing a specified objective function
    '''

    PK2_sdstore = []
    SIGMA_sdstore = []
    M_sdstore = []
 
    # Read in the data
    filename = input_file
    data, geometry, topology = XRT.parse_xdmf_output(filename)
    nqp = 8
    ninc = numpy.shape(data['time'])[0]
    nel = numpy.shape(data['cauchy_stress_0_0'])[0]
 
    # Read in the position information
    displacement, gradu, phi, gradphi = XRT.construct_degrees_of_freedom(data, nqp, nel)

    # Read in the stress information
    cauchy, symm, m = XRT.collect_stresses(data, nqp, nel)
    PK2, SIGMA, M = XRT.get_reference_configuration_stresses(data, nqp, nel)

    # Read in the strain information
    E, Ecal, Gamma, F, chi, grad_chi, estrain, h = XRT.compute_deformations(data, nqp, nel)

    # Get times
    times = numpy.unique(data['time'])
    ninc = len(times)
    print(f'times = {times}')
    
    # isolate element and qp
    cauchy = calibration_tools.isolate_element_and_qp(cauchy, '3x3', element, qp)
    symm = calibration_tools.isolate_element_and_qp(symm, '3x3', element, qp)
    PK2 = calibration_tools.isolate_element_and_qp(PK2, '3x3', element, qp)
    SIGMA = calibration_tools.isolate_element_and_qp(SIGMA, '3x3', element, qp)
    E = calibration_tools.isolate_element_and_qp(E, '3x3', element, qp)
    displacement = calibration_tools.isolate_element_and_qp(displacement, '3', element, qp)
    gradu = calibration_tools.isolate_element_and_qp(gradu, '3x3', element, qp)
    phi = calibration_tools.isolate_element_and_qp(phi, '3x3', element, qp)
    estrain = calibration_tools.isolate_element_and_qp(estrain, '3x3', element, qp)
    h = calibration_tools.isolate_element_and_qp(h, '3x3', element, qp)
    gradphi = calibration_tools.isolate_element_and_qp(gradphi, '3x3x3', element, qp)
    # reset nqp and qp since quantities are isolated
    nqp = 1
    qp = 0

    # store data for calibration
    Y = [PK2, SIGMA, M]
    inputs = [E, displacement, gradu, phi, gradphi, times]
    
    # get target nu from E
    # define time steps to calibrate against
    print(f'increment = {increment}')
    if increment == None:
        nu_inc = ninc - 1
    elif increment and (len(increment) == 1):
        nu_inc = int(increment[0])
    elif increment and (len(increment) > 1):
        nu_inc = int(increment[-1])
    else:
        print('something went wrong determining the increment for calculation Poisson ratio')
    nu_targ = -1*numpy.average([E[qp][nu_inc,0,0,0],E[qp][nu_inc,0,1,1]])/E[qp][nu_inc,0,2,2]

    # Estimate initial parameters
    if case == 1:
        param_est = calibration_tools.Isbuga_micrormorphic_elasticity_parameters(Emod, nu, L, case_1_override=True)
    else:
        param_est = calibration_tools.Isbuga_micrormorphic_elasticity_parameters(Emod, nu, L)
 
    # Define the elastic bounds
    upper = bound_half_width
    lower = -1*upper
    parameter_bounds = [[0.0, upper]] + [[lower, upper] for _ in range(6)] + [[lower, upper] for _ in range(11)]
    if len(elastic_parameter_ordering) != len(parameter_bounds):
        raise ValueError(f"The parameter bounds and the parameter names do not have the same length {len(parmaeter_bounds)} vs. {len(elastic_parameter_ordering)}")

    # calibrate!
    maxit = 2000
    num_workers = 8
    if case == 1:
        # Case 1 - calibrate just lambda and mu
        print(f'Target Poisson ratio = {nu_targ}')
        if numpy.isclose(nu_targ, 0.0, atol=1e-4):
            print(f'nu_targ is too close to zero to calibrate lambda')
            nu_targ = 0
            lamb_cal = False
        else:
            lamb_cal = True
        param_mask = [lamb_cal, True, False, False, False, False, False,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        #param_est = [59.25, 70.395]
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_1,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=1,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = calibrate_element.opti_options_1(list(res.x), Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 2:
        # Case 2 - calibrate first 7 parameters
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_2,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_2(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 3:
        # Case 3 - calibrate first 7 parameters and tau 7, do not use error for M in objective function
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 4:
        # Case 4 - calibrate all parameters simultaneously
        param_mask = [True, True, True, True, True, True, True,
                      True, True, True, True, True, True, True, True, True, True, True]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 5:
        # Case 5 - Same as case 3, but include M error in objective function
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 6:
        # Case 6 - Calibrate only the 11 higher order tau parmaeters
        param_mask = [False, False, False, False, False, False, False,
                      True, True, True, True, True, True, True, True, True, True, True]
        second_order_params = param_est[0:7]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, second_order_params, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_6(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, second_order_params, calibrate=False)
    elif case == 7:
        # Case 7 - Calibrate 7 second order parmaeters and read in existing higher order parameters which are fixed
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        # unpack previously calibrated parameters
        assert input_elastic_parameters != None, "input_elastic_parameters must be provided!"
        input_parameters, _ = calibration_tools.parse_input_parameters(input_elastic_parameters)
        third_order_params = input_parameters[10:21]
        print(f'input third_order_parameters = {third_order_params}')
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_7,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, third_order_params, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_7(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, third_order_params, calibrate=False)
    elif case == 8:
        # Case 8 - Calibrate all 18 parameters and read in existing higher order parameters which as an initial estimate with adjusted bounds
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        print('initial parameter estimation:')
        print(param_est)
        # unpack previously calibrated parameters
        assert input_elastic_parameters != None, "input_elastic_parameters must be provided!"
        input_parameters, _ = calibration_tools.parse_input_parameters(input_elastic_parameters)
        third_order_params = input_parameters[10:21]
        print(f'input third_order_parameters = {third_order_params}')
        param_est = numpy.hstack([param_est, third_order_params])
        third_order_bounds = [sorted([0.1*p, 10.*p]) for p in third_order_params]
        parameter_bounds = parameter_bounds + third_order_bounds
        res = scipy.optimize.differential_evolution(func=calibrate_element.opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = calibrate_element.opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    else:
        print('Select valid calibration case!')

    # Make a csv of all function evaluations and parameter sets
    if UQ_file:
        print(f'shape of Xstore = {numpy.shape(Xstore)}')
        print(f'shape of Lstore = {numpy.shape(Lstore)}')
        UQ_dict = handle_output_for_UQ(Xstore, Lstore, case)
        df = pandas.DataFrame(UQ_dict)
        df.to_csv(UQ_file, header=True, sep=',', index=False)

    # look at population energy info
    #population = res.population
    #energies = res.population_energies
    #print(f'size of population = {numpy.shape(population)}')
    #print(f'size of population_energies = {numpy.shape(energies)}')
    #print(f'population = \n {population}\n')
    #print(f'energies = \n {energies}\n')

    # plot resulting calibration
    if plot_file:
        print('plotting...')
        model_name=r'LinearElasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibration_tools.evaluate_model(inputs, params, model_name, calibrate_element.parameters_to_fparams, 0, element, nqp)
        PK2_sim = XRT.map_sim(PK2_sim, ninc)
        SIGMA_sim = XRT.map_sim(SIGMA_sim, ninc)
        M_sim = XRT.map_sim(M_sim, ninc, third_order=True)
        cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, inputs[2], inputs[3])

        calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}.PNG',
                                        element, nqp, 'E', 'S', increment=increment)
        calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}.PNG',
                                        element, nqp, '\mathcal{E}', '\Sigma', increment=increment)
        calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}.PNG',
                                                     element, nqp, increment=increment)
        calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}_bounded.PNG',
                                        element, nqp, 'E', 'S',increment=increment, find_bounds=True)
        calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}_bounded.PNG',
                                        element, nqp, '\mathcal{E}', '\Sigma', increment=increment, find_bounds=True)
        calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}_bounded.PNG',
                                                     element, nqp, increment=increment, find_bounds=True)
        calibration_tools.plot_stress_norm_calibration_comparison(PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim, E, Ecal, Gamma,
                                                                  f'{plot_file}_norms_case_{case}.PNG', nqp, increment=increment)
        # plot the full range of predictions if specific increments were speficied
        if increment:
            calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}_ALL.PNG',
                                            element, nqp, 'E', 'S')
            calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}_ALL.PNG',
                                            element, nqp, '\mathcal{E}', '\Sigma')
            calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}_ALL.PNG',
                                                         element, nqp)
            calibration_tools.plot_stress_norm_calibration_comparison(PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim, E, Ecal, Gamma,
                                                                      f'{plot_file}_norms_case_{case}_ALL.PNG', nqp)

    # output parameters
    output_filename = output_file
    output_dict = {}
    p = params
    output_dict['line 1'] = f"2 {p[0]} {p[1]}"
    output_dict['line 2'] = f"5 {p[2]} {p[3]} {p[4]} {p[5]} {p[6]}"
    output_dict['line 3'] = f"11 {p[7]} {p[8]} {p[9]} {p[10]} {p[11]} {p[12]} {p[13]} {p[14]} {p[15]} {p[16]} {p[17]}"
    output_dict['line 4'] = f"2 {p[3]} {p[6]}"
    output_dict['obj_func_value'] = f"{res.fun}"
    with open(output_filename, 'w') as f:
        yaml.dump(output_dict, f)
    return


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Calibrate micromorphic linear elasticity on a single filter domain (i.e. macroscale element) and quadrature point"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-i', '--input-file', type=str,
        help="The homogenized XDMF file output by the Micromorphic Filter")
    parser.add_argument('-o', '--output-file', type=str,
        help="The resulting list of parameters stored in a yaml file")
    parser.add_argument('--Emod', type=float,
        help="DNS elastic modulus, used for initial parameter estimation.")
    parser.add_argument('--nu', type=float,
        help="DNS Poisson's ratio, used for initial parameter estimation.")
    parser.add_argument('--L', type=float,
        help="DNS max dimension (width, height, depth, etc.), used for initial parameter estimation.")
    parser.add_argument('--element', type=int, default=0,
        help="The macro (filter) element to calibrate")
    parser.add_argument('--qp', type=int, default=0,
        help="The quadrature point of the macro (filter) element to calibrate")
    parser.add_argument('--increment', nargs="+", required=False, default=None,
        help="An optional list of one or more increments to perform calibration")
    parser.add_argument('--case', type=int, required=True,
        help="The calibration 'case'.\
              1: two parameter,\
              2: 7 parameter,\
              3: 7 parameter plus tau7 without error for M,\
              4: all 18 parameters,\
              5: 7 parameter plus tau7 with error for M,\
              6: 11 higher order parameters,\
              7: 7 parameters using fixed higher order parameters determined from case 6,\
              8: 7 parameters using initial guess and tighter bounds for higher order parameters determined from case 6")
    parser.add_argument('--plot-file', type=str, required=False, default=None,
        help="Optional root filename to for plotting results")
    parser.add_argument('--UQ-file', type=str, required=False,
        help='Optional csv filename to store function evaluations and parameter sets for UQ')
    parser.add_argument('--cal-norm', type=str, required=False, default='L1',
        help='The type of norm to use for calibration ("L1", "L2", or "L1-L2")')
    parser.add_argument('--bound-half-width', type=float, required=False, default=1.e5,
        help='The uniform parameter bound "half-width" to apply for all parameters to be calibrated.\
              Bounds for lambda will be [0., bound_half_width].\
              All other parameter bounds will be [-1*bound_half_width, bound_half_width]')
    parser.add_argument('--dev-norm-errors', type=str, required=False, default=False,
        help='Boolean whether to inclue deviatoric stress norms during calibration')
    parser.add_argument('--input-elastic-parameters', type=str, required=False, default=None,
        help='Yaml file containing previously calibrated elastic parameters')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    # Objective function evaluation lists
    Xstore = []
    Lstore = []

    args, unknown = parser.parse_known_args()
    sys.exit(calibrate_qp(input_file=args.input_file,
                          output_file=args.output_file,
                          Emod=args.Emod,
                          nu=args.nu,
                          L=args.L,
                          element=args.element,
                          qp=args.qp,
                          increment=args.increment,
                          case=args.case,
                          plot_file=args.plot_file,
                          UQ_file=args.UQ_file,
                          cal_norm=args.cal_norm,
                          bound_half_width=args.bound_half_width,
                          dev_norm_errors=calibrate_element.str2bool(args.dev_norm_errors),
                          input_elastic_parameters=args.input_elastic_parameters,
                          ))