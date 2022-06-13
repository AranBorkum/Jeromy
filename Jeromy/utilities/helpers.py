import os
import math
import numpy as np

def poisson_likelihood(exp, obs):

    output = 0
    if type(exp) is list and type(obs) is list:
        for i in range(len(list)):
            output += (exp[i] - obs[i] + obs[i]*np.log(obs[i]/exp[i]))
        return 2*output

    else:
        return 2 * (exp - obs + obs*np.log(obs/exp))

def gaussian(x, m, s, a):
    term1 = a / (s * np.sqrt(2*np.pi))
    term2 = (x-m)**2
    term3 = 2 * (s**2)
    return term1 * np.exp(-term2 / term3)

def add_in_quadrature(list):
    output = 0
    for i in list:
        output += i**2
    output = np.sqrt(output)
    return output

def directory_exists(path):
    return os.path.exists(path)

def ensure_directory_exists(path):
    if not directory_exists(path):
        os.makedirs(path)
        print(f"{path} created")

def get_systematics(config):
    sys_unc = config.systematic_uncertainties
    O15_tot_syst_unc     = [sys_unc["O15_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"], sys_unc["cross_section"]]
    F17_tot_syst_unc     = [sys_unc["F17_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"], sys_unc["cross_section"]]
    B8_CC_tot_syst_unc   = [sys_unc["B8_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"], sys_unc["cross_section"]]
    B8_ES_tot_syst_unc   = [sys_unc["B8_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"], sys_unc["cross_section"]]
    HEP_CC_tot_syst_unc  = [sys_unc["HEP_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"], sys_unc["cross_section"]]
    HEP_ES_tot_syst_unc  = [sys_unc["HEP_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    Radon_tot_syst_unc   = [sys_unc["radon"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    Neutron_tot_syst_unc = [sys_unc["neutron"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    Ar42_tot_syst_unc    = [sys_unc["argon_upper"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    
    O15_tot_syst_unc     = add_in_quadrature(O15_tot_syst_unc    )
    F17_tot_syst_unc     = add_in_quadrature(F17_tot_syst_unc    )
    B8_CC_tot_syst_unc   = add_in_quadrature(B8_CC_tot_syst_unc  )
    B8_ES_tot_syst_unc   = add_in_quadrature(B8_ES_tot_syst_unc  )
    HEP_CC_tot_syst_unc  = add_in_quadrature(HEP_CC_tot_syst_unc )
    HEP_ES_tot_syst_unc  = add_in_quadrature(HEP_ES_tot_syst_unc )
    Radon_tot_syst_unc   = add_in_quadrature(Radon_tot_syst_unc  )
    Neutron_tot_syst_unc = add_in_quadrature(Neutron_tot_syst_unc)
    Ar42_tot_syst_unc    = add_in_quadrature(Ar42_tot_syst_unc   )
    
    systematic_uncertainties = {"O15_syst"    : O15_tot_syst_unc,    
                                "F17_syst"    : F17_tot_syst_unc,    
                                "B8_CC_syst"  : B8_CC_tot_syst_unc,  
                                "B8_ES_syst"  : B8_ES_tot_syst_unc,  
                                "HEP_CC_syst" : HEP_CC_tot_syst_unc, 
                                "HEP_ES_syst" : HEP_ES_tot_syst_unc, 
                                "Radon_syst"  : Radon_tot_syst_unc,  
                                "Neutron_syst": Neutron_tot_syst_unc,
                                "Ar42_syst"   : Ar42_tot_syst_unc}   

    return systematic_uncertainties
        
    
    
def dump_systematics(config):
    sys_unc = config.systematic_uncertainties
    print("{: <15} {: <15}".format("component", "pecentage"))
    print("-"*30)
    for i in sys_unc: print("{: <15} {: <15}".format(i, sys_unc[i]))
    print()
    
    O15_tot_syst_unc     = [sys_unc["O15_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    F17_tot_syst_unc     = [sys_unc["F17_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    B8_CC_tot_syst_unc   = [sys_unc["B8_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    B8_ES_tot_syst_unc   = [sys_unc["B8_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    HEP_CC_tot_syst_unc  = [sys_unc["HEP_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    HEP_ES_tot_syst_unc  = [sys_unc["HEP_flux_e"], sys_unc["Sur_prob_e"],
                            sys_unc["trigger_e"], sys_unc["reconstruction"]]
    Radon_tot_syst_unc   = [sys_unc["radon"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    Neutron_tot_syst_unc = [sys_unc["neutron"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    Ar42_tot_syst_unc    = [sys_unc["argon_upper"], sys_unc["reconstruction"],
                            sys_unc["trigger_e"]]
    
    O15_tot_syst_unc     = add_in_quadrature(O15_tot_syst_unc    )
    F17_tot_syst_unc     = add_in_quadrature(F17_tot_syst_unc    )
    B8_CC_tot_syst_unc   = add_in_quadrature(B8_CC_tot_syst_unc  )
    B8_ES_tot_syst_unc   = add_in_quadrature(B8_ES_tot_syst_unc  )
    HEP_CC_tot_syst_unc  = add_in_quadrature(HEP_CC_tot_syst_unc )
    HEP_ES_tot_syst_unc  = add_in_quadrature(HEP_ES_tot_syst_unc )
    Radon_tot_syst_unc   = add_in_quadrature(Radon_tot_syst_unc  )
    Neutron_tot_syst_unc = add_in_quadrature(Neutron_tot_syst_unc)
    Ar42_tot_syst_unc    = add_in_quadrature(Ar42_tot_syst_unc   )

    print("{: <15} {: <15}".format("component", "percentage"                  ))
    print("-"*30)
    print("{: <15} {: <15}".format("O15"      , round(O15_tot_syst_unc    , 2)))
    print("{: <15} {: <15}".format("F17"      , round(F17_tot_syst_unc    , 2)))
    print("{: <15} {: <15}".format("B8 CC"    , round(B8_CC_tot_syst_unc  , 2)))
    print("{: <15} {: <15}".format("B8 ES"    , round(B8_ES_tot_syst_unc  , 2)))
    print("{: <15} {: <15}".format("HEP CC"   , round(HEP_CC_tot_syst_unc , 2)))
    print("{: <15} {: <15}".format("HEP ES"   , round(HEP_ES_tot_syst_unc , 2)))
    print("{: <15} {: <15}".format("Radon"    , round(Radon_tot_syst_unc  , 2)))
    print("{: <15} {: <15}".format("Neutron"  , round(Neutron_tot_syst_unc, 2)))
    print("{: <15} {: <15}".format("Argon"    , round(Ar42_tot_syst_unc   , 2)))


p_to_sigma = {0.682689492137086: 1,
              0.954499736103642: 2,
              0.997300203936740: 3,
              0.999936657516334: 4,
              0.999999426696856: 5,
              0.999999998026825: 6,
              0.999999999997440: 7}
    
def p_value_to_sigma(value):

    for i in p_to_sigma:
        if i < abs(1-value):
            return p_to_sigma[i]
        else:
            return 0

