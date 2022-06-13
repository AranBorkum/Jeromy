import matplotlib.pyplot as plt

import Jeromy
import Jeromy.IO.ImportDataHighZ as IDh
import Jeromy.IO.ImportDataLowZ as IDl;

def make_spectra(configuration):

    if configuration.metallicity == "high": module = IDh
    if configuration.metallicity == "low" : module = IDl
    
    # Make the CNO spectrum
    O15_yValues = [i*float(configuration.exposure) for i in module.O15.yValues]
    F17_yValues = [i*float(configuration.exposure) for i in module.F17.yValues]
    CNO_yValues = [i*float(configuration.exposure) for i in module.CNO.yValues]

    labels = [r"O15 $R_{evt, tot} = %.2E$" % sum(O15_yValues),
              r"F17 $R_{evt, tot} = %.2E$" % sum(F17_yValues),
              r"CNO $R_{evt, tot} = %.2E$" % sum(CNO_yValues)]
    
    
    plt.figure(f"CNO_spectrum_{configuration.name}", figsize=(16, 9))
    # plt.hist(module.O15.xValues, bins=module.O15.xValues, weights=O15_yValues, histtype="step", label=labels[0])
    # plt.hist(module.F17.xValues, bins=module.F17.xValues, weights=F17_yValues, histtype="step", label=labels[1])
    # plt.hist(module.CNO.xValues, bins=module.CNO.xValues, weights=CNO_yValues, histtype="step", label=labels[2])

    plt.hist(module.O15.xValues, bins=40, weights=O15_yValues, histtype="step", label=labels[0])
    plt.hist(module.F17.xValues, bins=40, weights=F17_yValues, histtype="step", label=labels[1])
    CNO_hist = plt.hist(module.CNO.xValues, bins=40, weights=CNO_yValues, histtype="step", label=labels[2])

    plt.yscale("log")
    plt.xlabel("Reconstructed neutrino energy [MeV]", fontsize=16)
    plt.ylabel(f"Event rate [evts / {float(configuration.exposure) * 10} kt-year / 1 MeV]", fontsize=16)
    plt.legend(loc="best")
    plt.savefig(f"{Jeromy.__jeromy_output__}/CNO_spectrum_{configuration.name}.png")
    plt.close()
    
    # Make the radiological spectra
    if int(configuration.background_reduction): reduction = 1/int(configuration.background_reduction)
    else:                                       reduction = 1
    
    Radon_yValues   = [i*float(configuration.exposure)*reduction for i in module.Radon        .yValues]
    Neutron_yValues = [i*float(configuration.exposure)*reduction for i in module.Neutron      .yValues]
    Ar42_yValues    = [i*float(configuration.exposure)*reduction for i in module.Ar42         .yValues]
    RadTot_yValues  = [i*float(configuration.exposure)*reduction for i in module.Radiologicals.yValues]
    labels = [r"Radon   $R_{evt, tot} = %.2E$" % sum(Radon_yValues  ),
              r"Neutron $R_{evt, tot} = %.2E$" % sum(Neutron_yValues),
              r"Ar42    $R_{evt, tot} = %.2E$" % sum(Ar42_yValues   ),
              r"Total   $R_{evt, tot} = %.2E$" % sum(RadTot_yValues )]              

    plt.figure(f"Radiological_spectrum_{configuration.name}", figsize=(16, 9))
    # plt.hist(module.Radon  .xValues, bins=module.Radon  .xValues, weights=Radon_yValues  , histtype="step", label=labels[0])
    # plt.hist(module.Neutron.xValues, bins=module.Neutron.xValues, weights=Neutron_yValues, histtype="step", label=labels[1])
    # plt.hist(module.Ar42   .xValues, bins=module.Ar42   .xValues, weights=Ar42_yValues   , histtype="step", label=labels[2])
    # plt.hist(module.Radiologicals.xValues, bins=module.Radon   .xValues,
    #          weights=RadTot_yValues   , histtype="step", label=labels[3])

    Radon_hist   = plt.hist(module.Radon  .xValues, bins=range(10), weights=Radon_yValues  , histtype="step", label=labels[0])
    Neutron_hist = plt.hist(module.Neutron.xValues, bins=range(10), weights=Neutron_yValues, histtype="step", label=labels[1])
    Ar42_hist    = plt.hist(module.Ar42   .xValues, bins=range(10), weights=Ar42_yValues   , histtype="step", label=labels[2])

    plt.yscale("log")
    plt.xlabel("Reconstructed neutrino energy [MeV]", fontsize=16)
    plt.ylabel(f"Event rate [evts / {float(configuration.exposure) * 10} kt-year / 1 MeV]", fontsize=16)
    plt.legend(loc="best")
    plt.xlim([0, 8])
    plt.savefig(f"{Jeromy.__jeromy_output__}/Radiological_spectrum_{configuration.name}.png")
    plt.close()
    
    # Make the other solar spectra
    B8_CC_yValues  = [i*float(configuration.exposure) for i in module.B8_CC .yValues]
    B8_ES_yValues  = [i*float(configuration.exposure) for i in module.B8_ES .yValues]
    B8_yValues     = [i*float(configuration.exposure) for i in module.B8    .yValues]
    HEP_CC_yValues = [i*float(configuration.exposure) for i in module.HEP_CC.yValues]
    HEP_ES_yValues = [i*float(configuration.exposure) for i in module.HEP_ES.yValues]
    HEP_yValues    = [i*float(configuration.exposure) for i in module.HEP   .yValues]

    labels = [r"Boron8 CC $R_{evt, tot} = %.2E$" % sum(B8_CC_yValues),
              r"Boron8 ES $R_{evt, tot} = %.2E$" % sum(B8_ES_yValues),
              r"Boron8    $R_{evt, tot} = %.2E$" % sum(B8_yValues   ),
              r"HEP CC $R_{evt, tot} = %.2E$" % sum(HEP_CC_yValues),
              r"HEP ES $R_{evt, tot} = %.2E$" % sum(HEP_ES_yValues),
              r"HEP    $R_{evt, tot} = %.2E$" % sum(HEP_yValues   )]

    plt.figure(f"Solar_neutrino_spectrum{configuration.name}", figsize=(16, 9))
    # plt.hist(module.B8_CC .xValues, bins=module.B8_CC .xValues, weights=B8_CC_yValues , histtype="step", label=labels[0])
    # plt.hist(module.B8_ES .xValues, bins=module.B8_ES .xValues, weights=B8_ES_yValues , histtype="step", label=labels[1])
    # plt.hist(module.B8    .xValues, bins=module.B8    .xValues, weights=B8_yValues    , histtype="step", label=labels[2])
    # plt.hist(module.HEP_CC.xValues, bins=module.HEP_CC.xValues, weights=HEP_CC_yValues, histtype="step", label=labels[3])
    # plt.hist(module.HEP_ES.xValues, bins=module.HEP_ES.xValues, weights=HEP_ES_yValues, histtype="step", label=labels[4])
    # plt.hist(module.HEP   .xValues, bins=module.HEP   .xValues, weights=HEP_yValues   , histtype="step", label=labels[5])
    plt.hist(module.B8_CC .xValues, bins=40, weights=B8_CC_yValues , histtype="step", label=labels[0])
    plt.hist(module.B8_ES .xValues, bins=40, weights=B8_ES_yValues , histtype="step", label=labels[1])
    B8_hist = plt.hist(module.B8    .xValues, bins=40, weights=B8_yValues    , histtype="step", label=labels[2])
    plt.hist(module.HEP_CC.xValues, bins=40, weights=HEP_CC_yValues, histtype="step", label=labels[3])
    plt.hist(module.HEP_ES.xValues, bins=40, weights=HEP_ES_yValues, histtype="step", label=labels[4])
    HEP_hist = plt.hist(module.HEP   .xValues, bins=40, weights=HEP_yValues   , histtype="step", label=labels[5])

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Reconstructed neutrino energy [MeV]", fontsize=16)
    plt.ylabel(f"Event rate [evts / {float(configuration.exposure) * 10} kt-year / 1 MeV]", fontsize=16)
    plt.legend(loc="best")
    plt.savefig(f"{Jeromy.__jeromy_output__}/Solar_neutrino_spectrum_{configuration.name}.png")
    plt.close()

    B8_map = {}
    HEP_map = {}
    BKG_map = {}
    Radon_map = {}
    Neutron_map = {}
    Ar42_map = {}
    CNO_map = {}

    for i in range(len(CNO_hist    [0])): CNO_map    [int(CNO_hist    [1][i])] = CNO_hist    [0][i]
    for i in range(len(B8_hist     [0])): B8_map     [int(B8_hist     [1][i])] = B8_hist     [0][i]
    for i in range(len(HEP_hist    [0])): HEP_map    [int(HEP_hist    [1][i])] = HEP_hist    [0][i]
    for i in range(len(Radon_hist  [0])): Radon_map  [int(Radon_hist  [1][i])] = Radon_hist  [0][i]
    for i in range(len(Neutron_hist[0])): Neutron_map[int(Neutron_hist[1][i])] = Neutron_hist[0][i]
    for i in range(len(Ar42_hist   [0])): Ar42_map   [int(Ar42_hist   [1][i])] = Ar42_hist   [0][i]


    print(CNO_map    , "\n")
    print(B8_map     , "\n")
    print(HEP_map    , "\n")
    print(Radon_map  , "\n")
    print(Neutron_map, "\n")
    print(Ar42_map   , "\n")


    Radiologicals_xValues = [i for i in Radon_map]
    Radiologicals_yValues = []
    Backgrounds_xValues = [i for i in B8_map]
    Backgrounds_yValues = []
    TotalSignal_xValues = [i for i in B8_map]
    TotalSignal_yValues = []

    for i in Radiologicals_xValues:
        value_b = 0
        if i in Radon_map:
            value_b += Radon_map[i]

        if i in Neutron_map:
            value_b += Neutron_map[i]

        if i in Ar42_map:
            value_b += Ar42_map[i]

        Radiologicals_yValues.append(value_b)
            
    
    for i in Backgrounds_xValues:
        value_b = 0
        value_t = 0

        if i in B8_map:
            value_b += B8_map[i]
            value_t += B8_map[i]
        
        if i in HEP_map:
            value_b += HEP_map[i]
            value_t += HEP_map[i]

        if i in Radon_map:
            value_b += Radon_map[i]
            value_t += Radon_map[i]

        if i in Neutron_map:
            value_b += Neutron_map[i]
            value_t += Neutron_map[i]

        if i in Ar42_map:
            value_b += Ar42_map[i]
            value_t += Ar42_map[i]

        if i in CNO_map:
            value_t += CNO_map[i]

            
        Backgrounds_yValues.append(value_b)
        TotalSignal_yValues.append(value_t)
    
    
    labels = [r"CNO          $R_{evt, tot} = %.2E$" % sum(CNO_yValues        ),
              r"Radiological $R_{evt, tot} = %.2E$" % sum(RadTot_yValues     ),
              r"Boron8       $R_{evt, tot} = %.2E$" % sum(B8_yValues         ),
              r"HEP          $R_{evt, tot} = %.2E$" % sum(HEP_yValues        ),
              r"Backgrounds  $R_{evt, tot} = %.2E$" % sum(Backgrounds_yValues),
              r"All Signal   $R_{evt, tot} = %.2E$" % sum(TotalSignal_yValues)]
    
        
    plt.figure(f"All_signal_spectrum{configuration.name}", figsize=(16, 9))

    plt.hist(module.CNO.xValues          , bins=40, weights=CNO_yValues, histtype="step", label=labels[0])
    plt.hist(Radiologicals_xValues       , bins=Radiologicals_xValues , weights=Radiologicals_yValues   , histtype="step", label=labels[1])
    plt.hist(module.B8    .xValues       , bins=40, weights=B8_yValues    , histtype="step", label=labels[2])
    plt.hist(module.HEP   .xValues       , bins=40, weights=HEP_yValues   , histtype="step", label=labels[3])
    plt.hist(Backgrounds_xValues         , bins=Backgrounds_xValues, weights=Backgrounds_yValues   , histtype="step", label=labels[4])
    plt.hist(TotalSignal_xValues         , bins=TotalSignal_xValues, weights=TotalSignal_yValues   , histtype="step", label=labels[5])

    plt.yscale("log")
    plt.xlim([0,8])
    plt.xlabel("Reconstructed neutrino energy [MeV]", fontsize=16)
    plt.ylabel(f"Event rate [evts / {float(configuration.exposure) * 10} kt-year / 1 MeV]", fontsize=16)
    plt.legend(loc="best")
    plt.savefig(f"{Jeromy.__jeromy_output__}/All_signal_spectrum_{configuration.name}.png")
    plt.close()
