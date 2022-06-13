
import matplotlib.pyplot as plt

import Jeromy
import Jeromy.statistics.run_chi2_tests as run_chi2_tests
import Jeromy.statistics.test_statistic as test_statistic
import Jeromy.utilities.helpers as h

def test_ensemble(configuration):

    baker_cousins = run_chi2_tests.Chi2Test(configuration = configuration,
                                            ntests        = 1000,
                                            low           = 1,
                                            high          = 3,
                                            test          = test_statistic.BakerCousinsChi,
                                            name          = "baker_cousins")

    baker_cousins_up = run_chi2_tests.Chi2Test(configuration = configuration,
                                               ntests        = 1000,
                                               low           = 1,
                                               high          = 3,
                                               test          = test_statistic.BakerCousinsChiUp,
                                               name          = "baker_cousins_up")

    baker_cousins_ll = run_chi2_tests.Chi2Test(configuration = configuration,
                                               ntests        = 1000,
                                               low           = 1,
                                               high          = 3,
                                               test          = test_statistic.BakerCousinsLL,
                                               name          = "baker_cousins_ll")

    neyman = run_chi2_tests.Chi2Test(configuration = configuration,
                                     ntests        = 1000,
                                     low           = 1,
                                     high          = 3,
                                     test          = test_statistic.Neyman,
                                     name          = "neyman")

    pearson = run_chi2_tests.Chi2Test(configuration = configuration,
                                      ntests        = 1000,
                                      low           = 1,
                                      high          = 3,
                                      test          = test_statistic.Pearson,
                                      name          = "pearson")

    extended_ll = run_chi2_tests.Chi2Test(configuration = configuration,
                                          ntests        = 1000,
                                          low           = 1,
                                          high          = 3,
                                          test          = test_statistic.ExtendedLL,
                                          name          = "extended_ll")

    
    print("Baker Cousins"   );baker_cousins_values    = baker_cousins.run_tests()
    print("Baker Cousins Up");baker_cousins_up_values = baker_cousins_up.run_tests()
    print("Baker Cousins LL");baker_cousins_ll_values = baker_cousins_ll.run_tests()
    print("Neyman"          );neyman_values           = neyman.run_tests()
    print("Pearson"         );pearson_values          = pearson.run_tests()
    print("Extended LL"     );extended_ll_values      = extended_ll.run_tests()

    h.ensure_directory_exists(f"{Jeromy.__jeromy_output__}/test_statistics_{configuration.name}")
    output_directory = f"{Jeromy.__jeromy_output__}/test_statistics_{configuration.name}"
    outfile = open(f"{output_directory}/test_statistics_{configuration.name}.txt", "w")
    outfile.write(f"baker_cousins_values,"   +
                  f"baker_cousins_up_values,"+
                  f"baker_cousins_ll_values,"+
                  f"neyman_values,"          +
                  f"pearson_values,"         +
                  f"extended_ll_values\n")

    for i in range(len(baker_cousins_values)):
        outfile.write(f"{baker_cousins_values   [i]},"+
                      f"{baker_cousins_up_values[i]},"+
                      f"{baker_cousins_ll_values[i]},"+
                      f"{neyman_values          [i]},"+
                      f"{pearson_values         [i]},"+
                      f"{extended_ll_values     [i]}\n")
    outfile.close()
                      
            
    
    plt.figure(r"$\chi^{2}$ for various methods", figsize=(16, 9))
    plt.subplot(2, 3, 1)
    plt.hist(baker_cousins_values)
    plt.title("Baker Cousins")

    plt.subplot(2, 3, 2)
    plt.hist(baker_cousins_up_values)
    plt.title("Baker Cousins Up")
    
    plt.subplot(2, 3, 3)
    plt.hist(baker_cousins_ll_values)
    plt.title("Baker Cousins Log Likelihood")

    plt.subplot(2, 3, 4)
    plt.hist(neyman_values)
    plt.title("Neyman")

    plt.subplot(2, 3, 5)
    plt.hist(pearson_values)
    plt.title("Pearson")

    plt.subplot(2, 3, 6)
    plt.hist(extended_ll_values)
    plt.title("Extended Log Likelihood")
    plt.savefig(f"{output_directory}/test_statistics_{configuration.name}.png")
    plt.close()
    
