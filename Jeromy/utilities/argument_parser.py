import argparse

def create_parser():
    """Argument parser for the main script
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-ch", "--config-help", required=False, action="store_true",
                        help="Print the available configurations")
    parser.add_argument("-c", "--configuration", required=False,
                        help="The desired configuration to run over")
    parser.add_argument("-mp", "--make-plots", required=False, action="store_true",
                        help="Use this flag to make the spectrum plots")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Use this flag to see the terminal outputs")
    parser.add_argument("-ds", "--dump-systematics", action="store_true",
                        help="Use flag to dump systematic uncertainties")
    parser.add_argument("-e", "--exposure",
                        help="The exposure in kt-years / 10 (1=10 kt-year, 0.1=1 kt-year)")
    parser.add_argument("-br", "--background-reduction",
                        help="divisor of the background, 1=no reduction, 1000=1/1000 origional rate")
    parser.add_argument("-s", "--systematics",
                        help="The name of the systematics table to use")
    parser.add_argument("-m", "--metallicity",
                        help="high or low depending on the metallicity you want to test")
    return parser.parse_args()
