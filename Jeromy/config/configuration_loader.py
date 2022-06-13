import yaml
import Jeromy
import Jeromy.IO.UncertaintyImporter as UncertaintyImporter


with open(Jeromy.__jeromy_home__ +"/config/my_config.yml", "r") as stream:
    configuration_loader = yaml.safe_load(stream)

def dump_configurations():
    for config in configuration_loader:
        print("--- " + config + " ---")
        ConfigurationLoader(config).print()
        print()
    
class ConfigurationLoader():

    def __init__(self, name="", systematics="", metallicity="", exposure="", background_reduction=""):
        if name != "":
            self.name                     = name
            systematic_uncertainties = configuration_loader[self.name]["systematic_uncertainties"]
            self.systematic_uncertainties = UncertaintyImporter.load_uncertainties(systematic_uncertainties)
            self.metallicity              = configuration_loader[self.name]["metallicity"]
            self.exposure                 = configuration_loader[self.name]["exposure"]
            self.background_reduction     = configuration_loader[self.name]["background_reduction"]
        
        else:
            self.name = name
            self.systematic_uncertainties = UncertaintyImporter.load_uncertainties(systematics)
            self.metallicity = metallicity
            self.exposure = exposure
            self.background_reduction = background_reduction

            
        if self.metallicity == "high":
            exec(open(f"{Jeromy.__jeromy_home__}/IO/ImportDataHighZ.py").read())
        if self.metallicity == "low":
            exec(open(f"{Jeromy.__jeromy_home__}/IO/ImportDataLowZ.py").read())

    def set_name(self, name): self.name = name
    def set_systematic_uncertainties(self, sys): self.systematic_uncertainties = UncertaintyImporter.load_uncertainties(sys)
    def set_metallicity(self, met): self.metallicity = met
    def set_exposure(self, exp): self.exposure = exp
    def set_background_reduction(self, bkg): self.background_reduction = bkg
            
        
    def print(self):
        print(self.name                    )
        print(self.systematic_uncertainties)
        print(self.metallicity             )
        print(self.exposure                )
        print(self.background_reduction    )

    
def load_configuration(configuration):
    return configuration_loader[configuration]
