import yaml
import Jeromy

with open(Jeromy.__jeromy_data__ +"/systematic_uncertainties.yml", "r") as stream:
    uncertaities_loaded = yaml.safe_load(stream)

def load_uncertainties(name):
    return uncertaities_loaded[name]
