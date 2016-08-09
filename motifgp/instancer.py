import os.path
import sys
from objectives import *


class Instancer(object):
    """
    Load configurations in a simple mapping format, Key=Value.
    Provides other functions to instance the configurations.
    """
    def __init__(self, config_path=None):
        self.config = self.load_config( config_path  ) 


    def load_config(self, config_path=None, prefix_key=False):
        """
        Given a config file path, load all configs from the mapping.
        Prefix key adds the config file basename as a prefix.
        """
        if config_path == None:
            return None

        prefix = ""
        if prefix_key:
            prefix = os.path.splitext(os.path.basename(config_path))[0] + "."
        
        configs = {}
        with open(config_path) as config_file:
            for line in config_file:
                if "#" == line[0]:
                    continue # Comment line
                elif "=" not in line:
                    print "Error one line:", line
                    print "Comments should start with # and other mapping should be Character=Instance."
                    continue
                else:
                    a,b = [x.strip() for x in line.split("=")]
                    configs[a] = prefix + b
        return configs
        
    def construct(self, key):
        """
        Construct module instance for a key in the configs.
        """
        return getattr(sys.modules[__name__], self.config[key])
