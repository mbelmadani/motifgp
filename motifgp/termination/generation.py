"""
    __author__ = "Manuel Belmadani"
    __copyright__ = "Copyright 2016, MotifGP"
    __credits__ = ["Manuel Belmadani"]

    __license__ = "LGPL v3.0"
    __version__ = "1.0.0"
    __maintainer__ = "Manuel Belmadani"
    __email__ = "mbelm006@uottawa.ca"
    __status__ = "Development"
    
"""
from termination.terminator import Terminator
import time

class GenerationTerminator(Terminator):
    """
    Terminator class for generation based termination
    """

    def __init__(self, limit=0, DEBUG=False):
        super(GenerationTerminator, self).__init__(DEBUG=DEBUG)
        self.limit = limit
        self.gen = 0
        
    def status(self, gen):
        """ Prints current status of the terminator """
        print "Generation:", gen,"|", str(t) + "/" + str(self.limit)

    def done(self):
        """ Returns true when termination has been reached """
        return self.terminate

    def update(self, gen=None):
        """
        If t is greater or equal to limit, set self.terminate to True.
        """
        self.ngen = gen
        if self.ngen >= self.limit:
            self.terminate = True
