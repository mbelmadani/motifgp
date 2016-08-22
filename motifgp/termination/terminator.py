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

from exceptions import NotImplementedError

class Terminator(object):

    """
    Generic condition termination class interface
    """

    def __init__(self, DEBUG=False):
        self.DEBUG = DEBUG
        self.terminate = False # Termination state

    def status(self):
        """ Prints current status of the terminator """
        pass

    def done(self):
        """ Returns true when termination has been reached """
        return self.terminate    

    def update(self):
        """ Checks and update termination status if the criteria is met. """
        raise NotImplementedError
