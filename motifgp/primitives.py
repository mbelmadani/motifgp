class PrimitivePSSM():
    """
    A class representing multiple positional-probability array values
    Basically this class is only here to prevent the GP grammar from creating
    recursive matrices.
    """

    def __init__(self,
                 positions):
        self.positions = positions

    def __add__(self, x):
        if len(self.positions) < 0:
            self.positions = [x]
        else:
            self.positions.append(x)
