###############################################################################
#
# Boundary layer parameter set
# This function passed the FT test
# This function is documented in the Leaf page
#
###############################################################################
""" Create a fixed boundary layer with 0 resistance """
BLFixed(FT) = LeafBLParaSetFixed{FT, 0.0}()
