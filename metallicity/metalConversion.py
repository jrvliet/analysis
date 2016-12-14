

def mrToZ(Z_cell):

    # mfToZ.py
    # 
    # Converts mass fractions to Z/Z_sun
    #
    # Z/Z_sun = (Z_m / X_h)_cell
    #           ----------------
    #           (Z_m / X_h)_sun
    #
    # where X_h + Y_he + Z_m = 1
    # Since the simulation does not track Helium, need to assume 
    # a value for r = Y/X
    

    # Solar values
    # Taken from Chris' Notes 
    X_sun = 0.70683
    Y_sun = 0.27431
    Z_sun = 0.0188
    
    r = 0.3347

    # Loop through cells
    X_cell = (1 - Z_cell) / (1 + r)
    Z = (Z_cell / X_cell) / (Z_sun / X_sun)

    return Z

