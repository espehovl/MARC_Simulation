# Minimal working example of use of the MARC class

# From the MARC file, we need the MARC and Ring classes
from MARC import MARC, Ring

# The MARC and Ring objects require the wavelength range input as numpy arrays
import numpy as np

# For simpler plotting
import matplotlib.pyplot as plt

# The wavelength range we are interested in
wavelengths = np.arange(1530e-9, 1560e-9, 1e-13)

# Create each add-drop ring you want in your MARC
# => Ring(<wavelength range>, <ring radius>, <angular separation>)
ring1 = Ring(wavelengths, 30e-6, 90)
ring2 = Ring(wavelengths, 45e-6, 135)
# ring3 = ...

# Create the MARC by letting the class know the wavelength range and the constituent rings
# All transmission and phase responses are calculated upon creating the MARC object
sensor = MARC(wavelengths, ring1, ring2)

plt.plot(wavelengths*1e9, sensor.intensity)
plt.xlabel("Wavelength, [nm]")
plt.ylabel("Intensity transmission")
plt.show()