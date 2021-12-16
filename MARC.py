# This package is converted from a Matlab script provided by Mukesh Yadav
# Translated and expanded upon by Espen Hovland, 2021,
# as part of my specialization project at NTNU.

import numpy as np
from fractions import Fraction

class Ring:
    """
    Add-drop ring resonator class. Calculates all relevant data when an instance is created.

    Parameters:
        wavelengths          (np.array): Wavelength sweep
        radius               (float):    Ring radius
        angular_separation   (float):    Angular separation of drop- and through-port waveguides
        coupling_coefficient (float):    Self-coupling coefficient of the input waveguide
        loss_coefficient     (float):    Round-trip loss coefficient of the ring
        n_eff                (float):    Effective refractive index of the ring waveguide

    Static methods:
        FSR(ring_radius, n_eff, lambda_0) -> FSR in [m]

    Available data (member variables):
        .r                   (float):    Radius of ring
        .ang_sep             (float):    Angular sep. of through- and drop-port
        .a                   (float):    Round-trip loss coefficient of the ring
        .n_eff               (float):    Effective refractive index
        .t1                  (float):    Self-coupling coefficient of input waveguide
        .t2                  (float):    Self-coupling coefficient of input waveguide
        .fsr                 (float):    Free spectral range
        .fsr_eff             (float):    Effective FSR (single-MARC)
        .dp_amplitude        (np.array): Drop-port amplitude response
        .tp_amplitude        (np.array): Through-port amplitude response
        .dp_intensity        (np.array): Drop-port intensity response
        .tp_intensity        (np.array): Through-port intensity response
        .dp_phase            (np.array): Drop-port phase response
        .tp_phase            (np.array): Through-port phase response
    """
    def __init__(self,
                 wavelengths,
                 radius:                float,
                 angular_separation:    float,
                 coupling_coefficient:  float = 0.95,
                 loss_coefficient:      float = 1, 
                 n_eff:                 float = 3.9 ) -> None:
        self.r       = radius               # [m]   Radius of ring
        self.ang_sep = angular_separation   # [deg] Angular sep. of drop- and through-port
        self.a       = loss_coefficient     # [1]   Round-trip loss coefficient
        self.n_eff   = n_eff                # [1]   Effective refractive index
        self.t1      = coupling_coefficient # [1]   Self-coupling coeff. of input waveguide
        self.t2      = self.t1 / self.a     # [1]   Self-coupling coeff. of drop waveguide

        # Round-trip phase shifts
        phi = 2 * np.pi * self.n_eff * 2 * np.pi * self.r / wavelengths 

        traversed = self.ang_sep / 360  # Proportion of ring travelled
        act_shift = phi * traversed     # Actual phase shift

        # Amplitude response at drop-port
        self.dp_amplitude = np.asarray(
            (-(np.sqrt(1-self.t1**2) * np.sqrt(1-self.t2**2) * np.power(self.a, traversed) 
            * np.exp(1j*act_shift)) / ( 1 - (self.t1 * self.t2 * self.a * np.exp(1j * phi)))), 
            dtype=complex
        )
        # Amplitude response at through-port
        self.tp_amplitude = np.asarray(
              (self.t1 - self.t2 * self.a * np.exp(1j * phi))
            / (1 - self.t1 * self.t2 * self.a * np.exp(1j * phi)),
            dtype=complex
        )

        # Intensity response
        self.dp_intensity = np.asarray(np.abs(self.dp_amplitude)**2, dtype=float) # Drop-port
        self.tp_intensity = np.asarray(np.abs(self.tp_amplitude)**2, dtype=float) # Through-port

        # Phase response
        self.dp_phase = np.unwrap(np.angle(self.dp_amplitude)) # Drop-port
        self.tp_phase = np.unwrap(np.angle(self.tp_amplitude)) # Through-port

        # Get the FSR, as if the center frequency is resonant frequency (approximation)
        self.fsr = self.get_FSR(self.r, self.n_eff, wavelengths[len(wavelengths)//2])

        # Calculate effective FSR (for single-ring MARCs)
        L    = 1 / traversed
        frac = Fraction.from_float(L).limit_denominator() # Find the reduced fraction N/M
        N    = frac.numerator
        self.fsr_eff: float = N * self.fsr # [nm]

        return

    @staticmethod
    def get_FSR(ring_radius: float, n_eff: float, lambda_0: float = 1550e-9) -> float:
        """
        Calculate the (approximate) FSR.
        Parameters:
            ring_radius (float): [m] Radius of ring
            n_eff       (float): [1] Effective refractive index of ring
            lambda_0    (float): [m] Resonance wavelength of ring. 
                                     Defaults to 1550 nm
        Returns:
            Calculated FSR [m]
        """
        return lambda_0**2 / (n_eff*2*np.pi*ring_radius)


#######################################################################################################
#######################################################################################################


class MZI:
    """ Mach-Zehnder interferometer base class """
    def __init__(self, initial_phase: float = 0) -> None:
        self.initial_phase = initial_phase # Phase imbalance of MZI

    @staticmethod
    def get_trans(amplitude, phi_1, phi_2) -> np.array:
        """
        The interference equation, for calculating the amplitude transmittance.

        Parameters:
            amplitude (np.array): The complex amplitude of the signal from the affected arm.
            phi_1     (np.array): The phase of the affected arm.
            phi_2     (np.array): The phase of the unaffected arm.

        Returns:
            t         (np.array): The transmittance of the MZ interferometer.
        """
        t = 0.5*np.exp(0.5j*(phi_1+phi_2+np.pi))*(np.abs(amplitude)*np.exp(0.5j*(phi_1-phi_2))
          + np.exp(-0.5j*(phi_1-phi_2)))
        return t

#######################################################################################################
#######################################################################################################


class MARC(MZI):
    """
    Multi-ring MARC class. Calculates all data upon creating an instance.

    Parameters:
        wavelengths   (np.array): Wavelength sweep
        rings         (Ring):     Rings to include in the MARC
        initial_phase (float):    Phase imbalance from MZI

    Available data (member variables):
        .num_rings    (int):      number of rings in MARC
        .rings        (list):     list of rings in device
        .amplitude    (np.array): output amplitude response of device
        .intensity    (np.array): output intensity response of device
        .tp_intensity (np.array): through-port intensity response of device
        .phase_output (np.array): output phase response of device
    """
    def __init__(self, wavelengths, *rings: Ring, initial_phase: float = 0) -> None:
        super().__init__(initial_phase=initial_phase)

        self.num_rings  = len(rings)
        self.rings      = [r for r in rings]

        amp_post_rings  = np.zeros(len(wavelengths), dtype="complex") # Amplitude resp. of all rings
        amp_tp          = np.ones(len(wavelengths), dtype="complex")  # Amplitude resp. of through-port

        for idx, ring in enumerate(self.rings):
            amp_ring = ring.dp_amplitude # Drop-port amplitude transmission of last ring
            for i in range(0, idx):
                amp_ring *= self.rings[i].tp_amplitude
            amp_post_rings += amp_ring

        for ring in self.rings:
            amp_tp *= ring.tp_amplitude

        amplitude = self.get_trans(amp_post_rings, np.unwrap(np.angle(amp_post_rings)), initial_phase)

        self.amplitude     = amplitude / np.max(np.abs(amplitude)) # Normalize the amplitude
        self.intensity     = np.asarray(np.abs(self.amplitude)**2, dtype=float)
        self.tp_intensity  = np.asarray(np.abs(amp_tp)**2, dtype=float)
        self.phase_output  = np.unwrap(np.angle(self.amplitude))

        return

#######################################################################################################
#######################################################################################################