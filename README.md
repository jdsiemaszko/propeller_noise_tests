A validation dataset for elements of propeller noise prediction models. The case mimics the experimental setup of [1] in the configuration D20L20 and covers the modelling elements in [2]. The two tests cover the potential interaction model for loading prediction (see [2,3]) and the propagation of sources to the acoustic far field via the acoustic analogy (see [4,5]):

1) Source Prediction via Potential Interaction:
    The test file contains model inputs: geometry, mean loading, and mean inflow velocity distribution
    The expected output is an array of blade loading harmonics (see equations 3 and 4 in [1]), strut loading harmonics (see equation 18), and a downwash profile along the rotation (equations 7-8).

2) Propagation via Acoustic Analogy:
    The test contains inputs same as above, additionally including the blade and strut loading harmonics and mean thickness-to-chord of the blade.
    The expected outputs are far-field acoustic pressure spectra for several radiating components:
        -rotor steady loading noise (eq. 2)
        -rotor unsteady loading noise (eq. 2)
        -rotor thickness noise
        -strut loading noise (eq. 9)
    The test assumes the "standard" implementation of the Hanson & Patrzych propagator model, rotating counterclockwise at the frequency Omega / 2 / pi measured in Hz.
    The rotor is assumed positioned at (0, 0, 0) in 3D, with the axis oriented along (0, 0, 1) with positive z taken as upstream. The zero-azimuth datum is (1, 0, 0). 
    
In the current dataset, the variables in the frequency domain are defined with the Fourier transform convention: \hat{f}_k = \Omega / 2 / \pi \int_0^{2 \pi / \Omega} (f(t) * exp(1j * k * \Omega * t)) dt. Mind the sign of the exponential term.

File contents:

    test_loading_prediction.h5
    │
    ├── inputs
    │   ├── Omega_rad_p_s
    │   ├── rho_kg_p_m3
    │   ├── sos_m_p_s
    │   ├── radius_inner_m
    │   ├── radius_outer_m
    │   ├── twist_inner_rad
    │   ├── twist_outer_rad
    │   ├── chord_inner_m
    │   ├── chord_outer_m
    │   ├── L_cylinder_m
    │   ├── D_cylinder_m
    │   ├── B
    │   ├── F_z_prime_N_p_m
    │   ├── F_phi_prime_N_p_m
    │   ├── F_r_prime_N_p_m
    │   └── U_inf_m_p_s
    │
    └── outputs
        ├── blade_loading_harmonics_N_p_m
        │   ├── real
        │   └── imag
        │
        ├── strut_loading_harmonics_N_p_m
        │   ├── real
        │   └── imag
        │
        ├── azimuth_rad
        ├── harmonics_k
        ├── frequency_k_Hz
        └── blade_downwash_m_p_s



    test_propagator.h5
    │
    ├── inputs
    │   ├── Omega_rad_p_s
    │   ├── rho_kg_p_m3
    │   ├── sos_m_p_s
    │   ├── radius_inner_m
    │   ├── radius_outer_m
    │   ├── twist_inner_rad
    │   ├── twist_outer_rad
    │   ├── chord_inner_m
    │   ├── chord_outer_m
    │   ├── t_c_inner
    │   ├── t_c_outer
    │   ├── L_cylinder_m
    │   ├── D_cylinder_m
    │   │
    │   ├── blade_loading_harmonics_N_p_m
    │   │   ├── real
    │   │   └── imag
    │   │
    │   ├── strut_loading_harmonics_N_p_m
    │   │   ├── real
    │   │   └── imag
    │   │
    │   └── B
    │
    └── outputs
        ├── p_loading_blade_steady_Pa
        │   ├── real
        │   └── imag
        │
        ├── p_loading_blade_unsteady_Pa
        │   ├── real
        │   └── imag
        │
        ├── p_loading_strut_Pa
        │   ├── real
        │   └── imag 
        │
        ├── p_thickness_blade_Pa
        │   ├── real
        │   └── imag
        │
        ├── observer_positions_m
        ├── harmonics_m
        ├── frequency_mB_Hz
        │
        ├── observer_radius_m
        ├── observer_polar_rad
        └── observer_azimuth_rad

Description of Variables:

inputs:

    Omega_rad_p_s : float - rotor rotational frequency in rad/s.
    rho_kg_p_m3 : float - ambient air density in kg/m^3.
    sos_m_p_s : float - ambient speed of sound in m/s.

    radius_inner_m : array of float of shape Nr - radial stations used in the discretization, in units m. Values represent the attachment point of the loading/thickness sources in part 2) and the radial stations assimilated to the 2D problem in part 1).

    radius_outer_m : array of float of shape Nr+1 - radial stations of the element edges in the discretization in units m. Related to radius_inner_m as: radius_inner_m = 1/2 * (radius_outer_m[1:] + radius_outer_m[:-1]). 

    chord_inner_m : array of float of shape Nr - chord length associated with radial stations radius_inner_m in units m.

    twist_inner_rad : array of float of shape Nr - twist of the propeller chordline w.r.t the propeller plane in units rad, associated with radial stations radius_inner_m.

    t_c_inner : array of float of shape Nr - mean thickness-to-chord ratio at the the radial station radius_inner_m, defined as: 1/c^2\int_0^c t(x)dx for t(x) representing the local thickness of the airfoil along the chord.

    chord_outer_m, twist_outer_rad, t_c_outer : arrays of float of shape Nr+1 - corresponding values of chord, twist, and t/c associated with radial stations radius_outer_m.

    L_cylinder_m : float - spacing between the strut and propeller plane, measured from the strut center in units m. See Fig. 6 in [2].

    D_cylinder_m : float - strut radius in units m.

    B : int - number of blades of the propeller.

    F_z_prime_N_p_m : array of float of shape Nr - distributed mean loading over the propeller blades in the axial direction, in units Newton per meter. Sign convention: positive upstream.

    F_phi_prime_N_p_m : array of float of shape Nr - distributed mean loading over the propeller blades in the azimuthal direction, in units Newton per meter. *Sign convention: positive opposite to the direction of rotation (drag is positive)*.

    F_r_prime_N_p_m : array of float of shape Nr - distributed mean loading over the propeller blades in the radial direction, in units Newton per meter. Set to zero in the current dataset. Sign convention: positive outwards.

    U_inf_m_p_s : array of float of shape (2, Nr) - Inflow velocity in the 2D potential flow problem, assumed unchanged between the propeller and strut stations, corresponding to the radial stations radius_inner_m, in units m/s. The sign convention is positive to the right (in the rotation direction), positive upwards (upstream), see Fig. 6 in [2]. x-component set to zero in the current dataset.

outputs:

    blade_loading_harmonics_N_p_m/real, imag : arrays of float of shape (3, Nk, Nr) - loading harmonics acting on a propeller blade, measured in N/m and separated between real and imaginary components. Axis 0 corresponds to the three force components, in order: radial (positive outwards), axial (positive upwards), azimuthal (positive opposite to the direction of travel). Axis 1 corresponds to the resolved harmonics of rotational frequency: k * Omega / 2 / pi measured in Hz. Harmonics k are provided in entry harmonics_k; element 0 along axis 1 corresponds to mean loading defined via F_z_prime_N_p_m, F_phi_prime_N_p_m, F_r_prime_N_p_m. Axis 2 corresponds to radial stations defined in radius_inner_m. Total loading harmonics are recovered as real + 1j * imag.

    strut_loading_harmonics_N_p_m/real, imag : arrays of float of shape (3, Nk, Nr) - loading harmonics over the strut; same conventions as above. In the current set, all odd harmonics are zero.

    harmonics_k : array of int of shape Nk - integers k defining the harmonics of the rotational frequency Omega corresponding to the loading harmonics blade_loading_harmonics_N_p_m/real, imag and strut_loading_harmonics_N_p_m/real, imag. the corresponding frequency in Hz is: f = k * Omega / 2 / pi, provided in frequency_k_Hz. 

    frequency_k_Hz : array of float of shape Nk - frequency associated with the loading harmonics defined as: k * Omega / 2 / pi, measured in Hz.

    blade_downwash_m_p_s : array of shape (Nr, Nt) - downwash acting on the propeller blade normal to the chord, in the time domain, measured in m/s. Axis 0 corresponds to radial stations radius_inner_m. Axis 1 corresponds to azimuth stations: phi = Omega * t defined in azimuth_rad. The time datum t=0 corresponds to the strut azimuth station.

    blade_downwash_harmonics_m_p_s : array of shape (Nk, Nr) - harmonics of downwash defined in blade_downwash_m_p_s computed with the same Fourier transform convention as defined in the header. Axis 1 corresponds to resolved harmonics in harmonics_k, axis 2 corresponds to radial stations radius_inner_m.

    azimuth_rad : array of shape Nt - azimuth stations phi = Omega * t corresponding to blade_downwash_m_p_s in radians.

    p_loading_blade_steady_Pa/real, imag : arrays of float of shape (Nx, Nm) - far-field pressure field of the propeller steady loading noise in units Pa at locations defined by observer_position_m and harmonics m of frequency Omega*B/2/pi defined by harmonics_m. Axis 0 corresponds to observer positions, axis 1 to harmonic number.

    p_loading_blade_unsteady_Pa/real, imag : arrays of float of shape (Nx, Nm) - far-field pressure field of the propeller unsteady loading noise in units Pa. Convention as above.

    p_loading_strut_Pa/real, imag : arrays of float of shape (Nx, Nm) - far-field pressure field of the strut loading noise in units Pa. Convention as above.

    p_thickness_blade_Pa/real, imag : arrays of float of shape (Nx, Nm) - far-field pressure field of the propeller thickness noise in units Pa. Convention as above.

    observer_position_m: array of floats of shape (3, Nx) - observer positions in the far field, in units m. Axis 0 corresponds to the 3 vector components x, y, z. Propeller is positioned at (0, 0, 0) with its axis along (0, 0, 1). The zero azimuth station is taken in the direction (1, 0, 0), corresponding to the strut azimuth.

    observer_radius_m, observer_polar_rad, observer_azimuth_rad : arrays of float of shape Nx - observer locations in spherical coordinates in the propeller frame, consistent with the model in [3]. Polar angle defined w.r.t to the propeller axis, increasing in the downstream direction. Azimuth defined counterclockwise (in the direction of rotation) w.r.t. the direction (1, 0, 0).

    harmonics_m : array of int of shape Nm - integers m defining the harmonics of the blade-passing frequency Omega*B corresponding to the far-field pressure harmonics p_loading_blade_steady_Pa, etc. The corresponding frequency in Hz is: m * B * Omega / 2 /pi, provided in frequency_mB_Hz.

    frequency_mB_Hz : array of float of shape Nm - frequency associated with the far-field pressure harmonics m defined as: m * B * Omega / 2 / pi, measured in Hz.

Bibliography

[1] Gojon, R., Parisot-Dupuis, H., Mellot, B., and Jardin, T., “Aeroacoustic radiation of low Reynolds number rotors in interaction with beams,” The Journal of the Acoustical Society of America, Vol. 154, No. 2, 8 2023, pp. 1248–1260.

[2]  Vella, E., Gojon, R., Parisot-Dupuis, H., Doué, N., Jardin, T., and Roger, M., “Mutual Interaction Noise in Rotor–strut Configuration,” AIAA Journal, Vol. 0, No. 0, 0, pp. 1–16.

[3] Wu, Y., Kingan, M. J., and Go, S. T., “Propeller-strut interaction tone noise,” Physics of Fluids, Vol. 34, No. 5, 5 2022.

[4] Hanson, D. B. and Parzych, D. J., “Theory for Noise of Propellers in Angular Inflow With Parametric Studies and Experimental Verification,” Tech. Rep. NASA Contractor Report 4499, NASA, 1993.

[5] Lowson, M. V., “Theoretical Analysis of Compressor Noise,” Tech. Rep. 1B, 01 1970.
