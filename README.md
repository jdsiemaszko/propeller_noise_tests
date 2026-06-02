This repository summarizes validation data for testing propeller noise prediction models covered in [1]. A total of 5 (2+3) tests are included for the source prediction and propagation.

1) Predictor:
    1.1) test_PIN_blade_loadings.h5

        inputs: 
        geometry parameters
        mean blade loading 
        radial inflow distribution

        outputs:
        blade loading harmonics

    1.2) test_PIN_beam_loadings.h5

        inputs: 
            geometry parameters
            mean blade loading 
            radial inflow distribution

        outputs:
            blade loading harmonics  

2) Propagator:
    2.1) test_blade_loading_noise.h5

        inputs:
            operating point
            ...
            blade loading harmonics...

        outputs:
            steady loading noise
            unsteady loading noise

    2.2) test_blade_thickness_noise.h5

        inputs:
            operating point
            ...
            blade thickness distribution (NACA 0012)

        outputs:
            blade thickness noise

    2.3) test_beam_loading_noise.h5

        inputs:

            beam loading harmonics

        outputs:

            beam loading noise