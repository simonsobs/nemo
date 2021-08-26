*** Test Cases ***

Quickstart clusters
    Run quickstart clusters

Quickstart sources
    Run quickstart sources

Make cluster sims
    Generate simulated cluster maps

Make source sims
    Generate simulated source maps
    
Recovered sim source amplitudes are unbiased
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Cross match     testsCache/smallTest_f090_inputCatalog.fits     testsCache/sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits
    Check recovered ratio   deltaT_c    deltaT_c    tolerance=0.01  SNRKey=SNR  SNRCut=5.0  plotFileName=plots/sim_sources_amplitudes.png

    
*** Keywords ***

Run quickstart clusters
    Set config              ../examples/quickstart/quickstart-clusters.yml
    Run nemo
    Cross match             testsCache/quickstart-clusters/quickstart-clusters_optimalCatalog.fits      testsCache/DR5_cluster-catalog_v1.1.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    tolerance=0.01  SNRKey=fixed_SNR  SNRCut=5.0  plotFileName=plots/quickstart-clusters.png
    Status should be        SUCCESS

Run quickstart sources
    Set config              ../examples/quickstart/quickstart-sources.yml
    Run nemo

Generate simulated cluster maps
    Setup quickstart
    Make sim f150   DR5_cluster-catalog_v1.1.fits   100.0    1234
    Make sim f090   DR5_cluster-catalog_v1.1.fits   100.0    1234

Generate simulated source maps
    Setup quickstart
    Make sim f150   pointsources-100   100.0    1234
    Make sim f090   pointsources-100   100.0    1234
    
Clean up
    Remove directory        testsCache/sim_ptsrc_f090  True
    Remove directory        testsCache/quickstart-clusters  True
    Remove directory        testsCache/quickstart-sources   True

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Setup quickstart
Suite Teardown              Clean up

*** Variables ***
