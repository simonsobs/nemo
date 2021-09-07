*** Test Cases ***

Quickstart clusters tutorial runs
    Run quickstart clusters
    Cross match             testsCache/quickstart-clusters/quickstart-clusters_optimalCatalog.fits      testsCache/DR5_cluster-catalog_v1.1.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    tolerance=0.01  SNRKey=fixed_SNR  SNRCut=5.0  plotFileName=plots/quickstart-clusters.png
    Status should be        SUCCESS

Quickstart sources tutorial runs
    Run quickstart sources

Cluster sim with nemoModel runs
    Generate simulated cluster maps

Source sim with nemoModel runs
    Generate simulated source maps
    
Recovered sim source amplitudes are unbiased
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Cross match     testsCache/smallTest_f090_inputCatalog.fits     testsCache/sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits
    Check recovered ratio   deltaT_c    deltaT_c    tolerance=0.02  SNRKey=SNR  SNRCut=5.0  plotFileName=plots/sim_ptsrc_amplitudes.png
    Status should be        SUCCESS

Recovered sim cluster amplitudes are unbiased
    Generate simulated cluster maps
    Set config      configs/sim_cl_MFMF.yml
    Run nemo
    Cross match     testsCache/DR5_cluster-catalog_v1.1.fits    testsCache/sim_cl_MFMF/sim_cl_MFMF_optimalCatalog.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    tolerance=0.02  SNRKey=fixed_SNR  SNRCut=5.0  plotFileName=plots/sim_cl_amplitudes.png
    Status should be        SUCCESS
    
End-to-end source recovery and subtraction
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Make signal only sim f090   sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits
    Subtract maps   testsCache/smallTest_f090.fits  testsCache/signal_model_only_f090.fits  testsCache/sources_diff.fits
    Subtract maps   testsCache/signal_free_f090.fits  testsCache/sources_diff.fits  testsCache/signal_free-sources_diff.fits
    Check map sigma     testsCache/signal_free-sources_diff.fits    100
    Status should be        SUCCESS

    
*** Keywords ***

Run quickstart clusters
    Set config              ../examples/quickstart/quickstart-clusters.yml
    Run nemo

Run quickstart sources
    Set config              ../examples/quickstart/quickstart-sources.yml
    Run nemo

Generate simulated cluster maps
    Setup quickstart
    Make sim f150   DR5_cluster-catalog_v1.1.fits   50.0    1234
    Make sim f090   DR5_cluster-catalog_v1.1.fits   50.0    1234

Generate simulated source maps
    Setup quickstart
    Make sim f150   pointsources-100   100.0    1234
    Make sim f090   pointsources-100   100.0    1234
    Make signal free sim f090   0.0   1234
    
Clean up
    Remove directory        testsCache/sim_cl_MFMF  True
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
