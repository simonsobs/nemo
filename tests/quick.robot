*** Test Cases ***

Quickstart clusters tutorial runs
    Run quickstart clusters
    Cross match             testsCache/quickstart-clusters/quickstart-clusters_optimalCatalog.fits      testsCache/DR5_cluster-catalog_v1.1.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    toleranceSigma=3.0    expectedRatio=0.94  errInKey=fixed_err_y_c  errOutKey=fixed_err_y_c  SNRKey=fixed_SNR  SNRCut=5.0  plotLabel=quickstart-clusters
    Status should be        SUCCESS

Quickstart sources tutorial runs
    Run quickstart sources

Mass estimation works
    Run quickstart clusters with Q
    Run nemo mass

Forced photometry using nemo works
    Set config              configs/quickstart-clusters-Q.yml
    Run nemo    DR5_cluster-catalog_v1.1.fits

Forced photometry using nemoMass works
    Run quickstart clusters with Q
    Run nemo mass   DR5_cluster-catalog_v1.1.fits   True

Source injection test runs
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo injection test

Quickstart multipass config runs
    Run quickstart multipass
    Cross match             testsCache/quickstart-multipass/quickstart-multipass_optimalCatalog.fits      testsCache/DR5_cluster-catalog_v1.1.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    toleranceSigma=3.0    expectedRatio=0.92  errInKey=fixed_err_y_c  errOutKey=fixed_err_y_c  SNRKey=fixed_SNR  SNRCut=5.0  plotLabel=quickstart-multipass
    Status should be        SUCCESS

Cluster sim with nemoModel runs
    Generate simulated cluster maps

Source sim with nemoModel runs
    Generate simulated source maps
    
Recovered sim source amplitudes are unbiased
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Cross match     testsCache/sim_f090_inputCatalog.fits     testsCache/sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits
    Check recovered ratio   deltaT_c    deltaT_c    toleranceSigma=3.0  SNRKey=SNR  SNRCut=7.0  plotLabel=sim_ptsrc_amplitudes
    Status should be        SUCCESS
    
End-to-end source recovery and subtraction
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Make signal only sim    sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits   f090    small
    Subtract maps   testsCache/sim_f090.fits  testsCache/signal_model_only_f090.fits  testsCache/sources_diff.fits
    Subtract maps   testsCache/signal_free_f090.fits  testsCache/sources_diff.fits  testsCache/signal_free-sources_diff.fits
    Check map sigma     testsCache/signal_free-sources_diff.fits    100
    Status should be        SUCCESS

End-to-end A10 cluster modeling and subtraction
    End-to-end cluster modeling and subtraction     A10

End-to-end B12 cluster modeling and subtraction
    End-to-end cluster modeling and subtraction     B12
    
    
*** Keywords ***

Run quickstart clusters
    Set config              ../examples/quickstart/quickstart-clusters.yml
    Run nemo

Run quickstart clusters with Q
    Set config              configs/quickstart-clusters-Q.yml
    Run nemo

Run quickstart sources
    Set config              ../examples/quickstart/quickstart-sources.yml
    Run nemo

Run quickstart multipass
    Set config              configs/quickstart-multipass.yml
    Run nemo

Generate simulated cluster maps
    Setup quickstart
    Make sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large
    Make sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large

Generate simulated cluster maps in parallel
    [Arguments]   ${profile}
    Setup quickstart
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large  ${profile}
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large  ${profile}
    Make signal free sim    0.0   1234  f150    large
    Make signal free sim    0.0   1234  f090    large
    
Generate simulated source maps
    Setup quickstart
    Make sim    pointsources-100   100.0    1234    f150    small
    Make sim    pointsources-100   100.0    1234    f090    small
    Make signal free sim    0.0   1234  f090    small

Generate large simulated source maps
    Setup quickstart
    Make sim    pointsources-100   100.0    1234    f150    large
    Make sim    pointsources-100   100.0    1234    f090    large
    Make signal free sim    0.0   1234  f090    large

End-to-end cluster modeling and subtraction
    [Arguments]   ${profile}
    Generate simulated cluster maps in parallel     ${profile}
    Set config      configs/sim_cl_${profile}_MFMF_tiles.yml
    Run parallel nemo
    Cross match     testsCache/DR5_cluster-catalog_v1.1.fits    testsCache/sim_cl_${profile}_MFMF_tiles/sim_cl_${profile}_MFMF_tiles_optimalCatalog.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    toleranceSigma=3.0  errInKey=fixed_err_y_c  errOutKey=fixed_err_y_c   SNRKey=fixed_SNR  SNRCut=5.0  plotLabel=${profile}_sim_cl_amplitudes
    Status should be        SUCCESS
    Make parallel signal only sim    sim_cl_${profile}_MFMF_tiles/sim_cl_${profile}_MFMF_tiles_optimalCatalog.fits  f090    large   ${profile}
    Subtract maps   testsCache/sim_f090.fits  testsCache/signal_model_only_f090.fits  testsCache/cl_${profile}_f090_diff.fits
    Subtract maps   testsCache/signal_free_f090.fits  testsCache/cl_${profile}_f090_diff.fits  testsCache/signal_free-cl_${profile}_f090_diff.fits
    Check map sigma     testsCache/signal_free-cl_${profile}_f090_diff.fits    50
    Status should be        SUCCESS
    [Teardown]  Clean up cluster sim    ${profile}

Clean up cluster sim
    [Arguments]   ${profile}
    Remove directory        testsCache/sim_cl_${profile}_MFMF_tiles     True

Clean up
    Remove directory        testsCache/sim_cl_MFMF                      True
    Remove directory        testsCache/sim_ptsrc_f090                   True
    Remove directory        testsCache/quickstart-clusters              True
    Remove directory        testsCache/quickstart-sources               True
    Remove directory        testsCache/quickstart-multipass             True
    Remove directory        testsCache/quickstart-multipass-Q           True
    Remove file             testsCache/*forcedCatalog*
    Remove file             testsCache/*_mass.fits

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Setup quickstart
Suite Teardown              Clean up


*** Variables ***
