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

Nemo runs in parallel
    Generate simulated A10 cluster maps in parallel
    Set config      configs/sim_cl_MFMF_tiles.yml
    Run parallel nemo
    
Recovered sim source amplitudes are unbiased
    Generate simulated source maps
    Set config      configs/sim_ptsrc_f090.yml
    Run nemo
    Cross match     testsCache/sim_f090_inputCatalog.fits     testsCache/sim_ptsrc_f090/sim_ptsrc_f090_optimalCatalog.fits
    Check recovered ratio   deltaT_c    deltaT_c    tolerance=0.02  SNRKey=SNR  SNRCut=5.0  plotFileName=plots/sim_ptsrc_amplitudes.png
    Status should be        SUCCESS

Recovered sim cluster amplitudes are unbiased
    Generate simulated cluster maps
    Set config      configs/sim_cl_MFMF.yml
    Run nemo
    Set config      configs/sim_cl_MFMF_pass2.yml
    Run nemo
    Cross match     testsCache/DR5_cluster-catalog_v1.1.fits    testsCache/sim_cl_MFMF/sim_cl_MFMF_optimalCatalog.fits
    Check recovered ratio   fixed_y_c    fixed_y_c    tolerance=0.02  SNRKey=fixed_SNR  SNRCut=5.0  plotFileName=plots/sim_cl_amplitudes.png
    Status should be        SUCCESS

#Recovered sim cluster amplitudes are unbiased - two pass
#    Generate simulated cluster maps
#    Set config      configs/sim_cl_MFMF.yml
#    Run nemo
#    Set config      configs/sim_cl_MFMF_pass2.yml
#    Run nemo
#    Cross match     testsCache/DR5_cluster-catalog_v1.1.fits    testsCache/sim_cl_MFMF_pass2/sim_cl_MFMF_pass2_optimalCatalog.fits
#    Check recovered ratio   fixed_y_c    fixed_y_c    tolerance=0.02  SNRKey=fixed_SNR  SNRCut=5.0  plotFileName=plots/sim_cl_amplitudes_pass2.png
#    Status should be        SUCCESS
    
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
    Generate simulated cluster maps in parallel     A10
    Set config      configs/sim_cl_MFMF_tiles.yml
    Run parallel nemo
    Make parallel signal only sim    sim_cl_MFMF_tiles/sim_cl_MFMF_tiles_optimalCatalog.fits   f090    large
    Subtract maps   testsCache/sim_f090.fits  testsCache/signal_model_only_f090.fits  testsCache/cl_f090_A10_diff.fits
    Subtract maps   testsCache/signal_free_f090.fits  testsCache/cl_f090_A10_diff.fits  testsCache/signal_free-cl_f090_A10_diff.fits
    Check map sigma     testsCache/signal_free-cl_f090_A10_diff.fits    50
    Status should be        SUCCESS

#    Make parallel sim signal only sim    sim_cl_MFMF/sim_cl_MFMF_optimalCatalog.fits   f150    large

    
*** Keywords ***

Run quickstart clusters
    Set config              ../examples/quickstart/quickstart-clusters.yml
    Run nemo

Run quickstart sources
    Set config              ../examples/quickstart/quickstart-sources.yml
    Run nemo

Generate simulated cluster maps
    Setup quickstart
    Make sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large
    Make sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large

Generate simulated cluster maps in parallel
    Arguments   ${profile}
    Setup quickstart
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large  ${profile}
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large  ${profile}
    Make signal free sim    0.0   1234  f150    large
    Make signal free sim    0.0   1234  f090    large

Generate simulated A10 cluster maps in parallel
    Setup quickstart
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large  A10
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large  A10

Generate simulated B12 cluster maps in parallel
    Setup quickstart
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f150    large  B12
    Make parallel sim    DR5_cluster-catalog_v1.1.fits   50.0    1234    f090    large  B12
    
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
    
Clean up
    Remove directory        testsCache/sim_cl_MFMF  True
    Remove directory        testsCache/tileDir_auto_1.0_surveyMask  True
    Remove directory        testsCache/tileDir_auto_1.0_sim_f090.fits   True
    Remove directory        testsCache/tileDir_auto_1.0_sim_f150.fits   True
    Remove directory        testsCache/sim_cl_MFMF_tiles  True
    #Remove directory        testsCache/sim_cl_MFMF_pass2  True
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
