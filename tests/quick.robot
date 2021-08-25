*** Test Cases ***

Quickstart clusters
    Run quickstart clusters

Quickstart sources
    Run quickstart sources
    
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

Run ACTPol clusters pipeline - real space filter
    Set config              configs/equD56_quick.yml
    Run nemo
    Run nemo mass

Run ACTPol clusters pipeline - Fourier space filter
    Set config              configs/equD56_quick_fourier.yml
    Run nemo
    Run nemo mass

Clean up
    Remove directory        testsCache/quickstart-clusters  True
    Remove directory        testsCache/quickstart-sources   True

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Setup quickstart
Suite Teardown              Clean up

*** Variables ***
