*** Test Cases ***

Recover published ACTPol cluster masses - real space filter
    Run ACTPol clusters pipeline - real space filter
    Cross match             testsCache/E-D56Clusters.fits   testsCache/equD56_quick/equD56_quick_mass.fits
    Check recovered ratio   M500cUPP    M500c   tolerance=0.03    SNRCut=4.0  plotFileName=plots/M500Recovery_E-D56.png
    Status should be        SUCCESS

Recover unbiased masses from mock catalog - real space filter
    Run ACTPol clusters pipeline - real space filter
    Run nemo selfn          ../examples/ACT-DR3-clusters/equD56_quick.yml
    Run nemo mock           testsCache/equD56_quick/selFn  testsCache/equD56_quick/mocks
    Run nemo mass           testsCache/equD56_quick/mocks/mockCatalog_1.fits
    Cross match             testsCache/equD56_quick/mocks/mockCatalog_1.fits   testsCache/equD56_quick/mocks/mockCatalog_1_mass.fits
    Check recovered ratio   true_M500c   M500c  tolerance=0.02    SNRCut=5.0  plotFileName=plots/M500Recovery_mock.png
    Status should be        SUCCESS

Recover published ACTPol cluster masses - Fourier space filter
    Run ACTPol clusters pipeline - Fourier space filter
    Cross match             testsCache/E-D56Clusters.fits   testsCache/equD56_quick_fourier/equD56_quick_fourier_mass.fits
    Check recovered ratio   M500cUPP    M500c   tolerance=0.03    SNRCut=4.0  plotFileName=plots/M500Recovery_E-D56_fourier.png
    Status should be        SUCCESS
    
*** Keywords ***

Run ACTPol clusters pipeline - real space filter
    Set config              ../examples/ACT-DR3-clusters/equD56_quick.yml
    Run nemo
    Run nemo mass

Run ACTPol clusters pipeline - Fourier space filter
    Set config              configs/equD56_quick_fourier.yml
    Run nemo
    Run nemo mass

Clean up
#    Remove directory        testsCache/equD56_quick    True
    Remove directory        testsCache/equD56_quick_fourier    True

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Setup equD56
Suite Teardown              Clean up

*** Variables ***
