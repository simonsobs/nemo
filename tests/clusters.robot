*** Test Cases ***

Recover published ACTPol cluster masses - real space filter
    Run ACTPol clusters pipeline - real space filter
    Cross match             ../examples/equD56/E-D56Clusters.fits   configs/equD56_quick/equD56_quick_M500.fits
    Check recovered ratio   M500cUPP    M500    tolerance=0.03    SNRCut=4.0  plotFileName=plots/M500Recovery_E-D56.png
    Status should be        SUCCESS

Recover unbiased masses from mock catalog - real space filter
    Set config              configs/equD56_quick.yml
    Run nemo selfn
    Run nemo mass
    Run nemo mock
    Run nemo mass           configs/equD56_quick/mocks/mockCatalog_1.fits
    Cross match             configs/equD56_quick/mocks/mockCatalog_1.fits   configs/equD56_quick/mocks/mockCatalog_1_M500.fits
    Check recovered ratio   true_M500   M500    tolerance=0.02    SNRCut=5.0  plotFileName=plots/M500Recovery_mock.png
    Status should be        SUCCESS

Recover published ACTPol cluster masses - Fourier space filter
    Run ACTPol clusters pipeline - Fourier space filter
    Cross match             ../examples/equD56/E-D56Clusters.fits   configs/equD56_quick_fourier/equD56_quick_fourier_M500.fits
    Check recovered ratio   M500cUPP    M500    tolerance=0.03    SNRCut=4.0  plotFileName=plots/M500Recovery_E-D56_fourier.png
    Status should be        SUCCESS
    
*** Keywords ***

Run ACTPol clusters pipeline - real space filter
    Set config              configs/equD56_quick.yml
    Run nemo
    Run nemo mass

Run ACTPol clusters pipeline - Fourier space filter
    Set config              configs/equD56_quick_fourier.yml
    Run nemo
    Run nemo mass

Clean up
    Remove directory        configs/equD56_quick    True
    Remove directory        configs/equD56_quick_fourier    True

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Setup equD56
Suite Teardown              Clean up

*** Variables ***
