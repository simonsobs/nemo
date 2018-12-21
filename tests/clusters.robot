*** Test Cases ***

Recover published ACTPol cluster masses
    Cross match             ../examples/equD56/E-D56Clusters.fits   configs/equD56_quick/equD56_quick_M500.fits
    Check recovered ratio   M500cUPP    M500    tolerance=0.03    SNRCut=4.0  plotFileName=plots/M500Recovery_E-D56.png
    Status should be        SUCCESS

Recover unbiased masses from mock catalog
    Set config              configs/equD56_quick.yml
    Run nemo selfn
    Run nemo mock
    Run nemo mass           configs/equD56_quick/mocks/mockCatalog_1.fits
    Cross match             configs/equD56_quick/mocks/mockCatalog_1.fits   configs/equD56_quick/mocks/mockCatalog_1_M500.fits
    Check recovered ratio   true_M500   M500    tolerance=0.02    SNRCut=5.0  plotFileName=plots/M500Recovery_mock.png
    Status should be        SUCCESS
    
*** Keywords ***

Run ACTPol clusters pipeline
    Set config              configs/equD56_quick.yml
    Run nemo
    Run nemo mass

Clean up
    Remove directory        configs/equD56_quick    True

*** Settings ***
#Documentation              To be added here
Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Run ACTPol clusters pipeline
Suite Teardown              Clean up

*** Variables ***
