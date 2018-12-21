*** Test Cases ***
    
Recovered source amplitudes are unbiased
    Cross match             inputSourcesCatalog.fits    configs/PSTest_E-D56/PSTest_E-D56_optimalCatalog.fits
    Check recovered ratio   I   deltaT_c    tolerance=0.01    SNRKey=SNR    SNRCut=5.0  plotFileName=plots/amplitudeRecovery.png
    Status should be        SUCCESS

Recovered source positions are accurate
    Cross match             inputSourcesCatalog.fits    configs/PSTest_E-D56/PSTest_E-D56_optimalCatalog.fits
    Check recovered positions   toleranceArcsec=12.0    SNRMax=10.0     plotFileName=plots/positionRecovery.png
    Status should be        SUCCESS
    
*** Keywords ***

Find injected sources
    Set config              configs/PSTest_E-D56.yml
    Inject sources
    Run nemo

Clean up
    Remove file             sourceInjectedMap.fits
    Remove file             inputSourcesCatalog.fits
    Remove directory        configs/PSTest_E-D56    True

*** Settings ***

Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Find injected sources
Suite Teardown              Clean up

*** Variables ***

