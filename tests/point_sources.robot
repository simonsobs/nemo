*** Test Cases ***
    
Recovered source amplitudes are unbiased
    Cross match             testsCache/inputSourcesCatalog.fits    testsCache/PSTest_E-D56/PSTest_E-D56_optimalCatalog.fits
    Check recovered ratio   deltaT_c   deltaT_c    tolerance=0.01    SNRKey=SNR    SNRCut=5.0  plotFileName=plots/amplitudeRecovery.png
    Status should be        SUCCESS

Recovered source positions are accurate
    Cross match             testsCache/inputSourcesCatalog.fits    testsCache/PSTest_E-D56/PSTest_E-D56_optimalCatalog.fits
    Check recovered positions   toleranceArcsec=12.0    SNRMax=10.0     plotFileName=plots/positionRecovery.png
    Status should be        SUCCESS
    
Recover published 2008 survey source fluxes - real space filter
    Setup south2008
    Set config              configs/PSTest_south2008.yml
    Run nemo
    Cross match             testsCache/act_source_catalog_AR1_2008.fits     testsCache/PSTest_south2008/PSTest_south2008_optimalCatalog.fits
    Check recovered ratio   fluxJy    fluxJy    tolerance=0.01  SNRKey=SNR  SNRCut=5.0  plotFileName=plots/amplitudeRecovery_M11.png
    Status should be        SUCCESS

Recover published 2008 survey source fluxes - Fourier space filter
    Setup south2008
    Set config              configs/PSTest_south2008_fourier.yml
    Run nemo
    Cross match             testsCache/act_source_catalog_AR1_2008.fits     testsCache/PSTest_south2008_fourier/PSTest_south2008_fourier_optimalCatalog.fits
    Check recovered ratio   fluxJy    fluxJy    tolerance=0.01  SNRKey=SNR  SNRCut=5.0  plotFileName=plots/amplitudeRecovery_M11_fourier.png
    Status should be        SUCCESS
    
*** Keywords ***

Find injected sources
    Setup equD56
    Set config              configs/PSTest_E-D56.yml
    Inject sources using nemo   testsCache/sourceInjectedMap.fits   testsCache/inputSourcesCatalog.fits
    Run nemo

Clean up
    Remove file             testsCache/sourceInjectedMap.fits
    Remove file             testsCache/inputSourcesCatalog.fits
    Remove directory        testsCache/PSTest_E-D56    True
    Remove directory        testsCache/PSTest_south2008    True
    Remove directory        testsCache/PSTest_south2008_fourier    True

*** Settings ***

Library                     OperatingSystem
Library                     lib/NemoTests.py
Suite Setup                 Find injected sources
Suite Teardown              Clean up

*** Variables ***

