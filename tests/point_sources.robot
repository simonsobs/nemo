*** Test Cases ***

Recovered source amplitudes are unbiased
    check recovered amplitudes
    status should be    SUCCESS

Recovered source positions are accurate
    check recovered positions
    status should be    SUCCESS
    
*** Keywords ***

Find injected sources
    inject sources
    run nemo

Clean up
    Remove file         sourceInjectedMap.fits
    Remove file         inputSourcesCatalog.fits
    Remove directory    configs/PSTest_E-D56    True

*** Settings ***

Library             OperatingSystem
Library             lib/PointSourceTests.py
Suite Setup         Find injected sources
Suite Teardown      Clean up

*** Variables ***
