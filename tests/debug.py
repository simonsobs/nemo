"""

For debugging test module(s) under lib/

"""

from lib import NemoTests

n=NemoTests.NemoTests()

# 2008 stuff
n.setup_south2008()

# Point sources
#n.cross_match("inputSourcesCatalog.fits", "configs/PSTest_E-D56/PSTest_E-D56_optimalCatalog.fits")
#n.check_recovered_ratio("I", "deltaT_c", tolerance=0.01, SNRKey="SNR", SNRCut=5.0, plotFileName="plots/amplitudeRecovery.png")
#n.check_recovered_positions(plotFileName="plots/positionRecovery.png")

#c.run_nemo_mass(catalogFileName = "configs/equD56_quick/mocks/mockCatalog_1.fits")
