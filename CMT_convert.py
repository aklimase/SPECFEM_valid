from pyrocko import moment_tensor as mtm

#magnitude = 3.818# Magnitude of the earthquake
magnitude = 4.0 # Magnitude of the earthquake

#5.984e+14 N-m

m0 = mtm.magnitude_to_moment(magnitude)  # convert the mag to moment
print m0

strike = 297
dip = 38
rake = 90
mt = mtm.MomentTensor(strike=strike, dip=dip, rake=rake, scalar_moment=m0)

m6 = [mt.mnn, mt.mee, mt.mdd, mt.mne, mt.mnd, mt.med]  # The six MT components
print(m6)

#Nm to dyne cm = 10^7
m6 = [mt.mnn*10.**7, mt.mee*10.**7, mt.mdd*10.**7, mt.mne*10.**7, mt.mnd*10.**7, mt.med*10.**7]  # The six MT components
print m6