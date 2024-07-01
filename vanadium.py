# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

import os
import scipy.optimize

instrument = 'TOPAZ'
output_directory = '/SNS/TOPAZ/shared/Vanadium/2022C_1202_AG/'

pixels_to_mask = ['0-18','237-255']
tubes_to_mask = ['0-18','237-255']

detector_calibration = '2022B/TOPAZ_2022B_AG.DetCal'
tube_calibration = None

bkg_ipts = 31189
bkg_no = 46756
bkg_scale = 1.0

ipts = 31189
run_nos = [46757]

k_min, k_max = 1.8, 18
tof_min, tof_max = None, None

calibration_directory = '/SNS/{}/shared/calibration/'.format(instrument)

file_directory = '/SNS/{}/IPTS-{}/nexus/'
file_name = '{}_{}.nxs.h5'

LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=instrument)

CreateGroupingWorkspace(InputWorkspace=instrument,
                        GroupDetectorsBy='bank',
                        OutputWorkspace='group')

if pixels_to_mask is not None:
    pixels_to_mask = ','.join([str(pixels[0])+'-'+str(pixels[-1]) if type(pixels) is list else str(pixels) for pixels in pixels_to_mask])
    MaskBTP(Workspace=instrument, Pixel=pixels_to_mask)

if tubes_to_mask is not None:
    tubes_to_mask = ','.join([str(tubes[0])+'-'+str(tubes[-1]) if type(tubes) is list else str(tubes) for tubes in tubes_to_mask])
    MaskBTP(Workspace=instrument, Tube=tubes_to_mask)

ExtractMask(InputWorkspace=instrument,
            OutputWorkspace='mask')

if tube_calibration is not None:

    LoadNexus(Filename=os.path.join(calibration_directory, tube_calibration),
              OutputWorkspace='tube_table')
    ApplyCalibration(Workspace=instrument, CalibrationTable='tube_table')

if detector_calibration is not None:

    ext = os.path.splitext(detector_calibration)[1]
    if ext == '.xml':
        LoadParameterFile(Workspace=instrument,
                          Filename=os.path.join(calibration_directory, detector_calibration))
    else:
        LoadIsawDetCal(InputWorkspace=instrument,
                       Filename=os.path.join(calibration_directory, detector_calibration))

files_to_load = '+'.join([os.path.join(file_directory.format(instrument,ipts), file_name.format(instrument,run_no)) for run_no in run_nos])
bkg_file = os.path.join(file_directory.format(instrument,bkg_ipts), file_name.format(instrument,bkg_no)) 

rebin_param = '{},{},{}'.format(k_min,k_max,k_max)

Load(Filename=files_to_load,
     OutputWorkspace='van',
     LoadMonitors=True,
     FilterByTofMin=tof_min,
     FilterByTofMax=tof_max)

SetSample(InputWorkspace='van',
          Geometry={'Shape': 'Sphere', 'Radius': 0.2,'Center': [0.,0.,0.]},
          Material={'ChemicalFormula': 'V', 'UnitCellVolume': 27.642, 'ZParameter': 2.})

# vanadium = CrystalStructure('3.0278 3.0278 3.0278', 'I m -3 m', 'V 0 0 0 1.0 0.00605')

NormaliseByCurrent(InputWorkspace='van',
                   OutputWorkspace='van')

Load(Filename=bkg_file,
     OutputWorkspace='bkg',
     FilterByTofMin=tof_min,
     FilterByTofMax=tof_max)

NormaliseByCurrent(InputWorkspace='bkg',
                   OutputWorkspace='bkg')

mtd['bkg'] *= bkg_scale

Minus(LHSWorkspace='van', RHSWorkspace='bkg', OutputWorkspace='van')

MaskDetectors(Workspace='van', MaskedWorkspace='mask')
ConvertUnits(InputWorkspace='van', OutputWorkspace='van', Target='Momentum')
CropWorkspace(InputWorkspace='van', OutputWorkspace='van', XMin=k_min, XMax=k_max)
Rebin(InputWorkspace='van', OutputWorkspace='van', Params='{},{},{}'.format(k_min,(k_max-k_min)/500, k_max))

ConvertUnits(InputWorkspace='van', OutputWorkspace='van_corr', Target='Wavelength')
#MonteCarloAbsorption(InputWorkspace='van_corr', OutputWorkspace='abs_corr')
AbsorptionCorrection(InputWorkspace='van_corr', OutputWorkspace='abs_corr')
#MultipleScatteringCorrection(InputWorkspace='van_corr', OutputWorkspace='mult_corr')
Divide(LHSWorkspace='van_corr', RHSWorkspace='abs_corr', OutputWorkspace='van_corr')
Divide(LHSWorkspace='van_corr', RHSWorkspace='mult_corr_sampleOnly', OutputWorkspace='van_corr')

ConvertUnits(InputWorkspace='van_corr', OutputWorkspace='van_corr', Target='Momentum')
ConvertUnits(InputWorkspace='abs_corr', OutputWorkspace='abs_corr', Target='Momentum')

Rebin(InputWorkspace='van_corr',
      OutputWorkspace='sa',
      Params=rebin_param,
      PreserveEvents=False)

GroupDetectors(InputWorkspace='van_corr',
               CopyGroupingFromWorkspace='group',
               OutputWorkspace='van_corr')

GroupDetectors(InputWorkspace='van',
               CopyGroupingFromWorkspace='group',
               OutputWorkspace='van')

RemoveMaskedSpectra(InputWorkspace='van', MaskedWorkspace='van', OutputWorkspace='van')
RemoveMaskedSpectra(InputWorkspace='van_corr', MaskedWorkspace='van_corr', OutputWorkspace='van_corr')

def flux(k, phi_0, lambda_0, C_epi, C_fast, alpha, beta):
    therm = phi_0*(4*np.pi**2)/k**2*np.exp(-(4*np.pi**2)/(k**2*lambda_0**2))
    epi = C_epi*k**3/(2*np.pi)**3
    fast = C_fast*(2*np.pi)**alpha/k**alpha*np.exp(-(2*np.pi*beta)/k)
    return therm+epi+fast

def residuals(params, k, y, e):
    return (flux(k, *params)-y)/e

y = mtd['van'].extractY()
x = mtd['van'].extractX()
e = mtd['van'].extractE()

k = (x[:,:-1]+x[:,1:])/2

y_corr = mtd['van_corr'].extractY()
x_corr = mtd['van_corr'].extractX()
e_corr = mtd['van_corr'].extractE()

y_mult, y_abs = mtd['mult_corr_sampleOnly'].extractY(), mtd['abs_corr'].extractY()
y_abs = mtd['abs_corr'].extractY()
y_fact = y_mult*y_abs

k_corr = (x_corr[:,:-1]+x_corr[:,1:])/2

for i in range(k.shape[0]):

    fig, ax = plt.subplots(2, 2, sharex='col', sharey='row')

    x0 = [0, 0, 0, 0, 0, 0]
    args = (k[i,:], y[i,:], e[i,:])

    sol = scipy.optimize.least_squares(residuals,
                                       x0=x0,
                                       args=args,
                                       loss='soft_l1')

    phi_0, lambda_0, C_epi, C_fast, alpha, beta = sol.x

    fit = flux(k[i,:], *sol.x)

    ax[0,0].errorbar(k[i,:], y[i,:]/y[i,:].max(), e[i,:]/y[i,:].max(), fmt='.', label='orig')
    ax[0,0].plot(k[i,:], fit/y[i,:].max(), '-', label='orig')
    ax[0,1].errorbar(2*np.pi/k[i,:], y[i,:]/y[i,:].max(), e[i,:]/y[i,:].max(), fmt='.', label='orig')
    ax[0,1].plot(2*np.pi/k[i,:], fit/y[i,:].max(), '-', label='orig')

    x0 = [0, 0, 0, 0, 0, 0]
    args = (k_corr[i,:], y_corr[i,:], e_corr[i,:])

    sol = scipy.optimize.least_squares(residuals,
                                       x0=x0,
                                       args=args,
                                       loss='soft_l1')

    phi_0, lambda_0, C_epi, C_fast, alpha, beta = sol.x

    fit = flux(k_corr[i,:], *sol.x)

    ax[0,0].errorbar(k_corr[i,:], y_corr[i,:]/y_corr[i,:].max(), e_corr[i,:]/y_corr[i,:].max(), fmt='.', label='corr')
    ax[0,0].plot(k_corr[i,:], fit/y_corr[i,:].max(), '-', label='corr')
    ax[0,1].errorbar(2*np.pi/k_corr[i,:], y_corr[i,:]/y_corr[i,:].max(), e_corr[i,:]/y_corr[i,:].max(), fmt='.', label='corr')
    ax[0,1].plot(2*np.pi/k_corr[i,:], fit/y_corr[i,:].max(), '-', label='corr')

    ax[1,0].plot(k_corr[i,:], 1/y_abs[i,:]*y_abs[i,:].max(), '-', label='abs')
    #ax[1,0].plot(k_corr[i,:], 1/y_mult[i,:]*y_mult[i,:].max(), '-', label='mult')
    #ax[1,0].plot(k_corr[i,:], 1/y_fact[i,:]*y_fact[i,:].max(), '-', label='corr')

    ax[1,1].plot(2*np.pi/k_corr[i,:], 1/y_abs[i,:]*y_abs[i,:].max(), '-', label='abs')
    #ax[1,1].plot(2*np.pi/k_corr[i,:], 1/y_mult[i,:]*y_mult[i,:].max(), '-', label='mult')
    #ax[1,1].plot(2*np.pi/k_corr[i,:], 1/y_fact[i,:]*y_fact[i,:].max(), '-', label='corr')

    ax[1,0].legend()
    ax[1,1].legend()

    ax[1,0].set_xlabel('$k$ [$\AA^{-1}$]')
    ax[1,1].set_xlabel('$\lambda$ [$\AA$]')
    ax[0,0].axvline(2*np.pi/0.75, color='k', linestyle='--')
    ax[1,0].axvline(2*np.pi/0.75, color='k', linestyle='--')
    ax[0,1].axvline(0.75, color='k', linestyle='--')
    ax[1,1].axvline(0.75, color='k', linestyle='--')
    ax[0,0].legend()
    ax[0,1].legend()
    fig.show()

plt.close('all')

Rebin(InputWorkspace='van_corr',
      OutputWorkspace='van_flux',
      Params=rebin_param)

flux = mtd['van_flux']
for i in range(flux.getNumberHistograms()):
    el=flux.getSpectrum(i)
    if flux.readY(i)[0] > 0:
        el.divide(flux.readY(i)[0],flux.readE(i)[0])

SortEvents(InputWorkspace='van_flux', SortBy='X Value')

IntegrateFlux(InputWorkspace='van_flux', OutputWorkspace='flux', NPoints=1000)

SaveNexus(InputWorkspace='sa', Filename=os.path.join(output_directory, 'solid_angle.nxs'))
SaveNexus(InputWorkspace='flux', Filename=os.path.join(output_directory, 'flux.nxs'))