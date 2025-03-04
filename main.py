from spectra_processing import *

if __name__ == '__main__':
 
    savefold = r'D:\PhD\oceanview\Spectra_stars\\'
    # savefold = r'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/'
    path = r'D:\PhD\oceanview\16012025\\'
    # path = r'/media/astronomer/Transcend/PhD/oceanview/16012025/'
    path1 = r'D:\PhD\oceanview\Asiago23052024\\'
    # path1 = r'/media/astronomer/Transcend/PhD/oceanview/Asiago23052024/'
    path2 = r'D:\PhD\oceanview\17012025\\'
    # path2 = r'/media/astronomer/Transcend/PhD/oceanview/17012025/'
    path_at = r'D:\PhD\oceanview\\'
    # path_at = r'/media/astronomer/Transcend/PhD/oceanview/'
    file_at = 'FG_UGA_LGA_UCA_LCA_UEA_LEA_Attenuation_Data_2.xlsx'

    spectral_lines = {
    "He I": [402.6, 447.3, 492., 587.6,  667.8],#447.1
    "Hα": [656.],
    "Hβ": [486.1],
    "Hγ": [434.0], #434.0
    "Hδ": [410.2],
    "He II": [420.0, 454.1],
    "Fe I": [430.9, 438.4, 495.8, 466.8, 527.0],
    # "Ca I": [420.8],
    # "Fe II": [516.90],
    "Mg I": [518.0],
    # "Na I D1": [589.00],
    # "Na I D2": [589.60],
    # "Ca II H": [396.85],
    # "Ca II K": [393.37],
    "Ca II IR 1": [849.80],
    "Ca II IR 2": [854.20],
    "Ca II IR 3": [866.20],
    "[O I] 1": [630.0],
    "[O I] 2": [636.4],
    # "O II": [501.5],
    # "C II": [426.7],
    "Ti I": [626.5],
    "Si II": [412.8, 634.7, 637.1],
    "Mg II": [448.1],
    # "O I": [627.7, 686.7, 759.4, 777.5, 822.7, 898.8], #most of these are telluric
    "TiO": [476.1, 495.4, 516.7],
    "CH4": [509., 543., 619., 576., 597., 668., 727.67, 863.],
    "NH3": [648.]
    }

    telluric_lines = {"TL":
    [687.,
     718.5,
     719.4,
     725.0,
     759.3,  # O2 (760.0 nm region)
     760.5,
     761.0,  # O2 (760.5 nm region)
     762.0,  # O2 (762.0 nm region)
     764.0,  # H2O (763.0-764.0 nm region)
     790.0,  # H2O (790.0region)
     701.5,  # H2O
     820.5,  # H2O (820.0-821.0 nm region)
     822.2,  # H2O (822.0 nm region)
     892.,  #
     935.0,  # H2O (934.0-935.0 nm region)
     940.0,  # H2O (940.0 nm region)
     942.0,  # H2O (942.0 nm region)
     940.5,  # H2O (940.5 nm region)
     943.0,  # H2O (943.0 nm region)
     946.0,  # H2O (946.0 nm region)
     953.0,  # H2O (953.0 nm region)
     960.0,  # H2O (960.0 nm region)
]}

    threshold = 0.02

    #=======================================================================================#
    #                                  Zeta Cas Spectrophometric Standard                   #
    #=======================================================================================#
    #Reference Data Read
    # ref_spectra = LoadData(path_at+'zetaCas_standard.txt',skiprows=0)
    # ref_spectra['wavelength']/=10
    #
    # #######
    #
    # spectra = LoadAllSpectra(path, 'zetaCas_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'zetaCas_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # #flat = LoadAllSpectra('/Users/taniamachado/Documents/PhD/Oceanview/Asiago_Jan16-18/Flats_16012025/', 'flat_16s_HDX017711')
    #
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # #flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # #flat_av = AverageSpectra(flat_subtr)
    #
    #
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = Build_Instrument_Response(sum_spectra, ref_spectra, minwavelength = 370, maxwavelength = 1100, plot=True)
    # print(IRF.shape)
    # IRF.to_csv(path_at + 'IRF.csv', index = False)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # epsilon = 1e-10
#
#     wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 370, 1100)
#     intensity_truncated = intensity_truncated/(IRF['intensity'].values + epsilon)
#
#
#     spec_norm = FitContinuum(10, wavelength_truncated, intensity_truncated, image = True)
    # spec_norm = BlackBodyFit(10000, wavelength_truncated, intensity_truncated, image = True)
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'ζ Cas, ΗR153 [m(V) = 3.66]', save = True, path = savefold, name = "zetaCas_spectrum_lines.png")

    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], AddSpectra(spectra), 'zetaCas', 'zetaCas raw data')
    # PlotData(axs[0], background_av, 'zetaCas', 'background')
    # PlotData(axs[0], dark_av, 'zetaCas', 'dark')
    # PlotData(axs[1], sum_spectra, 'zetaCas, x [m(V) = x]', 'zetaCas processed data')
    # plt.tight_layout()
    # # plt.savefig(savefold+'zetaCas.png', dpi = 600)

    # plt.figure()
    # plt.plot(wavelength_truncated, spec_norm)
    # plt.show()

    #=======================================================================================#
    #                                  HD93521                 #
    #=======================================================================================#
    # Reference Data Read
    # ref_spectra = LoadData(path_at+'fhd93521.dat.txt',skiprows=0)
    # ref_spectra['wavelength']/=10
    #
    # #######
    #
    # spectra = LoadAllSpectra(path2, 'HD93521_400umfiber_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path2, 'HD93521_400umfiber_background_HDX017711')
    #
    # # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    #
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # #flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # #flat_av = AverageSpectra(flat_subtr)
    #
    #
    # background_subtr = ProcessDataEach(background, None, None, None, subtract_background = False, subtract_dark = False, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, background_subtr_av, None, None, subtract_background = True, subtract_dark = False, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = Build_Instrument_Response(sum_spectra, ref_spectra, minwavelength = 380, maxwavelength = 1100, plot=True)
    # print(IRF.shape)
    # IRF.to_csv(path_at + 'IRF_HD93521.csv', index = False)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values + epsilon)
    #
    #
    # spec_norm = FitContinuum(12, wavelength_truncated, intensity_truncated, image = True)
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'HD93521 [m(V) = 7.03]', save = False, path = savefold, name = "HD93521_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], AddSpectra(spectra), 'HD93521', 'HD93521 raw data')
    # PlotData(axs[0], background_av, 'HD93521', 'background')
    # PlotData(axs[1], sum_spectra, 'HD93521 , x [m(V) = x]', 'HD93521 summed data')
    # plt.tight_layout()
    #
    # # plt.figure()
    # # plt.plot(wavelength_truncated, spec_norm)
    # plt.show()

    #=======================================================================================#
    #                                  Castor                                               #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'Castor_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Castor_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_5s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_5s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(16, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Castor, α Geminorum [m(V) = 1.58]', save = True, path = savefold, name = "Castor_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Castor', 'Castor raw data')
    # PlotData(axs[0], background_av, 'Castor', 'background')
    # PlotData(axs[0], dark_av, 'Castor', 'dark')
    # PlotData(axs[0], flat_av, 'Castor', 'flat')
    # # PlotData(axs[1], IRF, 'Capella', 'IRF')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Castor processed data')
    # plt.tight_layout()
    #
    # plt.show()

    # #=======================================================================================#
    # #                                  Gamma Uma                                              #
    # #=======================================================================================#
    #
    # spectra = LoadAllSpectra(path, 'gammaUma_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'gammaUma_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_1s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_1s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(17, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Phecda, γ Ursae Majoris [m(V) = 2.44]', save = False, path = savefold, name = "gammaUma_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'gammaUma', 'gammaUma raw data')
    # PlotData(axs[0], background_av, 'gammaUma', 'background')
    # PlotData(axs[0], dark_av, 'gammaUma', 'dark')
    # PlotData(axs[0], flat_av, 'gammaUmar', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'gammaUma processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Jupiter                                              #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'Jupiter_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Jupiter_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_2s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_2s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 385, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 385)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(25, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Jupiter [m(V) = -2.7]', save = False, path = savefold, name = "Jupiter_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Jupiter', 'Jupiter raw data')
    # PlotData(axs[0], background_av, 'Jupiter', 'background')
    # PlotData(axs[0], dark_av, 'Jupiter', 'dark')
    # PlotData(axs[0], flat_av, 'Jupiter', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Jupiter processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Ganymede                                             #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'Jupiter_Ganymede_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Jupiter_Ganymede_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_16s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(25, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Ganymede [m(V) = 4.6]', save = True, path = savefold, name = "Ganymede_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Ganymede', 'Ganymede raw data')
    # PlotData(axs[0], background_av, 'Ganymede', 'background')
    # PlotData(axs[0], dark_av, 'Ganymede', 'dark')
    # PlotData(axs[0], flat_av, 'Ganymede', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Ganymede processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                 Mars                                            #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'Mars_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Mars_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_3s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_3s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(25, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Mars [m(V) = -2.3]', save = True, path = savefold, name = "Mars_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Mars', 'Mars raw data')
    # PlotData(axs[0], background_av, 'Mars', 'background')
    # PlotData(axs[0], dark_av, 'Mars', 'dark')
    # PlotData(axs[0], flat_av, 'Mars', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Mars processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Uranus                                             #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'Uranus_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Uranus_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_16s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(38, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Uranus [m(V) = 5.7]', save = True, path = savefold, name = "Uranus_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Uranus', 'Uranus raw data')
    # PlotData(axs[0], background_av, 'Uranus', 'background')
    # PlotData(axs[0], dark_av, 'Uranus', 'dark')
    # PlotData(axs[0], flat_av, 'Uranus', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Uranus processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Tau Aurigae                                             #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, 'TauAurigae_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'TauAurigae_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_16s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    # # remove artefact
    # artifact_indices = np.where((wavelength_truncated > 1030) & (wavelength_truncated < 1032))[0]
    # wavelength_truncated = np.delete(wavelength_truncated, artifact_indices)
    # intensity_truncated = np.delete(intensity_truncated, artifact_indices)
    #
    # spec_norm = FitContinuum(18, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'τ Aurigae [m(V) = 4.5]', save = True, path = savefold, name = "TauAurigae_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'TauAurigae', 'TauAurigae raw data')
    # PlotData(axs[0], background_av, 'TauAurigae', 'background')
    # PlotData(axs[0], dark_av, 'TauAurigae', 'dark')
    # PlotData(axs[0], flat_av, 'TauAurigae', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'TauAurigae processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Eta Hydrae                                             #
    #=======================================================================================#

    #Reference Data Read
    # ref_spectra = LoadData(os.path.join(path_at, 'fhr3454.dat.txt'), skiprows=0)
    # # print(f"Final path: <{os.path.join(path_at, 'fhr3454.dat')}")
    # # with open(os.path.join(path_at, 'fhr3454.dat'), 'r') as f:
    # #     print("File opened successfully!")
    # ref_spectra['wavelength']/=10
    # spectra = LoadAllSpectra(path2, 'etaHya_400umfiber_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path2, 'etaHya_400umfiber_background_HDX017711')
    #
    # # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_16s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, None, None, subtract_background = False, subtract_dark = False, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, None, None, subtract_background = False, subtract_dark = False, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, None, None, subtract_background = True, subtract_dark = False, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1090)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # # IRF = IRF[(IRF['wavelength'] >= 380)]
    # IRF = Build_Instrument_Response(sum_spectra, ref_spectra, minwavelength = 380, maxwavelength = 1090, plot=True)
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/np.max(intensity_truncated)/(IRF['intensity'].values)
    #
    # # remove artefact
    # artifact_indices = np.where((wavelength_truncated > 1030) & (wavelength_truncated < 1032))[0]
    # wavelength_truncated = np.delete(wavelength_truncated, artifact_indices)
    # intensity_truncated = np.delete(intensity_truncated, artifact_indices)
    #
    # spec_norm = FitContinuum(12, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'η Hydrae [m(V) = 4.3]', save = False, path = savefold, name = "etaHya_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'eta Hya', 'eta Hya raw data')
    # PlotData(axs[0], background_av, 'eta Hya', 'background')
    # # PlotData(axs[0], dark_av, 'eta Hydrae', 'dark')
    # PlotData(axs[0], flat_av, 'eta Hydrae', 'flat')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'eta Hydrae processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  63 Aurigae                                              #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path, '63Aurigae_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, '63Aurigae_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_16s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_16s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(18, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = '63Aurigae [m(V) = 4.89]', save = True, path = savefold, name = "63Aurigae_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], '63Aurigae', '63Aurigae raw data')
    # PlotData(axs[0], background_av, '63Aurigae', 'background')
    # PlotData(axs[0], dark_av, '63Aurigae', 'dark')
    # PlotData(axs[0], flat_av, '63Aurigae', 'flat')
    # # PlotData(axs[1], IRF, 'Capella', 'IRF')
    #
    #
    # PlotData(axs[1], sum_spectra, None, '63Aurigae processed data')
    # plt.tight_layout()
    #
    # plt.show()

    #=======================================================================================#
    #                                  Capella                                               #
    #=======================================================================================#
    
    # path_aux = path+'Capella/'
    # spectra = LoadAllSpectra(path, 'Capella_data_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path, 'Capella_background_HDX017711')
    #
    # dark = LoadAllSpectra(path, 'dark_1s_HDX017711')
    # ###print(spectra[0])
    # ###print(background[0])
    # flat = LoadAllSpectra(path, 'flat_1s_HDX017711')
    #
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # flat_subtr = ProcessDataEach(flat, background_av, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    #
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ####print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 1100)
    # print(min(wavelength_truncated))
    # print(intensity_truncated.shape)
    # IRF = IRF[(IRF['wavelength'] >= 380)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    #
    # spec_norm = FitContinuum(18, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Capella, α Aurigae [m(V) = 0.08]', save = False, path = savefold, name = "Capella_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'Capella', 'Capella raw data')
    # PlotData(axs[0], background_av, 'Capella', 'background')
    # PlotData(axs[0], dark_av, 'Capella', 'dark')
    # PlotData(axs[0], flat_av, 'Capella', 'flat')
    # # PlotData(axs[1], IRF, 'Capella', 'IRF')
    #
    #
    # PlotData(axs[1], sum_spectra, None, 'Capella processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Capella.png', dpi = 600)
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Capella_spectrum.png', dpi = 600)
    # plt.show()

    #=======================================================================================#
    #                                  HR6212                                               #
    #=======================================================================================#

    # spectra = LoadAllSpectra(path1, 'HR6212_spectrum_HDX017711')
    # ##print(spectra)
    #
    # background = LoadAllSpectra(path1, 'HR6212_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_HDX017711')
    # #print(spectra[0])
    # #print(background[0])
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # #print(background_av)
    # dark_av = AverageSpectra(dark)
    # #print(dark_av)
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ##print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 940)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 940)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # spec_norm = FitContinuum(22, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'HR6212, ζ Herculis [m(V) = 2.81]', save = True, path = savefold, name = "HR6212_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'HR6212', 'HR6212 raw data')
    # PlotData(axs[0], background_av, 'HR6212', 'background')
    # PlotData(axs[0], dark_av, 'HR6212, ζ Herculis [m(V) = 2.81]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'HR6212 processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/HR6212.png', dpi = 600)
    # plt.show()

    #=======================================================================================#
    #                                  Arcturus                                               #
    #=======================================================================================#
    #

    # spectra = LoadAllSpectra(path1, 'Arcturus_spectrum_HDX017711')
    # # print(spectra)
    #
    # background = LoadAllSpectra(path1, 'Arcturus_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_3s_HDX017711')
    #
    # flat = LoadAllSpectra(path1, 'flat_HDX017711')
    # ##print(spectra[0])
    # ###print(background[0])
    # fiber_at = LoadAttenuation(path_at, file_at)
    #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ##print(background_av)
    # dark_av = AverageSpectra(dark)
    #
    # flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    # ##print(dark_av)
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(data_atcor, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ##print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 940)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 940)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # # spec_norm = BlackBodyFit(4200, wavelength_truncated, intensity_truncated, True)
    # spec_norm = FitContinuum(16, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Arcturus, α Boötis [m(V) = -0.05]', save = False, path = savefold, name = "Arcturus_spectrum_lines.png")
    # #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # # PlotData(axs[0], spectra[0], 'Arcturus', 'Arcturus raw data')
    # PlotData(axs[0], background_av, 'Arcturus', 'background')
    # PlotData(axs[0], flat_av, 'Arcturus', 'flat')
    # PlotData(axs[0], dark_av, 'Arcturus, α Boötis [m(V) = -0.05]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Arcturus processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Arcturus.png', dpi = 600)
    #
    # plt.show()
    #
    # #=======================================================================================#
    # #                                  Deneb                                              #
    # #=======================================================================================#
    #
    # spectra = LoadAllSpectra(path1, 'Deneb_spectrum_HDX017711')
    # #print(spectra)
    #
    # background = LoadAllSpectra(path1, 'Deneb_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_HDX017711')
    #
    # flat = LoadAllSpectra(path1, 'flat_HDX017711')
    #
    # ###print(background[0])
    # fiber_at = LoadAttenuation(path_at, file_at)
    # #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # #print(background_av)
    # dark_av = AverageSpectra(dark)
    # #print(dark_av)
    #
    # flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # flat_av = AverageSpectra(flat_subtr)
    # ##print(spectra[0])
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(data_atcor, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ###print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 930)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 930)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # spec_norm = FitContinuum(34, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Deneb, α Cygni [m(V) = 1.25]', save = False, path = savefold, name = "Deneb_spectrum_lines.png")
    # #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], None, 'Deneb raw data')
    # PlotData(axs[0], background_av, None, 'background')
    # PlotData(axs[0], dark_av, 'Deneb, α Cygni [m(V) = 1.25]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Deneb processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Deneb.png', dpi = 600)
    #
    # plt.show()
    #
    # #=======================================================================================#
    # #                                  alpha Oph                                              #
    # #=======================================================================================#
    # # #
    # spectra = LoadAllSpectra(path1, 'alphaOph_spectrum_HDX017711')
    # #print(spectra)
    #
    # background = LoadAllSpectra(path1, 'alphaOph_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_HDX017711')
    # ##print(spectra[0])
    # #print(background[0])
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # #print(background_av)
    # dark_av = AverageSpectra(dark)
    # #print(dark_av)
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # #print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 930)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 930)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # spec_norm = FitContinuum(13, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Rasalhague, α Ophiuchi [m(V) = 2.07]', save = False, path = savefold, name = "Rasalhague_spectrum_lines.png")
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], None, 'Rasalhague raw data')
    # PlotData(axs[0], background_av, None, 'background')
    # PlotData(axs[0], dark_av, 'Rasalhague, α Ophiuchi [m(V) = 2.07]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Rasalhague processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/alphaOph.png', dpi = 600)
    # plt.show()
    #
    # #=======================================================================================#
    # #                                  Alkaid                                             #
    # #=======================================================================================#
    #
    # spectra = LoadAllSpectra(path1, 'Alkaid_spectrum_HDX017711')
    # ###print(spectra)
    #
    # background = LoadAllSpectra(path1, 'Alkaid_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_HDX017711')
    # ##print(spectra[0])
    # ###print(background[0])
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # background_av = AverageSpectra(background)
    # ###print(background_av)
    # dark_av = AverageSpectra(dark)
    # ###print(dark_av)
    # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, None, subtract_background = True, subtract_dark = True, divide_flat = False)
    # ###print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 960)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 960)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Alkaid, η Ursae Majoris [m(V) = 1.86]', save = False, path = savefold, name = "Alkaid_spectrum_lines.png")
    #
    # #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], None, 'Alkaid raw data')
    # PlotData(axs[0], background_av, None, 'background')
    # PlotData(axs[0], dark_av, 'Alkaid, η Ursae Majoris [m(V) = 1.86]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Alkaid processed data')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Alkaid.png', dpi = 600)
    # plt.show()
    # #
    # #=======================================================================================#
    # #                                  delta Sge  (no background data)                      #
    # #=======================================================================================#
    #
    # spectra = LoadAllSpectra(path1, 'deltaSag_spectrum_HDX017711')
    # #print(spectra)
    #
    # # background = LoadAllSpectra(path1, 'deltaSag_background_HDX017711')
    #
    # dark = LoadAllSpectra(path1, 'dark_HDX017711')
    # # print(spectra[0])
    # #print(background[0])
    # fiber_at = LoadAttenuation(path_at, file_at)
    # # #
    # data_atcor = CorrectAttenuation(spectra, fiber_at)
    #
    # sum_spectra_raw = AddSpectra(spectra)
    # # background_av = AverageSpectra(background)
    # #print(background_av)
    # dark_av = AverageSpectra(dark)
    # #print(dark_av)
    # # background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # # background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(data_atcor, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    # #print(f'subtr: \n {subtr}')
    #
    # sum_spectra = AddSpectra(subtr)
    # IRF = pd.read_csv(path_at+"IRF.csv")
    # epsilon = 1e-10
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 380, 940)
    # IRF = IRF[(IRF['wavelength'] >= 380) & (IRF['wavelength'] <= 940)]
    # print(IRF['wavelength'].min())
    # interpolate = interp1d(IRF['wavelength'], IRF['intensity'], kind='linear', fill_value="extrapolate")
    # IRF['intensity'] = interpolate(wavelength_truncated)
    # # irf_norm = IRF/np.median(IRF)
    # intensity_truncated = intensity_truncated/(IRF['intensity'].values)
    #
    # spec_norm = FitContinuum(14, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'δ Sagittae [m(V) = 4.21 (3.82)]', save = False, path = savefold, name = "DeltaSge_spectrum_lines.png")
    # #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], None, 'delta Sag raw data')
    # # PlotData(axs[0], background_av, None, 'background')
    # PlotData(axs[0], dark_av, 'Kaus Media, δ Sagittarii [m(V) = 2.7]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'delta Sag processed data (only dark subtracted)')
    # plt.tight_layout()
    # # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/deltaSag.png', dpi = 600)
    # plt.show()

