
import pandas as pd  
import numpy as np


'''
los leemos de uno en uno para controlar el formato 

'''
################################################
# filtro J 
manufact = 'J 1250-160nm BARR ED626-1 data.xls'
cold = "filterJ_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### MKO J band 1250/160nm (1170-1330nm)\n"]
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Cold Test Witness',skiprows=11)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['wl(nm)'],df1['T%.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'H 1635-290nm BARR ED561-1 data.xls'
cold = "filterH_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["## MKO H band 1635/290nm (1490-1780nm) - ED561 \n"]
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Cold Test Witness',skiprows=11)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['wl(nm)'],df1['T%.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro Ks 
manufact = 'Ks 2150-320nm BARR ED562-2 data.xls'
cold = "filterKs_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### MKO K band 2150/320nm (1990-2310nm)- ED562 \n"]
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Cold Test Witness',skiprows=11)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['wl(nm)'],df1['T%.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro Ks 
manufact = 'Y 1020-100nm BARR ED655-1 data.xls'
cold = "filterY_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### MKO Y band 1020/100nm (970-1070nm)  \n"]
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Cold Test Witness',skiprows=11)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['wl(nm)'],df1['T%.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro 
manufact = 'UNAM 1257nm-14nm WO#110452 Data Pack.xls'
cold = "filterNB1p26_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### NB1.26 1257/14nm 1.5inch\n"]
comments.append("### FWHM: 14nm ± 0.2% (2.514nm)\n")
comments.append("### CWL: 1257nm ± 0.2% (2.514nm)\n")
comments.append("### Peak %T: >80%")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro 
manufact = 'UNAM 1282nm-20nm WO#110453 Data Pack.xls'
cold = "filterNB1p28_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["###\n"]
comments.append("### CWL: 1282nm ± 0.2% (3.287nm)\n")
comments.append("### FWHM: 20nm ± 0.2% (3.287nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 1643.5nm-25nm WO#109227 Data Pack.xls'
cold = "filterNB1p64_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 1643.5nm ± 0.2% (3.287nm)\n")
comments.append("### FWHM: 25nm ± 0.2% (3.287nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 1710nm-26nm WO#109228 Data Pack.xls'
cold = "filterNB1p71_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 1710nm ± 0.2% (3.42nm)\n")
comments.append("### FWHM: 26nm ± 0.2% (3.42nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2093nm-30nm WO#109229 Data Pack.xls'
cold = "filterNB2p09_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2093nm ± 0.2% (4.186nm)\n")
comments.append("### FWHM: 30nm ± 0.2% (4.186nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2122nm-32nm WO#109230 Data Pack.xls'
cold = "filterNB2p12_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2122nm ± 0.2% (4.244nm)\n")
comments.append("### FWHM: 32nm ± 0.2% (4.244nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 9'],df1['Witness.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2166nm-32nm WO#109232 Data Pack.xls'
cold = "filterNB2p17_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2166nm ± 0.2% (4.332nm) \n")
comments.append("### FWHM: 32nm ± 0.2% (4.332nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2266nm-27nm WO#109233 Data Pack.xls'
cold = "filterNB2p27_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2266nm ± 0.2% (4.533nm)\n")
comments.append("### FWHM: 27nm ± 0.2% (4.533nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2420nm-60nm WO#110766 Data Pack.xls'
cold = "filterNB2p42_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2420nm ± 0.2% (4.84nm)\n")
comments.append("### FWHM: 60nm ± 0.2% (4.84nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM 2480nm-60nm WO#110768 Data Pack.xls'
cold = "filterNB2p48_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### \n"]
comments.append("### CWL: 2480nm ± 0.2% (4.9nm)\n")
comments.append("### FWHM: 60nm ± 0.2% (4.9nm)\n")
comments.append("### Peak %T: >80%\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 10'],df1['Data.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM ND2 1000nm-2500nm WO#109236 Data Pack.xls'
cold = "filterND2_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### ND2 1000nm-2500nm 1.5inch\n"]
comments.append("### \n")
comments.append("### OD2 ± 0.2 Avg from 1000-2500nm\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=4)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 7'],df1['Piece 1.2'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM WB 1450-2500nm WO#109329 Data Pack.xls'
cold = "filterBF-HK_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### BF-HK 1450-2500nm 1.5inch\n"]
comments.append("### Cut-on: 1450nm ± 2% \n")
comments.append("### Out of Band Blocking: OD4 from 300-1250nm \n")
comments.append("### > 85% avg over the central 90% of the 80% BW\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=5)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 4'],df1['Witness.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM WB 1990-2310nm WO#110356 Data Pack.xls'
cold = "filterBF-K_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### Ks 1990-2310nm 1.5inch\n"]
comments.append("### Cut-on: 1990nm ± 0.5%\n")
comments.append("### Cut-off: 2310nm ± 0.5%\n")
comments.append("### > 80% avg over the central 90% of the 80% BW\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=5)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 4'],df1['Witness.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

################################################
# filtro H 
manufact = 'UNAM WB 900-1350nm WO#109242 Data Pack.xls'
cold = "filterBF-ZJ_cold.dat"
#filters=np.dtype({'manufact': ['col1', 'col2'], 'formats': ['i4','f4']})
#dtype([('col1', '<i4'), ('col2', '<f4')])
comments=["### BF-ZJ 900-1350nm 1.5inch\n"]
comments.append("### Cut-on: 900nm ± 2%\n")
comments.append("### Cut-off: 1350nm ± 2%\n")
comments.append("### > 85% avg over the central 90% of the 80% BW\n")
comments.append("## Wavelength[micron] Transmission[%]\n")
xl = pd.ExcelFile(manufact)
## leemos la hoja Cold Test Witness
df1 = xl.parse('Correlation Data',skiprows=5)

aguardar = ["{:.4f} {:.5f} \n".format(wave/1000.,ctrans) for wave,ctrans in zip(df1['Unnamed: 4'],df1['Witness.1'])] 
cold_file = open(cold, "w")
cold_file.writelines(comments)
cold_file.writelines(aguardar)
cold_file.close()

