#############################################################################################
#############################################################################################
#                                                                                           #
#                      COUNTS2SAR FUNCTION FOR 210PB PROCESSING                             #
#                                EVAN LAHR, 2022                                            #
#                                                                                           #
#                                                                                           #
#   'COUNTS_FNAME' must be a csv with column format:                                        #
#   'Detector Name', 'DepthMidpoint (cm)', '209Po decay (counts)', '210Po decay (counts)',  #
#   'TotalCounts (sec)','Counting_StartDate','Counting_StartTime','Plating_StartDate',      #
#   'Plating_StartTime','Pan (g)','WetSed+Pan (g)','DrySed+Pan (g)','WetChemSedWt (g)',     #
#   'SampleDepthInterval (cm)','VolFrac SiltClay'                                           #
#                                                                                           #
#    'BKG_FNAME' must be a csv with column format:                                          #
#    'Detector Name', 'Counts Po209', 'Counts Po210', 'Counting time (s)'                   #
#                                                                                           #
#    'SUPLVL' (int/float/double) is the background activity (dpm/g) of local sediments      #
#                                                                                           #
#############################################################################################
#############################################################################################


def counts2sar(counts_fname,bkg_fname,supLvl):
    import matplotlib.pyplot as plt
    from datetime import datetime
    import numpy as np
    import pandas as pd

    pd.set_option("precision", 20)

    #Define constants
    supportedLevel     = supLvl
    Collection_DateTime = datetime.strptime('10/01/2021', '%m/%d/%Y')

    lambda_210Pb_min   = 0.00000005914
    lambda_210Po_min   = 0.000003472848
    lambda_209Po_min   = 0.00000001292

    spike_volume       = 0.998
    spike_act_atCalibration = 12.0469862348134
    spike_calibrateDate = datetime.strptime('08/20/2018', '%m/%d/%Y')
    spike_error        = 0.4

    porewater_density  = 1.025
    particle_density   = 2.65
    porewater_saltFrac = 0.025

    pipette_error      = 0.003
    spike_error        = 0.033

    #read in a csv of background activity with the following column names
    #"Detector Name", "counts Po209", "counts Po210", "counting time (sec)"
    bkg = pd.read_csv(bkg_fname)
    bkg

    #read in a csv with the columns listed exactly as in "cts.columns"
    cts = pd.read_csv(counts_fname,float_precision='round_trip',header=None)
    cts.columns =['Detector Name', 'DepthMidpoint (cm)', '209Po decay (counts)', '210Po decay (counts)',
                  'TotalCounts (sec)','Counting_StartDate','Counting_StartTime','Plating_StartDate',
                  'Plating_StartTime','Pan (g)','WetSed+Pan (g)','DrySed+Pan (g)','WetChemSedWt (g)',
                  'SampleDepthInterval (cm)','VolFrac SiltClay']

    #convert csv dates to datetimes
    cts["Plating_DateTime"] = cts["Plating_StartDate"] +' '+ cts["Plating_StartTime"]
    cts['Plating_DateTime']=pd.to_datetime(cts['Plating_DateTime'], infer_datetime_format=True)
    cts["Counting_DateTime"] = cts["Counting_StartDate"] +' '+ cts["Counting_StartTime"]
    cts['Counting_DateTime']=pd.to_datetime(cts['Counting_DateTime'], infer_datetime_format=True)

    #CORRECT FOR BACKGROUND ACTIVITY
    a=[]
    for i in range(len(cts)):
        a.append(bkg.loc[bkg['Detector Name'] == cts['Detector Name'][i]].values[0])
    a=pd.DataFrame(a)
    a=a.drop(labels=[0], axis=1)
    a = a.rename(columns={1: 'BKG counts Pb209', 2: 'BKG counts Pb210', 3: 'BKG counts (sec)'})
    cts = pd.concat([a, cts], axis=1)
    cts['209Po decay minus bkg (counts)']=cts['209Po decay (counts)']-((cts['BKG counts (sec)']/60)*
                                              (cts['BKG counts Pb209']/(cts['BKG counts (sec)']/60)))
    cts['210Po decay minus bkg (counts)']=cts['210Po decay (counts)']-((cts['BKG counts (sec)']/60)*
                                              (cts['BKG counts Pb209']/(cts['BKG counts (sec)']/60)))

    #CALCULATE SAMPLE ACTIVITY
    cts['ElapsedTime_Plate2Count (min)'] = (cts['Counting_DateTime']-cts['Plating_DateTime'])/np.timedelta64(1,'m')
    cts['ElapsedTime_Collect2Plate (min)'] = (cts['Plating_DateTime']-Collection_DateTime)/np.timedelta64(1,'m')
    cts['ElapsedTime_SpikeCal2Count (min)']=(cts['Counting_DateTime']-spike_calibrateDate)/np.timedelta64(1,'m')
    cts['210Po_Decay_Count2Plate']=np.exp(lambda_210Po_min*cts['ElapsedTime_Plate2Count (min)'])
    cts['210Pb_Decay_Collect2Plate']=np.exp((-lambda_210Pb_min)*cts['ElapsedTime_Collect2Plate (min)'])
    cts['CorrAct_210Po_Spike (dpm/g)']=np.exp((-lambda_209Po_min)*cts['ElapsedTime_SpikeCal2Count (min)'])
    cts['WeightFrac_WaterAndSalt'] = ((cts['WetSed+Pan (g)']-cts['Pan (g)'])-(cts['DrySed+Pan (g)']-cts['Pan (g)']))/(cts['WetSed+Pan (g)']-cts['Pan (g)'])
    cts['WeightFrac_Sed'] = 1-cts['WeightFrac_WaterAndSalt']
    cts['WetChemSedWt_SaltCorrected (g)'] = cts['WetChemSedWt (g)']-((cts['WeightFrac_WaterAndSalt']*porewater_saltFrac)*(cts['WetChemSedWt (g)']/cts['WeightFrac_Sed']))
    cts['Raw activity (dpm/g)']=((((cts['210Po decay minus bkg (counts)'])*(cts['210Po_Decay_Count2Plate'])*(spike_volume))*
              ((cts['CorrAct_210Po_Spike (dpm/g)'])*(spike_act_atCalibration)))  /  ((cts['WetChemSedWt (g)'])
              *(cts['209Po decay minus bkg (counts)'])))  /  (cts['210Pb_Decay_Collect2Plate'])
    cts['salt/mud corrected activity (dpm/g)']=((cts['Raw activity (dpm/g)'])*(cts['WetChemSedWt (g)']/cts['WetChemSedWt_SaltCorrected (g)'])-supportedLevel)/cts['VolFrac SiltClay']

    #CALCULATE SEDIMENT BULK DENSITY AND POROSITY
    cts['VolFrac_water+salt']=(cts['WeightFrac_WaterAndSalt']/(1-porewater_saltFrac))*(1/porewater_density)
    cts['VolFrac_sed']       =(cts['WeightFrac_Sed']-(cts['WeightFrac_WaterAndSalt']*(porewater_saltFrac/(1-porewater_saltFrac))))  /  (particle_density)
    cts['Porosity (Uncorr.)']=(cts['WeightFrac_WaterAndSalt']*particle_density)/((cts['WeightFrac_WaterAndSalt']*particle_density)+(1-cts['WeightFrac_WaterAndSalt'])*porewater_density)
    cts['Porosity (Salt corr.)']=(cts['VolFrac_water+salt'])/(cts['VolFrac_water+salt']+cts['VolFrac_sed'])
    cts['BulkDensity_dry'] = (1-cts['Porosity (Salt corr.)'])*(particle_density)
    cts['BulkDensity_wet'] = (1-cts['Porosity (Salt corr.)'])*(particle_density)+(cts['Porosity (Salt corr.)'])*(porewater_density)

    #CALCULATE ERROR
    cts['Po209_error']=((cts['209Po decay (counts)'])**(1/2))/(cts['209Po decay (counts)'])
    cts['Po210_error']=((cts['210Po decay (counts)'])**(1/2))/(cts['210Po decay (counts)'])
    cts['Error_total'] = (((pipette_error)**2)+((spike_error)**2)+((cts['Po210_error'])**2)+((cts['Po209_error'])**2))**(1/2)
    cts['Error'] = cts['Raw activity (dpm/g)']*cts['Error_total']
    cts['Error_SaltCorr'] = (cts['Error']/cts['Raw activity (dpm/g)'])*(cts['salt/mud corrected activity (dpm/g)'])
    cts['Error_MudSaltCorr (Xdir_error)'] = cts['Error_SaltCorr']/cts['VolFrac SiltClay']

    #DROP REDUNDANT COLS
    cts=cts.drop(labels=["Counting_StartDate","Counting_StartTime","Plating_StartDate","Plating_StartTime"], axis=1)
    
    return cts