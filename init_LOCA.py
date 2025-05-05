def LOCA_getPrecip():
  #Initial Code by Colin Evans 
  #Returns LOCA precipitation to dir
  # =============================================================================
  # LOCA
  # =============================================================================
  
  # List of climate models used in the LOCA dataset. These models provide historical
  # and future climate projections downscaled to a high resolution.
  models = [
      "ACCESS1-0","ACCESS1-3","CCSM4","CESM1-BGC","CESM1-CAM5","CMCC-CM","CMCC-CMS",
      "CNRM-CM5","CSIRO-Mk3-6-0","CanESM2","EC-EARTH","FGOALS-g2","GFDL-CM3",
      "GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R","HadGEM2-AO","HadGEM2-CC","HadGEM2-ES",
      "IPSL-CM5A-LR","IPSL-CM5A-MR","MIROC-ESM","MIROC-ESM-CHEM","MIROC5",
      "MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","bcc-csm1-1","bcc-csm1-1-m",
      "inmcm4"
  ]
  
  # List of selected locations (longitude, latitude).
  # These are the coordinates of the stations where IDF curves will be computed.
  locs = [
      (-80.21448, 40.4846),  # Approx. Lat/Lon for Pittsburgh Int. Airport Station
  ]
  
  # Create an empty array to store precipitation data from all models and locations.
  all_precip = {}
  for scen in ['historical', 'rcp45', 'rcp85']:
      if scen == 'historical':
        all_precip[scen]= np.empty([len(models), len(locs), 20454])
      else:
        all_precip[scen]= np.empty([len(models), len(locs), 34333])
  
  # Loop through each climate model
  for i in range(len(models)):
    for scen in ['historical', 'rcp45', 'rcp85']:
      # Array to store precipitation data for each model
      if scen == 'historical':
        # Shape: (8, 8, 20454) - storing precipitation over a small 8x8 grid around each location
        precip = np.zeros([8, 8, 20454])
      else:
        precip = np.zeros([8, 8, 34333])
  
      # Loop through each station location
      for j in range(len(locs)):
          lon1 = locs[j][0] - (3 * 0.0625)  # Define western boundary of bounding box
          lat1 = locs[j][1] - (3 * 0.0625)  # Define southern boundary of bounding box
          lon2 = locs[j][0] + (3 * 0.0625)  # Define eastern boundary of bounding box
          lat2 = locs[j][1] + (3 * 0.0625)  # Define northern boundary of bounding box
          bbox = [lon1, lat1, lon2, lat2]  # Define bounding box (small region around station)
  
          # Define request parameters for the LOCA dataset
          if scen == 'historical':
            input_dict = {
                "bbox": bbox,  # Bounding box for spatial selection
                "sdate": "1950-01-01",  # Start date of historical period
                "edate": "2005-12-31",  # End date of historical period
                "grid": "loca:" + str(models[i]) + ":rcp85",  # LOCA model with RCP8.5 scenario
                "elems": [{"name": "pcpn", "interval": [0, 0, 1]}]  # Request daily precipitation
            }
          elif scen == 'rcp45':
            input_dict = {
                "bbox": bbox,  # Bounding box for spatial selection
                "sdate": "2006-01-01",  # Start date 
                "edate": "2099-12-31",  # End date 
                "grid": "loca:" + str(models[i]) + ":rcp45",  # LOCA model with RCP8.5 scenario
                "elems": [{"name": "pcpn", "interval": [0, 0, 1]}]  # Request daily precipitation
            }
          elif scen == 'rcp85':
            input_dict = {
                "bbox": bbox,  # Bounding box for spatial selection
                "sdate": "2006-01-01",  # Start date of historical period
                "edate": "2099-12-31",  # End date of historical period
                "grid": "loca:" + str(models[i]) + ":rcp85",  # LOCA model with RCP8.5 scenario
                "elems": [{"name": "pcpn", "interval": [0, 0, 1]}]  # Request daily precipitation
            }
  
          # Send request to the RCC-ACIS API to retrieve precipitation data
          req = requests.post("http://grid2.rcc-acis.org/GridData", json=input_dict)
          data_vals = req.json()  # Convert response to JSON format
          data = data_vals['data']  # Extract the precipitation data section
  
          # Loop through all available time steps and extract precipitation values
          for k in range(len(data)):
              precip[:, :, k] = np.array(data[k][1])  # Store precipitation for grid cells
  
          # Replace negative values (which indicate missing data) with NaN
          precip[precip < 0] = np.nan
  
          # Compute the mean precipitation over the 8x8 grid for each time step
          all_precip[scen][i, j, :] = np.nanmean(precip, axis=(0, 1))
  #Save All Precipitation to a pickle for future access
  with open('all_precip.pkl', 'wb') as f:
      pickle.dump(all_precip, f)
  return

def LOCA_IDF():
  # =============================================================================
  # LOCA - Compute IDF Thresholds for Each Locations
  # =============================================================================
  
  scenarios_loca = ['historical', 'rcp45', 'rcp85']
  
  index_levels = [models, scenarios_loca]  # Replace with your actual model and scenario values
  index_names = ['model', 'scenario']
  
  multi_index = pd.MultiIndex.from_product(index_levels, names=index_names)
  loca_annual_totals = pd.DataFrame(index=multi_index, columns=range(100))
  
  loca_thresholds_cp = {}
  loca_change_factor = {}
  loca_cp_idf = {}
  
  j=0 #Station Number
  
  for scen in ['historical', 'rcp45', 'rcp85']:
    i=0 #Model Counter
    loca_thresholds_cp[scen]={}
    loca_change_factor[scen] = {}
    loca_cp_idf[scen]={}
    for model in models:
      current_precip = all_precip[scen][i, j, :]
      #Log annual precipitation for each model, i, and each scenario, j
      annual_temp = annual_totals(current_precip)
      loca_annual_totals.loc[(model, scen), :] = np.pad(annual_temp, (0, 100-len(annual_temp)), 'constant', constant_values=np.nan)
      #Log empirical thresholds
      loca_thresholds_cp[scen][model] = calc_rp_values(current_precip, type = 'empirical')
      # Log Change Factors
      loca_change_factor[scen][model] = loca_thresholds_cp[scen][model]/loca_thresholds_cp['historical'][model]
      # Log Adjusted Atlas14 Curve
      loca_cp_idf[scen][model] = atlas14[['2', '5', '10', '25', '50', '100']]*np.nanmean(loca_change_factor[scen][model],axis=0)
  
      i=i+1
    with open('loca_thresholds_cp.pkl', 'wb') as f:
        pickle.dump(loca_thresholds_cp, f)
    with open('loca_change_factor.pkl', 'wb') as f:
        pickle.dump(loca_change_factor, f)
    with open('loca_cp_idf.pkl', 'wb') as f:
        pickle.dump(loca_cp_idf, f)
