# science_RODSEX/initial_data_processing
Processing raw ADV files for vorticity and other stats.

vort_data_1023.m
    Processes ring data for vorticity, potential vorticity, velocity, and pressure 
    time series (8Hz), and hourly u/v/p spectra saved in daily files.  
    Process hourly averages and standard deviation of vorticity, potential vorticity, 
    and velocity, hourly spectra and wave stats, and Lippmann rotation velocities. 
    
steve_data_0927.m and steve_data_1014_1015.m
    Process vorticity for Steve, text files:  [decimal day, vorticity]

vortex_lippmann_ratio.m
    Get Lippmann rotational velocities: hourly rms infragravity rotational from 
    excess uv variance as compared to expected wave pressure variance 
    (all depth corrected).