RCM_DICT = {
    'EUR-11': 
    { 
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-6',
    },
    'EUR-22': 
    {
        'REMO': 'GERICS_REMO2015',
    },
    'WAS-22': {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'EAS-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-6',
    },
    'CAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'SAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'NAM-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-6',
    },
    'AUS-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'AFR-22':
     {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
    'SEA-22':
    {
        'REMO': 'GERICS_REMO2015',
        'RegCM': 'ICTP_RegCM4-7',
    },
}


MODEL_DICT={
    'REMO' : dict(sftuf='orig_v3', orog='orog',sftlf='sftlf'),
    'RegCM' : dict(sftuf='', orog='orog',sftlf='sftlf'),
}

# Dictionary containing city locations and their respective domains
#location = {
#     'Mexico City' : dict(lon=-99.0833, lat=19.4667, domain = 'CAM-22'),
#     'Buenos Aires' : dict(lon=-58.416, lat=-34.559, domain = 'SAM-22'),
#     'New York' : dict(lon=-74.2261, lat=40.8858, domain = 'NAM-22'),
#     'Sydney' : dict(lon=151.01810, lat=-33.79170, domain = 'AUS-22'),
#     'Beijing' : dict(lon=116.41, lat=39.90, domain = 'EAS-22'),
#     'Tokyo' : dict(lon = 139.84, lat = 35.65, domain = 'EAS-22'),
#     'Jakarta' : dict(lon = 106.81, lat = -6.2, domain = 'SEA-22'), 
#     'Johannesburg' : dict(lon=28.183, lat=-25.733, domain = 'AFR-22'),
#     'Riyadh' : dict(lon=46.73300, lat=24.7000, domain = 'WAS-22'),
#     'Berlin' : dict(lon=13.4039, lat=52.4683, domain = 'EUR-11'),
#     'Paris' : dict(lon=  2.35, lat=48.85, domain = 'EUR-11'),
#     'London' : dict(lon= -0.13, lat=51.50, domain = 'EUR-11'),
#     'Madrid' : dict(lon= -3.70, lat=40.42, domain = 'EUR-11'),
#     'Los Angeles': dict(lon = -118.24, lat = 34.05, domain = 'NAM-22'),
#     'Montreal': dict(lon = -73.56, lat = 45.50, domain = 'NAM-22'),
#     'Chicago': dict(lon = -87.55, lat = 41.73, domain = 'NAM-22'),
#     'Bogota': dict(lon = -74.06, lat = 4.62, domain = 'SAM-22'),
#     'Baghdad': dict(lon = 44.40, lat = 33.34, domain = 'WAS-22'),
#     'Tehran': dict(lon = 51.42, lat = 35.69, domain = 'WAS-22'),
#     'Tashkent': dict(lon = 69.24, lat = 41.31, domain = 'WAS-22'),
#     'Cairo': dict(lon = 31.25, lat = 30.06, domain = 'AFR-22'),
#     'Delhi [New Delhi]': dict(lon = 77.22, lat = 28.64, domain = 'WAS-22'),
#    'Barcelona': dict(lon = 2.18, lat = 41.39, domain = 'EUR-11'),
#    'Rome': dict(lon =  12.51, lat = 41.89, domain = 'EUR-11'),
#    'Athens': dict(lon =   23.72, lat =  37.98, domain = 'EUR-11'),
#}

cities = {
     'Jakarta': dict(city = 'Jakarta', lon = 106.83, lat = -6.18, domain = 'SEA-22', domain_lim = 1, obs_im = 0.5),
     'Quezon City [Manila]': dict(city = 'Manila', lon = 120.98, lat = 14.60, domain = 'SEA-22', domain_lim = 1, obs_im = 0.5),
     'Mumbai': dict(city = 'MumbaiBombay', lon = 72.82, lat = 18.96, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Dhaka': dict(city = 'Dhaka', lon = 90.41, lat = 23.70, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Delhi [New Delhi]': dict(city = 'NewDelhi', lon = 77.22, lat = 28.64, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Tehran': dict(city = 'Tehran', lon = 51.42, lat = 35.69, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Cairo': dict(city = 'Cairo', lon = 31.25, lat = 30.06, domain = 'AFR-22', domain_lim = 1, obs_im = 0.5), 
     'New York' : dict(city = 'NewYork', lon=-74.2261, lat=40.8858, domain = 'NAM-22', domain_lim = 0.75, obs_im = 0.5), 
     'Tokyo' : dict(city = 'Tokyo', lon = 139.84, lat = 35.65, domain = 'EAS-22', domain_lim = 1, obs_im = 0.5),
     'Buenos Aires' : dict(city = 'BuenosAires', lon=-58.416, lat=-34.559, domain = 'SAM-22', domain_lim = 1, obs_im = 0.5),
     'São Paulo' : dict(city = 'SaoPaulo', lon=-46.633, lat=-23.550, domain = 'SAM-22', domain_lim = 1, obs_im = 0.5),
     'Los Angeles': dict(city = 'LosAngeles', lon = -118.24, lat = 34.05, domain = 'NAM-22', domain_lim = 1, obs_im = 0.5),
     'Beijing' : dict(city = 'Beijing', lon=116.41, lat=39.90, domain = 'EAS-22', domain_lim = 1, obs_im = 0.5), 
     'Chengdu' : dict(city = 'Chengdu', lon=104.07, lat=30.67, domain = 'EAS-22', domain_lim = 1, obs_im = 0.5), 
     'Mexico City' : dict(city = 'MexicoCity', lon=-99.0833, lat=19.4667, domain = 'CAM-22', domain_lim = 1, obs_im = 0.5),
     'Johannesburg' : dict(city = 'Johannesburg', lon=28.183, lat=-25.733, domain = 'AFR-22', domain_lim = 1, obs_im = 0.5),
     'Chicago': dict(city = 'Chicago', lon = -87.55, lat = 41.73, domain = 'NAM-22', domain_lim = 0.75, obs_im = 0.5),
     'Montreal': dict(city = 'Montreal', lon = -73.56, lat = 45.50, domain = 'NAM-22', domain_lim = 1, obs_im = 0.5),
     'Seoul': dict(city = 'Seoul', lon = 126.98, lat = 37.57, domain = 'EAS-22', domain_lim = 1, obs_im = 0.5),
     'Lima': dict(city = 'Lima', lon = -77.03, lat = -12.04, domain = 'SAM-22', domain_lim = 1, obs_im = 0.5),
     'Lagos': dict(city = 'Lagos', lon = 3.38, lat = 6.52, domain = 'AFR-22', domain_lim = 1, obs_im = 0.5),
     'Luanda': dict(city = 'Luanda', lon = 13.24, lat = -8.81, domain = 'AFR-22', domain_lim = 1, obs_im = 0.5),
     'Riyadh' : dict(city = 'Riyadh', lon=46.73300, lat=24.7000, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Tashkent': dict(city = 'Tashkent', lon = 69.24, lat = 41.31, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Sydney' : dict(city = 'Sydney', lon=151.01810, lat=-33.79170, domain = 'AUS-22', domain_lim = 1, obs_im = 0.5),
     'Shanghai': dict(city = 'Shanghai', lon = 121.47, lat = 31.23, domain = 'EAS-22', domain_lim = 1, obs_im = 0.5),
     'Baghdad': dict(city = 'Baghdad', lon = 44.40, lat = 33.34, domain = 'WAS-22', domain_lim = 1, obs_im = 0.5),
     'Khartoum': dict(city = 'Khartoum', lon = 32.53, lat = 15.58, domain = 'AFR-22', domain_lim = 1, obs_im = 0.5),
     'Santiago': dict(city = 'SantiagoChile', lon = -70.65, lat = -33.46, domain = 'SAM-22', domain_lim = 1, obs_im = 0.5),
     'Melbourne': dict(city = 'Melbourne', lon = 144.96, lat = -37.81, domain = 'AUS-22', domain_lim = 1, obs_im = 0.5),
     'Singapore': dict(city = 'Singapore', lon = 103.85, lat = 1.29, domain = 'SEA-22', domain_lim = 1, obs_im = 0.5),
     'Bogota': dict(city = 'Bogota', lon = -74.06, lat = 4.62, domain = 'SAM-22', domain_lim = 1, obs_im = 0.5),

     ## European cities within Global Analysis (for REMO on 0.22 grid res as well)
     'Paris' : dict(city = 'Paris', lon=  2.35, lat=48.85, domain = 'EUR-11', domain_lim = 0.5, obs_im = 0.5),
     'London' : dict(city = 'London', lon= -0.13, lat=51.50, domain = 'EUR-11', domain_lim = 0.5, obs_im = 0.5),
     'Moscow' : dict(city = 'Moscow', lon=37.62, lat=55.75, domain = 'EUR-11', domain_lim = 0.5, obs_im = 0.5),
     'Istanbul': dict(city = 'Istanbul', lon=28.98, lat=41.01, domain = 'EUR-11', domain_lim = 1, obs_im = 0.5),
     'Berlin' : dict(city = 'Berlin', lon=13.4039, lat=52.4683, domain = 'EUR-11', domain_lim = 0.5, obs_im = 0.5),
     'Naples': dict(city = 'Naples', lon = 14.27, lat = 40.85, domain = 'EUR-11', domain_lim = 0.5, obs_im = 0.5),

     ## Additional cities for 0.11 - 0.22 comparison 
     'Barcelona': dict(city = 'Barcelona', lon = 2.18, lat = 41.39, domain = 'EUR-11', domain_lim = 1, obs_im = 0.5),
     'Athens': dict(city = 'Athens', lon = 23.72, lat = 37.98, domain = 'EUR-11', domain_lim = 1, obs_im = 0.5),
     'Prague': dict(city = 'Prague', lon = 14.42, lat = 50.08, domain = 'EUR-11', domain_lim = 1, obs_im = 0.5),

}
