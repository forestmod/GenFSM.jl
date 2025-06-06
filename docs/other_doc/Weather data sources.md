# WorldClimate

zip of tif images

## Historic:
1960-2018, monthly
2.5 minutes
tmin,tmax,prec

## Future
monthly 2021-2100 by (not moving averages) 20 years periods
2.5 minutes
tmin,tmax,prec
23 models, 4 SSPs



TO Get:

Scenarios:
ssp126 SSP1-RCP2.6 climate as simulated by the GCMs.
ssp370 SSP3-RCP7 climate as simulated by the GCMs.
ssp585 SSP5-RCP8.5 climate as simulated by the GCMs.

Variables:
bio{1|19} (yearly), pr,tas,tasmax,tasmin (monthly)

Model:
MPI-ESM1-2-LR
(in second order: IPSL-CM6A-LR)

activity_id='ScenarioMIP', 
             table_id='Amon', 
             experiment_id='ssp585', 
             institution_id='MPI-M', 
             source_id='MPI-ESM1-2-LR', 
             member_id='r1i1p1f1',


Historis: 1999:2014 data from C/historic



### Chelsa availability

Historical
monthly pr/tas variables: 1980 - 2018
https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.zhdk.cloud.switch.ch%2Fchelsav2%2F&prefix=GLOBAL%2Fmonthly%2F
Or 
https://envicloud.wsl.ch/#/?bucket=https%3A%2F%2Fos.zhdk.cloud.switch.ch%2Fchelsav2%2F&prefix=%2F


https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/tasmax/CHELSA_tasmax_01_1999_V.2.1.tif
- 

tas/tas_min,tas_max in historical are in K/10 units, in scenarios are in K
pr in obs are in mm month/100, in the projection are in mm month (more precisely mm month are kg/m² month)


What I need ?

Observed data at high resolution level (2x2km) for EU for specific time frames (1999-2023)
I need the same variables in the future