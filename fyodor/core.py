from netCDF4 import Dataset                                             
import os                                                               
import numpy as np                                                    
import time
import matplotlib.pyplot as plt  
import matplotlib.patches as mpatches
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
# from astropy.utils.data import clear_download_cache
from astropy.constants import R_earth

__all__ = ['pwv']


def load_files():
    # Access working directory
    os.chdir('/Users/bmmorris/Downloads/')
    nc_filesT = glob.glob('*OR_ABI-L2-LVTPF*')
    nc_filesT = sorted(nc_filesT)
    nc_filesM = glob.glob('*OR_ABI-L2-LVMPF*')
    nc_filesM = sorted(nc_filesM)
    print('Number of Temperature files:', len(nc_filesT))
    print('Number of Moisture files:', len(nc_filesM))

    # Open the first file to retrieve earth parameters
    Proj_info = Dataset(nc_filesT[0], 'r')
    proj_info = Proj_info.variables['goes_imager_projection']
    lon_origin = proj_info.longitude_of_projection_origin
    H = proj_info.perspective_point_height + proj_info.semi_major_axis
    r_eq = proj_info.semi_major_axis
    r_pol = proj_info.semi_minor_axis

    # Retrieve pressure data
    P = Proj_info.variables['pressure'][:]
    Proj_info.close()

    # Retrieve time data
    g16_data_file = []
    t = []
    epoch = []
    date = []
    day = []

    for i in range(0, len(nc_filesT)):
        g16_data = nc_filesT[i]
        g16_data_file.append(g16_data)
        g16 = Dataset(g16_data_file[i], 'r')
        ttemp = g16.variables['t'][:]
        t.append(ttemp)
        epochtemp = 946728000 + int(t[i])
        epoch.append(epochtemp)
        datetemp = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(epoch[i]))
        date.append(datetemp)
        daytemp = time.strftime("%d-%m-%Y", time.gmtime(epoch[i]))
        day.append(daytemp)

    # Use astropy.time to keep format for target coordinates:
    times = Time(date, format='iso', scale='utc')

    # Barometric formula
    p0 = P[0]  # hPa
    R_D = 287  # Jkg-1K-1
    g = 9.81  # m/s2
    T_s = 288  # K

    H_sh = float((R_D * T_s) / g)  # Scale height
    h = H_sh * np.log((p0 / P)) * u.m

    e = np.sqrt((r_eq ** 2 - r_pol ** 2) / (r_eq ** 2))  # Eccentricity

    return (times, day, epoch, date, nc_filesT, nc_filesM, h, e, r_pol, r_eq, P,
            H, lon_origin, g16_data_file)


def pwv(location, Ra=None, Dec=None, P_minb=300, P_maxb=750,
        line_of_site='zenith'):
    """

    location :
    Ra :
    Dec :
    P_minb :
    P_maxb :
    P :
    line_of_site :  {"target", "zenith"}

    """
    params = load_files()
    times, day, epoch, date, nc_filesT, nc_filesM, h, e, r_pol, r_eq, P, H, lon_origin, g16_data_file = params

    # Pressure level boundaries
    P_minb = np.abs(P-P_minb).argmin()
    P_maxb = np.abs(P-P_maxb).argmin()

    # Convert from radian to degrees:
    raddeg = 180/np.pi

    Ra = float(Ra)
    # print('\n' 'Declination:')
    Dec = float(Dec)
    Sky = SkyCoord(ra=Ra*u.degree, dec=Dec*u.degree)
    Aa = Sky.transform_to(AltAz(obstime=times,
                                location=location.to_EarthLocation()))
    latt = Dec
    lont = Ra

    # print('Please wait...')

    # Computes PWV along line of sight
    if line_of_site == 'target':
        INDEX = np.ravel(np.where(Aa.alt.degree<30))
        INDEXP = np.ravel(np.where(Aa.alt.degree>30))

        # Keep time values corresponding to Alt above 30 degrees
        EPOCH = epoch
        for index in sorted(INDEX, reverse=True):
            del EPOCH[index]
        DATE = []
        for i in range(0, len(epoch)):
            DATEt = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(EPOCH[i]))
            DATE.append(DATEt)

        Alt = Aa.alt.rad
        Az = Aa.az.rad

        # Compute distance from location to projection point
        delta_x = []
        d_lat = []
        d_lon = []

        for i in range(0,len(Alt)):
            delta_xi = []
            for j in h:
                delta_xt = j/np.tan(Alt[i])
                delta_xt = delta_xt*u.m**-1
                delta_xi.append(delta_xt)
            delta_x.append(delta_xi)
        delta_x = np.array(delta_x)

        for i in range(0,len(Az)):
            d_latt = delta_x[i,]*np.cos(Az[i])
            d_lont = delta_x[i,]*np.sin(Az[i])
            d_lat.append(d_latt)
            d_lon.append(d_lont)
        d_lat = np.array(d_lat)
        d_lon = np.array(d_lon)

        # Compute latitude and longitude of projection points
        lat_proj = []
        lon_proj = []
        for i in range(0, len(Alt)):
            obs_latt = location.latdeg + raddeg*(np.arctan(d_lat[i,]/R_earth*u.m**1)*u.rad**-1)
            obs_lont = location.londeg + raddeg*(np.arctan(d_lon[i,]/R_earth*u.m**1)*u.rad**-1)
            lat_proj.append(obs_latt)
            lon_proj.append(obs_lont)
        lat_proj = np.array(lat_proj)
        lon_proj = np.array(lon_proj)
        rad = (np.pi)/180

        lat_proj_rad = rad*lat_proj
        lon_proj_rad = rad*lon_proj
        lambda_0 = rad*lon_origin

        #T ransform into scan angles
        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(lat_proj_rad))

        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H -r_c*np.cos(lat_origin)*np.cos(lon_proj_rad-lambda_0)
        s_y = -r_c*np.cos(lat_origin)*np.sin(lon_proj_rad-lambda_0)
        s_z = r_c*np.sin(lat_origin)

        s = np.sqrt(s_x**2+s_y**2+s_z**2)

        x = np.arcsin(-s_y/s)
        y = np.arctan(s_z/s_x)

        g16_data_fileT = []
        xscanT = []
        yscanT = []
        XT = []
        YT = []
        LVT = []

        # Retrieve Temperature data
        for i in INDEXP:
            g16_data = nc_filesT[i]
            g16_data_fileT.append(g16_data)
            g16 = Dataset(nc_filesT[i], 'r')
            xtemp = g16.variables['x'][:]
            xscanT.append(xtemp)
            ytemp = g16.variables['y'][:]
            yscanT.append(ytemp)

            LVTi = []
            Xi = []
            Yi = []
            for j in range(P_minb, P_maxb+1):
                Xtemp = np.abs( xtemp-x[i,j]).argmin()
                Xi.append(Xtemp)
                Ytemp = np.abs( ytemp-y[i,j]).argmin()
                Yi.append(Ytemp)
                LVTtemp = g16.variables['LVT'][Ytemp, Xtemp, j]
                LVTi.append(LVTtemp)
            LVT.append(LVTi)
            XT.append(Xi)
            YT.append(Yi)
        LVT = np.array(LVT)

        # Retrieve Relative humidity data
        g16_data_fileM = []
        xscanM = []
        yscanM = []
        XM = []
        YM = []
        LVM = []

        for i in INDEXP:
            g16_dataM = nc_filesM[i]
            g16_data_fileM.append(g16_dataM)
            g16M = Dataset(nc_filesM[i], 'r')
            xtempM = g16M.variables['x'][:]
            xscanM.append(xtempM)
            ytempM = g16M.variables['y'][:]
            yscanM.append(ytempM)

            LVMi = []
            Xi = []
            Yi = []
            for j in range(P_minb, P_maxb+1):
                XtempM = np.abs( xtempM-x[i,j]).argmin()
                Xi.append(XtempM)
                YtempM = np.abs( ytempM-y[i,j]).argmin()
                Yi.append(YtempM)
                LVMtemp = g16M.variables['LVM'][YtempM, XtempM, j]
                LVMi.append(LVMtemp)
            LVM.append(LVMi)
            XM.append(Xi)
            YM.append(Yi)
        LVM = np.array(LVM)

        P = 100*P
        LVT = LVT-273.15

        Pi = P[P_minb:P_maxb+1]

        # Constants needed and integrand
        rho_w = 1000 # kg/m**3
        g = 9.81 #m/s**2
        C = (-1)/(rho_w*g)
        ev = 100*6.11*LVM*10**((7.5*LVT)/(LVT+237.15)) # Partial water vapour pressure in Pa
        q = (0.622*ev)/(Pi-0.378*ev) #Specific humdity
        f = 1000*C*q #Complete integrand multiplied by 1000 to get the PWV in mm.

        # Numerical integration
        PWV = []
        for j in range(0, len(LVT)):
            integral = 0
            for i in range(1, len(Pi)):
                integral = integral +(Pi[i]-Pi[i-1])*((f[j,i]+f[j,i-1])/2)
            PWV.append(integral)

        PWV = np.asarray(PWV)

        print(PWV)

        PWV = np.nan_to_num(PWV)

        # Plot and save data
        fig = plt.figure(figsize=(8,5))
        ax=fig.add_subplot(111)
        ax.plot(DATE, PWV, 'bo', ms=4)
        plt.title('Precipitable Water Vapor along line of sight, {} on {}'
                  .format(location.site, day[0]), fontsize=20)
        plt.xticks(rotation='vertical', fontsize=16)
        plt.yticks(fontsize=16)
        ax.set_xlabel("Date", color="C0", fontsize =18)
        ax.set_ylabel("PWV (mm)", color="C0", fontsize =18)
        RA_patch = mpatches.Patch(color='white', label='RA: {} degrees'.format(Ra))
        Dec_patch = mpatches.Patch(color='white', label='Dec: {} degrees'.format(Dec))
        every_nth = 4
        for n, label in enumerate(ax.xaxis.get_ticklabels()):
            if n % every_nth != 0:
                label.set_visible(False)
        for n, label in enumerate(ax.xaxis.get_ticklines()):
            if n % every_nth != 0:
                label.set_visible(False)
        plt.tight_layout()
        plt.legend(handles=[RA_patch, Dec_patch], loc='lower right', fontsize=18)
        plt.show()
        fig.savefig('PWV_line_of_sight_{}_{}.png'.format(location.site, day[0]))

        np.savetxt('PWV_line_of_sight_{}_{}.csv'.format(location.site, day[0]),
                   np.column_stack((DATE, PWV)),
                   delimiter=',', fmt='%s', header='Time, PWV')

    # Computes PWV at zenith
    elif line_of_site == 'zenith':

        # Transform latitude and longitude into scan angles
        rad = (np.pi)/180
        lambda_0 = rad*lon_origin
        obs_lat_rad = rad*latt
        obs_lon_rad = rad*lont

        lat_origin = np.arctan(((r_pol**2)/(r_eq**2))*np.tan(obs_lat_rad))

        r_c = r_pol/(np.sqrt(1-(e**2)*(np.cos(lat_origin))**2))

        s_x = H -r_c*np.cos(lat_origin)*np.cos(obs_lon_rad-lambda_0)
        s_y = -r_c*np.cos(lat_origin)*np.sin(obs_lon_rad-lambda_0)
        s_z = r_c*np.sin(lat_origin)

        s = np.sqrt(s_x**2+s_y**2+s_z**2)

        x = np.arcsin(-s_y/s)
        y = np.arctan(s_z/s_x)

        xscanT = []
        yscanT = []

        # Retrieve Temperature data
        LVT = []
        for i in range(0, len(nc_filesT)):
            g16_data = nc_filesT[i]
            g16_data_file.append(g16_data)
            g16 = Dataset(g16_data_file[i], 'r')
            xtemp = g16.variables['x'][:]
            xscanT.append(xtemp)
            ytemp = g16.variables['y'][:]
            yscanT.append(ytemp)

            XT = []
            YT = []
            LVTi = []
            for j in range(0, len(P)):
                Xtemp = np.abs( xtemp-x).argmin()
                XT.append(Xtemp)
                Ytemp = np.abs( ytemp-y).argmin()
                YT.append(Ytemp)
                LVTtemp = g16.variables['LVT'][Ytemp, Xtemp, j]
                LVTi.append(LVTtemp)
            LVT.append(LVTi)

        LVT = np.array(LVT)

        # Retrieve Relative humidity data
        g16_data_fileM = []
        g16ncM = []
        xscanM = []
        yscanM = []
        LVM = []

        for i in range(0, len(nc_filesM)):
            g16_dataM = nc_filesM[i]
            g16_data_fileM.append(g16_dataM)
            g16M = Dataset(g16_data_fileM[i], 'r')
            xtempM = g16M.variables['x'][:]
            xscanM.append(xtempM)
            ytempM = g16M.variables['y'][:]
            yscanM.append(ytempM)

            XM = []
            YM = []
            LVMi = []
            for j in range(0, len(P)):
                XtempM = np.abs( xtempM-x).argmin()
                XM.append(XtempM)
                YtempM = np.abs( ytempM-y).argmin()
                YM.append(YtempM)
                LVMtemp = g16M.variables['LVM'][YtempM, XtempM, j]
                LVMi.append(LVMtemp)
            LVM.append(LVMi)

        LVM = np.array(LVM)

        # Change pressure units to Pa and Temperature to K
        P = 100*P
        LVT = LVT-273.15

        # Constants needed and integrand
        rho_w = 1000 # kg/m**3
        g = 9.81 #m/s**2
        C = (-1)/(rho_w*g)
        ev = 100*6.11*LVM*10**((7.5*LVT)/(LVT+237.15)) # Partial water vapour pressure in Pa
        #ev = 100*6.094*LVM*np.exp(17.625*LVT/(LVT+243.04))
        q = (0.622*ev)/(P-0.378*ev) #Specific humdity
        f = 1000*C*q #Complete integrand multiplied by 1000 to get the PWV in mm.

        # Numerical integration
        PWV = []
        for j in range(0, len(nc_filesT)):
            integral = 0
            for i in range(P_minb+1, P_maxb+1):
                integral = integral +(P[i]-P[i-1])*((f[j,i]+f[j,i-1])/2)
            PWV.append(integral)

        PWV = np.asarray(PWV)

        # print(PWV)

        PWV = np.nan_to_num(PWV)

        # Plot and save data
        fig = plt.figure(figsize=(8,5))
        ax=fig.add_subplot(111)
        ax.plot(date, PWV, 'bo', ms=4)
        plt.title('Precipitable Water Vapor at zenith, {} on {}'
                  .format(location.site, day[0]), fontsize=20)
        plt.xticks(rotation='vertical', fontsize=16)
        plt.yticks(fontsize=16)
        ax.set_xlabel("Date", color="C0", fontsize =18)
        ax.set_ylabel("PWV (mm)", color="C0", fontsize =18)
        every_nth = 4
        for n, label in enumerate(ax.xaxis.get_ticklabels()):
            if n % every_nth != 0:
                label.set_visible(False)
        for n, label in enumerate(ax.xaxis.get_ticklines()):
            if n % every_nth != 0:
                label.set_visible(False)
        plt.tight_layout()
        plt.show()
        fig.savefig('PWV_at_zenith_{}_{}.png'.format(location.site, day[0]))

        np.savetxt('PWV_at_zenith_{}_{}.csv'.format(location.site, day[0]),
                   np.column_stack((date, PWV)),
                   delimiter=',', fmt = '%s', header= 'Time, PWV')
