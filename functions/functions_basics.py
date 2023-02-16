def geo_avg(x,lat,dim=2):
    '''
    geo_avg: to calculate weighting average according to latitude
    input: 
        x: variable 
        lat: corresponding latittude
        dim: the order of the lat dimension, two cases: 2:[time,lev,lat,*lon],or 1:[time or lev, lat, *lon]
    output:
        result: 1d or 2d average result 
    '''
    import numpy as np
    s = x.shape
    if ((len(s)==4) & (dim==2)) or ((len(s)==3) & (dim==1)):
        x = np.nanmean(x,axis=-1)
    coslat = np.cos(lat/180*np.pi)
    s = x.shape
    if len(s)==3:
        result = np.nanmean(x*coslat[np.newaxis,np.newaxis,:],axis=-1)/np.nanmean(coslat)
    if len(s)==2:
        result = np.nanmean(x*coslat[np.newaxis,:],axis=-1)/np.nanmean(coslat)
    return result

def cal_anomaly(x):
    '''
    calculate anomaly of a numpy array 
    input: x: 1-d,2-d,3-d or 4d numpy array, !!! the first dimension must be month 
    output: x with seasonal cycle removed 
    '''
    import numpy as np
    s = x.shape
    n_time = s[0]
    monthly_mean = np.nanmean(x.reshape([n_time//12,12,*s[1:]]),axis=0).\
        reshape([1,12,*s[1:]]).repeat(len(x)//12,axis=0).reshape(s)
    return x-monthly_mean

def select_month(x,target_mon):
    '''
    select month or season from a monthly time series
    input: 
        x: array, 1,2,3,4 dimension
        target_mon: 
            1. number of month, from 1-12 
            2. name of month, e.g. Jan, Feb
            3. season name: DJF: 1,2,12; JJA: 6,7,8 SON: 9,10,11, MAM: 3,4,5
            4. phase name: dry: 1,2,3,12; wet: 6,7,8,9
    output: 
        array with month selected or seasonal mean 
    '''
    s = x.shape
    n_mon = s[0]
    if type(target_mon) != str:
        i_mon = [i for i in range(n_mon) if i%12 == target_mon-1]
        return x[i_mon]
    else:
        import numpy as np
        from datetime import datetime,timedelta
        mon_name_list = [(datetime(2000,1,1)+timedelta(days=31*i)).strftime("%b") for i in range(12)]
        mon_dict = {mon_name_list[i]:i for i in range(12)}
        season_dict = {'DJF':[0,1,11],'JJA':[5,6,7],'SON':[8,9,10],'MAM':[2,3,4]}
        phase_dict = {'dry':[0,1,2,11],'wet':[5,6,7,8]}
        
        if target_mon in mon_dict:
            i_mon = [i for i in range(n_mon) if i%12 == mon_dict[target_mon]]
            return x[i_mon]
        elif target_mon in season_dict:
            i_mon = [i for i in range(n_mon) if i%12 in season_dict[target_mon]]
            x_mon = x[i_mon]
            if target_mon == 'DJF':
                x_mon = np.append(np.nan,x_mon[:-1])
            return np.nanmean(x_mon.reshape([s[0]//12,3,*s[1:]]),axis=1)
        else:
            i_mon = [i for i in range(n_mon) if i%12 in phase_dict[target_mon]]
            x_mon = x[i_mon]
            if target_mon == 'dry':
                x_mon = np.append(np.nan,x_mon[:-1])
            return np.nanmean(x_mon.reshape([s[0]//12,4,*s[1:]]),axis=1)

def normalize(x):
    '''
    function to normalize data 
    '''
    import numpy as np
    return (x-np.nanmean(x))/np.nanstd(x)

def find_index(arr,target,method='nearest'):
    '''
    find an index of target value from amonotonous 1-d array arr
    '''
    import numpy as np
    if method == 'nearest':
        return (np.abs(arr - target)).argmin()
    else:
        if arr[1]<arr[0]:  ## if x is a decreasing array, reverse 
            arr = arr[::-1]  
        if method == 'higher':
            return np.where(arr>=target)[0][0]
        if method == 'lower':
            return np.where(arr<=target)[0][-1]
        
    
def moving_average(arr,n,method = 'nan'):
    '''
    calculate moving average values of 1-d array, and return an array with the same length 
    input:
        arr: 1-d array 
        n: moving window length 
        method:
            nan: fill in nan 
            avg: average from 0-1, 0-2, 0-3 ...
            diff: only use this when calculate annual mean, n = 13
    '''
    import numpy as np
    def moving_average_center(a, n) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    l1 = n//2-1
    l2 = n-l1
    l = len(arr)
    arr_new = np.zeros(l)
    if method == 'nan':
        arr_new[:l1] = np.nan
        arr_new[l1:l-l2+1] = moving_average_center(arr, n)
        arr_new[l-l2+1:] = np.nan
    if method == 'avg':
        for i in range(l1):
            arr_new[i] = np.nanmean(arr[:i+1])
        for i in range(l2):
            arr_new[-i-1] = np.nanmean(arr[-i-1:])
        arr_new[l1:l-l2+1] = moving_average_center(arr, n)
    if method == 'diff' and n==13:
        a2 = moving_average_center(arr, n)
        diff = (arr[l1:l-l2+1]-a2).reshape([(len(arr)-n+1)//12,12]).mean(axis=0)  # monthly mean difference between arr and running mean
        a1 = arr[:6] - diff[6:]
        a12 = np.append(a1,a2)
        a3 = arr[-6:] - diff[:6]
        arr_new = np.append(a12,a3)
    return arr_new

def convert_cftime_to_int(t):
    '''
    convert cftime to integer 
    input:
        arr: 1-d array 
        n: moving window length 
        method:
            nan: fill in nan 
            avg: average from 0-1, 0-2, 0-3 ...
            diff: only use this when calculate annual mean, n = 13
    '''
    from datetime import datetime
    return int(datetime.strftime(datetime.strptime(t.isoformat(),'%Y-%m-%dT%H:%M:%S'),
                 '%Y%m%d'))

def get_lat_lim(lat,lat_min,lat_max):
    '''
    calculate a range of latitude, in both hemispheres
    '''
    import numpy as np
    i_lat_n = np.where((lat>=lat_min) & (lat<=lat_max))[0]
    i_lat_s = np.where((lat<=-lat_min) & (lat>=-lat_max))[0]
    i_lats = [i_lat_s,i_lat_n]
    return i_lats
