def add_order(ax,num,x=0, y=0):
    '''
    add picture order
    input:
        ax: figure axis
        num: figure number, 0-->(a), 1-->(b)
        x, y: mannuly define order position
    '''
    x0,x1=ax.get_xlim()
    y0,y1=ax.get_ylim()
    text='('+chr(97+num)+')' 
    if x == 0: 
        x=x0-(x1-x0)*0.02
    if y == 0:
        y=y1+(y1-y0)*0.05
    ax.text(x,y,text,size=14,weight='semibold')
    

def print_lat(lat,dim=2):
    '''
    print from number to latitude
    '''
    if lat>0:
        return '{}°N'.format(lat)
    elif lat<0:
        return '{}°S'.format(-lat)
    else:
        return '0°'


def plot_p(x,y,p_values,threshold = 0.05,**kw):
    '''
    plot dots over the values where is significant (p_values<threshold)
    input:
        x: x dim
        y: y dim
        p_values: 2-d array of p_values
        threshold: threshold to plot out the significance
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    s = p_values.shape
    nx,ny = len(x),len(y)
    if s[0]!=nx:
        p_values = p_values.T
    sig_x = x[:,np.newaxis].repeat(ny,axis=1)
    sig_y = y[np.newaxis,:].repeat(nx,axis=0)
    sigs = np.zeros(p_values.shape)
    sigs[p_values<threshold] = 1
    i_sig = np.where(sigs.flatten())
    plt.plot(sig_x.flatten()[i_sig],sig_y.flatten()[i_sig],'.',
        ms=0.8,**kw)

    
def get_color(color_type='my',n = 10):
    '''
    return a list of color hex 
    ''' 
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    import numpy as np
    cm_dict = {'Reds':cm.Reds,'Blues':cm.Blues,'rainbow':cm.rainbow,'bwr':cm.bwr,
              'GnBu':cm.GnBu,'coolwarm':cm.coolwarm,'gist_ncar':cm.gist_ncar,
              'seismic':cm.seismic,'Pastel1':cm.Pastel1}
    if color_type=='default':
        return plt.rcParams['axes.prop_cycle'].by_key()['color']
    elif color_type=='like_default':
        return ['royalblue','darkorange','forestgreen','firebrick','mediumpurple']
    elif color_type=='my':
        return ['royalblue','firebrick','forestgreen','darkorange','mediumpurple']
    elif color_type in cm_dict:
        cmap = cm_dict[color_type]
        var = range(n)
        colors_array = cmap((var-np.min(var))/(np.nanmax(var)-np.nanmin(var)))
        rainbow = [colors.rgb2hex(i) for i in colors_array]
        return rainbow 

def add_z_axis(ticks=range(3,-10,-1),ylabel=True,fontsize=13):
    '''
    when the figure y axis is pressure levels, add altitude as a secondary yaxis 
    input: ticks: your desire ticks in pressure levels
    ''' 
    import numpy as np
    import matplotlib.pyplot as plt
    def p2z(p):
        # pressure to altitude estimation
        return -7*np.log(p/1000)
    def z2p(z):
        return 1000*np.exp(-z/7)

    secax = plt.gca().secondary_yaxis('right', functions=(p2z,z2p))
    secax.minorticks_off()
    secax.set_yticks([int(p2z(10**(i))) for i in ticks],
                    [str(int(p2z(10**(i)))) for i in ticks]);
    if ylabel:
        secax.set_ylabel('Altitude (km)',fontsize=fontsize)
        

def mon_arr(start_year,end_year):
    '''
    get an array of month string 
    input:
        start_year: start year
        end_year: end year of the data +1
    output:
        array of month string 
    '''
    import numpy as np
    return np.array(['{}-{}-01'.format(str((i+start_year)),
                str(j+1).zfill(2)) for i in range(end_year-start_year) for j in range(12)],
                    dtype='datetime64')