import matplotlib
matplotlib.use('Agg')

from pylab import *
rcParams['figure.facecolor']='white'
rcParams['font.size']=16
rcParams['legend.fontsize']=14

from matplotlib.backends.backend_pdf import PdfPages
pdf = PdfPages('ana_sim_plots.pdf')

energies = [1.e17,  3.e17,  1.e18, 3.e18,  1.e19, 3.e19,    1.e20, 3.e20, 1.e21]

#ANITA_liveTime = 
#Auger_liveTime = 
#IC_liveTime    = 
Auger    = [3.e15,  2.5e16, 1.e17,  2.e17,  3.e17,  4.e17,  4.e17, 4.e17, 4.e17]
IC_x_3   = [2.5e16, 5.5e16,   1.1e17, 2.5e17, 4.5e17, 7.5e17, 1.e18, 1.e18, 1.e18]
def results_plot(tag, num):
    A_Omega_start = []
    A_Omega_range = []
    A_Omega_exit  = []
    A_Omega_trig  = []

    for e in energies:
        e_str = '%1.0e'%e
        #e_str = e_str.replace('+','')
        #print e_str
        fnm = tag.replace('_X_', '_'+e_str+'_')
        print fnm
        f = np.load(fnm)
        A_Omega_start.append(f['A_Omega_start'])
        A_Omega_exit.append(f['A_Omega_exit'])
        A_Omega_range.append(f['A_Omega_range'])
        A_Omega_trig.append(f['A_Omega_trig'])
        print f['triggered_events'].shape
        if(f['triggered_events'].size ==0): continue

        #'''
        fig = figure(num, figsize=(12,10))
        subplot(223)
        plot(f['triggered_events'][:,0], 90.-f['triggered_events'][:,6], '.')
        xlabel('log10 ( Tau Exit Energy / eV )')
        ylabel('Tau Emergence Angle, deg')
        grid(True)
        xlim(16.,21.)
        ylim(0.,60.)
        subplot(221)
        hist(f['triggered_events'][:,0], bins=np.linspace(16.,21.,50))
        xlabel('log10 ( Tau Exit Energy / eV )')
        xlim(16.,21.)
        #ylabel('Tau Emergence Angle, deg')
        grid(True)
        suptitle(fnm)
        subplot(224)
        hist(90.-f['triggered_events'][:,6], bins=np.linspace(0.,60.,60))
        xlabel('Tau Emergence Angle, deg')
        xlim(0.,60.)
        grid(True)
        subplots_adjust(top=0.95)
        suptitle(fnm.split('/')[-1].split('.npz')[0])
        #pdf.savefig()
        fig.savefig('./ana_sim_figures/'+fnm.split('/')[-1].split('.npz')[0] + '.png')
        fig.clear()
        f.close()
        #'''


        '''
        figure()
        h,b = np.histogram(90.-f['triggered_events'][:,6], bins=np.linspace(-1.,90.,92))
        h = h.astype(float)
        h *= f['A_Omega_trig']/np.sum(h)
        plot(b[1:], h, drawstyle='steps', lw=2)
        grid(True)
        #hist(90.-f['triggered_events'][:,6], bins=np.linspace(0.,90.,91.))
        title(fnm)
        xlabel('Tau Exit Angle, deg')
        '''
        '''
        figure()
        h,b = np.histogram(f['triggered_events'][:,0], bins=np.linspace(15.,21.,61))
        h = h.astype(float)
        h *= f['A_Omega_trig']/np.sum(h)
        plot(b[1:], h, drawstyle='steps', lw=2)
        grid(True)
        #hist(90.-f['triggered_events'][:,6], bins=np.linspace(0.,90.,91.))
        title(fnm)
        xlabel('log10 ( Tau Exit Energy / eV )')
        '''


        

    A_Omega_start = np.array(A_Omega_start)
    A_Omega_range = np.array(A_Omega_range)
    A_Omega_exit  = np.array(A_Omega_exit)
    A_Omega_trig  = np.array(A_Omega_trig)

    fac = 1.e10 # km^2 in cm^2
    fac *= 17.*24*60*60. # 17 days in seconds
    '''
    figure()
    loglog(energies, A_Omega_start * fac, lw=2)
    plot(energies,   A_Omega_exit  * fac, lw=2)
    plot(energies,   A_Omega_range * fac, lw=2)
    plot(energies,   A_Omega_trig  * fac, lw=2)
    grid(True)
    '''
    return energies,  A_Omega_start * fac, A_Omega_exit  * fac, A_Omega_range * fac,  A_Omega_trig  * fac

#tag = 'ice_1.5km_alt_37km_X_regCS.npz'
#figure()
#energies_regCS, exp_start_regCS, exp_exit_regCS, exp_range_regCS, exp_trig_regCS = results_plot(tag)


colors = cm.hot(np.linspace(0., 1., 8))
model_count = 0
#for model in ['regCS_regEL', 'regCS_lowEL', 'lowCS_regEL', 'lowCS_lowEL']:
for model in ['lowCS_stdEL', 'lowCS_lowEL', 'midCS_stdEL', 'midCS_lowEL', 'uppCS_stdEL', 'uppCS_lowEL']:
    model_count+=1
    fig = figure(1000 + model_count,figsize=(6,10))
    ax = subplot(111)
    ax.set_yscale('log')
    ax.set_xscale('log')
    count = 0
    for th in [0.0, 1.0, 2.0, 3.0, 4.0]:
        tag = '/aurora_nobackup/eva/romerowo/nutau_outputs/20170508/ice_%1.1fkm_alt_37km_X_eV_%s.npz'%(th, model)
        energies, exp_start, exp_exit, exp_range, exp_trig = results_plot(tag, count+1)
        figure(1000 + model_count)
        loglog(energies, exp_trig, '-', lw=2, label='%1.1f km'%th, color=colors[count])
        count += 1

    #'''
    fig = figure(1000 + model_count)
    title(model.replace('_', ', '))
    if(model == 'midCS_stdEL'):
        plot(energies[:-2],  Auger[:-2],               'g--',  lw=4,  alpha=0.8, label='Auger (2015)')
        plot(energies[:-2],  np.array(IC_x_3[:-2]), 'b-',   lw=4,  alpha=0.8, label='IceCube (2016)' )
    grid(True, which='both')
    xlim(1.e17, 1.e21)
    ylim(1.e13, 1.e18)
    xlabel('Neutrino Energy, eV')
    ylabel('Exposure, cm$^2$ sr s')
    legend(loc=4, fontsize=14)
    subplots_adjust(bottom=0.125, left=0.2)
    pdf.savefig()
    #savefig('./ana_sim_figures/ana_sim_summary_%d.png'%model_count)
    #'''

pdf.close()
show()
exit()
#show()
figure()
#loglog(energies_regCS, exp_start_regCS, 'b:', lw=2, label='Geom <5$^\circ$ view angle')
#loglog(energies_lowCS, exp_start_lowCS, 'r:', lw=2)

#loglog(energies_regCS, exp_exit_regCS, 'b-.', lw=2, label=r'$\tau$ exit prob., reg CS')
loglog(energies_lowCS, exp_exit_lowCS, 'r-.', lw=2, label=r'$\tau$ exit prob., low CS')

#loglog(energies_regCS, exp_range_regCS, 'b--', lw=2, label='range cut, reg CS')
loglog(energies_lowCS, exp_range_lowCS, 'r--', lw=2, label='range cut, low CS')

#loglog(energies_regCS, exp_trig_regCS, 'b-', lw=2, label='Triggered, reg CS')
loglog(energies_lowCS, exp_trig_lowCS, 'r-', lw=2, label='Triggered, low CS')

plot(energies_lowCS[:-2],  Auger[:-2],               'k--',  lw=4,  alpha=0.35, label='Auger (2015)')
plot(energies_lowCS[:-2],  np.array(IC_x_3[:-2]), 'k-',   lw=4,  alpha=0.35, label='IceCube (2016)' )
grid(True)
xlim(1.e17, 1.e20)
xlabel('Neutrino Energy, eV')
ylabel('Exposure, cm$^2$ sr s')
legend(loc=0, fontsize=14)
subplots_adjust(bottom=0.125)
show()
