import numpy as np
import scipy.interpolate
import os
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.gridspec as gridspec
import ctypes
from fatiando.gravmag import transform
import os.path

class DBMpy(object):
    def __init__(self):
        print("You're now using DBMpy!\n\n")
    def cwind_txt(self, file_name, wsize, wdx):
        x,y,m = np.loadtxt(file_name, usecols=(0,1,2), unpack=True)
        w_size = wsize
        w_dx = wdx
        w = []
        wins = 0        
        for i in range(int(max(x)/w_dx)):
            x1 = min(x)+(i*w_dx)
            x2 = x1+w_size
            for j in range(int(max(y)/w_dx)):
                y1 = min(y)+(j*w_dx)
                y2 = y1+w_size
                if y1 <= max(y):
                    for k in range(len(m)):
                        if x1 <= x[k] <= x2 and y1 <= y[k] <= y2:
                            w.append([x[k],y[k],m[k]])
                if len(w) > 0:
                    if w[-1][0]-w[0][0] == w_size and w[-1][1]-w[0][1] == w_size:
                        wins+=1
                        arq_w = open('window_%d.txt'%wins, 'a')
                        for k in range(len(w)):
                            arq_w.write('%f, %f, %f\n'%(w[k][0],w[k][1],w[k][2]))
                        print("File 'window_%d.txt' was created!\n"%wins)
                        arq_w.close()
                    del w[:]
                    
    def check_winds(self):
        i = 1
        while True:
            try:
                nome = "window_%d.txt"%i
                x,y,z = np.loadtxt(nome, delimiter=',', usecols=(0,1,2),unpack=True)
                N = 100
                x_grid = np.linspace(x.min(), x.max(), N)
                y_grid = np.linspace(y.min(), y.max(), N)
                xi,yi = np.meshgrid(x_grid,y_grid)
                zi = scipy.interpolate.griddata((x, y), z,
                                        (xi,yi), method = 'cubic')
                i += 1
                if(any(np.isnan(zi)[0])):
                    os.rename(nome, 'BAD_'+nome)
                    print(nome + ' is not a good window for interpretation. File renamed!\n')
            except:
                print('Scan completed!')
                break
        
    def centroid_analysis(self, wi, wf):
        num_janela = wi
        ultima_janela = wf
        while num_janela <= ultima_janela:
            if os.path.isfile("window_%d.txt"%num_janela):
        
                nome = "window_%d.txt"%num_janela
                x,y,mag = np.loadtxt(nome, usecols=(0,1,2), unpack = True, delimiter=',')
                N = int(np.sqrt(len(mag)))
                shape = (N,N)
                x_grid = np.linspace(x.min(), x.max(), N)
                y_grid = np.linspace(y.min(), y.max(), N)
                xi,yi = np.meshgrid(x_grid,y_grid)
                zi = scipy.interpolate.griddata((x, y), mag, 
                                            (xi,yi), method='cubic')
                kx, ky, pds = transform.power_density_spectra(xi, yi, zi, shape)
                k_radial, pds_radial = transform.radial_average_spectrum(kx, ky, pds)
                dk = k_radial[1]-k_radial[0]
                gs = gridspec.GridSpec(2, 2)
                fig = plt.figure()
                ax1 = plt.subplot(gs[:,0])
                cm = ax1.imshow(zi, cmap = 'CMRmap', origin='lower',
                                interpolation = 'spline36', extent=[x.min(),x.max(),
                                                                    y.min(),y.max()])
                plt.colorbar(cm,orientation="horizontal",
                             pad = 0.1, label = "nT", fraction=0.1)
                ax1.set_xlabel("WE %d m"%(max(x)-min(x)))
                ax1.set_ylabel("SN %d m"%(max(y)-min(y)))
                plt.setp(ax1.get_xticklabels(), visible=False)
                plt.setp(ax1.get_yticklabels(), visible=False)
                ax1.tick_params(axis='both', which='both', length=0)
                plt.title("Anomalia magnetica - " + nome)
                ax2 = plt.subplot(gs[0, 1])
                ax2.set_title('Media radial do espectro de potencia de densidade')
                ptsZt = []
                for i,j in zip(k_radial,pds_radial):
                    ptZt = ax2.scatter(i/2*np.pi/10, np.log(np.sqrt(j))/(2*np.pi),
                          c = 'red', s = 7, picker = 1.5)
                    ptsZt.append(ptZt)
                plt.autoscale(False)
                ax2.set_ylabel('ln(sqr(P))/2pi')
                ax2.ticklabel_format(axis='x', style='sci', scilimits=(1, 1))
                plt.grid(ls = '--', lw = .3)
                ax2.set_xlim(ptsZt[0].get_offsets()[0][0]-dk,
                             ptsZt[-1].get_offsets()[0][0]+dk)
                #ax2.set_xlim(-2*10**-6,2*10**-4)
                ax3 = plt.subplot(gs[1, 1])
                ptsZ0 = []
                for i,j in zip(k_radial,pds_radial):
                    ptZ0 = ax3.scatter(i/2*np.pi/10, np.log(np.sqrt(j)/i)/(2*np.pi),
                         c = 'blue', s = 7, picker = 1.5)
                    ptsZ0.append(ptZ0)
                ax3.set_ylabel('(ln(sqr(P))/|k|)/2pi')
                ax3.set_xlabel('|k|/2pi (km-1)')
                ax3.ticklabel_format(axis='x', style='sci', scilimits=(1, 1))
                plt.grid(ls = '--', lw = .3)
                plt.autoscale(False)
                ax3.set_xlim(ptsZt[0].get_offsets()[0][0]-dk,
                             ptsZt[-1].get_offsets()[0][0]+dk)
               # ax3.set_xlim(-2*10**-6,2*10**-4)
                arts_zt, arts_z0, zt_x, zt_y, z0_x, z0_y = [],[],[],[],[],[]
                pts_calc_zt, pts_calc_z0 = {'x':[],'y':[]},{'x':[],'y':[]}
                zs = {'zt':None, 'z0':None}
                def on_pick(event):
                    for i in ptsZt:
                        if i == event.artist:
                            if len(arts_zt) == 0:
                                zt_x.append(i.get_offsets()[0][0])
                                zt_y.append(i.get_offsets()[0][1])
                                vx = i.get_offsets()[0][0]
                                lv = ax2.plot([vx,vx],[ax2.get_ylim()[0],ax2.get_ylim()[1]],
                                              c = 'red', lw = 1, ls = '--')[0]
                                arts_zt.append(lv)
                                fig.canvas.draw()
        
                            elif len(arts_zt) == 1:
                                vx = i.get_offsets()[0][0]
                                lv = ax2.plot([vx,vx],[ax2.get_ylim()[0],ax2.get_ylim()[1]],
                                              c = 'red', lw = 1, ls = '--')[0]
                                arts_zt.append(lv)
                                zt_x.append(i.get_offsets()[0][0])
                                zt_y.append(i.get_offsets()[0][1])
        
                                for pt in ptsZt:
                                    xpt = pt.get_offsets()[0][0]
                                    ypt = pt.get_offsets()[0][1]
                                    if zt_x[0] <= xpt <= zt_x[1]:
                                        pts_calc_zt['x'].append(xpt)
                                        pts_calc_zt['y'].append(ypt)
                                    
                                m, b = np.polyfit(pts_calc_zt['x'], pts_calc_zt['y'], 1)
                                l = ax2.plot(pts_calc_zt['x'], m*np.array(pts_calc_zt['x'])+b,
                                             c = 'black', lw = 1.5)[0]
                                arts_zt.append(l)
                                t = ax2.text(i.get_offsets()[0][0],i.get_offsets()[0][1],
                                            '  zt = %.2f m'%abs(m), fontsize=12)
                                arts_zt.append(t)
                                zs['zt'] = abs(m)                
                                fig.canvas.draw()
        
                            else:
                                for j in arts_zt:
                                    j.remove()
                                del arts_zt[:]
                                del zt_x[:]
                                del zt_y[:]
                                zs['zt'] = None
                                del pts_calc_zt['x'][:]
                                del pts_calc_zt['y'][:]
                                fig.canvas.draw()
                                
                    for i in ptsZ0:
                        if i == event.artist:
                            if len(arts_z0) == 0:
                                z0_x.append(i.get_offsets()[0][0])
                                z0_y.append(i.get_offsets()[0][1])
                                vx = i.get_offsets()[0][0]
                                lv = ax3.plot([vx,vx],[ax3.get_ylim()[0],ax3.get_ylim()[1]],
                                              c = 'blue', lw = 1, ls = '--')[0]
                                arts_z0.append(lv)
                                fig.canvas.draw()
        
                            elif len(arts_z0) == 1:
                                vx = i.get_offsets()[0][0]
                                lv = ax3.plot([vx,vx],[ax3.get_ylim()[0],ax3.get_ylim()[1]],
                                              c = 'blue', lw = 1, ls = '--')[0]
                                arts_z0.append(lv)
                                z0_x.append(i.get_offsets()[0][0])
                                z0_y.append(i.get_offsets()[0][1])
                                for pt in ptsZ0:
                                    try:
                                        xpt = pt.get_offsets()[0][0]
                                        ypt = pt.get_offsets()[0][1]
                                        if z0_x[0] <= xpt <= z0_x[1]:
                                            pts_calc_z0['x'].append(xpt)
                                            pts_calc_z0['y'].append(ypt)
                                    except:
                                        pass
                                    
                                m, b = np.polyfit(pts_calc_z0['x'], pts_calc_z0['y'], 1)
                                l = ax3.plot(pts_calc_z0['x'], m*np.array(pts_calc_z0['x'])+b,
                                             c = 'black', lw = 1.5)[0]
                                arts_z0.append(l)
                                t = ax3.text(i.get_offsets()[0][0],i.get_offsets()[0][1],
                                            '  z0 = %.2f m'%abs(m), fontsize=12)
                                arts_z0.append(t)
                                zs['z0'] = abs(m)                
                                fig.canvas.draw()
        
                            else:
                                for j in arts_z0:
                                    j.remove()
                                del arts_z0[:]
                                del z0_x[:]
                                del z0_y[:]
                                zs['z0'] = None
                                del pts_calc_z0['x'][:]
                                del pts_calc_z0['y'][:]
                                fig.canvas.draw()
        
                def on_key(event):
                    if zs['zt'] != None and zs['z0'] != None:
                        zb = 2*zs['z0']-zs['zt']
                        xmed = x.min()+(x.max()-x.min())/2
                        ymed = y.min()+(y.max()-y.min())/2
                        if event.key == 'd' or event.key == 'D':
                            print('zt = %.2f m\nz0 = %.2f m\nzb = %.2f m\nx = %.2f (m)\ny = %.2f (m)'%(zs['zt'],
                                zs['z0'],zb,xmed,ymed))
                            ctypes.windll.user32.MessageBoxW(None,
                                u'zt = %.2f m\nz0 = %.2f m\nzb = %.2f m\nx = %.2f (m)\ny = %.2f (m)'%(zs['zt'],
                                zs['z0'],zb,xmed,ymed), u'Curie depth', 0)
                        elif event.key == 's' or event.key == 'S':
                            arq = open('sup_curie.txt', 'a')
                            arq.write("%f, %f, %f, %f, %f, %s\n"%(xmed,
                                ymed,zb,zs['zt'],zs['z0'],nome))
                            arq.close()
                            print('''Foi adicionado ao arquivo sup_curie.txt os seguintes dados:
                 X (m), Y (m), zb (m), zt (m), z0 (m), nome da janela
                 %f, %f, %f, %f, %f, %s'''%(xmed,ymed,zb,zs['zt'],zs['z0'],nome))
        
                def on_press(event):
                    if event.inaxes == ax1:
                        if cm.get_cmap() == plt.get_cmap('jet'):
                            cm.set_cmap('CMRmap')
                        elif cm.get_cmap() == plt.get_cmap('CMRmap'):
                            cm.set_cmap('viridis')
                        elif cm.get_cmap() == plt.get_cmap('viridis'):
                            cm.set_cmap('hot')
                        elif cm.get_cmap() == plt.get_cmap('hot'):
                            cm.set_cmap('gist_gray')
                        elif cm.get_cmap() == plt.get_cmap('gist_gray'):
                            cm.set_cmap('jet')
                        fig.canvas.draw()
        
                mouse = fig.canvas.mpl_connect('button_press_event', on_press)
                key = fig.canvas.mpl_connect('key_press_event', on_key)
                cid = fig.canvas.mpl_connect('pick_event', on_pick)
                
                plt.show()
                num_janela+=1
            else:
                num_janela+=1
        
d = DBMpy()
#d.cwind_txt(file_name = "AM_DOC_1km_PLB_WGS84.XYZ", wsize = 250000, wdx = 100000)
#d.check_winds()
#d.centroid_analysis(wi=3,wf=198)














