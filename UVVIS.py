import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import seabreeze.spectrometers as sb
import time
import glob

#todo:
#running averages
#previous measurements
#irradiance
#setaverages
#setintegrationtime

class Measurement:
  def __init__(self, integrationtime = 4000, averages = 1):
    self.integrationtime = integrationtime
    self.averages = averages
    self.darkset = False
    self.lightset = False
    self.spec = sb.Spectrometer.from_serial_number()
    self.spec.integration_time_micros(self.integrationtime)
    self.wavelengths = self.spec.wavelengths()
    self.currentspec = self.spec.intensities()
    self.dark = self.spec.intensities()
    self.light = self.spec.intensities()
    self.mode = ''

  def setaverages(self,n):
    self.averages = n
    print('Averages set to: ' + str(n))

  def setintegrationtime(self):
    print('change integration time')

  def setlight(self):
    self.fig, (self.ax0) = plt.subplots(1,1)
    self.line0, = self.ax0.plot([], [], lw=1)
    self.line = [self.line0]
    self.ax0.set_xlim(min(self.wavelengths),max(self.wavelengths))
    self.ax0.set_ylim(-1000, 67000)
    self.ax0.set_ylabel('Intensity (count)')
    self.ax0.set_xlabel('Wavelength (nm)')
    self.ax0.set_title('LIGHT')
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.savelight)
    plt.show()

  def savelight(self,i):
    self.light = self.currentspec
    plt.close('all')
    self.lightset = True
    print('light set')

  def setdark(self):
    self.fig, (self.ax0) = plt.subplots(1,1)
    self.line0, = self.ax0.plot([], [], lw=1)
    self.line = [self.line0]
    self.ax0.set_xlim(min(self.wavelengths),max(self.wavelengths))
    self.ax0.set_ylim(-1000, 67000)
    self.ax0.set_ylabel('Intensity (count)')
    self.ax0.set_xlabel('Wavelength (nm)')
    self.ax0.set_title('DARK')
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.savedark)
    plt.show()

  def savedark(self,i):
    self.dark = self.currentspec
    plt.close('all')
    self.darkset = True
    print('dark set')

  def irradiance(self):
    print('irradiance')

  def absorbance(self,name):
    previous = glob.glob(self.name+'*')
    self.fig, (self.ax0, self.ax1) = plt.subplots(2,1,sharex=True)
    for measurement in previous:
      m = np.genfromtxt(measurement)
      self.ax1.plot(m[:,0], m[:,1], lw=1, alpha=0.3, color='black')
    self.line0, = self.ax0.plot([], [], lw=1)
    self.line1, = self.ax1.plot([], [], lw=1)
    self.line = [self.line0, self.line1]
    self.ax0.set_xlim(min(self.wavelengths),max(self.wavelengths))
    self.ax0.set_ylim(-1000, 67000)
    self.ax0.set_ylabel('Intensity (count)')
    self.ax0.set_title(name)
    self.ax1.set_ylim(-0.1, 3.5)
    self.ax1.set_ylabel('Absorbance')
    self.ax1.set_xlabel('Wavelength (nm)')
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.saveabsorbance)
    self.mode = 'a'
    self.name = name
    plt.show()

  def saveabsorbance(self,i):
    t = int(time.time())
    counter = len(glob.glob(self.name+'*'))
    header = 'Integration Time:\t'+str(self.integrationtime)+'\n'
    header = header+'Averages:\t'+str(self.averages)+'\n'
    header = header+'Unix-Time:\t'+str(t)+'\n\n'
    header = header+'Wavelength (nm)\tAbsorbance\tScope (counts)\tLight (counts)\tDark (counts)'
    np.savetxt(self.name+'_'+str(counter)+'.txt', np.transpose([self.wavelengths,np.log10((self.light-self.dark)/(self.currentspec-self.dark)),self.currentspec,self.light,self.dark]),header=header)
    plt.close('all')
    self.mode = ''

  def animate(self,i):
    spectrum = np.zeros(len(self.wavelengths))
    for k in range(self.averages):
      spectrum = spectrum + self.spec.intensities()
    self.currentspec = spectrum/self.averages
    if self.darkset:
      self.line[0].set_data(self.wavelengths, self.currentspec-self.dark)
    else:
      self.line[0].set_data(self.wavelengths, self.currentspec)
    if self.mode == 'a':
      self.line[1].set_data(self.wavelengths, np.log10((self.light-self.dark)/(self.currentspec-self.dark)))
    return self.line,

m = Measurement(4000,30)
m.setlight()
m.setdark()
windows = ['1_30nm','2_15nm','3_5nm']
for window in windows:
  m.absorbance(window)
