import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import seabreeze.spectrometers as sb
import time
import glob

#todo:
#running averages
#irradiance
#setaverages
#setintegrationtime

class Measurement:
  def __init__(self, integrationtime = 4000, averages = 1):
    self.integrationtime = integrationtime
    self.averages = averages
    # state variables to know in case integration time or averages changes
    self.darkset = False
    self.lightset = False
    # start up spectrometer
    try:
      self.spec = sb.Spectrometer.from_serial_number()
      self.spec.integration_time_micros(self.integrationtime)
      self.wavelengths = self.spec.wavelengths()
      self.currentspec = self.spec.intensities()
      self.dark = self.spec.intensities()
      self.light = self.spec.intensities()
    except:
      print('Error starting spectrometer')
    self.mode = ''

  def setaverages(self,n):
    self.averages = n
    print('Averages set to: ' + str(n))

  def setintegrationtime(self,t=None,auto=False):
    if t:
      self.integrationtime = t
      print('Integration time set to: ' + str(t))
    elif auto:
      #todo
      print('do auto int time')
    else:
      #todo
      print('do manual int time')

  def setlight(self):
    # generate a plot to set the light spectrum with just the scope
    self.fig, (self.ax0) = plt.subplots(1,1)
    # generate empty line for live plotting
    self.line0, = self.ax0.plot([], [], lw=1)
    self.line = [self.line0]
    # aesthetics ...
    self.ax0.set_xlim(min(self.wavelengths),max(self.wavelengths))
    self.ax0.set_ylim(-1000, 67000)
    self.ax0.set_ylabel('Intensity (count)')
    self.ax0.set_xlabel('Wavelength (nm)')
    self.ax0.set_title('LIGHT')
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.savelight)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.savelight)
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
    # now start the animation and connect action events
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.savedark)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.savedark)
    plt.show()

  def savedark(self,i):
    # set dark, close plot
    self.dark = self.currentspec
    plt.close('all')
    self.darkset = True
    print('dark set')

  def irradiance(self):
    # todo
    print('irradiance')

  def absorbance(self,name):
    # set name for correct filename
    self.name = name
    # generate new figure with two axis for scope and absorbance
    self.fig, (self.ax0, self.ax1) = plt.subplots(2,1,sharex=True)
    # look for previously saved spectra and plot them into the absorbance graph
    previous = glob.glob(self.name+'*')
    for measurement in previous:
      m = np.genfromtxt(measurement)
      self.ax1.plot(m[:,0], m[:,1], lw=1, alpha=0.3, color='black')
    # plot empty lines for live plotting
    self.line0, = self.ax0.plot([], [], lw=1)
    self.line1, = self.ax1.plot([], [], lw=1)
    self.line = [self.line0, self.line1]
    # some aesthetics ...
    self.ax0.set_xlim(min(self.wavelengths),max(self.wavelengths))
    self.ax0.set_ylim(-1000, 67000)
    self.ax0.set_ylabel('Intensity (count)')
    self.ax0.set_title(name)
    self.ax1.set_ylim(-0.1, 2.5)
    self.ax1.set_ylabel('Absorbance')
    self.ax1.set_xlabel('Wavelength (nm)')
    # now start the animation and connect action events
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.saveabsorbance)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.saveabsorbance)
    # set absorbance mode so that animation function knows what to plot
    self.mode = 'a'
    plt.show()

  def saveabsorbance(self,i):
    # get time of measurement
    t = int(time.time())
    # get measurement number using previous files
    counter = len(glob.glob(self.name+'*'))
    # save important information as header
    header = 'Integration Time:\t'+str(self.integrationtime)+'\n'
    header = header+'Averages:\t'+str(self.averages)+'\n'
    header = header+'Unix-Time:\t'+str(t)+'\n\n'
    header = header+'Wavelength (nm)\tAbsorbance\tScope (counts)\tLight (counts)\tDark (counts)'
    # save all kinds of spectra to file and close plot, reset mode
    np.savetxt(self.name+'_'+str(counter)+'.txt', np.transpose([self.wavelengths,np.log10((self.light-self.dark)/(self.currentspec-self.dark)),self.currentspec,self.light,self.dark]),header=header)
    plt.close('all')
    self.mode = ''

  def animate(self,i):
    # new empty spectrum
    spectrum = np.zeros(len(self.wavelengths))
    # collect spectra
    for k in range(self.averages):
      spectrum = spectrum + self.spec.intensities()
    # calculate average
    self.currentspec = spectrum/self.averages
    if self.darkset:
      # if dark spectrum is set: subtract that
      self.line[0].set_data(self.wavelengths, self.currentspec-self.dark)
    else:
      # else plot just spectrum
      self.line[0].set_data(self.wavelengths, self.currentspec)
    if self.mode == 'a':
      # for absorbance mode plot absorbance
      self.line[1].set_data(self.wavelengths, np.log10((self.light-self.dark)/(self.currentspec-self.dark)))
    if self.mode == 'i':
      # todo: irradiance plot
      pass
    return self.line,

m = Measurement(4000,30)
m.setlight()
m.setdark()
windows = ['1_30nm','2_15nm','3_5nm']
for window in windows:
  m.absorbance(window)
