import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import seabreeze.spectrometers as sb
import time
import glob

#todo:
#irradiance
#setintegrationtime

class Measurement:
  def __init__(self, integrationtime = 4000, averages = 1):
    # state variables to know in case integration time or averages changes
    self.darkset = False
    self.lightset = False
    self.spectra = []
    # start up spectrometer
    try:
      self.spec = sb.Spectrometer.from_serial_number()
      self.spec.integration_time_micros(self.integrationtime)
      self.wavelengths = self.spec.wavelengths()
      self.currentspec = self.spec.intensities()
      self.dark = self.currentspec
      self.light = self.currentspec
    except:
      print('Error starting spectrometer')
      # generate some fake data anyways
      self.wavelengths = np.linspace(200,1000,3648)
      self.currentspec = np.random.normal(5000,50,3648)
      self.dark = self.currentspec
      self.light = self.currentspec
    try:
      self.calibration = np.genfromtxt('flame_calibration.txt')
    except:
      print('Error: could not load calibration data')
    try:
      self.am0 = np.genfromtxt('AM0.tsv')
    except:
      print('Error: could not load AM0 data')
    self.setintegrationtime(integrationtime)
    self.setaverages(averages)
    self.mode = ''

  def setaverages(self,n):
    self.averages = n
    # fill spectra safe
    self.spectra = []
    for i in range(self.averages):
      self.spectra.append(self.light)
    self.darkset = False
    self.lightset = False
    print('Averages set to: ' + str(n))

  def setintegrationtime(self,t=None,auto=False):
    if t:
      self.integrationtime = t
      print('Integration time set to: ' + str(t))
    elif auto:
      self.integrationtime = 4000
      for i in range(5):
        print('Current Integration Time: '+str(self.integrationtime))
        self.spec.integration_time_micros(self.integrationtime)
        time.sleep(self.integrationtime*1e-6 + 1)
        self.currentspec = self.spec.intensities()
        print('Maximum: '+str(np.max(self.currentspec)))
        self.integrationtime = int(self.integrationtime*55000/np.max(self.currentspec))
        self.integrationtime = np.max([self.integrationtime, 4000])
    else:
      #todo
      print('do manual int time')
    self.darkset = False
    self.lightset = False

  def scope(self):
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
    self.ax0.set_title('SCOPE')
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    plt.show()

  def setlight(self):
    self.mode = 'l'
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
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.save)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.save)
    plt.show()

  def setdark(self):
    self.mode = 'd'
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
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.save)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.save)
    plt.show()

  def irradiance(self,name):
    # check if dark is set, else warn and break execution
    if not self.darkset:
      print('Dark Spectrum is not set')
      return
    self.mode = 'i'
    # set name for correct filename
    self.name = name
    # generate new figure with two axis for scope and absorbance
    self.fig, (self.ax0, self.ax1) = plt.subplots(2,1,sharex=True)
    #plot am0
    self.ax1.plot(self.am0[:,0]*1000, self.am0[:,1]/1000, lw=1, alpha=0.6, color='orange')
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
    self.ax1.set_ylim(-0.3, 2.3)
    self.ax1.set_ylabel('Irradiance (W/m$^2$/nm)')
    self.ax1.set_xlabel('Wavelength (nm)')
    # now start the animation and connect action events
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.save)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.save)
    # set absorbance mode so that animation function knows what to plot
    plt.show()

  def absorbance(self,name):
    # check if dark and light are set, else warn and break execution
    if not self.darkset:
      print('Dark Spectrum is not set')
      return
    if not self.lightset:
      print('Light Spectrum is not set')
      return
    self.mode = 'a'
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
    self.ax1.set_ylim(-0.1, 1.5)
    self.ax1.set_ylabel('Absorbance')
    self.ax1.set_xlabel('Wavelength (nm)')
    # now start the animation and connect action events
    self.anim = ani.FuncAnimation(self.fig, self.animate)
    self.cid = self.fig.canvas.mpl_connect('button_press_event', self.save)
    self.cid = self.fig.canvas.mpl_connect('key_press_event', self.save)
    # set absorbance mode so that animation function knows what to plot
    plt.show()

  def save(self,i):
    if self.mode == 'd':
      # set dark
      self.dark = self.currentspec
      self.darkset = True
      print('dark set')
    elif self.mode == 'l':
      # set light
      self.light = self.currentspec
      self.lightset = True
      print('light set')
    elif self.mode == 'a':
      # get time of measurement
      t = int(time.time())
      # get measurement number using previous files
      counter = len(glob.glob(self.name+'*'))
      # save important information as header
      header = 'Integration Time:\t'+str(self.integrationtime)+'\n'
      header = header+'Averages:\t'+str(self.averages)+'\n'
      header = header+'Unix-Time:\t'+str(t)+'\n\n'
      header = header+'Wavelength (nm)\tAbsorbance\tScope (counts)\tLight (counts)\tDark (counts)'
      # save all kinds of spectra to file
      np.savetxt(self.name+'_'+str(counter)+'.txt', np.transpose([self.wavelengths,np.log10((self.light-self.dark)/(self.currentspec-self.dark)),self.currentspec,self.light,self.dark]),header=header)
    elif self.mode == 'i':
      # get time of measurement
      t = int(time.time())
      # get measurement number using previous files
      counter = len(glob.glob(self.name+'*'))
      # save important information as header
      header = 'Integration Time:\t'+str(self.integrationtime)+'\n'
      header = header+'Averages:\t'+str(self.averages)+'\n'
      header = header+'Unix-Time:\t'+str(t)+'\n\n'
      header = header+'Wavelength (nm)\tIrradiance\tScope (counts)\tLight (counts)\tDark (counts)'
      # save all kinds of spectra to file
      np.savetxt(self.name+'_'+str(counter)+'.txt', np.transpose([self.wavelengths,self.calibration[:,1]*(self.currentspec-self.dark)*1000000/self.integrationtime,self.currentspec,self.light,self.dark]),header=header)
    #close plot, reset mode
    plt.close('all')
    self.mode = ''

  def animate(self,i):
    # new empty spectrum
    #spectrum = np.zeros(len(self.wavelengths))
    # collect spectra
    #for k in range(self.averages):
    #  spectrum = spectrum + self.spec.intensities()
    # calculate average
    self.spectra[:-1] = self.spectra[1:]
    try:
      self.spectra[-1] = self.spec.intensities()
    except:
      if self.mode == 'd':
        self.spectra[-1] = np.random.normal(5000,500,3648)
      elif self.mode == 'l':
        self.spectra[-1] = np.random.normal(40000,500,3648)
      elif self.mode == 'a':
        self.spectra[-1] = np.random.normal(20000,500,3648)
      elif self.mode == 'i':
        self.spectra[-1] = np.random.normal(7000,500,3648)
    self.currentspec = np.mean(np.array(self.spectra), axis=0)
    #self.currentspec = spectrum/self.averages
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
      self.line[1].set_data(self.wavelengths, self.calibration[:,1]*(self.currentspec-self.dark)*1000000/self.integrationtime)
    return self.line,

    def calibrationFilePath(self,filepath):
      try:
        self.calibration = np.genfromtxt(filepath)
      except:
        print('Error: could not load calibration data')

    def am0FilePath(self,filepath):
      try:
        self.am0 = np.genfromtxt(filepath)
      except:
        print('Error: could not load AM0 data')


#m = Measurement(4000,30)
#m.setlight()
#m.setdark()
#windows = ['1_30nm','2_15nm','3_5nm']
#for window in windows:
#  m.absorbance(window)