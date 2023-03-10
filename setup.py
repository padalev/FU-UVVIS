from setuptools import setup, find_packages

setup(name='FU-UVVIS',
      version='0.1',
      description='Spectroscopic utilities for Ocean Insight UVVIS spectrometers',
      author='Ruben Nitsche',
      author_email='mail@rubennitsche.eu',
      url='https://github.com/padalev/FU-UVVIS',
      packages=find_packages('UVVIS'),
     )