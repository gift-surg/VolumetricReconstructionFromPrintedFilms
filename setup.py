###
# \file setup.py
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       September 2017
#


from setuptools import setup

long_description = "A flexible framework to reconstruct the volumetric "
"representation from printed films. Is is based on semi-automatic slice "
"extraction, followed by automated slice-to-volume registration with "
"inter-slice transformation regularisation and slice intensity correction. "

setup(name='VolumetricReconstructionFromPrintedFilms',
      version='0.2rc1',
      description='Volumetric Reconstruction from Printed Films: '
      'Enabling 30 Year Longitudinal Analysis in MR Neuroimaging',
      long_description=long_description,
      url='https://github.com/gift-surg/VolumetricReconstructionFromPrintedFilms',
      author='Michael Ebner',
      author_email='michael.ebner.14@ucl.ac.uk',
      license='BSD-3-Clause',
      packages=['volumetricreconstructionfromprintedfilms'],
      install_requires=[
          'pysitk==0.1',
          'nsol==0.1',
          'simplereg==0.1',
          'niftymic==0.1',
          'matplotlib==1.4.3',
          'natsort==5.0.3',
          'numpy==1.13.1',
          'SimpleITK==1.0.1',
      ],
      zip_safe=False,
      keywords="mrireconstruction historicalmrfilmdata brainmri "
      "regularizedimageregistration totalvariationreconstruction "
      "longitudinal analysis",
      classifiers=[
          'Intended Audience :: Developers',
          'Intended Audience :: Healthcare Industry',
          'Intended Audience :: Science/Research',

          'License :: OSI Approved :: BSD License',

          'Topic :: Software Development :: Build Tools',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',

          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
      ],
      entry_points={
          'console_scripts': [
              'vrpf_extract_stack_from_films = volumetricreconstructionfromprintedfilms.application.extract_stack_from_films:main',
              'vrpf_correct_motion = volumetricreconstructionfromprintedfilms.application.correct_motion:main',
              'vrpf_reconstruct_volume = volumetricreconstructionfromprintedfilms.application.reconstruct_volume:main',
          ],
      },
      )
