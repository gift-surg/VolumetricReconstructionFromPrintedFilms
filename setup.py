###
# \file setup.py
#
# Install with symlink: 'pip install -e .'
# Changes to the source file will be immediately available to other users
# of the package
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       September 2017
#
# \see https://python-packaging.readthedocs.io/en/latest/minimal.html
# \see https://python-packaging-user-guide.readthedocs.io/tutorials/distributing-packages/


from setuptools import setup

long_description = ""

setup(name='VolumetricReconstructionFromPrintedFilms',
      version='0.1.dev1',
      description='Volumetric Reconstruction from Printed Films: '
      'Enabling 30 Year Longitudinal Analysis in MR Neuroimaging',
      long_description=long_description,
      url='https://github.com/gift-surg/VolumetricReconstructionFromPrintedFilms',
      author='Michael Ebner',
      author_email='michael.ebner.14@ucl.ac.uk',
      license='BSD-3-Clause',
      packages=['volumetricreconstructionfromprintedfilms'],
      install_requires=[
          'pysitk',
          'nsol',
          'simplereg',
          'niftymic',
      ],
      zip_safe=False,
      keywords="mrireconstruction historicalmrfilmdata brainmri "
      "regularizedimageregistration totalvariationreconstruction "
      "longitudinal analysis",
      classifiers=[
          'Development Status :: 3 - Alpha',

          'Intended Audience :: Developers',
          'Topic :: Software Development :: Build Tools',

          'License :: OSI Approved :: BSD License',

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
