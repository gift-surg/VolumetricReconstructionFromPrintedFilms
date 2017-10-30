# Volumetric Reconstruction from Printed Films

This is a Python-based open-source toolkit for research developed within the [GIFT-Surg][giftsurg] project to reconstruct a volumetric representation from printed brain MR films as described in [Ebner et al., 2017][citation].

The algorithm and software were developed by [Michael Ebner][mebner] at the [Translational Imaging Group][tig] in the [Centre for Medical Image Computing][cmic] at [University College London (UCL)][ucl]. Please note that currently **only Python 2** is supported.

If you have any questions or comments, please drop an email to `michael.ebner.14@ucl.ac.uk`.


## How it works

1. **Run the Semi-Automatic Slice Extraction Tool:**
A semi-automatic slice extraction tool is used to create a digital image stack from printed slices selected from the scanned brain MR films. It provides an initial digital 3D representation of acquired slices printed on a 2D film where the correct spatial position and dimension of each single slice needs to be recovered.

1. **Recover Meta-Data Information and Correct for In-plane Motion:**
A fully automatic volumetric reconstruction framework to estimate the lost meta-data information of each slice in the 3D space. It is based on a joint slice-to-volume affine registration with inter-slice 2D transformation regularisation and affine slice-intensity correction. Missing meta-data information is contributed by a longitudinal scan of the same subject.

1. **Recover the 3D Volume using Isotropic Total Variation Denoising:**
A final isotropic total variation in-plane deconvolution technique serves to revitalise the visual appearance of the reconstructed stack of printed slices.


## How to cite
If you use this software in your work, please cite [Ebner et al., 2017][citation].

* Ebner, M., Chung, K. K., Prados, F., Cardoso, M. J., Chard, D. T., Vercauteren, T., and Ourselin, S. (In press). Volumetric Reconstruction from Printed Films: Enabling 30 Year Longitudinal Analysis in MR Neuroimaging. NeuroImage.

## Installation

This toolkit is currently supported for **Python 2 only** and was tested on

* Mac OS X 10.10 and 10.12
* Ubuntu 14.04 and 16.04

It builds on a couple of additional libraries developed within the [GIFT-Surg][giftsurg] project including 
* [NiftyMIC][niftymic]
* [NSoL][nsol]
* [SimpleReg][simplereg]
* [PySiTK][pysitk]
* [ITK_NiftyMIC][itkniftymic]

whose installation requirements need to be met:

1. [Installation of ITK_NiftyMIC][itkniftymic]
1. [Installation of SimpleReg dependencies][simplereg-dependencies]

Afterwards, clone this repository via

* `git clone git@github.com:gift-surg/VolumetricReconstructionFromPrintedFilms.git`

where all remaining dependencies can be installed using `pip`:
* `pip install -r requirements.txt`
* `pip install -e .`

Note, that we suggest using the matplotlib version 1.4.3 as more recent versions (we tried matplotlib >= 2.0.0) may slow down the visualization performance of the semi-automatic slice extraction tool substantially.

Check installation via
* `python -m nose tests/installation_test.py`


## Example usage

Run the semi-automatic slice extraction tool to create a a digital image stack from printed slices:

```
vrpf_extract_stack_from_films \
--films path-to-film1.dcm ... path-to-filmN.dcm \
--stack path-to-extracted-stack.nii.gz \
--inplane-spacing initial-guess-in-plane-spacing \
--slice-thickness known-slice-thickness
```

The current version accepts DICOM (`.dcm`) and NIfTI (`.nii` or `.nii.gz`) as input film types.
A handbook on how to use the semi-automatic slice extraction tool can be found in `doc/`.

---

Recover the meta-data information and correct for in-plane motion of each individual slice by using a reference volume.

```
vrpf_correct_motion \
--stack path-to-extracted-stack.nii.gz \
--reference path-to-reference-volume.nii.gz \
--dir-output output-directory-for-motion-correction-results \
--dir-output-verbose output-directory-for-intermediate-results
```

Estimated transform parameters for both in-plane similarity and affine transforms are written to the output directory for each single slice in separate folders called `Similarity` and `Affine`, respectively. The obtained motion correction results are
used as input for 'vrpf_reconstruct_volume.py' which provides a volumetric
reconstructions in a subsequent step.

---

Based on the estimated transformations a volumetric representation is reconstructed. An additional total variation denoising step is performed for improved visual appearance.

An example call reads:
```
vrpf_reconstruct_volume \
--stack path-to-extracted-stack.nii.gz \
--reference path-to-reference-volume.nii.gz \
--dir-input path-to-motion-correction-results/Similarity \
--dir-output output-directory \
--regularization TV \
--alpha 0.005
```

## Licensing and Copyright
Copyright (c) 2017, [University College London][ucl].
This framework is made available as free open-source software under the [BSD-3-Clause License][bsd]. Other licenses may apply for dependencies.


## Funding
This work is partially funded by the UCL [Engineering and Physical Sciences Research Council (EPSRC)][epsrc] Centre for Doctoral Training in Medical Imaging (EP/L016478/1), the Innovative Engineering for Health award ([Wellcome Trust][wellcometrust] [WT101957] and [EPSRC][epsrc] [NS/A000027/1]), the [Multiple Sclerosis Society of Great Britain and Northern Ireland][mssociety] (grant references 20 and 984) and supported by researchers at the [National Institute for Health Research][nihr] [University College London Hospitals (UCLH)][uclh] Biomedical Research Centre. FP is supported by the [Guarantors of Brain][guarantors].


[citation]: https://www.sciencedirect.com/science/article/pii/S1053811917308042
[mebner]: http://cmictig.cs.ucl.ac.uk/people/phd-students/michael-ebner
[tig]: http://cmictig.cs.ucl.ac.uk
[bsd]: https://opensource.org/licenses/BSD-3-Clause
[giftsurg]: http://www.gift-surg.ac.uk
[cmic]: http://cmic.cs.ucl.ac.uk
[guarantors]: https://guarantorsofbrain.org/
[ucl]: http://www.ucl.ac.uk
[uclh]: http://www.uclh.nhs.uk
[epsrc]: http://www.epsrc.ac.uk
[wellcometrust]: http://www.wellcome.ac.uk
[mssociety]: https://www.mssociety.org.uk/
[nihr]: http://www.nihr.ac.uk/research
[itkniftymic]: https://github.com/gift-surg/ITK_NiftyMIC/wikis/home
[niftymic]: https://github.com/gift-surg/NiftyMIC
[nsol]: https://github.com/gift-surg/NSoL
[simplereg]: https://github.com/gift-surg/SimpleReg
[simplereg-dependencies]: https://github.com/gift-surg/SimpleReg/wikis/simplereg-dependencies
[pysitk]: https://github.com/gift-surg/PySiTK
