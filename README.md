# Volumetric Reconstruction from Printed MR Films 

This is a research-focused toolkit developed within the [GIFT-Surg](http://www.gift-surg.ac.uk/) project to reconstruct a volumetric representation from printed brain MR films as described in [[Ebner2017]](https://www.journals.elsevier.com/neuroimage).

If you have any questions or comments (or find bugs), please drop me an email to `michael.ebner.14@ucl.ac.uk`.

## How it works

1. **Run the Semi-Automatic Slice Extraction Tool:**
A semi-automatic slice extraction tool is used to create a digital image stack from printed slices selected from the scanned brain MR films. It provides an initial digital 3D representation of acquired slices printed on a 2D film where the correct spatial position and dimension of each single slice needs to be recovered.

1. **Recover Meta-Data Information and Correct for In-plane Motion:**
A fully automatic volumetric reconstruction framework to estimate the lost meta-data information of each slice in the 3D space. It is based on a joint slice-to-volume a ne registration with inter-slice 2D transformation regularisation and affine slice-intensity correction. Missing meta-data information is contributed by a longitudinal scan of the same subject.

1. **Recover the 3D Volume using Isotropic Total Variation Denoising:**
A final isotropic total variation in-plane deconvolution technique serves to revitalise the visual appearance of the reconstructed stack of printed slices.


## Installation

This toolkit depends on a variety of software packages developed within the [GIFT-Surg](http://www.gift-surg.ac.uk/) including
* [VolumetricReconstruction](https://cmiclab.cs.ucl.ac.uk/mebner/VolumetricReconstruction)
* [NumericalSolver](https://cmiclab.cs.ucl.ac.uk/mebner/NumericalSolver)
* [RegistrationTools](https://cmiclab.cs.ucl.ac.uk/mebner/RegistrationTools)
* [PythonHelper](https://cmiclab.cs.ucl.ac.uk/mebner/PythonHelper)

All dependencies can be installed using `pip` by running
* `pip install -r requirements.txt`
* `pip install -e .`

Note, that
* the current version relies on the deconvolution operator as implemented in the [VolumetricReconstruction](https://cmiclab.cs.ucl.ac.uk/mebner/VolumetricReconstruction) toolkit. It therefore relies on the ITK version as described there which needs to be installed separately.
* we suggest using the matplotlib version 1.4.3 as more recent versions (we tried matplotlib >= 2.0.0) may slow down the visualization performance of the semi-automatic slice extraction tool substantially.


## Example usage

Run the semi-automatic slice extraction tool to create a a digital image stack from printed slices:

* `python bin/extractStackFromFilms \
--films path-to-film1.dcm ... path-to-filmN.dcm \
--stack path-to-extracted-stack.nii.gz \
--inplane-spacing initial-guess-in-plane-spacing \
--slice-thickness known-slice-thickness
`

A handbook on how to use the semi-automatic slice extraction tool can be found in `doc/`.

---

Recover the meta-data information and correct for in-plane motion of each individual slice by using a reference volume.

* `python bin/correctMotion.py \
--stack path-to-extracted-stack.nii.gz \
--reference path-to-reference-volume.nii.gz \
--dir-output output-directory-for-motion-correction-results \
--dir-output-verbose output-directory-for-intermediate-results
`

Estimated transform parameters for both in-plane similarity and affine transforms are written to the output directory for each single slice. The obtained motion correction results are
used as input for 'reconstructVolume.py' which provides a volumetric
reconstructions in a subsequent step.

---

Based on the estimated transformations a volumetric representation is reconstructed. An additional total variation denoising step is performed for improved visual appearance.

An example call reads:
* `python bin/reconstructVolume.py \
--stack path-to-extracted-stack.nii.gz \
--reference path-to-reference-volume.nii.gz \
--dir-input path-to-motion-correction-results/Similarity \
--dir-output output-directory \
--regularization TV \
--alpha 0.005
`

# References

[[Ebner2017]](https://www.journals.elsevier.com/neuroimage) Ebner, M., Chung, K. K., Prados, F., Cardoso, M. J., Chard, D. T., Vercauteren, T., and Ourselin, S. (In press). Volumetric Reconstruction from Printed Films: Enabling 30 Year Longitudinal Analysis in MR Neuroimaging. NeuroImage.

## License
This framework is licensed under the [MIT license ![MIT](https://raw.githubusercontent.com/legacy-icons/license-icons/master/dist/32x32/mit.png)](http://opensource.org/licenses/MIT)
