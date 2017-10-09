##
# \file utilities.py
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       November 2016
#


# Import libraries
import os
import re
import numpy as np
import SimpleITK as sitk
from natsort import natsorted


# Import modules
import pythonhelper.SimpleITKHelper as sitkh
import volumetricreconstruction.base.Stack as st

from volumetricreconstruction.definitions import DIR_TMP


##
# Gets the updated affine transforms. At least one input parameter must be a
# list of transforms!
# \date       2016-11-02 23:38:58+0000
#
# \param      transforms_outer  Either as np.array containing list or
#                               sitk.Transform object
# \param      transforms_inner  Either as np.array containing list or
#                               sitk.Transform object
#
# \return     The updated affine transforms with elements transforms_outer[i]
# \f$ \circ\f$ transforms_inner[i]
#
def get_updated_affine_transforms(transforms_outer, transforms_inner):

    # In case transforms_inner is not an array
    try:
        len(transforms_inner)
    except:
        transforms_inner = np.array([transforms_inner])

    # In case transforms_outer is not an array
    try:
        len(transforms_outer)
    except:
        transforms_outer = np.array([transforms_outer])

    # Make len(transforms_inner) = len(transforms_outer)
    if len(transforms_inner) is 1:
        transforms_inner = [transforms_inner[0]] * len(transforms_outer)

    if len(transforms_outer) is 1:
        transforms_outer = [transforms_outer[0]] * len(transforms_inner)

    # Prepare output
    composite_transforms = [None] * len(transforms_inner)

    # Apply transform
    for i in range(0, len(transforms_inner)):
        composite_transforms[i] = sitkh.get_composite_sitk_affine_transform(
            transforms_outer[i], transforms_inner[i])

    return composite_transforms


def get_left_right_mirrored_stack(image_sitk, dir_tmp=DIR_TMP):

    filename_original = os.path.join(dir_tmp, "original.nii.gz")
    filename_swapped = os.path.join(dir_tmp, "original_swapped.nii.gz")
    sitk.WriteImage(image_sitk, filename_original)

    cmd = "fslswapdim "
    cmd += filename_original + " -x y z "
    cmd += filename_swapped
    ph.execute_command(cmd)

    image_sitk = sitk.ReadImage(filename_swapped)
    return image_sitk


def write_results_motion_correction(directory,
                                    filename,
                                    stack0,
                                    stack_corrected,
                                    slice_transforms,
                                    reference_image,
                                    suffix_stack_corrected="_corrected",
                                    suffix_transforms="_slicetransforms_",
                                    ):

    print("Write results after motion correction")
    stack0.write(directory=directory,
                 filename=filename, write_mask=True)
    stack_corrected.write(
        directory=directory,
        filename=filename+suffix_stack_corrected,
        write_mask=True,
        write_slices=True)
    reference_image.write(
        directory=directory,
        filename="Ref_" + reference_image.get_filename(),
        write_mask=False,
        write_slices=False)

    for i in range(0, len(slice_transforms)):
        sitk.WriteTransform(
            slice_transforms[i],
            os.path.join(directory,
                         filename + suffix_transforms + str(i) + ".tfm"))


def read_results_motion_correction(directory,
                                   suffix_stack_corrected="_corrected",
                                   suffix_transforms="_slicetransforms_",
                                   ):

    # if directory[-1] is not "/":
    #     directory += "/"

    p = re.compile("(.*)" + suffix_transforms + "[0-9]+[.]tfm")
    slice_transforms = [p.match(f).group(0)
                        for f in os.listdir(directory) if p.match(f)]
    slice_transforms = natsorted(slice_transforms, key=lambda y: y.lower())
    slice_transforms_sitk = [None]*len(slice_transforms)

    # Required to define resampling grid for reconstruction
    stack_corrected = st.Stack.from_slice_filenames(
        dir_input=directory,
        prefix_stack=p.match(slice_transforms[0]).group(
            1)+suffix_stack_corrected,
        suffix_mask="_mask")

    for i in range(0, len(slice_transforms_sitk)):
        slice_transforms_sitk[i] = sitk.ReadTransform(
            os.path.join(directory, slice_transforms[i]))
        slice_transforms_sitk[i] = sitk.AffineTransform(
            slice_transforms_sitk[i])

    return slice_transforms_sitk, stack_corrected

# def read_results_motion_correction(directory, filename):

#     if directory[-1] is not "/":
#         directory += "/"

#     stack0 = st.Stack.from_filename(
#         os.path.join(directory, filename + ".nii.gz"))
#     stack_corrected = st.Stack.from_slice_filenames(
#         dir_input=directory,
#         prefix_stack=filename+"_corrected",
#         suffix_mask="_mask")

#     slice_transforms_sitk = [None]*stack_corrected.sitk.GetDepth()
#     for i in range(0, len(slice_transforms_sitk)):
#         slice_transforms_sitk[i] = sitk.ReadTransform(
#             os.path.join(directory,
#                          filename + "_slicetransforms_" + str(i) + ".tfm"))
#         slice_transforms_sitk[i] = sitk.AffineTransform(
#             slice_transforms_sitk[i])

#     # Read PD reference used for motion correction
#     p = re.compile("Ref[_](.*)[.]nii.*")
#     filename_reference = [p.match(f).group(1)
#                           for f in os.listdir(directory) if p.match(f)][0]

#     reference_image = st.Stack.from_filename(
#         os.path.join(directory, "Ref_" + filename_reference + ".nii.gz"))
#     reference_image.set_filename(filename_reference)

#     return stack0, stack_corrected, slice_transforms_sitk, reference_image
