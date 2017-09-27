#!/usr/bin/python

# \file
#
#  \author Michael Ebner (michael.ebner.14@ucl.ac.uk)
#  \date Aug 2016


import os
import argparse
import numpy as np
import SimpleITK as sitk

import pythonhelper.SimpleITKHelper as sitkh
import pythonhelper.PythonHelper as ph
import volumetricreconstruction.base.Stack as st
import volumetricreconstruction.preprocessing.BrainStripping as bs
import volumetricreconstruction.preprocessing.IntensityCorrection as ic
import volumetricreconstruction.registration.RegistrationSimpleITK as regsitk
import volumetricreconstruction.registration.RegistrationCppITK as regitk
import volumetricreconstruction.registration.NiftyReg as regniftyreg
import volumetricreconstruction.registration.IntraStackRegistration as intrareg

import volumetricreconstructionfromprintedmrfilms.InputArgparser as inargs
import volumetricreconstructionfromprintedmrfilms.utilities as utils
import volumetricreconstructionfromprintedmrfilms.ScanExtractor as se


if __name__ == '__main__':

    time_start = ph.start_timing()

    input_parser = inargs.InputArgparser(
        description="Run the semi-automatic slice extraction tool to create a "
        "digital image stack from historical slices selected from the scanned "
        "brain MR films. It provides an initial digital 3D representation of "
        "acquired slices printed on a 2D film where the correct spatial "
        "position and dimension of each single slice needs to be recovered "
        "in subsequent steps by using "
        "'correctMotion.py' and 'reconstructVolume.py', respectively.",
        prog="python " + os.path.basename(__file__),
    )
    input_parser.add_films(required=True)
    input_parser.add_stack(required=True)
    input_parser.add_verbose(default=True)
    input_parser.add_inplane_spacing(default=0.14)
    input_parser.add_slice_thickness(default=5.)

    args = input_parser.parse_args()
    input_parser.print_arguments(args)

    if (".").join(os.path.basename(args.stack).split(".")[1:]) \
            not in ["nii", "nii.gz"]:
        raise IOError(
            "Output image (--stack) must be of type 'nii' or 'nii.gz'")

    if args.verbose:
        dir_output_verbose = os.path.dirname(args.stack)
    else:
        dir_output_verbose = None

    # Extract the stack semi-automatically from the specified MRI films
    scan_extractor = se.ScanExtractor(
        filenames=args.films,
        dir_output_verbose=dir_output_verbose)
    scan_extractor.run_semiautomatic_image_extraction()
    image_sitk = scan_extractor.get_sitk_stack_of_extracted_scans()

    # Update meta-data
    spacing = np.array([args.inplane_spacing,
                        args.inplane_spacing,
                        args.slice_thickness])
    image_sitk.SetSpacing(spacing)
    ph.print_info("Stack spacing updated")

    sitk.WriteImage(image_sitk, args.stack)
    ph.print_info("Extracted image stack written to '%s'" % (args.stack))

    elapsed_time = ph.stop_timing(time_start)
    ph.print_title("Summary Semi-automatic Slice Extraction")
    ph.print_info("Elapsed time: %s" % elapsed_time)
