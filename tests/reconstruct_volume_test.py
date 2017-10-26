##
# \file reconstruct_volume_test.py
#  \brief  Class to test semi-automatic brain extraction tool
#
#  \author Michael Ebner (michael.ebner.14@ucl.ac.uk)
#  \date October 2017


import os
import unittest
import SimpleITK as sitk

import pysitk.python_helper as ph

from volumetricreconstructionfromprintedfilms.definitions import \
    DIR_TMP, DIR_TEST


class ReconstructVolumeTest(unittest.TestCase):

    def test_reconstruct_volume(self):

        path_to_stack = os.path.join(
            DIR_TEST, "SemiAutomaticExtraction_A5208316-B0702861-5yr.nii.gz")
        path_to_reference = os.path.join(
            DIR_TEST, "Reference_A5208316-B0702861-20yr-0-PD.nii.gz")
        dir_input = os.path.join(
            DIR_TEST, "A5208316-B0702861-5yr_motion_correction", "Affine")
        dir_output = DIR_TMP

        # Choose course resolution for quick processing
        resolution_processing = 4
        resolution_reconstruction = 4

        cmd_args = []
        cmd_args.append("--stack %s" % path_to_stack)
        cmd_args.append("--reference %s" % path_to_reference)
        cmd_args.append("--dir-input %s" % dir_input)
        cmd_args.append("--dir-output %s" % dir_output)
        cmd_args.append("--resolution-processing %s" %
                        resolution_processing)
        cmd_args.append("--resolution-reconstruction %s" %
                        resolution_reconstruction)

        cmd = "vrpf_reconstruct_volume %s" % (" ").join(cmd_args)
        ph.execute_command(cmd)
