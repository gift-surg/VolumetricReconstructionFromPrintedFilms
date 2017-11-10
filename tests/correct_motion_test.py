##
# \file correct_motion_test.py
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


class CorrectMotionTest(unittest.TestCase):

    def test_correct_motion(self):

        path_to_stack = os.path.join(
            DIR_TEST, "SemiAutomaticExtraction_A5208316-B0702861-5yr.nii.gz")
        path_to_reference = os.path.join(
            DIR_TEST, "Reference_A5208316-B0702861-20yr-0-PD.nii.gz")
        dir_output = DIR_TMP

        # Choose parameters for quick processing
        factor_downsampling = 20
        iter_max = 5

        cmd_args = []
        cmd_args.append("--stack %s" % path_to_stack)
        cmd_args.append("--reference %s" % path_to_reference)
        cmd_args.append("--dir-output %s" % dir_output)
        cmd_args.append("--dir-output-verbose %s" % dir_output)
        cmd_args.append("--factor-downsampling %s" % factor_downsampling)
        cmd_args.append("--iter-max %s" % iter_max)

        cmd = "vrpf_correct_motion %s" % (" ").join(cmd_args)
        self.assertEqual(ph.execute_command(cmd), 0)
