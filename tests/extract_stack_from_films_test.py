##
# \file extract_stack_from_films_test.py
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


class ExtractStackFromFilmsTest(unittest.TestCase):

    def test_extract_stack_from_films(self):

        films = [
            "Film1_A5208316-B0702861-5yr.dcm",
            "Film2_A5208316-B0702861-5yr.nii.gz",
        ]
        path_to_films = [os.path.join(DIR_TEST, f) for f in films]
        path_to_stack = os.path.join(DIR_TMP, "extracted_stack.nii.gz")
        cmd_args = []
        cmd_args.append("--films %s" % (" ").join(path_to_films))
        cmd_args.append("--stack %s" % path_to_stack)
        cmd_args.append("--inplane-spacing 0.14")
        cmd_args.append("--slice-thickness 5")
        cmd = "vrpf_extract_stack_from_films %s" % (" ").join(cmd_args)
        self.assertEqual(ph.execute_command(cmd), 0)
