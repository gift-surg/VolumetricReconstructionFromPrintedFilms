#!/usr/bin/python

# \file ScanExtractor.py
#  \brief Extract scans semi-automatically from MR films
#
#  \author Michael Ebner (michael.ebner.14@ucl.ac.uk)
#  \date Nov 2016


# Import libraries
import SimpleITK as sitk
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import pythonhelper.PythonHelper as ph
import volumetricreconstructionfromprintedfilms.utilities.FigureEventHandling as feh

##
#       Extract scans semi-automatically from MR films
# \date       2016-11-02 00:17:44+0000
#


class ScanExtractor(object):

    ##
    # Constructor
    # \date       2016-09-19 11:00:37+0100
    #
    # \param      self                        The object
    # \param      filenames                   Filenames of all MR films of same
    #                                         scan (including extension)
    # \param      dir_output_verbose          Output directory to store results
    # \param      selection_window_offset     The offset for the class
    #                                         FigureEventHandling
    # \param      selection_window_dimension  The length for the class
    #                                         FigureEventHandling
    # \param      bookmark_default_integer    The bookmark default integer
    # \return     Stack of extracted slices consisting of all images on
    #             number_mr_films MR films as sitk.Image object
    #
    def __init__(self,
                 filenames,
                 dir_output_verbose=None,
                 selection_window_offset=np.array([100, -900]),
                 selection_window_dimension=np.array([1350, 1700]),
                 bookmark_default_integer=3,
                 ):

        self._filenames = filenames
        self._dir_output_verbose = dir_output_verbose

        self._selection_window_offset = selection_window_offset
        self._selection_window_dimension = selection_window_dimension
        self._bookmark_default_integer = bookmark_default_integer

        self._stack_sitk = None

    ##
    # Sets the bookmark default integer value to choose predefined
    # dimensions for selection window
    # \date       2017-01-26 16:33:08+0000
    #
    # \param      self                      The object
    # \param      bookmark_default_integer  The bookmark default integer
    #
    def set_bookmark_default_integer(self, bookmark_default_integer):
        self._bookmark_default_integer = bookmark_default_integer

    def get_bookmark_default_integer(self):
        return self._bookmark_default_integer

    ##
    # Get semi-automatically extracted stack of scans
    # \date       2016-11-02 01:23:39+0000
    #
    # \param      self  The object
    #
    # \return     The stack of extracted scans as sitk.Image object
    #
    def get_sitk_stack_of_extracted_scans(self):
        if self._stack_sitk is None:
            raise RuntimeError(
                "Stack has not been extracted yet. "
                "Run 'run_semiautomatic_image_extraction' first.")
        return sitk.Image(self._stack_sitk)

    ##
    #       Run semi-automatic pipeline to extract the scans from the film
    # \date       2016-11-02 00:23:58+0000
    #
    # \param  dir_input                   The dir input
    # \param  timepoint                   The timepoint
    # \param  number_mr_films             The number mr films
    # \param  selection_window_offset     The selection window offset
    # \param  selection_window_dimension  The selection window dimension
    # \param  self._dir_output_verbose                  The dir output
    #
    # \return     { description_of_the_return_value }
    #
    def run_semiautomatic_image_extraction(self):

        stack_nda_list = []
        self._partial_stack_sitk = []

        for i, filename in enumerate(self._filenames):

            ph.print_info("Open MRI film %s (%d/%d) ... " %
                          (filename, i+1, len(self._filenames)))

            # Read image
            image = sitk.ReadImage(filename)

            # Convert to data array
            nda = sitk.GetArrayFromImage(image).squeeze()

            plt.imshow(nda, cmap="Greys_r")

            # Save image as png
            if self._dir_output_verbose is not None:
                ph.create_directory(self._dir_output_verbose)
                filename_without_ext = os.path.basename(
                    filename).split(".")[0]
                filename_out = os.path.join(self._dir_output_verbose,
                                            filename_without_ext + ".png")
                plt.savefig(filename_out, dpi=400)
                ph.print_info("File written to " + filename_out)

            # Instantiate object to mark coordinates and feed with initial
            # offset and length
            figure_event_handling = feh.FigureEventHandling(
                nda, title=filename)

            figure_event_handling.set_offset(self._selection_window_offset)
            figure_event_handling.set_length(self._selection_window_dimension)

            # Mark coordinates to extract images and stack them
            figure_event_handling.extract_slices_semiautomatically()

            # Get coordinates, offset and length of selected windows, i.e.
            # slices
            coordinates = figure_event_handling.get_coordinates()
            offset = figure_event_handling.get_offset()
            length = figure_event_handling.get_length()

            # Update member attributes
            self._selection_window_offset = offset
            self._selection_window_dimension = length

            # Possible that no frame was selected (i.e. no appropriate image on
            # film)
            if len(coordinates) > 0:
                # Get stacked array of slices
                partial_stack_nda = \
                    self._get_stacked_slices_data_array_from_MR_film(
                        nda, coordinates, offset, length)

                # Append partial stack array to list
                stack_nda_list.append(partial_stack_nda)

                # Convert to sitk.Image
                self._partial_stack_sitk.append(
                    sitk.GetImageFromArray(partial_stack_nda))

            # Show image
            # sitk.Show(self._partial_stack_sitk[i])

        if len(stack_nda_list) > 0:
            # Get one entire data array
            stack_nda = \
                self._get_combined_stacked_slices_data_array_from_MR_films(
                    stack_nda_list)

            # Create sitk.Image
            self._stack_sitk = sitk.GetImageFromArray(stack_nda)

    ##
    #       Get the data array of the stacked slices of one MRI film
    # \date       2016-09-19 11:07:14+0100
    #
    # \param  nda          2D data array representing the MRI film
    # \param  coordinates  The coordinates of selected points of all selected
    #                          regions, i.e. single scans
    # \param  offset       The offset describing the north-west corner of each
    #                          selected region in nda
    # \param  length       The length describing the length in x and y of each
    #                          single scan
    #
    # \return     The stacked slices of partial stack as 3D data array
    #
    def _get_stacked_slices_data_array_from_MR_film(self,
                                                    nda,
                                                    coordinates,
                                                    offset,
                                                    length):
        N_slices = len(coordinates)
        stack_nda = np.zeros((N_slices, length[1], length[0]))

        for i in range(0, N_slices):
            stack_nda[i, :, :] = nda[
                coordinates[i][1]+offset[1]:
                coordinates[i][1]+offset[1]+length[1],
                coordinates[i][0]+offset[0]:
                coordinates[i][0]+offset[0]+length[0]
            ]

        return stack_nda

    ##
    #       Gets the combined stack of all slices from all MR films to one
    #             patient
    # \date       2016-09-19 11:15:46+0100
    #
    # \param  stack_nda_list  List of 3D data arrays comprising all the MR
    #                             films
    # \return     3D data array comprising all 2D slices compound to one stack.
    #
    def _get_combined_stacked_slices_data_array_from_MR_films(self,
                                                              stack_nda_list):

        N_slices = 0
        for i in range(0, len(stack_nda_list)):
            N_slices += stack_nda_list[i].shape[0]

        nda_stack = np.zeros((
            N_slices,
            stack_nda_list[-1].shape[1],
            stack_nda_list[-1].shape[2]))

        i_min = 0
        i_max = 0
        for i in range(0, len(stack_nda_list)):
            i_max += stack_nda_list[i].shape[0]
            nda_stack[i_min:i_max, :, :] = stack_nda_list[i][:, :, :]
            i_min += stack_nda_list[i].shape[0]

        return nda_stack
