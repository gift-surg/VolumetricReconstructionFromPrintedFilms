##
# \file InputArparser.py
# \brief      Class holding a collection of possible arguments to parse for
#             scripts
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       September 2017
#

import argparse

import pysitk.python_helper as ph

ALLOWED_FILE_EXTENSIONS = ["nii", "nii.gz"]

# Allowed input file types
FILE_EXTENSIONS = "(" + (", ").join(ALLOWED_FILE_EXTENSIONS) + ")"


##
# Class holding a collection of possible arguments to parse for scripts
# \date       2017-08-07 01:26:11+0100
#
class InputArgparser(object):

    def __init__(self,
                 description=None,
                 prog=None,
                 epilog="Author: Michael Ebner (michael.ebner.14@ucl.ac.uk)",
                 ):

        kwargs = {}
        if description is not None:
            kwargs['description'] = description
        if prog is not None:
            kwargs['prog'] = prog
        if epilog is not None:
            kwargs['epilog'] = epilog

        self._parser = argparse.ArgumentParser(**kwargs)

    def get_parser(self):
        return self._parser

    def parse_args(self):
        return self._parser.parse_args()

    def print_arguments(self, args, title="Input Parameters:"):
        ph.print_title(title)
        for arg in sorted(vars(args)):
            ph.print_info("%s: " % (arg), newline=False)
            print(getattr(args, arg))

    def add_films(
        self,
        option_string="--films",
        type=str,
        nargs="+",
        help="Path to films %s." % ("(.dcm)"),
        required=False,
    ):
        self._add_argument(dict(locals()))

    def add_reference(
        self,
        option_string="--reference",
        type=str,
        help="Path to reference image %s." % (FILE_EXTENSIONS),
        required=False,
    ):
        self._add_argument(dict(locals()))

    def add_stack(
        self,
        option_string="--stack",
        type=str,
        help="Path to semi-automatically extracted image %s." % (
            FILE_EXTENSIONS),
        required=False,
    ):
        self._add_argument(dict(locals()))

    def add_dir_output(
        self,
        option_string="--dir-output",
        type=str,
        help="Output directory.",
        required=False,
        default=None,
    ):
        self._add_argument(dict(locals()))

    def add_dir_input(
        self,
        option_string="--dir-input",
        type=str,
        help="Input directory holding the to obtained motion correction "
        "results obtained by 'correct_motion.py', "
        "e.g. 'dir-motion-correction/Similarity' or "
        "'dir-motion-correction/Affine'.",
        required=False,
        default=None,
    ):
        self._add_argument(dict(locals()))

    def add_dir_output_verbose(
        self,
        option_string="--dir-output-verbose",
        type=str,
        help="Output directory to store all intermediate results.",
        required=False,
        default=None,
    ):
        self._add_argument(dict(locals()))

    def add_iter_max(
        self,
        option_string="--iter-max",
        type=int,
        help="Number of maximum iterations for the numerical solver.",
        default=10,
    ):
        self._add_argument(dict(locals()))

    def add_inplane_spacing(
        self,
        option_string="--inplane-spacing",
        type=float,
        help="Set in-plane spacing for semi-automatically extracted stack.",
        default=1.,
    ):
        self._add_argument(dict(locals()))

    def add_slice_thickness(
        self,
        option_string="--slice-thickness",
        type=float,
        help="Set slice-thickness for semi-automatically extracted stack.",
        default=1.,
    ):
        self._add_argument(dict(locals()))

    def add_factor_inplane_spacing(
        self,
        option_string="--factor-inplane-spacing",
        type=float,
        help="Factor to multiply in-plane spacing of the semi-automatically "
        "extracted stack with",
        default=1.,
    ):
        self._add_argument(dict(locals()))

    def add_factor_downsampling(
        self,
        option_string="--factor-downsampling",
        type=int,
        help="Factor which used to downsample the in-plane resolution of the "
        "semi-automatically extracted stack. This will speed-up the "
        "computations but might affect the accuracy.",
        default=10,
    ):
        self._add_argument(dict(locals()))

    def add_regularization(
        self,
        option_string="--regularization",
        type=str,
        help="Type of regularization g(x) to solve the denoising problem "
        "min_x [||Ax-y||_2^2 + alpha g(x). "
        "Possible choices are 'TK0', 'TK1' or 'TV' for zeroth/first order "
        "Tikhonov  or total variation regularization, respectively."
        "I.e. "
        "g(x) = ||x||_2^2 for 'TK0', "
        "g(x) = ||Dx||_2^2 for 'TK1', "
        "or "
        "g(x) = ||Dx||_1 for 'TV'. ",
        required=False,
        default="TV",
    ):
        self._add_argument(dict(locals()))

    def add_resolution_processing(
        self,
        option_string="--resolution-processing",
        help="In-plane resolution used for processing the images.",
        type=float,
        default=0.25,
    ):
        self._add_argument(dict(locals()))

    def add_resolution_reconstruction(
        self,
        option_string="--resolution-reconstruction",
        help="In-plane resolution for the final volumetric reconstruction.",
        type=float,
        default=1,
    ):
        self._add_argument(dict(locals()))

    def add_alpha(
        self,
        option_string="--alpha",
        type=float,
        help="Regularization parameter alpha to solve minimization problem "
        "min_x [||Ax-y||_2^2 + alpha g(x)] with g denoting the "
        "regularization term.",
        default=0.03,
    ):
        self._add_argument(dict(locals()))

    def add_sigma(
        self,
        option_string="--sigma",
        type=float,
        help="Standard deviation to define blurring operator A ",
        default=0.25,
    ):
        self._add_argument(dict(locals()))

    def add_rho(
        self,
        option_string="--rho",
        type=float,
        help="Regularization parameter for augmented Lagrangian term required "
        "by ADMM approach for TV regularization",
        default=0.5,
    ):
        self._add_argument(dict(locals()))

    def add_iterations(
        self,
        option_string="--iterations",
        type=int,
        help="Number of ADMM iterations.",
        default=10,
    ):
        self._add_argument(dict(locals()))

    def add_verbose(
        self,
        option_string="--verbose",
        type=int,
        help="Turn on/off verbose output.",
        default=1,
    ):
        self._add_argument(dict(locals()))

    ##
    # Adds an argument to argument parser.
    #
    # Rationale: Make interface as generic as possible so that function call
    # works regardless the name of the desired option
    # \date       2017-08-06 21:54:51+0100
    #
    # \param      self     The object
    # \param      allvars  all variables set at respective function call as
    #                      dictionary
    #
    def _add_argument(self, allvars):

        # Skip variable 'self'
        allvars.pop('self')

        # Get name of argument to add
        option_string = allvars.pop('option_string')

        # Build dictionary for additional, optional parameters
        kwargs = {}
        for key, value in allvars.iteritems():
            kwargs[key] = value

        # Add information on default value in case provided
        if 'default' in kwargs.keys():

            txt_default = " [default: %s]" % (str(kwargs['default']))

            # Case where 'required' key is given:
            if 'required' in kwargs.keys():

                # Only add information in case argument is not mandatory to
                # parse
                if kwargs['default'] is not None and not kwargs['required']:
                    kwargs['help'] += txt_default

            # Case where no such field was provided
            else:
                if kwargs['default'] is not None:
                    kwargs['help'] += txt_default

        # Add argument with its options
        self._parser.add_argument(option_string, **kwargs)
