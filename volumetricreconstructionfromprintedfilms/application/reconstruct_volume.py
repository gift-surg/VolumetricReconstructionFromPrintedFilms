##
# \file reconstruct_volume.py
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       Nov 2016
#


import SimpleITK as sitk
import numpy as np
import os

import niftymic.base.stack as st
import niftymic.reconstruction.admm_solver as admm
import niftymic.reconstruction.tikhonov_solver as tk
import niftymic.utilities.brain_stripping as bs
import niftymic.utilities.intensity_correction as ic
import niftymic.base.psf as psf
import pysitk.python_helper as ph
import pysitk.simple_itk_helper as sitkh
import volumetricreconstructionfromprintedfilms.utilities.input_argparser as inargs
import volumetricreconstructionfromprintedfilms.utilities.utilities as utils


# noinspection PyPep8Naming
def main():
    time_start = ph.start_timing()

    input_parser = inargs.InputArgparser(
        description="Based on the estimated transformations obtained by "
                    "'correct_motion.py' a volumetric representation is "
                    "reconstructed. An additional in-plane denoising "
                    "step is performed for improved visual appearance",
    )
    input_parser.add_stack(required=True)
    input_parser.add_reference(required=True)
    input_parser.add_dir_input(required=True)
    input_parser.add_dir_output(
        required=True,
        help="Output directory to store volumetric "
             "reconstruction results")
    input_parser.add_regularization(default="TV")
    input_parser.add_alpha(
        default=0.003  # TV
        # default=0.03  # TK1
    )
    input_parser.add_option(
        option_string="--intensity-correction",
        type=int,
        default=1,
        help="Turn on/off intensity correction.")
    input_parser.add_option(
        option_string="--additional-frame",
        type=float,
        default=0,
        help="Specify extra surrounding frame (in mm) to be added in the "
        "in-plane dimensions for the reconstructed image. "
        "The slice-thickness remains unchanged.")
    input_parser.add_rho(default=0.5)
    input_parser.add_iterations(default=10)
    input_parser.add_iter_max(default=10)
    input_parser.add_sigma2(default=0.25)
    input_parser.add_resolution_processing(default=0.25)
    input_parser.add_resolution_reconstruction(default=1.)
    input_parser.add_verbose(default=False)

    args = input_parser.parse_args()
    input_parser.print_arguments(args)

    # ---------------------------------------------------------------------
    # Read reference image
    ph.print_title("Read Data")
    ph.print_info("Read reference image")
    reference_image = st.Stack.from_filename(args.reference)

    # ---------------------------------------------------------------------
    # Read motion correction results
    ph.print_info("Read motion correction results")
    slice_transforms_sitk, stack_corrected = \
        utils.read_results_motion_correction(args.dir_input)
    stack0 = st.Stack.from_filename(args.stack)

    # Extract filename without filename extension
    filename_stack = os.path.basename(args.stack).split(".")[0]

    # ---------------------------------------------------------------------
    # Define resampling and reconstruction grids
    ph.print_info("Define resampling and reconstruction grids")
    slice_thickness = stack0.sitk.GetSpacing()[2]
    resampling_grid_sitk = 0 * sitkh.get_downsampled_sitk_image(
        stack_corrected.sitk,
        new_spacing=(args.resolution_processing,
                     args.resolution_processing,
                     slice_thickness)
    )

    # Get enlarged FOV
    resampling_grid_sitk = sitkh.get_altered_field_of_view_sitk_image(
        resampling_grid_sitk,
        boundary_i=args.additional_frame,
        boundary_j=args.additional_frame,
        boundary_k=0,
        unit="mm")

    recon_grid_sitk = 0 * sitkh.get_downsampled_sitk_image(
        resampling_grid_sitk,
        new_spacing=(args.resolution_reconstruction,
                     args.resolution_reconstruction,
                     slice_thickness))
    recon_grid = st.Stack.from_sitk_image(recon_grid_sitk)

    # ---------------------------------------------------------------------
    # Correct for motion
    ph.print_title("Correct for motion")
    stack = st.Stack.from_stack(stack0)
    stack.update_motion_correction_of_slices(slice_transforms_sitk)

    # Default pixel value for resampling
    # Rationale: Due to high background noise, a zero pixel value would not be
    # suitable
    # noinspection PyTypeChecker
    default_pixel_value = np.percentile(
        np.array(sitk.GetArrayFromImage(stack.sitk)), 0.1)

    # verbose:
    if args.verbose:
        sitkh.show_stacks([
            stack.get_resampled_stack_from_slices(
                resampling_grid=recon_grid_sitk,
                interpolator="BSpline",
                default_pixel_value=default_pixel_value),
            reference_image.get_resampled_stack(
                resampling_grid=recon_grid_sitk,
                interpolator="BSpline")
        ])

    # Write results
    i = len(slice_transforms_sitk) / 2  # select midslice trafo for alignment
    stack0.update_motion_correction(slice_transforms_sitk[i])
    stack_naivelyscaled_recon_grid = \
        stack0.get_resampled_stack_from_slices(
            resampling_grid=recon_grid_sitk,
            interpolator="BSpline",
            default_pixel_value=default_pixel_value)
    stack_naivelyscaled_recon_grid.write(
        args.dir_output, filename_stack + "_recon-space")

    stack_motioncorrected_recon_grid = \
        stack.get_resampled_stack_from_slices(
            resampling_grid=recon_grid_sitk,
            interpolator="BSpline",
            default_pixel_value=default_pixel_value)
    stack_motioncorrected_recon_grid.write(
        args.dir_output, filename_stack + "_motion-corrected")

    # ---------------------------------------------------------------------
    # Get brain mask for reference image
    ph.print_title("Get brain mask for reference image")
    brain_stripping = bs.BrainStripping()
    brain_stripping.set_input_image_sitk(reference_image.sitk)
    brain_stripping.run()
    reference_image_sitk_mask = brain_stripping.get_brain_mask_sitk()
    reference_image = st.Stack.from_sitk_image(
        reference_image.sitk,
        reference_image.get_filename(),
        reference_image_sitk_mask)

    # ---------------------------------------------------------------------
    # Resampling to processing and reconstruction grids
    ph.print_title("Resample motion-corrected stack to processing grid")
    # noinspection PyTypeChecker
    default_pixel_value = np.percentile(
        np.array(sitk.GetArrayFromImage(stack.sitk)), 0.1)
    stack_resampled = stack.get_resampled_stack_from_slices(
        resampling_grid=resampling_grid_sitk,
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    reference_image_resampled = reference_image.get_resampled_stack(
        resampling_grid=resampling_grid_sitk, interpolator="BSpline")
    stack_resampled.set_filename(filename_stack + "_motion-corrected")

    ph.print_title("Resample original stack to reconstruction grid")
    stack0_resampled = st.Stack.from_stack(stack0)
    stack0_resampled = stack0_resampled.get_resampled_stack_from_slices(
        resampling_grid=recon_grid_sitk,
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    stack0_resampled.set_filename(filename_stack)

    # ---------------------------------------------------------------------
    if args.intensity_correction:
        # Perform intensity correction
        # sitkh.show_stacks([stack0_resampled, stack_resampled], title=["0","1"])
        ph.print_title("Perform intensity correction")
        intensity_correction = ic.IntensityCorrection(
            stack=stack_resampled,
            reference=reference_image_resampled,
            use_reference_mask=True,
            use_verbose=True)
        intensity_correction.set_additional_stack(stack0_resampled)
        intensity_correction.use_individual_slice_correction(False)
        intensity_correction.run_affine_intensity_correction()
        # intensity_correction.use_individual_slice_correction(False)
        intensity_correction.run_lower_percentile_capping_of_stack(
            percentile=25)
        # intensity_correction.use_individual_slice_correction(True)
        intensity_correction.run_linear_intensity_correction()
        # noinspection PyPep8Naming
        stack_intensityCorrected = intensity_correction.get_intensity_corrected_stack()
        stack0_intensityCorrected = \
            intensity_correction.get_intensity_corrected_additional_stack()

        stack_intensityCorrected.set_filename(
            filename_stack + "_motion-corrected-ic")
        stack0_intensityCorrected.set_filename(filename_stack + "_recon-space-ic")
    else:
        stack_intensityCorrected = stack_resampled
        stack0_intensityCorrected = stack0_resampled

    # Write results
    if args.intensity_correction:
        stack_naivelyscaledic_recon_grid = stack0_intensityCorrected
        stack_naivelyscaledic_recon_grid.write(
            args.dir_output, filename_stack + "_recon-space-ic")

        stack_motioncorrectedic_recon_grid = \
            stack_intensityCorrected.get_resampled_stack_from_slices(
                resampling_grid=recon_grid_sitk, interpolator="BSpline")
        stack_motioncorrectedic_recon_grid.write(
            args.dir_output, filename_stack + "_motion-corrected-ic")

    # verbose:
    if args.verbose:
        tmp = [
            stack_naivelyscaled_recon_grid,
            stack_motioncorrected_recon_grid,
        ]
        if args.intensity_correction:
            tmp.append(stack_naivelyscaledic_recon_grid)
            tmp.append(stack_motioncorrectedic_recon_grid)
        tmp.append(reference_image_resampled)
        sitkh.show_stacks(tmp)
    # ---------------------------------------------------------------------
    # Extract mask from reference
    ph.print_title("Extract mask from reference")
    stack_masked = stack_intensityCorrected

    # ---------------------------------------------------------------------
    # Perform SR step
    ph.print_title("Perform In-plane Deconvolution Step")

    # Deconvolution only in-plane
    if args.sigma2 < 0:
        # Estimate PSF automatically given original in-plane resolution
        cov = psf.PSF().get_gaussian_psf_covariance_matrix_from_spacing(
            stack0.sitk.GetSpacing())
        cov[2, 2] = 1e-5
    else:
        # Use predefined covariance
        cov = np.array([args.sigma2, args.sigma2, 1e-5])

    if args.regularization != "TV":
        volumetric_recon = tk.TikhonovSolver(
            stacks=[stack_masked],
            reconstruction=recon_grid,
            alpha=args.alpha,
            iter_max=args.iter_max,
            deconvolution_mode="predefined_covariance",
            predefined_covariance=cov,
        )

    else:
        # Initial value
        volumetric_recon = tk.TikhonovSolver(
            stacks=[stack_masked],
            reconstruction=recon_grid,
            alpha=0.02,
            iter_max=5,
            deconvolution_mode="predefined_covariance",
            predefined_covariance=cov,
        )
        volumetric_recon.run()
        HR_volume0 = volumetric_recon.get_reconstruction()

        volumetric_recon = admm.ADMMSolver(
            stacks=[stack_masked],
            reconstruction=HR_volume0,
            alpha=args.alpha,
            iter_max=args.iter_max,
            deconvolution_mode="predefined_covariance",
            predefined_covariance=cov,
            rho=args.rho,
            iterations=args.iterations,
        )

    volumetric_recon.run()

    stack_reconstructed = volumetric_recon.get_reconstruction()
    stack_reconstructed.set_filename(
        volumetric_recon.get_setting_specific_filename(prefix="recon_"))

    if args.verbose:
        tmp = [
            stack_naivelyscaled_recon_grid,
            stack_motioncorrected_recon_grid,
        ]
        if args.intensity_correction:
            tmp.append(stack_naivelyscaledic_recon_grid)
            tmp.append(stack_motioncorrectedic_recon_grid)
        tmp.append(stack_reconstructed)
        tmp.append(reference_image_resampled)
        sitkh.show_stacks(tmp)

    # Write results
    stack_reconstructed.write(
        args.dir_output,
        stack.get_filename() + "_" + stack_reconstructed.get_filename())
    elapsed_time = ph.stop_timing(time_start)

    ph.print_title("Summary Motion Correction")
    ph.print_info("Computational time: %s" % elapsed_time)

    return 0


if __name__ == '__main__':
    main()
