#!/usr/bin/python

# \file reconstructVolume.py
#  \brief
#
#  \author Michael Ebner (michael.ebner.14@ucl.ac.uk)
#  \date Nov 2016


# Import libraries
import os
import argparse
import numpy as np
import SimpleITK as sitk

import pythonhelper.SimpleITKHelper as sitkh
import pythonhelper.PythonHelper as ph
import volumetricreconstruction.base.Stack as st
import volumetricreconstruction.preprocessing.BrainStripping as bs
import volumetricreconstruction.preprocessing.IntensityCorrection as ic
import volumetricreconstruction.reconstruction.solver.TikhonovSolver as tk
import volumetricreconstruction.reconstruction.solver.ADMMSolver as admm

import volumetricreconstructionfromprintedmrfilms.utilities as utils


def get_parsed_input_line(
    verbose,
    minimizer,
    prefix_output,
    regularization,
    alpha,
    rho,
    iter_max,
    admm_iterations,
    sigma2,
    resolution_processing,
    resolution_reconstruction,
):

    parser = argparse.ArgumentParser(description="Run motion correction")

    parser.add_argument('--reference',
                        required=True,
                        type=str,
                        help="Path to reference image (*.nii.gz or *.nii)",
                        )
    parser.add_argument('--stack',
                        required=True,
                        type=str,
                        help="Path to naively stacked data "
                        "(*.nii.gz or *.nii)",
                        )
    parser.add_argument('--dir-input',
                        required=True,
                        type=str,
                        help="Input directory where final motion correction "
                        "results are stored (obtained by 'correctMotion.py')",
                        )
    parser.add_argument('--dir-output',
                        required=True,
                        type=str,
                        help="Output directory to store volumetric "
                        "reconstruction results",
                        )
    parser.add_argument('--regularization',
                        type=str,
                        help="Type of regularization for inverse problem. "
                        "Possible choices are 'TK0','TK1' or 'TV' for zeroth"
                        " or first order Tikhonov or isotropic total variation"
                        " regularization, respectively. I.e. "
                        "R(x) = ||x||^2 for 'TK0', "
                        "R(x) = ||Dx||^2 for 'TK1', "
                        "R(x) = TV(x) for 'TV'."
                        "[default: %s]"
                        % (regularization), default=regularization)
    parser.add_argument('--alpha',
                        type=float,
                        help="Regularization parameter alpha to solve the "
                        "inverse problem for each slice k:  argmin_x"
                        "[0.5 * ||y_k - A(sigma^2) x_k||^2 + alpha * R(x_k)]. "
                        "Recommendations for alpha: TK0, TK1: 0.05, TV: 5. "
                        "[default: %g]" % (alpha), default=alpha)
    parser.add_argument('--sigma2',
                        type=float,
                        help="Covariance for blurring operator A(\sigma^2) "
                        "[default: %g]" % (sigma2), default=sigma2)
    parser.add_argument('--rho',
                        type=float,
                        help="Regularization parameter rho for augmented "
                        "Lagrangian term used for ADMM"
                        "[default: %g]" % (rho), default=rho)
    parser.add_argument('--minimizer',
                        type=str,
                        help="Choice of minimizer used for the inverse problem"
                        " associated to the volumetric reconstruction step. "
                        "Possible choices are 'lsmr' or 'L-BFGS-B'. "
                        "[default: %s]" % (minimizer), default=minimizer)
    parser.add_argument('--iter-max',
                        type=int,
                        help="Number of maximum iterations for the numerical "
                        "solver. [default: %s]" % (iter_max), default=iter_max)
    parser.add_argument('--admm-iterations',
                        type=int,
                        help="Number of ADMM iterations. [default: %s]"
                        % (admm_iterations), default=admm_iterations)
    parser.add_argument('--resolution-processing',
                        type=float,
                        help="In-plane resolution used for processing scanned "
                        "images. [default: %g]"
                        % (resolution_processing),
                        default=resolution_processing)
    parser.add_argument('--resolution-reconstruction',
                        type=float,
                        help="In-plane resolution used for reconstructing the"
                        " final volume."
                        "[default: %g]"
                        % (resolution_reconstruction),
                        default=resolution_reconstruction)
    parser.add_argument('--prefix-output',
                        type=str,
                        help="Prefix for volumetric reconstruction output "
                        "filename. [default: %s]"
                        % (prefix_output), default=prefix_output)
    parser.add_argument('--verbose',
                        type=int,
                        help="Turn on/off verbose output. "
                        "[default %s]" % (verbose),
                        default=verbose,
                        )

    args = parser.parse_args()

    ph.print_title("Given Input")
    print("Chosen Parameters:")
    for arg in sorted(vars(args)):
        ph.print_info("%s: " % (arg), newline=False)
        print(getattr(args, arg))

    return args

if __name__ == '__main__':

    args = get_parsed_input_line(
        verbose=1,
        resolution_processing=0.25,
        resolution_reconstruction=1,
        minimizer="lsmr",
        prefix_output="recon_",
        iter_max=10,
        sigma2=0.25,
        # regularization="TK1",
        # alpha=0.3,    # TK1
        regularization="TV",
        alpha=5,    # TV
        rho=0.5,
        admm_iterations=10,
    )

    time_start = ph.start_timing()

    # Create output directory
    # ph.create_directory(args.dir_output)

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
    slice_thickness = stack0.sitk.GetSpacing()[-1]

    resampling_grid_sitk = 0 * sitkh.get_downsampled_sitk_image(
        stack_corrected.sitk,
        new_spacing=(args.resolution_processing,
                     args.resolution_processing,
                     slice_thickness)
    )

    # Get enlarged FOV
    resampling_grid_sitk = sitkh.get_altered_field_of_view_sitk_image(
        resampling_grid_sitk, unit="mm",
        boundary_i=5,
        boundary_j=5,
        boundary_k=0)

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

    # ---------------------------------------------------------------------
    # Get brain mask for reference image
    ph.print_title("Get brain mask for reference image")
    brain_stripping = bs.BrainStripping()
    brain_stripping.set_input_image_sitk(reference_image.sitk)
    brain_stripping.run_stripping()
    reference_image_sitk_mask = brain_stripping.get_brain_mask_sitk()
    reference_image = st.Stack.from_sitk_image(
        reference_image.sitk,
        reference_image.get_filename(),
        reference_image_sitk_mask)

    # ---------------------------------------------------------------------
    # Resampling to processing and reconstruction grids
    ph.print_title("Resample motion-corrected stack to processing grid")
    default_pixel_value = np.percentile(
        np.array(sitk.GetArrayFromImage(stack.sitk)), 0.1)
    stack_resampled = stack.get_resampled_stack_from_slices(
        resampling_grid=resampling_grid_sitk,
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    reference_image_resampled = reference_image.get_resampled_stack(
        resampling_grid=resampling_grid_sitk, interpolator="BSpline")
    stack_resampled.set_filename(filename_stack+"_motion-corrected")

    ph.print_title("Resample original stack to reconstruction grid")
    stack0_resampled = st.Stack.from_stack(stack0)
    i = len(slice_transforms_sitk)/2
    stack0_resampled.update_motion_correction(slice_transforms_sitk[i])
    stack0_resampled = stack0_resampled.get_resampled_stack_from_slices(
        resampling_grid=recon_grid_sitk,
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    stack0_resampled.set_filename(filename_stack)

    # Write results
    stack_naivelyscaled_recon_grid = \
        stack0_resampled.get_resampled_stack_from_slices(
            resampling_grid=recon_grid_sitk, interpolator="BSpline")
    stack_naivelyscaled_recon_grid.write(
        args.dir_output, filename_stack+"_scaled")

    stack_motioncorrected_recon_grid = \
        stack_resampled.get_resampled_stack_from_slices(
            resampling_grid=recon_grid_sitk, interpolator="BSpline")
    stack_motioncorrected_recon_grid.write(
        args.dir_output, filename_stack+"_motion-corrected")

    # ---------------------------------------------------------------------
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
    intensity_correction.run_lower_percentile_capping_of_stack(percentile=25)
    # intensity_correction.use_individual_slice_correction(True)
    intensity_correction.run_linear_intensity_correction()
    stack_intensityCorrected = intensity_correction.get_intensity_corrected_stack()
    stack0_intensityCorrected = \
        intensity_correction.get_intensity_corrected_additional_stack()

    stack_intensityCorrected.set_filename(
        filename_stack+"_motion-corrected-ic")
    stack0_intensityCorrected.set_filename(filename_stack+"_scaled-ic")

    # Write results
    stack_naivelyscaledic_recon_grid = stack0_intensityCorrected
    stack_naivelyscaledic_recon_grid.write(
        args.dir_output, filename_stack+"_scaled-ic")

    stack_motioncorrectedic_recon_grid = \
        stack_intensityCorrected.get_resampled_stack_from_slices(
            resampling_grid=recon_grid_sitk, interpolator="BSpline")
    stack_motioncorrectedic_recon_grid.write(
        args.dir_output, filename_stack+"_motion-corrected-ic")

    # verbose:
    if args.verbose:

        # On reconstruction grid:
        sitkh.show_stacks(
            [stack_naivelyscaled_recon_grid,
             stack_motioncorrected_recon_grid,
             stack_naivelyscaledic_recon_grid,
             stack_motioncorrectedic_recon_grid,
             reference_image_resampled],
            # label=[
            # "NaivelyScaled",
            # "MotionCorrected",
            # "NaivelyScaledIC",
            # "MotionCorrectedIC",
            # "Reference"],
        )

    # ---------------------------------------------------------------------
    # Extract mask from reference
    ph.print_title("Extract mask from reference")
    stack_masked = stack_intensityCorrected

    # ---------------------------------------------------------------------
    # Perform SR step
    ph.print_title("Perform SR step")

    # Deconvolution only in-plane
    cov = np.array([args.sigma2, args.sigma2, 1e-5])

    if args.regularization != "TV":
        volumetric_recon = tk.TikhonovSolver(
            stacks=[stack_masked],
            reconstruction=recon_grid,
            alpha=args.alpha,
            iter_max=args.iter_max,
            minimizer=args.minimizer,
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
            minimizer="lsmr",
            x_scale=1,
            deconvolution_mode="predefined_covariance",
            predefined_covariance=cov,
        )
        volumetric_recon.run_reconstruction()
        HR_volume0 = volumetric_recon.get_reconstruction()

        volumetric_recon = admm.ADMMSolver(
            stacks=[stack_masked],
            reconstruction=HR_volume0,
            alpha=args.alpha,
            iter_max=args.iter_max,
            minimizer=args.minimizer,
            x_scale=1,
            deconvolution_mode="predefined_covariance",
            predefined_covariance=cov,
            rho=args.rho,
            iterations=args.admm_iterations,
        )

    volumetric_recon.run_reconstruction()
    volumetric_recon.compute_statistics()

    stack_reconstructed = volumetric_recon.get_reconstruction()
    stack_reconstructed.set_filename(
        volumetric_recon.get_setting_specific_filename(prefix=args.prefix_output))

    if args.verbose:
        sitkh.show_stacks(
            [stack_naivelyscaled_recon_grid,
             stack_motioncorrected_recon_grid,
             stack_naivelyscaledic_recon_grid,
             stack_motioncorrectedic_recon_grid,
             stack_reconstructed,
             reference_image_resampled],
            label=["NaivelyScaled",
                   "MotionCorrected",
                   "NaivelyScaledIC",
                   "MotionCorrectedIC",
                   "Recon",
                   "Reference"],
        )

    # Write results
    stack_reconstructed.write(
        args.dir_output,
        stack.get_filename() + "_" + stack_reconstructed.get_filename())
    elapsed_time = ph.stop_timing(time_start)

    ph.print_title("Summary Motion Correction")
    ph.print_info("Computational time: %s" % elapsed_time)
