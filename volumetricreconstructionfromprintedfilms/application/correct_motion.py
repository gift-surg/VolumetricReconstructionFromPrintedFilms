##
# \file correct_motion.py
#
# \author     Michael Ebner (michael.ebner.14@ucl.ac.uk)
# \date       Aug 2016
#


import SimpleITK as sitk
import numpy as np
import os

import niftymic.base.stack as st
import niftymic.registration.cpp_itk_registration as regitk
import niftymic.registration.intra_stack_registration as intrareg
import niftymic.registration.niftyreg as regniftyreg
import niftymic.utilities.brain_stripping as bs
import niftymic.utilities.intensity_correction as ic
import pysitk.python_helper as ph
import pysitk.simple_itk_helper as sitkh
import volumetricreconstructionfromprintedfilms.utilities.input_argparser as inargs
import volumetricreconstructionfromprintedfilms.utilities.utilities as utils


def main():

    time_start = ph.start_timing()

    input_parser = inargs.InputArgparser(
        description="Run motion correction pipeline to estimate meta-data "
        "information and physical position of each slice in the 3D space. "
        "Estimated transform parameters for both in-plane similarity and "
        "affine transforms are written to the output directory for each "
        "single slice. "
        "The obtained motion correction results are used as input for "
        "'reconstruct_volume.py' which provides a volumetric "
        "reconstructions in a subsequent step. "

        "Intermediate results can be stored by using the optional "
        "(but recommended) argument 'dir-output-verbose'. ",
    )
    input_parser.add_stack(required=True)
    input_parser.add_reference(required=True)
    input_parser.add_dir_output(
        required=True,
        help="Output directory to store motion correction results for both "
        "similarity and affine in-plane transformations.",)
    input_parser.add_iter_max(default=20)
    input_parser.add_verbose(default=0)
    input_parser.add_dir_output_verbose(required=False)
    input_parser.add_factor_downsampling(default=2)
    input_parser.add_factor_inplane_spacing(default=1)
    input_parser.add_option(
        option_string="--rigid-only",
        type=int,
        default=0,
        help="Turn on/off the use of in-plane rigid motion correction only. "
        "Intensity correction is not performed then. Slice motion correction "
        "transforms are written to folder 'Rigid' in the output folder."
    )

    args = input_parser.parse_args()
    input_parser.print_arguments(args)

    if args.dir_output_verbose is None and args.verbose is True:
        raise IOError(
            "Provide --dir-output-verbose option in case you want to "
            "run verbose")

    # Extract filenames from given path without filename extension
    filename_reference = os.path.basename(args.reference).split(".")[0]
    filename_stack = os.path.basename(args.stack).split(".")[0]

    # ---------------------------------------------------------------------
    # Read reference image
    ph.print_title("Read Data")
    ph.print_info("Read reference image")
    reference_image_sitk = sitk.ReadImage(args.reference, sitk.sitkFloat64)

    # ---------------------------------------------------------------------
    # Read original stack of slices
    ph.print_info("Read original stack of slices")
    stack_sitk = sitk.ReadImage(args.stack, sitk.sitkFloat64)

    # Define original stack
    stack0 = st.Stack.from_sitk_image(stack_sitk, filename_stack)

    # ---------------------------------------------------------------------
    # Scale original stack with initial in-plane scale estimate
    spacing = np.array(stack_sitk.GetSpacing())
    spacing[0:-1] *= args.factor_inplane_spacing
    stack_sitk.SetSpacing(spacing)
    stack = st.Stack.from_sitk_image(stack_sitk, filename_stack)

    # ---------------------------------------------------------------------
    # Compensation transform for scaling w.r.t. to original stack0

    # Create scaling matrix + get info from image
    scaling_matrix = np.diag(
        [args.factor_inplane_spacing, args.factor_inplane_spacing, 1])
    direction_matrix = np.array(stack_sitk.GetDirection()).reshape(3, 3)
    origin = np.array(stack_sitk.GetOrigin())

    # Affine transform to correct for scaling in image header
    matrix = direction_matrix.dot(scaling_matrix).dot(
        direction_matrix.transpose())
    translation = origin - matrix.dot(origin)
    scaling_transform_sitk = sitk.AffineTransform(3)
    scaling_transform_sitk.SetMatrix(matrix.flatten())
    scaling_transform_sitk.SetTranslation(translation)

    # Create scaling transform incorporating that info
    slice_transforms_sitk = [scaling_transform_sitk] * \
        stack.get_number_of_slices()

    # ---------------------------------------------------------------------
    # Get brain mask for reference image
    ph.print_title("Get brain mask for reference image")
    brain_stripping = bs.BrainStripping()
    brain_stripping.set_input_image_sitk(reference_image_sitk)
    brain_stripping.run()
    reference_image_sitk_mask = brain_stripping.get_mask_around_skull(
        dilate_radius=0, erode_radius=0)
    reference_image = st.Stack.from_sitk_image(
        reference_image_sitk, filename_reference, reference_image_sitk_mask)

    # ---------------------------------------------------------------------
    # Downsample image
    ph.print_title("Downsample stack image")
    # noinspection PyTypeChecker
    default_pixel_value = np.percentile(
        np.array(sitk.GetArrayFromImage(stack.sitk)), 0.1)
    factor_downsampling = args.factor_downsampling / \
        float(args.factor_inplane_spacing)
    stack_sitk_downsampled = sitkh.get_downsampled_sitk_image(
        stack_sitk,
        downsampling_factors=(factor_downsampling,
                              factor_downsampling, 1),
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    ph.print_info("Default pixel value for resampling: %.2f"
                  % default_pixel_value)
    ph.print_info("Downsampling factor (corrected by in-plane spacing): %g"
                  % factor_downsampling)

    # ---------------------------------------------------------------------
    # Skull mask stripping
    ph.print_title("Skull mask stripping for stack image")
    brain_stripping = bs.BrainStripping(
        compute_brain_mask=True, compute_skull_image=True)
    brain_stripping.set_input_image_sitk(stack_sitk_downsampled)
    brain_stripping.run()
    # stack_sitk_mask_downsampled = brain_stripping.get_brain_mask_sitk()
    stack_sitk_mask_downsampled = brain_stripping.get_mask_around_skull(
        dilate_radius=0, erode_radius=0)
    stack_downsampled = st.Stack.from_sitk_image(
        stack_sitk_downsampled,
        filename_stack + "_downsampled",
        stack_sitk_mask_downsampled)

    # Counter to write the output images in a consecutive sequence
    ctr = [-1]

    # Write result
    if args.dir_output_verbose is not None:

        filename_suffix = "_downsampled" + \
            str(factor_downsampling).replace(".", "p")
        stack_downsampled.set_filename(filename_stack + filename_suffix)
        stack_downsampled.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=True)

    # ---------------------------------------------------------------------
    # Self in-plane rigid registration
    if not args.rigid_only:
        ph.print_title("Self in-plane rigid registration")
        inplane_registration = intrareg.IntraStackRegistration(
            stack=stack_downsampled,
            transform_initializer_type="moments",
            # use_stack_mask=True,
            optimizer_iter_max=args.iter_max,
            optimizer_loss="linear",
            use_verbose=True,
        )
        inplane_registration.run()
        inplane_registration.print_statistics()
        stack_selfRigidInplane = inplane_registration.get_corrected_stack()

        # Get slice transforms
        slice_transforms_sitk_update = \
            inplane_registration.get_slice_transforms_sitk()
        slice_transforms_sitk = utils.get_updated_affine_transforms(
            slice_transforms_sitk_update, slice_transforms_sitk)

        # Write result
        if args.dir_output_verbose is not None:
            filename_suffix = "_selfinplane"
            stack_selfRigidInplane.set_filename(
                filename_stack + filename_suffix)
            tmp = stack_selfRigidInplane.get_resampled_stack_from_slices(
                interpolator="BSpline", default_pixel_value=default_pixel_value)
            tmp.write(
                directory=args.dir_output_verbose,
                filename=filename_stack + "_" +
                str(ph.add_one(ctr)) + filename_suffix,
                write_mask=False)

        if args.verbose:
            sitkh.show_stacks(
                [stack_downsampled,
                 stack_selfRigidInplane.get_resampled_stack_from_slices(
                     resampling_grid=stack_downsampled.sitk,
                     interpolator="BSpline",
                     default_pixel_value=default_pixel_value),
                 ],
                segmentation=stack_downsampled
            )

    else:
        stack_selfRigidInplane = stack_downsampled

    # ---------------------------------------------------------------------
    # Skull mask stripping
    # ph.print_title("Skull mask stripping (for reference alignment)")

    tmp_fixed = stack_selfRigidInplane.get_resampled_stack_from_slices(
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)

    # ---------------------------------------------------------------------
    # Rigid registration to reference
    ph.print_title("Rigid registration to reference")
    # noinspection PyTypeChecker
    default_pixel_value = np.max((np.percentile(
        np.array(sitk.GetArrayFromImage(stack_downsampled.sitk)), 0.1), 0))

    registration_method = regniftyreg.RegAladin(
        fixed=tmp_fixed,
        moving=reference_image,
        use_verbose=False,
        # use_fixed_mask=True,
        # use_moving_mask=True,
    )
    registration_method.run()
    stack_rigidToReference = registration_method.get_transformed_fixed()

    # Get slice transforms
    slice_transforms_sitk_update = \
        registration_method.get_registration_transform_sitk()
    slice_transforms_sitk = utils.get_updated_affine_transforms(
        slice_transforms_sitk_update, slice_transforms_sitk)

    # Write result
    if args.dir_output_verbose is not None:
        filename_suffix = "_RigidToReference"
        stack_rigidToReference.set_filename(filename_stack + filename_suffix)
        tmp = stack_rigidToReference.get_resampled_stack_from_slices(
            interpolator="BSpline")
        tmp.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=False,
            # write_mask=True,
        )

    if args.verbose:
        sitkh.show_stacks(
            [stack_rigidToReference.get_resampled_stack_from_slices(
                resampling_grid=reference_image.sitk,
                interpolator="BSpline",
                default_pixel_value=default_pixel_value),
             reference_image])

    if args.rigid_only:
        # ---------------------------------------------------------------------
        # Inplane 2D rigid registration to reference
        ph.print_title("Inplane 2D rigid registration to reference")
        inplane_registration = intrareg.IntraStackRegistration(
            stack=stack_rigidToReference,
            reference=reference_image.get_resampled_stack(
                resampling_grid=stack_rigidToReference.sitk,
                interpolator="BSpline"))
        inplane_registration.use_verbose(True)
        inplane_registration.use_stack_mask_reference_fit_term(True)
        inplane_registration.use_stack_mask_neighbour_fit_term(True)
        inplane_registration.use_reference_mask(True)
        inplane_registration.set_transform_type("rigid")
        # inplane_registration.set_image_transform_reference_fit_term("gradient_magnitude")
        # inplane_registration.set_image_transform_reference_fit_term("partial_derivative")
        inplane_registration.set_intensity_correction_initializer_type(None)
        inplane_registration.set_intensity_correction_type_reference_fit(
            "affine")
        inplane_registration.set_intensity_correction_type_slice_neighbour_fit(
            "affine")
        inplane_registration.set_optimizer_iter_max(args.iter_max)
        inplane_registration.set_alpha_reference(10)
        inplane_registration.set_alpha_neighbour(1)
        inplane_registration.set_alpha_parameter(1e6)
        inplane_registration.set_optimizer_loss("soft_l1")
        inplane_registration.set_transform_initializer_type("identity")
        inplane_registration.run()
        inplane_registration.print_statistics()
        stack_RigidInplane = inplane_registration.get_corrected_stack()

        filename_suffix = "_inplane2Drigid"
        stack_RigidInplane.set_filename(filename_stack + filename_suffix)

        # Get slice transforms
        slice_transforms_sitk_update = \
            inplane_registration.get_slice_transforms_sitk()
        slice_transforms_sitk = utils.get_updated_affine_transforms(
            slice_transforms_sitk_update, slice_transforms_sitk)

        # Write result
        if args.dir_output_verbose is not None:
            tmp = stack_RigidInplane.get_resampled_stack_from_slices(
                interpolator="BSpline")
            tmp.write(
                directory=args.dir_output_verbose,
                filename=filename_stack + "_" +
                str(ph.add_one(ctr)) + filename_suffix,
                write_mask=False)

        if args.verbose:
            sitkh.show_stacks([
                stack_rigidToReference.get_resampled_stack_from_slices(
                    resampling_grid=stack_rigidToReference.sitk,
                    interpolator="BSpline",
                    default_pixel_value=default_pixel_value),
                stack_RigidInplane.get_resampled_stack_from_slices(
                    resampling_grid=stack_rigidToReference.sitk,
                    interpolator="BSpline",
                    default_pixel_value=default_pixel_value),
                reference_image,
            ]
            )

        # ---------------------------------------------------------------------
        # Write results: In-plane 2D Similar
        stack_final = stack_RigidInplane
        utils.write_results_motion_correction(
            os.path.join(args.dir_output, "Rigid"),
            filename_stack,
            stack0,
            stack_final,
            slice_transforms_sitk,
            reference_image)

        elapsed_time = ph.stop_timing(time_start)
        ph.print_title("Summary Motion Correction")
        ph.print_info("Computational time: %s" % elapsed_time)

        return 0

    # ---------------------------------------------------------------------
    # Inplane 3D similarity registration to reference
    ph.print_title("Inplane 3D similarity registration to reference")
    tmp_fixed = stack_rigidToReference.get_resampled_stack_from_slices(
        interpolator="BSpline",
        default_pixel_value=default_pixel_value)
    registration_itk = regitk.CppItkRegistration(
        fixed=tmp_fixed,
        moving=reference_image,
        registration_type="InplaneSimilarity",
        # use_fixed_mask=True,
        # use_moving_mask=True,
        use_multiresolution_framework=True,
        interpolator="Linear",
        scales_estimator="PhysicalShift",
        # scales_estimator="Jacobian",
        metric="Correlation",
        # metric="MattesMutualInformation",
        use_verbose=False,
    )
    registration_itk.run()
    stack_inplane3DSimilar = \
        registration_itk.get_stack_with_similarity_inplane_transformed_slices(
            stack_rigidToReference)

    # Get uniform in-plane scaling factor
    inplane_scale_3D = registration_itk.get_parameters()[6]
    ph.print_info("inplane_scale_3D = %g" % inplane_scale_3D)

    # Get all affine transforms to keep track of corrections
    # (Important: Affine transforms incorporate scaling!)
    slice_transforms_sitk_update = \
        registration_itk.get_registration_transform_sitk()
    slice_transforms_sitk = utils.get_updated_affine_transforms(
        slice_transforms_sitk_update, slice_transforms_sitk)

    # Write result
    if args.dir_output_verbose is not None:
        filename_suffix = "_inplane3Dsimilar"
        stack_inplane3DSimilar.set_filename(filename_stack + filename_suffix)
        tmp = stack_inplane3DSimilar.get_resampled_stack_from_slices(
            interpolator="BSpline")
        tmp.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=False)

    if args.verbose:
        sitkh.show_stacks([
            stack_inplane3DSimilar.get_resampled_stack_from_slices(
                resampling_grid=reference_image.sitk,
                interpolator="BSpline",
                default_pixel_value=default_pixel_value,
            ),
            reference_image])

    # ---------------------------------------------------------------------
    # Resampling of reference to stack grid
    ph.print_title("Resampling of reference to stack grid")
    reference_image_downsampled = reference_image.get_resampled_stack(
        resampling_grid=stack_inplane3DSimilar.sitk, interpolator="BSpline")
    if args.dir_output_verbose is not None:
        reference_image_downsampled.write(
            directory=args.dir_output_verbose, write_mask=True)

    # ---------------------------------------------------------------------
    # Skull mask stripping
    ph.print_title("Skull mask stripping")
    stack_tmp_sitk = stack_inplane3DSimilar.get_resampled_stack_from_slices(
        interpolator="BSpline", default_pixel_value=default_pixel_value).sitk
    brain_stripping = bs.BrainStripping(
        compute_brain_mask=True, compute_skull_image=True)
    brain_stripping.set_input_image_sitk(stack_tmp_sitk)
    brain_stripping.run()
    stack_sitk_mask_downsampled = brain_stripping.get_mask_around_skull(
        dilate_radius=10, erode_radius=0)
    stack_tmp = st.Stack.from_sitk_image(
        stack_tmp_sitk,
        filename_stack + "_inplane3Dsimilar",
        stack_sitk_mask_downsampled)

    # ---------------------------------------------------------------------
    # Perform intensity correction
    ph.print_title("Perform intensity correction")
    intensity_correction = ic.IntensityCorrection(
        stack=stack_tmp, reference=reference_image_downsampled, use_verbose=1)
    intensity_correction.use_individual_slice_correction(False)
    intensity_correction.run_affine_intensity_correction()
    intensity_correction.run_lower_percentile_capping_of_stack(percentile=20)
    intensity_correction.use_individual_slice_correction(False)
    intensity_correction.run_linear_intensity_correction()
    stack_intensityCorrected = \
        intensity_correction.get_intensity_corrected_stack()

    # Write result
    if args.dir_output_verbose is not None:
        filename_suffix = "_intensityCorrected"
        stack_intensityCorrected.set_filename(filename_stack + filename_suffix)
        tmp = stack_intensityCorrected.get_resampled_stack_from_slices(
            interpolator="BSpline")
        tmp.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=False)

    if args.verbose:
        sitkh.show_stacks(
            [stack_intensityCorrected.get_resampled_stack_from_slices(
                resampling_grid=reference_image_downsampled.sitk,
                interpolator="BSpline"),
             reference_image_downsampled],
            segmentation=reference_image_downsampled)

    # ---------------------------------------------------------------------
    # Inplane 2D similarity registration to reference
    ph.print_title("Inplane 2D similarity registration to reference")
    inplane_registration = intrareg.IntraStackRegistration(
        stack=stack_intensityCorrected,
        reference=reference_image_downsampled)
    inplane_registration.use_verbose(True)
    inplane_registration.use_stack_mask_reference_fit_term(True)
    inplane_registration.use_stack_mask_neighbour_fit_term(True)
    inplane_registration.use_reference_mask(True)
    inplane_registration.set_transform_type("similarity")
    # inplane_registration.set_image_transform_reference_fit_term("gradient_magnitude")
    # inplane_registration.set_image_transform_reference_fit_term("partial_derivative")
    inplane_registration.set_intensity_correction_initializer_type(None)
    inplane_registration.set_intensity_correction_type_reference_fit("affine")
    inplane_registration.set_intensity_correction_type_slice_neighbour_fit(
        "affine")
    inplane_registration.set_optimizer_iter_max(args.iter_max)
    inplane_registration.set_alpha_reference(10)
    inplane_registration.set_alpha_neighbour(1)
    inplane_registration.set_alpha_parameter(1e6)
    inplane_registration.set_optimizer_loss("soft_l1")
    # inplane_registration.set_optimizer_loss("linear")
    # inplane_registration.set_optimizer_loss("huber")

    # Variant A: Identity initialization
    inplane_registration.set_transform_initializer_type("identity")
    inplane_registration.run()
    inplane_registration.print_statistics()
    final_cost = inplane_registration.get_final_cost()
    slice_transforms_sitk_update = \
        inplane_registration.get_slice_transforms_sitk()
    stack_inplane2Dsimilar = inplane_registration.get_corrected_stack()
    filename_suffix = inplane_registration.get_setting_specific_filename()
    stack_inplane2Dsimilar.set_filename(filename_stack + filename_suffix)

    # sitkh.show_stacks(
    #     [stack_inplane3DSimilar.get_resampled_stack_from_slices(
    #         interpolator="BSpline"),
    #      stack_inplane2Dsimilar.get_resampled_stack_from_slices(
    #         resampling_grid=stack_inplane3DSimilar.sitk, interpolator="BSpline"),
    #      reference_image_downsampled
    #      ],
    #     segmentation=reference_image_downsampled
    # )

    # Variant B: Moments initialization
    inplane_registration.set_transform_initializer_type("moments")
    inplane_registration.run()
    inplane_registration.print_statistics()

    # Pick best fit between variant A and B
    if inplane_registration.get_final_cost() < final_cost:
        print("In-plane registration based on moments has smaller cost " +
              "(%.3e < %.3e)" % (inplane_registration.get_final_cost(),
                                 final_cost))
        slice_transforms_sitk_update = \
            inplane_registration.get_slice_transforms_sitk()
        stack_inplane2Dsimilar = inplane_registration.get_corrected_stack()
        filename_suffix = inplane_registration.get_setting_specific_filename()
        stack_inplane2Dsimilar.set_filename(filename_stack + filename_suffix)
    else:
        print("In-plane registration based on identity has smaller cost " +
              "(%.3e < %.3e)" % (final_cost,
                                 inplane_registration.get_final_cost()))

    if args.verbose:
        sitkh.show_stacks(
            [stack_inplane3DSimilar.get_resampled_stack_from_slices(
                interpolator="BSpline"),
             stack_inplane2Dsimilar.get_resampled_stack_from_slices(
                resampling_grid=stack_inplane3DSimilar.sitk,
                interpolator="BSpline"),
             reference_image_downsampled
             ],
            segmentation=reference_image_downsampled
        )

    # Get slice transforms
    slice_transforms_sitk = utils.get_updated_affine_transforms(
        slice_transforms_sitk_update, slice_transforms_sitk)

    # Write result
    if args.dir_output_verbose is not None:
        tmp = stack_inplane2Dsimilar.get_resampled_stack_from_slices(
            interpolator="BSpline")
        tmp.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=False)

    # ---------------------------------------------------------------------
    # Write results: In-plane 2D Similar
    stack_final = stack_inplane2Dsimilar
    utils.write_results_motion_correction(
        os.path.join(args.dir_output, "Similarity"),
        filename_stack,
        stack0,
        stack_final,
        slice_transforms_sitk,
        reference_image)

    # ---------------------------------------------------------------------
    # Inplane affine registration to reference
    print("Inplane affine registration to reference")
    inplane_registration = intrareg.IntraStackRegistration(
        stack_inplane2Dsimilar, reference_image_downsampled)
    inplane_registration.set_transform_initializer_type("identity")
    # inplane_registration.use_parameter_normalization(True)
    inplane_registration.use_verbose(True)
    # inplane_registration.use_stack_mask(True)
    inplane_registration.use_stack_mask_reference_fit_term(True)
    inplane_registration.use_stack_mask_neighbour_fit_term(True)
    inplane_registration.use_reference_mask(True)
    inplane_registration.set_transform_type("affine")
    # inplane_registration.set_image_transform_reference_fit_term("gradient_magnitude")
    # inplane_registration.set_image_transform_reference_fit_term("partial_derivative")
    inplane_registration.set_intensity_correction_initializer_type(None)
    inplane_registration.set_intensity_correction_type_reference_fit("affine")
    inplane_registration.set_intensity_correction_type_slice_neighbour_fit(
        "affine")
    inplane_registration.set_optimizer_iter_max(args.iter_max)
    inplane_registration.set_alpha_reference(10)
    inplane_registration.set_alpha_neighbour(1)
    inplane_registration.set_alpha_parameter(1e3)
    inplane_registration.set_optimizer_loss("soft_l1")
    inplane_registration.run()
    inplane_registration.print_statistics()
    stack_inplane2Daffine = inplane_registration.get_corrected_stack()

    filename_suffix = inplane_registration.get_setting_specific_filename()
    stack_inplane2Daffine.set_filename(filename_stack + filename_suffix)

    if args.verbose:
        sitkh.show_stacks(
            [stack_inplane3DSimilar.get_resampled_stack_from_slices(
                resampling_grid=stack_inplane3DSimilar.sitk,
                interpolator="BSpline"),
             stack_inplane2Dsimilar.get_resampled_stack_from_slices(
                 resampling_grid=stack_inplane3DSimilar.sitk,
                 interpolator="BSpline"),
             stack_inplane2Daffine.get_resampled_stack_from_slices(
                 resampling_grid=stack_inplane3DSimilar.sitk,
                 interpolator="BSpline"),
             reference_image_downsampled],
            segmentation=reference_image_downsampled
        )

    # Get slice transforms
    slice_transforms_sitk_update = \
        inplane_registration.get_slice_transforms_sitk()
    slice_transforms_sitk = utils.get_updated_affine_transforms(
        slice_transforms_sitk_update, slice_transforms_sitk)

    # Write result
    if args.dir_output_verbose is not None:
        tmp = stack_inplane2Daffine.get_resampled_stack_from_slices(
            interpolator="BSpline")
        tmp.write(
            directory=args.dir_output_verbose,
            filename=filename_stack + "_" +
            str(ph.add_one(ctr)) + filename_suffix,
            write_mask=False)

    # ---------------------------------------------------------------------
    # Write results: In-plane 2D Affine
    stack_final = stack_inplane2Daffine
    utils.write_results_motion_correction(
        os.path.join(args.dir_output, "Affine"),
        filename_stack,
        stack0,
        stack_final,
        slice_transforms_sitk,
        reference_image)

    elapsed_time = ph.stop_timing(time_start)
    ph.print_title("Summary Motion Correction")
    ph.print_info("Computational time: %s" % elapsed_time)

    return 0


if __name__ == '__main__':
    main()
