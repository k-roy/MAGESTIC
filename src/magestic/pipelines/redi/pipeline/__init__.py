"""magestic.pipelines.redi.pipeline — step orchestrators (step_a..step_i)."""

from . import (
    step_a_trim_counts,
    step_b_combine_keys,
    step_c_process_bcs,
    step_d_purity,
    step_e_annotate,
    step_f_cross_array,
    step_g_wgs_check,
    step_h_plot,
    step_i_select,
)

__all__ = [
    "step_a_trim_counts",
    "step_b_combine_keys",
    "step_c_process_bcs",
    "step_d_purity",
    "step_e_annotate",
    "step_f_cross_array",
    "step_g_wgs_check",
    "step_h_plot",
    "step_i_select",
]
