# CUT&RUN

| Stage | Script | Purpose |
| --- | --- | --- |
| 04 | `04_align_and_normalize.py` | Align mouse and yeast reads, filter BAMs, and generate normalized bigWigs. |
| 05 | `05_call_peaks_and_define_binding.py` | Call and filter peaks, then define promoter-bound genes. |
| 06 | `06_render_figures.py` | Render CUT&RUN figures and generate promoter-priority peak-to-gene assignments. |

The peak-to-gene assignments from Step 06 are used by the optional Step 07
shared-peak analysis. CUT&RUN figure generation itself is unchanged by the
Step 07 `--shared-peaks` option.
