# Models

Installed with the package and loaded by default (`rlmhelix.student.load_head`, `rlmhelix.Calibrator.load`).

| File | Description | Size | SHA-256 |
|---|---|---|---|
| `student_bilinear_apc.pt` | Distilled student head (bilinear, 2,167,040 parameters) on the frozen RiNALMo-giga backbone; epoch 13, validation pixel-wise Pearson correlation of APC maps 0.512. A checkpoint with `config` and `state_dict` (plus the metadata `source`, `epoch`, `validation_apc_pixel_pearson` and `training`) that loads with `torch.load(..., weights_only=True)`. | 8.7 MB | 993cdf533436383f66b49525c0c121b43dbc8f8f686252c714210f54cce81583 |
| `calibrator_gev.json` | Conditional GEV calibrator of the RLMHelix score: feature means and standard deviations, location and log-scale coefficients (25 features + intercept), shape ξ = −0.1845; fitted on 570,000 dinucleotide-preserving shuffles. | 2.9 kB | 0924e6760e330db63e864be0562616784936729fb0fd3d9fbfdf80e06e00965f |

The RNA language model (RiNALMo-giga, `multimolecule/rinalmo-giga`) is downloaded from Hugging Face on first use and is
not redistributed here. It is loaded at the revision used for the paper, `2a71f6f98fb41dd2e6542a5e131d3778111d1468`
(`rlmhelix.dependency.BACKBONE_REVISION`; the student head was trained on its hidden states, and the genome scan and the
benchmarks used it). `Student(revision=None)` and `TeacherModel(revision=None)` load the latest revision instead.
