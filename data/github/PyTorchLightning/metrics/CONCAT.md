# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**Note: we move fast, but still we preserve 0.1 version (one feature release) back compatibility.**


## [UnReleased] - 2022-MM-DD

### Added

- Added new image metric `SpectralAngleMapper` ([#885](https://github.com/PyTorchLightning/metrics/pull/885))


- Added `CoverageError` to classification metrics ([#787](https://github.com/PyTorchLightning/metrics/pull/787))


- Added `LabelRankingAveragePrecision` and `LabelRankingLoss` ([#787](https://github.com/PyTorchLightning/metrics/pull/787))


- Added support for `MetricCollection` in `MetricTracker` ([#718](https://github.com/PyTorchLightning/metrics/pull/718))


- Added new image metric `ErrorRelativeGlobalDimensionlessSynthesis` ([#894](https://github.com/PyTorchLightning/metrics/pull/894))


- Added new image metric `UniversalImageQualityIndex` ([#824](https://github.com/PyTorchLightning/metrics/pull/824))


- Added support for 3D image and uniform kernel in `StructuralSimilarityIndexMeasure` ([#818](https://github.com/PyTorchLightning/metrics/pull/818))


- Added smart update of `MetricCollection` ([#709](https://github.com/PyTorchLightning/metrics/pull/709))


- Added `ClasswiseWrapper` for better logging of classification metrics with multiple output values ([#832](https://github.com/PyTorchLightning/metrics/pull/832))


- Added `**kwargs` argument for passing additional arguments to base class ([#833](https://github.com/PyTorchLightning/metrics/pull/833))


- Added negative `ignore_index` for the Accuracy metric ([#362](https://github.com/PyTorchLightning/metrics/pull/362))


- Added `SpectralDistortionIndex` metric to image domain ([#873](https://github.com/PyTorchLightning/metrics/pull/873))


- Added `adaptive_k` for the `RetrievalPrecision` metric ([#910](https://github.com/PyTorchLightning/metrics/pull/910))


- Added `reset_real_features` argument image quality assesment metrics ([#722](https://github.com/PyTorchLightning/metrics/pull/722))


- Added new keyword argument `compute_on_cpu` to all metrics ([#867](https://github.com/PyTorchLightning/metrics/pull/867))


- Added `WeightedMeanAbsolutePercentageError` to regression package ([#948](https://github.com/PyTorchLightning/metrics/pull/948))


### Changed

- Made `num_classes` in `jaccard_index` a required argument ([#853](https://github.com/PyTorchLightning/metrics/pull/853), [#914](https://github.com/PyTorchLightning/metrics/pull/914))


- Added normalizer, tokenizer to ROUGE metric ([#838](https://github.com/PyTorchLightning/metrics/pull/838))


- Improved shape checking of `permutation_invariant_training` ([#864](https://github.com/PyTorchLightning/metrics/pull/864))


- Allowed reduction `None` ([#891](https://github.com/PyTorchLightning/metrics/pull/891))


- `MetricTracker.best_metric` will now give a warning when computing on metric that do not have a best ([#913](https://github.com/PyTorchLightning/metrics/pull/913))


### Deprecated

- Deprecated argument `compute_on_step` ([#792](https://github.com/PyTorchLightning/metrics/pull/792))


- Deprecated passing in `dist_sync_on_step`, `process_group`, `dist_sync_fn` direct argument ([#833](https://github.com/PyTorchLightning/metrics/pull/833))


### Removed

- Removed support for versions of Lightning lower than v1.5 ([#788](https://github.com/PyTorchLightning/metrics/pull/788))


- Removed deprecated functions, and warnings in Text ([#773](https://github.com/PyTorchLightning/metrics/pull/773))
  * `functional.wer`
  * `WER`
- Removed deprecated functions and warnings in Image ([#796](https://github.com/PyTorchLightning/metrics/pull/796))
  * `functional.ssim`
  * `functional.psnr`
  * `SSIM`
  * `PSNR`


- Removed deprecated functions, and warnings in classification and regression ([#806](https://github.com/PyTorchLightning/metrics/pull/806))
  * `FBeta`
  * `F1`
  * `Hinge`
  * `IoU`
  * `functional.iou`
  * `MatthewsCorrcoef`
  * `PearsonCorrcoef`
  * `SpearmanCorrcoef`
  * `functional.fbeta`
  * `functional.f1`
  * `functional.hinge`



- Removed deprecated functions, and warnings in detection and pairwise ([#804](https://github.com/PyTorchLightning/metrics/pull/804))
  * `MAP`
  * `functional.pairwise.manhatten`


- Removed deprecated functions, and warnings in Audio ([#805](https://github.com/PyTorchLightning/metrics/pull/805))
  * `PESQ`
  * `PIT`
  * `SDR`
  * `SNR`
  * `STOI`
  * `functional.audio.pesq`
  * `functional.audio.pit`
  * `functional.audio.sdr`
  * `functional.audio.snr`
  * `functional.audio.stoi`
  * `functional.audio.si_sdr`
  * `functional.audio.si_snr`


### Fixed

- Fixed device mismatch for `MAP` metric in specific cases ([#950](https://github.com/PyTorchLightning/metrics/pull/950))

- Improved testing speed ([#820](https://github.com/PyTorchLightning/metrics/pull/820))


- Fixed compatibility of `ClasswiseWrapper` with the `prefix` argument of `MetricCollection` ([#843](https://github.com/PyTorchLightning/metrics/pull/843))


- Fixed `BestScore` on GPU ([#912](https://github.com/PyTorchLightning/metrics/pull/912))

- Fixed Lsum computation for `ROUGEScore` ([#944](https://github.com/PyTorchLightning/metrics/pull/944))


## [0.7.3] - 2022-03-23

### Fixed

- Fixed unsafe log operation in `TweedieDeviace` for power=1 ([#847](https://github.com/PyTorchLightning/metrics/pull/847))
- Fixed bug in MAP metric related to either no ground truth or no predictions ([#884](https://github.com/PyTorchLightning/metrics/pull/884))
- Fixed `ConfusionMatrix`, `AUROC` and `AveragePrecision` on GPU when running in deterministic mode ([#900](https://github.com/PyTorchLightning/metrics/pull/900))
- Fixed NaN or Inf results returned by `signal_distortion_ratio` ([#899](https://github.com/PyTorchLightning/metrics/pull/899))
- Fixed memory leak when using `update` method with tensor where `requires_grad=True` ([#902](https://github.com/PyTorchLightning/metrics/pull/902))


## [0.7.2] - 2022-02-10

### Fixed

- Minor patches in JOSS paper.


## [0.7.1] - 2022-02-03

### Changed

- Used `torch.bucketize` in calibration error when `torch>1.8` for faster computations ([#769](https://github.com/PyTorchLightning/metrics/pull/769))
- Improve mAP performance ([#742](https://github.com/PyTorchLightning/metrics/pull/742))

### Fixed

- Fixed check for available modules ([#772](https://github.com/PyTorchLightning/metrics/pull/772))
- Fixed Matthews correlation coefficient when the denominator is 0 ([#781](https://github.com/PyTorchLightning/metrics/pull/781))


## [0.7.0] - 2022-01-17

### Added

- Added NLP metrics:
  - `MatchErrorRate` ([#619](https://github.com/PyTorchLightning/metrics/pull/619))
  - `WordInfoLost` and `WordInfoPreserved` ([#630](https://github.com/PyTorchLightning/metrics/pull/630))
  - `SQuAD` ([#623](https://github.com/PyTorchLightning/metrics/pull/623))
  - `CHRFScore` ([#641](https://github.com/PyTorchLightning/metrics/pull/641))
  - `TranslationEditRate` ([#646](https://github.com/PyTorchLightning/metrics/pull/646))
  - `ExtendedEditDistance` ([#668](https://github.com/PyTorchLightning/metrics/pull/668))
- Added `MultiScaleSSIM` into image metrics ([#679](https://github.com/PyTorchLightning/metrics/pull/679))
- Added Signal to Distortion Ratio (`SDR`) to audio package ([#565](https://github.com/PyTorchLightning/metrics/pull/565))
- Added `MinMaxMetric` to wrappers ([#556](https://github.com/PyTorchLightning/metrics/pull/556))
- Added `ignore_index` to retrieval metrics ([#676](https://github.com/PyTorchLightning/metrics/pull/676))
- Added support for multi references in `ROUGEScore` ([#680](https://github.com/PyTorchLightning/metrics/pull/680))
- Added a default VSCode devcontainer configuration ([#621](https://github.com/PyTorchLightning/metrics/pull/621))

### Changed

- Scalar metrics will now consistently have additional dimensions squeezed ([#622](https://github.com/PyTorchLightning/metrics/pull/622))
- Metrics having third party dependencies removed from global import ([#463](https://github.com/PyTorchLightning/metrics/pull/463))
- Untokenized for `BLEUScore` input stay consistent with all the other text metrics ([#640](https://github.com/PyTorchLightning/metrics/pull/640))
- Arguments reordered for `TER`, `BLEUScore`, `SacreBLEUScore`, `CHRFScore` now expect input order as predictions first and target second ([#696](https://github.com/PyTorchLightning/metrics/pull/696))
- Changed dtype of metric state from `torch.float` to `torch.long` in `ConfusionMatrix` to accommodate larger values ([#715](https://github.com/PyTorchLightning/metrics/pull/715))
- Unify `preds`, `target` input argument's naming across all text metrics ([#723](https://github.com/PyTorchLightning/metrics/pull/723), [#727](https://github.com/PyTorchLightning/metrics/pull/727))
  * `bert`, `bleu`, `chrf`, `sacre_bleu`, `wip`, `wil`, `cer`, `ter`, `wer`, `mer`, `rouge`, `squad`

### Deprecated

- Renamed IoU -> Jaccard Index ([#662](https://github.com/PyTorchLightning/metrics/pull/662))
- Renamed text WER metric ([#714](https://github.com/PyTorchLightning/metrics/pull/714))
  * `functional.wer` -> `functional.word_error_rate`
  * `WER` -> `WordErrorRate`
- Renamed correlation coefficient classes: ([#710](https://github.com/PyTorchLightning/metrics/pull/710))
  * `MatthewsCorrcoef` -> `MatthewsCorrCoef`
  * `PearsonCorrcoef` -> `PearsonCorrCoef`
  * `SpearmanCorrcoef` -> `SpearmanCorrCoef`
- Renamed audio STOI metric: ([#753](https://github.com/PyTorchLightning/metrics/pull/753), [#758](https://github.com/PyTorchLightning/metrics/pull/758))
  * `audio.STOI` to `audio.ShortTimeObjectiveIntelligibility`
  * `functional.audio.stoi` to `functional.audio.short_time_objective_intelligibility`
- Renamed audio PESQ metrics: ([#751](https://github.com/PyTorchLightning/metrics/pull/751))
  * `functional.audio.pesq` -> `functional.audio.perceptual_evaluation_speech_quality`
  * `audio.PESQ` -> `audio.PerceptualEvaluationSpeechQuality`
- Renamed audio SDR metrics: ([#711](https://github.com/PyTorchLightning/metrics/pull/711))
  * `functional.sdr` -> `functional.signal_distortion_ratio`
  * `functional.si_sdr` -> `functional.scale_invariant_signal_distortion_ratio`
  * `SDR` -> `SignalDistortionRatio`
  * `SI_SDR` -> `ScaleInvariantSignalDistortionRatio`
- Renamed audio SNR metrics: ([#712](https://github.com/PyTorchLightning/metrics/pull/712))
  * `functional.snr` -> `functional.signal_distortion_ratio`
  * `functional.si_snr` -> `functional.scale_invariant_signal_noise_ratio`
  * `SNR` -> `SignalNoiseRatio`
  * `SI_SNR` -> `ScaleInvariantSignalNoiseRatio`
- Renamed F-score metrics: ([#731](https://github.com/PyTorchLightning/metrics/pull/731), [#740](https://github.com/PyTorchLightning/metrics/pull/740))
  * `functional.f1` ->  `functional.f1_score`
  * `F1` ->  `F1Score`
  * `functional.fbeta` ->  `functional.fbeta_score`
  * `FBeta` ->  `FBetaScore`
- Renamed Hinge metric: ([#734](https://github.com/PyTorchLightning/metrics/pull/734))
  * `functional.hinge` ->  `functional.hinge_loss`
  * `Hinge` ->  `HingeLoss`
- Renamed image PSNR metrics ([#732](https://github.com/PyTorchLightning/metrics/pull/732))
  * `functional.psnr` -> `functional.peak_signal_noise_ratio`
  * `PSNR` -> `PeakSignalNoiseRatio`
- Renamed image PIT metric: ([#737](https://github.com/PyTorchLightning/metrics/pull/737))
  * `functional.pit` ->  `functional.permutation_invariant_training`
  * `PIT` ->  `PermutationInvariantTraining`
- Renamed image SSIM metric: ([#747](https://github.com/PyTorchLightning/metrics/pull/747))
  * `functional.ssim` ->  `functional.scale_invariant_signal_noise_ratio`
  * `SSIM` ->  `StructuralSimilarityIndexMeasure`
- Renamed detection `MAP` to `MeanAveragePrecision` metric ([#754](https://github.com/PyTorchLightning/metrics/pull/754))
- Renamed Fidelity & LPIPS image metric: ([#752](https://github.com/PyTorchLightning/metrics/pull/752))
  * `image.FID` ->  `image.FrechetInceptionDistance`
  * `image.KID` ->  `image.KernelInceptionDistance`
  * `image.LPIPS` ->  `image.LearnedPerceptualImagePatchSimilarity`

### Removed

- Removed `embedding_similarity` metric ([#638](https://github.com/PyTorchLightning/metrics/pull/638))
- Removed argument `concatenate_texts` from `wer` metric ([#638](https://github.com/PyTorchLightning/metrics/pull/638))
- Removed arguments `newline_sep` and `decimal_places` from `rouge` metric ([#638](https://github.com/PyTorchLightning/metrics/pull/638))

### Fixed

- Fixed MetricCollection kwargs filtering when no `kwargs` are present in update signature ([#707](https://github.com/PyTorchLightning/metrics/pull/707))


## [0.6.2] - 2021-12-15

### Fixed

- Fixed `torch.sort` currently does not support bool `dtype` on CUDA ([#665](https://github.com/PyTorchLightning/metrics/pull/665))
- Fixed mAP properly checks if ground truths are empty ([#684](https://github.com/PyTorchLightning/metrics/pull/684))
- Fixed initialization of tensors to be on correct device for `MAP` metric ([#673](https://github.com/PyTorchLightning/metrics/pull/673))


## [0.6.1] - 2021-12-06

### Changed

- Migrate MAP metrics from pycocotools to PyTorch ([#632](https://github.com/PyTorchLightning/metrics/pull/632))
- Use `torch.topk` instead of `torch.argsort` in retrieval precision for speedup ([#627](https://github.com/PyTorchLightning/metrics/pull/627))

### Fixed

- Fix empty predictions in MAP metric ([#594](https://github.com/PyTorchLightning/metrics/pull/594), [#610](https://github.com/PyTorchLightning/metrics/pull/610), [#624](https://github.com/PyTorchLightning/metrics/pull/624))
- Fix edge case of AUROC with `average=weighted` on GPU ([#606](https://github.com/PyTorchLightning/metrics/pull/606))
- Fixed `forward` in compositional metrics ([#645](https://github.com/PyTorchLightning/metrics/pull/645))


## [0.6.0] - 2021-10-28

### Added

- Added audio metrics:
  - Perceptual Evaluation of Speech Quality (PESQ) ([#353](https://github.com/PyTorchLightning/metrics/pull/353))
  - Short-Time Objective Intelligibility (STOI) ([#353](https://github.com/PyTorchLightning/metrics/pull/353))
- Added Information retrieval metrics:
  - `RetrievalRPrecision` ([#577](https://github.com/PyTorchLightning/metrics/pull/577))
  - `RetrievalHitRate` ([#576](https://github.com/PyTorchLightning/metrics/pull/576))
- Added NLP metrics:
  - `SacreBLEUScore` ([#546](https://github.com/PyTorchLightning/metrics/pull/546))
  - `CharErrorRate` ([#575](https://github.com/PyTorchLightning/metrics/pull/575))
- Added other metrics:
  - Tweedie Deviance Score ([#499](https://github.com/PyTorchLightning/metrics/pull/499))
  - Learned Perceptual Image Patch Similarity (LPIPS) ([#431](https://github.com/PyTorchLightning/metrics/pull/431))
- Added `MAP` (mean average precision) metric to new detection package ([#467](https://github.com/PyTorchLightning/metrics/pull/467))
- Added support for float targets in `nDCG` metric ([#437](https://github.com/PyTorchLightning/metrics/pull/437))
- Added `average` argument to `AveragePrecision` metric for reducing multi-label and multi-class problems ([#477](https://github.com/PyTorchLightning/metrics/pull/477))
- Added `MultioutputWrapper` ([#510](https://github.com/PyTorchLightning/metrics/pull/510))
- Added metric sweeping:
  - `higher_is_better` as constant attribute ([#544](https://github.com/PyTorchLightning/metrics/pull/544))
  - `higher_is_better` to rest of codebase ([#584](https://github.com/PyTorchLightning/metrics/pull/584))
- Added simple aggregation metrics: `SumMetric`, `MeanMetric`, `CatMetric`, `MinMetric`, `MaxMetric` ([#506](https://github.com/PyTorchLightning/metrics/pull/506))
- Added pairwise submodule with metrics ([#553](https://github.com/PyTorchLightning/metrics/pull/553))
  - `pairwise_cosine_similarity`
  - `pairwise_euclidean_distance`
  - `pairwise_linear_similarity`
  - `pairwise_manhatten_distance`

### Changed

- `AveragePrecision` will now as default output the `macro` average for multilabel and multiclass problems ([#477](https://github.com/PyTorchLightning/metrics/pull/477))
- `half`, `double`, `float` will no longer change the dtype of the metric states. Use `metric.set_dtype` instead ([#493](https://github.com/PyTorchLightning/metrics/pull/493))
- Renamed `AverageMeter` to `MeanMetric` ([#506](https://github.com/PyTorchLightning/metrics/pull/506))
- Changed `is_differentiable` from property to a constant attribute ([#551](https://github.com/PyTorchLightning/metrics/pull/551))
- `ROC` and `AUROC` will no longer throw an error when either the positive or negative class is missing. Instead return 0 score and give a warning

### Deprecated

- Deprecated  `functional.self_supervised.embedding_similarity` in favour of new pairwise submodule

### Removed

- Removed `dtype` property ([#493](https://github.com/PyTorchLightning/metrics/pull/493))

### Fixed

- Fixed bug in `F1` with `average='macro'` and `ignore_index!=None` ([#495](https://github.com/PyTorchLightning/metrics/pull/495))
- Fixed bug in `pit` by using the returned first result to initialize device and type ([#533](https://github.com/PyTorchLightning/metrics/pull/533))
- Fixed `SSIM` metric using too much memory ([#539](https://github.com/PyTorchLightning/metrics/pull/539))
- Fixed bug where `device` property was not properly update when metric was a child of a module (#542)

## [0.5.1] - 2021-08-30

### Added

- Added `device` and `dtype` properties ([#462](https://github.com/PyTorchLightning/metrics/pull/462))
- Added `TextTester` class for robustly testing text metrics ([#450](https://github.com/PyTorchLightning/metrics/pull/450))

### Changed

- Added support for float targets in `nDCG` metric ([#437](https://github.com/PyTorchLightning/metrics/pull/437))

### Removed

- Removed `rouge-score` as dependency for text package ([#443](https://github.com/PyTorchLightning/metrics/pull/443))
- Removed `jiwer` as dependency for text package ([#446](https://github.com/PyTorchLightning/metrics/pull/446))
- Removed `bert-score` as dependency for text package ([#473](https://github.com/PyTorchLightning/metrics/pull/473))

### Fixed

- Fixed ranking of samples in `SpearmanCorrCoef` metric ([#448](https://github.com/PyTorchLightning/metrics/pull/448))
- Fixed bug where compositional metrics where unable to sync because of type mismatch ([#454](https://github.com/PyTorchLightning/metrics/pull/454))
- Fixed metric hashing ([#478](https://github.com/PyTorchLightning/metrics/pull/478))
- Fixed `BootStrapper` metrics not working on GPU ([#462](https://github.com/PyTorchLightning/metrics/pull/462))
- Fixed the semantic ordering of kernel height and width in `SSIM` metric ([#474](https://github.com/PyTorchLightning/metrics/pull/474))


## [0.5.0] - 2021-08-09

### Added

- Added **Text-related (NLP) metrics**:
  - Word Error Rate (WER) ([#383](https://github.com/PyTorchLightning/metrics/pull/383))
  - ROUGE ([#399](https://github.com/PyTorchLightning/metrics/pull/399))
  - BERT score ([#424](https://github.com/PyTorchLightning/metrics/pull/424))
  - BLUE score ([#360](https://github.com/PyTorchLightning/metrics/pull/360))
- Added `MetricTracker` wrapper metric for keeping track of the same metric over multiple epochs ([#238](https://github.com/PyTorchLightning/metrics/pull/238))
- Added other metrics:
  - Symmetric Mean Absolute Percentage error (SMAPE) ([#375](https://github.com/PyTorchLightning/metrics/pull/375))
  - Calibration error ([#394](https://github.com/PyTorchLightning/metrics/pull/394))
  - Permutation Invariant Training (PIT) ([#384](https://github.com/PyTorchLightning/metrics/pull/384))
- Added support in `nDCG` metric for target with values larger than 1 ([#349](https://github.com/PyTorchLightning/metrics/pull/349))
- Added support for negative targets in `nDCG` metric ([#378](https://github.com/PyTorchLightning/metrics/pull/378))
- Added `None` as reduction option in `CosineSimilarity` metric ([#400](https://github.com/PyTorchLightning/metrics/pull/400))
- Allowed passing labels in (n_samples, n_classes) to `AveragePrecision` ([#386](https://github.com/PyTorchLightning/metrics/pull/386))

### Changed

- Moved `psnr` and `ssim` from `functional.regression.*` to `functional.image.*` ([#382](https://github.com/PyTorchLightning/metrics/pull/382))
- Moved `image_gradient` from `functional.image_gradients` to `functional.image.gradients` ([#381](https://github.com/PyTorchLightning/metrics/pull/381))
- Moved `R2Score` from `regression.r2score` to `regression.r2` ([#371](https://github.com/PyTorchLightning/metrics/pull/371))
- Pearson metric now only store 6 statistics instead of all predictions and targets ([#380](https://github.com/PyTorchLightning/metrics/pull/380))
- Use `torch.argmax` instead of `torch.topk` when `k=1` for better performance ([#419](https://github.com/PyTorchLightning/metrics/pull/419))
- Moved check for number of samples in R2 score to support single sample updating ([#426](https://github.com/PyTorchLightning/metrics/pull/426))

### Deprecated

- Rename `r2score` >> `r2_score` and `kldivergence` >> `kl_divergence` in `functional` ([#371](https://github.com/PyTorchLightning/metrics/pull/371))
- Moved `bleu_score` from `functional.nlp` to `functional.text.bleu` ([#360](https://github.com/PyTorchLightning/metrics/pull/360))

### Removed

- Removed restriction that `threshold` has to be in (0,1) range to support logit input (
    [#351](https://github.com/PyTorchLightning/metrics/pull/351)
    [#401](https://github.com/PyTorchLightning/metrics/pull/401))
- Removed restriction that `preds` could not be bigger than `num_classes` to support logit input ([#357](https://github.com/PyTorchLightning/metrics/pull/357))
- Removed module `regression.psnr` and `regression.ssim` ([#382](https://github.com/PyTorchLightning/metrics/pull/382)):
- Removed ([#379](https://github.com/PyTorchLightning/metrics/pull/379)):
    * function `functional.mean_relative_error`
    * `num_thresholds` argument in `BinnedPrecisionRecallCurve`

### Fixed

- Fixed bug where classification metrics with `average='macro'` would lead to wrong result if a class was missing ([#303](https://github.com/PyTorchLightning/metrics/pull/303))
- Fixed `weighted`, `multi-class` AUROC computation to allow for 0 observations of some class, as contribution to final AUROC is 0 ([#376](https://github.com/PyTorchLightning/metrics/pull/376))
- Fixed that `_forward_cache` and `_computed` attributes are also moved to the correct device if metric is moved ([#413](https://github.com/PyTorchLightning/metrics/pull/413))
- Fixed calculation in `IoU` metric when using `ignore_index` argument ([#328](https://github.com/PyTorchLightning/metrics/pull/328))


## [0.4.1] - 2021-07-05

### Changed

- Extend typing ([#330](https://github.com/PyTorchLightning/metrics/pull/330),
    [#332](https://github.com/PyTorchLightning/metrics/pull/332),
    [#333](https://github.com/PyTorchLightning/metrics/pull/333),
    [#335](https://github.com/PyTorchLightning/metrics/pull/335),
    [#314](https://github.com/PyTorchLightning/metrics/pull/314))

### Fixed

- Fixed DDP by `is_sync` logic to `Metric` ([#339](https://github.com/PyTorchLightning/metrics/pull/339))


## [0.4.0] - 2021-06-29

### Added

- Added **Image-related metrics**:
  - Fréchet inception distance (FID) ([#213](https://github.com/PyTorchLightning/metrics/pull/213))
  - Kernel Inception Distance (KID) ([#301](https://github.com/PyTorchLightning/metrics/pull/301))
  - Inception Score ([#299](https://github.com/PyTorchLightning/metrics/pull/299))
  - KL divergence ([#247](https://github.com/PyTorchLightning/metrics/pull/247))
- Added **Audio metrics**: SNR, SI_SDR, SI_SNR ([#292](https://github.com/PyTorchLightning/metrics/pull/292))
- Added other metrics:
  - Cosine Similarity ([#305](https://github.com/PyTorchLightning/metrics/pull/305))
  - Specificity ([#210](https://github.com/PyTorchLightning/metrics/pull/210))
  - Mean Absolute Percentage error (MAPE) ([#248](https://github.com/PyTorchLightning/metrics/pull/248))
- Added `add_metrics` method to `MetricCollection` for adding additional metrics after initialization ([#221](https://github.com/PyTorchLightning/metrics/pull/221))
- Added pre-gather reduction in the case of `dist_reduce_fx="cat"` to reduce communication cost ([#217](https://github.com/PyTorchLightning/metrics/pull/217))
- Added better error message for `AUROC` when `num_classes` is not provided for multiclass input ([#244](https://github.com/PyTorchLightning/metrics/pull/244))
- Added support for unnormalized scores (e.g. logits) in `Accuracy`, `Precision`, `Recall`, `FBeta`, `F1`, `StatScore`, `Hamming`, `ConfusionMatrix` metrics ([#200](https://github.com/PyTorchLightning/metrics/pull/200))
- Added `squared` argument to `MeanSquaredError` for computing `RMSE` ([#249](https://github.com/PyTorchLightning/metrics/pull/249))
- Added `is_differentiable` property to `ConfusionMatrix`, `F1`, `FBeta`, `Hamming`, `Hinge`, `IOU`, `MatthewsCorrcoef`, `Precision`, `Recall`, `PrecisionRecallCurve`, `ROC`, `StatScores` ([#253](https://github.com/PyTorchLightning/metrics/pull/253))
- Added `sync` and `sync_context` methods for manually controlling when metric states are synced ([#302](https://github.com/PyTorchLightning/metrics/pull/302))

### Changed

- Forward cache is reset when `reset` method is called ([#260](https://github.com/PyTorchLightning/metrics/pull/260))
- Improved per-class metric handling for imbalanced datasets for `precision`, `recall`, `precision_recall`, `fbeta`, `f1`, `accuracy`, and `specificity` ([#204](https://github.com/PyTorchLightning/metrics/pull/204))
- Decorated `torch.jit.unused` to `MetricCollection` forward ([#307](https://github.com/PyTorchLightning/metrics/pull/307))
- Renamed `thresholds` argument to binned metrics for manually controlling the thresholds ([#322](https://github.com/PyTorchLightning/metrics/pull/322))
- Extend typing ([#324](https://github.com/PyTorchLightning/metrics/pull/324),
    [#326](https://github.com/PyTorchLightning/metrics/pull/326),
    [#327](https://github.com/PyTorchLightning/metrics/pull/327))

### Deprecated

- Deprecated `functional.mean_relative_error`, use `functional.mean_absolute_percentage_error` ([#248](https://github.com/PyTorchLightning/metrics/pull/248))
- Deprecated `num_thresholds` argument in `BinnedPrecisionRecallCurve` ([#322](https://github.com/PyTorchLightning/metrics/pull/322))

### Removed

- Removed argument `is_multiclass` ([#319](https://github.com/PyTorchLightning/metrics/pull/319))

### Fixed

- AUC can also support more dimensional inputs when all but one dimension are of size 1 ([#242](https://github.com/PyTorchLightning/metrics/pull/242))
- Fixed `dtype` of modular metrics after reset has been called ([#243](https://github.com/PyTorchLightning/metrics/pull/243))
- Fixed calculation in `matthews_corrcoef` to correctly match formula ([#321](https://github.com/PyTorchLightning/metrics/pull/321))

## [0.3.2] - 2021-05-10

### Added

- Added `is_differentiable` property:
    * To `AUC`, `AUROC`, `CohenKappa` and `AveragePrecision` ([#178](https://github.com/PyTorchLightning/metrics/pull/178))
    * To `PearsonCorrCoef`, `SpearmanCorrcoef`, `R2Score` and `ExplainedVariance` ([#225](https://github.com/PyTorchLightning/metrics/pull/225))

### Changed

- `MetricCollection` should return metrics with prefix on `items()`, `keys()` ([#209](https://github.com/PyTorchLightning/metrics/pull/209))
- Calling `compute` before `update` will now give warning ([#164](https://github.com/PyTorchLightning/metrics/pull/164))

### Removed

- Removed `numpy` as direct dependency ([#212](https://github.com/PyTorchLightning/metrics/pull/212))

### Fixed

- Fixed auc calculation and add tests ([#197](https://github.com/PyTorchLightning/metrics/pull/197))
- Fixed loading persisted metric states using `load_state_dict()` ([#202](https://github.com/PyTorchLightning/metrics/pull/202))
- Fixed `PSNR` not working with `DDP` ([#214](https://github.com/PyTorchLightning/metrics/pull/214))
- Fixed metric calculation with unequal batch sizes ([#220](https://github.com/PyTorchLightning/metrics/pull/220))
- Fixed metric concatenation for list states for zero-dim input ([#229](https://github.com/PyTorchLightning/metrics/pull/229))
- Fixed numerical instability in `AUROC` metric for large input ([#230](https://github.com/PyTorchLightning/metrics/pull/230))

## [0.3.1] - 2021-04-21

- Cleaning remaining inconsistency and fix PL develop integration (
    [#191](https://github.com/PyTorchLightning/metrics/pull/191),
    [#192](https://github.com/PyTorchLightning/metrics/pull/192),
    [#193](https://github.com/PyTorchLightning/metrics/pull/193),
    [#194](https://github.com/PyTorchLightning/metrics/pull/194)
)


## [0.3.0] - 2021-04-20

### Added

- Added `BootStrapper` to easily calculate confidence intervals for metrics ([#101](https://github.com/PyTorchLightning/metrics/pull/101))
- Added Binned metrics  ([#128](https://github.com/PyTorchLightning/metrics/pull/128))
- Added metrics for Information Retrieval ([(PL^5032)](https://github.com/PyTorchLightning/pytorch-lightning/pull/5032)):
    * `RetrievalMAP` ([PL^5032](https://github.com/PyTorchLightning/pytorch-lightning/pull/5032))
    * `RetrievalMRR` ([#119](https://github.com/PyTorchLightning/metrics/pull/119))
    * `RetrievalPrecision` ([#139](https://github.com/PyTorchLightning/metrics/pull/139))
    * `RetrievalRecall` ([#146](https://github.com/PyTorchLightning/metrics/pull/146))
    * `RetrievalNormalizedDCG` ([#160](https://github.com/PyTorchLightning/metrics/pull/160))
    * `RetrievalFallOut` ([#161](https://github.com/PyTorchLightning/metrics/pull/161))
- Added other metrics:
    * `CohenKappa` ([#69](https://github.com/PyTorchLightning/metrics/pull/69))
    * `MatthewsCorrcoef` ([#98](https://github.com/PyTorchLightning/metrics/pull/98))
    * `PearsonCorrcoef` ([#157](https://github.com/PyTorchLightning/metrics/pull/157))
    * `SpearmanCorrcoef` ([#158](https://github.com/PyTorchLightning/metrics/pull/158))
    * `Hinge` ([#120](https://github.com/PyTorchLightning/metrics/pull/120))
- Added `average='micro'` as an option in AUROC for multilabel problems ([#110](https://github.com/PyTorchLightning/metrics/pull/110))
- Added multilabel support to `ROC` metric ([#114](https://github.com/PyTorchLightning/metrics/pull/114))
- Added testing for `half` precision ([#77](https://github.com/PyTorchLightning/metrics/pull/77),
    [#135](https://github.com/PyTorchLightning/metrics/pull/135)
)
- Added `AverageMeter` for ad-hoc averages of values ([#138](https://github.com/PyTorchLightning/metrics/pull/138))
- Added `prefix` argument to `MetricCollection` ([#70](https://github.com/PyTorchLightning/metrics/pull/70))
- Added `__getitem__` as metric arithmetic operation ([#142](https://github.com/PyTorchLightning/metrics/pull/142))
- Added property `is_differentiable` to metrics and test for differentiability ([#154](https://github.com/PyTorchLightning/metrics/pull/154))
- Added support for `average`, `ignore_index` and `mdmc_average` in `Accuracy` metric ([#166](https://github.com/PyTorchLightning/metrics/pull/166))
- Added `postfix` arg to `MetricCollection` ([#188](https://github.com/PyTorchLightning/metrics/pull/188))

### Changed

- Changed `ExplainedVariance` from storing all preds/targets to tracking 5 statistics ([#68](https://github.com/PyTorchLightning/metrics/pull/68))
- Changed behaviour of `confusionmatrix` for multilabel data to better match `multilabel_confusion_matrix` from sklearn ([#134](https://github.com/PyTorchLightning/metrics/pull/134))
- Updated FBeta arguments ([#111](https://github.com/PyTorchLightning/metrics/pull/111))
- Changed `reset` method to use `detach.clone()` instead of `deepcopy` when resetting to default ([#163](https://github.com/PyTorchLightning/metrics/pull/163))
- Metrics passed as dict to `MetricCollection` will now always be in deterministic order ([#173](https://github.com/PyTorchLightning/metrics/pull/173))
- Allowed `MetricCollection` pass metrics as arguments ([#176](https://github.com/PyTorchLightning/metrics/pull/176))

### Deprecated

- Rename argument `is_multiclass` -> `multiclass` ([#162](https://github.com/PyTorchLightning/metrics/pull/162))

### Removed

- Prune remaining deprecated ([#92](https://github.com/PyTorchLightning/metrics/pull/92))

### Fixed

- Fixed when `_stable_1d_sort` to work when `n>=N` ([PL^6177](https://github.com/PyTorchLightning/pytorch-lightning/pull/6177))
- Fixed `_computed` attribute not being correctly reset ([#147](https://github.com/PyTorchLightning/metrics/pull/147))
- Fixed to Blau score ([#165](https://github.com/PyTorchLightning/metrics/pull/165))
- Fixed backwards compatibility for logging with older version of pytorch-lightning ([#182](https://github.com/PyTorchLightning/metrics/pull/182))


## [0.2.0] - 2021-03-12

### Changed

- Decoupled PL dependency ([#13](https://github.com/PyTorchLightning/metrics/pull/13))
- Refactored functional - mimic the module-like structure: classification, regression, etc. ([#16](https://github.com/PyTorchLightning/metrics/pull/16))
- Refactored utilities -  split to topics/submodules ([#14](https://github.com/PyTorchLightning/metrics/pull/14))
- Refactored `MetricCollection` ([#19](https://github.com/PyTorchLightning/metrics/pull/19))

### Removed

- Removed deprecated metrics from PL base ([#12](https://github.com/PyTorchLightning/metrics/pull/12),
    [#15](https://github.com/PyTorchLightning/metrics/pull/15))



## [0.1.0] - 2021-02-22

- Added `Accuracy` metric now generalizes to Top-k accuracy for (multi-dimensional) multi-class inputs using the `top_k` parameter ([PL^4838](https://github.com/PyTorchLightning/pytorch-lightning/pull/4838))
- Added `Accuracy` metric now enables the computation of subset accuracy for multi-label or multi-dimensional multi-class inputs with the `subset_accuracy` parameter ([PL^4838](https://github.com/PyTorchLightning/pytorch-lightning/pull/4838))
- Added `HammingDistance` metric to compute the hamming distance (loss) ([PL^4838](https://github.com/PyTorchLightning/pytorch-lightning/pull/4838))
- Added `StatScores` metric to compute the number of true positives, false positives, true negatives and false negatives ([PL^4839](https://github.com/PyTorchLightning/pytorch-lightning/pull/4839))
- Added `R2Score` metric ([PL^5241](https://github.com/PyTorchLightning/pytorch-lightning/pull/5241))
- Added `MetricCollection` ([PL^4318](https://github.com/PyTorchLightning/pytorch-lightning/pull/4318))
- Added `.clone()` method to metrics ([PL^4318](https://github.com/PyTorchLightning/pytorch-lightning/pull/4318))
- Added `IoU` class interface ([PL^4704](https://github.com/PyTorchLightning/pytorch-lightning/pull/4704))
- The `Recall` and `Precision` metrics (and their functional counterparts `recall` and `precision`) can now be generalized to Recall@K and Precision@K with the use of `top_k` parameter ([PL^4842](https://github.com/PyTorchLightning/pytorch-lightning/pull/4842))
- Added compositional metrics ([PL^5464](https://github.com/PyTorchLightning/pytorch-lightning/pull/5464))
- Added AUC/AUROC class interface ([PL^5479](https://github.com/PyTorchLightning/pytorch-lightning/pull/5479))
- Added `QuantizationAwareTraining` callback ([PL^5706](https://github.com/PyTorchLightning/pytorch-lightning/pull/5706))
- Added `ConfusionMatrix` class interface ([PL^4348](https://github.com/PyTorchLightning/pytorch-lightning/pull/4348))
- Added multiclass AUROC metric ([PL^4236](https://github.com/PyTorchLightning/pytorch-lightning/pull/4236))
- Added `PrecisionRecallCurve, ROC, AveragePrecision` class metric ([PL^4549](https://github.com/PyTorchLightning/pytorch-lightning/pull/4549))
- Classification metrics overhaul ([PL^4837](https://github.com/PyTorchLightning/pytorch-lightning/pull/4837))
- Added `F1` class metric ([PL^4656](https://github.com/PyTorchLightning/pytorch-lightning/pull/4656))
- Added metrics aggregation in Horovod and fixed early stopping ([PL^3775](https://github.com/PyTorchLightning/pytorch-lightning/pull/3775))
- Added `persistent(mode)` method to metrics, to enable and disable metric states being added to `state_dict` ([PL^4482](https://github.com/PyTorchLightning/pytorch-lightning/pull/4482))
- Added unification of regression metrics ([PL^4166](https://github.com/PyTorchLightning/pytorch-lightning/pull/4166))
- Added persistent flag to `Metric.add_state` ([PL^4195](https://github.com/PyTorchLightning/pytorch-lightning/pull/4195))
- Added classification metrics ([PL^4043](https://github.com/PyTorchLightning/pytorch-lightning/pull/4043))
- Added new Metrics API. ([PL^3868](https://github.com/PyTorchLightning/pytorch-lightning/pull/3868), [PL^3921](https://github.com/PyTorchLightning/pytorch-lightning/pull/3921))
- Added EMB similarity ([PL^3349](https://github.com/PyTorchLightning/pytorch-lightning/pull/3349))
- Added SSIM metrics ([PL^2671](https://github.com/PyTorchLightning/pytorch-lightning/pull/2671))
- Added BLEU metrics ([PL^2535](https://github.com/PyTorchLightning/pytorch-lightning/pull/2535))
<div align="center">

<img src="docs/source/_static/images/logo.png" width="400px">

**Machine learning metrics for distributed, scalable PyTorch applications.**

______________________________________________________________________

<p align="center">
  <a href="#what-is-torchmetrics">What is Torchmetrics</a> •
  <a href="#implementing-your-own-module-metric">Implementing a metric</a> •
  <a href="#build-in-metrics">Built-in metrics</a> •
  <a href="https://torchmetrics.readthedocs.io/en/stable/">Docs</a> •
  <a href="#community">Community</a> •
  <a href="#license">License</a>
</p>

______________________________________________________________________

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/torchmetrics)](https://pypi.org/project/torchmetrics/)
[![PyPI Status](https://badge.fury.io/py/torchmetrics.svg)](https://badge.fury.io/py/torchmetrics)
[![PyPI Status](https://pepy.tech/badge/torchmetrics)](https://pepy.tech/project/torchmetrics)
[![Conda](https://img.shields.io/conda/v/conda-forge/torchmetrics?label=conda&color=success)](https://anaconda.org/conda-forge/torchmetrics)
![Conda](https://img.shields.io/conda/dn/conda-forge/torchmetrics)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/PytorchLightning/metrics/blob/master/LICENSE)

[![CI testing - base](https://github.com/PyTorchLightning/metrics/actions/workflows/ci_test-base.yml/badge.svg?branch=master&event=push)](https://github.com/PyTorchLightning/metrics/actions/workflows/ci_test-base.yml)
[![PyTorch & Conda](https://github.com/PyTorchLightning/metrics/actions/workflows/ci_test-conda.yml/badge.svg?branch=master&event=push)](https://github.com/PyTorchLightning/metrics/actions/workflows/ci_test-conda.yml)
[![Build Status](https://dev.azure.com/PytorchLightning/Metrics/_apis/build/status/PyTorchLightning.metrics?branchName=master)](https://dev.azure.com/PytorchLightning/Metrics/_build/latest?definitionId=3&branchName=master)
[![codecov](https://codecov.io/gh/PyTorchLightning/metrics/branch/master/graph/badge.svg?token=NER6LPI3HS)](https://codecov.io/gh/PyTorchLightning/metrics)

[![Slack](https://img.shields.io/badge/slack-chat-green.svg?logo=slack)](https://join.slack.com/t/pytorch-lightning/shared_invite/zt-12iz3cds1-uyyyBYJLiaL2bqVmMN7n~A)
[![Documentation Status](https://readthedocs.org/projects/torchmetrics/badge/?version=latest)](https://torchmetrics.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5844769.svg)](https://doi.org/10.5281/zenodo.5844769)
[![JOSS status](https://joss.theoj.org/papers/561d9bb59b400158bc8204e2639dca43/status.svg)](https://joss.theoj.org/papers/561d9bb59b400158bc8204e2639dca43)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/PyTorchLightning/metrics/master.svg)](https://results.pre-commit.ci/latest/github/PyTorchLightning/metrics/master)

______________________________________________________________________

</div>

## Installation

Simple installation from PyPI

```bash
pip install torchmetrics
```

<details>
  <summary>Other installations</summary>

Install using conda

```bash
conda install -c conda-forge torchmetrics
```

Pip from source

```bash
# with git
pip install git+https://github.com/PytorchLightning/metrics.git@release/latest
```

Pip from archive

```bash
pip install https://github.com/PyTorchLightning/metrics/archive/refs/heads/release/latest.zip
```

Extra dependencies for specialized metrics:

```bash
pip install torchmetrics[audio]
pip install torchmetrics[image]
pip install torchmetrics[text]
pip install torchmetrics[all]  # install all of the above
```

Install latest developer version

```bash
pip install https://github.com/PyTorchLightning/metrics/archive/master.zip
```

</details>

______________________________________________________________________

## What is Torchmetrics

TorchMetrics is a collection of 80+ PyTorch metrics implementations and an easy-to-use API to create custom metrics. It offers:

- A standardized interface to increase reproducibility
- Reduces boilerplate
- Automatic accumulation over batches
- Metrics optimized for distributed-training
- Automatic synchronization between multiple devices

You can use TorchMetrics with any PyTorch model or with [PyTorch Lightning](https://pytorch-lightning.readthedocs.io/en/stable/) to enjoy additional features such as:

- Module metrics are automatically placed on the correct device.
- Native support for logging metrics in Lightning to reduce even more boilerplate.

## Using TorchMetrics

### Module metrics

The [module-based metrics](https://pytorchlightning.github.io/metrics/references/modules.html) contain internal metric states (similar to the parameters of the PyTorch module) that automate accumulation and synchronization across devices!

- Automatic accumulation over multiple batches
- Automatic synchronization between multiple devices
- Metric arithmetic

**This can be run on CPU, single GPU or multi-GPUs!**

For the single GPU/CPU case:

```python
import torch

# import our library
import torchmetrics

# initialize metric
metric = torchmetrics.Accuracy()

# move the metric to device you want computations to take place
device = "cuda" if torch.cuda.is_available() else "cpu"
metric.to(device)

n_batches = 10
for i in range(n_batches):
    # simulate a classification problem
    preds = torch.randn(10, 5).softmax(dim=-1).to(device)
    target = torch.randint(5, (10,)).to(device)

    # metric on current batch
    acc = metric(preds, target)
    print(f"Accuracy on batch {i}: {acc}")

# metric on all batches using custom accumulation
acc = metric.compute()
print(f"Accuracy on all data: {acc}")
```

Module metric usage remains the same when using multiple GPUs or multiple nodes.

<details>
  <summary>Example using DDP</summary>

<!--phmdoctest-mark.skip-->

```python
import os
import torch
import torch.distributed as dist
import torch.multiprocessing as mp
from torch import nn
from torch.nn.parallel import DistributedDataParallel as DDP
import torchmetrics


def metric_ddp(rank, world_size):
    os.environ["MASTER_ADDR"] = "localhost"
    os.environ["MASTER_PORT"] = "12355"

    # create default process group
    dist.init_process_group("gloo", rank=rank, world_size=world_size)

    # initialize model
    metric = torchmetrics.Accuracy()

    # define a model and append your metric to it
    # this allows metric states to be placed on correct accelerators when
    # .to(device) is called on the model
    model = nn.Linear(10, 10)
    model.metric = metric
    model = model.to(rank)

    # initialize DDP
    model = DDP(model, device_ids=[rank])

    n_epochs = 5
    # this shows iteration over multiple training epochs
    for n in range(n_epochs):

        # this will be replaced by a DataLoader with a DistributedSampler
        n_batches = 10
        for i in range(n_batches):
            # simulate a classification problem
            preds = torch.randn(10, 5).softmax(dim=-1)
            target = torch.randint(5, (10,))

            # metric on current batch
            acc = metric(preds, target)
            if rank == 0:  # print only for rank 0
                print(f"Accuracy on batch {i}: {acc}")

        # metric on all batches and all accelerators using custom accumulation
        # accuracy is same across both accelerators
        acc = metric.compute()
        print(f"Accuracy on all data: {acc}, accelerator rank: {rank}")

        # Reseting internal state such that metric ready for new data
        metric.reset()

    # cleanup
    dist.destroy_process_group()


if __name__ == "__main__":
    world_size = 2  # number of gpus to parallize over
    mp.spawn(metric_ddp, args=(world_size,), nprocs=world_size, join=True)
```

</details>

### Implementing your own Module metric

Implementing your own metric is as easy as subclassing an [`torch.nn.Module`](https://pytorch.org/docs/stable/generated/torch.nn.Module.html). Simply, subclass `torchmetrics.Metric`
and implement the following methods:

```python
import torch
from torchmetrics import Metric


class MyAccuracy(Metric):
    def __init__(self, dist_sync_on_step=False):
        # call `self.add_state`for every internal state that is needed for the metrics computations
        # dist_reduce_fx indicates the function that should be used to reduce
        # state from multiple processes
        super().__init__(dist_sync_on_step=dist_sync_on_step)

        self.add_state("correct", default=torch.tensor(0), dist_reduce_fx="sum")
        self.add_state("total", default=torch.tensor(0), dist_reduce_fx="sum")

    def update(self, preds: torch.Tensor, target: torch.Tensor):
        # update metric states
        preds, target = self._input_format(preds, target)
        assert preds.shape == target.shape

        self.correct += torch.sum(preds == target)
        self.total += target.numel()

    def compute(self):
        # compute final result
        return self.correct.float() / self.total
```

### Functional metrics

Similar to [`torch.nn`](https://pytorch.org/docs/stable/nn.html), most metrics have both a [module-based](https://torchmetrics.readthedocs.io/en/latest/references/modules.html) and a [functional](https://torchmetrics.readthedocs.io/en/latest/references/functional.html) version.
The functional versions are simple python functions that as input take [torch.tensors](https://pytorch.org/docs/stable/tensors.html) and return the corresponding metric as a [torch.tensor](https://pytorch.org/docs/stable/tensors.html).

```python
import torch

# import our library
import torchmetrics

# simulate a classification problem
preds = torch.randn(10, 5).softmax(dim=-1)
target = torch.randint(5, (10,))

acc = torchmetrics.functional.accuracy(preds, target)
```

### Covered domains and example metrics

We currently have implemented metrics within the following domains:

- Audio (
  [ScaleInvariantSignalDistortionRatio](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#ScaleInvariantSignalDistortionRatio),
  [ScaleInvariantSignalNoiseRatio](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#ScaleInvariantSignalNoiseRatio),
  [SignalNoiseRatio](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#SignalNoiseRatio)
  and [few more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#audio-metrics)
  )
- Classification (
  [Accuracy](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#accuracy),
  [F1Score](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#f1score),
  [AUROC](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#auroc)
  and [many more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#classification-metrics)
  )
- Detection (
  [MeanAveragePrecision](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#meanaverageprecision))
- Information Retrieval (
  [RetrievalMAP](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#retrievalmap),
  [RetrievalMRR](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#retrievalmrr),
  [RetrievalNormalizedDCG](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#retrievalnormalizeddcg)
  and [few more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#retrieval)
  )
- Image (
  [FrechetInceptionDistance](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#FrechetInceptionDistance),
  [KernelInceptionDistance](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#KernelInceptionDistance),
  [StructuralSimilarityIndexMeasure](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#StructuralSimilarityIndexMeasure)
  and [many more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#image-metrics)
  )
- Regression (
  [ExplainedVariance](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#explainedvariance),
  [PearsonCorrCoef](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#pearsoncorrcoef),
  [R2Score](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#r2score)
  and [many more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#regression-metrics)
  )
- Text (
  [BleuScore](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#bleuscore),
  [RougeScore](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#rougescore),
  [WordErrorRate](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#WordErrorRate)
  and [many more](https://torchmetrics.readthedocs.io/en/latest/references/modules.html#text)
  )

In total torchmetrics contains 80+ metrics!

## Contribute!

The lightning + torchmetric team is hard at work adding even more metrics.
But we're looking for incredible contributors like you to submit new metrics
and improve existing ones!

Join our [Slack](https://join.slack.com/t/pytorch-lightning/shared_invite/zt-12iz3cds1-uyyyBYJLiaL2bqVmMN7n~A)
to get help becoming a contributor!

## Community

For help or questions, join our huge community on [Slack](https://join.slack.com/t/pytorch-lightning/shared_invite/zt-12iz3cds1-uyyyBYJLiaL2bqVmMN7n~A)!

## Citation

We’re excited to continue the strong legacy of open source software and have been inspired
over the years by Caffe, Theano, Keras, PyTorch, torchbearer, ignite, sklearn and fast.ai.

If you want to cite this framework feel free to use GitHub's built-in citation option to generate a bibtex or APA-Style citation based on [this file](https://github.com/PyTorchLightning/metrics/blob/master/CITATION.cff)(but only if you loved it 😊).

## License

Please observe the Apache 2.0 license that is listed in this repository.
In addition, the Lightning framework is Patent Pending.
---
title: TorchMetrics - Measuring Reproducibility in PyTorch
tags:
  - python
  - deep learning
  - pytorch
authors:
  - name: Nicki Skafte Detlefsen
    affiliation: '1,2'
    orcid: 0000-0002-8133-682X
  - name: Jiri Borovec
    affiliation: '1'
    orcid: 0000-0001-7437-824X
  - name: Justus Schock
    affiliation: '1,3'
    orcid: 0000-0003-0512-3053
  - name: Ananya Harsh Jha
    affiliation: '1'
  - name: Teddy Koker
    affiliation: '1'
  - name: Luca Di Liello
    affiliation: '4'
  - name: Daniel Stancl
    affiliation: '5'
  - name: Changsheng Quan
    affiliation: '6'
  - name: Maxim Grechkin
    affiliation: '7'
  - name: William Falcon
    affiliation: '1,8'
affiliations:
  - name: Grid AI Labs
    index: 1
  - name: Technical University of Denmark
    index: 2
  - name: University Hospital Düsseldorf
    index: 3
  - name: University of Trento
    index: 4
  - name: Charles University
    index: 5
  - name: Zhejiang University
    index: 6
  - name: Independent Researcher
    index: 7
  - name: New York University
    index: 8
date: 08 Dec 2021
bibliography: paper.bib
---

# Summary

A main problem with reproducing machine learning publications is the variance of metric implementations across papers. A lack of standardization leads to different behavior in mechanisms such as checkpointing, learning rate schedulers or early stopping, that will influence the reported results. For example, a complex metric such as Fréchet inception distance (FID) for synthetic image quality evaluation [@fid] will differ based on the specific interpolation method used.

There have been a few attempts at tackling the reproducibility issues. Papers With Code [@papers_with_code] links research code with its corresponding paper. Similarly, arXiv [@arxiv] recently added a code and data section that links both official and community code to papers. However, these methods rely on the paper code to be made publicly accessible which is not always possible. Our approach is to provide the de-facto reference implementation for metrics. This approach enables proprietary work to still be comparable as long as they’ve used our reference implementations.

We introduce TorchMetrics, a general-purpose metrics package that covers a wide variety of tasks and domains used in the machine learning community. TorchMetrics provides standard classification and regression metrics; and domain-specific metrics for audio, computer vision, natural language processing, and information retrieval. Our process for adding a new metric is as follows, first we integrate a well-tested and established third-party library. Once we’ve verified the implementations and written tests for them, we re-implement them in native PyTorch [@pytorch] to enable hardware acceleration and remove any bottlenecks in inter-device transfer.

# Statement of need

Currently, there is no standard, widely-adopted metrics library for native PyTorch. Some native PyTorch libraries support domain-specific metrics such as Transformers [@transformers] for calculating NLP-specific metrics. However, no library exists that covers multiple domains. PyTorch users, therefore, often rely on non-PyTorch packages such as Scikit-learn [@scikit_learn] for computing even simple metrics such as accuracy, F1, or AUROC metrics.

However, while Scikit-learn is considered the gold standard for computing metrics in regression and classification, it relies on the core assumption that all predictions and targets are available simultaneously. This contradicts the typical workflow in a modern deep learning training/evaluation loop where data comes in batches. Therefore, the metric needs to be calculated in an online fashion. It is important to note that, in general, it is not possible to calculate a global metric as its average or sum of the metric calculated per batch.

TorchMetrics solves this problem by introducing stateful metrics that can calculate metric values on a stream of data alongside the classical functional and stateless metrics provided by other packages like Scikit-learn. We do this with an effortless `update` and `compute` interface, well known from packages such as Keras [@keras]. The `update` function takes in a batch of predictions and targets and updates the internal state. For example, for a metric such as accuracy, the internal states are simply the number of correctly classified samples and the total observed number of samples. When all batches have been passed to the `update` method, the `compute` method can get the accumulated accuracy over all the batches. In addition to `update` and `compute`, each metric also has a `forward` method (as any other `torch.nn.Module`) that can be used to both get the metric on the current batch of data and accumulate global state. This enables the user to get fine-grained info about the metric on the individual batch and the global metric of how well their model is doing.

```python
# Minimal example showcasing the TorchMetrics interface
import torch
from torch import tensor, Tensor
# base class all modular metrics inherit from
from torchmetrics import Metric

class Accuracy(Metric):
    def __init__(self):
        super().__init__()
        # `self.add_state` defines the states of the metric
        #  that should be accumulated and will automatically
        #  be synchronized between devices
        self.add_state("correct", default=tensor(0), dist_reduce_fx="sum")
        self.add_state("total", default=tensor(0), dist_reduce_fx="sum")

    def update(self, preds: Tensor, target: Tensor) -> None:
        # update takes `preds` and `target` and accumulate the current
        # stream of data into the global states for later
        self.correct += torch.sum(preds == target)
        self.total += target.numel()

    def compute(self) -> Tensor:
        # compute takes the accumulated states
        # and returns the final metric value
        return self.correct / self.total
```

Another core feature of TorchMetrics is its ability to scale to multiple devices seamlessly. Modern deep learning models are often trained on hundreds of devices such as GPUs or TPUs (see @large_example1; @large_example2 for examples). This scale introduces the need to synchronize metrics across machines to get the correct value during training and evaluation. In distributed environments, TorchMetrics automatically accumulates across devices before reporting the calculated metric to the user.

In addition to stateful metrics (called modular metrics in TorchMetrics), we also support a functional interface that works similar to Scikit-learn. This interface provides simple Python functions that take as input PyTorch Tensors and return the corresponding metric as a PyTorch Tensor. These can be used when metrics are evaluated on single devices, and no accumulation is needed, making them very fast to compute.

TorchMetrics exhibits high test coverage on the various configurations, including all three major OS platforms (Linux, macOS, and Windows), and various Python, CUDA, and PyTorch versions. We test both minimum and latest package requirements for all combinations of OS and Python versions and include additional tests for each PyTorch version from 1.3 up to future development versions. On every pull request and merge to master, we run a full test suite. All standard tests run on CPU. In addition, we run all tests on a multi-GPU setting which reflects realistic Deep Learning workloads. For usability, we have auto-generated HTML documentation (hosted at [readthedocs](https://torchmetrics.readthedocs.io/en/stable/)) from the source code which updates in real-time with new merged pull requests.

TorchMetrics is released under the Apache 2.0 license. The source code is available at https://github.com/PyTorchLightning/metrics.

# Acknowledgement

The TorchMetrics team thanks Thomas Chaton, Ethan Harris, Carlos Mocholí, Sean Narenthiran, Adrian Wälchli, and Ananth Subramaniam for contributing ideas, participating in discussions on API design, and completing Pull Request reviews. We also thank all of our open-source contributors for reporting and resolving issues with this package. We are grateful to the PyTorch Lightning team for their ongoing and dedicated support of this project, and Grid.ai for providing computing resources and cloud credits needed to run our Continuos Integrations.

# References
## What does this PR do?

Fixes #\<issue_number>

## Before submitting

- [ ] Was this **discussed/approved** via a Github issue? (no need for typos and docs improvements)
- [ ] Did you read the [contributor guideline](https://github.com/PyTorchLightning/metrics/blob/master/.github/CONTRIBUTING.md), Pull Request section?
- [ ] Did you make sure to **update the docs**?
- [ ] Did you write any new **necessary tests**?

## PR review

Anyone in the community is free to review the PR once the tests have passed.
If we didn't discuss your PR in Github issues there's a high chance it will not be merged.

## Did you have fun?

Make sure you had fun coding 🙃
# Contributing

Welcome to the Torchmetrics community! We're building largest collection of native pytorch metrics, with the
goal of reducing boilerplate and increasing reproducibility.

## Contribution Types

We are always looking for help implementing new features or fixing bugs.

### Bug Fixes:

1. If you find a bug please submit a github issue.

   - Make sure the title explains the issue.
   - Describe your setup, what you are trying to do, expected vs. actual behaviour. Please add configs and code samples.
   - Add details on how to reproduce the issue - a minimal test case is always best, colab is also great.
     Note, that the sample code shall be minimal and if needed with publicly available data.

1. Try to fix it or recommend a solution. We highly recommend to use test-driven approach:

   - Convert your minimal code example to a unit/integration test with assert on expected results.
   - Start by debugging the issue... You can run just this particular test in your IDE and draft a fix.
   - Verify that your test case fails on the master branch and only passes with the fix applied.

1. Submit a PR!

_**Note**, even if you do not find the solution, sending a PR with a test covering the issue is a valid contribution and we can
help you or finish it with you :\]_

### New Features:

1. Submit a github issue - describe what is the motivation of such feature (adding the use case or an example is helpful).

1. Let's discuss to determine the feature scope.

1. Submit a PR! We recommend test driven approach to adding new features as well:

   - Write a test for the functionality you want to add.
   - Write the functional code until the test passes.

1. Add/update the relevant tests!

- [This PR](https://github.com/PyTorchLightning/metrics/pull/98) is a good example for adding a new metric

### Test cases:

Want to keep Torchmetrics healthy? Love seeing those green tests? So do we! How to we keep it that way?
We write tests! We value tests contribution even more than new features. One of the core values of torchmetrics
is that our users can trust our metric implementation. We can only guarantee this if our metrics are well tested.

______________________________________________________________________

## Guidelines

### Developments scripts

To build the documentation locally, simply execute the following commands from project root (only for Unix):

- `make clean` cleans repo from temp/generated files
- `make docs` builds documentation under _docs/build/html_
- `make test` runs all project's tests with coverage

### Original code

All added or edited code shall be the own original work of the particular contributor.
If you use some third-party implementation, all such blocks/functions/modules shall be properly referred and if
possible also agreed by code's author. For example - `This code is inspired from http://...`.
In case you adding new dependencies, make sure that they are compatible with the actual Torchmetrics license
(ie. dependencies should be _at least_ as permissive as the Torchmetrics license).

### Coding Style

1. Use f-strings for output formation (except logging when we stay with lazy `logging.info("Hello %s!", name)`.
1. You can use `pre-commit` to make sure your code style is correct.

### Documentation

We are using Sphinx with Napoleon extension.
Moreover, we set Google style to follow with type convention.

- [Napoleon formatting with Google style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
- [ReStructured Text (reST)](https://docs.pylonsproject.org/projects/docs-style-guide/)
- [Paragraph-level markup](https://www.sphinx-doc.org/en/1.5/markup/para.html)

See following short example of a sample function taking one position string and optional

```python
from typing import Optional


def my_func(param_a: int, param_b: Optional[float] = None) -> str:
    """Sample function.

    Args:
        param_a: first parameter
        param_b: second parameter

    Return:
        sum of both numbers

    Example:
        Sample doctest example...
        >>> my_func(1, 2)
        3

    .. note:: If you want to add something.
    """
    p = param_b if param_b else 0
    return str(param_a + p)
```

When updating the docs make sure to build them first locally and visually inspect the html files (in the browser) for
formatting errors. In certain cases, a missing blank line or a wrong indent can lead to a broken layout.
Run these commands

```bash
make docs
```

and open `docs/build/html/index.html` in your browser.

Notes:

- You need to have LaTeX installed for rendering math equations. You can for example install TeXLive by doing one of the following:
  - on Ubuntu (Linux) run `apt-get install texlive` or otherwise follow the instructions on the TeXLive website
  - use the [RTD docker image](https://hub.docker.com/r/readthedocs/build)
- with PL used class meta you need to use python 3.7 or higher

When you send a PR the continuous integration will run tests and build the docs.

### Testing

**Local:** Testing your work locally will help you speed up the process since it allows you to focus on particular (failing) test-cases.
To setup a local development environment, install both local and test dependencies:

```bash
python -m pip install -r requirements/test.txt
python -m pip install pre-commit
```

You can run the full test-case in your terminal via this make script:

```bash
make test
# or natively
python -m pytest torchmetrics tests
```

Note: if your computer does not have multi-GPU nor TPU these tests are skipped.

**GitHub Actions:** For convenience, you can also use your own GHActions building which will be triggered with each commit.
This is useful if you do not test against all required dependency versions.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual attention or
  advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or electronic
  address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at waf2107@columbia.edu. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

[homepage]: https://www.contributor-covenant.org
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''
---

## 🚀 Feature

<!-- A clear and concise description of the feature proposal -->

### Motivation

<!-- Please outline the motivation for the proposal. Is your feature request related to a problem? e.g., I'm always frustrated when [...]. If this is related to another GitHub issue, please link here too -->

### Pitch

<!-- A clear and concise description of what you want to happen. -->

### Alternatives

<!-- A clear and concise description of any alternative solutions or features you've considered, if any. -->

### Additional context

<!-- Add any other context or screenshots about the feature request here. -->
---
name: Typos and doc fixes
about: Typos and doc fixes
title: ''
labels: documentation
assignees: ''
---

## 📚 Documentation

For typos and doc fixes, please go ahead and:

1. Create an issue.
1. Fix the typo.
1. Submit a PR.

Thanks!
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug / fix, help wanted
assignees: ''
---

## 🐛 Bug

<!-- A clear and concise description of what the bug is. -->

### To Reproduce

Steps to reproduce the behavior...

<!-- If you have a code sample, error messages, stack traces, please provide it here as well -->

#### Code sample

<!-- Ideally attach a minimal code sample to reproduce the decried issue.
Minimal means having the shortest code but still preserving the bug. -->

### Expected behavior

<!-- A clear and concise description of what you expected to happen. -->

### Environment

- OS (e.g., Linux):
- Python & PyTorch Version (e.g., 1.0):
- How you installed PyTorch (`conda`, `pip`, build command if you used source):
- Any other relevant information:

### Additional context

<!-- Add any other context about the problem here. -->

.. _scikit-learn's implementation of SMAPE: https://github.com/scikit-learn/scikit-learn/blob/15a949460/sklearn/metrics/_regression.py#L197
.. _scikit-learn's implementation of MAPE: https://github.com/scikit-learn/scikit-learn/blob/15a949460/sklearn/metrics/_regression.py#L197
.. _Mean Average Precision: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Mean_average_precision
.. _Fall-out: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Fall-out
.. _Normalized Discounted Cumulative Gain: https://en.wikipedia.org/wiki/Discounted_cumulative_gain
.. _IR Precision: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Precision
.. _IR R-Precision: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#R-precision
.. _IR Recall: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Recall
.. _Accuracy: https://en.wikipedia.org/wiki/Accuracy_and_precision
.. _SMAPE: https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error
.. _SNR: https://en.wikipedia.org/wiki/Signal-to-noise_ratio
.. _ROC AUC: https://en.wikipedia.org/wiki/Receiver_operating_characteristic#Further_interpretations
.. _Cohen's kappa score: https://en.wikipedia.org/wiki/Cohen%27s_kappa
.. _scikit-learn's implementation of confusion matrix: https://scikit-learn.org/stable/modules/model_evaluation.html#confusion-matrix
.. _confusion matrix gets calculated per label: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.multilabel_confusion_matrix.html
.. _F-score: https://en.wikipedia.org/wiki/F-score
.. _Hamming distance: https://en.wikipedia.org/wiki/Hamming_distance
.. _Hinge loss: https://en.wikipedia.org/wiki/Hinge_loss
.. _Jaccard index: https://en.wikipedia.org/wiki/Jaccard_index
.. _KL divergence: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
.. _Matthews correlation coefficient: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
.. _Precision: https://en.wikipedia.org/wiki/Precision_and_recall
.. _Recall: https://en.wikipedia.org/wiki/Precision_and_recall
.. _Specificity: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
.. _Type I and Type II errors: https://en.wikipedia.org/wiki/Type_I_and_type_II_errors
.. _confusion matrix: https://en.wikipedia.org/wiki/Confusion_matrix#Table_of_confusion
.. _sklearn averaging methods: https://scikit-learn.org/stable/modules/model_evaluation.html#multiclass-and-multilabel-classification
.. _Cosine Similarity: https://en.wikipedia.org/wiki/Cosine_similarity
.. _spearmans rank correlation coefficient: https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
.. _WordErrorRate: https://en.wikipedia.org/wiki/Word_error_rate
.. _FID: https://en.wikipedia.org/wiki/Fr%C3%A9chet_inception_distance
.. _mean-squared-error: https://en.wikipedia.org/wiki/Mean_squared_error
.. _SSIM: https://en.wikipedia.org/wiki/Structural_similarity
.. _explained variance: https://en.wikipedia.org/wiki/Explained_variation
.. _IR Average precision: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Average_precision
.. _IR Fall-out: https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Fall-out
.. _MAPE implementation returns: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mean_absolute_percentage_error.html
.. _mean squared logarithmic error: https://scikit-learn.org/stable/modules/model_evaluation.html#mean-squared-log-error
.. _LPIPS: https://arxiv.org/abs/1801.03924
.. _Mean-Average-Precision (mAP) and Mean-Average-Recall (mAR): https://jonathan-hui.medium.com/map-mean-average-precision-for-object-detection-45c121a31173
.. _Tweedie Deviance Score: https://en.wikipedia.org/wiki/Tweedie_distribution#The_Tweedie_deviance
.. _Permutation Invariant Training of Deep Models: https://ieeexplore.ieee.org/document/7952154
.. _Computes the Top-label Calibration Error: https://arxiv.org/pdf/1909.10155.pdf
.. _Gradient Computation of Image: https://en.wikipedia.org/wiki/Image_gradient
.. _R2 Score_Coefficient Determination: https://en.wikipedia.org/wiki/Coefficient_of_determination
.. _Rank of element tensor: https://github.com/scipy/scipy/blob/v1.6.2/scipy/stats/stats.py#L4140-L4303
.. _Mean Reciprocal Rank: https://en.wikipedia.org/wiki/Mean_reciprocal_rank
.. _BERT_score: https://github.com/Tiiiger/bert_score/blob/master/bert_score/utils.py
.. _Bert_score Evaluating Text Generation: https://arxiv.org/abs/1904.09675
.. _BLEU score: https://en.wikipedia.org/wiki/BLEU
.. _BLEU: http://www.aclweb.org/anthology/P02-1040.pdf
.. _Machine Translation Evolution: https://aclanthology.org/P04-1077.pdf
.. _Rouge score_Text Normalizition: https://github.com/google-research/google-research/blob/master/rouge/tokenize.py
.. _Calculate Rouge Score: https://en.wikipedia.org/wiki/ROUGE_(metric)
.. _Rouge Detail: https://aclanthology.org/W04-1013/
.. _Square Root of a Positive Definite Matrix: https://github.com/steveli/pytorch-sqrtm/blob/master/sqrtm.py
.. _Fid Score: https://github.com/photosynthesis-team/piq/blob/master/piq/fid.py
.. _Rethinking the Inception Architecture for ComputerVision: https://arxiv.org/abs/1512.00567
.. _GANs Trained by a Two Time-Scale: https://arxiv.org/abs/1706.08500
.. _Improved Techniques for Training GANs: https://arxiv.org/abs/1606.03498
.. _KID Score: https://github.com/toshas/torch-fidelity/blob/v0.3.0/torch_fidelity/metric_kid.py
.. _Demystifying MMD GANs: https://arxiv.org/abs/1801.01401
.. _Computes Peak Signal-to-Noise Ratio: https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
.. _Turn a Metric into a Bootstrapped: https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
.. _Metric Test for Reset: https://github.com/PyTorchLightning/pytorch-lightning/pull/7055
.. _Computes Mean Absolute Error: https://en.wikipedia.org/wiki/Mean_absolute_error
.. _Mean Absolute Percentage Error: https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
.. _mean squared error: https://en.wikipedia.org/wiki/Mean_squared_error
.. _Aggregate the statistics from multiple devices: https://stackoverflow.com/questions/68395368/estimate-running-correlation-on-multiple-nodes
.. _Pearson Correlation Coefficient: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
.. _Python ROUGE Implementation: https://pypi.org/project/rouge-score/
.. _Scikit_Learn-Ranking.py: https: //github.com/scikit-learn/scikit-learn/blob/master/sklearn/metrics/_ranking.py
.. _Verified Uncertainty Calibration: https://arxiv.org/abs/1909.10155
.. _SQuAD Metric: https://arxiv.org/pdf/1606.05250.pdf
.. _chrF score: https://aclanthology.org/W15-3049.pdf
.. _chrF++ score: https://aclanthology.org/W17-4770.pdf
.. _TER: https://aclanthology.org/2006.amta-papers.25.pdf
.. _ExtendedEditDistance: https://aclanthology.org/W19-5359.pdf
.. _MultiScaleSSIM: https://ece.uwaterloo.ca/~z70wang/publications/msssim.pdf
.. _UniversalImageQualityIndex: https://ieeexplore.ieee.org/document/995823
.. _SpectralDistortionIndex: https://www.semanticscholar.org/paper/Multispectral-and-panchromatic-data-fusion-without-Alparone-Aiazzi/b6db12e3785326577cb95fd743fecbf5bc66c7c9
.. _WMAPE: https://en.wikipedia.org/wiki/WMAPE
.. _governance:

TorchMetrics Governance
#######################

This document describes governance processes we follow in developing TorchMetrics.

Persons of Interest
*******************

Leads
-----
- Nicki Skafte (`skaftenicki <https://github.com/SkafteNicki>`_)
- Jirka Borovec (`Borda <https://github.com/Borda>`_)
- Justus Schock (`justusschock <https://github.com/justusschock>`_)


Core Maintainers
----------------
- Luca Di Liello (`lucadiliello <https://github.com/lucadiliello>`_)
- Daniel Stancl (`stancld <https://github.com/stancld>`_)
- Maxim Grechkin (`maximsch2 <https://github.com/maximsch2>`_)
- Changsheng Quan (`quancs <https://github.com/quancs>`_)


Alumni
------
- Ananya Harsh Jha (`ananyahjha93 <https://github.com/ananyahjha93>`_)
- Teddy Koker (`teddykoker <https://github.com/teddykoker>`_)


Releases
********

We release a new minor version (e.g., 0.5.0) every few months and bugfix releases if needed.
The minor versions contain new features, API changes, deprecations, removals, potential backward-incompatible
changes and also all previous bugfixes included in any bugfix release. With every release, we publish a changelog
where we list additions, removals, changed functionality and fixes.

Project Management and Decision Making
**************************************

The decision what goes into a release is governed by the :ref:`staff contributors and leaders <governance>` of
TorchMetrics development. Whenever possible, discussion happens publicly on GitHub and includes the whole community.
When a consensus is reached, staff and core contributors assign milestones and labels to the issue and/or pull request
and start tracking the development. It is possible that priorities change over time.

Commits to the project are exclusively to be added by pull requests on GitHub and anyone in the community is welcome to review them.
However, reviews submitted by
`code owners <https://github.com/PyTorchLightning/metrics/blob/master/.github/CODEOWNERS>`_
have higher weight and it is necessary to get the approval of code owners before a pull request can be merged.
Additional requirements may apply case by case.

API Evolution
*************

Torchmetrics development is driven by research and best practices in a rapidly developing field of AI and machine
learning. Change is inevitable and when it happens, the Torchmetric team is committed to minimizing user friction and
maximizing ease of transition from one version to the next. We take backward compatibility and reproducibility very
seriously.

For API removal, renaming or other forms of backward-incompatible changes, the procedure is:

#. A deprecation process is initiated at version X, producing warning messages at runtime and in the documentation.
#. Calls to the deprecated API remain unchanged in their function during the deprecation phase.
#. One minor versions in the future at version X+1 the breaking change takes effect.

The "X+1" rule is a recommendation and not a strict requirement. Longer deprecation cycles may apply for some cases.
.. TorchMetrics documentation master file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TorchMetrics documentation
==========================

.. include:: pages/brief_intro.rst


More reading
============

.. toctree::
   :maxdepth: 2
   :name: guide
   :caption: User Guide

   pages/quickstart
   pages/overview
   pages/implement
   pages/lightning

.. toctree::
   :maxdepth: 3
   :name: metrics
   :caption: Metrics API references

   references/modules
   references/functional

.. toctree::
   :maxdepth: 1
   :name: community
   :caption: Community

   governance
   generated/CODE_OF_CONDUCT.md
   generated/CONTRIBUTING.md
   generated/CHANGELOG.md

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. role:: hidden
    :class: hidden-section

.. include:: ../links.rst

##################
Functional metrics
##################

*****
Audio
*****

perceptual_evaluation_speech_quality [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.audio.pesq.perceptual_evaluation_speech_quality


permutation_invariant_training [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.permutation_invariant_training
    :noindex:


signal_distortion_ratio [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.signal_distortion_ratio
    :noindex:


scale_invariant_signal_distortion_ratio [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.scale_invariant_signal_distortion_ratio
    :noindex:


scale_invariant_signal_noise_ratio [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.scale_invariant_signal_noise_ratio
    :noindex:


signal_noise_ratio [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.signal_noise_ratio
    :noindex:


short_time_objective_intelligibility [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.audio.stoi.short_time_objective_intelligibility
    :noindex:


**************
Classification
**************

accuracy [func]
~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.accuracy
    :noindex:


auc [func]
~~~~~~~~~~

.. autofunction:: torchmetrics.functional.auc
    :noindex:


auroc [func]
~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.auroc
    :noindex:


average_precision [func]
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.average_precision
    :noindex:


calibration_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.calibration_error
    :noindex:


cohen_kappa [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.cohen_kappa
    :noindex:


confusion_matrix [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.confusion_matrix
    :noindex:


coverage_error [func]
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.coverage_error
    :noindex:


dice_score [func]
~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.dice_score
    :noindex:


f1_score [func]
~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.f1_score
    :noindex:


fbeta_score [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.fbeta_score
    :noindex:


hamming_distance [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.hamming_distance
    :noindex:


hinge_loss [func]
~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.hinge_loss
    :noindex:


jaccard_index [func]
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.jaccard_index
    :noindex:


kl_divergence [func]
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.kl_divergence
    :noindex:


label_ranking_average_precision [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.label_ranking_average_precision
    :noindex:


label_ranking_loss [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.label_ranking_loss
    :noindex:


matthews_corrcoef [func]
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.matthews_corrcoef
    :noindex:


roc [func]
~~~~~~~~~~

.. autofunction:: torchmetrics.functional.roc
    :noindex:


precision [func]
~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.precision
    :noindex:


precision_recall [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.precision_recall
    :noindex:


precision_recall_curve [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.precision_recall_curve
    :noindex:


recall [func]
~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.recall
    :noindex:


select_topk [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.utilities.data.select_topk
    :noindex:


specificity [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.specificity
    :noindex:


stat_scores [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.stat_scores
    :noindex:


to_categorical [func]
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.utilities.data.to_categorical
    :noindex:


to_onehot [func]
~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.utilities.data.to_onehot
    :noindex:


*****
Image
*****


error_relative_global_dimensionless_synthesis [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.error_relative_global_dimensionless_synthesis
    :noindex:


image_gradients [func]
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.image_gradients
    :noindex:


structural_similarity_index_measure [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.structural_similarity_index_measure
    :noindex:


multiscale_structural_similarity_index_measure [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.multiscale_structural_similarity_index_measure
    :noindex:


peak_signal_noise_ratio [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.peak_signal_noise_ratio
    :noindex:


spectral_angle_mapper [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.spectral_angle_mapper
    :noindex:


spectral_distortion_index [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.spectral_distortion_index


universal_image_quality_index [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.universal_image_quality_index
    :noindex:


**********
Regression
**********


cosine_similarity [func]
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.cosine_similarity
    :noindex:


explained_variance [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.explained_variance
    :noindex:


mean_absolute_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.mean_absolute_error
    :noindex:


mean_absolute_percentage_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.mean_absolute_percentage_error
    :noindex:


mean_squared_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.mean_squared_error
    :noindex:


mean_squared_log_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.mean_squared_log_error
    :noindex:


pearson_corrcoef [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.pearson_corrcoef
    :noindex:


r2_score [func]
~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.r2_score
    :noindex:


spearman_corrcoef [func]
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.spearman_corrcoef
    :noindex:


symmetric_mean_absolute_percentage_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.symmetric_mean_absolute_percentage_error
    :noindex:


tweedie_deviance_score [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.tweedie_deviance_score
    :noindex:


weighted_mean_absolute_percentage_error [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.weighted_mean_absolute_percentage_error
    :noindex:


****************
Pairwise Metrics
****************

pairwise_cosine_similarity [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.pairwise_cosine_similarity
    :noindex:


pairwise_euclidean_distance [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.pairwise_euclidean_distance
    :noindex:


pairwise_linear_similarity [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.pairwise_linear_similarity
    :noindex:


pairwise_manhattan_distance [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.pairwise_manhattan_distance
    :noindex:


*********
Retrieval
*********

retrieval_average_precision [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_average_precision
    :noindex:


retrieval_reciprocal_rank [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_reciprocal_rank
    :noindex:


retrieval_precision [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_precision
    :noindex:


retrieval_r_precision [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_r_precision
    :noindex:


retrieval_recall [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_recall
    :noindex:


retrieval_fall_out [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_fall_out
    :noindex:


retrieval_normalized_dcg [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_normalized_dcg
    :noindex:


retrieval_hit_rate [func]
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.retrieval_hit_rate
    :noindex:

****
Text
****

bert_score [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.text.bert.bert_score

bleu_score [func]
~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.bleu_score
    :noindex:

char_error_rate [func]
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.char_error_rate
    :noindex:

chrf_score [func]
~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.chrf_score
    :noindex:

extended_edit_distance [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.extended_edit_distance
    :noindex:

match_error_rate [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.match_error_rate
    :noindex:

rouge_score [func]
~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.text.rouge.rouge_score
    :noindex:

sacre_bleu_score [func]
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.sacre_bleu_score
    :noindex:

squad [func]
~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.squad
    :noindex:

translation_edit_rate [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.translation_edit_rate
    :noindex:

word_error_rate [func]
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.word_error_rate
    :noindex:

word_information_lost [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.word_information_lost
    :noindex:

word_information_preserved [func]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: torchmetrics.functional.word_information_preserved
    :noindex:
##############
Module metrics
##############

.. include:: ../links.rst

**********
Base class
**********

The base ``Metric`` class is an abstract base class that are used as the building block for all other Module
metrics.

.. autoclass:: torchmetrics.Metric
    :noindex:


*****************
Basic Aggregation
*****************

Torchmetrics comes with a number of metrics for aggregation of basic statistics: mean, max, min etc. of
either tensors or native python floats.

CatMetric
~~~~~~~~~

.. autoclass:: torchmetrics.CatMetric
    :noindex:

MaxMetric
~~~~~~~~~

.. autoclass:: torchmetrics.MaxMetric
    :noindex:

MeanMetric
~~~~~~~~~~

.. autoclass:: torchmetrics.MeanMetric
    :noindex:

MinMetric
~~~~~~~~~

.. autoclass:: torchmetrics.MinMetric
    :noindex:

SumMetric
~~~~~~~~~

.. autoclass:: torchmetrics.SumMetric
    :noindex:

*****
Audio
*****

For the purposes of audio metrics, inputs (predictions, targets) must have the same size.
If the input is 1D tensors the output will be a scalar. If the input is multi-dimensional with shape ``[...,time]``
the metric will be computed over the ``time`` dimension.

.. doctest::

    >>> import torch
    >>> from torchmetrics import SignalNoiseRatio
    >>> target = torch.tensor([3.0, -0.5, 2.0, 7.0])
    >>> preds = torch.tensor([2.5, 0.0, 2.0, 8.0])
    >>> snr = SignalNoiseRatio()
    >>> snr(preds, target)
    tensor(16.1805)

PerceptualEvaluationSpeechQuality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.audio.pesq.PerceptualEvaluationSpeechQuality

PermutationInvariantTraining
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.PermutationInvariantTraining
    :noindex:

SignalDistortionRatio
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SignalDistortionRatio
    :noindex:

ScaleInvariantSignalDistortionRatio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ScaleInvariantSignalDistortionRatio
    :noindex:

ScaleInvariantSignalNoiseRatio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ScaleInvariantSignalNoiseRatio
    :noindex:

SignalNoiseRatio
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SignalNoiseRatio
    :noindex:

ShortTimeObjectiveIntelligibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.audio.stoi.ShortTimeObjectiveIntelligibility
    :noindex:


**************
Classification
**************

Input types
~~~~~~~~~~~

For the purposes of classification metrics, inputs (predictions and targets) are split
into these categories (``N`` stands for the batch size and ``C`` for number of classes):

.. csv-table:: \*dtype ``binary`` means integers that are either 0 or 1
    :header: "Type", "preds shape", "preds dtype", "target shape", "target dtype"
    :widths: 20, 10, 10, 10, 10

    "Binary", "(N,)", "``float``", "(N,)", "``binary``\*"
    "Multi-class", "(N,)", "``int``", "(N,)", "``int``"
    "Multi-class with logits or probabilities", "(N, C)", "``float``", "(N,)", "``int``"
    "Multi-label", "(N, ...)", "``float``", "(N, ...)", "``binary``\*"
    "Multi-dimensional multi-class", "(N, ...)", "``int``", "(N, ...)", "``int``"
    "Multi-dimensional multi-class with logits or probabilities", "(N, C, ...)", "``float``", "(N, ...)", "``int``"

.. note::
    All dimensions of size 1 (except ``N``) are "squeezed out" at the beginning, so
    that, for example, a tensor of shape ``(N, 1)`` is treated as ``(N, )``.

When predictions or targets are integers, it is assumed that class labels start at 0, i.e.
the possible class labels are 0, 1, 2, 3, etc. Below are some examples of different input types

.. testcode::

    # Binary inputs
    binary_preds  = torch.tensor([0.6, 0.1, 0.9])
    binary_target = torch.tensor([1, 0, 2])

    # Multi-class inputs
    mc_preds  = torch.tensor([0, 2, 1])
    mc_target = torch.tensor([0, 1, 2])

    # Multi-class inputs with probabilities
    mc_preds_probs  = torch.tensor([[0.8, 0.2, 0], [0.1, 0.2, 0.7], [0.3, 0.6, 0.1]])
    mc_target_probs = torch.tensor([0, 1, 2])

    # Multi-label inputs
    ml_preds  = torch.tensor([[0.2, 0.8, 0.9], [0.5, 0.6, 0.1], [0.3, 0.1, 0.1]])
    ml_target = torch.tensor([[0, 1, 1], [1, 0, 0], [0, 0, 0]])


Using the multiclass parameter
------------------------------

In some cases, you might have inputs which appear to be (multi-dimensional) multi-class
but are actually binary/multi-label - for example, if both predictions and targets are
integer (binary) tensors. Or it could be the other way around, you want to treat
binary/multi-label inputs as 2-class (multi-dimensional) multi-class inputs.

For these cases, the metrics where this distinction would make a difference, expose the
``multiclass`` argument. Let's see how this is used on the example of
:class:`~torchmetrics.StatScores` metric.

First, let's consider the case with label predictions with 2 classes, which we want to
treat as binary.

.. testcode::

   from torchmetrics.functional import stat_scores

   # These inputs are supposed to be binary, but appear as multi-class
   preds  = torch.tensor([0, 1, 0])
   target = torch.tensor([1, 1, 0])

As you can see below, by default the inputs are treated
as multi-class. We can set ``multiclass=False`` to treat the inputs as binary -
which is the same as converting the predictions to float beforehand.

.. doctest::

    >>> stat_scores(preds, target, reduce='macro', num_classes=2)
    tensor([[1, 1, 1, 0, 1],
            [1, 0, 1, 1, 2]])
    >>> stat_scores(preds, target, reduce='macro', num_classes=1, multiclass=False)
    tensor([[1, 0, 1, 1, 2]])
    >>> stat_scores(preds.float(), target, reduce='macro', num_classes=1)
    tensor([[1, 0, 1, 1, 2]])

Next, consider the opposite example: inputs are binary (as predictions are probabilities),
but we would like to treat them as 2-class multi-class, to obtain the metric for both classes.

.. testcode::

   preds  = torch.tensor([0.2, 0.7, 0.3])
   target = torch.tensor([1, 1, 0])

In this case we can set ``multiclass=True``, to treat the inputs as multi-class.

.. doctest::

    >>> stat_scores(preds, target, reduce='macro', num_classes=1)
    tensor([[1, 0, 1, 1, 2]])
    >>> stat_scores(preds, target, reduce='macro', num_classes=2, multiclass=True)
    tensor([[1, 1, 1, 0, 1],
            [1, 0, 1, 1, 2]])

Accuracy
~~~~~~~~

.. autoclass:: torchmetrics.Accuracy
    :noindex:

AveragePrecision
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.AveragePrecision
    :noindex:

AUC
~~~

.. autoclass:: torchmetrics.AUC
    :noindex:

AUROC
~~~~~

.. autoclass:: torchmetrics.AUROC
    :noindex:

BinnedAveragePrecision
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.BinnedAveragePrecision
    :noindex:

BinnedPrecisionRecallCurve
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.BinnedPrecisionRecallCurve
    :noindex:

BinnedRecallAtFixedPrecision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.BinnedRecallAtFixedPrecision
    :noindex:

CalibrationError
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.CalibrationError
    :noindex:

CohenKappa
~~~~~~~~~~

.. autoclass:: torchmetrics.CohenKappa
    :noindex:

ConfusionMatrix
~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ConfusionMatrix
    :noindex:

CoverageError
~~~~~~~~~~~~~

.. autoclass:: torchmetrics.CoverageError
    :noindex:

F1Score
~~~~~~~

.. autoclass:: torchmetrics.F1Score
    :noindex:

FBetaScore
~~~~~~~~~~

.. autoclass:: torchmetrics.FBetaScore
    :noindex:

HammingDistance
~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.HammingDistance
    :noindex:

HingeLoss
~~~~~~~~~

.. autoclass:: torchmetrics.HingeLoss
    :noindex:

JaccardIndex
~~~~~~~~~~~~

.. autoclass:: torchmetrics.JaccardIndex
    :noindex:

KLDivergence
~~~~~~~~~~~~

.. autoclass:: torchmetrics.KLDivergence
    :noindex:

LabelRankingAveragePrecision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.LabelRankingAveragePrecision
    :noindex:

LabelRankingLoss
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.LabelRankingLoss
    :noindex:

MatthewsCorrCoef
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MatthewsCorrCoef
    :noindex:

Precision
~~~~~~~~~

.. autoclass:: torchmetrics.Precision
    :noindex:

PrecisionRecallCurve
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.PrecisionRecallCurve
    :noindex:

Recall
~~~~~~

.. autoclass:: torchmetrics.Recall
    :noindex:


ROC
~~~

.. autoclass:: torchmetrics.ROC
    :noindex:


Specificity
~~~~~~~~~~~

.. autoclass:: torchmetrics.Specificity
    :noindex:


StatScores
~~~~~~~~~~

.. autoclass:: torchmetrics.StatScores
    :noindex:

*****
Image
*****

Image quality metrics can be used to access the quality of synthetic generated images from machine
learning algorithms such as `Generative Adverserial Networks (GANs) <https://en.wikipedia.org/wiki/Generative_adversarial_network>`_.

ErrorRelativeGlobalDimensionlessSynthesis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.image.ergas.ErrorRelativeGlobalDimensionlessSynthesis
    :noindex:

FrechetInceptionDistance
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.image.fid.FrechetInceptionDistance
    :noindex:

InceptionScore
~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.image.inception.InceptionScore
    :noindex:

KernelInceptionDistance
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.image.kid.KernelInceptionDistance
    :noindex:

LearnedPerceptualImagePatchSimilarity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.image.lpip.LearnedPerceptualImagePatchSimilarity
    :noindex:

MultiScaleStructuralSimilarityIndexMeasure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MultiScaleStructuralSimilarityIndexMeasure
    :noindex:

PeakSignalNoiseRatio
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.PeakSignalNoiseRatio
    :noindex:


SpectralAngleMapper
~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SpectralAngleMapper
    :noindex:


SpectralDistortionIndex
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SpectralDistortionIndex
    :noindex:


StructuralSimilarityIndexMeasure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.StructuralSimilarityIndexMeasure
    :noindex:


UniversalImageQualityIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.UniversalImageQualityIndex
    :noindex:

*********
Detection
*********

Object detection metrics can be used to evaluate the predicted detections with given groundtruth detections on images.

MeanAveragePrecision
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.detection.mean_ap.MeanAveragePrecision
    :noindex:

**********
Regression
**********

CosineSimilarity
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.CosineSimilarity
    :noindex:


ExplainedVariance
~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ExplainedVariance
    :noindex:


MeanAbsoluteError
~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MeanAbsoluteError
    :noindex:


MeanAbsolutePercentageError
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MeanAbsolutePercentageError
    :noindex:


MeanSquaredError
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MeanSquaredError
    :noindex:


MeanSquaredLogError
~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MeanSquaredLogError
    :noindex:


PearsonCorrCoef
~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.PearsonCorrCoef
    :noindex:


R2Score
~~~~~~~

.. autoclass:: torchmetrics.R2Score
    :noindex:


SpearmanCorrCoef
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SpearmanCorrCoef
    :noindex:

SymmetricMeanAbsolutePercentageError
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SymmetricMeanAbsolutePercentageError
    :noindex:


TweedieDevianceScore
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.TweedieDevianceScore
    :noindex:


WeightedMeanAbsolutePercentageError
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.WeightedMeanAbsolutePercentageError
    :noindex:


*********
Retrieval
*********

Input details
~~~~~~~~~~~~~

For the purposes of retrieval metrics, inputs (indexes, predictions and targets) must have the same size
(``N`` stands for the batch size) and the following types:

.. csv-table::
    :header: "indexes shape", "indexes dtype", "preds shape", "preds dtype", "target shape", "target dtype"
    :widths: 10, 10, 10, 10, 10, 10

    "(N,...)", "``long``", "(N,...)", "``float``", "(N,...)", "``long`` or ``bool``"

.. note::
    All dimensions are flattened at the beginning, so
    that, for example, a tensor of shape ``(N, M)`` is treated as ``(N * M, )``.

In Information Retrieval you have a query that is compared with a variable number of documents. For each pair ``(Q_i, D_j)``,
a score is computed that measures the relevance of document ``D`` w.r.t. query ``Q``. Documents are then sorted by score
and you hope that relevant documents are scored higher. ``target`` contains the labels for the documents (relevant or not).

Since a query may be compared with a variable number of documents, we use ``indexes`` to keep track of which scores belong to
the set of pairs ``(Q_i, D_j)`` having the same query ``Q_i``.

.. note::
    `Retrieval` metrics are only intended to be used globally. This means that the average of the metric over each batch can be quite different
    from the metric computed on the whole dataset. For this reason, we suggest to compute the metric only when all the examples
    has been provided to the metric. When using `Pytorch Lightning`, we suggest to use ``on_step=False``
    and ``on_epoch=True`` in ``self.log`` or to place the metric calculation in ``training_epoch_end``, ``validation_epoch_end`` or ``test_epoch_end``.

.. doctest::

    >>> from torchmetrics import RetrievalMAP
    >>> # functional version works on a single query at a time
    >>> from torchmetrics.functional import retrieval_average_precision

    >>> # the first query was compared with two documents, the second with three
    >>> indexes = torch.tensor([0, 0, 1, 1, 1])
    >>> preds = torch.tensor([0.8, -0.4, 1.0, 1.4, 0.0])
    >>> target = torch.tensor([0, 1, 0, 1, 1])

    >>> rmap = RetrievalMAP() # or some other retrieval metric
    >>> rmap(preds, target, indexes=indexes)
    tensor(0.6667)

    >>> # the previous instruction is roughly equivalent to
    >>> res = []
    >>> # iterate over indexes of first and second query
    >>> for indexes in ([0, 1], [2, 3, 4]):
    ...     res.append(retrieval_average_precision(preds[indexes], target[indexes]))
    >>> torch.stack(res).mean()
    tensor(0.6667)


RetrievalMAP
~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalMAP
    :noindex:


RetrievalMRR
~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalMRR
    :noindex:


RetrievalPrecision
~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalPrecision
    :noindex:


RetrievalRPrecision
~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalRPrecision
    :noindex:


RetrievalRecall
~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalRecall
    :noindex:


RetrievalFallOut
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalFallOut
    :noindex:


RetrievalNormalizedDCG
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalNormalizedDCG
    :noindex:


RetrievalHitRate
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.RetrievalHitRate
    :noindex:

****
Text
****

BERTScore
~~~~~~~~~~

.. autoclass:: torchmetrics.text.bert.BERTScore
    :noindex:

BLEUScore
~~~~~~~~~

.. autoclass:: torchmetrics.BLEUScore
    :noindex:

CharErrorRate
~~~~~~~~~~~~~

.. autoclass:: torchmetrics.CharErrorRate
    :noindex:

CHRFScore
~~~~~~~~~

.. autoclass:: torchmetrics.CHRFScore
    :noindex:

ExtendedEditDistance
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ExtendedEditDistance
    :noindex:

MatchErrorRate
~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MatchErrorRate
    :noindex:

ROUGEScore
~~~~~~~~~~

.. autoclass:: torchmetrics.text.rouge.ROUGEScore
    :noindex:

SacreBLEUScore
~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.SacreBLEUScore
    :noindex:

SQuAD
~~~~~

.. autoclass:: torchmetrics.SQuAD
    :noindex:

TranslationEditRate
~~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.TranslationEditRate
    :noindex:

WordErrorRate
~~~~~~~~~~~~~

.. autoclass:: torchmetrics.WordErrorRate
    :noindex:

WordInfoLost
~~~~~~~~~~~~

.. autoclass:: torchmetrics.WordInfoLost
    :noindex:

WordInfoPreserved
~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.WordInfoPreserved
    :noindex:


********
Wrappers
********

Modular wrapper metrics are not metrics in themself, but instead take a metric and alter the internal logic
of the base metric.

BootStrapper
~~~~~~~~~~~~

.. autoclass:: torchmetrics.BootStrapper
    :noindex:

ClasswiseWrapper
~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.ClasswiseWrapper
    :noindex:

MetricTracker
~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MetricTracker
    :noindex:

MinMaxMetric
~~~~~~~~~~~~

.. autoclass:: torchmetrics.MinMaxMetric
    :noindex:

MultioutputWrapper
~~~~~~~~~~~~~~~~~~

.. autoclass:: torchmetrics.MultioutputWrapper
    :noindex:
.. role:: hidden
    :class: hidden-section
.. currentmodule:: {{ module }}


{{ name | underline}}

.. autoclass:: {{ name }}
    :members:


..
  autogenerated from source/_templates/classtemplate.rst
  note it does not have :inherited-members:
{{ name | escape | underline }}

.. currentmodule:: {{ fullname }}

{% block functions %}
{% if functions %}
.. rubric:: Functions

.. autosummary::
    :nosignatures:
{% for item in functions %}
    {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block classes %}
{% if classes %}
.. rubric:: Classes

.. autosummary::
    :nosignatures:
{% for item in classes %}
    {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block exceptions %}
{% if exceptions %}
.. rubric:: Exceptions

.. autosummary::
    :nosignatures:
{% for item in exceptions %}
    {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

.. automodule:: {{ fullname }}
.. _implement:

*********************
Implementing a Metric
*********************

To implement your own custom metric, subclass the base :class:`~torchmetrics.Metric` class and implement the following methods:

- ``__init__()``: Each state variable should be called using ``self.add_state(...)``.
- ``update()``: Any code needed to update the state given any inputs to the metric.
- ``compute()``: Computes a final value from the state of the metric.

We provide the remaining interface, such as ``reset()`` that will make sure to correctly reset all metric
states that have been added using ``add_state``. You should therefore not implement ``reset()`` yourself.
Additionally, adding metric states with ``add_state`` will make sure that states are correctly synchronized
in distributed settings (DDP). To see how metric states are synchronized across distributed processes,
refer to ``add_state()`` docs from the base ``Metric`` class.

Example implementation:

.. testcode::

    from torchmetrics import Metric

    class MyAccuracy(Metric):
        def __init__(self, dist_sync_on_step=False):
            super().__init__(dist_sync_on_step=dist_sync_on_step)

            self.add_state("correct", default=torch.tensor(0), dist_reduce_fx="sum")
            self.add_state("total", default=torch.tensor(0), dist_reduce_fx="sum")

        def update(self, preds: torch.Tensor, target: torch.Tensor):
            preds, target = self._input_format(preds, target)
            assert preds.shape == target.shape

            self.correct += torch.sum(preds == target)
            self.total += target.numel()

        def compute(self):
            return self.correct.float() / self.total


Internal implementation details
-------------------------------

This section briefly describes how metrics work internally. We encourage looking at the source code for more info.
Internally, TorchMetrics wraps the user defined ``update()`` and ``compute()`` method. We do this to automatically
synchronize and reduce metric states across multiple devices. More precisely, calling ``update()`` does the
following internally:

1. Clears computed cache.
2. Calls user-defined ``update()``.

Similarly, calling ``compute()`` does the following internally:

1. Syncs metric states between processes.
2. Reduce gathered metric states.
3. Calls the user defined ``compute()`` method on the gathered metric states.
4. Cache computed result.

From a user's standpoint this has one important side-effect: computed results are cached. This means that no
matter how many times ``compute`` is called after one and another, it will continue to return the same result.
The cache is first emptied on the next call to ``update``.

``forward`` serves the dual purpose of both returning the metric on the current data and updating the internal
metric state for accumulating over multiple batches. The ``forward()`` method achieves this by combining calls
to ``update`` and ``compute`` in the following way:

1. Calls ``update()`` to update the global metric state (for accumulation over multiple batches)
2. Caches the global state.
3. Calls ``reset()`` to clear global metric state.
4. Calls ``update()`` to update local metric state.
5. Calls ``compute()`` to calculate metric for current batch.
6. Restores the global state.

This procedure has the consequence of calling the user defined ``update`` **twice** during a single
forward call (one to update global statistics and one for getting the batch statistics).


---------

.. autoclass:: torchmetrics.Metric
    :members:

Contributing your metric to Torchmetrics
----------------------------------------

Wanting to contribute the metric you have implemented? Great, we are always open to adding more metrics to ``torchmetrics``
as long as they serve a general purpose. However, to keep all our metrics consistent we request that the implementation
and tests gets formatted in the following way:

1. Start by reading our `contribution guidelines <https://torchmetrics.readthedocs.io//en/latest/generated/CONTRIBUTING.html>`_.
2. First implement the functional backend. This takes cares of all the logic that goes into the metric. The code should
   be put into a single file placed under ``torchmetrics/functional/"domain"/"new_metric".py`` where ``domain`` is the type of
   metric (classification, regression, nlp etc) and ``new_metric`` is the name of the metric. In this file, there should be the
   following three functions:

  1. ``_new_metric_update(...)``: everything that has to do with type/shape checking and all logic required before distributed syncing need to go here.
  2. ``_new_metric_compute(...)``: all remaining logic.
  3. ``new_metric(...)``: essentially wraps the ``_update`` and ``_compute`` private functions into one public function that
     makes up the functional interface for the metric.

  .. note::
     The `functional accuracy <https://github.com/PyTorchLightning/metrics/blob/master/torchmetrics/functional/classification/accuracy.py>`_
     metric is a great example of this division of logic.

3. In a corresponding file placed in ``torchmetrics/"domain"/"new_metric".py`` create the module interface:

  1. Create a new module metric by subclassing ``torchmetrics.Metric``.
  2. In the ``__init__`` of the module call ``self.add_state`` for as many metric states are needed for the metric to
     proper accumulate metric statistics.
  3. The module interface should essentially call the private ``_new_metric_update(...)`` in its `update` method and similarly the
     ``_new_metric_compute(...)`` function in its ``compute``. No logic should really be implemented in the module interface.
     We do this to not have duplicate code to maintain.

   .. note::
     The module `Accuracy <https://github.com/PyTorchLightning/metrics/blob/master/torchmetrics/classification/accuracy.py>`_
     metric that corresponds to the above functional example showcases these steps.

4. Remember to add binding to the different relevant ``__init__`` files.

5. Testing is key to keeping ``torchmetrics`` trustworthy. This is why we have a very rigid testing protocol. This means
   that we in most cases require the metric to be tested against some other common framework (``sklearn``, ``scipy`` etc).

  1. Create a testing file in ``tests/"domain"/test_"new_metric".py``. Only one file is needed as it is intended to test
     both the functional and module interface.
  2. In that file, start by defining a number of test inputs that your metric should be evaluated on.
  3. Create a testclass ``class NewMetric(MetricTester)`` that inherits from ``tests.helpers.testers.MetricTester``.
     This testclass should essentially implement the ``test_"new_metric"_class`` and ``test_"new_metric"_fn`` methods that
     respectively tests the module interface and the functional interface.
  4. The testclass should be parameterized (using ``@pytest.mark.parametrize``) by the different test inputs defined initially.
     Additionally, the ``test_"new_metric"_class`` method should also be parameterized with an ``ddp`` parameter such that it gets
     tested in a distributed setting. If your metric has additional parameters, then make sure to also parameterize these
     such that different combinations of inputs and parameters gets tested.
  5. (optional) If your metric raises any exception, please add tests that showcase this.

  .. note::
    The `test file for accuracy <https://github.com/PyTorchLightning/metrics/blob/master/tests/classification/test_accuracy.py>`_ metric
    shows how to implement such tests.

If you only can figure out part of the steps, do not fear to send a PR. We will much rather receive working
metrics that are not formatted exactly like our codebase, than not receiving any. Formatting can always be applied.
We will gladly guide and/or help implement the remaining :]
TorchMetrics is a collection of Machine learning metrics for distributed, scalable PyTorch models and an easy-to-use API to create custom metrics. It offers the following benefits:

* Optimized for distributed-training
* A standardized interface to increase reproducibility
* Reduces Boilerplate
* Distributed-training compatible
* Rigorously tested
* Automatic accumulation over batches
* Automatic synchronization between multiple devices

You can use TorchMetrics in any PyTorch model, or with in `PyTorch Lightning <https://pytorch-lightning.readthedocs.io/en/stable/>`_ to enjoy additional features:

* This means that your data will always be placed on the same device as your metrics.
* Native support for logging metrics in Lightning to reduce even more boilerplate.

Using TorchMetrics
******************

Module metrics
~~~~~~~~~~~~~~

.. testcode::

    import torch
    import torchmetrics

    # initialize metric
    metric = torchmetrics.Accuracy()

    n_batches = 10
    for i in range(n_batches):
        # simulate a classification problem
        preds = torch.randn(10, 5).softmax(dim=-1)
        target = torch.randint(5, (10,))
        # metric on current batch
        acc = metric(preds, target)
        print(f"Accuracy on batch {i}: {acc}")

    # metric on all batches using custom accumulation
    acc = metric.compute()
    print(f"Accuracy on all data: {acc}")

.. testoutput::
   :hide:
   :options: +ELLIPSIS, +NORMALIZE_WHITESPACE

    Accuracy on batch ...

Module metric usage remains the same when using multiple GPUs or multiple nodes.


Functional metrics
~~~~~~~~~~~~~~~~~~

.. testcode::

    import torch
    import torchmetrics

    # simulate a classification problem
    preds = torch.randn(10, 5).softmax(dim=-1)
    target = torch.randint(5, (10,))

    acc = torchmetrics.functional.accuracy(preds, target)


Implementing a metric
~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    class MyAccuracy(Metric):
        def __init__(self, dist_sync_on_step=False):
            # call `self.add_state`for every internal state that is needed for the metrics computations
            # dist_reduce_fx indicates the function that should be used to reduce
            # state from multiple processes
            super().__init__(dist_sync_on_step=dist_sync_on_step)

            self.add_state("correct", default=torch.tensor(0), dist_reduce_fx="sum")
            self.add_state("total", default=torch.tensor(0), dist_reduce_fx="sum")

        def update(self, preds: torch.Tensor, target: torch.Tensor):
            # update metric states
            preds, target = self._input_format(preds, target)
            assert preds.shape == target.shape

            self.correct += torch.sum(preds == target)
            self.total += target.numel()

        def compute(self):
            # compute final result
            return self.correct.float() / self.total
.. testsetup:: *

    import torch
    from pytorch_lightning.core.lightning import LightningModule

########
Overview
########

The ``torchmetrics`` is a Metrics API created for easy metric development and usage in
PyTorch and PyTorch Lightning. It is rigorously tested for all edge cases and includes a growing list of
common metric implementations.

The metrics API provides ``update()``, ``compute()``, ``reset()`` functions to the user. The metric :ref:`base class <references/modules:base class>` inherits
:class:`torch.nn.Module` which allows us to call ``metric(...)`` directly. The ``forward()`` method of the base ``Metric`` class
serves the dual purpose of calling ``update()`` on its input and simultaneously returning the value of the metric over the
provided input.

These metrics work with DDP in PyTorch and PyTorch Lightning by default. When ``.compute()`` is called in
distributed mode, the internal state of each metric is synced and reduced across each process, so that the
logic present in ``.compute()`` is applied to state information from all processes.

This metrics API is independent of PyTorch Lightning. Metrics can directly be used in PyTorch as shown in the example:

.. code-block:: python

    from torchmetrics.classification import Accuracy

    train_accuracy = Accuracy()
    valid_accuracy = Accuracy()

    for epoch in range(epochs):
        for x, y in train_data:
            y_hat = model(x)

            # training step accuracy
            batch_acc = train_accuracy(y_hat, y)
            print(f"Accuracy of batch{i} is {batch_acc}")

        for x, y in valid_data:
            y_hat = model(x)
            valid_accuracy.update(y_hat, y)

        # total accuracy over all training batches
        total_train_accuracy = train_accuracy.compute()

        # total accuracy over all validation batches
        total_valid_accuracy = valid_accuracy.compute()

        print(f"Training acc for epoch {epoch}: {total_train_accuracy}")
        print(f"Validation acc for epoch {epoch}: {total_valid_accuracy}")

        # Reset metric states after each epoch
        train_accuracy.reset()
        valid_accuracy.reset()

.. note::

    Metrics contain internal states that keep track of the data seen so far.
    Do not mix metric states across training, validation and testing.
    It is highly recommended to re-initialize the metric per mode as
    shown in the examples above.

.. note::

    Metric states are **not** added to the models ``state_dict`` by default.
    To change this, after initializing the metric, the method ``.persistent(mode)`` can
    be used to enable (``mode=True``) or disable (``mode=False``) this behaviour.


*******************
Metrics and devices
*******************

Metrics are simple subclasses of :class:`~torch.nn.Module` and their metric states behave
similar to buffers and parameters of modules. This means that metrics states should
be moved to the same device as the input of the metric:

.. code-block:: python

    from torchmetrics import Accuracy

    target = torch.tensor([1, 1, 0, 0], device=torch.device("cuda", 0))
    preds = torch.tensor([0, 1, 0, 0], device=torch.device("cuda", 0))

    # Metric states are always initialized on cpu, and needs to be moved to
    # the correct device
    confmat = Accuracy(num_classes=2).to(torch.device("cuda", 0))
    out = confmat(preds, target)
    print(out.device) # cuda:0

However, when **properly defined** inside a :class:`~torch.nn.Module` or
:class:`~pytorch_lightning.core.lightning.LightningModule` the metric will be be automatically move
to the same device as the the module when using ``.to(device)``.  Being
**properly defined** means that the metric is correctly identified as a child module of the
model (check ``.children()`` attribute of the model). Therefore, metrics cannot be placed
in native python ``list`` and ``dict``, as they will not be correctly identified
as child modules. Instead of ``list`` use :class:`~torch.nn.ModuleList` and instead of
``dict`` use :class:`~torch.nn.ModuleDict`. Furthermore, when working with multiple metrics
the native `MetricCollection`_ module can also be used to wrap multiple metrics.

.. testcode::

    from torchmetrics import Accuracy, MetricCollection

    class MyModule(torch.nn.Module):
        def __init__(self):
            ...
            # valid ways metrics will be identified as child modules
            self.metric1 = Accuracy()
            self.metric2 = nn.ModuleList(Accuracy())
            self.metric3 = nn.ModuleDict({'accuracy': Accuracy()})
            self.metric4 = MetricCollection([Accuracy()]) # torchmetrics build-in collection class

        def forward(self, batch):
            data, target = batch
            preds = self(data)
            ...
            val1 = self.metric1(preds, target)
            val2 = self.metric2[0](preds, target)
            val3 = self.metric3['accuracy'](preds, target)
            val4 = self.metric4(preds, target)

You can always check which device the metric is located on using the `.device` property.

Metrics in Dataparallel (DP) mode
=================================

When using metrics in `Dataparallel (DP) <https://pytorch.org/docs/stable/generated/torch.nn.DataParallel.html#torch.nn.DataParallel>`_
mode, one should be aware DP will both create and clean-up replicas of Metric objects during a single forward pass.
This has the consequence, that the metric state of the replicas will as default be destroyed before we can sync
them. It is therefore recommended, when using metrics in DP mode, to initialize them with ``dist_sync_on_step=True``
such that metric states are synchonized between the main process and the replicas before they are destroyed.

Addtionally, if metrics are used together with a `LightningModule` the metric update/logging should be done
in the ``<mode>_step_end`` method (where ``<mode>`` is either ``training``, ``validation`` or ``test``), else
it will lead to wrong accumulation. In practice do the following:

.. testcode::

    def training_step(self, batch, batch_idx):
        data, target = batch
        preds = self(data)
        ...
        return {'loss': loss, 'preds': preds, 'target': target}

    def training_step_end(self, outputs):
        #update and log
        self.metric(outputs['preds'], outputs['target'])
        self.log('metric', self.metric)

Metrics in Distributed Data Parallel (DDP) mode
===============================================

When using metrics in `Distributed Data Parallel (DDP) <https://pytorch.org/docs/stable/generated/torch.nn.parallel.DistributedDataParallel.html>`_
mode, one should be aware that DDP will add additional samples to your dataset if the size of your dataset is
not equally divisible by ``batch_size * num_processors``. The added samples will always be replicates of datapoints
already in your dataset. This is done to secure an equal load for all processes. However, this has the consequence
that the calculated metric value will be sligtly bias towards those replicated samples, leading to a wrong result.

During training and/or validation this may not be important, however it is highly recommended when evaluating
the test dataset to only run on a single gpu or use a `join <https://pytorch.org/docs/stable/_modules/torch/nn/parallel/distributed.html#DistributedDataParallel.join>`_
context in conjunction with DDP to prevent this behaviour.

****************************
Metrics and 16-bit precision
****************************

Most metrics in our collection can be used with 16-bit precision (``torch.half``) tensors. However, we have found
the following limitations:

* In general ``pytorch`` had better support for 16-bit precision much earlier on GPU than CPU. Therefore, we
  recommend that anyone that want to use metrics with half precision on CPU, upgrade to atleast pytorch v1.6
  where support for operations such as addition, subtraction, multiplication ect. was added.
* Some metrics does not work at all in half precision on CPU. We have explicitly stated this in their docstring,
  but they are also listed below:

  - :ref:`references/modules:PeakSignalNoiseRatio` and :ref:`references/functional:peak_signal_noise_ratio [func]`
  - :ref:`references/modules:StructuralSimilarityIndexMeasure` and :ref:`references/functional:structural_similarity_index_measure [func]`
  - :ref:`references/modules:KLDivergence` and :ref:`references/functional:kl_divergence [func]`

You can always check the precision/dtype of the metric by checking the `.dtype` property.

******************
Metric Arithmetics
******************

Metrics support most of python built-in operators for arithmetic, logic and bitwise operations.

For example for a metric that should return the sum of two different metrics, implementing a new metric is an
overhead that is not necessary. It can now be done with:

.. code-block:: python

    first_metric = MyFirstMetric()
    second_metric = MySecondMetric()

    new_metric = first_metric + second_metric

``new_metric.update(*args, **kwargs)`` now calls update of ``first_metric`` and ``second_metric``. It forwards
all positional arguments but forwards only the keyword arguments that are available in respective metric's update
declaration. Similarly ``new_metric.compute()`` now calls compute of ``first_metric`` and ``second_metric`` and
adds the results up. It is important to note that all implemented operations always returns a new metric object. This means
that the line ``first_metric == second_metric`` will not return a bool indicating if ``first_metric`` and ``second_metric``
is the same metric, but will return a new metric that checks if the ``first_metric.compute() == second_metric.compute()``.

This pattern is implemented for the following operators (with ``a`` being metrics and ``b`` being metrics, tensors, integer or floats):

* Addition (``a + b``)
* Bitwise AND (``a & b``)
* Equality (``a == b``)
* Floordivision (``a // b``)
* Greater Equal (``a >= b``)
* Greater (``a > b``)
* Less Equal (``a <= b``)
* Less (``a < b``)
* Matrix Multiplication (``a @ b``)
* Modulo (``a % b``)
* Multiplication (``a * b``)
* Inequality (``a != b``)
* Bitwise OR (``a | b``)
* Power (``a ** b``)
* Subtraction (``a - b``)
* True Division (``a / b``)
* Bitwise XOR (``a ^ b``)
* Absolute Value (``abs(a)``)
* Inversion (``~a``)
* Negative Value (``neg(a)``)
* Positive Value (``pos(a)``)
* Indexing (``a[0]``)

.. note::

    Some of these operations are only fully supported from Pytorch v1.4 and onwards, explicitly we found:
    ``add``, ``mul``, ``rmatmul``, ``rsub``, ``rmod``


****************
MetricCollection
****************

In many cases it is beneficial to evaluate the model output by multiple metrics.
In this case the ``MetricCollection`` class may come in handy. It accepts a sequence
of metrics and wraps these into a single callable metric class, with the same
interface as any other metric.

Example:

.. testcode::

    from torchmetrics import MetricCollection, Accuracy, Precision, Recall
    target = torch.tensor([0, 2, 0, 2, 0, 1, 0, 2])
    preds = torch.tensor([2, 1, 2, 0, 1, 2, 2, 2])
    metric_collection = MetricCollection([
        Accuracy(),
        Precision(num_classes=3, average='macro'),
        Recall(num_classes=3, average='macro')
    ])
    print(metric_collection(preds, target))

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    {'Accuracy': tensor(0.1250),
     'Precision': tensor(0.0667),
     'Recall': tensor(0.1111)}

Similarly it can also reduce the amount of code required to log multiple metrics
inside your LightningModule

.. testcode::

    from torchmetrics import Accuracy, MetricCollection, Precision, Recall

    class MyModule(LightningModule):
        def __init__(self):
            metrics = MetricCollection([Accuracy(), Precision(), Recall()])
            self.train_metrics = metrics.clone(prefix='train_')
            self.valid_metrics = metrics.clone(prefix='val_')

        def training_step(self, batch, batch_idx):
            logits = self(x)
            # ...
            output = self.train_metrics(logits, y)
            # use log_dict instead of log
            # metrics are logged with keys: train_Accuracy, train_Precision and train_Recall
            self.log_dict(output)

        def validation_step(self, batch, batch_idx):
            logits = self(x)
            # ...
            output = self.valid_metrics(logits, y)
            # use log_dict instead of log
            # metrics are logged with keys: val_Accuracy, val_Precision and val_Recall
            self.log_dict(output)

.. note::

    `MetricCollection` as default assumes that all the metrics in the collection
    have the same call signature. If this is not the case, input that should be
    given to different metrics can given as keyword arguments to the collection.

An additional advantage of using the ``MetricCollection`` object is that it will
automatically try to reduce the computations needed by finding groups of metrics
that share the same underlying metric state. If such a group of metrics is found only one
of them is actually updated and the updated state will be broadcasted to the rest
of the metrics within the group. In the example above, this will lead to a 2x-3x lower computational
cost compared to disabling this feature. However, this speedup comes with a fixed cost upfront, where
the state-groups have to be determined after the first update. This overhead can be significantly higher then gains speed-up for very
a low number of steps (approx. up to 100) but still leads to an overall speedup for everything beyond that.
In case the groups are known beforehand, these can also be set manually to avoid this extra cost of the
dynamic search. See the *compute_groups* argument in the class docs below for more information on this topic.

.. autoclass:: torchmetrics.MetricCollection
    :noindex:


****************************
Module vs Functional Metrics
****************************

The functional metrics follow the simple paradigm input in, output out.
This means they don't provide any advanced mechanisms for syncing across DDP nodes or aggregation over batches.
They simply compute the metric value based on the given inputs.

Also, the integration within other parts of PyTorch Lightning will never be as tight as with the Module-based interface.
If you look for just computing the values, the functional metrics are the way to go.
However, if you are looking for the best integration and user experience, please consider also using the Module interface.


*****************************
Metrics and differentiability
*****************************

Metrics support backpropagation, if all computations involved in the metric calculation
are differentiable. All modular metrics have a property that determines if a metric is
differentiable or not.

However, note that the cached state is detached from the computational
graph and cannot be back-propagated. Not doing this would mean storing the computational
graph for each update call, which can lead to out-of-memory errors.
In practise this means that:

.. code-block:: python

    metric = MyMetric()
    val = metric(pred, target) # this value can be back-propagated
    val = metric.compute() # this value cannot be back-propagated

A functional metric is differentiable if its corresponding modular metric is differentiable.

.. _Metric kwargs:

************************
Advanced metric settings
************************

The following is a list of additional arguments that can be given to any metric class (in the ``**kwargs`` argument)
that will alter how metric states are stored and synced.

If you are running metrics on GPU and are encountering that you are running out of GPU VRAM then the following
argument can help:

- ``compute_on_cpu`` will automatically move the metric states to cpu after calling ``update``, making sure that
  GPU memory is not filling up. The consequence will be that the ``compute`` method will be called on CPU instead
  of GPU. Only applies to metric states that are lists.

If you are running in a distributed environment, ``TorchMetrics`` will automatically take care of the distributed
synchronization for you. However, the following three keyword arguments can be given to any metric class for
further control over the distributed aggregation:

- ``dist_sync_on_step``: This argument is ``bool`` that indicates if the metric should syncronize between
  different devices every time ``forward`` is called. Setting this to ``True`` is in general not recommended
  as syncronization is an expensive operation to do after each batch.

- ``process_group``: By default we syncronize across the *world* i.e. all proceses being computed on. You
  can provide an ``torch._C._distributed_c10d.ProcessGroup`` in this argument to specify exactly what
  devices should be syncronized over.

- ``dist_sync_fn``: By default we use :func:`torch.distributed.all_gather` to perform the synchronization between
  devices. Provide another callable function for this argument to perform custom distributed synchronization.
###########
Quick Start
###########

TorchMetrics is a collection of 80+ PyTorch metrics implementations and an easy-to-use API to create custom metrics. It offers:

* A standardized interface to increase reproducibility
* Reduces Boilerplate
* Distributed-training compatible
* Rigorously tested
* Automatic accumulation over batches
* Automatic synchronization between multiple devices

You can use TorchMetrics in any PyTorch model, or within `PyTorch Lightning <https://pytorch-lightning.readthedocs.io/en/stable/>`_ to enjoy additional features:

* This means that your data will always be placed on the same device as your metrics.
* Native support for logging metrics in Lightning to reduce even more boilerplate.

Install
*******

You can install TorchMetrics using pip or conda:

.. code-block:: bash

    # Python Package Index (PyPI)
    pip install torchmetrics
    # Conda
    conda install -c conda-forge torchmetrics

Eventually if there is a missing PyTorch wheel for your OS or Python version you can simply compile `PyTorch from source <https://github.com/pytorch/pytorch>`_:

.. code-block:: bash

    # Optional if you do not need compile GPU support
    export USE_CUDA=0  # just to keep it simple
    # you can install the latest state from master
    pip install git+https://github.com/pytorch/pytorch.git
    # OR set a particular PyTorch release
    pip install git+https://github.com/pytorch/pytorch.git@<release-tag>
    # and finalize with installing TorchMetrics
    pip install torchmetrics


Using TorchMetrics
******************

Functional metrics
~~~~~~~~~~~~~~~~~~

Similar to `torch.nn <https://pytorch.org/docs/stable/nn>`_, most metrics have both a class-based and a functional version.
The functional versions implement the basic operations required for computing each metric.
They are simple python functions that as input take `torch.tensors <https://pytorch.org/docs/stable/tensors.html>`_
and return the corresponding metric as a ``torch.tensor``.
The code-snippet below shows a simple example for calculating the accuracy using the functional interface:

.. testcode::

    import torch
    # import our library
    import torchmetrics

    # simulate a classification problem
    preds = torch.randn(10, 5).softmax(dim=-1)
    target = torch.randint(5, (10,))

    acc = torchmetrics.functional.accuracy(preds, target)

Module metrics
~~~~~~~~~~~~~~

Nearly all functional metrics have a corresponding class-based metric that calls it a functional counterpart underneath. The class-based metrics are characterized by having one or more internal metrics states (similar to the parameters of the PyTorch module) that allow them to offer additional functionalities:

* Accumulation of multiple batches
* Automatic synchronization between multiple devices
* Metric arithmetic

The code below shows how to use the class-based interface:

.. testcode::

    import torch
    # import our library
    import torchmetrics

    # initialize metric
    metric = torchmetrics.Accuracy()

    n_batches = 10
    for i in range(n_batches):
        # simulate a classification problem
        preds = torch.randn(10, 5).softmax(dim=-1)
        target = torch.randint(5, (10,))
        # metric on current batch
        acc = metric(preds, target)
        print(f"Accuracy on batch {i}: {acc}")

    # metric on all batches using custom accumulation
    acc = metric.compute()
    print(f"Accuracy on all data: {acc}")

    # Reseting internal state such that metric ready for new data
    metric.reset()

.. testoutput::
   :hide:
   :options: +ELLIPSIS, +NORMALIZE_WHITESPACE

    Accuracy on batch ...


Implementing your own metric
****************************

Implementing your own metric is as easy as subclassing a :class:`torch.nn.Module`. Simply, subclass :class:`~torchmetrics.Metric` and do the following:

1. Implement ``__init__`` where you call ``self.add_state`` for every internal state that is needed for the metrics computations
2. Implement ``update`` method, where all logic that is necessary for updating metric states go
3. Implement ``compute`` method, where the final metric computations happens

For practical examples and more info about implementing a metric, please see this :ref:`page <implement>`.


Development Environment
~~~~~~~~~~~~~~~~~~~~~~~

TorchMetrics provides a `Devcontainer <https://code.visualstudio.com/docs/remote/containers>`_ configuration for `Visual Studio Code <https://code.visualstudio.com/>`_ to use a `Docker container <https://www.docker.com/>`_ as a pre-configured development environment.
This avoids struggles setting up a development environment and makes them reproducible and consistent.
Please follow the `installation instructions <https://code.visualstudio.com/docs/remote/containers#_installation>`_ and make yourself familiar with the `container tutorials <https://code.visualstudio.com/docs/remote/containers-tutorial>`_ if you want to use them.
In order to use GPUs, you can enable them within the ``.devcontainer/devcontainer.json`` file.
.. testsetup:: *

    import torch
    from torch.nn import Module
    from pytorch_lightning.core.lightning import LightningModule
    from torchmetrics import Metric

#################################
TorchMetrics in PyTorch Lightning
#################################

TorchMetrics was originally created as part of `PyTorch Lightning <https://github.com/PyTorchLightning/pytorch-lightning>`_, a powerful deep learning research
framework designed for scaling models without boilerplate.

.. note::

    TorchMetrics always offers compatibility with the last 2 major PyTorch Lightning versions, but we recommend to always keep both frameworks
    up-to-date for the best experience.

While TorchMetrics was built to be used with native PyTorch, using TorchMetrics with Lightning offers additional benefits:

* Modular metrics are automatically placed on the correct device when properly defined inside a LightningModule.
  This means that your data will always be placed on the same device as your metrics. No need to call ``.to(device)`` anymore!
* Native support for logging metrics in Lightning using
  `self.log <https://pytorch-lightning.readthedocs.io/en/stable/extensions/logging.html#logging-from-a-lightningmodule>`_ inside
  your LightningModule.
* The ``.reset()`` method of the metric will automatically be called at the end of an epoch.

The example below shows how to use a metric in your `LightningModule <https://pytorch-lightning.readthedocs.io/en/stable/common/lightning_module.html>`_:

.. testcode:: python

    class MyModel(LightningModule):

        def __init__(self):
            ...
            self.accuracy = torchmetrics.Accuracy()

        def training_step(self, batch, batch_idx):
            x, y = batch
            preds = self(x)
            ...
            # log step metric
            self.accuracy(preds, y)
            self.log('train_acc_step', self.accuracy)
            ...

        def training_epoch_end(self, outs):
            # log epoch metric
            self.log('train_acc_epoch', self.accuracy)

Metric logging in Lightning happens through the ``self.log`` or ``self.log_dict`` method. Both methods only support the logging of *scalar-tensors*.
While the vast majority of metrics in torchmetrics returns a scalar tensor, some metrics such as :class:`~torchmetrics.ConfusionMatrix`, :class:`~torchmetrics.ROC`,
:class:`~torchmetrics.MeanAveragePrecision`, :class:`~torchmetrics.ROUGEScore` return outputs that are non-scalar tensors (often dicts or list of tensors) and should therefore be
dealt with separately. For info about the return type and shape please look at the documentation for the ``compute`` method for each metric you want to log.

********************
Logging TorchMetrics
********************

Logging metrics can be done in two ways: either logging the metric object directly or the computed metric values. When :class:`~torchmetrics.Metric` objects, which return a scalar tensor
are logged directly in Lightning using the LightningModule `self.log <https://pytorch-lightning.readthedocs.io/en/stable/extensions/logging.html#logging-from-a-lightningmodule>`_ method,
Lightning will log the metric based on ``on_step`` and ``on_epoch`` flags present in ``self.log(...)``. If ``on_epoch`` is True, the logger automatically logs the end of epoch metric
value by calling ``.compute()``.

.. note::

    ``sync_dist``, ``sync_dist_op``, ``sync_dist_group``, ``reduce_fx`` and ``tbptt_reduce_fx``
    flags from ``self.log(...)`` don't affect the metric logging in any manner. The metric class
    contains its own distributed synchronization logic.

    This however is only true for metrics that inherit the base class ``Metric``,
    and thus the functional metric API provides no support for in-built distributed synchronization
    or reduction functions.


.. testcode:: python

    class MyModule(LightningModule):

        def __init__(self):
            ...
            self.train_acc = torchmetrics.Accuracy()
            self.valid_acc = torchmetrics.Accuracy()

        def training_step(self, batch, batch_idx):
            x, y = batch
            preds = self(x)
            ...
            self.train_acc(preds, y)
            self.log('train_acc', self.train_acc, on_step=True, on_epoch=False)

        def validation_step(self, batch, batch_idx):
            logits = self(x)
            ...
            self.valid_acc(logits, y)
            self.log('valid_acc', self.valid_acc, on_step=True, on_epoch=True)

As an alternative to logging the metric object and letting Lightning take care of when to reset the metric etc. you can also manually log the output
of the metrics.

.. testcode:: python

    class MyModule(LightningModule):

        def __init__(self):
            ...
            self.train_acc = torchmetrics.Accuracy()
            self.valid_acc = torchmetrics.Accuracy()

        def training_step(self, batch, batch_idx):
            x, y = batch
            preds = self(x)
            ...
            batch_value = self.train_acc(preds, y)
            self.log('train_acc_step', batch_value)

        def training_epoch_end(self, outputs):
            self.train_acc.reset()

        def validation_step(self, batch, batch_idx):
            logits = self(x)
            ...
            self.valid_acc.update(logits, y)

        def validation_epoch_end(self, outputs):
            self.log('valid_acc_epoch', self.valid_acc.compute())
            self.valid_acc.reset()

Note that logging metrics this way will require you to manually reset the metrics at the end of the epoch yourself. In general, we recommend logging
the metric object to make sure that metrics are correctly computed and reset. Additionally, we highly recommend that the two ways of logging are not
mixed as it can lead to wrong results.

.. note::

    When using any Modular metric, calling ``self.metric(...)`` or ``self.metric.forward(...)`` serves the dual purpose of calling ``self.metric.update()``
    on its input and simultaneously returning the metric value over the provided input. So if you are logging a metric *only* on epoch-level (as in the
    example above), it is recommended to call ``self.metric.update()`` directly to avoid the extra computation.

    .. testcode:: python

        class MyModule(LightningModule):

            def __init__(self):
                ...
                self.valid_acc = torchmetrics.Accuracy()

            def validation_step(self, batch, batch_idx):
                logits = self(x)
                ...
                self.valid_acc.update(logits, y)
                self.log('valid_acc', self.valid_acc, on_step=True, on_epoch=True)


***************
Common Pitfalls
***************

The following contains a list of pitfalls to be aware of:

* If using metrics in data parallel mode (dp), the metric update/logging should be done
  in the ``<mode>_step_end`` method (where ``<mode>`` is either ``training``, ``validation``
  or ``test``). This is because ``dp`` split the batches during the forward pass and metric states are destroyed after each forward pass, thus leading to wrong accumulation. In practice do the following:

.. testcode:: python

    class MyModule(LightningModule):

        def training_step(self, batch, batch_idx):
            data, target = batch
            preds = self(data)
            # ...
            return {'loss': loss, 'preds': preds, 'target': target}

        def training_step_end(self, outputs):
            # update and log
            self.metric(outputs['preds'], outputs['target'])
            self.log('metric', self.metric)

* Modular metrics contain internal states that should belong to only one DataLoader. In case you are using multiple DataLoaders,
  it is recommended to initialize a separate modular metric instances for each DataLoader and use them separately. The same holds
  for using seperate metrics for training, validation and testing.

.. testcode:: python

    class MyModule(LightningModule):

        def __init__(self):
            ...
            self.val_acc = nn.ModuleList([torchmetrics.Accuracy() for _ in range(2)])

        def val_dataloader(self):
            return [DataLoader(...), DataLoader(...)]

        def validation_step(self, batch, batch_idx, dataloader_idx):
            x, y = batch
            preds = self(x)
            ...
            self.val_acc[dataloader_idx](preds, y)
            self.log('val_acc', self.val_acc[dataloader_idx])

* Mixing the two logging methods by calling ``self.log("val", self.metric)`` in ``{training}/{val}/{test}_step`` method and
  then calling ``self.log("val", self.metric.compute())`` in the corresponding ``{training}/{val}/{test}_epoch_end`` method.
  Because the object is logged in the first case, Lightning will reset the metric before calling the second line leading to
  errors or nonsense results.

* Calling ``self.log("val", self.metric(preds, target))`` with the intention of logging the metric object. Because
  ``self.metric(preds, target)`` corresponds to calling the forward method, this will return a tensor and not the
  metric object. Such logging will be wrong in this case. Instead it is important to seperate into seperate lines:

.. testcode:: python

    def training_step(self, batch, batch_idx):
        x, y = batch
        preds = self(x)
        ...
        # log step metric
        self.accuracy(preds, y)  # compute metrics
        self.log('train_acc_step', self.accuracy)  # log metric object
