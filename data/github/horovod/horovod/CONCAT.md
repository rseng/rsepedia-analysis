# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased] - YYYY-MM-DD

### Added
- Added Elastic keyword parameters to RayExecutor API: This API supports both static(non-elastic) and elastic horovod jobs. This resolves issue:
[#3190](https://github.com/horovod/horovod/issues/3190).

- TensorFlow: Added in-place broadcasting of variables. ([#3128](https://github.com/horovod/horovod/pull/3128))

- Added support for resurrecting blacklisted hosts. ([#3319](https://github.com/horovod/horovod/pull/3319))

### Changed

- Moved to CMake version 3.13 with first-class CUDA language support and re-enabled parallelized builds. ([#3261](https://github.com/horovod/horovod/pull/3261))

### Deprecated
- Deprecated ElasticRayExecutor APIs in favor of the new RayExecutor API for issue: [#3190](https://github.com/horovod/horovod/issues/3190).
### Removed

### Fixed

- fix the example of pytorch_lightning_mnist.py ([#3245](https://github.com/horovod/horovod/pull/3245))

- Call _setup in remote trainers to point to the correct shared lib path ([#3258](https://github.com/horovod/horovod/pull/3258))
## [v0.23.0] - 2021-10-06

### Added

- Added process sets to concurrently run collective operations on subsets of Horovod processes in TensorFlow, PyTorch, and MXNet. ([#2839](https://github.com/horovod/horovod/pull/2839), [#3042](https://github.com/horovod/horovod/pull/3042), [#3043](https://github.com/horovod/horovod/pull/3043), [#3054](https://github.com/horovod/horovod/pull/3054), [#3083](https://github.com/horovod/horovod/pull/3083), [#3090](https://github.com/horovod/horovod/pull/3090))

- Added XLA support for Allreduce via `tf.function(jit_compile=True)`. ([#3053](https://github.com/horovod/horovod/pull/3053))

- Added fused buffer scaling and unpack/pack kernels on GPU. ([#2973](https://github.com/horovod/horovod/pull/2973))

- Added support for NCCL on CUDA 11.4. ([#3182](https://github.com/horovod/horovod/issues/3182))

- Added fp16 compression for MXNet. ([#2987](https://github.com/horovod/horovod/issues/2987))

- Added terminate_on_nan flag to Spark Lightning estimator. ([#3088](https://github.com/horovod/horovod/issues/3088))

- Added barrier() API to torch module to support simple synchronization among ranks and to achieve parity with PyTorch DDP and similar frameworks. [#3139](https://github.com/horovod/horovod/pull/3139)

- Added params for customizing Tensorboard callback. ([#3153](https://github.com/horovod/horovod/issues/3153))

- Added `hvd.cross_rank()` for keras. ([#3008](https://github.com/horovod/horovod/issues/3008))

- Added barrier() API to torch module to support simple synchronization among ranks and to achieve parity with PyTorch DDP and similar frameworks. [#3139](https://github.com/horovod/horovod/pull/3139)

### Changed

- Implemented more asynchronous dependency handling on GPU. ([#2963](https://github.com/horovod/horovod/pull/2963))

- Ray: RayExecutor will now use the current placement group instead of always creating a new one. ([#3134](https://github.com/horovod/horovod/pull/3134))

- Lightning: turned off shuffling for validation dataset. ([#2974](https://github.com/horovod/horovod/pull/2974))

- Ray: RayExecutor will use the current placement group if one exists. ([#3134](https://github.com/horovod/horovod/pull/3134))

- Extended `hvd.join()` to return the last rank that joined. ([#3097](https://github.com/horovod/horovod/pull/3097)
### Deprecated

### Removed

- Spark/Keras: remove bare Keras support. ([#3191](https://github.com/horovod/horovod/pull/3191))

### Fixed

- Fix Horovod develop/editable install mode and incremental builds. ([#3074](https://github.com/horovod/horovod/pull/3074))

- Estimator/Lightning: use lightning datamodule. ([#3084](https://github.com/horovod/horovod/pull/3084))

- Fix Horovod Spark StringType and numpy type mapping issue. ([#3146](https://github.com/horovod/horovod/pull/3146))

- Fixed error in Keras LearningRateScheduler. ([#3135](https://github.com/horovod/horovod/pull/3135))

- Fixed bug in Lightning Profiler on Ray. ([#3122](https://github.com/horovod/horovod/pull/3122))

- Fixed torch op lazy release to prevent OOM in elastic training. ([#3110](https://github.com/horovod/horovod/pull/3110))

- Lightning: Fixed usage of the checkpoint callback. ([#3186](https://github.com/horovod/horovod/pull/3186))

- Fixed MPICH support to use Intel MPI's implementation. ([#3148](https://github.com/horovod/horovod/pull/3148))


- Fixed race condition in PyTorch async dataloader. ([#3120](https://github.com/horovod/horovod/pull/3120))

- Keras: Fixed learning rate scheduler. ([#3142](https://github.com/horovod/horovod/pull/3142), [#3135](https://github.com/horovod/horovod/pull/3135))

## [v0.22.1] - 2021-06-10

### Added

- Estimator: added support for loading data from S3, GCS, ADLS, and other remote filesystems. ([#2927](https://github.com/horovod/horovod/issues/2927))

- Estimator: added custom Spark data loader interface. ([#2938](https://github.com/horovod/horovod/issues/2923))

- LightningEstimator: added support to supply a logger and associated parameter to control the frequency of logging. ([#2926](https://github.com/horovod/horovod/pull/2926))

- Estimator: added check to ensure all ranks have the same device type. ([#2942](https://github.com/horovod/horovod/pull/2942))

### Changed

- Changed behavior from using TensorBoardLogger to now using it as a fallback if a logger is not supplied. ([#2926](https://github.com/horovod/horovod/pull/2926))

- Ray: disabled capturing child tasks in placement group. ([#2920](https://github.com/horovod/horovod/pull/2920))

### Fixed

- Fixed `hvd.tensorflow.keras.Compression`, accidentally removed in v0.22.0. ([#2945](https://github.com/horovod/horovod/pull/2945))

- TorchEstimator: fixed usage of `validation_steps` in place of `validation_steps_per_epoch`. ([#2918](https://github.com/horovod/horovod/pull/2918))

- TensorFlow: fixed C++ API for TF v2.6.0. ([#2932](https://github.com/horovod/horovod/pull/2932))

- PyTorch: fixed `sparse_allreduce_async` for PyTorch v0.10.0. ([#2965](https://github.com/horovod/horovod/pull/2965))

## [v0.22.0] - 2021-05-18

### Added

- Added pytorch_lightning spark estimator which enables training pytorch_lightning models. ([#2713](https://github.com/horovod/horovod/pull/2713))

- Added NVTX tracing hooks for profiling with Nsight Systems. ([#2723](https://github.com/horovod/horovod/pull/2723))

- Added a generic `num_workers` API for ``RayExecutor`` ([#2870](https://github.com/horovod/horovod/pull/2870))

- Supports Ray Client without code changes. ([#2882](https://github.com/horovod/horovod/pull/2882))

- Supports inmemory cache option for Keras Estimator. ([#2896](https://github.com/horovod/horovod/pull/2896))

- Added FP16 support for GPU tensor in mxnet. ([#2915](https://github.com/horovod/horovod/pull/2915))

- Added response caching for allgather operations. ([#2872](https://github.com/horovod/horovod/pull/2872))

- Estimator: add petastorm reader_pool_type into constructor ([#2903](https://github.com/horovod/horovod/pull/2903))

### Changed

- Changed `alltoall` to return the received splits as a second return value if non-uniform splits are sent. ([#2631](https://github.com/horovod/horovod/pull/2631))

- Changed ``RayExecutor`` to use [Ray Placement Groups](https://docs.ray.io/en/master/placement-group.html) for worker colocation. ([#2824](https://github.com/horovod/horovod/pull/2824))

- Changed ``Inmemory dataloader`` usage for Torch Estimator with petastorm v0.11.0 release. ([#2896](https://github.com/horovod/horovod/pull/2896))

### Fixed

- Changed RayExecutor to use Ray node ID to enable multi-container:single-host setups. ([#2883](https://github.com/horovod/horovod/pull/2882))

- Support sparse gradients aggregation in TF1 Keras. ([#2879](https://github.com/horovod/horovod/pull/2879))

- Respect `global_step` parameter for LegacyOptimizers when aggregating gradients.  ([#2879](https://github.com/horovod/horovod/pull/2879))

- Fixed compatibility with PyTorch 1.9.0. ([#2829](https://github.com/horovod/horovod/pull/2829))

## [v0.21.3] - 2021-02-15

### Added

- Add `groups` parameter in `DistributedOptimizer` for custom allreduce groups. ([#2523](https://github.com/horovod/horovod/pull/2523))

### Removed

- Removed `num_groups` parameter in `DistributedOptimizer`, replaced with `groups`. ([#2523](https://github.com/horovod/horovod/pull/2523))

### Fixed

- Fixed worker desynchronization deadlock issue in TensorFlow 2.4. ([#2647](https://github.com/horovod/horovod/pull/2647))

- Deduped Keras `LearningRateWarmupCallback` log after gradual learning rate warmup. ([#2661](https://github.com/horovod/horovod/pull/2661))

## [v0.21.2] - 2021-02-08

### Added

- Added support for Intel(R) MPI in horovodrun. ([#2374](https://github.com/horovod/horovod/pull/2374))

- Add support for callbacks in Ray Elastic Executor. ([#2639](https://github.com/horovod/horovod/pull/2639))

- Added forwarding of stdout/stderr captured to driver over Gloo. ([#2646](https://github.com/horovod/horovod/pull/2646))

### Fixed

- Fixed broadcast_optimizer_state to handle NoneType params for PyTorch 1.8. ([#2624](https://github.com/horovod/horovod/pull/2624))

- Fixed `local_rank` support for Ray. ([#2596](https://github.com/horovod/horovod/pull/2596))

- Fixed DL estimators to obtain the output df schema without sampling the input. ([#2611](https://github.com/horovod/horovod/pull/2611))

- Fixed wrong default for horovod.tensorflow.keras.allreduce average ([#2627](https://github.com/horovod/horovod/pull/2627))

## [v0.21.1] - 2021-01-06

### Added

- Added in-memory dataset caching param to `TorchEstimator`. ([#2434](https://github.com/horovod/horovod/pull/2434))

- Added `val_batch_size` param to the Estimator API. ([#2505](https://github.com/horovod/horovod/pull/2505))

- Added support for TorchScript modules when using `TorchEstimator`. ([#2494](https://github.com/horovod/horovod/pull/2494))

### Changed

- Migrated to oneCCL aligned with oneAPI specification v1.0. ([#2513](https://github.com/horovod/horovod/pull/2513))

- Added knob to set cache hint for oneCCL allreduce. ([#2560](https://github.com/horovod/horovod/pull/2560))

- Renamed `horovodrun` arg `--ccl-bgt-affinity` to `--thread-affinity`. ([#2562](https://github.com/horovod/horovod/pull/2562))

- Changed default build parallelism from `-j8` to `-j1` to address potential race condition. ([#2572](https://github.com/horovod/horovod/pull/2572))

### Fixed

- Fixed building Horovod for ROCm PyTorch with newer hipify script. ([#2360](https://github.com/horovod/horovod/pull/2360))

- Fixed "Executable class" support for Ray. ([#2510](https://github.com/horovod/horovod/pull/2510))

- Fixed TorchEstimator returning model without switching to eval mode. ([#2517](https://github.com/horovod/horovod/pull/2517))

- Remove ssh reliance for Ray elastic training. ([#2528](https://github.com/horovod/horovod/pull/2528))

- Fixed error handling for changing framework without reinstalling horovod. ([#2529](https://github.com/horovod/horovod/pull/2529))

- Fixed "Intermediate path does not exist" error with DBFSLocalStore. ([#2526](https://github.com/horovod/horovod/pull/2526))

- Avoid synchronization if workers are only shrinked in elastic mode. ([#2514](https://github.com/horovod/horovod/pull/2514))

- Fixed Ray resource test. ([#2575](https://github.com/horovod/horovod/pull/2575))

- Fixed usage of env variable `HOROVOD_GLOO_TIMEOUT_SECONDS` with `horovodrun`. ([#2571](https://github.com/horovod/horovod/pull/2571))

## [v0.21.0] - 2020-11-23

### Added

- Added support for backward_passes_per_step > 1 for TF Keras graph mode. ([#2346](https://github.com/horovod/horovod/pull/2346))

- Added support for backward_passes_per_step > 1 for TF Keras eager execution. ([#2371](https://github.com/horovod/horovod/pull/2371))

- Added support for backward_passes_per_step > 1 for TF LegacyOptimizer in graph mode. ([#2401](https://github.com/horovod/horovod/pull/2401))

- Added grouped allreduce to enable more efficient tensor fusion and deterministic training. ([#2453](https://github.com/horovod/horovod/pull/2453))

- Add support for specifying `op` and `compression` in `horovod.tensorflow.keras.allreduce()`. ([#2423](https://github.com/horovod/horovod/pull/2423))

- Adding support for batched D2D memcopy kernel on GPU. ([#2435](https://github.com/horovod/horovod/pull/2435))

- Added schema inference in Spark Estimator without sampling. ([#2373](https://github.com/horovod/horovod/pull/2373))

- Added `Store.create("dbfs:/")` mapping to `DBFSLocalStore("/dbfs/...")`. ([#2376](https://github.com/horovod/horovod/pull/2376))

### Changed

- Changed Keras callbacks to require parameter `initial_lr` of `LearningRateScheduleCallback` and `LearningRateWarmupCallback`. ([#2459](https://github.com/horovod/horovod/pull/2459))

- Changed default cycle time from 5ms to 1ms and fusion threshold from 64MB to 128MB. ([#2468](https://github.com/horovod/horovod/pull/2468))

### Fixed

- Fixed support for TensorFlow v2.4.0. ([#2381](https://github.com/horovod/horovod/pull/2381))

- Fixed averaging using CUDA half2 implementation one element half buffers. ([#2375](https://github.com/horovod/horovod/pull/2375))

- Fixed `HOROVOD_THREAD_AFFINITY` when using oneCCL. ([#2350](https://github.com/horovod/horovod/pull/2350))

- Added timeout to SSH check in horovodrun to prevent hanging. ([#2448](https://github.com/horovod/horovod/pull/2448))

- Added `HOROVOD_GLOO_TIMEOUT_SECONDS` value to error messages. ([#2436](https://github.com/horovod/horovod/pull/2436))

- Fixed race condition in dynamic timeline API. ([#2341](https://github.com/horovod/horovod/pull/2341))

- Fixed --log-hide-timestamp to apply to driver logs with Gloo. ([#2388](https://github.com/horovod/horovod/pull/2388))

- Fixed the search order of Eigen and Flatbuffers paths. ([#2473](https://github.com/horovod/horovod/pull/2473))

- Fixed type checks in `TorchEstimator` to correctly use `isinstance()`. ([#2480](https://github.com/horovod/horovod/pull/2480))

## [0.20.3] - 2020-10-01

### Added

- Added Elastic Ray integration. ([#2291](https://github.com/horovod/horovod/pull/2291))

### Changed

- Removed dependency on SSH access for Ray. ([#2275](https://github.com/horovod/horovod/pull/2275))

## [0.20.2] - 2020-09-25

### Fixed

- Fixed building Horovod without HOROVOD_WITHOUT_MXNET when MXNet is not installed. ([#2334](https://github.com/horovod/horovod/pull/2334))

## [0.20.1] - 2020-09-25

### Added

- Added Databricks storage `DBFSLocalStore` and support for GPU-aware scheduling to horovod.spark Estimator. ([#2234](https://github.com/horovod/horovod/pull/2234))

- Added ElasticSampler and PyTorch Elastic ImageNet example. ([#2297](https://github.com/horovod/horovod/pull/2297))

- Added ability to dynamically start and stop timeline programmatically. ([#2215](https://github.com/horovod/horovod/pull/2215))

- Added support for Gloo on macOS. ([#2254](https://github.com/horovod/horovod/pull/2254))

- Exposed name argument to TensorFlow allreduce operation. ([#2325](https://github.com/horovod/horovod/pull/2325))

- Added option to strip outer name scope from Horovod ops in TensorFlow. ([#2328](https://github.com/horovod/horovod/pull/2328))

### Fixed

- Fixed usage of VERBOSE=1 when setting custom MAKEFLAGS. ([#2239](https://github.com/horovod/horovod/pull/2239))

- Fixed bugs in Keras Elastic Callback classes. ([#2289](https://github.com/horovod/horovod/pull/2289))

- Fixed RelWithDebInfo build and made it the default with -03 optimizations. ([#2305](https://github.com/horovod/horovod/pull/2305))

- Fixed usage of tf.cond in TensorFlow alltoall gradient. ([#2327](https://github.com/horovod/horovod/pull/2327))

- Fixed allreduce averaging for TF IndexedSlices in ROCm path. ([#2279](https://github.com/horovod/horovod/pull/2279))

- Include stdexcept to handle certain compiler / frameworks that don't include it already. ([#2238](https://github.com/horovod/horovod/pull/2238))

- Fixed Debug builds by setting compiler options based on CMake build type. ([#2263](https://github.com/horovod/horovod/pull/2263))

- Skipped launching zero-sized send/recvs for NCCLAlltoall. ([#2273](https://github.com/horovod/horovod/pull/2273))

- Fixed missing run in tf keras elastic mode. ([#2272](https://github.com/horovod/horovod/pull/2272))

- Fixed loss function in TensorFlow2 elastic synthetic benchmark. ([#2265](https://github.com/horovod/horovod/pull/2265))

- Fixed usage of HOROVOD_MIXED_INSTALL env var in alltoall tests. ([#2266](https://github.com/horovod/horovod/pull/2266))

- Removed keras requirement from Ray example. ([#2262](https://github.com/horovod/horovod/pull/2262))

## [0.20.0] - 2020-09-02

### Added

- Added bare-metal elastic mode implementation to enable auto-scaling and fault tolerance. ([#1849](https://github.com/horovod/horovod/pull/1849))

- Added Elastic Horovod support for Spark auto-scaling. ([#1956](https://github.com/horovod/horovod/pull/1956))

- Added All-to-All operation for TensorFlow, PyTorch, and MXNet. ([#2143](https://github.com/horovod/horovod/pull/2143))

- Added support for `gradient_predivide_factor` and averaging in Horovod backend. ([#1949](https://github.com/horovod/horovod/pull/1949))

- Added NCCL implementation of the allgather operation. ([#1952](https://github.com/horovod/horovod/pull/1952))

- Added `HOROVOD_GPU_OPERATIONS` installation variable to simplify enabling NCCL support for all GPU operations. ([#1960](https://github.com/horovod/horovod/pull/1960))

- Added TensorFlow implementation of `SyncBatchNormalization` layer. ([#2075](https://github.com/horovod/horovod/pull/2075))

- Added `hvd.is_initialized()` method. ([#2020](https://github.com/horovod/horovod/pull/2020))

- Added `hvd.allgather_object` function for TensorFlow, PyTorch, and MXNet. ([#2166](https://github.com/horovod/horovod/pull/2166))

- Added `hvd.broadcast_object` function for MXNet. ([#2122](https://github.com/horovod/horovod/pull/2122))

- Added `label_shapes` parameter to KerasEstimator and TorchEstimator. ([#2140](https://github.com/horovod/horovod/pull/2140))

- Added optional `modelCheckPoint` callback to KerasEstimator params. ([#2124](https://github.com/horovod/horovod/pull/2124))

- Added `ssh_identity_file` argument to `horovodrun`. ([#2201](https://github.com/horovod/horovod/pull/2201))

- Added support for `horovodrun` on `kubeflow/mpi-job`. ([#2199](https://github.com/horovod/horovod/pull/2199))

- Added Ray integration. ([#2218](https://github.com/horovod/horovod/pull/2218))

### Changed

- Moved `horovod.run.runner.run` to `horovod.run`. ([#2099](https://github.com/horovod/horovod/pull/2099))

- HOROVOD_THREAD_AFFINITY accepts multiple values, one for every Horovod rank ([#2131](https://github.com/horovod/horovod/pull/2131))

- Migrated build system for native libraries to CMake ([#2009](https://github.com/horovod/horovod/pull/2009))

### Deprecated

- HOROVOD_CCL_BGT_AFFINITY is deprected. Use HOROVOD_THREAD_AFFINITY instead ([#2131](https://github.com/horovod/horovod/pull/2131))

### Removed

- Dropped support for Python 2. ([#1954](https://github.com/horovod/horovod/pull/1954))

- Dropped support for TensorFlow < 1.15. ([#2169](https://github.com/horovod/horovod/pull/2169))

- Dropped support for PyTorch < 1.2. ([#2086](https://github.com/horovod/horovod/pull/2086))

### Fixed

- Fixed MXNet allgather implementation to correctly handle resizing the output buffer. ([#2092](https://github.com/horovod/horovod/pull/2092))

- Fixed Keras Spark Estimator incompatibility with TensorFlow 1.15 due to `tf.autograph`. ([#2069](https://github.com/horovod/horovod/pull/2069))

- Fixed API compatibility with PyTorch 1.6. ([#2051](https://github.com/horovod/horovod/pull/2051))

- Fixed Keras API compatibility with TensorFlow 2.4.0. ([#2178](https://github.com/horovod/horovod/pull/2178))

- Fixed allgather gradient for TensorFlow 2 in cases where the tensor shape is not known during graph construction. ([#2121](https://github.com/horovod/horovod/pull/2121))

- Fixed running using Gloo with an imbalanced number of workers per host. ([#2212](https://github.com/horovod/horovod/pull/2212))
## Governance of the Horovod Project

Horovod is a graduated projected within the [LF AI & Data Foundation](https://lfaidata.foundation/).

### Charter

You can find Horovod Charter [here](https://wiki.lfai.foundation/download/attachments/7733301/Horovod%20Project%20Technical%20Charter%2012-22-2018%20FINAL.pdf?version=1&modificationDate=1558389484000&api=v2)

### Technical Steering Committee

Horovod development is governed by the Horovod Technical Steering Committee (TSC). The TSC consists of voting and
non-voting members, in addition to a chairman responsible for running TSC meetings, setting the meeting agenda, and
calling votes on proposals.

Current chairman of the Horovod TSC:
* [Travis Addair](https://github.com/tgaddair) - Predibase

Current voting members of the Horovod TSC:
* [Alex Sergeev](https://github.com/alsrgv) - Carbon Robotics
* [Travis Addair](https://github.com/tgaddair) - Predibase
* [Can Karakus](https://github.com/karakusc) - Amazon
* [Josh Romero](https://github.com/romerojosh) - NVIDIA
* [Nicolas Castet](https://github.com/nvcastet) - NVIDIA
* [Enrico Minack](https://github.com/EnricoMi) - G-Research
* [Xu Ning](https://github.com/thuningxu) - Uber
* [Todd Mytkowicz](https://github.com/klipto) - Microsoft

Current non-voting members of the Horovod TSC:
* [Leonard Lausen](https://github.com/leezu) - Amazon
* [Jonathan Dekhtiar](https://github.com/DEKHTIARJonathan) - NVIDIA
* [Richard Liaw](https://github.com/richardliaw) - Anyscale
* [Neil Conway](https://github.com/neilconway) - Determined AI, HPE
* [Min Cai](https://github.com/mincai) - Uber
* [Chongxiao Cao](https://github.com/chongxiaoc) - Uber
* [Max Gerlach](https://github.com/maxhgerlach) - DeepL
* [Ryan Beethe](https://github.com/rb-determined-ai) - Determined AI, HPE
* [Abin Shahab](https://github.com/ashahab) - LinkedIn
* [TJ Xu](https://github.com/Tixxx) - Uber

Emeritus members of the Horovod TSC:
* [Lin Yuan](https://github.com/apeforest)
* [Haibin Lin](https://github.com/eric-haibin-lin)
* [Yuxi Hu](https://github.com/yuxihu)
* [Emad Barsoum](https://github.com/ebarsoum)
* [Aaron Harlap](https://github.com/aaron276h)
* [Jaliya Ekanayake](https://github.com/jaliyae)
* [Kaarthik Sivashanmugam](https://github.com/skaarthik)
* [Armand McQueen](https://github.com/armandmcqueen)

Non-voting members of the TSC ("maintainers") have commit access to the Horovod GitHub repository, and take part in the
standing TSC meetings and mailing lists. Emeritus members are no longer active maintainers of the project, but are
welcome to participate in any TSC meeting.

The Horovod TSC meets monthly and publishes meeting notes via a [mailing list](https://lists.lfai.foundation/g/horovod-tsc).
This mailing list can also be utilized to reach out to the TSC.  Major decisions regarding the technical directions of
the Horovod project will be brought before the TSC for discussion, with an accompanying proposal document termed an RFC
(Request for Comments).

Technical decisions made by the TSC should be unanimous, with each voting and non-voting member either agreeing to the
proposal or abstaining for it to pass.  If consensus cannot be reached, then the proposal is to be put to a vote
among the voting members of the TSC, at which point a majority of the voting TSC must agree to the proposal for it to
pass.

Decisions to add or change members of the TSC in either a voting or non-voting capacity are handled the same as other
proposals (without an RFC): an attempt is made to reach a unanimous decision among the entire TSC, followed by a vote
among the voting members if no consensus can be reached.
# Security Policy

## Reporting a Vulnerability

Please report security vulnerabilities to horovod-security@lists.lfaidata.foundation.

Anyone can post to this mailing list, however, only active maintainers of the Horovod project will be able to read to the message.
# Code of Conduct

Horovod is a project hosted under the LF AI Foundation. As such, we would like to urge you to please be mindful of and adhere to the Linux Foundation’s [Code of Conduct](https://lfprojects.org/policies/code-of-conduct/) when contributing to the Horovod.

If you have any questions or concerns, please email info@lfai.foundation.

Thank you. 
## Contributing to Horovod

**Thanks for taking the time to contribute!**

Refer to the following guidelines to contribute new functionality or bug fixes to Horovod:
1. Use [autopep8](https://github.com/hhatto/autopep8) to format the Python code.
2. Use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to format C++ code.
3. Add unit tests for any new code you write.
4. Run unit tests in both CPU and GPU environments.

### Code of Conduct

Please be mindful of and adhere to the Linux Foundation's
[Code of Conduct](https://lfprojects.org/policies/code-of-conduct) when contributing to Horovod.
## Checklist before submitting

- [ ] Did you read the [contributor guide](https://github.com/horovod/horovod/blob/master/CONTRIBUTING.md)?
- [ ] Did you update the docs?
- [ ] Did you write any tests to validate this change?  
- [ ] Did you update the [CHANGELOG](https://github.com/horovod/horovod/blob/master/CHANGELOG.md), if this change affects users?

## Description

Fixes # (issue).

## Review process to land 

1. All tests and other checks must succeed.
2. At least one member of the [technical steering committee](https://github.com/horovod/horovod/blob/master/CONTRIBUTING.md) must review and approve.
3. If any member of the technical steering committee requests changes, they must be addressed.
---
name: Bug report
about: For questions about using Horovod or getting it to work in you environment, use https://github.com/horovod/horovod/discussions
title: ''
labels: bug
assignees: ''

---

**Environment:**
1. Framework: (TensorFlow, Keras, PyTorch, MXNet)
2. Framework version:
3. Horovod version:
4. MPI version:
5. CUDA version:
6. NCCL version:
7. Python version:
8. Spark / PySpark version:
9. Ray version:
10. OS and version:
11. GCC version:
12. CMake version:

**Checklist:**
1. Did you search issues to find if somebody asked this question before?
2. If your question is about hang, did you read [this doc](https://github.com/horovod/horovod/blob/master/docs/running.rst)?
3. If your question is about docker, did you read [this doc](https://github.com/horovod/horovod/blob/master/docs/docker.rst)?
4. Did you check if you question is answered in the [troubleshooting guide](https://github.com/horovod/horovod/blob/master/docs/troubleshooting.rst)?

**Bug report:**
Please describe erroneous behavior you're observing and steps to reproduce it.

---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Horovod Docker Images

Often installing Horovod on bare metal can be difficult if your environment is not setup
correctly with CUDA, MPI, G++, CMake, etc. These Docker images are provided to simplify
the onboarding process for new users, and can serve as a starting point for building your
own runtime environment.

## Repositories

Separate images are provided for different Horovod configurations, and are published
to separate repos in DockerHub.

* `horovod/horovod` Horovod built with CUDA support and packaged with the latest stable TensorFlow, PyTorch, MXNet, 
  and Spark releases
* `horovod/horovod-cpu` Horovod built for CPU training and packaged with the latest stable TensorFlow, PyTorch, MXNet, 
  and Spark releases
* `horovod/horovod-ray` Horoovd built with CUDA support from the latest 
  [ray-project/ray:nightly-gpu](https://github.com/ray-project/ray) and packaged with the latest stable 
  TensorFlow and PyTorch releases

## Image Tags

* `master` - built from Horovod's `master` branch
* `nightly` - nightly build of Horovod
* `sha-<commit point>` - version of Horovod at designated git sha1 7-character commit point

## Building Custom Images

Build arguments are provided to allow the user to build Horovod against custom versions of various frameworks,
including:

* `TENSORFLOW_VERSION` - version of `tensorflow` pip package to install
* `PYTORCH_VERSION` - version of `torch` pip package to install
* `PYTORCH_LIGHTNING_VERSION` - version of `pytorch_lightning` pip package to install
* `TORCHVISION_VERSION` - version of `torchvision` pip package to install
* `MXNET_VERSION` - version of `mxnet` pip package to install
* `CUDNN_VERSION` - version of `libcudnn` apt package to install (only for `horovod` image)
* `NCCL_VERSION` - version of `libnccl` apt package to install (only for `horovod` image)
* `CUDA_DOCKER_VERSION` - tag of the `nvidia/cuda` image to build from (only for `horovod` image)
* `RAY_DOCKER_VERSION` - tag of the `rayproject/ray` GPU image to build from (only for `horovod-ray` image)

Building the Docker images should be run from the root Horovod directory. For example:

```
export DOCKER_BUILDKIT=1
docker build \
    --build-arg TENSORFLOW_VERSION=2.3.1 \
    --build-arg PYTORCH_VERSION=1.7.0+cu110 \
    -f docker/horovod/Dockerfile .
```

## Running Containers

See the [Horovod in Docker](../docs/docker.rst) documentation for guidance on running these Docker images, and
[Horovod on Ray](../docs/ray.rst) for usage with Ray.
.. raw:: html

    <p align="center"><img src="https://user-images.githubusercontent.com/16640218/34506318-84d0c06c-efe0-11e7-8831-0425772ed8f2.png" alt="Logo" width="200"/></p>
    <br/>

Horovod
=======

.. raw:: html

   <div align="center">

.. image:: https://badge.fury.io/py/horovod.svg
   :target: https://badge.fury.io/py/horovod
   :alt: PyPI Version

.. image:: https://badge.buildkite.com/6f976bc161c69d9960fc00de01b69deb6199b25680a09e5e26.svg?branch=master
   :target: https://buildkite.com/horovod/horovod
   :alt: Build Status

.. image:: https://readthedocs.org/projects/horovod/badge/?version=latest
   :target: https://horovod.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/badge/slack-chat-green.svg?logo=slack
   :target: https://forms.gle/cPGvty5hp31tGfg79
   :alt: Slack

.. raw:: html

   </div>

.. raw:: html

   <div align="center">

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :target: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :alt: License

.. image:: https://app.fossa.com/api/projects/git%2Bgithub.com%2Fhorovod%2Fhorovod.svg?type=shield
   :target: https://app.fossa.com/projects/git%2Bgithub.com%2Fhorovod%2Fhorovod?ref=badge_shield
   :alt: FOSSA Status

.. image:: https://bestpractices.coreinfrastructure.org/projects/2373/badge
   :target: https://bestpractices.coreinfrastructure.org/projects/2373
   :alt: CII Best Practices

.. image:: https://pepy.tech/badge/horovod
   :target: https://pepy.tech/project/horovod
   :alt: Downloads

.. raw:: html

   </div>

.. inclusion-marker-start-do-not-remove

|

Horovod is a distributed deep learning training framework for TensorFlow, Keras, PyTorch, and Apache MXNet.
The goal of Horovod is to make distributed deep learning fast and easy to use.


.. raw:: html

   <p><img src="https://raw.githubusercontent.com/lfai/artwork/master/lfaidata-assets/lfaidata-project-badge/graduate/color/lfaidata-project-badge-graduate-color.png" alt="LF AI & Data" width="200"/></p>


Horovod is hosted by the `LF AI & Data Foundation <https://lfdl.io>`_ (LF AI & Data). If you are a company that is deeply
committed to using open source technologies in artificial intelligence, machine, and deep learning, and want to support
the communities of open source projects in these domains, consider joining the LF AI & Data Foundation. For details
about who's involved and how Horovod plays a role, read the Linux Foundation `announcement <https://lfdl.io/press/2018/12/13/lf-deep-learning-welcomes-horovod-distributed-training-framework-as-newest-project/>`_.

|

.. contents::

|

Documentation
-------------

- `Latest Release <https://horovod.readthedocs.io/en/stable>`_
- `master <https://horovod.readthedocs.io/en/latest>`_

|

Why Horovod?
------------
The primary motivation for this project is to make it easy to take a single-GPU training script and successfully scale
it to train across many GPUs in parallel. This has two aspects:

1. How much modification does one have to make to a program to make it distributed, and how easy is it to run it?
2. How much faster would it run in distributed mode?

Internally at Uber we found the MPI model to be much more straightforward and require far less code changes than previous
solutions such as Distributed TensorFlow with parameter servers. Once a training script has been written for scale with
Horovod, it can run on a single-GPU, multiple-GPUs, or even multiple hosts without any further code changes.
See the `Usage <#usage>`__ section for more details.

In addition to being easy to use, Horovod is fast. Below is a chart representing the benchmark that was done on 128
servers with 4 Pascal GPUs each connected by RoCE-capable 25 Gbit/s network:

.. image:: https://user-images.githubusercontent.com/16640218/38965607-bf5c46ca-4332-11e8-895a-b9c137e86013.png
   :alt: 512-GPU Benchmark

Horovod achieves 90% scaling efficiency for both Inception V3 and ResNet-101, and 68% scaling efficiency for VGG-16.
See `Benchmarks <docs/benchmarks.rst>`_ to find out how to reproduce these numbers.

While installing MPI and NCCL itself may seem like an extra hassle, it only needs to be done once by the team dealing
with infrastructure, while everyone else in the company who builds the models can enjoy the simplicity of training them at
scale.


Install
-------
To install Horovod on Linux or macOS:

1. Install `CMake <https://cmake.org/install/>`__

.. raw:: html

    <p/>

2. If you've installed TensorFlow from `PyPI <https://pypi.org/project/tensorflow>`__, make sure that ``g++-5`` or above is installed.

   If you've installed PyTorch from `PyPI <https://pypi.org/project/torch>`__, make sure that ``g++-5`` or above is installed.

   If you've installed either package from `Conda <https://conda.io>`_, make sure that the ``gxx_linux-64`` Conda package is installed.

.. raw:: html

    <p/>

3. Install the ``horovod`` pip package.

   To run on CPUs:

   .. code-block:: bash

      $ pip install horovod

   To run on GPUs with NCCL:

   .. code-block:: bash

      $ HOROVOD_GPU_OPERATIONS=NCCL pip install horovod

For more details on installing Horovod with GPU support, read `Horovod on GPU <docs/gpus.rst>`_.

For the full list of Horovod installation options, read the `Installation Guide <docs/install.rst>`_.

If you want to use MPI, read `Horovod with MPI <docs/mpi.rst>`_.

If you want to use Conda, read `Building a Conda environment with GPU support for Horovod <docs/conda.rst>`_.

If you want to use Docker, read `Horovod in Docker <docs/docker.rst>`_.

To compile Horovod from source, follow the instructions in the `Contributor Guide <docs/contributors.rst>`_.


Concepts
--------
Horovod core principles are based on `MPI <http://mpi-forum.org/>`_ concepts such as *size*, *rank*,
*local rank*, **allreduce**, **allgather**, **broadcast**, and **alltoall**. See `this page <docs/concepts.rst>`_
for more details.

Supported frameworks
--------------------
See these pages for Horovod examples and best practices:

- `Horovod with TensorFlow <docs/tensorflow.rst>`_
- `Horovod with XLA in Tensorflow <xla.rst>`_
- `Horovod with Keras <docs/keras.rst>`_
- `Horovod with PyTorch <docs/pytorch.rst>`_
- `Horovod with MXNet <docs/mxnet.rst>`_

Usage
-----

To use Horovod, make the following additions to your program:

1. Run ``hvd.init()`` to initialize Horovod.

.. raw:: html

    <p/>

2. Pin each GPU to a single process to avoid resource contention.

   With the typical setup of one GPU per process, set this to *local rank*. The first process on
   the server will be allocated the first GPU, the second process will be allocated the second GPU, and so forth.

.. raw:: html

    <p/>


3. Scale the learning rate by the number of workers.

   Effective batch size in synchronous distributed training is scaled by the number of workers.
   An increase in learning rate compensates for the increased batch size.

.. raw:: html

    <p/>


4. Wrap the optimizer in ``hvd.DistributedOptimizer``.

   The distributed optimizer delegates gradient computation to the original optimizer, averages gradients using **allreduce** or **allgather**, and then applies those averaged gradients.

.. raw:: html

    <p/>


5. Broadcast the initial variable states from rank 0 to all other processes.

   This is necessary to ensure consistent initialization of all workers when training is started with random weights or restored from a checkpoint.

.. raw:: html

    <p/>


6. Modify your code to save checkpoints only on worker 0 to prevent other workers from corrupting them.

.. raw:: html

    <p/>


Example using TensorFlow v1 (see the `examples <https://github.com/horovod/horovod/blob/master/examples/>`_ directory for full training examples):

.. code-block:: python

    import tensorflow as tf
    import horovod.tensorflow as hvd


    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    config = tf.ConfigProto()
    config.gpu_options.visible_device_list = str(hvd.local_rank())

    # Build model...
    loss = ...
    opt = tf.train.AdagradOptimizer(0.01 * hvd.size())

    # Add Horovod Distributed Optimizer
    opt = hvd.DistributedOptimizer(opt)

    # Add hook to broadcast variables from rank 0 to all other processes during
    # initialization.
    hooks = [hvd.BroadcastGlobalVariablesHook(0)]

    # Make training operation
    train_op = opt.minimize(loss)

    # Save checkpoints only on worker 0 to prevent other workers from corrupting them.
    checkpoint_dir = '/tmp/train_logs' if hvd.rank() == 0 else None

    # The MonitoredTrainingSession takes care of session initialization,
    # restoring from a checkpoint, saving to a checkpoint, and closing when done
    # or an error occurs.
    with tf.train.MonitoredTrainingSession(checkpoint_dir=checkpoint_dir,
                                           config=config,
                                           hooks=hooks) as mon_sess:
      while not mon_sess.should_stop():
        # Perform synchronous training.
        mon_sess.run(train_op)


Running Horovod
---------------
The example commands below show how to run distributed training.
See `Run Horovod <docs/running.rst>`_ for more details, including RoCE/InfiniBand tweaks and tips for dealing with hangs.

1. To run on a machine with 4 GPUs:

   .. code-block:: bash

        $ horovodrun -np 4 -H localhost:4 python train.py

2. To run on 4 machines with 4 GPUs each:

   .. code-block:: bash

       $ horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 python train.py

3. To run using Open MPI without the ``horovodrun`` wrapper, see `Running Horovod with Open MPI <docs/mpi.rst>`_.

4. To run in Docker, see `Horovod in Docker <docs/docker.rst>`_.

5. To run on Kubernetes, see `Kubeflow MPI Operator <https://github.com/kubeflow/mpi-operator/>`_, `Helm Chart <https://github.com/kubernetes/charts/tree/master/stable/horovod/>`_, `FfDL <https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/>`_, and `Polyaxon <https://docs.polyaxon.com/integrations/horovod/>`_.

6. To run on Spark, see `Horovod on Spark <docs/spark.rst>`_.

7. To run on Ray, see `Horovod on Ray <docs/ray.rst>`_.

8. To run in Singularity, see `Singularity <https://github.com/sylabs/examples/tree/master/machinelearning/horovod>`_.

9. To run in a LSF HPC cluster (e.g. Summit), see `LSF <docs/lsf.rst>`_.

10. To run on Hadoop Yarn, see `TonY <https://github.com/linkedin/TonY/>`_.

Gloo
----
`Gloo <https://github.com/facebookincubator/gloo>`_ is an open source collective communications library developed by Facebook.

Gloo comes included with Horovod, and allows users to run Horovod without requiring MPI to be installed.

For environments that have support both MPI and Gloo, you can choose to use Gloo at runtime by passing the ``--gloo`` argument to ``horovodrun``:

.. code-block:: bash

     $ horovodrun --gloo -np 2 python train.py

mpi4py
------
Horovod supports mixing and matching Horovod collectives with other MPI libraries, such as `mpi4py <https://mpi4py.scipy.org>`_,
provided that the MPI was built with multi-threading support.

You can check for MPI multi-threading support by querying the ``hvd.mpi_threads_supported()`` function.

.. code-block:: python

    import horovod.tensorflow as hvd

    # Initialize Horovod
    hvd.init()

    # Verify that MPI multi-threading is supported.
    assert hvd.mpi_threads_supported()

    from mpi4py import MPI
    assert hvd.size() == MPI.COMM_WORLD.Get_size()

You can also initialize Horovod with an `mpi4py` sub-communicator, in which case each sub-communicator
will run an independent Horovod training.

.. code-block:: python

    from mpi4py import MPI
    import horovod.tensorflow as hvd

    # Split COMM_WORLD into subcommunicators
    subcomm = MPI.COMM_WORLD.Split(color=MPI.COMM_WORLD.rank % 2,
                                   key=MPI.COMM_WORLD.rank)

    # Initialize Horovod
    hvd.init(comm=subcomm)

    print('COMM_WORLD rank: %d, Horovod rank: %d' % (MPI.COMM_WORLD.rank, hvd.rank()))


Inference
---------
Learn how to optimize your model for inference and remove Horovod operations from the graph `here <docs/inference.rst>`_.


Tensor Fusion
-------------
One of the unique things about Horovod is its ability to interleave communication and computation coupled with the ability
to batch small **allreduce** operations, which results in improved performance. We call this batching feature Tensor Fusion.

See `here <docs/tensor-fusion.rst>`__ for full details and tweaking instructions.


Horovod Timeline
----------------
Horovod has the ability to record the timeline of its activity, called Horovod Timeline.

.. image:: https://user-images.githubusercontent.com/16640218/29735271-9e148da0-89ac-11e7-9ae0-11d7a099ac89.png
   :alt: Horovod Timeline

Use Horovod timeline to analyze Horovod performance.
See `here <docs/timeline.rst>`__ for full details and usage instructions.


Automated Performance Tuning
----------------------------
Selecting the right values to efficiently make use of Tensor Fusion and other advanced Horovod features can involve
a good amount of trial and error. We provide a system to automate this performance optimization process called
**autotuning**, which you can enable with a single command line argument to ``horovodrun``.

See `here <docs/autotune.rst>`__ for full details and usage instructions.


Horovod Process Sets
--------------------
Horovod allows you to concurrently run distinct collective operations in different groups of processes taking part in
one distributed training. Set up ``hvd.process_set`` objects to make use of this capability.

See `Process Sets <docs/process_set.rst>`__ for detailed instructions.


Guides
------
1. Run distributed training in Microsoft Azure using `Batch AI and Horovod <https://github.com/Azure/BatchAI/tree/master/recipes/Horovod>`_.
2. `Distributed model training using Horovod <https://spell.ml/blog/distributed-model-training-using-horovod-XvqEGRUAACgAa5th>`_.

Send us links to any user guides you want to publish on this site

Troubleshooting
---------------
See `Troubleshooting <docs/troubleshooting.rst>`_ and submit a `ticket <https://github.com/horovod/horovod/issues/new>`_
if you can't find an answer.


Citation
--------
Please cite Horovod in your publications if it helps your research:

::

    @article{sergeev2018horovod,
      Author = {Alexander Sergeev and Mike Del Balso},
      Journal = {arXiv preprint arXiv:1802.05799},
      Title = {Horovod: fast and easy distributed deep learning in {TensorFlow}},
      Year = {2018}
    }


Publications
------------
1. Sergeev, A., Del Balso, M. (2017) *Meet Horovod: Uber’s Open Source Distributed Deep Learning Framework for TensorFlow*.
Retrieved from `https://eng.uber.com/horovod/ <https://eng.uber.com/horovod/>`_

2. Sergeev, A. (2017) *Horovod - Distributed TensorFlow Made Easy*. Retrieved from
`https://www.slideshare.net/AlexanderSergeev4/horovod-distributed-tensorflow-made-easy <https://www.slideshare.net/AlexanderSergeev4/horovod-distributed-tensorflow-made-easy>`_

3. Sergeev, A., Del Balso, M. (2018) *Horovod: fast and easy distributed deep learning in TensorFlow*. Retrieved from
`arXiv:1802.05799 <https://arxiv.org/abs/1802.05799>`_


References
----------
The Horovod source code was based off the Baidu `tensorflow-allreduce <https://github.com/baidu-research/tensorflow-allreduce>`_
repository written by Andrew Gibiansky and Joel Hestness. Their original work is described in the article
`Bringing HPC Techniques to Deep Learning <http://andrew.gibiansky.com/blog/machine-learning/baidu-allreduce/>`_.

Getting Involved
----------------
- `Community Slack <https://forms.gle/cPGvty5hp31tGfg79>`_ for collaboration and discussion
- `Horovod Announce <https://lists.lfai.foundation/g/horovod-announce>`_ for updates on the project
- `Horovod Technical-Discuss <https://lists.lfai.foundation/g/horovod-technical-discuss>`_ for public discussion
- `Horovod Security <https://lists.lfai.foundation/g/horovod-security>`_ to report security vulnerabilities


.. inclusion-marker-end-do-not-remove
   Place contents above here if they should also appear in read-the-docs.
   Contents below are already part of the read-the-docs table of contents.

.. inclusion-marker-start-do-not-remove


Concepts
========

Horovod core principles are based on the `MPI <http://mpi-forum.org/>`_ concepts *size*, *rank*,
*local rank*, *allreduce*, *allgather*, *broadcast*, and *alltoall*. These are best explained by example. Say we launched
a training script on 4 servers, each having 4 GPUs. If we launched one copy of the script per GPU:

* *Size* would be the number of processes, in this case, 16.

* *Rank* would be the unique process ID from 0 to 15 (*size* - 1).

* *Local rank* would be the unique process ID within the server from 0 to 3.

* *Allreduce* is an operation that aggregates data among multiple processes and distributes results back to them.  *Allreduce* is used to average dense tensors.  Here's an illustration from the `MPI Tutorial <http://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/>`__:

.. image:: http://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/mpi_allreduce_1.png
   :alt: Allreduce Illustration

* *Allgather* is an operation that gathers data from all processes on every process.  *Allgather* is used to collect values of sparse tensors.  Here's an illustration from the `MPI Tutorial <http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/>`__:

.. image:: http://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/allgather.png
   :alt: Allgather Illustration


* *Broadcast* is an operation that broadcasts data from one process, identified by root rank, onto every other process. Here's an illustration from the `MPI Tutorial <http://mpitutorial.com/tutorials/mpi-broadcast-and-collective-communication/>`__:

    .. image:: http://mpitutorial.com/tutorials/mpi-broadcast-and-collective-communication/broadcast_pattern.png
       :alt: Broadcast Illustration


* *Alltoall* is an operation to exchange data between all processes.  *Alltoall* may be useful to implement neural networks with advanced architectures that span multiple devices.


.. inclusion-marker-end-do-not-remove
.. include:: ./hyperparameter_search.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./mpi.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove.. include:: ./install.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove.. inclusion-marker-start-do-not-remove

Autotune: Automated Performance Tuning
======================================

Horovod comes with several adjustable "knobs" that can affect runtime performance, including
``--fusion-threshold-mb`` and ``--cycle-time-ms`` (tensor fusion), ``--cache-capacity`` (response cache), and
hierarchical collective algorithms ``--hierarchical-allreduce`` and ``--hierarchical-allgather``.

Determining the best combination of these values to maximize performance (minimize time to convergence) can be a
matter of trial-and-error, as many factors including model complexity, network bandwidth, GPU memory, etc. can all
affect inputs per second throughput during training.

Horovod provides a mechanism to automate the process of selecting the best values for these "knobs" called
**autotuning**. The Horovod autotuning system uses
`Bayesian optimization <https://en.wikipedia.org/wiki/Bayesian_optimization>`_ to intelligently search through the
space of parameter combinations during training. This feature can be enabled by setting the ``--autotune`` flag for
``horovodrun``:

.. code-block:: bash

    $ horovodrun -np 4 --autotune python train.py

When autotuning is enabled, Horovod will spend the first steps / epochs of training experimenting with different
parameter values and collecting metrics on performance (measured in bytes allreduced / allgathered per unit of time).
Once the experiment reaches convergence, or a set number of samples have been collected, the system will record the best
combination of parameters discovered and continue to use them for the duration of training.

A log of all parameter combinations explored (and the best values selected) can be recorded by providing
the ``--autotune-log-file`` option to ``horovodrun``:

.. code-block:: bash

    $ horovodrun -np 4 --autotune --autotune-log-file /tmp/autotune_log.csv python train.py

By logging the best parameters to a file, you can opt to set the best parameters discovered on the command line
instead of re-running autotuning if training is paused and later resumed.

Note that some configurable parameters, like tensor compression, are not included as part of the autotuning process
because they can affect model convergence. The purpose of autotuning at this time is entirely to improve scaling
efficiency without making any tradeoffs on model performance.


Constant Parameters
-------------------

Sometimes you may wish to hold certain values constant and only tune the unspecified parameters. This can be
accomplished by explicitly setting those values on the command line or in the config file provided
by ``--config-file``:

.. code-block:: bash

    $ horovodrun -np 4 --autotune --cache-capacity 1024 --no-hierarchical-allgather python train.py

In the above example, parameters ``cache-capacity`` and ``hierarchical-allgather`` will not be adjusted by
autotuning.


Advanced Autotuning
-------------------

Enabling autotuning imposes a tradeoff between degraded performance during the early phases of training in exchange for
better performance later on. As such, it's generally recommended to use autotuning in situations where training is both
expected to take a long time (many epochs on very large datasets) and where scaling efficiency has been found lacking
using the default settings.

You can tune the autotuning system itself to change the number of warmup samples (discarded samples at the beginning),
steps per sample, and maximum samples:

.. code-block:: bash

    $ horovodrun -np 4 --autotune \
    --autotune-warmup-samples 5 --autotune-steps-per-sample 20 --autotune-bayes-opt-max-samples 40 \
    python train.py

Increasing these values will generally improve the accuracy of the autotuning process at the cost of greater time
spent in the autotuning process with degraded performance.

Finally, for those familiar with the underlying theory of Bayesian optimization and Gaussian processes, you can tune
the noise regularization term (alpha) to account for variance in your network or other system resources:

.. code-block:: bash

    $ horovodrun -np 4 --autotune --autotune-gaussian-process-noise 0.75 python train.py

.. inclusion-marker-end-do-not-remove
Horovod with MXNet
==================
Horovod supports Apache MXNet and regular TensorFlow in similar ways.

See full training `MNIST <https://github.com/horovod/horovod/blob/master/examples/mxnet/mxnet_mnist.py>`__ and `ImageNet <https://github.com/horovod/horovod/blob/master/examples/mxnet/mxnet_imagenet_resnet50.py>`__ examples.
The script below provides a simple skeleton of code block based on the Apache MXNet Gluon API.

.. code-block:: python

    import mxnet as mx
    import horovod.mxnet as hvd
    from mxnet import autograd

    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank
    context = mx.gpu(hvd.local_rank())
    num_workers = hvd.size()

    # Build model
    model = ...
    model.hybridize()

    # Create optimizer
    optimizer_params = ...
    opt = mx.optimizer.create('sgd', **optimizer_params)

    # Initialize parameters
    model.initialize(initializer, ctx=context)

    # Fetch and broadcast parameters
    params = model.collect_params()
    if params is not None:
        hvd.broadcast_parameters(params, root_rank=0)

    # Create DistributedTrainer, a subclass of gluon.Trainer
    trainer = hvd.DistributedTrainer(params, opt)

    # Create loss function
    loss_fn = ...

    # Train model
    for epoch in range(num_epoch):
        train_data.reset()
        for nbatch, batch in enumerate(train_data, start=1):
            data = batch.data[0].as_in_context(context)
            label = batch.label[0].as_in_context(context)
            with autograd.record():
                output = model(data.astype(dtype, copy=False))
                loss = loss_fn(output, label)
            loss.backward()
            trainer.step(batch_size)



.. NOTE:: Some MXNet versions do not work with Horovod:

    - MXNet 1.4.0 and earlier have `GCC incompatibility issues <https://github.com/horovod/horovod/issues/884>`__. Use MXNet 1.4.1 or later with Horovod 0.16.2 or later to avoid these incompatibilities.
    - MXNet 1.5.1, 1.6.0, 1.7.0, and 1.7.0.post1 are missing MKLDNN headers, so they do not work with Horovod. Use 1.5.1.post0, 1.6.0.post0, and 1.7.0.post0, respectively.
    - MXNet 1.6.0.post0 and 1.7.0.post0 are only available as mxnet-cu101 and mxnet-cu102.
.. include:: ./tensor-fusion.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Distributed Hyperparameter Search
=================================

Horovod's data parallelism training capabilities allow you to scale out and speed up the workload of training a deep learning model. However, simply using 2x more workers does not necessarily mean the model will obtain the same accuracy in 2x less time.

To address this, you often need to re-tune hyperparameters when training at scale, as many hyperparameters exhibit different behaviors at larger scales.

Horovod offers a `Ray Tune`_ integration to enable parallel hyperparameter tuning with distributed training.

.. image:: media/tune.png
    :align: center
    :scale: 20%

`Ray Tune`_ is an industry standard tool for distributed hyperparameter tuning. `Ray Tune`_ includes the latest hyperparameter search algorithms, integrates with TensorBoard and other analysis libraries, and natively supports distributed training. The `Ray Tune`_ + Horovod integration leverages the underlying Ray framework to provide a scalable and comprehensive hyperparameter tuning setup.

**By the end of this guide, you will learn:**

* How to set up `Ray Tune`_ and Horovod to tune your hyperparameters
* Typical hyperparameters to configure for distributed training

Horovod + Ray Tune
------------------

Leverage `Ray Tune`_ with Horovod to combine distributed hyperparameter tuning with distributed training. Here is an example demonstrating basic usage:

.. code-block:: python

    import horovod.torch as hvd
    from ray import tune
    import time

    def training_function(config: Dict):
        hvd.init()
        for i in range(config["epochs"]):
            time.sleep(1)
            model = Model(learning_rate=config["lr"])
            tune.report(test=1, rank=hvd.rank())

    trainable = DistributedTrainableCreator(
            training_function, num_slots=2, use_gpu=use_gpu)
    analysis = tune.run(
        trainable,
        num_samples=2,
        config={
            "epochs": tune.grid_search([1, 2, 3]),
            "lr": tune.grid_search([0.1, 0.2, 0.3]),
        }
    )
    print(analysis.best_config)

Basic setup
-----------

Use Ray Tune's `DistributedTrainableCreator`_ function to adapt your Horovod training function to be compatible with Ray Tune.

`DistributedTrainableCreator`_ exposes ``num_hosts``, ``num_slots``, ``use_gpu``, and ``num_cpus_per_slot``. Use these parameters to specify the resource allocation of a single "trial" (or "Trainable") which itself can be a distributed training job.


.. code-block:: python

    # Each training job will use 2 GPUs.
    trainable = DistributedTrainableCreator(
        training_function, num_slots=2, use_gpu=True)

The training function itself must do three things:

1. It must adhere to the `Tune Function API signature <https://docs.ray.io/en/latest/tune/api_docs/trainable.html#function-api>`__.
2. Its body must include a ``horovod.init()`` call.
3. It must call ``tune.report`` (`docs <https://docs.ray.io/en/latest/tune/api_docs/trainable.html#tune-report-tune-checkpoint-function-api>`__) during training, typically called iteratively at the end of every epoch.


Optimization of hyperparameters
-------------------------------

`Ray Tune`_ is able to orchestrate complex computational patterns with the `Ray Actor API <https://docs.ray.io/en/latest/actors.html>`__. For hyperparameter tuning, `Ray Tune`_ is able to conduct `parallel bayesian optimization <https://docs.ray.io/en/latest/tune/api_docs/suggestion.html>`__ and `Population Based Training <https://docs.ray.io/en/latest/tune/api_docs/schedulers.html>`__ on a group of distributed models.

You may need to implement model checkpointing. The rest of the optimization process can be configured with a couple lines of code.

.. code-block:: python

    from ray import tune
    from ray.tune.suggest.bayesopt import BayesOptSearch
    from ray.tune.suggest import ConcurrencyLimiter

    def training_function(config):
        ...

    algo = BayesOptSearch()
    algo = ConcurrencyLimiter(algo, max_concurrent=4)
    results = tune.run(
        training_function,
        config={"lr": tune.uniform(0.001, 0.1)},
        name="horovod",
        metric="mean_loss",
        mode="min",
        search_alg=algo)

    print(results.best_config)

**Search Space**

Tune has a native interface for `specifying search spaces <https://docs.ray.io/en/master/tune/api_docs/search_space.html#tune-search-space>`__. You can specify the search space via ``tune.run(config=...)``.

Thereby, either use the ``tune.grid_search`` primitive to specify an axis of a grid search...

.. code-block:: python

    tune.run(
        trainable,
        config={"bar": tune.grid_search([True, False])})


... or one of the random sampling primitives to specify distributions:

.. code-block:: python

    tune.run(
        trainable,
        config={
            "param1": tune.choice([True, False]),
            "bar": tune.uniform(0, 10),
            "alpha": tune.sample_from(lambda _: np.random.uniform(100) ** 2),
            "const": "hello"  # It is also ok to specify constant values.
        })

Read more about Tune's `Search Space API <https://docs.ray.io/en/master/tune/api_docs/search_space.html#tune-search-space>`__.

**Analyzing Results**

``tune.run`` returns an `Analysis <https://docs.ray.io/en/master/tune/api_docs/analysis.html>`__ object which has methods for analyzing your training.

.. code-block:: python

    analysis = tune.run(trainable, search_alg=algo, stop={"training_iteration": 20})

    best_trial = analysis.best_trial  # Get best trial
    best_config = analysis.best_config  # Get best trial's hyperparameters
    best_logdir = analysis.best_logdir  # Get best trial's logdir
    best_checkpoint = analysis.best_checkpoint  # Get best trial's best checkpoint
    best_result = analysis.best_result  # Get best trial's last results
    best_result_df = analysis.best_result_df  # Get best result as pandas dataframe


Set up a tuning cluster
-----------------------

Leverage `Ray Tune`_ with Horovod on a laptop, single machine with multiple GPUs, or across multiple machines. To run on a single machine, execute your Python script as-is (for example, `horovod_simple.py <https://docs.ray.io/en/latest/tune/examples/horovod_simple.html>`__, assuming Ray and Horovod are installed properly):

.. code-block:: bash

    python horovod_simple.py


To leverage a distributed hyperparameter tuning setup with `Ray Tune`_ + Horovod, install Ray and set up a `Ray cluster <https://docs.ray.io/en/latest/cluster/index.html>`__. Start a Ray cluster with the `Ray Cluster Launcher <https://docs.ray.io/en/latest/cluster/launcher.html>`__ or manually.

Below, we’ll use the `Ray Cluster Launcher <https://docs.ray.io/en/latest/cluster/launcher.html>`__, but you can start Ray on any list of nodes, on any cluster manager or cloud provider.

First, specify a configuration file. Below we have an example of using AWS EC2, but you can launch the cluster on any cloud provider:

.. code-block:: yaml

    # ray_cluster.yaml
    cluster_name: horovod-cluster
    provider: {type: aws, region: us-west-2}
    auth: {ssh_user: ubuntu}
    min_workers: 3
    max_workers: 3

    # Deep Learning AMI (Ubuntu) Version 21.0
    head_node: {InstanceType: p3.2xlarge, ImageId: ami-0b294f219d14e6a82}
    worker_nodes: {
        InstanceType: p3.2xlarge, ImageId: ami-0b294f219d14e6a82}
    setup_commands: # Set up each node.
        - HOROVOD_WITH_GLOO=1 HOROVOD_GPU_OPERATIONS=NCCL pip install horovod[ray]

Run ``ray up ray_cluster.yaml``, and a cluster of 4 nodes (1 head node + 3 worker nodes) will be automatically started with Ray.

.. code-block:: bash


    [6/6] Starting the Ray runtime
    Did not find any active Ray processes.
    Shared connection to 34.217.192.11 closed.
    Local node IP: 172.31.43.22
    2020-11-04 04:24:33,882 INFO services.py:1106 -- View the Ray dashboard at http://localhost:8265

    --------------------
    Ray runtime started.
    --------------------

    Next steps
      To connect to this Ray runtime from another node, run
        ray start --address='172.31.43.22:6379' --redis-password='5241590000000000'

      Alternatively, use the following Python code:
        import ray
        ray.init(address='auto', _redis_password='5241590000000000')

      If connection fails, check your firewall settings and network configuration.

      To terminate the Ray runtime, run
        ray stop
    Shared connection to 34.217.192.11 closed.
      New status: up-to-date

    Useful commands
      Monitor autoscaling with
        ray exec ~/dev/cfgs/check-autoscaler.yaml 'tail -n 100 -f /tmp/ray/session_latest/logs/monitor*'
      Connect to a terminal on the cluster head:
        ray attach ~/dev/cfgs/check-autoscaler.yaml
      Get a remote shell to the cluster manually:
        ssh -o IdentitiesOnly=yes -i ~/.ssh/ray-autoscaler_2_us-west-2.pem ubuntu@34.217.192.11

After the cluster is up, you can ssh into the head node and run your Tune script there.

Implementation (underneath the hood)
------------------------------------

Underneath the hood, `Ray Tune`_ will launch multiple "`trials <https://docs.ray.io/en/latest/tune/key-concepts.html#tune-run-and-trials>`__" in parallel. Each of these trials reference a `set of Ray actors <https://docs.ray.io/en/latest/actors.html>`__. For each trial, there will be 1 “coordinator actor,” and this coordinator actor will manage N training actors. One basic assumption of this implementation is that all sub-workers of a trial will be placed evenly across different machines.

.. image:: media/tune-horovod.jpg

Training actors will each hold a copy of the model and will create a communication group for Horovod allreduce. Training will execute on each actor, reporting intermediate metrics back to Tune.

This API requires Gloo as the underlying communication primitive. Be sure to install Horovod with ``HOROVOD_WITH_GLOO`` `enabled <https://horovod.readthedocs.io/en/stable/install_include.html#gloo>`__.


Common Hyperparameters
----------------------

We will cover a couple common hyperparameters that you may need to re-tune at scale:

1. Batch Size
2. Learning Rate schedules
3. Optimizers

Parameter: Batch size
~~~~~~~~~~~~~~~~~~~~~

By using data parallelism, it is necessary to scale the batch size along with workers to avoid reducing the per-worker workload and maximizing worker efficiency. However, increasing batch size can easily cause generalization issues (see this `Facebook Imagenet Training paper <https://research.fb.com/wp-content/uploads/2017/06/imagenet1kin1h5.pdf>`__ for more details).

**What are common solutions?**

* Linear scaling of learning rates: When the minibatch size is multiplied by k, multiply the learning rate by k.
* Dynamically adjusting batch size over the course of training:

  - One of the original papers presents a simple baseline of increasing the batch size over time
  - `ABSA provides a way <https://openreview.net/pdf?id=H1lnJ2Rqt7>`__ to leverage second order information to guide the batch size over time
  - `Gradient noise scale <https://openai.com/blog/science-of-ai/>`__ can be calculated to guide the increase of batch size over time

To leverage a dynamically changing batch size in training, you should either:

* Leverage `gradient accumulation <https://gist.github.com/thomwolf/ac7a7da6b1888c2eeac8ac8b9b05d3d3>`__
* Implement your own TrialScheduler to dynamically change the number of workers (coming soon)

Parameter: Learning rate schedules (warmup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As noted in this `Facebook Imagenet Training paper <https://research.fb.com/wp-content/uploads/2017/06/imagenet1kin1h5.pdf>`__, the linear scaling rule breaks down when the network is rapidly changing, which commonly occurs in early stages of training. This issue can be addressed with a "warmup," which is a strategy of using less aggressive learning rates at the start of training.

**What are common solutions?**

Goyal et al. (2017) proposes a warm-up schedule, where training usually starts with a small learning rate, and gradually increased to match a larger target learning rate. After the warm-up period (usually a few epochs), a regular learning rate schedule is used ("multi-steps", polynomial decay etc). Thus, there are generally three parameters for warmup schedule:

* Length of warmup (number of epochs)
* Starting learning rate
* Peak learning rate

Parameter: Optimizers
~~~~~~~~~~~~~~~~~~~~~

Optimizers are algorithms/methods that are used to update network weights iteratively. Common optimizers in deep learning include Adam, RMSProp, and SGD with momentum.

In large scale learning, naive approaches to optimizing and updating neural network weights can lead to poor generalization or decreased performance. For example,  Alexnet on Imagenet using standard SGD with momentum (and a warmup scheme) will stop scaling after ``B=2K``.

**What are common solutions?**

* `LARS <https://arxiv.org/pdf/1708.03888.pdf>`__ calculates a local learning rate per layer at each optimization step. It normalizes the gradient magnitude of each layer and instead uses a user-set coefficient and magnitude of the layer weights to compute the learning rate. The original paper for LARS presents performance improvements for training AlexNet with large batch sizes.
* `LAMB <https://towardsdatascience.com/an-intuitive-understanding-of-the-lamb-optimizer-46f8c0ae4866>`__ stands for “Layer-wise Adaptive Moments optimizer for Batch training.” It makes a few small changes to LARS. In spirit, it is “combining the ADAM optimizer with layer-wise scaling of LARS”. The original motivation of the LAMB work is because LARS did not work well for attention-based architectures such as BERT.

.. _`Ray Tune`: https://docs.ray.io/en/latest/tune/

.. _`DistributedTrainableCreator`: https://docs.ray.io/en/latest/tune/api_docs/integration.html#horovod-tune-integration-horovod

.. inclusion-marker-end-do-not-remove
.. include:: ./lsf.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove


Contributor Guide
=================

This guide covers the process of contributing to Horovod as a developer.


Environment Setup
-----------------

Clone the repository locally:

.. code-block:: bash

    $ git clone --recursive https://github.com/horovod/horovod.git

Develop within a virtual environment to avoid dependency issues:

.. code-block:: bash

    $ python3 -m venv env
    $ . env/bin/activate

We recommend installing package versions that match with those under test in
`Buildkite <https://github.com/horovod/horovod/blob/master/.buildkite/gen-pipeline.sh>`__.
The following versions are recommended (see default versions defined through :code:`ARG` in
`Dockerfile.test.cpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.cpu>`__ and
`Dockerfile.test.gpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.gpu>`__ file.

You can find all other non-Python packages that need to be installed on your system for Horovod to build
in the `Dockerfile.test.cpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.cpu>`__ and
`Dockerfile.test.gpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.gpu>`__ files.
Specifically, see all :code:`RUN apt-get install` lines.

Build and Install
-----------------

From *inside* the Horovod root directory, install Horovod in develop/editable mode:

.. code-block:: bash

    $ HOROVOD_WITH_PYTORCH=1 HOROVOD_WITH_TENSORFLOW=1 pip install -v -e .

Set ``HOROVOD_WITHOUT_[FRAMEWORK]=1`` to disable building Horovod plugins for that framework.
This is useful when you’re testing a feature of one framework in particular and wish to save time.

Set ``HOROVOD_WITH_[FRAMEWORK]=1`` to generate an error if the Horovod plugin for that framework failed to build.

Set ``HOROVOD_DEBUG=1`` for a debug build with checked assertions, disabled compiler optimizations etc.

Other environmental variables can be found in the `install documentation <https://github.com/horovod/horovod/blob/master/docs/install.rst#environment-variables>`__.

You can install optional dependencies defined in `setup.py <https://github.com/horovod/horovod/blob/master/setup.py>`__ by adding brackets
at the end of the command line e.g. ``[test]`` for test dependencies.
If you have not installed specific DL frameworks yet, add ``[dev]`` to install the CPU version of all supported DL frameworks.

In develop mode, you can edit the Horovod source directly in the repo folder. For Python code, the changes will take effect
immediately. For **C++/CUDA code**, the ``... pip install -v -e .`` command needs to be invoked again to perform an incremental build.

Testing
-------

Horovod has unit tests for all frameworks under ``test/parallel``. These should be invoked via ``horovodrun`` or
``mpirun`` and each test script may require to be run independently from the other test scripts:

.. code-block:: bash

    $ cd test/parallel
    $ horovodrun -np 2 pytest -v test_tensorflow.py
    $ horovodrun -np 2 pytest -v test_torch.py
    # ...

    # Or to run all framework tests:
    $ cd test/parallel
    $ ls -1 test_*.py | xargs -n 1 horovodrun -np 2 pytest -v

Moreover, there are integration tests and non-parallelized tests to be run directly via ``pytest``:

.. code-block:: bash

    $ cd test/integration
    $ pytest -v

    $ cd test/single
    $ pytest -v

**Note:** You will need PySpark and Java to run the Spark tests.

**IMPORTANT:** Some tests contain GPU-only codepaths that will be skipped if running without GPU support or, in some
cases, if there are fewer than four GPUs installed.


Continuous Integration
----------------------

Horovod uses `Buildkite <https://buildkite.com/horovod/horovod>`__ for continuous integration in AWS running on both
Intel CPU hardware and NVIDIA GPUs (with NCCL).  Tests are run once per night on master automatically, and on each
commit to a remote branch.

Buildkite test configurations are defined in
`docker-compose.test.yml <https://github.com/horovod/horovod/blob/master/docker-compose.test.yml>`__.  Each test
configuration defines a Docker image that is built from either
`Docker.test.cpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.cpu>`__ (for CPU tests) or
`Docker.test.gpu <https://github.com/horovod/horovod/blob/master/Dockerfile.test.gpu>`__ (for GPU tests).

Individual tests are run on each configuration as defined in
`gen-pipeline.sh <https://github.com/horovod/horovod/blob/master/.buildkite/gen-pipeline.sh>`__.  Every test
configuration needs to also be defined here in order to be run at test time.  Each time ``run_test`` is called
a new test artifact will be generated in Buildkite that either succeeds or fails depending on exit code.

In our AWS configuration, GPU tests are run with 4 GPUs per container. Most tests are run with 2 worker processes
each, however, model parallelism require 2 GPUs per worker, requiring 4 GPUs total.


Documentation
-------------

The Horovod documentation is published to https://horovod.readthedocs.io/.

Those HTML pages can be rendered from ``.rst`` files located in the `docs` directory.
You need to set up Sphinx before you compile the documentation the first time:

.. code-block:: bash

    $ cd docs
    $ pip install -r requirements.txt
    $ make clean

Then you can build the HTML pages and open ``docs/_build/html/index.html``:

.. code-block:: bash

    $ cd docs
    $ make html
    $ open _build/html/index.html

Sphinx can render the documentation in many other formats. Type ``make`` to get a list of available formats.


Adding Custom Operations
------------------------

Operations in Horovod are used to transform Tensors across workers.  Horovod currently supports operations that
implement Broadcast, Allreduce, and Allgather interfaces.  Gradients in Horovod are aggregated through
Allreduce operations (with the exception of sparse gradients, which use Allgather).

All data transfer operations are implemented in the
`horovod/common/ops <https://github.com/horovod/horovod/tree/master/horovod/common/ops>`__ directory.  Implementations
are organized by the collective communication library used to perform the operation (e.g.,
`mpi_operations.cc <https://github.com/horovod/horovod/blob/master/horovod/common/ops/mpi_operations.cc>`__ for MPI).

To create a new custom operation, start by defining a new class that inherits from the base operation, in the file
corresponding to the library you'll use to implement the operation:

.. code-block:: c++

    class CustomAllreduce : public AllreduceOp {
    public:
      CustomAllreduce(MPIContext* mpi_context, HorovodGlobalState* global_state);

      virtual ~CustomAllreduce() = default;

      Status Execute(std::vector<TensorTableEntry>& entries, const Response& response) override;

      bool Enabled(const ParameterManager& parameter_manager,
                   const std::vector<TensorTableEntry>& entries,
                   const Response& response) const override;

The ``Execute`` member function is responsible for performing the operation on a list of Tensors. The ``entries``
parameter provides access to all the Tensor buffers and metadata that need to be processed,
and the ``response`` parameter contains additional metadata including which devices are being used by different ranks.

``Enabled`` should return true if your operation can be performed on the given Tensor entries subject to the
current parameter settings and response metadata.

Once you've written the implementation for your operation, add it to the ``OperationManager`` in the
``CreateOperationManager`` function of
`operations.cc <https://github.com/horovod/horovod/blob/master/horovod/common/operations.cc>`__.  Because more than one
operation may be *enabled* at a time, but only one will be performed on a given vector of Tensor entries, consider the
order of your operation in the ``OperationManager`` vector before adding it in.

The first operations in the vector will be checked before those at the end, and the first operation that is *enabled*
will be performed. Broadly, the order of operations should be:

1. Custom operations that trigger based on parameters configured at runtime (e.g., ``NCCLHierarchicalAllreduce``).
2. Accelerated operations that take advantage of specialized hardware where available (e.g., ``NCCLAllreduce``).
3. Default operations that can run using standard CPUs and host memory (e.g., ``MPIAllreduce``).

Most custom operations that require preconditions such as runtime flags will fall into the first category.


Adding Compression Algorithms
-----------------------------

Gradient compression is used to reduce the amount of data sent over the network during an Allreduce operation.  Such
compression algorithms are implemented per framework (TensorFlow, PyTorch, MXNet, etc.) in
``horovod/[framework]/compression.py``
(see: `TensorFlow <https://github.com/horovod/horovod/blob/master/horovod/tensorflow/compression.py>`__,
`PyTorch <https://github.com/horovod/horovod/blob/master/horovod/torch/compression.py>`__).

To implement a new compression algorithm, first add a new class inheriting from ``Compressor``:

.. code-block:: python

    class CustomCompressor(Compressor):
        @staticmethod
        def compress(tensor):
            # do something here ...
            return tensor_compressed, ctx

        @staticmethod
        def decompress(tensor, ctx):
            # do something here ...
            return tensor_decompressed

The ``compress`` method takes a Tensor gradient and returns it in its compressed form, along with any additional context
necessary to decompress the tensor back to its original form.  Similarly, ``decompress`` takes in a compressed tensor
with its context and returns a decompressed tensor.  Compression can be done in pure Python, or in C++ using a custom
op (e.g., in `mpi_ops.cc <https://github.com/horovod/horovod/blob/master/horovod/tensorflow/mpi_ops.cc>`__ for
TensorFlow).

Once implemented, add your ``Compressor`` subclass to the ``Compressor`` class, which emulates an enumeration API:

.. code-block:: python

    class Compression(object):
        # ...

        custom = CustomCompressor

Finally, you can start using your new compressor by passing it to the ``DistributedOptimizer``:

.. code-block:: python

    opt = hvd.DistributedOptimizer(opt, compression=hvd.Compression.custom)


Horovod on Spark
----------------

The ``horovod.spark`` package makes it easy to run Horovod jobs in Spark clusters. The following section
outlines how Horovod orchestrates Spark and MPI.

Your Horovod job becomes the Spark driver and creates ``num_proc`` tasks on the Spark cluster (``horovod.spark._make_spark_thread``).
Each task runs ``horovod.spark._task_fn`` that registers with the driver, so that the driver knows when all
tasks are up and which IP and port they are running at. They also send their host hash, a string that
is treated by MPI as a hostname.

**Note:** Horovod expects all tasks to run at the same time, so your cluster has to provide at least ``num_proc`` cores to your Horovod job.
There can be multiple cores per executor, so an executor can process multiple tasks. Hosts can also have multiple executors.

The driver signals all tasks that all other tasks are up running. Each task continues initialisation
and then waits for the RPC to terminate.

After signalling all tasks are up, the driver runs ``mpi_run`` to launch the Python function in those tasks (RPC).
Usually, MPI connects to the hosts via SSH, but this would not allow to launch the Python function inside the Spark executors.
Therefore, MPI connects to each executor by invoking the ``horovod.spark.driver.mpirun_rsh`` method to "remote shell"
into the executors. This method communicates with the task that has the smallest index per host hash.
This task executes the ``orted`` command provided by MPI.
This way, a single ``orted`` process runs per executor, even if the executor has multiple cores / tasks.
MPI then uses `orted` to launch the Python function for that executor.
There will be one Python function running per core in each executor inside the first task.
All other tasks with the same host hash wait for the first task to terminate.

The following diagram illustrates this process:

.. image:: _static/spark-mpi.png


Elastic Horovod on Spark
------------------------

Elastic Horovod on Spark has a few constraints:

- each host has at most a single slot, which simplifies auto-scaling on Spark
  - for this the host hash includes the index of the task
  - this dis-allows shared memory across tasks running on the same host
  - see "Host Hash" below.


Host Hash
~~~~~~~~~

The host hash represents a single unit of processing power that shares memory. Usually, this is a regular host.
In scenarios where YARN is used to allocate cores for your Spark job, memory allocation is only shared within an executor.
There can be multiple executors running for your Horovod job on the same host, but they have each limited memory allocation.
Hence each executor gets its own host hash.

If you require each Python function to run in their own task process within a Spark executor,
then the index of the task has to become part of the host hash as well. This has only been shown useful
for Elastic Horovod on Spark, but there only for simplification.


Release Process
---------------

This section applies to contributors with permissions to release new versions of Horovod to the public.


Version Bump
~~~~~~~~~~~~

Make a PR that changes ``__version__ in horovod/__init__.py``.  Example:
`#1352 <https://github.com/horovod/horovod/pull/1352>`_.


Tag
~~~

.. code-block:: bash

    $ git tag -a v0.18.0 -m "Horovodrun config file, bugfixes"
    $ git push origin v0.18.0

Create Release
~~~~~~~~~~~~~~

Follow the GitHub instructions for `Creating a Release <https://docs.github.com/en/github/administering-a-repository/releasing-projects-on-github/managing-releases-in-a-repository#creating-a-release>`_.

Once the release has been created, this will trigger a workflow that uploads the Horovod source distribution to `PyPI <https://pypi.org>`_ automatically using `Twine <https://pypi.org/project/twine>`_.

After the workflow completes, verify that the latest version of Horovod is now available:

.. code-block:: bash

    $ pip install --upgrade horovod

.. inclusion-marker-end-do-not-remove
Horovod with TensorFlow
=======================
To use Horovod with TensorFlow, make the following modifications to your training script:

1. Run ``hvd.init()``.

.. raw:: html

    <p/>

2. Pin each GPU to a single process.

   With the typical setup of one GPU per process, set this to *local rank*. The first process on
   the server will be allocated the first GPU, the second process will be allocated the second GPU, and so forth.

   For **TensorFlow v1**:

   .. code-block:: python

       config = tf.ConfigProto()
       config.gpu_options.visible_device_list = str(hvd.local_rank())

   For **TensorFlow v2**:

   .. code-block:: python

       gpus = tf.config.experimental.list_physical_devices('GPU')
       for gpu in gpus:
           tf.config.experimental.set_memory_growth(gpu, True)
       if gpus:
           tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

.. raw:: html

    <p/>


3. Scale the learning rate by the number of workers.

   Effective batch size in synchronous distributed training is scaled by the number of workers.
   An increase in learning rate compensates for the increased batch size.

.. raw:: html

    <p/>


4. Wrap the optimizer in ``hvd.DistributedOptimizer``.

   The distributed optimizer delegates gradient computation to the original optimizer, averages gradients using **allreduce** or **allgather**, and then applies those averaged gradients.

   For **TensorFlow v2**, when using a ``tf.GradientTape``, wrap the tape in ``hvd.DistributedGradientTape`` instead of wrapping the optimizer.

.. raw:: html

    <p/>


5. Broadcast the initial variable states from rank 0 to all other processes.

   This is necessary to ensure consistent initialization of all workers when training is started with random weights or restored from a checkpoint.

   For **TensorFlow v1**, add ``hvd.BroadcastGlobalVariablesHook(0)`` when using a ``MonitoredTrainingSession``.
   When not using ``MonitoredTrainingSession``, execute the ``hvd.broadcast_global_variables`` op after global variables have been initialized.

   For **TensorFlow v2**, use ``hvd.broadcast_variables`` after models and optimizers have been initialized.

.. raw:: html

    <p/>


6. Modify your code to save checkpoints only on worker 0 to prevent other workers from corrupting them.

   For **TensorFlow v1**, accomplish this by passing ``checkpoint_dir=None`` to ``tf.train.MonitoredTrainingSession`` if ``hvd.rank() != 0``.

   For **TensorFlow v2**, construct a ``tf.train.Checkpoint`` and only call ``checkpoint.save()`` when ``hvd.rank() == 0``.

.. raw:: html

    <p/>


TensorFlow v1 Example (see the `examples <https://github.com/horovod/horovod/blob/master/examples/>`_ directory for full training examples):

.. code-block:: python

    import tensorflow as tf
    import horovod.tensorflow as hvd


    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    config = tf.ConfigProto()
    config.gpu_options.visible_device_list = str(hvd.local_rank())

    # Build model...
    loss = ...
    opt = tf.train.AdagradOptimizer(0.01 * hvd.size())

    # Add Horovod Distributed Optimizer
    opt = hvd.DistributedOptimizer(opt)

    # Add hook to broadcast variables from rank 0 to all other processes during
    # initialization.
    hooks = [hvd.BroadcastGlobalVariablesHook(0)]

    # Make training operation
    train_op = opt.minimize(loss)

    # Save checkpoints only on worker 0 to prevent other workers from corrupting them.
    checkpoint_dir = '/tmp/train_logs' if hvd.rank() == 0 else None

    # The MonitoredTrainingSession takes care of session initialization,
    # restoring from a checkpoint, saving to a checkpoint, and closing when done
    # or an error occurs.
    with tf.train.MonitoredTrainingSession(checkpoint_dir=checkpoint_dir,
                                           config=config,
                                           hooks=hooks) as mon_sess:
      while not mon_sess.should_stop():
        # Perform synchronous training.
        mon_sess.run(train_op)

TensorFlow v2 Example (from the `MNIST <https://github.com/horovod/horovod/blob/master/examples/tensorflow2/tensorflow2_mnist.py>`_ example):

.. code-block:: python

    import tensorflow as tf
    import horovod.tensorflow as hvd

    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    gpus = tf.config.experimental.list_physical_devices('GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    if gpus:
        tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

    # Build model and dataset
    dataset = ...
    model = ...
    loss = tf.losses.SparseCategoricalCrossentropy()
    opt = tf.optimizers.Adam(0.001 * hvd.size())

    checkpoint_dir = './checkpoints'
    checkpoint = tf.train.Checkpoint(model=model, optimizer=opt)

    @tf.function
    def training_step(images, labels, first_batch):
        with tf.GradientTape() as tape:
            probs = mnist_model(images, training=True)
            loss_value = loss(labels, probs)

        # Horovod: add Horovod Distributed GradientTape.
        tape = hvd.DistributedGradientTape(tape)

        grads = tape.gradient(loss_value, mnist_model.trainable_variables)
        opt.apply_gradients(zip(grads, mnist_model.trainable_variables))

        # Horovod: broadcast initial variable states from rank 0 to all other processes.
        # This is necessary to ensure consistent initialization of all workers when
        # training is started with random weights or restored from a checkpoint.
        #
        # Note: broadcast should be done after the first gradient step to ensure optimizer
        # initialization.
        if first_batch:
            hvd.broadcast_variables(mnist_model.variables, root_rank=0)
            hvd.broadcast_variables(opt.variables(), root_rank=0)

        return loss_value

    # Horovod: adjust number of steps based on number of GPUs.
    for batch, (images, labels) in enumerate(dataset.take(10000 // hvd.size())):
        loss_value = training_step(images, labels, batch == 0)

        if batch % 10 == 0 and hvd.local_rank() == 0:
            print('Step #%d\tLoss: %.6f' % (batch, loss_value))

    # Horovod: save checkpoints only on worker 0 to prevent other workers from
    # corrupting it.
    if hvd.rank() == 0:
        checkpoint.save(checkpoint_dir)
.. include:: ./ray.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./oneccl.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove.. inclusion-marker-start-do-not-remove

AdaSum with Horovod
===================

The Adaptive Summation, or AdaSum, is a novel algorithm for improving distributed data parallel training of Deep Learning models. This improvement can be seen in different ways: reducing the number steps to achieve the same accuracy in some cases and allowing you to scale to more training workers without penalizing learning rate and convergence stability.
AdaSum can be used with Horovod and PyTorch/TensorFlow. 

|

.. Contents::

|


Introduction to the AdaSum Algorithm
------------------------------------


Scaling DNN training to many GPUs always comes at a convergence degradation. This is because with larger batch sizes, gradients are averaged and the learning rate per example is smaller. To address this, learning rate is usually scaled up, but this can lead to divergence of model parameters. AdaSum addresses these two issues without introducing any hyperparameter.

Suppose there are two almost-parallel gradients from two different GPUs, g1 and g2, and they need to be reduced as shown in the figure below. The two common practices for reductions are g1+g2, the gray vector, or (g1+g2)/2, the green vector. g1+g2 may cause divergence of the model since it is effectively moving in the direction of g1 or g2 by two times the magnitude of g1 or g2. Therefore, generally (g1+g2)/2 is safer and more desired. Note that (g1+g2)/2 penalizes both the components g1 and g2 equally.

.. image:: media/abc4d31f19a315321553564e2225615b.png

Now consider the two orthogonal gradients g1 and g2 in the figure below. Since g1 and g2 are in two different dimensions and independent of each other, g1+g2 may not cause divergence.

.. image:: media/173cffdbdc89620287996ac28ca4a9ae.png

Finally, consider the third scenario where g1 and g2 are neither parallel nor orthogonal as shown in the figure below. In such a case, where taking the sum might cause a divergence, AdaSum controls the effect of the overall gradient update by subtracting half of g1’s projection on g2(pink vector) from g2, subtracting half of g2’s projection on g1 (orange vector) from g1, and summing the two components together.

.. image:: media/afa201b07a16fb29829dd8390aa0cc07.png

.. image:: media/d9b318cc2d8c16fe4ade2fa73ad83ec6.png

This formula reduces to a sum when g1 and g2 are orthogonal and an average when g1 and g2 are parallel.

This idea extends to many gradients as well. Suppose there are 2\^n gradients coming from 2\^n different GPUs. AdaSum inductively takes pairs of gradients and reduces them using the method above until all of them are reduced into one gradient. Thus, AdaSum needs the number of nodes to be a power of 2 in the current implementation.


The Distributed Optimizer for AdaSum
------------------------------------


AdaSum uses the Distributed AdaSum Optimizer to update the weights of the model after each step. In the usual data-parallel training scenario, the gradients are calculated independently by backpropagating on all the nodes, doing a reduce (averaging the gradients) so that all the nodes now have the same gradients, and then updating the weights of the model.

The distributed optimizer for AdaSum first obtains the local gradients from the backpropagation step from the current local mini batch. Instead of performing the reduce at this point, it applies the optimization function to the local gradients to perform the weight update. Then, the delta, which is the difference in the weights before and after the update is obtained, which is then reduced instead of the gradients. Once all the workers have the same delta, the weight update step is then performed as the sum of the initial weights and delta.

Since the nature of AdaSum requires it to operate on the full magnitude of the gradient, the newly added distributed optimizer uses the difference in magnitude of weights between before and after the optimizer performs a step to deliver a more accurate estimation.


Installation and Usage Instructions
-----------------------------------


AdaSum can be used and experimented with Horovod and Pytorch/TensorFlow.

In addition, there are two options of using AdaSum with Horovod: with Message Passing Interface (MPI) and with `NCCL <https://developer.nvidia.com/nccl>`_. 
Any valid implementation of MPI can be used, but AdaSum has been tested with `OpenMPI <https://www.open-mpi.org/>`_ and `IntelMPI <https://software.intel.com/en-us/mpi-library>`_.

Setting up the environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

Below are the requirements for running Horovod with AdaSum:

-   cuda >= 6.0

-   OpenMPI >= 3.0

-   NCCL >= 2.0

-   Pytorch >= 1.2.0 OR

-   Tensorflow >= 1.11.0, < 2.0

-   Horovod >= 0.18.2

*Using NCCL:*

If the **HOROVOD_GPU_OPERATIONS=NCCL** flag is used to compile Horovod, NCCL is used instead. In this case, NCCL will be used for intra-node communication, and AdaSum will be used for inter-node communication.

Modes of Operation
------------------

Adasum can be used in the following ways depending on the hardware setup available.

Pure CPU
^^^^^^^^

When dealing with a hardware setup of multiple nodes, each node having worker GPUs that are not connected by a high speed interconnect like `NVLink <https://www.nvidia.com/en-us/data-center/nvlink/>`_, where the communication happens through the CPU, AdaSum through MPI can be used for both intra-node and inter-node communication. In this case, all of the AdaSum ops are performed on the CPU.

If the hardware setup allows for a different mode like Ring or Hierarchical to be used, those must be used instead to get the highest performance benefit.

.. image:: media/7220c70747b40ab58fce2dc246958218.png

Ring
^^^^

On specifically configured machines (`DGX1 <https://www.nvidia.com/en-us/data-center/dgx-1/>`_ nodes with 8 GPUs each), the Ring mode can be used instead of the pure CPU mode. This mode is identical to the pure CPU mode for inter-node communication, but is able to do intra-node communication without going through the CPU. It does this by utilizing CUDA-aware MPI (OpenMPI built with `UCX <https://www.openucx.org/>`_ support) in order to allow direct GPU to GPU communication within nodes. This results in identical convergence benefits to pure CPU mode, but much better throughput on nodes that support it.

Ring mode is currently supported only on **DGX1** nodes having 8 GPUs each.

.. image:: media/4920a765a77fa6eeca28c5aceaa405ec.png

Hierarchical
^^^^^^^^^^^^

In cases where the hardware does not support Ring mode, but throughput higher than that of the pure CPU mode is desired, the hierarchical mode can be used instead.

The hierarchical mode functions similar to the Ring mode, except for using NCCL to do regular averaging intra-node, instead of using CUDA-aware MPI to do an AdaSum-like ring. Note that hierarchical also works on any hardware configuration, and is not limited to DGX1s.

In practice, hierarchical yields the best throughput, but lowers the convergence benefits of AdaSum due to some of the ops being regular averaging. As a rule of thumb, typically the convergence benefit degradation is insignificant on clusters with large numbers of nodes (\>=8), as in that case there are enough inter-node AdaSum ops being performed. This is the ideal Hierarchical scenario.

The other reason to use Hierarchical even on smaller clusters is when Ring mode is not supported, and CPU mode throughput is simply too low to be viable. Note that in these cases the convergence benefits compared to not using AdaSum at all might be minor.

The learning rate that should be used is equal to the best learning rate for a single worker (GPU) scaled by the number of GPUs locally on a node. On very large clusters, scaling this even more by another factor of 1.5-2.0x might give better results but is not guaranteed and should be tried only if scaling by just the local size is not sufficient for good convergence

.. image:: media/a254d38d0e56319c0507a16ea09df959.png

Modification to the Code
------------------------

A new distributed optimizer has been added to both TensorFlow and Pytorch to support the AdaSum algorithm.

An optional parameter **op** has been added to DistributedOptimizer and allreduce API for users to specify which operation to perform.
When **op=hvd.AdaSum** is specified, the new optimizer will be used.

AdaSum is highly effective in scaling to large batch sizes. The **backward_passes_per_step** parameter of the DistributedOptimizer can be used for gradient accumulation in order to scale to larger effective batch sizes without being limited by GPU memory.

TensorFlow
^^^^^^^^^^

-   DistributedOptimizer

.. code-block:: python

    opt = tf.train.AdamOptimizer(0.001)
    opt = hvd.DistributedOptimizer(opt, backward_passes_per_step=5, op=hvd.AdaSum)

-   Allreduce

.. code-block:: python
    
    hvd.allreduce(tensor, op=hvd.AdaSum)

PyTorch
^^^^^^^

-   DistributedOptimizer

.. code-block:: python

    optimizer = optim.SGD(model.parameters(), lr=args.lr, momentum=args.momentum)
    optimizer = hvd.DistributedOptimizer(optimizer, named_parameters=model.named_parameters(), compression=compression, backward_passes_per_step = 5, op=hvd.AdaSum)

-   Allreduce

.. code-block:: python

    hvd.allreduce(tensor, op=hvd.AdaSum)

Case Studies
------------


Square and Cubic optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**A simple case study to understand AdaSum’s behavior**

In order to understand the behavior and potential benefits of AdaSum as compared to Averaging, consider a simple experiment in squared optimization using AdaSum. Here, the goal is to estimate the coefficients of a polynomial of degree 2. The features are generated by randomly sampling a uniform distribution, and scaling by a factor of x_max which can be specified. This sets the complexity of the data that is used to estimate the coefficients. Additionally, the learning rate and the op to be used for Allreduce can be specified as well. The true label is calculated with the original true coefficients, without adding any noise.

In order to estimate the coefficients, Stochastic Gradient Descent is used. The training is stopped once the gradients are zero for two consecutive runs. This optimization can be run over a range of learning rates, number of workers and data range (set by x_max). This can also be modified to a cubic optimization problem.

This experiment can be run through the jupyter notebook `adasum_bench.ipynb <../examples/adasum/adasum_bench.ipynb>`_, with the models being defined in `adasum_small_model.py <../examples/adasum/adasum_small_model.py>`_.

On running experiments with a different number of workers, we can draw the following conclusions for this simple scenario with plain SGD as the optimizer:
 
-   **On the number of steps for convergence:** For the same problem, AdaSum achieves the same accuracy (100% in this case) in lower number of steps as compared to averaging. Depending on the complexity of the problem, this reduction can be anywhere up to 50% for less complex square parameter optimization.



-   **On scaling learning rate for higher number of workers**: For traditional averaging, when the number of workers is increased with local batch size the	same, this increases the global batch size, causing a higher smoothing effect on the gradients. To increase the speed of convergence, it is recommended that the learning rate be scaled up by the number of workers as	recommended in the paper `Accurate, Large Minibatch SGD: Training ImageNet	in 1 Hour <https://arxiv.org/abs/1706.02677>`_.

 **From this example, we see that with AdaSum, the LR need not be scaled linearly with the number of workers, but a better scaling factor would be 2-2.5.**


-   **On using LR decay**: With AdaSum, we see that a form of regularization effect already takes place over the gradients. As the training progresses, the magnitude of the gradients reduces, simulating the same effect as that of decaying the learning rate. Although some decay might be necessary for training more complex models, this result must be kept in mind as the same extent of decay might not be necessary.


MNIST
^^^^^

**Higher accuracy with the same number of steps**

Here, we test the applicability of the observations from the simple cubic optimization problem to training MNIST with AdaSum. By scaling the best learning rate for a single worker case by 2.5 while using AdaSum with higher number of nodes, we see that we consistently get better accuracy with the same number of steps as compared to averaging.


|

Key Takeaways
-------------

|

-   AdaSum ensures correct convergence behavior even with large effective batch sizes.

-   As the number of ranks scales up, the learning rate does not need to be scaled linearly if using CPU to do AdaSum reduction. A good scaling factor would be between 2\-2.5 over the best learning rate for a single worker.

-   If the HOROVOD_GPU_OPERATIONS=NCCL flag is used to compile Horovod, the learning rate that should be used is equal to the best learning rate for a single	worker (GPU) scaled by the number of GPUs locally on a node. On very large	clusters, scaling this even more by another factor of 1.5\-2.0x might give	better results but is not guaranteed and should be tried only if scaling by just the local size is not sufficient for good convergence.

-   Pytorch training in fp16 format is not yet supported. Integration of Apex	into the new optimizer to enabled full mixed precision training with AdaSum in Pytorch is a work in progress.

-   When HOROVOD_GPU_OPERATIONS=NCCL flag is used to compile Horovod and training	is run on a single node, only averaging through NCCL library is used to	perform reductions and no AdaSum algorithm will take place in this configuration.

.. inclusion-marker-end-do-not-remove
.. include:: ./inference.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./running.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./gpus.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./timeline.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./concepts.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove

.. inclusion-marker-start-do-not-remove


Benchmarks
==========


.. image:: https://user-images.githubusercontent.com/16640218/38965607-bf5c46ca-4332-11e8-895a-b9c137e86013.png
   :alt: 512-GPU Benchmark


The above benchmark was done on 128 servers with 4 Pascal GPUs each connected by a RoCE-capable 25 Gbit/s network. Horovod
achieves 90% scaling efficiency for both Inception V3 and ResNet-101, and 68% scaling efficiency for VGG-16.

To reproduce the benchmarks:

1. Install Horovod using the instructions provided on the `Horovod on GPU <https://github.com/horovod/horovod/blob/master/docs/gpus.rst>`__ page.

2. Clone `https://github.com/tensorflow/benchmarks <https://github.com/tensorflow/benchmarks>`__

.. code-block:: bash

    $ git clone https://github.com/tensorflow/benchmarks
    $ cd benchmarks


3. Run the benchmark. Examples below are for Open MPI.

.. code-block:: bash

    $ horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 \
        python scripts/tf_cnn_benchmarks/tf_cnn_benchmarks.py \
            --model resnet101 \
            --batch_size 64 \
            --variable_update horovod


4. At the end of the run, you will see the number of images processed per second:

.. code-block:: bash

    total images/sec: 1656.82


Real data benchmarks
~~~~~~~~~~~~~~~~~~~~
The benchmark instructions above are for the synthetic data benchmark.

To run the benchmark on a real data, you need to download the `ImageNet dataset <http://image-net.org/download-images>`__
and convert it using the TFRecord `preprocessing script <https://github.com/tensorflow/models/blob/master/research/slim/datasets/download_and_convert_imagenet.sh>`__.

Now, simply add ``--data_dir /path/to/imagenet/tfrecords --data_name imagenet --num_batches=2000`` to your training command:

.. code-block:: bash

    $ horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 \
        python scripts/tf_cnn_benchmarks/tf_cnn_benchmarks.py \
            --model resnet101 \
            --batch_size 64 \
            --variable_update horovod \
            --data_dir /path/to/imagenet/tfrecords \
            --data_name imagenet \
            --num_batches=2000


Horovod synthetic benchmarks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Horovod also comes with out-of-the-box benchmarking support for
`TensorFlow v1 <https://github.com/horovod/horovod/blob/master/examples/tensorflow/tensorflow_synthetic_benchmark.py>`__,
`TensorFlow v2 <https://github.com/horovod/horovod/blob/master/examples/tensorflow2/tensorflow2_synthetic_benchmark.py>`__, and
`PyTorch <https://github.com/horovod/horovod/blob/master/examples/pytorch/pytorch_synthetic_benchmark.py>`__.

These benchmarks allow you to measure Horovod's performance and scalability in your environment, as well as try advanced
Horovod features like gradient compression:

.. code-block:: bash

    $ horovodrun -np 4 -H server1:2,server2:2 \
        python --fp16-allreduce tensorflow2_synthetic_benchmark.py

When diagnosing performance issues, we recommend running these synthetic benchmarks first to ensure that the issues are
not originating from the training script itself.

.. inclusion-marker-end-do-not-remove
Horovod with PyTorch
====================
To use Horovod with PyTorch, make the following modifications to your training script:

1. Run ``hvd.init()``.

.. raw:: html

    <p/>

2. Pin each GPU to a single process.

   With the typical setup of one GPU per process, set this to *local rank*. The first process on
   the server will be allocated the first GPU, the second process will be allocated the second GPU, and so forth.

   .. code-block:: python

       if torch.cuda.is_available():
           torch.cuda.set_device(hvd.local_rank())

.. raw:: html

    <p/>

3. Scale the learning rate by the number of workers.

   Effective batch size in synchronous distributed training is scaled by the number of workers.
   An increase in learning rate compensates for the increased batch size.

.. raw:: html

    <p/>

4. Wrap the optimizer in ``hvd.DistributedOptimizer``.

   The distributed optimizer delegates gradient computation to the original optimizer, averages gradients using *allreduce* or *allgather*, and then applies those averaged gradients.

.. raw:: html

    <p/>

5. Broadcast the initial variable states from rank 0 to all other processes:

   .. code-block:: python

       hvd.broadcast_parameters(model.state_dict(), root_rank=0)
       hvd.broadcast_optimizer_state(optimizer, root_rank=0)

   This is necessary to ensure consistent initialization of all workers when training is started with random weights or restored from a checkpoint.

.. raw:: html

    <p/>

6. Modify your code to save checkpoints only on worker 0 to prevent other workers from corrupting them.

   Accomplish this by guarding model checkpointing code with ``hvd.rank() != 0``.

.. raw:: html

    <p/>

Example (also see a full training `example <https://github.com/horovod/horovod/blob/master/examples/pytorch/pytorch_mnist.py>`__):

.. code-block:: python

    import torch
    import horovod.torch as hvd

    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    torch.cuda.set_device(hvd.local_rank())

    # Define dataset...
    train_dataset = ...

    # Partition dataset among workers using DistributedSampler
    train_sampler = torch.utils.data.distributed.DistributedSampler(
        train_dataset, num_replicas=hvd.size(), rank=hvd.rank())

    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=..., sampler=train_sampler)

    # Build model...
    model = ...
    model.cuda()

    optimizer = optim.SGD(model.parameters())

    # Add Horovod Distributed Optimizer
    optimizer = hvd.DistributedOptimizer(optimizer, named_parameters=model.named_parameters())

    # Broadcast parameters from rank 0 to all other processes.
    hvd.broadcast_parameters(model.state_dict(), root_rank=0)

    for epoch in range(100):
       for batch_idx, (data, target) in enumerate(train_loader):
           optimizer.zero_grad()
           output = model(data)
           loss = F.nll_loss(output, target)
           loss.backward()
           optimizer.step()
           if batch_idx % args.log_interval == 0:
               print('Train Epoch: {} [{}/{}]\tLoss: {}'.format(
                   epoch, batch_idx * len(data), len(train_sampler), loss.item()))


.. NOTE:: PyTorch GPU support requires NCCL 2.2 or later. It also works with NCCL 2.1.15 if you are not using RoCE or InfiniBand.


PyTorch Lightning
-----------------

Horovod is supported as a distributed backend in `PyTorch Lightning <https://github.com/PyTorchLightning/pytorch-lightning>`_ from v0.7.4 and above.

With PyTorch Lightning, distributed training using Horovod requires only a single line code change to your existing training script:

.. code-block:: python

    # train Horovod on GPU (number of GPUs / machines provided on command-line)
    trainer = pl.Trainer(accelerator='horovod', gpus=1)

    # train Horovod on CPU (number of processes / machines provided on command-line)
    trainer = pl.Trainer(accelerator='horovod')

May need to change parameter "accelerator" name to "distributed_backend" in some older version of pytorch_lightning.

Start the training job and specify the number of workers on the command line as you normally would when using Horovod:

.. code-block:: bash

    # run training with 4 GPUs on a single machine
    $ horovodrun -np 4 python train.py

    # run training with 8 GPUs on two machines (4 GPUs each)
    $ horovodrun -np 8 -H hostname1:4,hostname2:4 python train.py

You can find an example of using pytorch lightning trainer with horovod backend in `pytorch_lightning_mnist.py 
<../examples/pytorch/pytorch_lightning_mnist.py>`__

See the PyTorch Lightning `docs <https://pytorch-lightning.readthedocs.io/en/stable/multi_gpu.html#horovod>`_ for more details.

A Pytorch-Lightning based spark estimator is also added, example is in `pytorch_lightning_spark_mnist.py <../examples/spark/pytorch/pytorch_lightning_spark_mnist.py>`__
.. include:: ./process_set.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove.. inclusion-marker-start-do-not-remove

Inference
=========

What about inference? Inference may be done outside of the Python script that was used to train the model. If you do this, it
will not have references to the Horovod library.

To run inference on a checkpoint generated by the Horovod-enabled training script you should optimize the graph and only
keep operations necessary for a forward pass through model. The `Optimize for Inference <https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/tools/optimize_for_inference.py>`__
script from the TensorFlow repository will do that for you.

If you want to convert your checkpoint to `Frozen Graph <https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/tools/freeze_graph.py>`__,
you should do so after doing the optimization described above, otherwise the `Freeze Graph <https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/tools/freeze_graph.py>`__
script will fail to load Horovod op:

.. code-block:: bash

    ValueError: No op named HorovodAllreduce in defined operations.


.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Horovod with Intel(R) oneCCL
============================
To use Horovod with the Intel(R) oneAPI Collective Communications Library (oneCCL), follow the steps below.

1. Install `oneCCL <https://github.com/intel/oneccl>`_.

To install oneCCL, follow `these steps <https://github.com/intel/oneccl/blob/master/README.md>`_.

Source ``setvars.sh`` to start using oneCCL.

.. code-block:: bash

    source <install_dir>/env/setvars.sh

2. Set ``HOROVOD_CPU_OPERATIONS`` variable
    
.. code-block:: bash

    export HOROVOD_CPU_OPERATIONS=CCL

3. Install Horovod from source code

.. code-block:: bash

    python setup.py build
    python setup.py install

or via pip 

.. code-block:: bash
    
    pip install horovod

Advanced settings
*****************

Affinity
--------

You can specify the affinity for Horovod background thread with the ``HOROVOD_THREAD_AFFINITY`` environment variable.
See the instructions below.

Set Horovod background thread affinity according to the rule - if there is N Horovod processes per node, this variable should contain all the values for every local process using comma as a separator:

.. code-block:: bash
    
    export HOROVOD_THREAD_AFFINITY=c0,c1,...,c(N-1)

where c0,...,c(N-1) are core IDs to pin background threads from local processes.


Set the number of oneCCL workers:

.. code-block:: bash
    
    export CCL_WORKER_COUNT=X

where X is the number of oneCCL worker threads (workers) per process you'd like to dedicate to drive communication.


Set oneCCL workers affinity automatically:

.. code-block:: bash

    export CCL_WORKER_AFFINITY=auto

This is default mode. The exact core IDs will depend from process launcher used.

Set oneCCL workers affinity explicitly:

.. code-block:: bash

    export CCL_WORKER_AFFINITY=c0,c1,..,c(X-1)

where c0,c1,..,c(X-1) are core IDs dedicated to local oneCCL workers, i.e. X = ``CCL_WORKER_COUNT`` * Number of processes per node.

Please refer to `Execution of Communication Operations <https://oneapi-src.github.io/oneCCL/operation_execution.html>`_ for more information.


For example, we have 2 nodes and each node has 2 sockets: socket0 CPUs: 0-17,36-53 and socket1 CPUs: 18-35,54-71. We dedicate the last two cores of each socket for 2 oneCCL workers and pin Horovod background thread to one of the hyper-thread cores of oneCCL workers's cores. All these cores are excluded from Intel MPI pinning using ``I_MPI_PIN_PROCESSOR_EXCLUDE_LIST`` to dedicate them to oneCCL and Horovod tasks only, thus avoiding the conflict with framework's computational threads.

.. code-block:: bash
    
    export CCL_WORKER_COUNT=2
    export CCL_WORKER_AFFINITY="16,17,34,35"
    export HOROVOD_THREAD_AFFINITY="53,71"
    export I_MPI_PIN_DOMAIN=socket
    export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST="16,17,34,35,52,53,70,71"

    mpirun -n 4 -ppn 2 -hostfile hosts python ./run_example.py


Caching
-------

Set cache hint for oneCCL operations:

.. code-block:: bash
    
    export HOROVOD_CCL_CACHE=0|1

Available for ``allreduce`` only yet. Disabled by default.

Please refer to `Caching of Communication Operations <https://oneapi-src.github.io/oneCCL/operation_caching.html>`_ for more information.

.. inclusion-marker-end-do-not-remove
Horovod with XLA in Tensorflow
===============================

Basic usage
-----------

XLA Horovod ops can be enabled by setting ``HOROVOD_ENABLE_XLA_OPS = 1`` by controlling the registration of the ops to Tensorflow/XLA.

There are two main ways to enable XLA and they could work with Horovod in different ways:

For **Explicit compilation with tf.function(jit_compile=True)**:

.. code-block:: python

    os.environ["HOROVOD_ENABLE_XLA_OPS"] = "1"

     @tf.function(jit_compile=True)
     def compiled_hvd_allreduce(self, dtype, dim):
         tensor = self.random_uniform(
             [17] * dim, -100, 100, dtype=dtype)
         summed = hvd.allreduce(tensor, average=False)
         return summed

In this way, all the ops in the ``compiled_hvd_allreduce`` function are lowered into XLA per the compilation requirement. If the XLA Horovod ops are not enabled, XLA will report compilation errors.


For **Auto-clustering**:

Auto-clustering is a convenient way to use XLA by simply setting ``TF_XLA_FLAGS=--tf_xla_auto_jit=2`` and the XLA JIT automatically selects ops in the Tensorflow graph to be lowered into XLA. In this mode, enabling XLA Horovod ops is optional, because the auto-clustering can work even if the Horovod ops are left to be run by Tensorflow (devices) while only parts of the graphs are lowered onto XLA (devices).

List of supported XLA Horovod ops
---------------------------------

The supported op list is:

``HorovodAllreduce``

.. inclusion-marker-start-do-not-remove

Horovod on GPU
==============


To use Horovod on GPU, read the options below and see which one applies to you best.

Have GPUs?
~~~~~~~~~~
In most situations, using NCCL 2 will significantly improve performance over the CPU version.  NCCL 2 provides the **allreduce**
operation optimized for NVIDIA GPUs and a variety of networking devices, such as RoCE or InfiniBand.

1. Install `NCCL 2 <https://developer.nvidia.com/nccl>`__ following `these steps <http://docs.nvidia.com/deeplearning/sdk/nccl-install-guide/index.html>`__.

   If you have installed NCCL 2 using the ``nccl-<version>.txz`` package, you should add the library path to ``LD_LIBRARY_PATH``
   environment variable or register it in ``/etc/ld.so.conf``.

   .. code-block:: bash

       $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/nccl-<version>/lib


2. (Optional) If you're using an NVIDIA Tesla GPU and NIC with GPUDirect RDMA support, you can further speed up NCCL 2
by installing an `nv_peer_memory <https://github.com/Mellanox/nv_peer_memory>`__ driver.

   `GPUDirect <https://developer.nvidia.com/gpudirect>`__ allows GPUs to transfer memory among each other without CPU
   involvement, which significantly reduces latency and load on CPU.  NCCL 2 is able to use GPUDirect automatically for
   **allreduce** operation if it detects it.

3. Install `Open MPI <https://www.open-mpi.org/>`__ or another MPI implementation following `these steps <https://www.open-mpi.org/faq/?category=building#easy-build>`__.

   **Note**: Open MPI 3.1.3 has an issue that may cause hangs.  The recommended fix is to downgrade to Open MPI 3.1.2 or upgrade to Open MPI 4.0.0.

4. If you installed TensorFlow from `PyPI <https://pypi.org/project/tensorflow>`__, make sure that ``g++-5`` or above is installed.

   If you installed PyTorch from `PyPI <https://pypi.org/project/torch>`__, make sure that ``g++-5`` or above is installed.

   If you installed either package from `Conda <https://conda.io>`_, make sure that the ``gxx_linux-64`` Conda package is installed.

5. Install the ``horovod`` pip package.

   If you installed NCCL 2 using the ``nccl-<version>.txz`` package, you should specify the path to NCCL 2 using the ``HOROVOD_NCCL_HOME``
   environment variable.

   .. code-block:: bash

       $ HOROVOD_NCCL_HOME=/usr/local/nccl-<version> HOROVOD_GPU_OPERATIONS=NCCL pip install --no-cache-dir horovod


   If you installed NCCL 2 using the Ubuntu package, you can run:

   .. code-block:: bash

       $ HOROVOD_GPU_OPERATIONS=NCCL pip install --no-cache-dir horovod
   
   If you installed NCCL 2 using the `CentOS / RHEL package <https://docs.nvidia.com/deeplearning/sdk/nccl-install-guide/index.html#rhel_centos>`__, you can run:

   .. code-block:: bash

       $ HOROVOD_NCCL_INCLUDE=/usr/include HOROVOD_NCCL_LIB=/usr/lib64 HOROVOD_GPU_OPERATIONS=NCCL pip install --no-cache-dir horovod


**Note**: Some models with a high computation to communication ratio benefit from doing allreduce on CPU, even if a
GPU version is available. To force allreduce to happen on CPU, pass ``device_dense='/cpu:0'`` to ``hvd.DistributedOptimizer``:

.. code-block:: python

    opt = hvd.DistributedOptimizer(opt, device_dense='/cpu:0')


Advanced: Have a proprietary MPI implementation with GPU support optimized for your network?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This section is only relevant if you have a proprietary MPI implementation with GPU support, i.e. not Open MPI or MPICH.
Most users should follow one of the sections above.

If your MPI vendor's implementation of *allreduce* operation on GPU is faster than NCCL 2, you can configure Horovod to
use it instead:

.. code-block:: bash

    $ HOROVOD_GPU_ALLREDUCE=MPI pip install --no-cache-dir horovod


Additionally, if your MPI vendor's implementation supports *allgather* and *broadcast* operations on GPU, you can
configure Horovod to use them as well:

.. code-block:: bash

    $ HOROVOD_GPU_OPERATIONS=MPI pip install --no-cache-dir horovod


**Note**: Allgather allocates an output tensor which is proportionate to the number of processes participating in the
training.  If you find yourself running out of GPU memory, you can force allgather to happen on CPU by passing
``device_sparse='/cpu:0'`` to ``hvd.DistributedOptimizer``:

.. code-block:: python

    opt = hvd.DistributedOptimizer(opt, device_sparse='/cpu:0')


.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Troubleshooting
===============


Import TensorFlow failed during installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Is TensorFlow installed?

If you see the error message below, it means that TensorFlow is not installed.  Please install TensorFlow before installing Horovod.

.. code-block:: bash

    error: import tensorflow failed, is it installed?

    Traceback (most recent call last):
      File "/tmp/pip-OfE_YX-build/setup.py", line 29, in fully_define_extension
        import tensorflow as tf
    ImportError: No module named tensorflow


2. Are the CUDA libraries available?

If you see the error message below, it means that TensorFlow cannot be loaded.
If you're installing Horovod into a container on a machine without GPUs, you may use CUDA stub drivers to work around the issue.

.. code-block:: bash

    error: import tensorflow failed, is it installed?

    Traceback (most recent call last):
      File "/tmp/pip-41aCq9-build/setup.py", line 29, in fully_define_extension
        import tensorflow as tf
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/__init__.py", line 24, in <module>
        from tensorflow.python import *
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/python/__init__.py", line 49, in <module>
        from tensorflow.python import pywrap_tensorflow
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/python/pywrap_tensorflow.py", line 52, in <module>
        raise ImportError(msg)
    ImportError: Traceback (most recent call last):
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/python/pywrap_tensorflow.py", line 41, in <module>
        from tensorflow.python.pywrap_tensorflow_internal import *
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/python/pywrap_tensorflow_internal.py", line 28, in <module>
        _pywrap_tensorflow_internal = swig_import_helper()
      File "/usr/local/lib/python2.7/dist-packages/tensorflow/python/pywrap_tensorflow_internal.py", line 24, in swig_import_helper
        _mod = imp.load_module('_pywrap_tensorflow_internal', fp, pathname, description)
    ImportError: libcuda.so.1: cannot open shared object file: No such file or directory


To use CUDA stub drivers:

.. code-block:: bash

    # temporary add stub drivers to ld.so.cache
    $ ldconfig /usr/local/cuda/lib64/stubs

    # install Horovod, add other HOROVOD_* environment variables as necessary
    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install --no-cache-dir horovod

    # revert to standard libraries
    $ ldconfig


MPI is not found during installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Is MPI in PATH?

If you see the error message below, it means ``mpicxx`` was not found in PATH. Typically ``mpicxx`` is located in the same directory as ``mpirun``.
Add a directory containing ``mpicxx`` to PATH before installing Horovod.

.. code-block:: bash

    error: mpicxx -show failed, is mpicxx in $PATH?

    Traceback (most recent call last):
      File "/tmp/pip-dQ6A7a-build/setup.py", line 70, in get_mpi_flags
        ['mpicxx', '-show'], universal_newlines=True).strip()
      File "/usr/lib/python2.7/subprocess.py", line 566, in check_output
        process = Popen(stdout=PIPE, *popenargs, **kwargs)
      File "/usr/lib/python2.7/subprocess.py", line 710, in __init__
        errread, errwrite)
      File "/usr/lib/python2.7/subprocess.py", line 1335, in _execute_child
        raise child_exception
    OSError: [Errno 2] No such file or directory


To use custom MPI directory:

.. code-block:: bash

    $ export PATH=$PATH:/path/to/mpi/bin
    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install --no-cache-dir horovod


2. Are MPI libraries added to ``$LD_LIBRARY_PATH`` or ``ld.so.conf``?

If you see the error message below, it means ``mpicxx`` was not able to load some of the MPI libraries. If you recently
installed MPI, make sure that the path to MPI libraries is present the ``$LD_LIBRARY_PATH`` environment variable, or in the
``/etc/ld.so.conf`` file.

.. code-block:: bash

    mpicxx: error while loading shared libraries: libopen-pal.so.40: cannot open shared object file: No such file or directory
    error: mpicxx -show failed (see error below), is MPI in $PATH?

    Traceback (most recent call last):
    File "/tmp/pip-build-wrtVwH/horovod/setup.py", line 107, in get_mpi_flags
    shlex.split(show_command), universal_newlines=True).strip()
    File "/usr/lib/python2.7/subprocess.py", line 574, in check_output
    raise CalledProcessError(retcode, cmd, output=output)
    CalledProcessError: Command '['mpicxx', '-show']' returned non-zero exit status 127


If you have installed MPI in a user directory, you can add the MPI library directory to ``$LD_LIBRARY_PATH``:

.. code-block:: bash

    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/mpi/lib


If you have installed MPI in a non-standard system location (i.e. not ``/usr`` or ``/usr/local``), you should add it to the
``/etc/ld.so.conf`` file:

.. code-block:: bash

    $ echo /path/to/mpi/lib | sudo tee -a /etc/ld.so.conf


Additionally, if you have installed MPI in a system location, you should run ``sudo ldconfig`` after installation to
register libraries in the cache:

.. code-block:: bash

    $ sudo ldconfig


Error during installation: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below, it means that your MPI is likely outdated. We recommend installing
`Open MPI >=4.0.0 <https://www.open-mpi.org/faq/?category=building#easy-build>`__.

**Note**: Prior to installing a new version of Open MPI, don't forget to remove your existing MPI installation.

.. code-block:: bash

    horovod/tensorflow/mpi_ops.cc: In function ‘void horovod::tensorflow::{anonymous}::PerformOperation(horovod::tensorflow::{anonymous}::TensorTable&, horovod::tensorflow::MPIResponse)’:
    horovod/tensorflow/mpi_ops.cc:802:79: # error: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]
                                      recvcounts, displcmnts, dtype, MPI_COMM_WORLD);
                                                                                   ^
    In file included from horovod/tensorflow/mpi_ops.cc:38:0:
    /usr/anaconda2/include/mpi.h:633:5: error:   initializing argument 1 of ‘int MPI_Allgatherv(void*, int, MPI_Datatype, void*, int*, int*, MPI_Datatype, MPI_Comm)’ [-fpermissive]
     int MPI_Allgatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
         ^
    horovod/tensorflow/mpi_ops.cc:1102:45: error: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]
                                   MPI_COMM_WORLD))
                                                 ^


Error during installation: fatal error: pyconfig.h: No such file or directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below, it means that you need to install Python headers.

.. code-block:: bash

    build/horovod/torch/mpi_lib/_mpi_lib.c:22:24: fatal error: pyconfig.h: No such file or directory
     #  include <pyconfig.h>
                            ^
    compilation terminated.


You can do this by installing a ``python-dev`` or ``python3-dev`` package.  For example, on a Debian or Ubuntu system:

.. code-block:: bash

    $ sudo apt-get install python-dev


NCCL 2 is not found during installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below, it means NCCL 2 was not found in the standard libraries location. If you have a directory
where you installed NCCL 2 which has both ``include`` and ``lib`` directories containing ``nccl.h`` and ``libnccl.so``
respectively, you can pass it via ``HOROVOD_NCCL_HOME`` environment variable. Otherwise you can specify them separately
via ``HOROVOD_NCCL_INCLUDE`` and ``HOROVOD_NCCL_LIB`` environment variables.

.. code-block:: bash

    build/temp.linux-x86_64-2.7/test_compile/test_nccl.cc:1:18: fatal error: nccl.h: No such file or directory
     #include <nccl.h>
                      ^
    compilation terminated.
    error: NCCL 2.0 library or its later version was not found (see error above).
    Please specify correct NCCL location via HOROVOD_NCCL_HOME environment variable or combination of HOROVOD_NCCL_INCLUDE and HOROVOD_NCCL_LIB environment variables.

    HOROVOD_NCCL_HOME - path where NCCL include and lib directories can be found
    HOROVOD_NCCL_INCLUDE - path to NCCL include directory
    HOROVOD_NCCL_LIB - path to NCCL lib directory


For example:

.. code-block:: bash

    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install --no-cache-dir horovod


Or:

.. code-block:: bash

    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_INCLUDE=/path/to/nccl/include HOROVOD_NCCL_LIB=/path/to/nccl/lib pip install --no-cache-dir horovod


Pip install: no such option: --no-cache-dir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below, it means that your version of pip is out of date. You can remove the ``--no-cache-dir`` flag
since your version of pip does not do caching. The ``--no-cache-dir`` flag is added to all examples to ensure that when you
change Horovod compilation flags, it will be rebuilt from source and not just reinstalled from the pip cache, which is
modern pip's `default behavior <https://pip.pypa.io/en/stable/reference/pip_install/#caching>`__.

.. code-block:: bash

    $ pip install --no-cache-dir horovod

    Usage:
      pip install [options] <requirement specifier> ...
      pip install [options] -r <requirements file> ...
      pip install [options] [-e] <vcs project url> ...
      pip install [options] [-e] <local project path> ...
      pip install [options] <archive url/path> ...

    no such option: --no-cache-dir


For example:

.. code-block:: bash

    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install --no-cache-dir horovod


ncclAllReduce failed: invalid data type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below during the training, it means that Horovod was linked to the wrong version of NCCL
library.

.. code-block:: bash

    UnknownError (see above for traceback): ncclAllReduce failed: invalid data type
             [[Node: DistributedMomentumOptimizer_Allreduce/HorovodAllreduce_gradients_AddN_2_0 = HorovodAllreduce[T=DT_FLOAT, _device="/job:localhost/replica:0/task:0/device:GPU:0"](gradients/AddN_2)]]
             [[Node: train_op/_653 = _Recv[client_terminated=false, recv_device="/job:localhost/replica:0/task:0/device:CPU:0", send_device="/job:localhost/replica:0/task:0/device:GPU:0", send_device_incarnation=1, tensor_name="edge_1601_train_op", tensor_type=DT_FLOAT, _device="/job:localhost/replica:0/task:0/device:CPU:
    0"]()]]


If you're using Anaconda or Miniconda, you most likely have the ``nccl`` package installed. The solution is to remove
the package and reinstall Horovod:

.. code-block:: bash

    $ conda remove nccl
    $ pip uninstall -y horovod
    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install --no-cache-dir horovod


transport/p2p.cu:431 WARN failed to open CUDA IPC handle : 30 unknown error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below during the training with ``-x NCCL_DEBUG=INFO``, it likely means that multiple servers
share the same ``hostname``.

.. code-block:: bash

    node1:22671:22795 [1] transport/p2p.cu:431 WARN failed to open CUDA IPC handle : 30 unknown error


MPI and NCCL rely on hostnames to distinguish between servers, so you should make sure that every server has a unique
hostname.

Running out of memory
~~~~~~~~~~~~~~~~~~~~~
If you notice that your program is running out of GPU memory and multiple processes
are being placed on the same GPU, it's likely that your program (or its dependencies)
create a ``tf.Session`` that does not use the ``config`` that pins specific GPU.

If possible, track down the part of program that uses these additional tf.Sessions and pass the same configuration.

Alternatively, you can place following snippet in the beginning of your program to ask TensorFlow
to minimize the amount of memory it will pre-allocate on each GPU:

.. code-block:: python

    small_cfg = tf.ConfigProto()
    small_cfg.gpu_options.allow_growth = True
    with tf.Session(config=small_cfg):
        pass


As a last resort, you can **replace** setting ``config.gpu_options.visible_device_list``
with different code:

.. code-block:: python

    # Pin GPU to be used
    import os
    os.environ['CUDA_VISIBLE_DEVICES'] = str(hvd.local_rank())


**Note**: Setting ``CUDA_VISIBLE_DEVICES`` is incompatible with ``config.gpu_options.visible_device_list``.

Setting ``CUDA_VISIBLE_DEVICES`` has additional disadvantage for GPU version - CUDA will not be able to use IPC, which
will likely cause NCCL and MPI to fail.  In order to disable IPC in NCCL and MPI and allow it to fallback to shared
memory, use:
* ``export NCCL_P2P_DISABLE=1`` for NCCL.
* ``--mca btl_smcuda_use_cuda_ipc 0`` flag for OpenMPI and similar flags for other vendors.

libcudart.so.X.Y: cannot open shared object file: No such file or directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you notice that your program crashes with a ``libcudart.so.X.Y: cannot open shared object file: No such file or directory`` error, it's likely that your framework and Horovod were build with different versions of CUDA.

To build Horovod with a specific CUDA version, use the ``HOROVOD_CUDA_HOME`` environment variable during installation:

.. code-block:: bash

    $ pip uninstall -y horovod
    $ HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl HOROVOD_CUDA_HOME=/path/to/cuda pip install --no-cache-dir horovod


FORCE-TERMINATE AT Data unpack would read past end of buffer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below during the training, it's likely that you have a wrong version of ``hwloc`` installed in your system.

.. code-block:: bash

    --------------------------------------------------------------------------
    An internal error has occurred in ORTE:

    [[25215,0],1] FORCE-TERMINATE AT Data unpack would read past end of buffer:-26 - error grpcomm_direct.c(359)

    This is something that should be reported to the developers.
    --------------------------------------------------------------------------
    [future5.stanford.edu:12508] [[25215,0],1] ORTE_ERROR_LOG: Data unpack would read past end of buffer in file grpcomm_direct.c at line 355


Purge ``hwloc`` from your system:

.. code-block:: bash

    $ apt purge hwloc-nox libhwloc-dev libhwloc-plugins libhwloc5


After ``hwloc`` is purged, `re-install Open MPI <https://www.open-mpi.org/faq/?category=building#easy-build>`__.

See `this issue <https://github.com/open-mpi/ompi/issues/4437>`__ for more details.

segmentation fault with tensorflow 1.14 or higher mentioning `hwloc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are using TensorFlow 1.14 or 1.15 and are getting a segmentation fault, check whether it mentions `hwloc`:

    ...
    Signal: Segmentation fault (11)
    Signal code: Address not mapped (1)
    Failing at address: 0x99
    [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x3ef20)[0x7f309d34ff20]
    [ 1] /usr/lib/x86_64-linux-gnu/libopen-pal.so.20(opal_hwloc_base_free_topology+0x76)[0x7f3042871ca6]
    ...

If it does, this could be a conflict with the `hwloc` symbols explorted from TensorFlow.

To fix this, locate your hwloc library with `ldconfig -p | grep libhwloc.so`, and then set `LD_PRELOAD`. For example:

    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libhwloc.so python -c 'import horovod.tensorflow as hvd; hvd.init()'

See [this issue](https://github.com/horovod/horovod/issues/1123) for more information.

bash: orted: command not found
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you see the error message below during the training, it's likely that Open MPI cannot find one of its components in PATH.

.. code-block:: bash

    bash: orted: command not found
    --------------------------------------------------------------------------
    ORTE was unable to reliably start one or more daemons.
    This usually is caused by:

    * not finding the required libraries and/or binaries on
      one or more nodes. Please check your PATH and LD_LIBRARY_PATH
      settings, or configure OMPI with --enable-orterun-prefix-by-default

    * lack of authority to execute on one or more specified nodes.
      Please verify your allocation and authorities.

    * the inability to write startup files into /tmp (--tmpdir/orte_tmpdir_base).
      Please check with your sys admin to determine the correct location to use.

    *  compilation of the orted with dynamic libraries when static are required
      (e.g., on Cray). Please check your configure cmd line and consider using
      one of the contrib/platform definitions for your system type.

    * an inability to create a connection back to mpirun due to a
      lack of common network interfaces and/or no route found between
      them. Please check network connectivity (including firewalls
      and network routing requirements).
    --------------------------------------------------------------------------


We recommended reinstalling Open MPI with the ``--enable-orterun-prefix-by-default`` flag, like so:

.. code-block:: bash

    $ wget https://www.open-mpi.org/software/ompi/v4.0/downloads/openmpi-4.0.0.tar.gz
    $ tar zxf openmpi-4.0.0.tar.gz
    $ cd openmpi-4.0.0
    $ ./configure --enable-orterun-prefix-by-default
    $ make -j $(nproc) all
    $ make install
    $ ldconfig


.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Elastic Horovod
===============


Elastic training enables Horovod to scale up and down the number of workers dynamically at runtime, without
requiring a restart or resuming from checkpoints saved to durable storage. With elastic training, workers can come
and go from the Horovod job without interrupting the training process.


When to use elastic training
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- You are running an `autoscaling <https://en.wikipedia.org/wiki/Autoscaling>`__ job that may acquire more resources for training over time.
- Your job is running on preemptable or spot instances that may come and go with little warning.
- Your nodes are unreliable and you want your job to continue training if some of the hosts fail.


Requirements
~~~~~~~~~~~~

- TensorFlow >= 1.15 or PyTorch >= 1.0
- Horovod >= 0.20.0 with Gloo support (install Horovod using ``HOROVOD_WITH_GLOO=1`` to ensure it is installed)
- A way to discover available hosts at runtime


Modifying the training script with State Synchronization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The biggest difference when moving from normal distributed training to elastic training is the need to track and synchronize
state among the workers as workers are added or removed from the job.

To enable elastic training, make the following changes to your training script:

1. Wrap your main training process (everything following initialization) in a function decorated with ``hvd.elastic.run``.

   The first argument to this decorated function should be an instance of ``hvd.elastic.State``.  Before executing the
   decorated function, this state object will be synchronized across workers.  This ensures that workers that were
   newly added, as well as workers that might have inconsistent state, all share the same state before training begins.

   Because the sync function uses collective ops, and upon worker add the active workers will not reset from before this
   function, *no Horovod collective ops (broadcast, allreduce, allgather, etc.) can be called before this function*.

2. Place all variables that need to be kept in sync between worker replicas (model parameters, optimizer state, epoch and batch numbers, etc.) into a ``hvd.elastic.State`` object.

   Standard state implementations are provided for TensorFlow, Keras, and PyTorch.  However, it may be necessary in some cases to override
   the base ``hvd.elastic.State`` object to handle broadcasting custom types.

3. Periodically call ``state.commit()`` to backup a copy of your state in memory.

   This is useful to prevent corrupted state in the event that a worker fails unexpectedly. For example, if training fails
   in the middle of a parameter update, some gradient updates may have applied while others were still being allreduced.  When this
   happens, a ``HorovodInternalError`` will be raised, and all parameters will be restored to the values at the time of the last commit.

   Because commits can be expensive (as the model size increases), there is a tradeoff between the per-batch processing time
   and how far the training process needs to rollback in the event of a failure.  For example, if you commit once every 10
   batches, you reduce the amount of copying by a factor of 10. But if a failure occurs, you may need to redo up to 10
   previously processed batches.

   Elastic Horovod can avoid these rollbacks by performing what we call a *graceful removal* of a worker. If the driver
   process discovers that a host has been made available or flagged for removal, it will push a notification to the workers.
   The next time ``state.commit()`` or the more lightweight ``state.check_host_updates()`` is called, a ``HostsUpdatedInterrupt``
   will be raised.  This event is handled similar to the ``HorovodInternalError``, except that parameter state will not be
   restored to the last commit.

   In general, if your hardware is generally reliable, and your orchestration system gives the driver ample warning
   when a host is scheduled to be removed from the job, then you can safely call ``state.commit()`` on a reduced frequency,
   and call ``state.check_host_updates()`` at the end of each batch instead.

4. Register callbacks with the ``hvd.elastic.State`` object to respond to changes in the worker membership in the job.

   For example, rescaling the learning rate with the new world size or repartitioning the dataset would commonly be done
   through these callbacks.

   Callbacks are called after Horovod has reinitialized, but before state is synchronized across the workers.

The reset process following a ``HorovodInternalError`` (failure) or ``HostsUpdatedInterrupt`` (add/remove request) is as follows:

1. Catch exception within the ``hvd.elastic.run`` decorator.
2. Restore last committed state if ``HorovodInternalError`` was raised.
3. Reinitialize Horovod context performing a new round of rendezvous.
4. Synchronize state among the workers by broadcasting from the new worker-0.
5. Resume training by executing the underlying training function.

During rendezvous, older workers will take priority in being assigned worker-0 status to ensure that the state that
is broadcast is up to date.


Elastic TensorFlow
~~~~~~~~~~~~~~~~~~

TensorFlow v1 Example:

.. code-block:: python
    :emphasize-lines: 17,18,23,29,32,33

    import tensorflow as tf
    import horovod.tensorflow as hvd

    hvd.init()

    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.gpu_options.visible_device_list = str(hvd.local_rank())

    dataset = ...
    model = ...

    lr = tf.Variable(base_lr * hvd.size())
    optimizer = tf.train.GradientDescentOptimizer(lr)
    optimizer = hvd.DistributedOptimizer(optimizer)

    @hvd.elastic.run
    def train(state, train_one_batch):
        for state.epoch in range(state.epoch, epochs):
            for state.batch in range(state.batch, batches_per_epoch):
                train_one_batch()
                if state.batch % batches_per_commit == 0:
                    state.commit()
            state.batch = 0

    with tf.Session(config=config) as session:
        session.run(tf.global_variables_initializer())

        def on_state_reset():
            lr.load(base_lr * hvd.size(), session)

        state = hvd.elastic.TensorFlowState(session=session, batch=0, epoch=0)
        state.register_reset_callbacks([on_state_reset])

        train_opt = optimizer.minimize(loss)
        train(state, lambda: session.run(train_opt))

TensorFlow v2 Example:

.. code-block:: python
    :emphasize-lines: 33,34,40,43,46,47

    import tensorflow as tf
    import horovod.tensorflow as hvd

    hvd.init()

    gpus = tf.config.experimental.list_physical_devices('GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    if gpus:
        tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

    dataset = ...
    model = ...

    optimizer = tf.optimizers.Adam(lr * hvd.size())

    @tf.function
    def train_one_batch(data, target, allreduce=True):
        with tf.GradientTape() as tape:
            probs = model(data, training=True)
            loss = tf.losses.categorical_crossentropy(target, probs)

        if allreduce:
            tape = hvd.DistributedGradientTape(tape)

        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))

    # Initialize model and optimizer state so we can synchronize across workers
    data, target = get_random_batch()
    train_one_batch(data, target, allreduce=False)

    @hvd.elastic.run
    def train(state):
        for state.epoch in range(state.epoch, epochs):
            for state.batch in range(state.batch, batches_per_epoch):
                data, target = get_random_batch()
                train_one_batch(data, target)
                if state.batch % batches_per_commit == 0:
                    state.commit()
            state.batch = 0

    def on_state_reset():
        optimizer.lr.assign(lr * hvd.size())

    state = hvd.elastic.TensorFlowKerasState(model, optimizer, batch=0, epoch=0)
    state.register_reset_callbacks([on_state_reset])
    train(state)


Elastic Keras
~~~~~~~~~~~~~

.. code-block:: python
    :emphasize-lines: 21,24,25,28,29,30,36,37

    import tensorflow as tf
    import horovod.tensorflow.keras as hvd

    hvd.init()

    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.gpu_options.visible_device_list = str(hvd.local_rank())
    tf.keras.backend.set_session(tf.Session(config=config))

    dataset = ...
    model = ...

    opt = keras.optimizers.Adadelta(lr * hvd.size())
    opt = hvd.DistributedOptimizer(opt)

    model.compile(loss=keras.losses.sparse_categorical_crossentropy,
                  optimizer=opt,
                  metrics=['accuracy'])

    def on_state_reset():
        tf.keras.backend.set_value(model.optimizer.lr, lr * hvd.size())

    state = hvd.elastic.KerasState(model, batch=100, epoch=0)
    state.register_reset_callbacks([on_state_reset])

    callbacks = [
        hvd.elastic.CommitStateCallback(state),
        hvd.elastic.UpdateBatchStateCallback(state),
        hvd.elastic.UpdateEpochStateCallback(state),
    ]

    if hvd.rank() == 0:
        callbacks.append(keras.callbacks.ModelCheckpoint('./checkpoint-{epoch}.h5'))

    @hvd.elastic.run
    def train(state):
        model.fit(dataset,
                  steps_per_epoch=500 // hvd.size(),
                  callbacks=callbacks,
                  epochs=epochs - state.epoch,
                  verbose=1 if hvd.rank() == 0 else 0)

    train(state)


Elastic PyTorch
~~~~~~~~~~~~~~~

.. code-block:: python
    :emphasize-lines: 14,15,28,31,36,37

    import torch
    import horovod.torch as hvd

    hvd.init()

    torch.cuda.set_device(hvd.local_rank())

    dataset = ...
    model = ...

    optimizer = optim.SGD(model.parameters(), lr * hvd.size())
    optimizer = hvd.DistributedOptimizer(optimizer)

    @hvd.elastic.run
    def train(state):
        batch_offset = state.batch
        for state.epoch in range(state.epoch, epochs):
            for state.batch in range(state.batch, batches_per_epoch):
                data, target = get_random_batch()

                optimizer.zero_grad()
                output = model(data)
                loss = F.nll_loss(output, target)
                loss.backward()
                optimizer.step()

                if state.batch % batches_per_commit == 0:
                    state.commit()
            state.batch = 0

    def on_state_reset():
        # adjust learning rate on reset
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr * hvd.size()

    state = hvd.elastic.TorchState(model, optimizer, batch=0, epoch=0)
    state.register_reset_callbacks([on_state_reset])
    train(state)


Running with horovodrun
~~~~~~~~~~~~~~~~~~~~~~~

Elastic training jobs are started using the ``horovodrun`` command line tool. The major difference when launching
elastic jobs is that hosts are not specified explicitly, but instead **discovered** at runtime.  The most general way
to allow Horovod to discover available hosts is to provide a ``--host-discovery-script`` when launching the job:

.. code-block:: bash

    $ horovodrun -np 8 --host-discovery-script discover_hosts.sh python train.py

The host discovery script must have user executable permissions, and return one host with its available slots per line
of the form: ``<hostname>:<slots>``.  For example:

.. code-block:: bash

    $ ./discover_hosts.sh
    host-1:4
    host-2:4
    host-3:4

If the host discovery scripts fails to execute (due to a permissions issue) or otherwise returns a non-zero exit code
the first time it is called, the training process will fail immediately. However, subsequent errors will result in
retries until the job times-out (due to failure to discover a sufficient number of slots).

Your discovery script may omit the ``:<slots>`` if you explicitly specify the number of slots per host as an argument:

.. code-block:: bash

    $ horovodrun -np 8 --host-discovery-script discover_hosts.sh --slots 4 python train.py

The elastic training job will not start until at least ``-np`` slots are available for running worker processes.

You can additionally specify the minimum and maximum number of processes to run with during the job:

.. code-block:: bash

    $ horovodrun -np 8 --min-np 4 --max-np 12 --host-discovery-script discover_hosts.sh python train.py

If the number of available slots falls below ``--min-np`` (due to host failure, preemption, etc.), then the job will
pause waiting for more hosts to become available or until ``HOROVOD_ELASTIC_TIMEOUT`` (default: 600 seconds) has
elapsed.  If unspecified, minimum np defaults to ``-np``.

The maximum np can be used to cap the number of processes (to prevent over-utilizing available resources) and to serve
as a reference point for learning rate scales and data partitions (in cases where these need to be held constant
regardless of the current number of workers).  If unspecified, maximum np also defaults to ``-np``.

Instances that fail will be added to a blacklist, as they may have faulty hardware. Hosts will remain in blacklist for a configured cooldown period.
After the cooldown period ends, the hosts will be whitelisted back. This is to account for transient failures, and cases where the same host
is added back to a job.
Cooldown periods can be configured with the ``--blacklist-cooldown-range`` parameter like this:

.. code-block:: bash

    $ horovodrun -np 8 --blacklist-cooldown-range 10 100 --min-np 4 --max-np 12 --host-discovery-script discover_hosts.py python train.py

The above example configures the minimum cooldown period to 10 seconds and the maximum cooldown period to 100 seconds.
The intial cooldown period would be 10 seconds. For repeat failures the cooldown period would grow with an exponential
backoff delay (with a constant exponent of 2): 10s, 20s, 40s, and so on. However, the maximum cooldown period would be
capped at 100 seconds, regardless of failure count. A random backoff fraction of the cooldown lower limit is added
to the cooldown delay.
The default behavior is to have no cooldown period, and blacklisted hosts would remain in blacklist.

Ranks that fail repeatedly will result in job failure, as it may be the case that the training process cannot make progress.


Running on Ray
~~~~~~~~~~~~~~

Running an elastic training script with Ray is simple and provides additional benefits to existing Horovod Elastic functionality:

* You can execute training from interactive Python environments (i.e., a Jupyter notebook)
* You can automatically leverage Ray's autoscaler to add/remove spot instances on AWS/GCP/Azure/Kubernetes.


To use elastic training with Ray:

.. code-block:: python

    import horovod.torch as hvd

    # Put the Horovod concepts into a single function
    # This function will be serialized with Cloudpickle
    def training_fn():
        hvd.init()
        model = Model()
        torch.cuda.set_device(hvd.local_rank())

        @hvd.elastic.run
        def train(state):
            for state.epoch in range(state.epoch, epochs):
                ...
                state.commit()


        state = hvd.elastic.TorchState(model, optimizer, batch=0, epoch=0)
        state.register_reset_callbacks([on_state_reset])
        train(state)
        return


    from horovod.ray import ElasticRayExecutor
    import ray

    ray.init()  # or ray.init(address="auto") if on a Ray cluster

    settings = ElasticRayExecutor.create_settings(verbose=True)
    executor = ElasticRayExecutor(settings, use_gpu=True, cpus_per_slot=2)
    executor.start()
    executor.run(training_fn)


Running on Spark
~~~~~~~~~~~~~~~~

Current constraints:

- `max_np` and `min_np` are `None` or equal to `num_np`, i.e. no auto-scaling, only fault tolerant


Practical Considerations: Consistent training
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With workers frequently being added and removed from the training process, it creates the possibility for learning
rates, numbers of partitions, and other parameters that vary with the number of workers to hurt model convergence if
not properly handled.

Learning rate will need to be rescaled via callback when using gradient averaging.  Using Adasum, no adjustment will
need to be made assuming that local size remains the same.

If using random sampling to read data, then no repartitioning need occur. For the time being, this is the recommended
strategy to simplify elastic training configuration.

If using dataset partitioning, callbacks may be used to repartition dataset as necessary, skipping already processed
data. Care needs to be taken when partitioning the data to ensure that data is not processed more than once. As such,
the preferred approach is to keep the number of partitions constant (from ``hvd.max_size()``), but redistribute
partitions and use local gradient aggregation to keep total batch size constant.

.. inclusion-marker-end-do-not-remove
.. raw:: html

    <p align="center"><img src="https://user-images.githubusercontent.com/16640218/34506318-84d0c06c-efe0-11e7-8831-0425772ed8f2.png" alt="Logo" width="200"/></p>
    <br/>

Horovod
=======

.. raw:: html

   <div align="center">

.. image:: https://badge.fury.io/py/horovod.svg
   :target: https://badge.fury.io/py/horovod
   :alt: PyPI Version

.. image:: https://badge.buildkite.com/6f976bc161c69d9960fc00de01b69deb6199b25680a09e5e26.svg?branch=master
   :target: https://buildkite.com/horovod/horovod
   :alt: Build Status

.. image:: https://readthedocs.org/projects/horovod/badge/?version=latest
   :target: https://horovod.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/badge/slack-chat-green.svg?logo=slack
   :target: https://forms.gle/cPGvty5hp31tGfg79
   :alt: Slack

.. raw:: html

   </div>

.. raw:: html

   <div align="center">

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :target: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :alt: License

.. image:: https://app.fossa.com/api/projects/git%2Bgithub.com%2Fhorovod%2Fhorovod.svg?type=shield
   :target: https://app.fossa.com/projects/git%2Bgithub.com%2Fhorovod%2Fhorovod?ref=badge_shield
   :alt: FOSSA Status

.. image:: https://bestpractices.coreinfrastructure.org/projects/2373/badge
   :target: https://bestpractices.coreinfrastructure.org/projects/2373
   :alt: CII Best Practices

.. image:: https://pepy.tech/badge/horovod
   :target: https://pepy.tech/project/horovod
   :alt: Downloads

.. raw:: html

   </div>

.. inclusion-marker-start-do-not-remove

|

Horovod is a distributed deep learning training framework for TensorFlow, Keras, PyTorch, and Apache MXNet.
The goal of Horovod is to make distributed deep learning fast and easy to use.


.. raw:: html

   <p><img src="https://raw.githubusercontent.com/lfai/artwork/master/lfaidata-assets/lfaidata-project-badge/graduate/color/lfaidata-project-badge-graduate-color.png" alt="LF AI & Data" width="200"/></p>


Horovod is hosted by the `LF AI & Data Foundation <https://lfdl.io>`_ (LF AI & Data). If you are a company that is deeply
committed to using open source technologies in artificial intelligence, machine, and deep learning, and want to support
the communities of open source projects in these domains, consider joining the LF AI & Data Foundation. For details
about who's involved and how Horovod plays a role, read the Linux Foundation `announcement <https://lfdl.io/press/2018/12/13/lf-deep-learning-welcomes-horovod-distributed-training-framework-as-newest-project/>`_.

|

.. contents::

|

Why Horovod?
------------
The primary motivation for this project is to make it easy to take a single-GPU training script and successfully scale
it to train across many GPUs in parallel. This has two aspects:

1. How much modification does one have to make to a program to make it distributed, and how easy is it to run it?
2. How much faster would it run in distributed mode?

Internally at Uber we found the MPI model to be much more straightforward and require far less code changes than previous
solutions such as Distributed TensorFlow with parameter servers. Once a training script has been written for scale with
Horovod, it can run on a single-GPU, multiple-GPUs, or even multiple hosts without any further code changes.
See the `Usage <#usage>`__ section for more details.

In addition to being easy to use, Horovod is fast. Below is a chart representing the benchmark that was done on 128
servers with 4 Pascal GPUs each connected by RoCE-capable 25 Gbit/s network:

.. image:: https://user-images.githubusercontent.com/16640218/38965607-bf5c46ca-4332-11e8-895a-b9c137e86013.png
   :alt: 512-GPU Benchmark

Horovod achieves 90% scaling efficiency for both Inception V3 and ResNet-101, and 68% scaling efficiency for VGG-16.
See `Benchmarks <benchmarks.rst>`_ to find out how to reproduce these numbers.

While installing MPI and NCCL itself may seem like an extra hassle, it only needs to be done once by the team dealing
with infrastructure, while everyone else in the company who builds the models can enjoy the simplicity of training them at
scale.


Install
-------
To install Horovod on Linux or macOS:

1. Install `CMake <https://cmake.org/install/>`__

.. raw:: html

    <p/>

2. If you've installed TensorFlow from `PyPI <https://pypi.org/project/tensorflow>`__, make sure that ``g++-5`` or above is installed.

   If you've installed PyTorch from `PyPI <https://pypi.org/project/torch>`__, make sure that ``g++-5`` or above is installed.

   If you've installed either package from `Conda <https://conda.io>`_, make sure that the ``gxx_linux-64`` Conda package is installed.

.. raw:: html

    <p/>

3. Install the ``horovod`` pip package.

   To run on CPUs:

   .. code-block:: bash

      $ pip install horovod

   To run on GPUs with NCCL:

   .. code-block:: bash

      $ HOROVOD_GPU_OPERATIONS=NCCL pip install horovod

For more details on installing Horovod with GPU support, read `Horovod on GPU <gpus.rst>`_.

For the full list of Horovod installation options, read the `Installation Guide <install.rst>`_.

If you want to use MPI, read `Horovod with MPI <mpi.rst>`_.

If you want to use Conda, read `Building a Conda environment with GPU support for Horovod <conda.rst>`_.

If you want to use Docker, read `Horovod in Docker <docker.rst>`_.

To compile Horovod from source, follow the instructions in the `Contributor Guide <contributors.rst>`_.


Concepts
--------
Horovod core principles are based on `MPI <http://mpi-forum.org/>`_ concepts such as *size*, *rank*,
*local rank*, **allreduce**, **allgather**, **broadcast**, and **alltoall**. See `this page <concepts.rst>`_
for more details.

Supported frameworks
--------------------
See these pages for Horovod examples and best practices:

- `Horovod with TensorFlow <tensorflow.rst>`_
- `Horovod with XLA in Tensorflow <xla.rst>`_
- `Horovod with Keras <keras.rst>`_
- `Horovod with PyTorch <pytorch.rst>`_
- `Horovod with MXNet <mxnet.rst>`_

Usage
-----

To use Horovod, make the following additions to your program:

1. Run ``hvd.init()`` to initialize Horovod.

.. raw:: html

    <p/>

2. Pin each GPU to a single process to avoid resource contention.

   With the typical setup of one GPU per process, set this to *local rank*. The first process on
   the server will be allocated the first GPU, the second process will be allocated the second GPU, and so forth.

.. raw:: html

    <p/>


3. Scale the learning rate by the number of workers.

   Effective batch size in synchronous distributed training is scaled by the number of workers.
   An increase in learning rate compensates for the increased batch size.

.. raw:: html

    <p/>


4. Wrap the optimizer in ``hvd.DistributedOptimizer``.

   The distributed optimizer delegates gradient computation to the original optimizer, averages gradients using **allreduce** or **allgather**, and then applies those averaged gradients.

.. raw:: html

    <p/>


5. Broadcast the initial variable states from rank 0 to all other processes.

   This is necessary to ensure consistent initialization of all workers when training is started with random weights or restored from a checkpoint.

.. raw:: html

    <p/>


6. Modify your code to save checkpoints only on worker 0 to prevent other workers from corrupting them.

.. raw:: html

    <p/>


Example using TensorFlow v1 (see the `examples <https://github.com/horovod/horovod/blob/master/examples/>`_ directory for full training examples):

.. code-block:: python

    import tensorflow as tf
    import horovod.tensorflow as hvd


    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    config = tf.ConfigProto()
    config.gpu_options.visible_device_list = str(hvd.local_rank())

    # Build model...
    loss = ...
    opt = tf.train.AdagradOptimizer(0.01 * hvd.size())

    # Add Horovod Distributed Optimizer
    opt = hvd.DistributedOptimizer(opt)

    # Add hook to broadcast variables from rank 0 to all other processes during
    # initialization.
    hooks = [hvd.BroadcastGlobalVariablesHook(0)]

    # Make training operation
    train_op = opt.minimize(loss)

    # Save checkpoints only on worker 0 to prevent other workers from corrupting them.
    checkpoint_dir = '/tmp/train_logs' if hvd.rank() == 0 else None

    # The MonitoredTrainingSession takes care of session initialization,
    # restoring from a checkpoint, saving to a checkpoint, and closing when done
    # or an error occurs.
    with tf.train.MonitoredTrainingSession(checkpoint_dir=checkpoint_dir,
                                           config=config,
                                           hooks=hooks) as mon_sess:
      while not mon_sess.should_stop():
        # Perform synchronous training.
        mon_sess.run(train_op)


Running Horovod
---------------
The example commands below show how to run distributed training.
See `Run Horovod <running.rst>`_ for more details, including RoCE/InfiniBand tweaks and tips for dealing with hangs.

1. To run on a machine with 4 GPUs:

   .. code-block:: bash

        $ horovodrun -np 4 -H localhost:4 python train.py

2. To run on 4 machines with 4 GPUs each:

   .. code-block:: bash

       $ horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 python train.py

3. To run using Open MPI without the ``horovodrun`` wrapper, see `Running Horovod with Open MPI <mpi.rst>`_.

4. To run in Docker, see `Horovod in Docker <docker.rst>`_.

5. To run on Kubernetes, see `Kubeflow MPI Operator <https://github.com/kubeflow/mpi-operator/>`_, `Helm Chart <https://github.com/kubernetes/charts/tree/master/stable/horovod/>`_, `FfDL <https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/>`_, and `Polyaxon <https://docs.polyaxon.com/integrations/horovod/>`_.

6. To run on Spark, see `Horovod on Spark <spark.rst>`_.

7. To run on Ray, see `Horovod on Ray <ray.rst>`_.

8. To run in Singularity, see `Singularity <https://github.com/sylabs/examples/tree/master/machinelearning/horovod>`_.

9. To run in a LSF HPC cluster (e.g. Summit), see `LSF <lsf.rst>`_.

10. To run on Hadoop Yarn, see `TonY <https://github.com/linkedin/TonY/>`_.

Gloo
----
`Gloo <https://github.com/facebookincubator/gloo>`_ is an open source collective communications library developed by Facebook.

Gloo comes included with Horovod, and allows users to run Horovod without requiring MPI to be installed.

For environments that have support both MPI and Gloo, you can choose to use Gloo at runtime by passing the ``--gloo`` argument to ``horovodrun``:

.. code-block:: bash

     $ horovodrun --gloo -np 2 python train.py

mpi4py
------
Horovod supports mixing and matching Horovod collectives with other MPI libraries, such as `mpi4py <https://mpi4py.scipy.org>`_,
provided that the MPI was built with multi-threading support.

You can check for MPI multi-threading support by querying the ``hvd.mpi_threads_supported()`` function.

.. code-block:: python

    import horovod.tensorflow as hvd

    # Initialize Horovod
    hvd.init()

    # Verify that MPI multi-threading is supported.
    assert hvd.mpi_threads_supported()

    from mpi4py import MPI
    assert hvd.size() == MPI.COMM_WORLD.Get_size()

You can also initialize Horovod with an `mpi4py` sub-communicator, in which case each sub-communicator
will run an independent Horovod training.

.. code-block:: python

    from mpi4py import MPI
    import horovod.tensorflow as hvd

    # Split COMM_WORLD into subcommunicators
    subcomm = MPI.COMM_WORLD.Split(color=MPI.COMM_WORLD.rank % 2,
                                   key=MPI.COMM_WORLD.rank)

    # Initialize Horovod
    hvd.init(comm=subcomm)

    print('COMM_WORLD rank: %d, Horovod rank: %d' % (MPI.COMM_WORLD.rank, hvd.rank()))


Inference
---------
Learn how to optimize your model for inference and remove Horovod operations from the graph `here <inference.rst>`_.


Tensor Fusion
-------------
One of the unique things about Horovod is its ability to interleave communication and computation coupled with the ability
to batch small **allreduce** operations, which results in improved performance. We call this batching feature Tensor Fusion.

See `here <tensor-fusion.rst>`__ for full details and tweaking instructions.


Horovod Timeline
----------------
Horovod has the ability to record the timeline of its activity, called Horovod Timeline.

.. image:: https://user-images.githubusercontent.com/16640218/29735271-9e148da0-89ac-11e7-9ae0-11d7a099ac89.png
   :alt: Horovod Timeline

Use Horovod timeline to analyze Horovod performance.
See `here <timeline.rst>`__ for full details and usage instructions.


Automated Performance Tuning
----------------------------
Selecting the right values to efficiently make use of Tensor Fusion and other advanced Horovod features can involve
a good amount of trial and error. We provide a system to automate this performance optimization process called
**autotuning**, which you can enable with a single command line argument to ``horovodrun``.

See `here <autotune.rst>`__ for full details and usage instructions.


Horovod Process Sets
--------------------
Horovod allows you to concurrently run distinct collective operations in different groups of processes taking part in
one distributed training. Set up ``hvd.process_set`` objects to make use of this capability.

See `Process Sets <process_set.rst>`__ for detailed instructions.


Guides
------
1. Run distributed training in Microsoft Azure using `Batch AI and Horovod <https://github.com/Azure/BatchAI/tree/master/recipes/Horovod>`_.
2. `Distributed model training using Horovod <https://spell.ml/blog/distributed-model-training-using-horovod-XvqEGRUAACgAa5th>`_.

Send us links to any user guides you want to publish on this site

Troubleshooting
---------------
See `Troubleshooting <troubleshooting.rst>`_ and submit a `ticket <https://github.com/horovod/horovod/issues/new>`_
if you can't find an answer.


Citation
--------
Please cite Horovod in your publications if it helps your research:

::

    @article{sergeev2018horovod,
      Author = {Alexander Sergeev and Mike Del Balso},
      Journal = {arXiv preprint arXiv:1802.05799},
      Title = {Horovod: fast and easy distributed deep learning in {TensorFlow}},
      Year = {2018}
    }


Publications
------------
1. Sergeev, A., Del Balso, M. (2017) *Meet Horovod: Uber’s Open Source Distributed Deep Learning Framework for TensorFlow*.
Retrieved from `https://eng.uber.com/horovod/ <https://eng.uber.com/horovod/>`_

2. Sergeev, A. (2017) *Horovod - Distributed TensorFlow Made Easy*. Retrieved from
`https://www.slideshare.net/AlexanderSergeev4/horovod-distributed-tensorflow-made-easy <https://www.slideshare.net/AlexanderSergeev4/horovod-distributed-tensorflow-made-easy>`_

3. Sergeev, A., Del Balso, M. (2018) *Horovod: fast and easy distributed deep learning in TensorFlow*. Retrieved from
`arXiv:1802.05799 <https://arxiv.org/abs/1802.05799>`_


References
----------
The Horovod source code was based off the Baidu `tensorflow-allreduce <https://github.com/baidu-research/tensorflow-allreduce>`_
repository written by Andrew Gibiansky and Joel Hestness. Their original work is described in the article
`Bringing HPC Techniques to Deep Learning <http://andrew.gibiansky.com/blog/machine-learning/baidu-allreduce/>`_.

Getting Involved
----------------
- `Community Slack <https://forms.gle/cPGvty5hp31tGfg79>`_ for collaboration and discussion
- `Horovod Announce <https://lists.lfai.foundation/g/horovod-announce>`_ for updates on the project
- `Horovod Technical-Discuss <https://lists.lfai.foundation/g/horovod-technical-discuss>`_ for public discussion
- `Horovod Security <https://lists.lfai.foundation/g/horovod-security>`_ to report security vulnerabilities


.. inclusion-marker-end-do-not-remove
   Place contents above here if they should also appear in read-the-docs.
   Contents below are already part of the read-the-docs table of contents.
.. include:: ./conda.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove.. inclusion-marker-start-do-not-remove

Horovod with MPI
================

MPI can be used as an alternative to Gloo for coordinating work between processes in Horovod. When using NCCL, performance
will be similar between the two, but if you are doing CPU training, there are noticeable performance benefits to using MPI.

First install `Open MPI <https://www.open-mpi.org/>`_ or another MPI implementation. Learn how to install Open MPI `on this page <https://www.open-mpi.org/faq/?category=building#easy-build>`_.

**Note**: Open MPI 3.1.3 has an issue that may cause hangs. The recommended fix is to downgrade to Open MPI 3.1.2 or upgrade to Open MPI 4.0.0.

mpirun
------

``horovodrun`` introduces a convenient, Open MPI-based wrapper for running Horovod scripts.

In some cases it is desirable to have fine-grained control over options passed to Open MPI.  This page describes
running Horovod training directly using Open MPI.

1. Run on a machine with 4 GPUs:

   .. code-block:: bash

       horovodrun -np 4 python train.py

   Equivalent Open MPI command:

   .. code-block:: bash

       mpirun -np 4 \
           -bind-to none -map-by slot \
           -x NCCL_DEBUG=INFO -x LD_LIBRARY_PATH -x PATH \
           -mca pml ob1 -mca btl ^openib \
           python train.py

2. Run on 4 machines with 4 GPUs each:

   .. code-block:: bash

      horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 python train.py

   Equivalent Open MPI command:

   .. code-block:: bash

       mpirun -np 16 \
           -H server1:4,server2:4,server3:4,server4:4 \
           -bind-to none -map-by slot \
           -x NCCL_DEBUG=INFO -x LD_LIBRARY_PATH -x PATH \
           -mca pml ob1 -mca btl ^openib \
           python train.py

Starting with the Open MPI 3, it's important to add the ``-bind-to none`` and ``-map-by slot`` arguments.
``-bind-to none`` specifies Open MPI to not bind a training process to a single CPU core (which would hurt performance).
``-map-by slot`` allows you to have a mixture of different NUMA configurations because the default behavior is to bind
to the socket.

The ``-mca pml ob1`` and ``-mca btl ^openib`` flags force the use of TCP for MPI communication.  This avoids many
multiprocessing issues that Open MPI has with RDMA which typically results in segmentation faults.  Using TCP for MPI
does not have noticeable performance impact since most of the heavy communication is done by NCCL, which will use RDMA
via RoCE or InfiniBand if they're available (see `Horovod on GPU <gpus.rst>`_).  Notable exceptions from this rule are
models that heavily use ``hvd.broadcast()`` and ``hvd.allgather()`` operations.  To make those operations use RDMA,
read the `Open MPI with RDMA <#open-mpi-with-rdma>`_ section below.

With the ``-x`` option you can specify (``-x NCCL_DEBUG=INFO``) or copy (``-x LD_LIBRARY_PATH``) an environment variable to
all the workers.

Custom SSH ports
~~~~~~~~~~~~~~~~

Specify custom SSH ports with ``-mca plm_rsh_args "-p <port>"`` as follows:

.. code-block:: bash

    mpirun -np 16 \
        -H server1:4,server2:4,server3:4,server4:4 \
        -bind-to none -map-by slot \
        -mca plm_rsh_args "-p 12345"
        -x NCCL_DEBUG=INFO -x LD_LIBRARY_PATH -x PATH \
        -mca pml ob1 -mca btl ^openib \
        python train.py

This is frequently useful in the case of `running Horovod in Docker environment <docker.rst>`_.

Open MPI with RDMA
~~~~~~~~~~~~~~~~~~

As noted above, using TCP for MPI communication does not have any significant effects on performance in the majority of
cases. Models that make heavy use of ``hvd.broadcast()`` and ``hvd.allgather()`` operations are exceptions to that rule.

Default Open MPI ``openib`` BTL that provides RDMA functionality does not work well with MPI multithreading.  In order
to use RDMA with ``openib``, multithreading must be disabled via the ``-x HOROVOD_MPI_THREADS_DISABLE=1`` option.  See the
example below:

.. code-block:: bash

    mpirun -np 16 \
        -H server1:4,server2:4,server3:4,server4:4 \
        -bind-to none -map-by slot \
        -x NCCL_DEBUG=INFO -x LD_LIBRARY_PATH -x HOROVOD_MPI_THREADS_DISABLE=1 -x PATH \
        -mca pml ob1 \
        python train.py

Other MPI RDMA implementations may or may not benefit from disabling multithreading, so please consult vendor
documentation.

Horovod Parameter Knobs
~~~~~~~~~~~~~~~~~~~~~~~

Many of the configurable parameters available as command line arguments to ``horovodrun`` can be used with ``mpirun``
through the use of environment variables.

Tensor Fusion:

.. code-block:: bash

    $ mpirun -x HOROVOD_FUSION_THRESHOLD=33554432 -x HOROVOD_CYCLE_TIME=3.5 ... python train.py

Timeline:

.. code-block:: bash

    $ mpirun -x HOROVOD_TIMELINE=/path/to/timeline.json -x HOROVOD_TIMELINE_MARK_CYCLES=1 ... python train.py

Autotuning:

.. code-block:: bash

    $ mpirun -x HOROVOD_AUTOTUNE=1 -x HOROVOD_AUTOTUNE_LOG=/tmp/autotune_log.csv ... python train.py

Note that when using ``horovodrun``, any command line arguments will override values set in the environment.

Hangs due to non-routed network interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Having network interfaces that are not routed can cause Open MPI to hang. An example of such interface is ``docker0``.

If you see non-routed interfaces (like ``docker0``) in the output of ``ifconfig``, you should tell Open MPI and NCCL to not
use them via the ``-mca btl_tcp_if_exclude <interface>[,<interface>]`` and ``NCCL_SOCKET_IFNAME=^<interface>[,<interface>]``
parameters.

.. code-block:: bash

    ifconfig

Produces output like this::

    docker0   Link encap:Ethernet  HWaddr 02:42:2d:17:ea:66
              inet addr:172.17.0.1  Bcast:0.0.0.0  Mask:255.255.0.0
              UP BROADCAST MULTICAST  MTU:1500  Metric:1
              RX packets:0 errors:0 dropped:0 overruns:0 frame:0
              TX packets:0 errors:0 dropped:0 overruns:0 carrier:0
              collisions:0 txqueuelen:0
              RX bytes:0 (0.0 B)  TX bytes:0 (0.0 B)
    eth0      Link encap:Ethernet  HWaddr 24:8a:07:b3:7d:8b
              inet addr:10.0.0.1  Bcast:10.0.0.255  Mask:255.255.255.0
              UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
              RX packets:900002410 errors:0 dropped:405 overruns:0 frame:0
              TX packets:1521598641 errors:0 dropped:0 overruns:0 carrier:0
              collisions:0 txqueuelen:1000
              RX bytes:376184431726 (350.3 GiB)  TX bytes:954933846124 (889.3 GiB)
    eth1      Link encap:Ethernet  HWaddr 24:8a:07:b3:7d:8a
              inet addr:192.168.0.1  Bcast:192.168.0.255  Mask:255.255.255.0
              UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
              RX packets:2410141 errors:0 dropped:0 overruns:0 frame:0
              TX packets:2312177 errors:0 dropped:0 overruns:0 carrier:0
              collisions:0 txqueuelen:1000
              RX bytes:698398061 (666.0 MiB)  TX bytes:458504418 (437.2 MiB)
    lo        Link encap:Local Loopback
              inet addr:127.0.0.1  Mask:255.0.0.0
              inet6 addr: ::1/128 Scope:Host
              UP LOOPBACK RUNNING  MTU:65536  Metric:1
              RX packets:497075633 errors:0 dropped:0 overruns:0 frame:0
              TX packets:497075633 errors:0 dropped:0 overruns:0 carrier:0
              collisions:0 txqueuelen:1
              RX bytes:72680421398 (67.6 GiB)  TX bytes:72680421398 (67.6 GiB)

Example ``mpirun`` command with ``lo`` and ``docker0`` interfaces excluded:

.. code-block:: bash

    mpirun -np 16 \
        -H server1:4,server2:4,server3:4,server4:4 \
        -bind-to none -map-by slot \
        -x NCCL_DEBUG=INFO -x LD_LIBRARY_PATH -x PATH \
        -x NCCL_SOCKET_IFNAME=^lo,docker0 \
        -mca pml ob1 -mca btl ^openib \
        -mca btl_tcp_if_exclude lo,docker0 \
        python train.py

.. inclusion-marker-end-do-not-remove
.. include:: ./spark.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Horovod on Ray
==============

``horovod.ray`` allows users to leverage Horovod on `a Ray cluster <https://docs.ray.io/en/latest/cluster/index.html>`_.

Currently, the Ray + Horovod integration provides a :ref:`RayExecutor API <horovod_ray_api>`.

.. note:: The Ray + Horovod integration currently only supports a Gloo backend.

Installation
------------

Use the extra ``[ray]`` option to install Ray along with Horovod.

.. code-block:: bash

    $ HOROVOD_WITH_GLOO=1 ... pip install 'horovod[ray]'

See the Ray documentation for `advanced installation instructions <https://docs.ray.io/en/latest/installation.html>`_.


Horovod Ray Executor
--------------------

The Horovod Ray integration offers a ``RayExecutor`` abstraction (:ref:`docs <horovod_ray_api>`),
which is a wrapper over a group of `Ray actors (stateful processes) <https://docs.ray.io/en/latest/walkthrough.html#remote-classes-actors>`_.

.. code-block:: python

    from horovod.ray import RayExecutor

    # Start the Ray cluster or attach to an existing Ray cluster
    ray.init()

    # Start num_workers actors on the cluster
    executor = RayExecutor(
        setting, num_workers=num_workers, use_gpu=True)

    # This will launch `num_workers` actors on the Ray Cluster.
    executor.start()


All actors will be part of the Horovod ring, so ``RayExecutor`` invocations will be able to support arbitrary Horovod collective operations.

Note that there is an implicit assumption on the cluster being homogenous in shape (i.e., all machines have the same number of slots available). This is simply
an implementation detail and is not a fundamental limitation.

To actually execute a function, you can run the following:

.. code-block:: python

    # Using the stateless `run` method, a function can take in any args or kwargs
    def simple_fn():
        hvd.init()
        print("hvd rank", hvd.rank())
        return hvd.rank()

    # Execute the function on all workers at once
    result = executor.run(simple_fn)
    # Check that the rank of all workers is unique
    assert len(set(result)) == hosts * num_slots

    executor.shutdown()


Stateful Execution
~~~~~~~~~~~~~~~~~~

A unique feature of Ray is its support for `stateful Actors <https://docs.ray.io/en/latest/walkthrough.html#remote-classes-actors>`_. This means that you can start arbitrary Python classes on each worker, easily supporting operations and calls where data is cached in memory.

.. code-block:: python

    import torch
    from horovod.torch import hvd
    from horovod.ray import RayExecutor

    class MyModel:
        def __init__(self, learning_rate):
            self.model = NeuralNet()
            optimizer = torch.optim.SGD(
                self.model.parameters(),
                lr=learning_rate,
            )
            self.optimizer = hvd.DistributedOptimizer(optimizer)

        def get_weights(self):
            return dict(self.model.parameters())

        def train(self):
            return self._train(self.model, self.optimizer)


    ray.init()
    executor = RayExecutor(...)
    executor.start(executable_cls=MyModel)

    # Run 5 training steps
    for i in range(5):
        # Stateful `execute` method takes the current worker executable as a parameter
        executor.execute(lambda worker: worker.train())

    # Obtain the trained weights from each model replica
    result = executor.execute(lambda worker: worker.get_weights())

    # `result` will be N copies of the model weights
    assert all(isinstance(res, dict) for res in result)

Elastic Ray Executor
--------------------

Ray also supports `elastic execution <elastic.rst>`_ via :ref:`the RayExecutor <horovod_ray_api>`. Similar to default Horovod, the difference between the non-elastic and elastic versions of Ray is that the hosts and number of workers is dynamically determined at runtime.

You must first set up `a Ray cluster`_. Ray clusters can support autoscaling for any cloud provider (AWS, GCP, Azure).

.. code-block:: bash

    # First, run `pip install boto3` and `aws configure`
    #
    # Create or update the cluster. When the command finishes, it will print
    # out the command that can be used to SSH into the cluster head node.
    $ ray up ray/python/ray/autoscaler/aws/example-full.yaml

After you have a Ray cluster setup, you will need to move parts of your existing elastic Horovod training script into a training function. Specifically,
the instantiation of your model and the invocation of the ``hvd.elastic.run`` call should be done inside this function.

.. code-block:: python

    import horovod.torch as hvd

    # Put the Horovod concepts into a single function
    # This function will be serialized with Cloudpickle
    def training_fn():
        hvd.init()
        model = Model()
        torch.cuda.set_device(hvd.local_rank())

        @hvd.elastic.run
        def train(state):
            for state.epoch in range(state.epoch, epochs):
                ...
                state.commit()


        state = hvd.elastic.TorchState(model, optimizer, batch=0, epoch=0)
        state.register_reset_callbacks([on_state_reset])
        train(state)
        return

You can then attach to the underlying Ray cluster and execute the training function:

.. code-block:: python

    import ray
    from horovod.ray import RayExecutor

    ray.init(address="auto")  # attach to the Ray cluster
    settings = RayExecutor.create_settings(verbose=True)
    executor = RayExecutor(
        settings, min_workers=1, use_gpu=True, cpus_per_slot=2)
    executor.start()
    executor.run(training_fn)

Ray will automatically start remote actors which execute ``training_fn`` on nodes as they become available. Note that ``executor.run`` call will terminate whenever any one of the training functions terminates successfully, or if all workers fail.

AWS: Cluster Launcher
---------------------

You can also easily leverage the `Ray cluster launcher <https://docs.ray.io/en/latest/cluster/launcher.html>`_ to spin up cloud instances.

.. code-block:: yaml

    # Save as `ray_cluster.yaml`

    cluster_name: horovod-cluster
    provider: {type: aws, region: us-west-2}
    auth: {ssh_user: ubuntu}
    min_workers: 3
    max_workers: 3

    # Deep Learning AMI (Ubuntu) Version 21.0
    head_node: {InstanceType: p3.2xlarge, ImageId: ami-0b294f219d14e6a82}
    worker_nodes: {InstanceType: p3.2xlarge, ImageId: ami-0b294f219d14e6a82}
    setup_commands: # Set up each node.
        - HOROVOD_WITH_GLOO=1 HOROVOD_GPU_OPERATIONS=NCCL pip install horovod[ray]

You can start the specified Ray cluster and monitor its status with:

.. code-block:: bash

    $ ray up ray_cluster.yaml  # starts the head node
    $ ray monitor ray_cluster.yaml  # wait for worker nodes

Then, in your python script, make sure you add ``ray.init(address="auto")`` to connect
to the distributed Ray cluster.

.. code-block:: diff

    -ray.init()
    +ray.init(address="auto")

Then you can execute Ray scripts on the cluster:

.. code-block:: bash

    $ ray submit ray_cluster.yaml <your_script.py>

    # the above is is equivalent to
    $ ray attach ray_cluster.yaml  # ssh
    ubuntu@ip-172-31-24-53:~$ python <your_script.py>

.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Tensor Fusion
=============

One of the unique things about Horovod is its ability to interleave communication and computation coupled with the ability
to batch small **allreduce** operations, which results in improved performance. We call this batching feature Tensor Fusion.

Tensor Fusion works by attempting to combine all the tensors that are ready to be reduced at given moment of time into
one reduction operation. The algorithm of Tensor Fusion is as follows:

1. Determine which tensors are ready to be reduced. Select first few tensors that fit in ``HOROVOD_FUSION_THRESHOLD`` bytes and have the same data type.
2. Allocate fusion buffer of size ``HOROVOD_FUSION_THRESHOLD`` if it was not allocated before. Default fusion buffer size is 128 MB.
3. Copy data of selected tensors into the fusion buffer.
4. Execute the **allreduce** operation on the fusion buffer.
5. Copy data from the fusion buffer into the output tensors.
6. Repeat until there are no more tensors to reduce in this cycle.

The fusion buffer size can be adjusted using the ``--fusion-threshold-mb`` command line argument to ``horovodrun``:

.. code-block:: bash

    $ horovodrun -np 4 --fusion-threshold-mb 32 python train.py

Setting ``--fusion-threshold-mb`` to zero disables Tensor Fusion:

.. code-block:: bash

    $ horovodrun -np 4 --fusion-threshold-mb 0 python train.py

You can tweak time between cycles (defined in milliseconds) using the ``--cycle-time-ms`` command line argument:

.. code-block:: bash

    $ horovodrun -np 4 --cycle-time-ms 3.5 python train.py

.. inclusion-marker-end-do-not-remove
Horovod with Keras
==================
Horovod supports Keras and regular TensorFlow in similar ways. To use Horovod with Keras, make the following modifications to your training script:

1. Run ``hvd.init()``.

.. raw:: html

    <p/>

2. Pin each GPU to a single process.

   With the typical setup of one GPU per process, set this to *local rank*. The first process on
   the server will be allocated the first GPU, the second process will be allocated the second GPU, and so forth.

   For **TensorFlow v1**:

   .. code-block:: python

       config = tf.ConfigProto()
       config.gpu_options.visible_device_list = str(hvd.local_rank())
       K.set_session(tf.Session(config=config))

   For **TensorFlow v2**:

   .. code-block:: python

       gpus = tf.config.experimental.list_physical_devices('GPU')
       for gpu in gpus:
           tf.config.experimental.set_memory_growth(gpu, True)
       if gpus:
           tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

.. raw:: html

    <p/>

3. Scale the learning rate by the number of workers.

   Effective batch size in synchronous distributed training is scaled by the number of workers.
   An increase in learning rate compensates for the increased batch size.

.. raw:: html

    <p/>

4. Wrap the optimizer in ``hvd.DistributedOptimizer``.

   The distributed optimizer delegates gradient computation to the original optimizer, averages gradients using *allreduce* or *allgather*, and then applies those averaged gradients.

.. raw:: html

    <p/>

5. Add ``hvd.callbacks.BroadcastGlobalVariablesCallback(0)`` to broadcast initial variable states from rank 0 to all other processes.

   This is necessary to ensure consistent initialization of all workers when training is started with random weights or restored from a checkpoint.

.. raw:: html

    <p/>

6. Modify your code to save checkpoints only on worker 0 to prevent other workers from corrupting them.

   Accomplish this by guarding model checkpointing code with ``hvd.rank() != 0``.

.. raw:: html

    <p/>

.. NOTE:: - Keras 2.0.9 has a `known issue <https://github.com/fchollet/keras/issues/8353>`_ that makes each worker allocate all GPUs on the server, instead of the GPU assigned by the *local rank*. If you have multiple GPUs per server, upgrade to Keras 2.1.2 or downgrade to Keras 2.0.8.

          - To use ``keras`` bundled with ``tensorflow`` you must use ``from tensorflow import keras`` instead of ``import keras`` and ``import horovod.tensorflow.keras as hvd`` instead of ``import horovod.keras as hvd`` in the import statements.

See full training `simple <https://github.com/horovod/horovod/blob/master/examples/keras/keras_mnist.py>`_ (shown below) and `advanced <https://github.com/horovod/horovod/blob/master/examples/keras/keras_mnist_advanced.py>`_ examples.


.. code-block:: python

    from __future__ import print_function
    import keras
    from keras.datasets import mnist
    from keras.models import Sequential
    from keras.layers import Dense, Dropout, Flatten
    from keras.layers import Conv2D, MaxPooling2D
    from keras import backend as K
    import math
    import tensorflow as tf
    import horovod.keras as hvd

    # Horovod: initialize Horovod.
    hvd.init()

    # Horovod: pin GPU to be used to process local rank (one GPU per process)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.gpu_options.visible_device_list = str(hvd.local_rank())
    K.set_session(tf.Session(config=config))

    batch_size = 128
    num_classes = 10

    # Horovod: adjust number of epochs based on number of GPUs.
    epochs = int(math.ceil(12.0 / hvd.size()))

    # Input image dimensions
    img_rows, img_cols = 28, 28

    # The data, shuffled and split between train and test sets
    (x_train, y_train), (x_test, y_test) = mnist.load_data()

    if K.image_data_format() == 'channels_first':
        x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
        x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)
        input_shape = (1, img_rows, img_cols)
    else:
        x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
        x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
        input_shape = (img_rows, img_cols, 1)

    x_train = x_train.astype('float32')
    x_test = x_test.astype('float32')
    x_train /= 255
    x_test /= 255
    print('x_train shape:', x_train.shape)
    print(x_train.shape[0], 'train samples')
    print(x_test.shape[0], 'test samples')

    # Convert class vectors to binary class matrices
    y_train = keras.utils.to_categorical(y_train, num_classes)
    y_test = keras.utils.to_categorical(y_test, num_classes)

    model = Sequential()
    model.add(Conv2D(32, kernel_size=(3, 3),
                    activation='relu',
                    input_shape=input_shape))
    model.add(Conv2D(64, (3, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(num_classes, activation='softmax'))

    # Horovod: adjust learning rate based on number of GPUs.
    opt = keras.optimizers.Adadelta(1.0 * hvd.size())

    # Horovod: add Horovod Distributed Optimizer.
    opt = hvd.DistributedOptimizer(opt)

    model.compile(loss=keras.losses.categorical_crossentropy,
                  optimizer=opt,
                  metrics=['accuracy'])

    callbacks = [
        # Horovod: broadcast initial variable states from rank 0 to all other processes.
        # This is necessary to ensure consistent initialization of all workers when
        # training is started with random weights or restored from a checkpoint.
        hvd.callbacks.BroadcastGlobalVariablesCallback(0),
    ]

    # Horovod: save checkpoints only on worker 0 to prevent other workers from corrupting them.
    if hvd.rank() == 0:
        callbacks.append(keras.callbacks.ModelCheckpoint('./checkpoint-{epoch}.h5'))

    model.fit(x_train, y_train,
              batch_size=batch_size,
              callbacks=callbacks,
              epochs=epochs,
              verbose=1,
              validation_data=(x_test, y_test))
    score = model.evaluate(x_test, y_test, verbose=0)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])

TensorFlow v2 Keras Example (from the `MNIST <https://github.com/horovod/horovod/blob/master/examples/tensorflow2/tensorflow2_keras_mnist.py>`_ example):

.. code-block:: python

    import tensorflow as tf
    import horovod.tensorflow.keras as hvd

    # Initialize Horovod
    hvd.init()

    # Pin GPU to be used to process local rank (one GPU per process)
    gpus = tf.config.experimental.list_physical_devices('GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    if gpus:
        tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

    # Build model and dataset
    dataset = ...
    model = ...
    opt = tf.optimizers.Adam(0.001 * hvd.size())

    # Horovod: add Horovod DistributedOptimizer.
    opt = hvd.DistributedOptimizer(opt)

    # Horovod: Specify `experimental_run_tf_function=False` to ensure TensorFlow
    # uses hvd.DistributedOptimizer() to compute gradients.
    mnist_model.compile(loss=tf.losses.SparseCategoricalCrossentropy(),
                        optimizer=opt,
                        metrics=['accuracy'],
                        experimental_run_tf_function=False)

    callbacks = [
        # Horovod: broadcast initial variable states from rank 0 to all other processes.
        # This is necessary to ensure consistent initialization of all workers when
        # training is started with random weights or restored from a checkpoint.
        hvd.callbacks.BroadcastGlobalVariablesCallback(0),
    ]

    # Horovod: save checkpoints only on worker 0 to prevent other workers from corrupting them.
    if hvd.rank() == 0:
        callbacks.append(keras.callbacks.ModelCheckpoint('./checkpoint-{epoch}.h5'))

    model.fit(dataset,
              steps_per_epoch=500 // hvd.size(),
              callbacks=callbacks,
              epochs=24,
              verbose=1 if hvd.rank() == 0 else 0)
.. include:: ./contributors.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Horovod on Spark
================

The ``horovod.spark`` package provides a convenient wrapper around Horovod that makes running distributed training
jobs in Spark clusters easy.

In situations where training data originates from Spark, this enables
a tight model design loop in which data processing, model training, and
model evaluation are all done in Spark.

We provide two APIs for running Horovod on Spark: a high level **Estimator** API and a lower level **Run** API. Both
use the same underlying mechanism to launch Horovod on Spark executors, but the Estimator API abstracts the data
processing (from Spark DataFrames to deep learning datasets), model training loop, model checkpointing, metrics
collection, and distributed training.

We recommend using Horovod Spark Estimators if you:

* Are using Keras (``tf.keras`` or ``keras``) or PyTorch for training.
* Want to train directly on a Spark DataFrame from ``pyspark``.
* Are using a standard gradient descent optimization process as your training loop.

If for whatever reason the Estimator API does not meet your needs, the Run API offers more fine-grained control.

Installation
------------

When installing Horovod for usage with Spark, use the extra ``[spark]`` to install all Spark dependencies as well:

.. code-block:: bash

    $ ... pip install horovod[spark]

Note that Horovod Spark Estimators require the following:

*  ``horovod >= 0.19.0``
*  ``pyspark >= 2.3.2``

Not included in the list of dependencies by ``[spark]`` are deep learning frameworks (TensorFlow or PyTorch).
Horovod Spark Estimators additionally require at least one of these combinations:

*  ``tensorflow-gpu >= 1.12.0`` or ``tensorflow >= 1.12.0`` (for ``KerasEstimator``)
*  ``torch >= 1.0.0`` and ``tensorboard >= 1.14.0`` (for ``TorchEstimator``)
*  ``torch >= 1.4.0`` and ``pytorch_lightning >= 1.3.8`` (for ``LightningEstimator``)


Horovod Spark Estimators
~~~~~~~~~~~~~~~~~~~~~~~~
Horovod Spark Estimators allow you to train your deep neural network directly on an existing Spark DataFrame,
leveraging Horovod’s ability to scale across multiple workers, without any specialized code for distributed training:

.. code-block:: python

    from tensorflow import keras
    import tensorflow as tf
    import horovod.spark.keras as hvd

    model = keras.models.Sequential()
        .add(keras.layers.Dense(8, input_dim=2))
        .add(keras.layers.Activation('tanh'))
        .add(keras.layers.Dense(1))
        .add(keras.layers.Activation('sigmoid'))

    # NOTE: unscaled learning rate
    optimizer = keras.optimizers.SGD(lr=0.1)
    loss = 'binary_crossentropy'

    store = HDFSStore('/user/username/experiments')
    keras_estimator = hvd.KerasEstimator(
        num_proc=4,
        store=store,
        model=model,
        optimizer=optimizer,
        loss=loss,
        feature_cols=['features'],
        label_cols=['y'],
        batch_size=32,
        epochs=10)


    keras_model = keras_estimator.fit(train_df) \
        .setOutputCols(['predict'])
    predict_df = keras_model.transform(test_df)

The Estimator hides the complexity of gluing Spark DataFrames to a deep learning training script, reading data into a
format interpretable by the training framework, and distributing the training using Horovod.  The user only needs to
provide a Keras or PyTorch model, and the Estimator will do the work of fitting it to the DataFrame.

After training, the Estimator returns a Transformer representation of the trained model.  The model transformer can
be used like any Spark ML transformer to make predictions on an input DataFrame, writing them as new columns in the
output DataFrame.

Estimators can be used to track experiment history through model checkpointing, hot-start retraining, and metric
logging (for Tensorboard) using the Estimator ``Store`` abstraction.  Stores are used for persisting all training
artifacts including intermediate representations of the training data.  Horovod natively supports stores for HDFS
and local filesystems.

`Petastorm <https://github.com/uber/petastorm/blob/master/petastorm/pytorch.py#L258>` based data loader is used by default, but user can define a custom data loader by override the `BaseDataLoader` interface. A async data loader mixin can also added on top of the data loader.

End-to-end example
------------------
`keras_spark_rossmann_estimator.py script <../examples/spark/keras/keras_spark_rossmann_estimator.py>`__ provides
an example of end-to-end data preparation and training of a model for the
`Rossmann Store Sales <https://www.kaggle.com/c/rossmann-store-sales>`__ Kaggle
competition. It is inspired by an article `An Introduction to Deep Learning for Tabular Data <https://www.fast.ai/2018/04/29/categorical-embeddings/>`__
and leverages the code of the notebook referenced in the article. The example is split into three parts:

#. The first part performs complicated data preprocessing over an initial set of CSV files provided by the competition and gathered by the community.
#. The second part defines a Keras model and performs a distributed training of the model using Horovod on Spark.
#. The third part performs prediction using the best model and creates a submission file.

To run the example, be sure to install Horovod with ``[spark]``, then:

.. code-block:: bash

    $ wget https://raw.githubusercontent.com/horovod/horovod/master/examples/spark/keras/keras_spark_rossmann_estimator.py
    $ wget http://files.fast.ai/part2/lesson14/rossmann.tgz
    $ tar zxvf rossmann.tgz
    $ python keras_spark_rossmann_estimator.py

For pytorch, you can check `pytorch_lightning_spark_mnist.py script <../examples/spark/pytorch/pytorch_lightning_spark_mnist.py>`__ for how to use use lightning estimator with horovod backend to train mnist model on spark.

Training on existing Parquet datasets
-------------------------------------

If your data is already in the Parquet format and you wish to train on it with Horovod Spark Estimators, you
can do so without needing to reprocess the data in Spark. Using `Estimator.fit_on_parquet()`, you can train directly
on an existing Parquet dataset:

.. code-block:: python

    store = HDFSStore(train_path='/user/username/training_dataset', val_path='/user/username/val_dataset')
    keras_estimator = hvd.KerasEstimator(
        num_proc=4,
        store=store,
        model=model,
        optimizer=optimizer,
        loss=loss,
        feature_cols=['features'],
        label_cols=['y'],
        batch_size=32,
        epochs=10)

    keras_model = keras_estimator.fit_on_parquet()

The resulting ``keras_model`` can then be used the same way as any Spark Transformer, or you can extract the underlying
Keras model and use it outside of Spark:

.. code-block:: python

    model = keras_model.getModel()
    pred = model.predict([np.ones([1, 2], dtype=np.float32)])

This approach will work on datasets created using ``horovod.spark.common.util.prepare_data``. It will also work with
any Parquet file that contains no Spark user-defined data types (like ``DenseVector`` or ``SparseVector``).  It's
recommended to use ``prepare_data`` to ensure the data is properly prepared for training even if you have an existing
dataset in Parquet format.  Using ``prepare_data`` allows you to properly partition the dataset for the number of
training processes you intend to use, as well as compress large sparse data columns:

.. code-block:: python

    store = HDFSStore(train_path='/user/username/training_dataset', val_path='/user/username/val_dataset')
    with util.prepare_data(num_processes=4,
                           store=store,
                           df=df,
                           feature_columns=['features'],
                           label_columns=['y'],
                           validation=0.1,
                           compress_sparse=True):
        keras_estimator = hvd.KerasEstimator(
            num_proc=4,
            store=store,
            model=model,
            optimizer=optimizer,
            loss=loss,
            feature_cols=['features'],
            label_cols=['y'],
            batch_size=32,
            epochs=10)

        keras_model = keras_estimator.fit_on_parquet()

Once the data has been prepared, you can reuse it in future Spark applications without needing to call
``util.prepare_data`` again.

Horovod Spark Run
~~~~~~~~~~~~~~~~~
You can also use Horovod on Spark to run the same code you would within an ordinary training script using any
framework supported by Horovod.  To do so, simply write your training logic within a function, then use
``horovod.spark.run`` to execute the function in parallel with MPI on top of Spark.

Because Horovod on Spark uses ``cloudpickle`` to send the training function to workers for execution, you can capture
local variables from your training script or notebook within the training function, similar to using a user-defined
function in PySpark.

A toy example of running a Horovod job in Spark is provided below:

.. code-block:: bash

    $ pyspark
    [PySpark welcome message]

    >>> def fn(magic_number):
    ...   import horovod.torch as hvd
    ...   hvd.init()
    ...   print('Hello, rank = %d, local_rank = %d, size = %d, local_size = %d, magic_number = %d' % (hvd.rank(), hvd.local_rank(), hvd.size(), hvd.local_size(), magic_number))
    ...   return hvd.rank()
    ...
    >>> import horovod.spark
    >>> horovod.spark.run(fn, args=(42,))
    Running 16 processes...
    [Stage 0:>                                                        (0 + 16) / 16]
    Hello, rank = 15, local_rank = 3, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 13, local_rank = 1, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 8, local_rank = 0, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 9, local_rank = 1, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 10, local_rank = 2, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 11, local_rank = 3, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 6, local_rank = 2, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 4, local_rank = 0, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 0, local_rank = 0, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 1, local_rank = 1, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 2, local_rank = 2, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 5, local_rank = 1, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 3, local_rank = 3, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 12, local_rank = 0, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 7, local_rank = 3, size = 16, local_size = 4, magic_number = 42
    Hello, rank = 14, local_rank = 2, size = 16, local_size = 4, magic_number = 42
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    >>>

A more complete example can be found in `keras_spark_rossmann_run.py <../examples/spark/keras/keras_spark_rossmann_run.py>`__, which
shows how you can use the low level ``horovod.spark.run`` API to train a model end-to-end in the following steps:

.. code-block:: bash

    $ wget https://raw.githubusercontent.com/horovod/horovod/master/examples/spark/keras/keras_spark_rossmann_run.py
    $ wget http://files.fast.ai/part2/lesson14/rossmann.tgz
    $ tar zxvf rossmann.tgz
    $ python keras_spark_rossmann_run.py


Spark cluster setup
~~~~~~~~~~~~~~~~~~~
As deep learning workloads tend to have very different resource requirements
from typical data processing workloads, there are certain considerations
for DL Spark cluster setup.

GPU training
------------
For GPU training, one approach is to set up a separate GPU Spark cluster
and configure each executor with ``# of CPU cores`` = ``# of GPUs``. This can
be accomplished in standalone mode as follows:

.. code-block:: bash

    $ echo "export SPARK_WORKER_CORES=<# of GPUs>" >> /path/to/spark/conf/spark-env.sh
    $ /path/to/spark/sbin/start-all.sh


This approach turns the ``spark.task.cpus`` setting to control # of GPUs
requested per process (defaults to 1).

The ongoing `SPARK-24615 <https://issues.apache.org/jira/browse/SPARK-24615>`__ effort aims to
introduce GPU-aware resource scheduling in future versions of Spark.

CPU training
------------
For CPU training, one approach is to specify the ``spark.task.cpus`` setting
during the training session creation:

.. code-block:: python

    conf = SparkConf().setAppName('training') \
        .setMaster('spark://training-cluster:7077') \
        .set('spark.task.cpus', '16')
    spark = SparkSession.builder.config(conf=conf).getOrCreate()


This approach allows you to reuse the same Spark cluster for data preparation
and training.

Security
--------
Horovod on Spark uses Open MPI to run the Horovod jobs in Spark, so
it's as secure as the Open MPI implementation itself.

Since Open MPI does not use encrypted communication and is capable of
launching new processes, it's recommended to **use network level
security to isolate Horovod jobs from potential attackers**.

Environment knobs
-----------------
* ``HOROVOD_SPARK_START_TIMEOUT`` - sets the default timeout for Spark tasks to spawn, register, and start running the code.  If executors for Spark tasks are scheduled on-demand and can take a long time to start, it may be useful to increase this timeout on a system level.

Horovod on Databricks
------------------------------
To run Horovod in Spark on Databricks, create a Store instance with a DBFS path in one of the following patterns:

* ``/dbfs/...``
* ``dbfs:/...``
* ``file:///dbfs/...``

.. code-block:: python

    store = Store.create(dbfs_path)
    # or explicitly using DBFSLocalStore
    store = DBFSLocalStore(dbfs_path)

The `DBFSLocalStore` uses Databricks File System (DBFS) local file APIs
(`AWS <https://docs.databricks.com/data/databricks-file-system.html#local-file-apis>`__ |
`Azure <https://docs.microsoft.com/en-us/azure/databricks/data/databricks-file-system#--local-file-apis>`__)
as a store of intermediate data and training artifacts.

Databricks pre-configures GPU-aware scheduling on Databricks Runtime 7.0 ML GPU and above. See GPU scheduling instructions
(`AWS <https://docs.databricks.com/clusters/gpu.html#gpu-scheduling-1>`__ |
`Azure <https://docs.microsoft.com/en-us/azure/databricks/clusters/gpu#gpu-scheduling>`__)
for details.

With the Estimator API, horovod will launch ``# of tasks on each worker = # of GPUs on each worker``, and each task will
pin GPU to the assigned GPU from spark.

With the Run API, the function ``get_available_devices()`` from ``horovod.spark.task`` will return a list of assigned GPUs
for the spark task from which ``get_available_devices()`` is called.
See `keras_spark3_rossmann.py <../examples/spark/keras/keras_spark3_rossmann.py>`__ for an example of using
``get_available_devices()`` with the Run API.

.. inclusion-marker-end-do-not-remove
.. include:: ./autotune.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./docker.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove


Run Horovod
===========

This page includes examples for Open MPI that use ``horovodrun``. Check your
MPI documentation for arguments to the ``mpirun``
command on your system.

Typically one GPU will be allocated per process, so if a server has 4 GPUs,
you will run 4 processes. In ``horovodrun``,
the number of processes is specified with the ``-np`` flag.

To run on a machine with 4 GPUs:

.. code-block:: bash

    $ horovodrun -np 4 -H localhost:4 python train.py

To run on 4 machines with 4 GPUs each:

.. code-block:: bash

    $ horovodrun -np 16 -H server1:4,server2:4,server3:4,server4:4 python train.py

You can also specify host nodes in a host file. For example:

.. code-block:: bash

    $ cat myhostfile

    aa slots=2
    bb slots=2
    cc slots=2

This example lists the host names (aa, bb, and cc) and how many "slots" there
are for each.
Slots indicate how many processes can potentially execute on a node.
This format is the same as in
`mpirun command <https://www.open-mpi.org/doc/v4.0/man1/mpirun.1.php#toc6>`__.

To run on hosts specified in a hostfile:

.. code-block:: bash

    $ horovodrun -np 6 -hostfile myhostfile python train.py


Requirements
~~~~~~~~~~~~

Usage of ``horovodrun`` requires one of the following:

* Open MPI >= 2.X
* Spectrum MPI
* MPICH
* OpenRTE
* Gloo
* Intel(R) MPI

If you do not have MPI installed, you can run ``horovodrun`` using Gloo.  Gloo dependencies come with Horovod
automatically, and only require CMake to be available on your system at the time you install Horovod.

If you wish to use a different version of MPI, you may still be able to run Horovod using `mpirun <mpi.rst>`
directly.


Failures due to SSH issues
~~~~~~~~~~~~~~~~~~~~~~~~~~
The host where ``horovodrun`` is executed must be able to SSH to all other
hosts without any prompts.

If ``horovodrun`` fails with a permission error, verify that you can ssh to
every other server without entering a password or
answering questions like this:


``The authenticity of host '<hostname> (<ip address>)' can't be established.
RSA key fingerprint is xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx:xx.
Are you sure you want to continue connecting (yes/no)?``


To learn more about setting up passwordless authentication, see `this page <http://www.linuxproblem.org/art_9.html>`__.

To avoid ``The authenticity of host '<hostname> (<ip address>)' can't be
established`` prompts, add all the hosts to
the ``~/.ssh/known_hosts`` file using ``ssh-keyscan``:

.. code-block:: bash

    $ ssh-keyscan -t rsa,dsa server1 server2 > ~/.ssh/known_hosts


Advanced: Run Horovod with Open MPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In some advanced cases you might want fine-grained control over options passed to Open MPI.
To learn how to run Horovod training directly using Open MPI,
read `Run Horovod with Open MPI <mpi.rst>`_.

Run Horovod with Intel(R) MPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``horovodrun`` automatically converts some parameters to the format supported by Intel(R) MPI ``mpirun``. The set of allowed options includes ``-np``, ``-H`` and
ssh arguments (-p, -i). Intel(R) MPI ``mpirun`` does not support MCA parameters, but you can set some of the options via `environment variables <https://software.intel.com/content/www/us/en/develop/documentation/mpi-developer-reference-linux/environment-variable-reference.html>`__.
For additional information refer to `Intel(R) MPI official documentation <https://software.intel.com/content/www/us/en/develop/documentation/mpi-developer-reference-linux/top/command-reference/mpiexec-hydra/global-options.html>`__.

.. inclusion-marker-end-do-not-remove
.. include:: ./elastic.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove


Horovod in LSF
==============

This page includes examples for running Horovod in a LSF cluster.
``horovodrun`` will automatically detect the host names and GPUs of your LSF job.
If the LSF cluster supports ``jsrun``, ``horovodrun`` will use it as launcher
otherwise it will default to ``mpirun``.

Inside a LSF batch file or in an interactive session, you just need to use:

.. code-block:: bash

    horovodrun python train.py

Here, Horovod will start a process per GPU on all the hosts of the LSF job.

You can also limit the run to a subset of the job resources. For example, using only 6 GPUs:

.. code-block:: bash

    horovodrun -np 6 python train.py

You can still pass extra arguments to ``horovodrun``. For example, to trigger CUDA-Aware MPI:

.. code-block:: bash

    horovodrun --mpi-args="-gpu" python train.py

.. inclusion-marker-end-do-not-remove
.. include:: ./troubleshooting.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. include:: ./adasum_user_guide.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
Horovod documentation
=====================
Horovod improves the speed, scale, and resource utilization of deep learning training.

Get started
-----------
Choose your deep learning framework to learn how to get started with Horovod.

.. raw:: html

    <button class="accordion">TensorFlow</button>
    <div class="panel">
      <p>To use Horovod with TensorFlow on your laptop:
         <ol>
            <li><a href="https://www.open-mpi.org/faq/?category=building#easy-build">Install Open MPI 3.1.2 or 4.0.0</a>, or another MPI implementation. </li>
            <li>
               If you've installed TensorFlow from <a href="https://pypi.org/project/tensorflow">PyPI</a>, make sure that <code>g++-5</code> or above is installed.<br/>
               If you've installed TensorFlow from <a href="https://conda.io">Conda</a>, make sure that the <code>gxx_linux-64</code> Conda package is installed.
            </li>
            <li>Install the Horovod pip package: <code>pip install horovod</code></li>
            <li>Read <a href="https://horovod.readthedocs.io/en/latest/tensorflow.html">Horovod with TensorFlow</a> for best practices and examples. </li>
         </ol>
         Or, use <a href="https://horovod.readthedocs.io/en/latest/gpus_include.html">Horovod on GPUs</a>, in <a href="https://horovod.readthedocs.io/en/latest/spark_include.html">Spark</a>, <a href="https://horovod.readthedocs.io/en/latest/docker_include.html">Docker</a>, <a href="https://github.com/sylabs/examples/tree/master/machinelearning/horovod">Singularity</a>, or Kubernetes (<a href="https://github.com/kubeflow/examples/tree/master/demos/yelp_demo/ks_app/vendor/kubeflow/mpi-job">Kubeflow</a>, <a href="https://github.com/kubeflow/mpi-operator/">MPI Operator</a>, <a href="https://github.com/helm/charts/tree/master/stable/horovod">Helm Chart</a>, and <a href="https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/">FfDL</a>).
      </p>
    </div>

    <button class="accordion">Keras</button>
    <div class="panel">
      <p>To use Horovod with Keras on your laptop:
         <ol>
            <li><a href="https://www.open-mpi.org/faq/?category=building#easy-build">Install Open MPI 3.1.2 or 4.0.0</a>, or another MPI implementation. </li>
            <li>
               If you've installed TensorFlow from <a href="https://pypi.org/project/tensorflow">PyPI</a>, make sure that <code>g++-5</code> or above is installed.<br/>
               If you've installed TensorFlow from <a href="https://conda.io">Conda</a>, make sure that the <code>gxx_linux-64</code> Conda package is installed.
            </li>
            <li>Install the Horovod pip package: <code>pip install horovod</code></li>
            <li>Read <a href="https://horovod.readthedocs.io/en/latest/keras.html">Horovod with Keras</a> for best practices and examples. </li>
         </ol>
         Or, use <a href="https://horovod.readthedocs.io/en/latest/gpus_include.html">Horovod on GPUs</a>, in <a href="https://horovod.readthedocs.io/en/latest/spark_include.html">Spark</a>, <a href="https://horovod.readthedocs.io/en/latest/docker_include.html">Docker</a>, <a href="https://github.com/sylabs/examples/tree/master/machinelearning/horovod">Singularity</a>, or Kubernetes (<a href="https://github.com/kubeflow/examples/tree/master/demos/yelp_demo/ks_app/vendor/kubeflow/mpi-job">Kubeflow</a>, <a href="https://github.com/kubeflow/mpi-operator/">MPI Operator</a>, <a href="https://github.com/helm/charts/tree/master/stable/horovod">Helm Chart</a>, and <a href="https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/">FfDL</a>).
      </p>
    </div>

    <button class="accordion">PyTorch</button>
    <div class="panel">
      <p>To use Horovod with PyTorch on your laptop:
         <ol>
            <li><a href="https://www.open-mpi.org/faq/?category=building#easy-build">Install Open MPI 3.1.2 or 4.0.0</a>, or another MPI implementation. </li>
            <li>
               If you've installed PyTorch from <a href="https://pypi.org/project/torch">PyPI</a>, make sure that <code>g++-5</code> or above is installed.<br/>
               If you've installed PyTorch from <a href="https://conda.io">Conda</a>, make sure that the <code>gxx_linux-64</code> Conda package is installed.
            </li>
            <li>Install the Horovod pip package: <code>pip install horovod</code></li>
            <li>Read <a href="https://horovod.readthedocs.io/en/latest/pytorch.html">Horovod with PyTorch</a> for best practices and examples. </li>
         </ol>
         Or, use <a href="https://horovod.readthedocs.io/en/latest/gpus_include.html">Horovod on GPUs</a>, in <a href="https://horovod.readthedocs.io/en/latest/spark_include.html">Spark</a>, <a href="https://horovod.readthedocs.io/en/latest/docker_include.html">Docker</a>, <a href="https://github.com/sylabs/examples/tree/master/machinelearning/horovod">Singularity</a>, or Kubernetes (<a href="https://github.com/kubeflow/examples/tree/master/demos/yelp_demo/ks_app/vendor/kubeflow/mpi-job">Kubeflow</a>, <a href="https://github.com/kubeflow/mpi-operator/">MPI Operator</a>, <a href="https://github.com/helm/charts/tree/master/stable/horovod">Helm Chart</a>, and <a href="https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/">FfDL</a>).
      </p>
    </div>

    <button class="accordion">Apache MXNet</button>
    <div class="panel">
      <p>To use Horovod with Apache MXNet on your laptop:
         <ol>
            <li><a href="https://www.open-mpi.org/faq/?category=building#easy-build">Install Open MPI 3.1.2 or 4.0.0</a>, or another MPI implementation. </li>
            <li>Install the Horovod pip package: <code>pip install horovod</code></li>
            <li>Read <a href="https://horovod.readthedocs.io/en/latest/mxnet.html">Horovod with MXNet</a> for best practices and examples. </li>
         </ol>
         Or, use <a href="https://horovod.readthedocs.io/en/latest/gpus_include.html">Horovod on GPUs</a>, in <a href="https://horovod.readthedocs.io/en/latest/spark_include.html">Spark</a>, <a href="https://horovod.readthedocs.io/en/latest/docker_include.html">Docker</a>, <a href="https://github.com/sylabs/examples/tree/master/machinelearning/horovod">Singularity</a>, or Kubernetes (<a href="https://github.com/kubeflow/examples/tree/master/demos/yelp_demo/ks_app/vendor/kubeflow/mpi-job">Kubeflow</a>, <a href="https://github.com/kubeflow/mpi-operator/">MPI Operator</a>, <a href="https://github.com/helm/charts/tree/master/stable/horovod">Helm Chart</a>, and <a href="https://github.com/IBM/FfDL/tree/master/etc/examples/horovod/">FfDL</a>).
      </p>
    </div>

    <script>
        var acc = document.getElementsByClassName("accordion");
        var i;

        for (i = 0; i < acc.length; i++) {
          acc[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var panel = this.nextElementSibling;
            if (panel.style.maxHeight){
              panel.style.maxHeight = null;
            } else {
              panel.style.maxHeight = panel.scrollHeight + "px";
            }
          });
        }
     </script>

Guides
------

.. toctree::
   :maxdepth: 2

   summary_include

   concepts_include

   install_include

   api

   tensorflow

   xla

   keras

   pytorch

   mxnet

   running_include

   elastic_include

   benchmarks_include

   inference_include

   gpus_include

   mpi_include

   oneccl_include

   conda_include

   docker_include

   spark_include

   ray_include

   lsf_include

   tensor-fusion_include

   adasum_user_guide_include

   timeline_include

   hyperparameter_search_include

   autotune_include

   process_set_include

   troubleshooting_include

   contributors_include



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. inclusion-marker-start-do-not-remove

Horovod in Docker
=================

To streamline the installation process, we have published reference `Dockerfiles <https://github.com/horovod/horovod/blob/master/docker>`__ so
you can get started with Horovod in minutes. These containers include Horovod `examples <https://github.com/horovod/horovod/tree/master/examples>`__ in the ``/examples``
directory.

Pre-built Docker containers with Horovod are available on `DockerHub <https://hub.docker.com/r/horovod/horovod>`__ for GPU, CPU, and `Ray <https://ray.io>`__.


Running on a single machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~
After the container is built, run it using `nvidia-docker <https://github.com/NVIDIA/nvidia-docker>`__.

**Note**: You can replace ``horovod:latest`` with the `specific <https://hub.docker.com/r/horovod/horovod/tags>`__ pre-build
Docker container with Horovod instead of building it by yourself.

.. code-block:: bash

    $ nvidia-docker run -it horovod:latest
    root@c278c88dd552:/examples# horovodrun -np 4 -H localhost:4 python keras_mnist_advanced.py


If you don't run your container in privileged mode, you may see the following message:

.. code-block:: bash

    [a8c9914754d2:00040] Read -1, expected 131072, errno = 1


You can ignore this message.


Running on multiple machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here we describe a simple example involving a shared filesystem ``/mnt/share`` using a common port number ``12345`` for the SSH
daemon that will be run on all the containers. ``/mnt/share/ssh`` would contain a typical ``id_rsa`` and ``authorized_keys``
pair that allows `passwordless authentication <http://www.linuxproblem.org/art_9.html>`__.

**Note**: These are not hard requirements but they make the example more concise. A shared filesystem can be replaced by ``rsyncing``
SSH configuration and code across machines, and a common SSH port can be replaced by machine-specific ports
defined in ``/root/.ssh/ssh_config`` file.

Primary worker:

.. code-block:: bash

    host1$ nvidia-docker run -it --network=host -v /mnt/share/ssh:/root/.ssh horovod:latest
    root@c278c88dd552:/examples# horovodrun -np 16 -H host1:4,host2:4,host3:4,host4:4 -p 12345 python keras_mnist_advanced.py


Secondary workers:

.. code-block:: bash

    host2$ nvidia-docker run -it --network=host -v /mnt/share/ssh:/root/.ssh horovod:latest \
        bash -c "/usr/sbin/sshd -p 12345; sleep infinity"


.. code-block:: bash

    host3$ nvidia-docker run -it --network=host -v /mnt/share/ssh:/root/.ssh horovod:latest \
        bash -c "/usr/sbin/sshd -p 12345; sleep infinity"


.. code-block:: bash

    host4$ nvidia-docker run -it --network=host -v /mnt/share/ssh:/root/.ssh horovod:latest \
        bash -c "/usr/sbin/sshd -p 12345; sleep infinity"


Adding Mellanox RDMA support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you have Mellanox NICs, we recommend that you mount your Mellanox devices (``/dev/infiniband``) in the container
and enable the IPC_LOCK capability for memory registration:

.. code-block:: bash

   $ nvidia-docker run -it --network=host -v /mnt/share/ssh:/root/.ssh --cap-add=IPC_LOCK --device=/dev/infiniband horovod:latest 
   root@c278c88dd552:/examples# ...


You need to specify these additional configuration options on primary and secondary workers.


Running containers with different ports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To run in situations without a common SSH port (e.g., multiple containers on the same host):

1. Configure your `~/.ssh/config <https://linuxize.com/post/using-the-ssh-config-file>`__ file to assign custom host names and ports for each container

   .. code-block:: bash

        Host host1
          HostName 192.168.1.10
          Port 1234

        Host host2
          HostName 192.168.1.10
          Port 2345 

2. Use ``horovodrun`` directly as though each container were a separate host with its own IP
   
   .. code-block:: bash

        $ horovodrun -np 8 -H host1:4,host2:4 python keras_mnist_advanced.py

.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Analyze Performance
===================

Horovod has the ability to record the timeline of its activity, called Horovod Timeline.

.. image:: https://user-images.githubusercontent.com/16640218/29735271-9e148da0-89ac-11e7-9ae0-11d7a099ac89.png
   :alt: Horovod Timeline


To record a Horovod Timeline, set the ``--timeline-filename`` command line argument to the location of the timeline
file to be created.  This file is only recorded on rank 0, but it contains information about activity of all workers.

.. code-block:: bash

    $ horovodrun -np 4 --timeline-filename /path/to/timeline.json python train.py


You can then open the timeline file using the ``chrome://tracing`` facility of the `Chrome <https://www.google.com/chrome/browser/>`__ browser.

In the example above, you can see few tensors being reduced. There are two major phases for each tensor reduction:

1. **Negotiation** - a phase when all workers send to rank 0 signal that they're ready to reduce the given tensor.

   * Each worker reporting readiness is represented by a tick under the ``NEGOTIATE_ALLREDUCE`` bar, so you can see which workers were early and which were late.

   * Immediately after negotiation, rank 0 sends all other workers signal to start reducing the tensor.

2. **Processing** - a phase when the operation actually happens. It is further subdivided into multiple sub-phases:

   * ``WAIT_FOR_DATA`` indicates time taken to wait for GPU to finish computing input to the **allreduce**, *allgather*, or **broadcast** operations. This happens because TensorFlow tries to smartly interleave scheduling and GPU computation. This is only applicable to situations where the Horovod operation is placed on GPU.

   * ``WAIT_FOR_OTHER_TENSOR_DATA`` indicates time taken to wait for GPU to finish computing other inputs for other operations that are part of the same fusion batch.

   * ``QUEUE`` happens when reduction is done with NCCL, and the previous NCCL operation did not finish yet.

   * ``MEMCPY_IN_FUSION_BUFFER`` and ``MEMCPY_OUT_FUSION_BUFFER`` indicate time taken to copy data into and out of the fusion buffer.

   * ``NCCL_ALLREDUCE``, ``MPI_ALLREDUCE``, ``MPI_ALLGATHER``, or ``MPI_BCAST`` indicate time taken to do the actual operation on GPU (or CPU) and highlights whether the operation was performed using NCCL or pure MPI.

   * In case of ``HOROVOD_HIERARCHICAL_ALLREDUCE=1``, ``NCCL_ALLREDUCE`` will become a sequence or a subsequence of ``NCCL_REDUCESCATTER``, ``NCCL_REDUCE``, ``MEMCPY_IN_HOST_BUFFER``, ``MPI_ALLREDUCE``, ``MEMCPY_OUT_HOST_BUFFER``, ``NCCL_ALLGATHER``, ``NCCL_BCAST``.

Adding cycle markers
~~~~~~~~~~~~~~~~~~~~
Horovod performs work in cycles.  These cycles are used to aid `Tensor Fusion <https://github.com/horovod/horovod/blob/master/docs/tensor-fusion.rst>`__. Horovod has the ability to record the moment when each cycle starts for debugging of Tensor Fusion.

.. image:: https://user-images.githubusercontent.com/16640218/51659458-64806100-1f5f-11e9-9a27-ba934ceec75f.png
   :alt: Cycle Markers


Since this information makes timeline view very crowded, it is not enabled by default. To add cycle markers to the timeline, set the ``--timeline-mark-cycles`` flag:

.. code-block:: bash

    $ horovodrun -np 4 --timeline-filename /path/to/timeline.json --timeline-mark-cycles python train.py

.. inclusion-marker-end-do-not-remove
.. include:: ./benchmarks.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Horovod Installation Guide
==========================

Requirements
------------

- GNU Linux or macOS
- Python >= 3.6
- `g++-5` or above, or another compiler supporting C++14
- CMake 3.13 or newer
- TensorFlow, PyTorch, or MXNet
- (Optional) MPI

For best performance on GPU:

- `NCCL 2 <https://developer.nvidia.com/nccl>`__

If Horovod in unable to find the CMake binary, you may need to set ``HOROVOD_CMAKE`` in your environment before
installing.

Horovod does not support Windows.

Frameworks
----------

You can build Horovod for TensorFlow, PyTorch, and MXNet. By default, Horovod will attempt to build
support for all of them. At least one must be enabled for Horovod to install successfully.

To ensure that framework dependencies are properly installed before attempting to install Horovod, append
extra arguments that identify the required frameworks:

.. code-block:: bash

    $ pip install horovod[tensorflow,keras,pytorch,mxnet,spark]

In addition to specifying framework requirements individually, you can require all frameworks collectively:

.. code-block:: bash

    $ pip install horovod[all-frameworks]

This is useful when building Horovod as part of a larger collection of dependencies at once, relying on the pip
compiler to determine the correct install order.

TensorFlow
~~~~~~~~~~

To ensure that Horovod is built with TensorFlow support enabled:

.. code-block:: bash

    $ HOROVOD_WITH_TENSORFLOW=1 pip install horovod[tensorflow]

To skip TensorFlow, set ``HOROVOD_WITHOUT_TENSORFLOW=1`` in your environment.

PyTorch
~~~~~~~

To ensure that Horovod is built with PyTorch support enabled:

.. code-block:: bash

    $ HOROVOD_WITH_PYTORCH=1 pip install horovod[pytorch]

To skip PyTorch, set ``HOROVOD_WITHOUT_PYTORCH=1`` in your environment.

MXNet
~~~~~

To ensure that Horovod is built with MXNet CPU support enabled:

.. code-block:: bash

    $ HOROVOD_WITH_MXNET=1 pip install horovod[mxnet]

Some MXNet versions do not work with Horovod:

- MXNet 1.4.0 and earlier have `GCC incompatibility issues <https://github.com/horovod/horovod/issues/884>`__. Use MXNet 1.4.1 or later with Horovod 0.16.2 or later to avoid these incompatibilities.
- MXNet 1.5.1, 1.6.0, 1.7.0, and 1.7.0.post1 are missing MKLDNN headers, so they do not work with Horovod. Use 1.5.1.post0, 1.6.0.post0, and 1.7.0.post0, respectively.
- MXNet 1.6.0.post0 and 1.7.0.post0 are only available as mxnet-cu101 and mxnet-cu102.

To skip MXNet, set ``HOROVOD_WITHOUT_MXNET=1`` in your environment.

Keras
~~~~~

Standalone Keras support is currently only available for the TensorFlow backend.

To ensure that Horovod is built with Keras support available:

.. code-block:: bash

    $ HOROVOD_WITH_TENSORFLOW=1 pip install horovod[tensorflow,keras]

There are no plugins built for Keras, but the TensorFlow plugin must be enabled in order to use Horovod with Keras.

Spark
~~~~~

Horovod can be used with Spark in combination with any of the frameworks above.

To ensure Horovod has all the necessary requirements in order to run on top of Spark:

.. code-block:: bash

    $ pip install horovod[spark]

Controllers
-----------

The controller is used for coordinating work between Horovod processes (determining which tensors to process). We
provide controller implementations for both MPI and Gloo. By default, Horovod will attempt to build support for both
of them. At least one must be enabled for Horovod to install successfully.

MPI
~~~

MPI is the original controller for Horovod.  It uses ``mpirun`` to launch worker processes (``horovodrun`` will use
``mpirun`` under the hood when using MPI).

To use Horovod with MPI, install `Open MPI <https://www.open-mpi.org/>`_ or another MPI implementation.
Learn how to install Open MPI `on this page <https://www.open-mpi.org/faq/?category=building#easy-build>`_.

**Note**: Open MPI 3.1.3 has an issue that may cause hangs. The recommended fix is to downgrade to Open MPI 3.1.2 or
upgrade to Open MPI 4.0.0.

* To force Horovod to install with MPI support, set ``HOROVOD_WITH_MPI=1`` in your environment.
* To force Horovod to skip building MPI support, set ``HOROVOD_WITHOUT_MPI=1``.

If both MPI and Gloo are enabled in your installation, then MPI will be the default controller.

Gloo
~~~~

Gloo is a more recent controller for Horovod that does not require additional dependencies besides CMake to install.

When used as a controller in combination with NCCL, Gloo performs almost identically to MPI on standard benchmarks.

* To force Horovod to install with Gloo support, set ``HOROVOD_WITH_GLOO=1`` in your environment.
* To force Horovod to skip building Gloo support, set ``HOROVOD_WITHOUT_GLOO=1``.

Gloo mode uses ``horovodrun`` to launch worker processes.

Gloo is required to use the elastic / fault tolerant API for Horovod.

**Note**: macOS users must install `libuv <https://github.com/libuv/libuv>`_ in order to use Gloo:

.. code-block:: bash

    $ brew install libuv

Tensor Operations
-----------------

For running on GPUs with optimal performance, we recommend installing Horovod with NCCL support following the
`Horovod on GPU <gpus.rst>`_ guide.

For tensor data on CPU, you can use MPI, Gloo, and Intel's oneCCL. By default, the framework used by your controller
will be used for CPU operations. You can override this by setting ``HOROVOD_CPU_OPERATIONS`` in your environment.

NCCL
~~~~

NCCL is supported for Allreduce, Allgather, Broadcast, and Alltoall operations.  You can enable these by setting
``HOROVOD_GPU_OPERATIONS=NCCL`` during installation.

NCCL operations are supported on both Nvidia (CUDA) and AMD (ROCm) GPUs. You can set ``HOROVOD_GPU`` in your
environment to specify building with CUDA or ROCm. CUDA will be assumed if not specified.

Note that Alltoall requires NCCL version >= 2.7.0.

MPI
~~~

When using an MPI controller, MPI will be used when NCCL is unavailable, or if tensors are placed in host memory prior
to the allreduce request. In cases where NCCL is unavailable, MPI has been shown to outperform Gloo for CPU tensor
operations.

MPI can also be used for GPU operations, but this is not recommended in most cases. See `Horovod on GPU <gpus.rst>`_ for
more details.

Gloo
~~~~

When using a Gloo controller, Gloo will be used in place of MPI for CPU operations by default.

oneCCL
~~~~~~

oneCCL is an Intel library for accelerated collective operations on CPU. See
`Horovod with Intel(R) oneCCL <oneccl.rst>`_ for more details.

Set ``HOROVOD_CPU_OPERATIONS=CCL`` to use oneCCL.


Check Build
-----------

After successfully installing Horovod, run:

.. code-block:: bash

    $ horovodrun --check-build

Every feature that was successfully enabled will be marked with an 'X'. If you intended to install Horovod with a
feature that is not listed as enabled, you can reinstall Horovod, setting the appropriate environment variables to
diagnose failures:

.. code-block:: bash

    $ pip uninstall horovod
    $ HOROVOD_WITH_...=1 pip install --no-cache-dir horovod

Installing Horovod with Conda (+pip)
------------------------------------

To use Conda to install PyTorch, TensorFlow, MXNet, Horovod, as well as GPU dependencies such as
NVIDIA CUDA Toolkit, cuDNN, NCCL, etc., see `Build a Conda Environment with GPU Support for Horovod <conda.rst>`_.

Environment Variables
---------------------

Optional environment variables that can be set to configure the installation process for Horovod.

Possible values are given in curly brackets: {}.

* ``HOROVOD_DEBUG`` - {1}. Install a debug build of Horovod with checked assertions, disabled compiler optimizations etc.
* ``HOROVOD_BUILD_ARCH_FLAGS`` - additional C++ compilation flags to pass in for your build architecture.
* ``HOROVOD_CUDA_HOME`` - path where CUDA include and lib directories can be found.
* ``HOROVOD_BUILD_CUDA_CC_LIST`` - List of compute capabilities to build Horovod CUDA kernels for (example: ``HOROVOD_BUILD_CUDA_CC_LIST=60,70,75``)
* ``HOROVOD_ROCM_HOME`` - path where ROCm include and lib directories can be found.
* ``HOROVOD_NCCL_HOME`` - path where NCCL include and lib directories can be found.
* ``HOROVOD_NCCL_INCLUDE`` - path to NCCL include directory.
* ``HOROVOD_NCCL_LIB`` - path to NCCL lib directory.
* ``HOROVOD_NCCL_LINK`` - {SHARED, STATIC}. Mode to link NCCL library. Defaults to STATIC for CUDA, SHARED for ROCm.
* ``HOROVOD_WITH_GLOO`` - {1}. Require that Horovod is built with Gloo support enabled.
* ``HOROVOD_WITHOUT_GLOO`` - {1}. Skip building with Gloo support.
* ``HOROVOD_WITH_MPI`` - {1}. Require that Horovod is built with MPI support enabled.
* ``HOROVOD_WITHOUT_MPI`` - {1}. Skip building with MPI support.
* ``HOROVOD_GPU`` - {CUDA, ROCM}. Framework to use for GPU operations.
* ``HOROVOD_GPU_OPERATIONS`` - {NCCL, MPI}. Framework to use for GPU tensor allreduce, allgather, and broadcast.
* ``HOROVOD_GPU_ALLREDUCE`` - {NCCL, MPI}. Framework to use for GPU tensor allreduce.
* ``HOROVOD_GPU_ALLGATHER`` - {NCCL, MPI}. Framework to use for GPU tensor allgather.
* ``HOROVOD_GPU_BROADCAST`` - {NCCL, MPI}. Framework to use for GPU tensor broadcast.
* ``HOROVOD_ALLOW_MIXED_GPU_IMPL`` - {1}. Allow Horovod to install with NCCL allreduce and MPI GPU allgather / broadcast.  Not recommended due to a possible deadlock.
* ``HOROVOD_CPU_OPERATIONS`` - {MPI, GLOO, CCL}. Framework to use for CPU tensor allreduce, allgather, and broadcast.
* ``HOROVOD_CMAKE`` - path to the CMake binary used to build Gloo (not required when using MPI).
* ``HOROVOD_WITH_TENSORFLOW`` - {1}. Require Horovod to install with TensorFlow support enabled.
* ``HOROVOD_WITHOUT_TENSORFLOW`` - {1}. Skip installing TensorFlow support.
* ``HOROVOD_WITH_PYTORCH`` - {1}. Require Horovod to install with PyTorch support enabled.
* ``HOROVOD_WITHOUT_PYTORCH`` - {1}. Skip installing PyTorch support.
* ``HOROVOD_WITH_MXNET`` - {1}. Require Horovod to install with MXNet support enabled.
* ``HOROVOD_WITHOUT_MXNET`` - {1}. Skip installing MXNet support.

.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Process Sets: Concurrently Running Collective Operations
========================================================

Most Horovod operations in TensorFlow, PyTorch, or MXNet feature a ``process_set`` argument: By setting up different
process sets you may have multiple subsets of the world of Horovod processes run distinct collective operations in
parallel. Besides Horovod's fundamental operations like ``hvd.allgather``, ``hvd.allreduce``, ``hvd.alltoall``,
``hvd.broadcast``, or ``hvd.grouped_allreduce``, also many high-level utility objects such as
``hvd.DistributedOptimizer`` come with support for process sets.

As an example consider building a Horovod model to be trained by four worker processes with two concurrent allreduce
operations on the "even" or "odd" subset.  In the following we will see three ways to configure Horovod to use an even
and an odd process set, offering you as much flexibility as you need. The code examples are presented for TensorFlow,
but the interface for the other supported frameworks is equivalent.

1) Static process sets
----------------------

.. code-block:: python

    # on all ranks
    even_set = hvd.ProcessSet([0,2])
    odd_set = hvd.ProcessSet([1,3])
    hvd.init(process_sets=[even_set, odd_set])

    for p in [hvd.global_process_set, even_set, odd_set]:
      print(p)
    # ProcessSet(process_set_id=0, ranks=[0, 1, 2, 3], mpi_comm=None)
    # ProcessSet(process_set_id=1, ranks=[0, 2], mpi_comm=None)
    # ProcessSet(process_set_id=2, ranks=[1, 3], mpi_comm=None)

    # on ranks 0 and 2
    result = hvd.allreduce(tensor_for_even_ranks, process_set=even_set)

    # on ranks 1 and 3
    result = hvd.allreduce(tensor_for_odd_ranks, process_set=odd_set)

Having initialized Horovod like this, the configuration of process sets cannot be changed without restarting the
program.  If you only use the default global process set (``hvd.global_process_set``), there is no impact on
performance.

2) Static process sets from MPI communicators
---------------------------------------------

.. code-block:: python

    # on all ranks
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    subcomm = MPI.COMM_WORLD.Split(color=MPI.COMM_WORLD.rank % 2,
                                   key=MPI.COMM_WORLD.rank)

    split_process_set = hvd.ProcessSet(subcomm)

    hvd.init(comm, process_sets=[split_process_set])

    for p in [hvd.global_process_set, split_process_set]:
        print(p)
    # ProcessSet(process_set_id=0, ranks=[0, 1, 2, 3], mpi_comm=<mpi4py.MPI.Intracomm object at 0x7fb817323dd0>)
    # ProcessSet(process_set_id=1, ranks=[0, 2], mpi_comm=<mpi4py.MPI.Intracomm object at 0x7fb87e2ddfb0>)
    ## (split_process_set differs by rank)

    # on ranks 0 and 2
    result = hvd.allreduce(tensor_for_even_ranks, process_set=split_process_set)

    # on ranks 1 and 3
    result = hvd.allreduce(tensor_for_odd_ranks, process_set=split_process_set)

If you are already using multiple MPI communicators in your distributed program, you can plug them right in.

3) Dynamic process sets
-----------------------

.. code-block:: python

    # on all ranks
    hvd.init(process_sets="dynamic")  # alternatively set HOROVOD_DYNAMIC_PROCESS_SETS=1
    even_set = hvd.add_process_set([0,2])
    odd_set = hvd.add_process_set([1,3])

    for p in [hvd.global_process_set, even_set, odd_set]:
      print(p)
    # ProcessSet(process_set_id=0, ranks=[0, 1, 2, 3], mpi_comm=None)
    # ProcessSet(process_set_id=1, ranks=[0, 2], mpi_comm=None)
    # ProcessSet(process_set_id=2, ranks=[1, 3], mpi_comm=None)

    # on ranks 0 and 2
    result = hvd.allreduce(tensor_for_even_ranks, process_set=even_set)

    # on ranks 1 and 3
    result = hvd.allreduce(tensor_for_odd_ranks, process_set=odd_set)

The most flexible setup is achieved with "dynamic" process sets.  Process sets can be registered and deregistered
dynamically at any time after initializing Horovod via ``hvd.add_process_set()`` and ``hvd.remove_process_set()``.
Calls to these functions must be made identically and in the same order by all processes.

Note that dynamic process sets come with some slight extra synchronization overhead.

.. inclusion-marker-end-do-not-remove
.. inclusion-marker-start-do-not-remove

Build a Conda Environment with GPU Support for Horovod
======================================================

In this section we describe how to build Conda environments for deep learning projects using 
Horovod to enable distributed training across multiple GPUs (either on the same node or 
spread across multuple nodes).

Installing the NVIDIA CUDA Toolkit
----------------------------------

Install `NVIDIA CUDA Toolkit 10.1`_ (`documentation`_) which is the most recent version of NVIDIA 
CUDA Toolkit supported by all three deep learning frameworks that are currently supported by 
Horovod.

Why not just use the ``cudatoolkit`` package?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically when installing PyTorch, TensorFlow, or Apache MXNet with GPU support using Conda, you 
add the appropriate version of the ``cudatoolkit`` package to your ``environment.yml`` file. 
Unfortunately, for the moment at least, the cudatoolkit packages available via Conda do not 
include the `NVIDIA CUDA Compiler (NVCC)`_, which is required in order to build Horovod extensions 
for PyTorch, TensorFlow, or MXNet.

What about the ``cudatoolkit-dev`` package?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While there are ``cudatoolkit-dev`` packages available from ``conda-forge`` that do include NVCC, 
we have had difficulty getting these packages to consistently install properly. Some of the 
available builds require manual intervention to accept license agreements, making these builds 
unsuitable for installing on remote systems (which is critical functionality). Other builds seems 
to work on Ubuntu but not on other flavors of Linux.

Despite this, we would encourage you to try adding ``cudatoolkit-dev`` to your ``environment.yml`` 
file and see what happens! The package is well maintained so perhaps it will become more stable in 
the future.

Use the ``nvcc_linux-64`` meta-package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most robust approach to obtain NVCC and still use Conda to manage all the other dependencies 
is to install the NVIDIA CUDA Toolkit on your system and then install a meta-package 
`nvcc_linux-64`_ from conda-forge, which configures your Conda environment to use the NVCC 
installed on the system together with the other CUDA Toolkit components installed inside the Conda 
environment.

The ``environment.yml`` file
----------------------------

We prefer to specify as many dependencies as possible in the Conda ``environment.yml`` file 
and only specify dependencies in ``requirements.txt`` for install via ``pip`` that are not 
available via Conda channels. Check the Horovod `installation guide`_ for details of required 
dependencies.

Channel Priority
^^^^^^^^^^^^^^^^

Use the recommended channel priorities. Note that ``conda-forge`` has priority over 
``defaults`` and ``pytorch`` has priority over ``conda-forge``. ::

    name: null

    channels:
    - pytorch
    - conda-forge
    - defaults

Dependencies
^^^^^^^^^^^^

There are a few things worth noting about the dependencies.

1. Even though you have installed the NVIDIA CUDA Toolkit manually, you should still use Conda to 
   manage the other required CUDA components such as ``cudnn`` and ``nccl`` (and the optional 
   ``cupti``).
2. Use two meta-packages, ``cxx-compiler`` and ``nvcc_linux-64``, to make sure that suitable C, 
   and C++ compilers are installed and that the resulting Conda environment is aware of the 
   manually installed CUDA Toolkit.
3. Horovod requires some controller library to coordinate work between the various Horovod 
   processes. Typically this will be some MPI implementation such as `OpenMPI`_. However, rather 
   than specifying the ``openmpi`` package directly, you should instead opt for `mpi4py`_ Conda 
   package which provides a CUDA-aware build of OpenMPI.
4. Horovod also support the `Gloo`_ collective communications library that can be used in place of 
   MPI. Include ``cmake`` to insure that the Horovod extensions for Gloo are built.

Below are the core required dependencies. The complete ``environment.yml`` file is available 
on `GitHub`_. ::

    dependencies:
    - bokeh=1.4
    - cmake=3.16 # insures that Gloo library extensions will be built
    - cudnn=7.6
    - cupti=10.1
    - cxx-compiler=1.0 # insures C and C++ compilers are available
    - jupyterlab=1.2
    - mpi4py=3.0 # installs cuda-aware openmpi
    - nccl=2.5
    - nodejs=13
    - nvcc_linux-64=10.1 # configures environment to be "cuda-aware"
    - pip=20.0
    - pip:
        - mxnet-cu101mkl==1.6.* # MXNET is installed prior to horovod
        - -r file:requirements.txt
    - python=3.7
    - pytorch=1.5
    - tensorboard=2.1
    - tensorflow-gpu=2.1
    - torchvision=0.6

The ``requirements.txt`` file
-----------------------------

The ``requirements.txt`` file is where all of the ``pip`` dependencies, including Horovod itself, 
are listed for installation. In addition to Horovod we typically will also use ``pip`` to install 
JupyterLab extensions to enable GPU and CPU resource monitoring via `jupyterlab-nvdashboard`_ and 
Tensorboard support via `jupyter-tensorboard`_. ::

    horovod==0.19.*
    jupyterlab-nvdashboard==0.2.*
    jupyter-tensorboard==0.2.*

    # make sure horovod is re-compiled if environment is re-built
    --no-binary=horovod

Note the use of the ``--no-binary`` option at the end of the file. Including this option ensures 
that Horovod will be re-built whenever the Conda environment is re-built.

Building the Conda environment
------------------------------

After adding any necessary dependencies that should be downloaded via Conda to the 
``environment.yml`` file and any dependencies that should be downloaded via ``pip`` to the 
``requirements.txt`` file, create the Conda environment in a sub-directory ``env`` of your 
project directory by running the following commands.

.. code-block:: bash

    $ export ENV_PREFIX=$PWD/env
    $ export HOROVOD_CUDA_HOME=$CUDA_HOME
    $ export HOROVOD_NCCL_HOME=$ENV_PREFIX
    $ export HOROVOD_GPU_OPERATIONS=NCCL
    $ conda env create --prefix $ENV_PREFIX --file environment.yml --force

By default Horovod will try and build extensions for all detected frameworks. See the 
documentation on `environment variables`_ for the details on additional environment variables that 
can be set prior to building Horovod.

Once the new environment has been created you can activate the environment with the following 
command.

.. code-block:: bash

    $ conda activate $ENV_PREFIX

The ``postBuild`` file
^^^^^^^^^^^^^^^^^^^^^^

If you wish to use any JupyterLab extensions included in the ``environment.yml`` and 
``requirements.txt`` files, then you may need to rebuild the JupyterLab application.

For simplicity, we typically include the instructions for re-building JupyterLab in a 
``postBuild`` script. Here is what this script looks like in the example Horovod environments.

.. code-block:: bash

    jupyter labextension install --no-build jupyterlab-nvdashboard 
    jupyter labextension install --no-build jupyterlab_tensorboard
    jupyter lab build

Use the following commands to source the ``postBuild`` script.

.. code-block:: bash

    $ conda activate $ENV_PREFIX # optional if environment already active
    $ . postBuild

Listing the contents of the Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To see the full list of packages installed into the environment, run the following command.

.. code-block:: bash

    $ conda activate $ENV_PREFIX # optional if environment already active
    $ conda list

Verifying the Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After building the Conda environment, check that Horovod has been built with support for the deep 
learning frameworks TensorFlow, PyTorch, Apache MXNet, and the contollers MPI and Gloo with the 
following command.

.. code-block:: bash

    $ conda activate $ENV_PREFIX # optional if environment already active
    $ horovodrun --check-build

You should see output similar to the following.::

    Horovod v0.19.4:
    Available Frameworks:
        [X] TensorFlow
        [X] PyTorch
        [X] MXNet
    Available Controllers:
        [X] MPI
        [X] Gloo
    Available Tensor Operations:
        [X] NCCL
        [ ] DDL
        [ ] CCL
        [X] MPI
        [X] Gloo

Wrapping it all up in a Bash script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We typically wrap these commands into a shell script ``create-conda-env.sh``. Running the shell 
script will set the Horovod build variables, create the Conda environment, activate the Conda 
environment, and build JupyterLab with any additional extensions.

.. code-block:: bash

    #!/bin/bash --login

    set -e
    
    export ENV_PREFIX=$PWD/env
    export HOROVOD_CUDA_HOME=$CUDA_HOME
    export HOROVOD_NCCL_HOME=$ENV_PREFIX
    export HOROVOD_GPU_OPERATIONS=NCCL
    conda env create --prefix $ENV_PREFIX --file environment.yml --force
    conda activate $ENV_PREFIX
    . postBuild

We recommend that you put scripts inside a ``bin`` directory in your project root directory. Run 
the script from the project root directory as follows.

.. code-block:: bash

    ./bin/create-conda-env.sh # assumes that $CUDA_HOME is set properly

Updating the Conda environment
------------------------------

If you add (remove) dependencies to (from) the ``environment.yml`` file or the 
``requirements.txt`` file after the environment has already been created, then you can 
re-create the environment with the following command.

.. code-block:: bash

    $ conda env create --prefix $ENV_PREFIX --file environment.yml --force

However, whenever we add (remove) dependencies we prefer to re-run the Bash script which will re-build 
both the Conda environment and JupyterLab.

.. code-block:: bash

    $ ./bin/create-conda-env.sh

.. _NVIDIA CUDA Toolkit 10.1: https://developer.nvidia.com/cuda-10.1-download-archive-update2
.. _documentation: https://docs.nvidia.com/cuda/archive/10.1/
.. _NVIDIA CUDA Compiler (NVCC): https://docs.nvidia.com/cuda/archive/10.1/cuda-compiler-driver-nvcc/index.html
.. _nvcc_linux-64: https://github.com/conda-forge/nvcc-feedstock
.. _installation guide: https://horovod.readthedocs.io/en/latest/install_include.html
.. _OpenMPI: https://www.open-mpi.org/
.. _mpi4py: https://mpi4py.readthedocs.io/en/stable/
.. _Gloo: https://github.com/facebookincubator/gloo
.. _GitHub: https://github.com/kaust-vislab/horovod-gpu-data-science-project
.. _jupyterlab-nvdashboard: https://github.com/rapidsai/jupyterlab-nvdashboard
.. _jupyter-tensorboard: https://github.com/lspvic/jupyter_tensorboard
.. _environment variables: https://horovod.readthedocs.io/en/latest/install_include.html#environment-variables

.. inclusion-marker-end-do-not-removeOverview
========

.. include:: ./summary.rst
   :start-after: inclusion-marker-start-do-not-remove
   :end-before: inclusion-marker-end-do-not-remove
API
===

horovod.tensorflow
------------------
.. automodule:: horovod.tensorflow
.. automodule:: horovod.tensorflow.elastic

horovod.tensorflow.keras
------------------------
.. automodule:: horovod.tensorflow.keras
.. automodule:: horovod.tensorflow.keras.callbacks
    :special-members: __init__
.. automodule:: horovod.tensorflow.keras.elastic

horovod.keras
-------------
.. automodule:: horovod.keras
.. automodule:: horovod.keras.callbacks
    :special-members: __init__
.. automodule:: horovod.keras.elastic

horovod.torch
-------------
.. automodule:: horovod.torch
.. automodule:: horovod.torch.elastic

horovod.mxnet
-------------
.. automodule:: horovod.mxnet

horovod.spark
-------------
.. automodule:: horovod.spark

horovod.spark.keras
-------------------
.. autoclass:: horovod.spark.keras.KerasEstimator
    :show-inheritance:
.. autoclass:: horovod.spark.keras.KerasModel
    :show-inheritance:

horovod.spark.torch
-------------------
.. autoclass:: horovod.spark.torch.TorchEstimator
    :show-inheritance:
.. autoclass:: horovod.spark.torch.TorchModel
    :show-inheritance:

horovod.spark.common
--------------------
.. autoclass:: horovod.spark.common.estimator.HorovodEstimator
.. autoclass:: horovod.spark.common.estimator.HorovodModel
.. automodule:: horovod.spark.common.backend
    :show-inheritance:
.. automodule:: horovod.spark.common.store
    :show-inheritance:

.. _horovod_ray_api:

horovod.ray
-----------
.. automodule:: horovod.ray

horovod.run
-------------
.. automodule:: horovod.run
