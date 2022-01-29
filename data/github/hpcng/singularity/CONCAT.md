# Singularity Support

Need help? We have several ways to reach us, depending on your preferences and needs.

## Documentation

If you haven't already, read our [documentation](https://singularity.hpcng.org/docs/).
These docs have common use cases, and might be helpful to browse before
submitting an issue.
You can contribute to our documentation by creating PRs against the
[singularity-admindocs](https://github.com/hpcng/singularity-admindocs) and
[singularity-userdocs](https://github.com/hpcng/singularity-userdocs) repositories.

## Github

For issues with code (and especially if you need to share debug output)
we recommend Github issues boards.

- [Singularity Issues](https://github.com/hpcng/singularity/issues): is
    recommended for most issues with the Singularity software.
- [User Documentation](https://github.com/hpcng/singularity-userdocs/issues)
    questions, feedback, and suggestions should go here. Feel free to
    create an issue on a board and additionally request updated content
    here.
- [Admin Documentation](https://github.com/hpcng/singularity-admindocs/issues)
    questions, feedback, and suggestions should go here. Feel free to create
    an issue on a board and additionally request updated content here.

Note that usage questions, or problems related to running a specific piece of
software in a container, are best asked on the Google Group, or Slack channel.
Questions in these venues will be seen by a greater number of users, who may
already know the answer!

### How do I ask for help?

After you identify a bug, you should search the respective issue board
for similar problems reported by other users. Another user may be facing
the same issue, and you can add a +1 (in message or icon) to indicate to
the maintainers that the issue is pressing for you as well. The squeaky
wheel gets the grease!

### How is time allocated to addressing issues?

While we wish we could address every issue, there are only so many hours
in the day. We rank issues based on the following questions:

1. How many users are affected?
1. Is there a proposed work-around?
1. In how many instances does the proposed work-around fail?

With these simple questions, we can ensure that work is directed and has
the maximum impact! However, if your issue doesn't seem to be getting
attention you can still move it along using some of the strategies
discussed below.

### What if my issue goes stale?

Issues can go stale for a number of reasons. In the bullets below, we
will review some of these reasons, along with strategies for managing
them:

1. *The issue needs a gentle reminder*. Try targeting a few people with
    a "`ping @username any thoughts about this?`" in the case that it
    was forgotten.
1. *Was your issue properly explained*? You are much more likely to get
    help when you give clear instructions for reproducing the issue, and
    show effort on your part to think about what the problem might be.
    If possible, try to come up with a way to reproduce the issue that
    does not involve a special environment or exotic hardware.
1. *Is there broad need*? It could be that your issue isn't having a big
    enough impact for other users to warrant the time for the small
    development team. In this case, you might try implementing a
    suggested fix, and then asking for help with the details.
1. *Is your issue scattered?* When many issues pile up on boards, it
    sometimes is the case that issues are duplicated. It's important to
    find these duplicates and merge them into one, because in finding
    the duplicate you find another user to talk to about the issue.
1. *Does your issue need to have scope?* The idea of scoping an issue
    means framing it with respect to other components of the software.
    For example, if you have a feature request to see metadata about an
    object, you might frame that in the context of container
    introspection, and suggest an addition to the software that fits
    with the "inspect" command. A very powerful thing to do would be to
    open up an issue that (not only discusses your specific addition)
    but also opens up discussion to the general community for "How we
    can do introspection" better. Then create a set of issues and add
    them to a [Github
    milestone](https://help.github.com/articles/about-milestones/).
    This kind of contribution is much more powerful than simply asking
    for something.

## Google Group

You can reach the community quickly by way of joining our [Google
Group](https://groups.google.com/a/lbl.gov/forum/#!forum/singularity).

## Slack

For real time support from the community, you can join our community on
slack at [https://hpcng.slack.com/](https://hpcng.slack.com/). Ping the
Google Group or one of the admins here to request to be added.

Is there something missing here you'd like to see? Please [let us
know](https://github.com/hpcng/singularity/issues).
# Copyright

Copyright (c) 2017-2018, Sylabs, Inc. All rights reserved.

Copyright (c) 2015-2017, Gregory M. Kurtzer. All rights reserved.

Copyright (c) 2016-2017, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.
# Singularity Changelog

## Changes since last release

- `--writable-tmpfs` can be used with `singularity build` to run
  the `%test` section of the build with a ephemeral tmpfs overlay,
  permitting tests that write to the container filesystem.
- `--compat` flag for actions is a new short-hand to enable a number of
  options that increase OCI/Docker compatibility. Infers `--containall,
  --no-init, --no-umask, --writable-tmpfs`. Does not use user, uts, or
  network namespaces as these may not be supported on many installations.
- `--no-https` now applies to connections made to library services specified
  in `--library://<hostname>/...` URIs.
- The experimental `--nvccli` flag will use `nvidia-container-cli` to setup the
  container for Nvidia GPU operation. Singularity will not bind GPU libraries
  itself. Environment variables that are used with Nvidia's `docker-nvidia`
  runtime to configure GPU visibility / driver capabilities & requirements are
  parsed by the `--nvccli` flag from the environment of the calling user. By
  default, the `compute` and `utility` GPU capabilities are configured. The `use
  nvidia-container-cli` option in `singularity.conf` can be set to `yes` to
  always use `nvidia-container-cli` when supported.
  `--nvccli` is not supported in the setuid workflow,
  and it requires being used in combination with `--writable` in user
  namespace mode.
  Please see documentation for more details.
- A new `--mount` flag and `SINGULARITY_MOUNT` environment variable can be used
  to specify bind mounts in
  `type=bind,source=<src>,destination=<dst>[,options...]` format. This improves
  CLI compatibility with other runtimes, and allows binding paths containing
  `:` and `,` characters (using CSV style escaping).
- Perform concurrent multi-part downloads for `library://` URIs. Uses 3
  concurrent downloads by default, and is configurable in `singularity.conf` or
  via environment variables.

### Changed defaults / behaviours

- Building Singularity from source requires go >=1.16. We now aim to support
  the two most recent stable versions of Go. This corresponds to the Go
  [Release Maintenance Policy](https://github.com/golang/go/wiki/Go-Release-Cycle#release-maintenance)
  and [Security Policy](https://golang.org/security), ensuring critical bug
  fixes and security patches are available for all supported language versions.
  However, rpm and debian packaging apply patches to support older native go
  installations.
- LABELs from Docker/OCI images are now inherited. This fixes a longstanding
  regression from Singularity 2.x. Note that you will now need to use
  `--force` in a build to override a label that already exists in the source
  Docker/OCI container.
- Instances are no longer created with an IPC namespace by default. An IPC
  namespace can be specified with the `-i|--ipc` flag.
- `--bind`, `--nv` and `--rocm` options for `build` command can't be set through
  environment variables `SINGULARITY_BIND`, `SINGULARITY_BINDPATH`, `SINGULARITY_NV`,
  `SINGULARITY_ROCM` anymore due to side effects reported by users in this
  [issue](https://github.com/hpcng/singularity/pull/6211), they must be explicitely
  requested via command line.
- `--nohttps` flag has been deprecated in favour of `--no-https`. The old flag
  is still accepted, but will display a deprecation warning.
- Removed `--nonet` flag, which was intended to disable networking for in-VM
  execution, but has no effect.
- Paths for `cryptsetup`, `go`, `ldconfig`, `mksquashfs`, `nvidia-container-cli`,
  `unsquashfs` are now found at build time by `mconfig` and written into
  `singularity.conf`. The path to these executables can be overridden by
  changing the value in `singularity.conf`. If the path for any of them other
  than `cryptsetup` or `ldconfig` is not set in `singularity.conf` then the
  executable will be found by searching `$PATH`.
- When calling `ldconfig` to find GPU libraries, singularity will *not* fall back
  to `/sbin/ldconfig` if the `ldconfig` on `$PATH` errors. If installing in a
  Guix/Nix on environment on top of a standard host distribution you *must* set
  `ldconfig path = /sbin/ldconfig` to use the host distribution `ldconfig` to
  find GPU libraries.
- Example log-plugin rewritten as a CLI callback that can log all commands
  executed, instead of only container execution, and has access to command
  arguments.
- The bundled reference CNI plugins are updated to v1.0.1. The `flannel` plugin
  is no longer included, as it is maintained as a separate plugin at:
  <https://github.com/flannel-io/cni-plugin>. If you use the flannel CNI plugin
  you should install it from this repository.
- `--nv` will not call `nvidia-container-cli` to find host libraries, unless
  the new experimental GPU setup flow that employs `nvidia-container-cli`
  for all GPU related operations is enabled (see above).
- If a container is run with `--nvccli` and `--contain`, only GPU devices
  specified via the `NVIDIA_VISIBLE_DEVICES` environment variable will be
  exposed within the container. Use `NVIDIA_VISIBLE_DEVICES=all` to access all
  GPUs inside a container run with `--nvccli`.
- Build `--bind` option allows to set multiple bind mount without specifying
  the `--bind` option for each bindings.
- The behaviour of the `allow container` directives in `singularity.conf` has
  been modified, to support more intuitive limitations on the usage of SIF and non-SIF
  container images. If you use these directives, _you may need to make changes
  to singularity.conf to preserve behaviour_.
  - A new `allow container sif` directive permits or denies usage of
    _unencrypted_ SIF images, irrespective of the filesystem(s) inside the SIF.
  - The `allow container encrypted` directive permits or denies usage of SIF
    images with an encrypted root filesystem.
  - The `allow container squashfs/extfs` directives in `singularity.conf`
    permit or deny usage of bare SquashFS and EXT image files only.
  - The effect of the `allow container dir` directive is unchanged.

## v3.8.4 - \[2021-11-09\]

### Bug fixes

- Fix the oras contexts to avoid hangs upon failed pushed to Harbor registry.

### Enhancements

- Added seccomp, cryptsetup, devscripts & correct go version test to
  debian packaging.

Additional changes include dependency updates for the SIF module (to v2.0.0),
and migration to maintained versions of other modules. There is no change to
functionality, on-disk SIF format etc.

## v3.8.3 - \[2021-09-07\]

### Bug fixes

- Fix regression introduced in 3.8.1 that caused bind mounts without a
  destination to be added twice.

## v3.8.2 - \[2021-08-31\]

### Bug fixes

- Fix regression when files `source`d from `%environment` contain `\` escaped
  shell builtins (fixes issue with `source` of conda profile.d script).
- The `oci` commands will operate on systems that use the v2 unified cgroups
  hierarchy.
- `singularity delete` will use the correct library service when the hostname
  is specified in the `library://` URI.
- `singularity build` will use the correct library service when the hostname
  is specified in the `library://` URI / definition file.
- Call `debootstrap` with correct Debian arch when it is not identical to the
  value of `runtime.GOARCH`. E.g. `ppc64el -> ppc64le`.
- When destination is ommitted in `%files` entry in definition file, ensure
  globbed files are copied to correct resolved path.
- Return an error if `--tokenfile` used for `remote login` to an OCI registry,
  as this is not supported.
- Ensure repeated `remote login` to same URI does not create duplicate entries
  in `~/.singularity/remote.yaml`.
- Properly escape single quotes in Docker `CMD` / `ENTRYPOINT` translation.
- Use host uid when choosing unsquashfs flags, to avoid selinux xattr errors
  with `--fakeroot` on non-EL/Fedora distributions with recent squashfs-tools.
- Updated the modified golang-x-crypto module with the latest upstream
  version.

## v3.8.1 - \[2021-08-12\]

### Bug Fixes

- Allow escaped `\$` in a SINGULARITYENV_ var to set a literal `$` in
  a container env var. Also allow escaped commas and colons in the
  source bind path.
- Handle absolute symlinks correctly in multi-stage build `%copy from`
  blocks.
- Fix incorrect reference in sandbox restrictive permissions warning.
- Prevent garbage collection from closing the container image file
  descriptor.
- Update to Arch Linux pacman.conf URL and remove file size verification.
- Avoid panic when mountinfo line has a blank field.

## v3.8.0 - \[2021-06-15\]

### Changed defaults / behaviours

> :warning: Go module was renamed from `github.com/sylabs/singularity` to `github.com/hpcng/singularity`

### New features / functionalities

- A new `overlay` command allows creation and addition of writable overlays.
- Administrators can allow named users/groups to use specific CNI network
  configurations. Managed by directives in `singularity.conf`.
- The `build` command now honors `--nv`, `--rocm`, and `--bind` flags,
  permitting builds that require GPU access or files bound in from the host.
- A library service hostname can be specified as the first component of a
  `library://` URL.
- Singularity is now relocatable for unprivileged installations only.

### Bug Fixes

- Respect http proxy server environment variables in key operations.
- When pushing SIF images to `oras://` endpoints, work around Harbor
  & GitLab failure to accept the `SifConfigMediaType`.
- Avoid a `setfsuid` compilation warning on some gcc versions.
- Fix a crash when silent/quiet log levels used on pulls from
  `shub://` and `http(s)://` URIs.
- Wait for dm device to appear when mounting an encrypted container
  rootfs.
- Accommodate ppc64le pageSize in TestCgroups and disable -race.
- Fix Debian packaging

### Testing / Development

Testing changes are not generally itemized. However, developers and contributors
should note that this release has modified the behavior of `make test` for ease
of use:

- `make test` runs limited unit and integration tests that will not
  require docker hub credentials.
- `make testall` runs the full unit/integration/e2e test suite that
  requires docker credentials to be set with `E2E_DOCKER_USERNAME`
  and `E2E_DOCKER_PASSWORD` environment variables.

## v3.7.4 - \[2021-05-26\]

### Security Related Fixes

- [CVE-2021-32635](https://github.com/hpcng/singularity/security/advisories/GHSA-jq42-hfch-42f3):
  Due to incorrect use of a default URL, singularity action commands
  (run/shell/exec) specifying a container using a library:// URI will
  always attempt to retrieve the container from the default remote
  endpoint (cloud.sylabs.io) rather than the configured remote
  endpoint.  An attacker may be able to push a malicious container to
  the default remote endpoint with a URI that is identical to the URI
  used by a victim with a non-default remote endpoint, thus executing
  the malicious container.

## v3.7.3 - \[2021-04-06\]

### Security Related Fixes

- [CVE-2021-29136](https://github.com/opencontainers/umoci/security/advisories/GHSA-9m95-8hx6-7p9v):
  A dependency used by Singularity to extract docker/OCI image layers can be
  tricked into modifying host files by creating a malicious layer that has a
  symlink with the name "." (or "/"), when running as root. This vulnerability
  affects a `singularity build` or `singularity pull` as root, from a docker or
  OCI source.

## v3.7.2 - \[2021-03-09\]

### Bug Fixes

- Fix progress bar display when source image size is unknown.
- Fix a memory usage / leak issue when building from an existing image file.
- Fix to allow use of `--library` flag to point push/pull at default cloud
  library when another remote is in use.
- Address false positive loop test errors, and an e2e test registry setup issue.

## v3.7.1 - \[2021-01-12\]

### Bug Fixes

- Accommodate /sys/fs/selinux mount changes on kernel 5.9+.
- Fix loop devices file descriptor leak when shared loop devices is enabled.
- Use MaxLoopDevices variable from config file in all appropriate locations.
- Use -buildmode=default (non pie) on ppc64le to prevent crashes when using
  plugins.
- Remove spurious warning in parseTokenSection()
- e2e test fixes for new kernels, new unsquashfs version.
- Show correct web URI for detached builds against alternate remotes.

### New features / functionalities

- The singularity binary is now relocatable when built without setuid support

## v3.7.0 - \[2020-11-24\]

### New features / functionalities

- Allow configuration of global custom keyservers, separate from remote
  endpoints.
- Add a new global keyring, for public keys only (used for ECL).
- The `remote login` command now supports authentication to Docker/OCI
  registries and custom keyservers.
- New `--exclusive` option for `remote use` allows admin to lock usage to a
  specific remote.
- A new `Fingerprints:` header in definition files will check that a SIF source
  image can be verified, and is signed with keys matching all specified
  fingerprints.
- Labels can be set dynamically from a build's `%post` section by setting them
  in the `SINGULARITY_LABELS` environment variable.
- New `build-arch` label is automatically set to the architecture of the host
  during a container build.
- New `-D/--description` flag for `singularity push` sets description for a
  library container image.
- `singularity remote status` shows validity of authentication token if set.
- `singularity push` reports quota usage and URL on successful push to a library
  server that supports this.
- A new `--no-mount` flag for actions allows a user to disable
  proc/sys/dev/devpts/home/tmp/hostfs/cwd mounts, even if they are enabled in
  `singularity.conf`.

### Changed defaults / behaviours

- When actions (run/shell/exec...) are used without `--fakeroot` the umask from
  the calling environment will be propagated into the container, so that files
  are created with expected permissions. Use the new `--no-umask` flag to return
  to the previous behaviour of setting a default 0022 umask.
- Container metadata, environment, scripts are recorded in a descriptor in
  builds to SIF files, and `inspect` will use this if present.
- The `--nv` flag for NVIDIA GPU support will not resolve libraries reported by
  `nvidia-container-cli` via the ld cache. Will instead respect absolute paths
  to libraries reported by the tool, and bind all versioned symlinks to them.
- General re-work of the `remote login` flow, adds prompts and token
  verification before replacing an existing authentication token.
- The Execution Control List (ECL) now verifies container fingerprints using the
  new global keyring. Previously all users would need relevant keys in their own
  keyring.
- The SIF layer mediatype for ORAS has been changed to
  `application/vnd.sylabs.sif.layer.v1.sif` reflecting the published
  [opencontainers/artifacts](https://github.com/opencontainers/artifacts/blob/master/artifact-authors.md#defining-layermediatypes)
  value.
- `SINGULARITY_BIND` has been restored as an environment variable set within a
  running container. It now reflects all user binds requested by the `-B/--bind`
  flag, as well as via `SINGULARITY_BIND[PATHS]`.
- `singularity search` now correctly searches for container images matching the
  host architecture by default. A new `--arch` flag allows searching for other
  architectures. A new results format gives more detail about container image
  results, while users and collections are no longer returned.

### Bug Fixes

- Support larger definition files, environments etc. by passing engine
  configuration in the environment vs. via socket buffer.
- Ensure `docker-daemon:` and other source operations respect
  `SINGULARITY_TMPDIR` for all temporary files.
- Support double quoted filenames in the `%files` section of build definitions.
- Correct `cache list` sizes to show KiB with powers of 1024, matching `du` etc.
- Don't fail on `enable fusemount=no` when no fuse mounts are needed.
- Pull OCI images to the correct requested location when the cache is disabled.
- Ensure `Singularity>` prompt is set when container has no environment script,
  or singularity is called through a wrapper script.
- Avoid build failures in `yum/dnf` operations against the 'setup' package on
  `RHEL/CentOS/Fedora` by ensuring staged `/etc/` files do not match distro
  default content.
- Failed binds to `/etc/hosts` and `/etc/localtime` in a container run with
  `--contain` are no longer fatal errors.
- Don't initialize the cache for actions where it is not required.
- Increase embedded shell interpreter timeout, to allow slow-running environment
  scripts to complete.
- Correct buffer handling for key import to allow import from STDIN.
- Reset environment to avoid `LD_LIBRARY_PATH` issues when resolving
  dependencies for the `unsquashfs` sandbox.
- Fall back to `/sbin/ldconfig` if `ldconfig` on `PATH` fails while resolving
  GPU libraries. Fixes problems on systems using Nix / Guix.
- Address issues caused by error code changes in `unsquashfs` version 4.4.
- Ensure `/dev/kfd` is bound into container for ROCm when `--rocm` is used with
  `--contain`.
- Tolerate comments on `%files` sections in build definition files.
- Fix a loop device file descriptor leak.

### Known Issues

- A change in Linux kernel 5.9 causes `--fakeroot` builds to fail with a
  `/sys/fs/selinux` remount error. This will be addressed in Singularity v3.7.1.

## v3.6.4 - \[2020-10-13\]

### Security related fixes

Singularity 3.6.4 addresses the following security issue.

- [CVE-2020-15229](https://github.com/hpcng/singularity/security/advisories/GHSA-7gcp-w6ww-2xv9):
  Due to insecure handling of path traversal and the lack of path sanitization
  within unsquashfs (a distribution provided utility used by Singularity), it is
  possible to overwrite/create files on the host filesystem during the
  extraction of a crafted squashfs filesystem. Affects unprivileged execution of
  SIF / SquashFS images, and image builds from SIF / SquashFS images.

### Bug Fixes

- Update scs-library-client to support `library://` backends using an 3rd party
  S3 object store that does not strictly conform to v4 signature spec.

## v3.6.3 - \[2020-09-15\]

### Security related fixes

Singularity 3.6.3 addresses the following security issues.

- [CVE-2020-25039](https://github.com/hpcng/singularity/security/advisories/GHSA-w6v2-qchm-grj7):
  When a Singularity action command (run, shell, exec) is run with the fakeroot
  or user namespace option, Singularity will extract a container image to a
  temporary sandbox directory. Due to insecure permissions on the temporary
  directory it is possible for any user with access to the system to read the
  contents of the image. Additionally, if the image contains a world-writable
  file or directory, it is possible for a user to inject arbitrary content into
  the running container.

- [CVE-2020-25040](https://github.com/hpcng/singularity/security/advisories/GHSA-jv9c-w74q-6762):
  When a Singularity command that results in a container build operation is
  executed, it is possible for a user with access to the system to read the
  contents of the image during the build. Additionally, if the image contains a
  world-writable file or directory, it is possible for a user to inject
  arbitrary content into the running build, which in certain circumstances may
  enable arbitrary code execution during the build and/or when the built
  container is run.

## Change defaults / behaviours

- The value for maximum number of loop devices in the config file is now used
  everywhere instead of redefining this value

### Bug Fixes

- Add CAP_MKNOD in capability bounding set of RPC to fix issue with cryptsetup
  when decrypting image from within a docker container.
- Fix decryption issue when using both IPC and PID namespaces.
- Fix unsupported builtins panic from shell interpreter and add umask support
  for definition file scripts.
- Do not load keyring in prepare_linux if ECL not enabled.
- Ensure sandbox option overrides remote build destination.

## v3.6.2 - \[2020-08-25\]

### New features / functionalities

- Add --force option to `singularity delete` for non-interactive workflows.

### Change defaults / behaviours

- Default to current architecture for `singularity delete`.

### Bug Fixes

- Respect current remote for `singularity delete` command.
- Allow `rw` as a (noop) bind option.
- Fix capability handling regression in overlay mount.
- Fix LD_LIBRARY_PATH environment override regression with `--nv/--rocm`.
- Fix environment variable duplication within singularity engine.
- Use `-user-xattrs` for unsquashfs to avoid error with rootless extraction
  using unsquashfs 3.4 (Ubuntu 20.04).
- Correct `--no-home` message for 3.6 CWD behavior.
- Don't fail if parent of cache dir not accessible.
- Fix tests for Go 1.15 Ctty handling.
- Fix additional issues with test images on ARM64.
- Fix FUSE e2e tests to use container ssh_config.

## v3.6.1 - \[2020-07-21\]

### New features / functionalities

- Support compilation with `FORTIFY_SOURCE=2` and build in `pie` mode with
  `fstack-protector` enabled (#5433).

### Bug Fixes

- Provide advisory message r.e. need for `upper` and `work` to exist in overlay
  images.
- Use squashfs mem and processor limits in squashfs gzip check.
- Ensure build destination path is not an empty string - do not overwrite CWD.
- Don't unset PATH when interpreting legacy /environment files.

## v3.6.0 - \[2020-07-14\]

### Security related fixes

Singularity 3.6.0 introduces a new signature format for SIF images, and changes
to the signing / verification code to address:

- [CVE-2020-13845](https://cve.mitre.org/cgi-bin/cvename.cgi?name=2020-13845) In
  Singularity 3.x versions below 3.6.0, issues allow the ECL to be bypassed by a
  malicious user.
- [CVE-2020-13846](https://cve.mitre.org/cgi-bin/cvename.cgi?name=2020-13846) In
  Singularity 3.5 the `--all / -a` option to `singularity verify` returns
  success even when some objects in a SIF container are not signed, or cannot be
  verified.
- [CVE-2020-13847](https://cve.mitre.org/cgi-bin/cvename.cgi?name=2020-13847) In
  Singularity 3.x versions below 3.6.0, Singularity's sign and verify commands
  do not sign metadata found in the global header or data object descriptors of
  a SIF file, allowing an attacker to cause unexpected behavior. A signed
  container may verify successfully, even when it has been modified in ways that
  could be exploited to cause malicious behavior.

Please see the published security advisories at
<https://github.com/hpcng/singularity/security/advisories> for full detail of
these security issues.

Note that the new signature format is necessarily incompatible with Singularity
\< 3.6.0 - e.g. Singularity 3.5.3 cannot verify containers signed by 3.6.0.

We thank Tru Huynh for a report that led to the review of, and changes to, the
signature implementation.

### New features / functionalities

- Singularity now supports the execution of minimal Docker/OCI containers that
  do not contain `/bin/sh`, e.g. `docker://hello-world`.
- A new cache structure is used that is concurrency safe on a filesystem that
  supports atomic rename. *If you downgrade to Singularity 3.5 or older after
  using 3.6 you will need to run `singularity cache clean`.*
- A plugin system rework adds new hook points that will allow the development of
  plugins that modify behavior of the runtime. An image driver concept is
  introduced for plugins to support new ways of handling image and overlay
  mounts. *Plugins built for \<=3.5 are not compatible with 3.6*.
- The `--bind` flag can now bind directories from a SIF or ext3 image into a
  container.
- The `--fusemount` feature to mount filesystems to a container via FUSE drivers
  is now a supported feature (previously an experimental hidden flag). This
  permits users to mount e.g. `sshfs` and `cvmfs` filesystems to the container
  at runtime.
- A new `-c/--config` flag allows an alternative `singularity.conf` to be
  specified by the `root` user, or all users in an unprivileged installation.
- A new `--env` flag allows container environment variables to be set via the
  Singularity command line.
- A new `--env-file` flag allows container environment variables to be set from
  a specified file.
- A new `--days` flag for `cache clean` allows removal of items older than a
  specified number of days. Replaces the `--name` flag which is not generally
  useful as the cache entries are stored by hash, not a friendly name.
- A new '--legacy-insecure' flag to `verify` allows verification of SIF
  signatures in the old, insecure format.
- A new '-l / --logs' flag for `instance list` that shows the paths to instance
  STDERR / STDOUT log files.
- The `--json` output of `instance list` now include paths to STDERR / STDOUT
  log files.

### Changed defaults / behaviours

- New signature format (see security fixes above).
- Environment variables prefixed with `SINGULARITYENV_` always take precedence
  over variables without `SINGULARITYENV_` prefix.
- The `%post` build section inherits environment variables from the base image.
- `%files from ...` will now follow symlinks for sources that are directly
  specified, or directly resolved from a glob pattern. It will not follow
  symlinks found through directory traversal. This mirrors Docker multi-stage
  COPY behaviour.
- Restored the CWD mount behaviour of v2, implying that CWD path is not
  recreated inside container and any symlinks in the CWD path are not resolved
  anymore to determine the destination path inside container.
- The `%test` build section is executed the same manner as
  `singularity test image`.
- `--fusemount` with the `container:` default directive will foreground the FUSE
  process. Use `container-daemon:` for previous behavior.
- Fixed spacing of `singularity instance list` to be dynamically changing based
  off of input lengths instead of fixed number of spaces to account for long
  instance names.

### Deprecated / removed commands

- Removed `--name` flag for `cache clean`; replaced with `--days`.
- Deprecate `-a / --all` option to `sign/verify` as new signature behavior makes
  this the default.

### Bug Fixes

- Don't try to mount `$HOME` when it is `/` (e.g. `nobody` user).
- Process `%appinstall` sections in order when building from a definition file.
- Ensure `SINGULARITY_CONTAINER`, `SINGULARITY_ENVIRONMENT` and the custom shell
  prompt are set inside a container.
- Honor insecure registry settings from `/etc/containers/registries.conf`.
- Fix `http_proxy` env var handling in `yum` bootstrap builds.
- Disable log colorization when output location is not a terminal.
- Check encryption keys are usable before beginning an encrypted build.
- Allow app names with non-alphanumeric characters.
- Use the `base` metapackage for arch bootstrap builds - arch no longer has a
  `base` group.
- Ensure library client messages are logged with `--debug`.
- Do not mount `$HOME` with `--fakeroot --contain`.
- Fall back to underlay automatically when using a sandbox on GPFS.
- Fix Ctrl-Z handling - propagation of signal.

## v3.5.3 - \[2020-02-18\]

### Changed defaults / behaviours

The following minor behaviour changes have been made in 3.5.3 to allow correct
operation on CRAY CLE6, and correct an issue with multi-stage image builds that
was blocking use by build systems such as Spack:

- Container action scripts are no longer bound in from `etc/actions.d` on the
  host. They are created dynamically and inserted at container startup.
- `%files from ...` will no longer follow symlinks when copying between stages
  in a multi stage build, as symlinks should be copied so that they resolve
  identically in later stages. Copying `%files` from the host will still
  maintain previous behavior of following links.

### Bug Fixes

- Bind additional CUDA 10.2 libs when using the `--nv` option without
  `nvidia-container-cli`.
- Fix an NVIDIA persistenced socket bind error with `--writable`.
- Add detection of ceph to allow workarounds that avoid issues with sandboxes on
  ceph filesystems.
- Ensure setgid is inherited during make install.
- Ensure the root directory of a build has owner write permissions, regardless
  of the permissions in the bootstrap source.
- Fix a regression in `%post` and `%test` to honor the `-c` option.
- Fix an issue running `%post` when a container doesn't have `/etc/resolv.conf`
  or `/etc/hosts` files.
- Fix an issue with UID detection on RHEL6 when running instances.
- Fix a logic error when a sandbox image is in an overlay incompatible location,
  and both overlay and underlay are disabled globally.
- Fix an issue causing user namespace to always be used when `allow-setuid=no`
  was configured in a setuid installation.
- Always allow key IDs and fingerprints to be specified with or without a `0x`
  prefix when using `singularity keys`
- Fix an issue preventing joining an instance started with `--boot`.
- Provide a useful error message if an invalid library:// path is provided.
- Bring in multi-part upload client functionality that will address large image
  upload / proxied upload issues with a future update to Sylabs cloud.

In addition, numerous improvements have been made to the test suites, allowing
them to pass cleanly on a range of kernel versions and distributions that are
not covered by the open-source CI runs.

## v3.5.2 - \[2019-12-17\]

### [Security related fix](https://cve.mitre.org/cgi-bin/cvename.cgi?name=2019-19724)

- 700 permissions are enforced on `$HOME/.singularity` and
  `SINGULARITY_CACHEDIR` directories (CVE-2019-19724). Many thanks to Stuart
  Barkley for reporting this issue.

### Bug Fixes

- Fixes an issue preventing use of `.docker/config` for docker registry
  authentication.

- Fixes the `run-help` command in the unprivileged workflow.

- Fixes a regression in the `inspect` command to support older image formats.

- Adds a workaround for an EL6 kernel bug regarding shared bind mounts.

- Fixes caching of http(s) sources with conflicting filenames.

- Fixes a fakeroot sandbox build error on certain filesystems, e.g. lustre,
  GPFS.

- Fixes a fakeroot build failure to a sandbox in $HOME.

- Fixes a fakeroot build failure from a bad def file section script location.

- Fixes container execution errors when CWD is a symlink.

- Provides a useful warning r.e. possible fakeroot build issues when seccomp
  support is not available.

- Fixes an issue where the `--disable-cache` option was not being honored.

- Deprecated `--groupid` flag for `sign` and `verify`; replaced with
  `--group-id`.

- Removed useless flag `--url` for `sign`.

## v3.5.1 - \[2019-12-05\]

### New features / functionalities

A single feature has been added in the bugfix release, with specific
functionality:

- A new option `allow container encrypted` can be set to `no` in
  `singularity.conf` to prevent execution of encrypted containers.

### Bug Fixes

This point release addresses the following issues:

- Fixes a disk space leak when building from docker-archive.
- Makes container process SIGABRT return the expected code.
- Fixes the `inspect` command in unprivileged workflow.
- Sets an appropriate default umask during build stages, to avoid issues with
  very restrictive user umasks.
- Fixes an issue with build script content being consumed from STDIN.
- Corrects the behaviour of underlay with non-empty / symlinked CWD and absolute
  symlink binds targets.
- Fixes execution of containers when binding BTRFS filesystems.
- Fixes build / check failures for MIPS & PPC64.
- Ensures file ownership maintained when building image from sandbox.
- Fixes a squashfs mount error on kernel 5.4.0 and above.
- Fixes an underlay fallback problem, which prevented use of sandboxes on lustre
  filesystems.

## v3.5.0 - \[2019-11-13\]

### New features / functionalities

- New support for AMD GPUs via `--rocm` option added to bind ROCm devices and
  libraries into containers.
- Plugins can now modify Singularity behaviour with two mutators: CLI and
  Runtime.
- Introduced the `config global` command to edit `singularity.conf` settings
  from the CLI.
- Introduced the `config fakeroot` command to setup `subuid` and `subgid`
  mappings for `--fakeroot` from the Singularity CLI.

### Changed defaults / behaviours

- Go 1.13 adopted.
- Vendored modules removed from the Git tree, will be included in release
  tarballs.
- Singularity will now fail with an error if a requested bind mount cannot be
  made.
  - This is beneficial to fail fast in workflows where a task may fail a long
    way downstream if a bind mount is unavailable.
  - Any unavailable bind mount sources must be removed from `singularity.conf`.
- Docker/OCI image extraction now faithfully respects layer permissions.
  - This may lead to sandboxes that cannot be removed without modifying
    permissions.
  - `--fix-perms` option added to preserve old behaviour when building
    sandboxes.
  - Discussion issue for this change at:
    <https://github.com/sylabs/singularity/issues/4671>
- `Singularity>` prompt is always set when entering shell in a container.
- The current `umask` will be honored when building a SIF file.
- `instance exec` processes acquire cgroups set on `instance start`
- `--fakeroot` supports uid/subgid ranges >65536
- `singularity version` now reports semver compliant version information.

### Deprecated / removed commands

- Deprecated `--id` flag for `sign` and `verify`; replaced with `--sif-id`.

## v3.4.2 - \[2019-10-08\]

- This point release addresses the following issues:
  - Sets workable permissions on OCI -> sandbox rootless builds
  - Fallback correctly to user namespace for non setuid installation
  - Correctly handle the starter-suid binary for non-root installs
  - Creates CACHEDIR if it doesn't exist
  - Set apex loglevel for umoci to match singularity loglevel

## v3.4.1 - \[2019-09-17\]

- This point release addresses the following issues:
  - Fixes an issue where a PID namespace was always being used
  - Fixes compilation on non 64-bit architectures
  - Allows fakeroot builds for zypper, pacstrap, and debootstrap
  - Correctly detects seccomp on OpenSUSE
  - Honors GO_MODFLAGS properly in the mconfig generated makefile
  - Passes the Mac hostname to the VM in MacOS Singularity builds
  - Handles temporary EAGAIN failures when setting up loop devices on recent
    kernels
  - Fixes excessive memory usage in singularity push

## v3.4.0 - \[2019-08-30\]

### New features / functionalities

- New support for building and running encrypted containers with RSA keys and
  passphrases
  - `--pem-path` option added to the `build` and action commands for RSA based
    encrypted containers
  - `--passphrase` option added to `build` and action commands for passphrase
    based encrypted containers
  - `SINGULARITY_ENCRYPTION_PEM_PATH` and `SINGULARITY_ENCRYPTION_PASSPHRASE`
    environment variables added to serve same functions as above
  - `--encrypt` option added to `build` command to build an encrypted container
    when environment variables contain a secret
- New `--disable-cache` flag prevents caching of downloaded containers
- Added support for multi-line variables in singularity def-files
- Added support for 'indexed' def-file variables (like arrays)
- Added support for SUSE SLE Products
- Added the def-file variables: product, user, regcode, productpgp, registerurl,
  modules, otherurl (indexed)
- Support multiple-architecture tags in the SCS library
- Added a `--dry-run` flag to `cache clean`
- Added a `SINGULARITY_SYPGPDIR` environment variable to specify the location of
  PGP key data
- Added a `--nonet` option to the action commands to disable networking when
  running with the `--vm` option
- Added a `--long-list` flag to the `key search` command to preserve
- Added experimental, hidden `--fusemount` flag to pass a command to mount a
  libfuse3 based file system within the container

### Changed defaults / behaviors

- Runtime now properly honors `SINGULARITY_DISABLE_CACHE` environment variable
- `remote add` command now automatically attempts to login and a `--no-login`
  flag is added to disable this behavior
- Using the `pull` command to download an unsigned container no longer produces
  an error code
- `cache clean` command now prompts user before cleaning when run without
  `--force` option and is more verbose
- Shortened the default output of the `key search` command

### Deprecated / removed commands

- The `--allow-unsigned` flag to `pull` has been deprecated and will be removed
  in the future

## v3.3.0 - \[2019-06-17\]

### Changed defaults / behaviors

- Remote login and status commands will now use the default remote if a remote
  name is not supplied
- Added Singularity hub (`shub`) cache support when using the `pull` command
- Clean cache in a safer way by only deleting the cache subdirectories
- Improvements to the `cache clean` command

### New features / functionalities

- new `oras` URI for pushing and pulling SIF files to and from supported OCI
  registries
- added the `--fakeroot` option to `build`, `exec`, `run`, `shell`, `test`, and
  `instance start` commands to run container in a new user namespace as uid 0
- added the `fakeroot` network type for use with the `--network` option
- `sif` command to allow for the inspection and manipulation of SIF files with
  the following subcommands
  - `add` Add a data object to a SIF file
  - `del` Delete a specified object descriptor and data from SIF file
  - `dump` Extract and output data objects from SIF files
  - `header` Display SIF global headers
  - `info` Display detailed information of object descriptors
  - `list` List object descriptors from SIF files
  - `new` Create a new empty SIF image file
  - `setprim` Set primary system partition

## v3.2.1 - \[2019-05-28\]

- This point release fixes the following bugs:
  - Allows users to join instances with non-suid workflow
  - Removes false warning when seccomp is disabled on the host
  - Fixes an issue in the terminal when piping output to commands
  - Binds NVIDIA persistenced socket when `--nv` is invoked

## v3.2.0 - \[2019-05-14\]

### [Security related fix](https://cve.mitre.org/cgi-bin/cvename.cgi?name=2019-11328)

- Instance files are now stored in user's home directory for privacy and many
  checks have been added to ensure that a user can't manipulate files to change
  `starter-suid` behavior when instances are joined (many thanks to Matthias
  Gerstner from the SUSE security team for finding and securely reporting this
  vulnerability)

### New features / functionalities

- Introduced a new basic framework for creating and managing plugins
- Added the ability to create containers through multi-stage builds
  - Definitions now require `Bootstrap` be the first parameter of header
- Created the concept of a Sylabs Cloud "remote" endpoint and added the ability
  for users and admins to set them through CLI and conf files
- Added caching for images from Singularity Hub
- Made it possible to compile Singularity outside of `$GOPATH`
- Added a json partition to SIF files for OCI configuration when building from
  an OCI source
- Full integration with Singularity desktop for MacOS code base

### New Commands

- Introduced the `plugin` command group for creating and managing plugins

  - `compile` Compile a singularity plugin
  - `disable` disable an installed singularity plugin
  - `enable` Enable an installed singularity plugin
  - `inspect` Inspect a singularity plugin (either an installed one or an image)
  - `install` Install a singularity plugin
  - `list` List installed singularity plugins
  - `uninstall` Uninstall removes the named plugin from the system

- Introduced the `remote` command group to support management of Singularity
  endpoints:

  - `add` Create a new Sylabs Cloud remote endpoint
  - `list` List all remote endpoints that are configured
  - `login` Log into a remote endpoint using an authentication token
  - `remove` Remove an existing Sylabs Cloud remote endpoint
  - `status` Check the status of the services at an endpoint
  - `use` Set a remote endpoint to be used by default

- Added to the `key` command group to improve PGP key management:

  - `export` Export a public or private key into a specific file
  - `import` Import a local key into the local keyring
  - `remove` Remove a local public key

- Added the `Stage: <name>` keyword to the definition file header and the
  `from <stage name>` option/argument pair to the `%files` section to support
  multistage builds

### Deprecated / removed commands

- The `--token/-t` option has been deprecated in favor of the
  `singularity remote` command group

### Changed defaults / behaviors

- Ask to confirm password on a newly generated PGP key
- Prompt to push a key to the KeyStore when generated
- Refuse to push an unsigned container unless overridden with
  `--allow-unauthenticated/-U` option
- Warn and prompt when pulling an unsigned container without the
  `--allow-unauthenticated/-U` option
- `Bootstrap` must now be the first field of every header because of parser
  requirements for multi-stage builds

## v3.1.1 - \[2019-04-02\]

### New Commands

- New hidden `buildcfg` command to display compile-time parameters
- Added support for `LDFLAGS`, `CFLAGS`, `CGO_` variables in build system
- Added `--nocolor` flag to Singularity client to disable color in logging

### Removed Commands

- `singularity capability <add/drop> --desc` has been removed
- `singularity capability list <--all/--group/--user>` flags have all been
  removed

### New features / functionalities

- The `--builder` flag to the `build` command implicitly sets `--remote`
- Repeated binds no longer cause Singularity to exit and fail, just warn instead
- Corrected typos and improved docstrings throughout
- Removed warning when CWD does not exist on the host system
- Added support to spec file for RPM building on SLES 11

## v3.1.0 - \[2019-02-22\]

### New Commands

- Introduced the `oci` command group to support a new OCI compliant variant of
  the Singularity runtime:
  - `attach` Attach console to a running container process
  - `create` Create a container from a bundle directory
  - `delete` Delete container
  - `exec` Execute a command within container
  - `kill` Kill a container
  - `mount` Mount create an OCI bundle from SIF image
  - `pause` Suspends all processes inside the container
  - `resume` Resumes all processes previously paused inside the container
  - `run` Create/start/attach/delete a container from a bundle directory
  - `start` Start container process
  - `state` Query state of a container
  - `umount` Umount delete bundle
  - `update` Update container cgroups resources
- Added `cache` command group to inspect and manage cached files
  - `clean` Clean your local Singularity cache
  - `list` List your local Singularity cache

### New features / functionalities

- Can now build CLI on darwin for limited functionality on Mac
- Added the `scratch` bootstrap agent to build from anything
- Reintroduced support for zypper bootstrap agent
- Added the ability to overwrite a new `singularity.conf` when building from RPM
  if desired
- Fixed several regressions and omissions in [SCIF](https://sci-f.github.io/)
  support
- Added caching for containers pulled/built from the
  [Container Library](https://cloud.sylabs.io/library)
- Changed `keys` command group to `key` (retained hidden `keys` command for
  backward compatibility)
- Created an `RPMPREFIX` variable to allow RPMs to be installed in custom
  locations
- Greatly expanded CI unit and end-to-end testing

## v3.0.3 - \[2019-01-21\]

- Bind paths in `singularity.conf` are properly parsed and applied at runtime
- Singularity runtime will properly fail if `singularity.conf` file is not owned
  by the root user
- Several improvements to RPM packaging including using golang from epel,
  improved support for Fedora, and avoiding overwriting conf file on new RPM
  install
- Unprivileged `--contain` option now properly mounts `devpts` on older kernels
- Uppercase proxy environment variables are now rightly respected
- Add http/https protocols for singularity run/pull commands
- Update to SIF 1.0.2
- Add _noPrompt_ parameter to `pkg/signing/Verify` function to enable silent
  verification

## v3.0.2 - \[2019-01-04\]

- Added the `--docker-login` flag to enable interactive authentication with
  docker registries
- Added support for pulling directly from HTTP and HTTPS
- Made minor improvements to RPM packaging and added basic support for alpine
  packaging
- The `$SINGULARITY_NOHTTPS`,`$SINGULARITY_TMPDIR`, and
  `$SINGULARITY_DOCKER_USERNAME`/`$SINGULARITY_DOCKER_PASSWORD` environment
  variables are now correctly respected
- Pulling from a private shub registry now works as expected
- Running a container with `--network="none"` no longer incorrectly fails with
  an error message
- Commands now correctly return 1 when incorrectly executed without arguments
- Progress bars no longer incorrectly display when running with `--quiet` or
  `--silent`
- Contents of `91-environment.sh` file are now displayed if appropriate when
  running `inspect --environment`

## v3.0.1 - \[2018-10-31\]

- Improved RPM packaging procedure via makeit
- Enhanced general stability of runtime

## v3.0.0 - \[2018-10-08\]

- Singularity is now written primarily in Go to bring better integration with
  the existing container ecosystem
- Added support for new URIs (`build` & `run/exec/shell/start`):
  - `library://` - Supports the
    [Sylabs.io Cloud Library](https://cloud.sylabs.io/library)
  - `docker-daemon:` - Supports images managed by the locally running docker
    daemon
  - `docker-archive:` - Supports archived docker images
  - `oci:` - Supports oci images
  - `oci-archive:` - Supports archived oci images
- Handling of `docker` & `oci` URIs/images now utilizes
  [containers/image](https://github.com/containers/image) to parse and convert
  those image types in a supported way
- Replaced `singularity instance.*` command group with `singularity instance *`
- The command `singularity help` now only provides help regarding the usage of
  the `singularity` command. To display an image's `help` message, use
  `singularity run-help <image path>` instead

### Removed Deprecated Commands

- Removed deprecated `singularity image.*` command group
- Removed deprecated `singularity create` command
- Removed deprecated `singularity bootstrap` command
- Removed deprecated `singularity mount` command
- Removed deprecated `singularity check` command

### New Commands

- Added `singularity run-help <image path>` command to output an image's `help`
  message
- Added `singularity sign <image path>` command to allow a user to
  cryptographically sign a SIF image
- Added `singularity verify <image path>` command to allow a user to verify a
  SIF image's cryptographic signatures
- Added `singularity keys` command to allow the management of `OpenPGP` key
  stores
- Added `singularity capability` command to allow fine grained control over the
  capabilities of running containers
- Added `singularity push` command to push images to the
  [Sylabs.io Cloud Library](https://cloud.sylabs.io/library)

### Changed Commands

#### Action Command Group (`run/shell/exec/instance start`)

- Added flags:
  - `--add-caps <string>`: Run the contained process with the specified
    capability set (requires root)
  - `--allow-setuid`: Allows setuid binaries to be mounted into the container
    (requires root)
  - `--apply-cgroups <path>`: Apply cgroups configuration from file to contained
    processes (requires root)
  - `--dns <string>`: Adds the comma separated list of DNS servers to the
    containers `resolv.conf` file
  - `--drop-caps <string>`: Drop the specified capabilities from the container
    (requires root)
  - `--fakeroot`: Run the container in a user namespace as `uid=0`. Requires a
    recent kernel to function properly
  - `--hostname <string>`: Set the hostname of the container
  - `--keep-privs`: Keep root user privilege inside the container (requires
    root)
  - `--network <string>`: Specify a list of comma separated network types
    ([CNI Plugins](https://github.com/containernetworking/cni)) to be present
    inside the container, each with its own dedicated interface in the container
  - `--network-args <string>`: Specify arguments to pass to CNI network plugins
    (set by `--network`)
  - `--no-privs`: Drop all privileges from root user inside the container
    (requires root)
  - `--security <string>`: Configure security features such as SELinux,
    Apparmor, Seccomp...
  - `--writable-tmpfs`: Run container with a `tmpfs` overlay
- The command `singularity instance start` now supports the `--boot` flag to
  boot the container via `/sbin/init`
- Changes to image mounting behavior:
  - All image formats are mounted as read only by default
  - `--writable` only works on images which can be mounted in read/write
    \[applicable to: `sandbox` and legacy `ext3` images\]
  - `--writable-tmpfs` runs the container with a writable `tmpfs`-based overlay
    \[applicable to: all image formats\]
  - `--overlay <string>` now specifies a list of `ext3`/`sandbox` images which
    are set as the containers overlay \[applicable to: all image formats\]

#### Build Command

- All images are now built as
  [Singularity Image Format (SIF)](https://www.sylabs.io/2018/03/sif-containing-your-containers/)
  images by default
- When building to a path that already exists, `singularity build` will now
  prompt the user if they wish to overwrite the file existing at the specified
  location
- The `-w|--writable` flag has been removed
- The `-F|--force` flag now overrides the interactive prompt and will always
  attempt to overwrite the file existing at the specified location
- The `-u|--update` flag has been added to support the workflow of running a
  definition file on top of an existing container \[implies `--sandbox`, only
  supports `sandbox` image types\]
- The `singularity build` command now supports the following flags for
  integration with the
  [Sylabs.io Cloud Library](https://cloud.sylabs.io/library):
  - `-r|--remote`: Build the image remotely on the Sylabs Remote Builder
    (currently unavailable)
  - `-d|--detached`: Detach from the `stdout` of the remote build \[requires
    `--remote`\]
  - `--builder <string>`: Specifies the URL of the remote builder to access
  - `--library <string>`: Specifies the URL of the
    [Sylabs.io Cloud Library](https://cloud.sylabs.io/library) to push the built
    image to when the build command destination is in the form
    `library://<reference>`
- The `bootstrap` keyword in the definition file now supports the following
  values:
  - `library`
  - `docker-daemon`
  - `docker-archive`
  - `oci`
  - `oci-archive`
- The `from` keyword in the definition file now correctly parses a `docker` URI
  which includes the `registry` and/or `namespace` components
- The `registry` and `namespace` keywords in the definition file are no longer
  supported. Instead, those values may all go into the `from` keyword
- Building from a tar archive of a `sandbox` no longer works
# Contributors to Singularity

## Project Lead

```text
- Gregory M. Kurtzer <gmk@resf.org>, <gmkurtzer@gmail.com>, <gmk@ctrliq.com>
```

## Maintainers

```text
- Cedric Clerget <cedric@ctrliq.com>, <cedric.clerget@univ-fcomte.fr>
- Ian Kaneshiro <ian@ctrliq.com>, <iankane@umich.edu>
- Krishna Muriki <kmuriki@lbl.gov>, <kmuriki@gmail.com>
- Dave Dykstra <dwd@fnal.gov>
```

## Contributors

```text
- Adam Hughes <adam@sylabs.io>, <stickmanica@gmail.com>
- Adam Simpson <asimpson@nvidia.com>, <adambsimpson@gmail.com>
- Afif Elghraoui <afif.elghraoui@nih.gov>
- Amanda Duffy <aduffy@lenovo.com>
- Ana Guerrero Lopez <aguerrero@suse.com>
- Ángel Bejarano <abejarano@ontropos.com>
- Apuã Paquola <apuapaquola@gmail.com>
- Aron Öfjörð Jóhannesson <aron1991@gmail.com>
- Bernard Li <bernardli@lbl.gov>
- Brian Bockelman <bbockelm@cse.unl.edu>
- Carl Madison <carl@sylabs.io>
- Chris Burr <christopher.burr@cern.ch>
- Chris Hollowell <hollowec@bnl.gov>
- Christian Neyers <foss@neyers.org>
- Daniele Tamino <daniele.tamino@gmail.com>
- Dave Godlove <d@sylabs.io>, <davidgodlove@gmail.com>
- Dave Love <d.love@liverpool.ac.uk>
- David Trudgian <david.trudgian@utsouthwestern.edu>,
  <david.trudgian@sylabs.io>, <dave@trudgian.net>
- Diana Langenbach <dcl@dcl.sh>
- Dimitri Papadopoulos <3234522+DimitriPapadopoulos@users.noreply.github.com>
- Divya Cote <divya.cote@gmail.com>
- Eduardo Arango <eduardo@sylabs.io>, <arangogutierrez@gmail.com>
- Egbert Eich <eich@suse.com>
- Eric Müller <mueller@kip.uni-heidelberg.de>
- Felix Abecassis <fabecassis@nvidia.com>
- Geoffroy Vallee <geoffroy@sylabs.io>, <geoffroy.vallee@gmail.com>
- George Hartzell <hartzell@alerce.com>
- Gert Hulselmans <gert.hulselmans@kuleuven.vib.be>
- Götz Waschk <goetz.waschk@desy.de>
- Hakon Enger <hakonenger@github.com>
- Hugo Meiland <hugo.meiland@microsoft.com>
- Jack Morrison <morrisonjc@ornl.gov>, <jack@rescale.com>
- Jacob Chappell <chappellind@gmail.com>, <jacob.chappell@uky.edu>
- Jarrod Johnson <jjohnson2@lenovo.com>
- Jason Stover <jms@sylabs.io>, <jason.stover@gmail.com>
- Jeff Kriske <jekriske@gmail.com>
- Jia Li <jiali@sylabs.io>
- Joana Chavez <joana@sylabs.io>, <j.chavezlavalle@gmail.com>
- Josef Hrabal <josef.hrabal@vsb.cz>
- Justin Cook <justin@sylabs.io>
- Justin Riley <justin_riley@harvard.edu>
- Krishna Muriki <kmuriki@lbl.gov>
- Kumar Sukhani <kumarsukhani@gmail.com>
- Kundan Kumar <iamkundankumar28@gmail.com>
- Maciej Sieczka <msieczka@sieczka.org>
- Marcelo Magallon <marcelo@sylabs.io>
- Mark Egan-Fuller <markeganfuller@googlemail.com>
- Matt Wiens <mwiens91@gmail.com>
- Michael Bauer <m@sylabs.io>, <bauerm@umich.edu>
- Michael Herzberg <michael@mherzberg.de>
- Michael Milton <ttmigueltt@gmail.com>
- Michael Moore <michael.moore@nuance.com>
- Mike Frisch <michael.frisch@sylabs.io>
- Mike Gray <mike@sylabs.io>
- Nathan Chou <nathan.chou@sylabs.io>, <choun@berkeley.edu>
- Nathan Lin <nathan.lin@yale.edu>
- Oleksandr Moskalenko <om@rc.ufl.edu>
- Oliver Breitwieser <obreitwi@kip.uni-heidelberg.de>, <oliver@breitwieser.eu>
- Oliver Freyermuth <freyermuth@physik.uni-bonn.de>
- Olivier Sallou <olivier.sallou@irisa.fr>
- Peter Steinbach <steinbach@scionics.de>
- Petr Votava <votava.petr@gene.com>
- Rafal Gumienny <rafal.gumienny@gmail.com>
- Ralph Castain <rhc@open-mpi.org>
- Rémy Dernat <remy.dernat@umontpellier.fr>
- Richard Neuboeck <hawk@tbi.univie.ac.at>
- Sasha Yakovtseva <sasha@sylabs.io>, <sashayakovtseva@gmail.com>
- Satish Chebrolu  <satish@sylabs.io>
- Shane Loretz <sloretz@openrobotics.org>, <shane.loretz@gmail.com>
- Shengjing Zhu <i@zhsj.me>
- Tarcisio Fedrizzi <tarcisio.fedrizzi@gmail.com>
- Thomas Hamel <hmlth@t-hamel.fr>
- Tim Wright <7im.Wright@protonmail.com>
- Tru Huynh <tru@pasteur.fr>
- Tyson Whitehead <twhitehead@gmail.com>
- Vanessa Sochat <vsochat@stanford.edu>
- Westley Kurtzer <westley@sylabs.io>, <westleyk@nym.hush.com>
- Yannick Cote <y@sylabs.io>, <yhcote@gmail.com>
- Yaroslav Halchenko <debian@onerussian.com>
- Onur Yılmaz <csonuryilmaz@gmail.com>
- Pranathi Locula <locula@deshaw.com>
- Pedro Alves Batista <pedro.pesquisapb@gmail.com>
```
# Singularity

IMPORTANT NOTE: Singularity is being renamed to
[Apptainer](https://apptainer.org).
This repository is now only for maintaining the 3.8 series and archiving
the history.
Submit all current issues and pull requests to
[https://github.com/apptainer/apptainer](https://github.com/apptainer/apptainer).

[![CI](https://github.com/hpcng/singularity/actions/workflows/ci.yml/badge.svg)](https://github.com/hpcng/singularity/actions/workflows/ci.yml)

- [Documentation](https://singularity.hpcng.org/docs/)
- [Support](#support)
- [Community Meetings / Minutes / Roadmap](https://drive.google.com/drive/u/0/folders/1npfBhIDxqeJIUHZ0tMeuHPvc_iB4T2B6)
- [Project License](LICENSE.md)
- [Guidelines for Contributing](CONTRIBUTING.md)
- [Code of Conduct](CODE_OF_CONDUCT.md)
- [Citation](#citing-singularity)

## What is Singularity?

Singularity is an open source container platform designed to be simple, fast,
and secure. Many container platforms are available, but Singularity is designed
for ease-of-use on shared systems and in high performance computing (HPC)
environments. It features:

- An immutable single-file container image format, supporting cryptographic
  signatures and encryption.
- Integration over isolation by default. Easily make use of GPUs, high speed
  networks, parallel filesystems on a cluster or server.
- Mobility of compute. The single file SIF container format is easy to transport
  and share.
- A simple, effective security model. You are the same user inside a container
  as outside, and cannot gain additional privilege on the host system by
  default.

Singularity is open source software, distributed under the [BSD License](LICENSE.md).

Check out [talks about Singularity](https://singularity.hpcng.org/talks)
and some [use cases of Singularity](https://singularity.hpcng.org/usecases)
on our website.

## Getting Started with Singularity

To install Singularity from source, see the [installation
instructions](INSTALL.md). For other installation options, see [our
guide](https://singularity.hpcng.org/admin-docs/master/installation.html).

System administrators can learn how to configure Singularity, and get an
overview of its architecture and security features in the [administrator
guide](https://singularity.hpcng.org/admin-docs/master/).

For users, see the [user guide](https://singularity.hpcng.org/user-docs/master/)
for details on how to run and build containers with Singularity.

## Contributing to Singularity

Community contributions are always greatly appreciated. To start developing
Singularity, check out the [guidelines for contributing](CONTRIBUTING.md).

Please note we have a [code of conduct](CODE_OF_CONDUCT.md). Please follow it in
all your interactions with the project members and users.

Our roadmap, other documents, and user/developer meeting information can be
found in the [singularity community page](https://singularity.hpcng.org/help).

We also welcome contributions to our [user
guide](https://github.com/hpcng/singularity-userdocs) and [admin
guide](https://github.com/hpcng/singularity-admindocs).

## Support

To get help with Singularity, check out the [Singularity
Help](https://singularity.hpcng.org/help) web page.

## Go Version Compatibility

Singularity aims to maintain support for the two most recent stable versions
of Go. This corresponds to the Go
[Release Maintenance
Policy](https://github.com/golang/go/wiki/Go-Release-Cycle#release-maintenance)
and [Security Policy](https://golang.org/security),
ensuring critical bug fixes and security patches are available for all
supported language versions.

## Citing Singularity

The Singularity software may be cited using our Zenodo DOI `10.5281/zenodo.1310023`:

> Singularity Developers (2021) Singularity. 10.5281/zenodo.1310023
> <https://doi.org/10.5281/zenodo.1310023>

This is an 'all versions' DOI for referencing Singularity in a manner that is
not version-specific. You may wish to reference the particular version of
Singularity used in your work. Zenodo creates a unique DOI for each release,
and these can be found in the 'Versions' sidebar on the [Zenodo record page](https://doi.org/10.5281/zenodo.1310023).

Please also consider citing the original publication describing Singularity:

> Kurtzer GM, Sochat V, Bauer MW (2017) Singularity: Scientific containers for
> mobility of compute. PLoS ONE 12(5): e0177459.
> <https://doi.org/10.1371/journal.pone.0177459>

## License

_Unless otherwise noted, this project is licensed under a 3-clause BSD license
found in the [license file](LICENSE.md)._
# Release Procedure

The release procedure below can be performed by a project member with
"maintainer" or higher privileges on the GitHub repository. It assumes
that you will be working in an up-to-date local clone of the GitHub
repository, where the `upstream` remote points to
`github.com/hpcng/singularity`.

## Prior to Release

1. Set a target date for the release candidate (if required) and release.
   Generally 2 weeks from RC -> release is appropriate for new 3.X.0 minor
   versions.
1. Aim to specifically discuss the release timeline and progress in community
   meetings at least 2 months prior to the scheduled date.
1. Use a GitHub milestone to track issues and PRs that will form part of the
   next release.
1. Ensure that the `CHANGELOG.md` is kept up-to-date on the `master` branch,
   with all relevant changes listed under a "Changes Since Last Release"
   section.
1. Monitor and merge dependabot updates, such that a release is made with as
   up-to-date versions of dependencies as possible. This lessens the burden in
   addressing patch release fixes that require dependency updates, as we use
   several dependencies that move quickly.

## Creating the Release Branch and Release Candidate

When a new 3.Y.0 minor version of Singularity is issued the release
process begins by branching, and then issuing a release candidate for
broader testing.

When a new 3.Y.Z patch release is issued, the branch will already be present,
and steps 1-2 should be skipped.

1. From a repository that is up-to-date with master, create a release
   branch e.g. `git checkout upstream/master -b release-3.8`.
1. Push the release branch to GitHub via `git push upstream release-3.8`.
1. Examine the GitHub branch protection rules, to extend them to the
   new release branch if needed.
1. Modify the `README.md`, `INSTALL.md`, `CHANGELOG.md` via PR against
   the release branch, so that they reflect the version to be released.
1. Apply an annotated tag via `git tag -a -m "Singularity v3.8.0
   Release Candidate 1" v3.8.0-rc.1`.
1. Push the tag via `git push upstream v3.8.0-rc.1`.
1. Create a tarball via `./mconfig --only-rpm -v && make dist`.
1. Test intallation from the tarball.
1. Compute the sha256sum of the tarball e.g. `sha256sum *.tar.gz > sha256sums`.
1. Create a GitHub release, marked as a 'pre-release', incorporating
   `CHANGELOG.md` information, and attaching the tarball and
   `sha256sums`.
1. Notify the community about the RC via the Google Group and Slack.

There will often be multiple release candidates issued prior to the final
release of a new 3.Y.0 minor version.

A small 3.Y.Z patch release may not require release candidates where the code
changes are contained, confirmed by the person reporting the bug(s), and well
covered by tests.

## Creating a Final Release

1. Ensure the user and admin documentation is up-to-date for the new
   version, branched, and tagged.
   - [User Docs](https://singularity.hpcng.org/user-docs/master/) can be
     edited [here](https://github.com/hpcng/singularity-userdocs)
   - [Admin Docs](https://singularity.hpcng.org/admin-docs/master/) can be
     edited [here](https://github.com/hpcng/singularity-admindocs)
1. Ensure the user and admin documentation has been deployed to the
   singularity.hpcng.org website.
1. Modify the `README.md`, `INSTALL.md`, `CHANGELOG.md` via PR against
   the release branch, so that they reflect the version to be released.
1. Apply an annotated tag via `git tag -a -m "Singularity v3.8.0" v3.8.0`.
1. Push the tag via `git push upstream v3.8.0-rc.1`.
1. Create a tarball via `./mconfig -v && make dist`.
1. Test intallation from the tarball.
1. Compute the sha256sum of the tarball e.g. `sha256sum *.tar.gz > sha256sums`.
1. Create a GitHub release, incorporating `CHANGELOG.md` information,
   and attaching the tarball and `sha256sums`.
1. Notify the community about the RC via the Google Group and Slack.

## After the Release

1. Create and merge a PR from the `release-3.x` branch into `master`, so that
   history from the RC process etc. is captured on `master`.
1. If the release is a new major/minor version, move the prior `release-3.x`
   branch to `vault/release-3.x`.
1. If the release is a new major/minor version, update the
   `.github/dependabot.yml` configuration so that dependabot is tracking the new
   stable release branch.
1. Start scheduling / setting up milestones etc. to track the next release!
# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of
experience, nationality, personal appearance, race, religion, or sexual identity
and orientation.

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

Project maintainers have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, or to ban temporarily or permanently any
contributor for other behaviors that they deem inappropriate, threatening,
offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project leader (gmkurtzer@gmail.com). All
complaints will be reviewed and investigated and will result in a
response that is deemed necessary and appropriate to the circumstances.
The project team is obligated to maintain confidentiality with regard to
the reporter of an incident. Further details of specific enforcement
policies may be posted separately.

Project maintainers, contributors and users who do not follow or enforce the
Code of Conduct in good faith may face temporary or permanent repercussions with
their involvement in the project as determined by the project's leader(s).

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing

## Contributor's Agreement

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the project, without imposing a
separate written license agreement for such Enhancements, then you hereby grant
the following license: a non-exclusive, royalty-free perpetual license to
install, use, modify, prepare derivative works, incorporate into other computer
software, distribute, and sublicense such enhancements or derivative works
thereof, in binary and source code form.

## Getting Started

When contributing to Singularity, it is important to properly communicate the
gist of the contribution. If it is a simple code or editorial fix, simply
explaining this within the GitHub Pull Request (PR) will suffice. But if this is
a larger fix or Enhancement, you are advised to first discuss the change with
the project leader or developers.

Please note we have a [code of conduct](CODE_OF_CONDUCT.md). Please follow it in
all your interactions with the project members and users.

## Pull Requests (PRs)

1. Essential bug fix PRs should be sent to both master and release branches.
1. Small bug fix and feature enhancement PRs should be sent to master only.
1. Follow the existing code style precedent, especially for C. For Go, you
   will mostly conform to the style and form enforced by the "go fmt" and
   "golint" tools for proper formatting.
1. For any new functionality, please write appropriate go tests that will run as
   part of the Continuous Integration (github workflow actions) system.
1. Make sure that the project's default copyright and header have been included
   in any new source files.
1. Make sure your code passes linting, by running `make check` before submitting
   the PR. We use `golangci-lint` as our linter. You may need to address linting
   errors by:
   - Running `gofumpt -w .` to format all `.go` files. We use
     [gofumpt](https://github.com/mvdan/gofumpt) instead of `gofmt` as it adds
     additional formatting rules which are helpful for clarity.
   - Leaving a function comment on **every** new exported function and package
     that your PR has introduced. To learn about how to properly comment Go
     code, read
     [this post on golang.org](https://golang.org/doc/effective_go.html#commentary)
1. Make sure you have locally tested using `make -C builddir test` and that all
   tests succeed before submitting the PR.
1. If possible, run `make -C builddir testall` locally, after setting the
   environment variables `E2E_DOCKER_USERNAME` and `E2E_DOCKER_PASSWORD`
   appropriately for an authorized Docker Hub account. This is required as
   Singularity's end-to-end tests perform many tests that build from or execute
   docker images. Our CI is authorized to run these tests if you cannot.
1. Ask yourself is the code human understandable? This can be accomplished via a
   clear code style as well as documentation and/or comments.
1. The pull request will be reviewed by others, and finally merged when all
   requirements are met.
1. The `CHANGELOG.md` must be updated for any of the following changes:
   - Renamed commands
   - Deprecated / removed commands
   - Changed defaults / behaviors
   - Backwards incompatible changes
   - New features / functionalities
1. PRs which introduce a new Go dependency to the project via `go get` and
   additions to `go.mod` should explain why the dependency is required.

## Documentation

There are a few places where documentation for the Singularity project lives.
The [changelog](CHANGELOG.md) is where PRs should include documentation if
necessary. When a new release is tagged, the
[user-docs](https://singularity.hpcng.org/user-docs/master/) and
[admin-docs](https://singularity.hpcng.org/admin-docs/master/) will be updated
using the contents of the `CHANGELOG.md` file as reference.

1. The [changelog](CHANGELOG.md) is a place to document **functional**
   differences between versions of Singularity. PRs which require
   documentation must update this file. This should be a document which can be
   used to explain what the new features of each version of Singularity are,
   and should **not** read like a commit log. Once a release is tagged (*e.g.
   v3.0.0*), a new top level section will be made titled **Changes Since
   vX.Y.Z** (*e.g. Changes Since v3.0.0*) where new changes will now be
   documented, leaving the previous section immutable.
1. The [README](README.md) is a place to document critical information for new
   users of Singularity. It should typically not change, but in the case where
   a change is necessary a PR may update it.
1. The [user-docs](https://singularity.hpcng.org/user-docs/master/) should
   document anything pertinent to the usage of Singularity.
1. The [admin-docs](https://singularity.hpcng.org/admin-docs/master/)
   document anything that is pertinent to a system administrator who manages a
   system with Singularity installed.
1. If necessary, changes to the message displayed when running
   `singularity help *` can be made by editing `docs/content.go`.
# Installing Singularity

Since you are reading this from the Singularity source code, it will be assumed
that you are building/compiling from source.

Singularity packages are available for various Linux distributions, but may not
always be up-to-date with the latest source release version.

For full instructions on installation, including building RPMs,
installing pre-built EPEL packages etc. please check the
[installation section of the admin guide](https://singularity.hpcng.org/admin-docs/master/installation.html).

## Install system dependencies

You must first install development tools and libraries to your host.

On Debian-based systems, including Ubuntu:

```sh
# Ensure repositories are up-to-date
sudo apt-get update
# Install debian packages for dependencies
sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    curl wget git
```

On CentOS/RHEL:

```sh
# Install basic tools for compiling
sudo yum groupinstall -y 'Development Tools'
# Ensure EPEL repository is available
sudo yum install -y epel-release
# Install RPM packages for dependencies
sudo yum install -y \
    libseccomp-devel \
    squashfs-tools \
    cryptsetup \
    wget git
```

## Install Go

Singularity is written in Go, and may require a newer version of Go than is
available in the repositories of your distribution. We recommend installing the
latest version of Go from the [official binaries](https://golang.org/dl/).

First, download the Go tar.gz archive to `/tmp`, then extract the archive to
`/usr/local`.

_**NOTE:** if you are updating Go from a older version, make sure you remove
`/usr/local/go` before reinstalling it._

```sh
export GOVERSION=1.17.3 OS=linux ARCH=amd64  # change this as you need

wget -O /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${GOVERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz
```

Finally, add `/usr/local/go/bin` to the `PATH` environment variable:

```sh
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc
```

## Install golangci-lint

If you will be making changes to the source code, and submitting PRs, you should
install `golangci-lint`, which is the linting tool used in the Singularity
project to ensure code consistency.

Every pull request must pass the `golangci-lint` checks, and these will be run
automatically before attempting to merge the code. If you are modifying
Singularity and contributing your changes to the repository, it's faster to run
these checks locally before uploading your pull request.

In order to download and install the latest version of `golangci-lint`, you can
run:

<!-- markdownlint-disable MD013 -->

```sh
curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $(go env GOPATH)/bin v1.43.0
```

<!-- markdownlint-enable MD013 -->

Add `$(go env GOPATH)` to the `PATH` environment variable:

```sh
echo 'export PATH=$PATH:$(go env GOPATH)/bin' >> ~/.bashrc
source ~/.bashrc
```

## Clone the repo

With the adoption of Go modules you no longer need to clone the Singularity
repository to a specific location.

Clone the repository with `git` in a location of your choice:

```sh
git clone https://github.com/hpcng/singularity.git
cd singularity
```

By default your clone will be on the `master` branch which is where development
of Singularity happens.
To build a specific version of Singularity, check out a
[release tag](https://github.com/hpcng/singularity/tags) before compiling,
for example:

```sh
git checkout v3.8.4
```

## Compiling Singularity

You can configure, build, and install Singularity using the following commands:

```sh
./mconfig
cd ./builddir
make
sudo make install
```

And that's it! Now you can check your Singularity version by running:

```sh
singularity --version
```

The `mconfig` command accepts options that can modify the build and installation
of Singularity. For example, to build in a different folder and to set the
install prefix to a different path:

```sh
./mconfig -b ./buildtree -p /usr/local
```

See the output of `./mconfig -h` for available options.

## Building & Installing from an RPM

On a RHEL / CentOS / Fedora machine you can build a Singularity into an rpm
package, and install it from the rpm. This is useful if you need to install
Singularity across multiple machines, or wish to manage all software via
`yum/dnf`.

To build the rpm, in addition to the
[dependencies](#install-system-dependencies),
install `rpm-build`, `wget`, and `golang`:

```sh
sudo yum install -y rpm-build wget golang
```

The rpm build can use the distribution or EPEL version of Go, even
though as of this writing that version is older than the default
minimum version of Go that Singularity requires.
This is because the rpm applies a source code patch to lower the minimum
required.

To build from a release source tarball do these commands:

<!-- markdownlint-disable MD013 -->

```sh
export VERSION=3.8.4  # this is the singularity version, change as you need

# Fetch the source
wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
# Build the rpm from the source tar.gz
rpmbuild -tb singularity-${VERSION}.tar.gz
# Install Singularity using the resulting rpm
sudo rpm -ivh ~/rpmbuild/RPMS/x86_64/singularity-${VERSION}-1.el7.x86_64.rpm
# (Optionally) Remove the build tree and source to save space
rm -rf ~/rpmbuild singularity-${VERSION}*.tar.gz
```

<!-- markdownlint-enable MD013 -->

Alternatively, to build an RPM from the latest master you can
[clone the repo as detailed above](#clone-the-repo).
Create the build configuration using the `--only-rpm` option of
`mconfig` if you're using the system's too-old golang installation,
to lower the minimum required version.
Then use the `rpm` make target to build Singularity as an rpm package:

<!-- markdownlint-disable MD013 -->

```sh
./mconfig --only-rpm
make -C builddir rpm
sudo rpm -ivh ~/rpmbuild/RPMS/x86_64/singularity-3.8.4*.x86_64.rpm # or whatever version you built
```

<!-- markdownlint-enable MD013 -->

By default, the rpm will be built so that Singularity is installed under
`/usr/local`.

To build an rpm with an alternative install prefix set RPMPREFIX on the make
step, for example:

```sh
make -C builddir rpm RPMPREFIX=/opt/singularity
```

For more information on installing/updating/uninstalling the RPM, check out our
[admin docs](https://singularity.hpcng.org/admin-docs/master/admin_quickstart.html).

## Debian Package

Additional information on how to build a Debian package can be found in [dist/debian/DEBIAN_PACKAGE.md](dist/debian/DEBIAN_PACKAGE.md).
# Creating the Debian package

## Preparation

As long as the debian directory is in the sub-directory you need to link
or copy it in the top directory.

In the top directory do this:

```sh
rm -rf debian
cp -r dist/debian .
```

Make sure all the dependencies are met. See the `INSTALL.md` for this.
Otherwise `debuild` will complain and quit.

## Configuration

To do some configuration for the build, some environment variables can
be used.

Due to the fact, that `debuild` filters out some variables, all the
configuration variables need to be prefixed by `DEB_`

### mconfig

See `mconfig --help` for details about the configuration options.

`export DEB_NOSUID=1`    adds --without-suid

`export DEB_NONETWORK=1` adds --without-network

`export DEB_NOSECCOMP=1` adds --without-seccomp

`export DEB_NOALL=1`     adds all of the above

To select a specific profile for `mconfig`.

__REMINDER:__
to build with seccomp you need to install `libseccomp-dev` package !

For real production environment us this configuration:

```sh
export DEB_SC_PROFILE=release-stripped
```

or if debugging is needed use this.

```sh
export DEB_SC_PROFILE=debug
```

In case a different build directory is needed:

```sh
export DEB_SC_BUILDDIR=builddir
```

### debchange

One way to update the changelog would be that the developer of singularity
update the Debian changelog on every commit. As this is double work, because
of the CHANGELOG.md in the top directory, the changelog is automatically
updated with the version of the source which is currently checked out.
Which means you can easily build Debian packages for all the different tagged
versions of the software. See `INSTALL.md` on how to checkout a specific
version.

Be aware, that `debchange` will complain about a lower version as the top in
the current changelog. Which means you have to cleanup the changelog if needed.
If you did not change anything in the debian directory manually, it might be easiest
to [start from scratch](#Preparation).
Be aware, that the Debian install directory as you see it now, might not be available
in older versions (branches, tags). Make sure you have a clean copy of the debian
directory before you switch to (checkout) an older version.

Usually `debchange` is configured by the environment variables
`DEBFULLNAME` and `EMAIL`. As `debuild` creates a clean environment it
filters out most of the environment variables. To set `DEBFULLNAME` for
the `debchange` command in the makefile, you have to set `DEB_FULLNAME`.
If these variables are not set, `debchange` will try to find appropriate
values from the system configuration. Usually by using the login name
and the domain-name.

```sh
export DEB_FULLNAME="Your Name"
export EMAIL="you@example.org"
```

## Building

As usual for creating a Debian package you can use `dpkg-buildpackage`
or `debuild` which is a kind of wrapper for the first and includes the start
of `lintian`, too.

```sh
dpkg-buildpackage --build=binary --no-sign
lintian --verbose --display-info --show-overrides
```

or all in one

```sh
debuild --build=binary --no-sign --lintian-opts --display-info --show-overrides
```

After successful build the Debian package can be found in the parent directory.

To clean up the temporary files created by `debuild` use the command:

```sh
dh clean
```

To cleanup the copy of the debian directory, make sure you saved your
changes (if any) and remove it.

```sh
rm -rf debian
```

For details on Debian package building see the man-page of `debuild` and
`dpkg-buildpackage` and `lintian`

## Debian Repository

In the current version this is by far not ready for using it in official
Debian Repositories.

This might change in future. I updated the old debian directory to make
it just work, for people needing it.

Any help is welcome to provide a Debian installer which can be used for
building a Debian package,
that can be used in official Debian Repositories.
### Version of Singularity:

What version of Singularity are you using? Run:

```
$ singularity --version

```


<!-- please include command-line output in a code block -->

### Expected behavior

What did you expect to see when you do...?


### Actual behavior

What actually happend? Why was it incorrect?



<!-- if this is a feature request, you can ignore this next part -->

### Steps to reproduce this behavior

How can others reproduce this issue/problem?


### What OS/distro are you running

```
$ cat /etc/os-release


```


### How did you install Singularity

Write here how you installed Singularity. Eg. RPM, source.


## Description of the Pull Request (PR):

Write your description of the PR here. Be sure to include as much background,
and details necessary for the reviewers to understand exactly what this is
fixing or enhancing.


### This fixes or addresses the following GitHub issues:

 - Fixes #


#### Before submitting a PR, make sure you have done the following:

- Read the [Guidelines for Contributing](https://github.com/hpcng/singularity/blob/master/CONTRIBUTING.md), and this PR conforms to the stated requirements.
- Added changes to the [CHANGELOG](https://github.com/hpcng/singularity/blob/master/CHANGELOG.md) if necessary according to the [Contribution Guidelines](https://github.com/hpcng/singularity/blob/master/CONTRIBUTING.md)
- Added tests to validate this PR, linted with `make check`  and tested this PR locally with a `make test`, and `make testall` if possible (see CONTRIBUTING.md).
- Based this PR against the appropriate branch according to the [Contribution Guidelines](https://github.com/hpcng/singularity/blob/master/CONTRIBUTING.md)
- Added myself as a contributor to the [Contributors File](https://github.com/hpcng/singularity/blob/master/CONTRIBUTORS.md)
# Examples

The example bootstrap definition files, each called `Singularity`, are located
in their respectively named folders in this directory. These files can be used
to create new container images on a variety of Linux distributions or become the
basis for customization to build reproducible containers for a specific purpose.
While many of these examples use core mirrors and OS distributions, keep in mind
that you can use a Docker bootstrap to create almost any of them.

## Contributing

If you have a specific scientific (or other) container, we suggest that you
consider [singularity hub](https://singularity-hub.org) to serve it. If you do
not intend to build or use the container, or want to provide a base template,
then you might also want to send a pull request to add it here.

### contrib

If you wish to contribute a definition file that does not fall within one of the
folders here, it should go into [contrib](contrib). In this case, please send a
pull request and contribute it to the examples/contribs directory with the
format being hyphen ('-') delimited of the following format:

1. Base distribution name and version if applicable (e.g. centos7 or
   ubuntu_trusty)
1. Target nomenclature that describes the container (e.g. tensorflow)
1. Any relevant version strings to the application or work-flow
1. Always end in .def

An example of this:

```text
examples/contrib/debian84-tensorflow-0.10.def
```

### base

If your contribution is more appropriate for one of the base or template
distributions, then please make a respective folder in the [examples](.)
directory, and name the definition file `Singularity`.
# Arch for Singularity

This bootstrap spec will generate an arch linux distribution using Singularity
2.3 (current development branch). Note that you can also just bootstrap a Docker
image:

If you want to move forward with the raw, old school, jeans and hard toes
bootstrap, here is what to do. I work on an Ubuntu machine, so I had to use a
Docker Arch Linux image to do this. This first part you should do on your local
machine (if not arch linux) is to use Docker to interactively work in an Arch
Linux image. If you don't want to copy paste the build spec file, you can use
`--volume` to mount a directory from your host to a folder in the image (I would
recommend `/tmp` or similar). Here we run the docker image:

```bash
docker run -it  --privileged pritunl/archlinux bash
```

```bash
pacman -S -y git autoconf libtool automake gcc python make sudo vim \
 arch-install-scripts wget
git clone https://github.com/hpcng/singularity
cd singularity
git checkout -b development
git pull origin development
./autogen.sh
./configure --prefix=/usr/local
```

You can add the [Singularity](Singularity) build spec here, or cd to where it is
if you have mounted a volume.

```bash
cd /tmp
singularity create arch.img
sudo singularity bootstrap arch.img Singularity
```

That should do the trick!
# Singularity example plugin

This directory contains an example CLI plugin for singularity. It demonstrates
how to add a command and flags.

## Building

In order to build the plugin you need a copy of code matching the version of
singularity that you wish to use. You can find the commit matching the
singularity binary by running:

```console
$ singularity version
3.1.1-723.g7998470e7
```

this means this version of singularity is _post_ 3.1.1 (but before the
next version after that one). The suffix .gXXXXXXXXX indicates the exact
commit in github.com/hpcng/singularity used to build this binary
(7998470e7 in this example).

Obtain a copy of the source code by running:

```sh
git clone https://github.com/hpcng/singularity.git
cd singularity
git checkout 7998470e7
```

Still from within that directory, run:

```sh
singularity plugin compile ./examples/plugins/cli-plugin
```

This will produce a file `./examples/plugins/cli-plugin/cli-plugin.sif`.

Currently there's a limitation regarding the location of the plugin code: it
must reside somewhere _inside_ the singularity source code tree.

## Installing

Once you have compiled the plugin into a SIF file, you can install it into the
correct singularity directory using the command:

```sh
singularity plugin install ./examples/plugins/cli-plugin/cli-plugin.sif
```

Singularity will automatically load the plugin code from now on.

## Other commands

You can query the list of installed plugins:

```console
$ singularity plugin list
ENABLED  NAME
    yes  sylabs.io/cli-plugin
```

Disable an installed plugin:

```sh
singularity plugin disable sylabs.io/cli-plugin
```

Enable a disabled plugin:

```sh
singularity plugin enable sylabs.io/cli-plugin
```

Uninstall an installed plugin:

```sh
singularity plugin uninstall sylabs.io/cli-plugin
```

And inspect a SIF file before installing:

```console
$ singularity plugin inspect examples/plugins/cli-plugin/cli-plugin.sif
Name: sylabs.io/cli-plugin
Description: This is a short test CLI plugin for Singularity
Author: Sylabs Team
Version: 0.1.0
```
# Build Singularity

## Summary

This is a build container that generates installable singularity packages for
singularity v3.X.X. The container will output a deb and rpm in the current
directory.

## Known Bugs

Some versions of singularity contain the character 'v', such as v3.0.0. The
container will have to be rebuilt with the following statement modified:

```sh
curl -L -o singularity-${VERSION}.tar.gz
    https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
```

## Usage

```sh
sudo singularity build build-singularity.sif build-singularity.def

./build-singularity.sif {version}

./build-singularity.sif 3.8.0
```
# Singularity SCI-F Apps

The Scientific Filesystem is well suited for Singularity containers to allow you
to build a container that has multiple entrypoints, along with modular environments,
libraries, and executables. Here we will review the basic building and using of a
Singularity container that implements SCIF. For more quick start tutorials, see
the [official documentation for SCIF](https://vsoch.github.io/scif/).

Build your image

```sh
sudo singularity build cowsay.simg Singularity.cowsay 
```

What apps are installed?

```console
$ singularity apps cowsay.simg
cowsay
fortune
lolcat
```

Ask for help for a specific app!

```console
$ singularity help --app fortune cowsay.simg
fortune is the best app
```

Run a particular app

```console
$ singularity run --app fortune cowsay.simg
When I reflect upon the number of disagreeable people who I know who have gone
to a better world, I am moved to lead a different life.
    -- Mark Twain, "Pudd'nhead Wilson's Calendar"
```

Inspect an app

```console
$ singularity inspect --app fortune cowsay.img 
{
    "SCIF_APPNAME": "fortune",
    "SCIF_APPSIZE": "1MB"
}
```
# Bootstrap Self

A self bootstrap means packaging the current operating system that you live on
into an image. Since we assume the root, you don't need to define a `From`. It
looks like this:

## Options

```text
Bootstrap: self
```

If you really wanted to specify some root, you could do this:

```text
Bootstrap: self
From: /
```

And we highly recommend that you exclude paths that you don't want added to the
tar. For example, Docker stores a lot of data in `/var`, so I chose to exclude
that, along with some of the applications in `/opt`:

```text
Bootstrap: self
Exclude: /var/lib/docker /home/vanessa /opt/*
```

## Build Example

so we could do the following with the specification build file in this folder:

```sh
singularity create --size 8000 container.img
sudo singularity bootstrap container.img Singularity
```
# MAKEIT

## Goals

- To generate **native Makefiles** for the system *makeit* is running on.
  To accomplish that, we transform a set of Makefile fragments and
  module config files into a non-recursive Makefile that uses features
  available to all reasonable versions of Make (GNU, BSD, SVR4, etc).

- If **makeit** can be called a build system, it should stay so small and
  platform independent that it could be included in each project that
  it helps to build.

- To include/install and setup *makeit* for your project take a look at the
  INSTALL.md file.

## Module (\*.mconf) Keywords

- **name** : name of the module, just a handle
- **prog** : name of a program to link
- **lib** : name of a library to create, without the **lib** prefix
- **data** : name of a data file (symbols, pictures, text, etc.) to embed
- **asrc** : list of (.S) assembly source files
- **csrc** : list of C source files to build
- **win_asrc** : windows only list of C source files to build
- **win_csrc** : windows only list of C source files to build
- **unix_asrc** : unix only list of C source files to build
- **unix_csrc** : unix only list of C source files to build
- **depends** : list of module **name**'s that a prog or a lib depends on
- **cflags** : list of CFLAGS to add for this module
- **ldflags** : list of LDFLAGS to add for this module
- **extralibs** : list of extra libs needed by the program (e.g., -lgcc)
- **cleanfiles** : list of extra files to remove when *make clean* is called

## Implementation

- POSIX portable tools mainly awk and sh with system commands
# End-to-End Testing

This package contains the end-to-end tests for `singularity`.

## Contributing

For this example, we're going to use a topic of `env` or
`environment variable tests`.

- Add your topic as a runtime-hook in `suite.go`.

```go
// RunE2ETests by functionality
t.Run("BUILD", imgbuild.RunE2ETests)
t.Run("ACTIONS", actions.RunE2ETests)
t.Run("ENV", env.RunE2ETests)
```

- Create a directory for your topic.

```sh
mkdir -p e2e/env
```

- Create a source file to include your topic's tests.

```sh
touch e2e/env/env.go
```

- Optionally create a source file to include helpers for your topic's test.

```sh
touch e2e/env/env_utils.go
```

- Add a package declaration to your topic's test file that matches what you put
  in `suite.go`

```go
package env
```

- Add a variable to store the testing settings in your topic's test file.

```go
import (
        "github.com/kelseyhightower/envconfig"
)

type testingEnv struct {
	// base env for running tests
	CmdPath     string `split_words:"true"`
	TestDir     string `split_words:"true"`
	RunDisabled bool   `default:"false"`
}

var testenv testingEnv
```

- Add a entry-point to your topic's test file that matches what you put in
  `suite.go`

```go
//RunE2ETests is the main func to trigger the test suite
func RunE2ETests(t *testing.T) {
	err := envconfig.Process("E2E", &testenv)
	if err != nil {
		t.Fatal(err.Error())
	}
}
```

- Create a test in your topic's test file as you normally would in `go`.

```go
func TestEnv(t *Testing.T) {
	...
}
```

- Run your test from your entry-point function using a `go` sub-test.

```go
//RunE2ETests is the main func to trigger the test suite
func RunE2ETests(t *testing.T) {
	err := envconfig.Process("E2E", &testenv)
        if err != nil {
        	t.Fatal(err.Error())
        }
        
	// Add tests
	t.Run("TestEnv", TestEnv)
}
```

- Example of what your topic's test file might look like:

```go
package env 

import (
	"github.com/kelseyhightower/envconfig"
)

type testingEnv struct {
	// base env for running tests
	CmdPath     string `split_words:"true"`
	TestDir     string `split_words:"true"`
	RunDisabled bool   `default:"false"`
}

var testenv testingEnv

func TestEnv(t *testing.T) {
	...
}

//RunE2ETests is the main func to trigger the test suite
func RunE2ETests(t *testing.T) {
	err := envconfig.Process("E2E", &testenv)
	if err != nil {
		t.Fatal(err.Error())
	}

	t.Run("TestEnv", TestEnv)
}
```

## Running

Test your topic using the `e2e` target in the `Makefile`. To avoid skipping
these tests (default), make sure you set the environment variable
`SINGULARITY_E2E` to `1`.

```sh
SINGULARITY_E2E=1 make -C builddir e2e-test
```

- Verify that your test was run by modifying the `Makefile` to add a verbose
  flag (`go test -v`) and re-running the previous `make` step.
# Testing `singularity help` content

This package contains the end-to-end tests for `singularity help`.

## Contributing new help tests

For this example, we're going to create a new test for
`singularity help inspect`.

- Add the help text to the `testdata/help` directory.

```sh
singularity help inspect > e2e/testdata/help/help-inspect.txt
```

- Add the help command to the `helpContentTests` struct in `help.go`

```go
var helpContentTests = []struct {
        cmds []string
}{
	...
	// singularity inspect
	{[]string{"help", "inspect"}},
	...
}	
```

## Updating existing help tests

For this example, we're going to update an existing test for
`singularity help inspect`.

- When a help test fails, we need to check why it failed.

  - Was the failure a result of an unintended change? If so, we open an issue.
  - Was the failure a result of an intended change? If so, we update the help
    text.

- Update the help text in the `testdata/help` directory.

```sh
singularity help inspect > e2e/testdata/help/help-inspect.txt
```

## Running the help tests

To verify this test, modify the `Makefile` to add both a verbose flag and a
filter flag (`go test -v -r helpContentTests`) and then run the tests.

```sh
SINGULARITY_E2E=1 make -C builddir e2e-test
```
