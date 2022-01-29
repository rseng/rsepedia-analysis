# Cromwell Change Log

## 75 Release Notes

### New `AwaitingCloudQuota` backend status

For Cloud Life Sciences v2beta only.

When a user's GCP project reaches a quota limit, Cromwell continues to submit jobs and Life Sciences acknowledges them as created even if the physical VM cannot yet start. Cromwell now detects this condition in the backend and reports `AwaitingCloudQuota`.

The status is informational and does not require any action. Users wishing to maximize throughput can use `AwaitingCloudQuota` as an indication they should check quota in Cloud Console and request a quota increase from GCP.

`AwaitingCloudQuota` will appear between the `Initializing` and `Running` backend statuses, and will be skipped if not applicable.

Now:

| Status in metadata |Quota normal| Quota delay          | Status meaning                                    |
|--------------------|----|----------------------|---------------------------------------------------|
| `executionStatus`    |`Running`| `Running`            | Job state Cromwell is requesting from the backend |
| `backendStatus`      |`Running`| `AwaitingCloudQuota` | Job state reported by backend                          |

Previously:

| Status in metadata |Quota normal|Quota delay| Status meaning                                            |
|--------------------|----|----|-----------------------------------------------------------|
| `executionStatus`    |`Running`|`Running`| Job state Cromwell is requesting from the backend |
| `backendStatus`      |`Running`|`Running`| Job state reported by backend |

## 73 Release Notes

### Workflow Restart Performance Improvements

Cromwell now allows for improved performance restarting large workflows through the use of a separate rate limiter for restart checks than the rate limiter used for starting new jobs.
The restart check rate limiter is pre-configured in Cromwell's bundled [reference.conf](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf); see the `job-restart-check-rate-control` stanza in that file for explanations of the various parameters if adjustments are desired.

## 71 Release Notes

### Bug Fixes

* Fixed an issue handling data in Google Cloud Storage buckets with requester pays enabled that could sometimes cause I/O to fail.

## 70 Release Notes

### CWL security fix [#6510](https://github.com/broadinstitute/cromwell/pull/6510)

Fixed an issue that could allow submission of an untrusted CWL file to initiate remote code execution. The vector was improper deserialization of the YAML source file.

CWL execution is enabled by default unless a `CWL` [stanza](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf#L460-L482) is present in the configuration that specifies `enabled: false`. Cromwell instances with CWL disabled were not affected. Consequently, users who wish to mitigate the vulnerability without upgrading Cromwell may do so via this config change.

- Thank you to [Bruno P. Kinoshita](https://github.com/kinow) who first found the issue in a different CWL project ([CVE-2021-41110](https://github.com/common-workflow-language/cwlviewer/security/advisories/GHSA-7g7j-f5g3-fqp7)) and [Michael R. Crusoe](https://github.com/mr-c) who suggested we investigate ours.

## 68 Release Notes

### Virtual Private Cloud

Previous Cromwell versions allowed PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the
information to Google Cloud project labels.

Cromwell now allows PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the network and
subnetwork name directly inside the `virtual-private-cloud` backend configuration. More info
[here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

## 67 Release Notes

### Configuration updates for improved scaling

Some configuration changes were introduced in Cromwell 67 to support improved scaling. See Cromwell's `reference.conf` for details on new parameters.

* I/O throttling moved from `io` to its own `io.throttle` stanza; config updates may be required if these values are currently being overridden in local deployments.

* The default `system.job-rate-control` has been changed from 50 per second to 20 per 10 seconds.

* New configuration parameters have been introduced for values which were previously hardcoded constants:
  * `system.file-hash-batch-size`, value updated from `100` to `50`.
  * `io.gcs.max-batch-size`, value stays the same at `100`.
  * `io.gcs.max-batch-duration`, value stays the same at `5 seconds`.

* New configuration parameters which should not require updating:
  * `io.command-backpressure-staleness`
  * `io.backpressure-extension-log-threshold`
  * `load-control.io-normal-window-minimum`
  * `load-control.io-normal-window-maximum`

* `io.nio.parallelism` was previously misspelled in `reference.conf` but not in Cromwell's configuration reading code. Only correct spellings of this configuration key had or will have effect.

## 66 Release Notes

### Google Artifact Registry Support
Cromwell now supports call caching when using Docker images hosted on
[Google Artifact Registry](https://cloud.google.com/artifact-registry).

### Google Image Repository Hashing Updates
The previously documented `docker.hash-lookup.gcr` configuration has been renamed to `docker.hash-lookup.google` and
now applies to both Google Container Registry (GCR) and Google Artifact Registry (GAR) repositories.
Support for the `docker.hash-lookup.gcr-api-queries-per-100-seconds` configuration key has been formally discontinued
and a bug preventing correct handling of `docker.hash-lookup...throttle` configuration has been fixed.
Please see Cromwell's bundled
[`reference.conf`](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf)
for more details.

## 65 Release Notes

* An additional set of metrics relating to metadata age were added.

### AMD Rome support on PAPI v2
On the PAPI v2 backends "AMD Rome" is now supported as a CPU platform. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#cpuplatform).

## 64 Release Notes

### Intel Cascade Lake support on PAPI v2

On the PAPI v2 backends "Intel Cascade Lake" is now supported as a CPU platform. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#cpuplatform).

## 63 Release Notes

### Removed refresh token authentication mode

Google Pipelines API v1 supported authentication with refresh tokens, while v2 of the API does not.

Now that v1 has been discontinued and shut down, this version of Cromwell removes support for refresh tokens.

## 62 Release Notes

### Downloading Access URLs

Added experimental support to download data during Google [Cloud Life Sciences](https://cloud.google.com/life-sciences)
jobs using [DRS
AccessURLs](https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.1.0/docs/#_accessurl).

## 61 Release Notes

### No labels update for Archived workflows

If **- and ONLY if -** you have metadata archiving turned on, then for a workflow whose metadata has been archived by Cromwell 
according to the lifecycle policy, Cromwell will no longer add new labels or update existing labels for this workflow 
coming through PATCH `/labels` endpoint.

## 60 Release Notes

### Java 11

As of this version, a distribution of Java 11 is required to run Cromwell. Cromwell is developed, tested, and
containerized using [AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/).

### Hybrid metadata storage ("carboniting") removed

Carboniting functionality has been removed from Cromwell. 
There will be no effect for customers who store metadata permanently in the relational database (most common),
and there will also be no effect for customers who use the in-memory database.

Breaking change only for customers who explicitly enabled `carbonite-metadata-service` in their configuration to split
metadata storage between a relational database and Google Cloud Storage. If you had previously enabled carboniting and 
deletion, any workflows marked as `ArchivedAndPurged` in your database will no longer be accessible via the Cromwell metadata API.

## 59 Release Notes

### Bug Fixes

* Fixed a pair of bugs that could cause workflows to fail unexpectedly with the errors "413 Request Entity Too Large"
  and "java.net.SocketTimeoutException: Read timed out" when accessing Google Cloud Storage.

## 58 Release Notes

Internal CI-related changes only.

## 57 Release Notes

### Breaking configuration change to reference disk support on PAPI v2

Beginning with Cromwell 57, reference disk manifests are now specified completely within Cromwell configuration
rather than through a level of indirection to a manifest file stored in GCS. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/backends/Google#reference-disk-support).

## 56 Release Notes

### Retry with More Memory as workflow option

The experimental memory retry feature gains per-workflow customization and includes breaking changes:
* The per-backend configuration key `<backend>.config.memory-retry.error-keys` has been removed and replaced 
with global key `system.memory-retry-error-keys`
* The per-backend configuration key `<backend>.config.memory-retry.multiplier` has been replaced with **workflow option** 
`memory_retry_multiplier`

More details can be found [here](https://cromwell.readthedocs.io/en/develop/wf_options/Overview.md#retry-with-more-memory-multiplier).

### Bug Fixes

* Fixed a bug that caused Cromwell to mark workflows as failed after a single `500`, `503`, or `504` error from Google Cloud Storage.
  * Cromwell will now retry these errors as designed.
  * The default retry count is `5` and may be customized with `system.io.number-of-attempts`. 

## 55 Release Notes

### Apple Silicon support statement

Users with access to the new Mac hardware should review [important information provided here](https://cromwell.readthedocs.io/en/stable/Releases).

### Bug Fixes

* Fixed a bug that prevented `read_json()` from working with arrays and primitives. The function now works as expected for all valid JSON data inputs. 
More information on JSON Type to WDL Type conversion can be found [here](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#mixed-read_jsonstringfile).

* Now retries HTTP 408 responses as well as HTTP 429 responses during DOS/DRS resolution requests.

* Fixed a bug that prevented the call caching diff endpoint from working with scatters in workflows with archived metadata.

### New Features

#### Reference disk support on PAPI v2

Cromwell now offers support for the use of reference disks on the PAPI v2 backend as an alternative to localizing
reference inputs. More details [here](https://cromwell.readthedocs.io/en/develop/backends/Google#reference-disk-support).

#### Docker image cache support on PAPI v2 lifesciences beta

Cromwell now offers support for the use of Docker image caches on the PAPI v2 lifesciences beta backend. More details [here](https://cromwell.readthedocs.io/en/develop/backends/Google#docker-image-cache-support).

#### Preemptible Recovery via Checkpointing

* Cromwell can now help tasks recover from preemption by allowing them to specify a 'checkpoint' file which will be restored
to the worker VM on the next attempt if the task is interrupted. More details [here](https://cromwell.readthedocs.io/en/develop/optimizations/CheckpointFiles)

## 54 Release Notes

### Bug Fixes

* Fixed a bug that prevented `write_json()` from working with arrays and primitives. The function now works as expected for `Boolean`, `String`, `Integer`, `Float`,
 `Pair[_, _]`, `Object`, `Map[_, _]` and `Array[_]` (including array of objects) type inputs. More information on WDL Type to JSON Type 
 conversion can be found [here](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#mixed-read_jsonstringfile).

### Spark backend support removal

Spark backend was not widely used and it was decided to remove it from the codebase in order to narrow the scope of Cromwell code. 

### Improved DRS Localizer logging

Error logging while localizing a DRS URI should now be more clear especially when there is a Requester Pays bucket involved.

### Per-backend hog factors
Cromwell now allows overriding system-level log factors on back-end level. First, Cromwell will try to use hog-factor 
defined in the backend config, and if it is not defined, it will default to using system-wide hog factor.
```conf
backend {
  providers {
    PAPIv2 {
      config {
        hog-factor: 2
      }
    }
  }
}
```
For more information about hog factors please see [this page](https://cromwell.readthedocs.io/en/develop/cromwell_features/HogFactors/).

### `martha_v2` Support Removed

Cromwell now only supports resolving DOS or DRS URIs through [Martha](https://github.com/broadinstitute/martha)'s most
recent metadata endpoint `martha_v3`, dropping support for Martha's previous metadata endpoint `martha_v2`. To switch to
the new version of Martha's metadata endpoint, update the `martha.url` found in the [filesystems
config](https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#overview) to point to `/martha_v3`. More
information on Martha's `martha_v3` request and response schema can be found
[here](https://github.com/broadinstitute/martha#martha-v3).

### DOS/DRS `localization_optional` Support

When running on a backend that supports `localization_optional: true` any DOS or DRS `File` values in the generated
command line will be substituted with the `gsUri` returned from Martha's `martha_v3` endpoint. More information on
`localization_optional` can be found [here](https://cromwell.readthedocs.io/en/stable/optimizations/FileLocalization/).

### DOS/DRS metadata retrieval retried by default

Attempts to retrieve DOS/DRS metadata from Martha will be retried by default. More information can be found
[here](https://cromwell.readthedocs.io/en/stable/optimizations/FileLocalization/).

## 53 Release Notes

### Martha v3 Support

Cromwell now supports resolving DRS URIs through Martha v3 (in addition to Martha v2). To switch to the new version of Martha, update the `martha.url` found in the [filesystems config](https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#overview) to
point to `/martha_v3`. More information on Martha v3 request and response schema can be found [here](https://github.com/broadinstitute/martha#martha-v3).

### Support for custom entrypoints on Docker images

Cromwell can now support docker images which have custom entrypoints in the PAPIv2 alpha and beta backends.

### Alpha support for WDL optional outputs on PAPI v2

* Alpha support for WDL optional output files on the PAPI v2 backend has been added, please see the
[documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#alpha-support-for-wdl-optional-outputs-on-papi-v2)
for known limitations.

### Monitoring Image Script

* Cromwell now supports an optional `monitoring_image_script` workflow option in addition to the existing
`monitoring_script` and `monitoring_image` options. For more information see the [Google Pipelines API Workflow Options
documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#google-pipelines-api-workflow-options).

## 52 Release Notes

### Documentation

Information on how to properly use the Singularity cache with Cromwell is now
provided in the [Cromwell Singularity documentation](
https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#singularity).

### Google library upgrade [(#5565)](https://github.com/broadinstitute/cromwell/pull/5565)

All previous versions of Cromwell shipped with Google Cloud Storage (GCS) libraries that are now deprecated and will [stop working in August 2020](https://developers.googleblog.com/2018/03/discontinuing-support-for-json-rpc-and.html). This release adopts updated libraries to ensure uninterrupted operation. The only user action required is upgrading Cromwell.   

### Bug fixes

* Fixed a bug that required Cromwell to be restarted in order to pick up DNS changes.
    * By default, the JVM caches DNS records with a TTL of infinity.
    * Cromwell now configures its JVM with a 3-minute TTL. This value can be customized by setting `system.dns-cache-ttl`.  
* Clarified an error message that Cromwell emits when the compute backend terminates a job of its own volition (as opposed to termination in response to an abort request from Cromwell)
    * Previously, the error read `The job was aborted from outside Cromwell`
    * The new error reads `The compute backend terminated the job. If this termination is unexpected, examine likely causes such as preemption, running out of disk or memory on the compute instance, or exceeding the backend's maximum job duration.` 

## 51 Release Notes

### Changes and Warnings

The configuration format for call cache blacklisting has been updated, please see the [call caching documentation](
https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching) for details.

### Bug fixes

* Fixed a bug where the `size(...)` function did not work correctly on files 
  from a shared filesystem if `size(...)` was called in the input section on a 
  relative path.
+ Fixed a bug where the `use_relative_output_paths` option would not preserve intermediate folders.

### New functionality

#### Call caching blacklisting improvements

Cromwell previously supported blacklisting GCS buckets containing cache hits which could not be copied for permissions 
reasons. Cromwell now adds support for blacklisting individual cache hits which could not be copied for any reason,
as well as grouping blacklist caches according to a workflow option key. More information available in the [
call caching documentation]( https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching). 

#### new xxh64 and fingerprint strategies for call caching

Existing call cache strategies `path` and `path+modtime` don't work when using docker on shared filesystems 
(SFS backend, i.e. not in cloud storage). The `file` (md5sum) strategy works, but uses a lot of resources.
Two faster strategies have been added for this use case: `xxh64` and 
`fingerprint`. `xxh64` is a lightweight hashing algorithm, `fingerprint` is a strategy designed to be very 
lightweight. Read more about it in the [call caching documentation](
https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching).

## 50 Release Notes

### Changes and Warnings

#### Metadata Archival Config Change

**Note:** Unless you have already opted-in to GCS-archival of metadata during its development, this change will not affect you.
Cromwell's metadata archival configuration has changed in a backwards incompatible way to increase consistency,
please see
[the updated documentation](https://cromwell.readthedocs.io/en/stable/Configuring#hybrid-metadata-storage-classic-carbonite) for details.

## 49 Release Notes

### Changes and Warnings

#### Job store database refactoring

The primary keys of Cromwell's job store tables have been refactored to use a `BIGINT` datatype in place of the previous
`INT` datatype. Cromwell will not be usable during the time the Liquibase migration for this refactor is running.
In the Google Cloud SQL with SSD environment this migration runs at a rate of approximately 40,000 `JOB_STORE_SIMPLETON_ENTRY`
rows per second. In deployments with millions or billions of `JOB_STORE_SIMPLETON_ENTRY` rows the migration may require
a significant amount of downtime so please plan accordingly. The following SQL could be used to estimate the number of
rows in this table:

```
SELECT table_rows FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'cromwell' AND table_name = 'JOB_STORE_SIMPLETON_ENTRY';
```

#### Execution Directory Layout (cache copies)

When an attempt to copy a cache result is made, you'll now see a `cacheCopy` directory in the call root directory. 
This prevents them clashing with the files staged to the same directory for attempt 1 if the cache copy fails (see also: Bug Fixes).

The directory layout used to be:

```
[...]/callRoot/
  - script [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - stdout [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - output.file [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - attempt-2/ [if attempt 1 fails]
    - script
    - stdout
    - output.file
```

but is now:

```
[...]/callRoot/
  - cacheCopy/
    - script
    - stdout
    - output.file
  - script [for attempt 1 if the cache copy fails]
  - stdout [for attempt 1 if the cache copy fails]
  - output.file [for attempt 1 if the cache copy fails]
  - attempt-2/ [if attempt 1 fails]
    - script
    - stdout
    - output.file
```

### New Functionality

#### Disable call-caching for tasks

It is now possible to indicate in a workflow that a task should not be call-cached. See details 
[here](https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks).

#### Delete Intermediate Outputs on PapiV2

* **Experimental:** When a new workflow option `delete_intermediate_output_files` is submitted with the workflow,
intermediate `File` objects will be deleted when the workflow completes. See the [Google Pipelines API Workflow Options
documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#google-pipelines-api-workflow-options)
for more information.

#### Metadata Archival Support

Cromwell 49 now offers the option to archive metadata to GCS and remove the equivalent metadata from relational
database storage. Please see 
[the documentation](https://cromwell.readthedocs.io/en/stable/Configuring#hybrid-metadata-storage-classic-carbonite) for more details. 

#### Adding support for Google Cloud Life Sciences v2beta
Cromwell now supports running workflows using Google Cloud Life Sciences v2beta API in addition to Google Cloud Genomics v2alpha1. 
More information about migration to the new API from v2alpha1 
[here](https://cromwell.readthedocs.io/en/stable/backends/Google#migration-from-google-cloud-genomics-v2alpha1-to-google-cloud-life-sciences-v2beta). 
* **Note** Google Cloud Life Sciences is the new name for newer versions of Google Cloud Genomics.
* **Note** Support for Google Cloud Genomics v2alpha1 will be removed in a future version of Cromwell. Advance notice will be provided.

### New Docs

#### Installation methods

Links to the conda package and docker container are now available in 
[the install documentation](https://cromwell.readthedocs.io/en/stable/Getting/).


### Bug Fixes

+ Fix a bug where zip files with directories could not be imported. 
  For example a zip with `a.wdl` and `b.wdl` could be imported but one with `sub_workflows/a.wdl` 
  and `imports/b.wdl` could not.
+ Fix a bug which sometimes allowed execution scripts copied by a failed cache-copy to be run instead
  of the attempt-1 script for a live job execution. 
  
## 48 Release Notes

### Womtool Graph for WDL 1.0

The `womtool graph` command now supports WDL 1.0 workflows. 
* **Note:** Generated graphs - including in WDL draft 2 - may look slightly different than they did in version 47.

### Documentation

+ Documented the use of a HSQLDB file-based database so users can try call-caching without needing a database server.
  Please checkout [the database documentation](https://cromwell.readthedocs.io/en/stable/Configuring#database).

## 47 Release Notes

### Retry with more memory on Papiv2 [(#5180)](https://github.com/broadinstitute/cromwell/pull/5180)

Cromwell now allows user defined retries. With `memory-retry` config you can specify an array of strings which when encountered in the `stderr` 
file by Cromwell, allows the task to be retried with multiplier factor mentioned in the config. More information [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### GCS Parallel Composite Upload Support

Cromwell 47 now supports GCS parallel composite uploads which can greatly improve delocalization performance.
This feature is turned off by default, it can be turned on by either a backend-level configuration setting or
on a per-workflow basis with workflow options. More details [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### Papi V2 Localization Using GCR [(#5200)](https://github.com/broadinstitute/cromwell/pull/5200)

The Docker image for the Google Cloud SDK was previously only [published on Docker
Hub](https://hub.docker.com/r/google/cloud-sdk). Now that the image is [publicly hosted in
GCR](http://gcr.io/google.com/cloudsdktool/cloud-sdk), Papi V2 jobs will localize inputs and delocalize outputs using
the GCR image.

## 46 Release Notes

### Nvidia GPU Driver Update

The default driver for Nvidia GPU's on Google Cloud has been updated from `390` to `418.87.00`.  A user may override this option at anytime by providing the `nvidiaDriverVersion` runtime attribute.  See the [Runtime Attribute description for GPUs](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#runtime-attribute-descriptions) for detailed information.

### Enhanced "error code 10" handling in PAPIv2

On Google Pipelines API v2, a worker VM that is preempted may emit a generic error message like
```
PAPI error code 10. The assigned worker has failed to complete the operation
```
instead of a preemption-specific message like
```
PAPI error code 14. Task was preempted for the 2nd time.
```
Cromwell 44 introduced special handling that detects both preemption indicators and re-runs the job consistent with the `preemptible` setting.

Cromwell 46 enhances this handling in response to user reports of possible continued issues.   

## 45 Release Notes

### Improved input and output transfer performance on PAPI v2

Cromwell now requires only a single PAPI "action" each for the entire localization or delocalization process, rather than two per file or directory.
This greatly increases execution speed for jobs with large numbers of input or output files.
In testing, total execution time for a call with 800 inputs improved from more than 70 minutes to less than 20 minutes.

### List dependencies flag in Womtool Command Line [(#5098)](https://github.com/broadinstitute/cromwell/pull/5098)

Womtool now outputs the list of files referenced in import statements using `-l` flag for `validate` command.
More info [here](https://cromwell.readthedocs.io/en/stable/WOMtool/)

### BCS backend new Features support

#### New docker registry
Alibaba Cloud Container Registry is now supported for the `docker` runtime attribute, and the previous `dockerTag` 
runtime attribute continues to be available for Alibaba Cloud OSS Registry.
#### Call caching
Cromwell now supports Call caching when using the BCS backend.
#### Workflow output glob
Globs can be used to define outputs for BCS backend.
#### NAS mount
Alibaba Cloud NAS is now supported for the `mounts` runtime attribute.

### Call Caching Failure Messages [(#5095)](https://github.com/broadinstitute/cromwell/pull/5095)

Call cache failures are no longer sent to the workflow metadata. Instead a limited number of call cache failure messages
will be sent to the workflow log. See [the Cromwell call caching
documentation](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching/) for more information on call
cache failure logging.

## 44 Release Notes

### Improved PAPI v2 Preemptible VM Support

In some cases PAPI v2 will report the preemption of a VM in a way that differs from PAPI v1. This novel means of reporting
preemption was not recognized by Cromwell's PAPI v2 backend and would result in preemptions being miscategorized as call failures.
Cromwell's PAPI v2 backend will now handle this type of preemption.

## 43 Release Notes

### Virtual Private Cloud with Subnetworks

Cromwell now allows PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the subnetwork key 
`subnetwork-label-key` inside `virtual-private-cloud` in backend configuration. More info [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### Call caching database refactoring

Cromwell's `CALL_CACHING_HASH_ENTRY` primary key has been refactored to use a `BIGINT` datatype in place of the previous
`INT` datatype. Cromwell will not be usable during the time the Liquibase migration for this refactor is running.
In the Google Cloud SQL with SSD environment this migration runs at a rate of approximately 100,000 `CALL_CACHING_HASH_ENTRY`
rows per second. In deployments with millions or billions of `CALL_CACHING_HASH_ENTRY` rows the migration may require  
a significant amount of downtime so please plan accordingly. The following SQL could be used to estimate the number of
rows in this table:

```
select max(CALL_CACHING_HASH_ENTRY_ID) from CALL_CACHING_HASH_ENTRY
```

### Stackdriver Instrumentation

Cromwell now supports sending metrics to [Google's Stackdriver API](https://cloud.google.com/monitoring/api/v3/). 
Learn more on how to configure [here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

### BigQuery in PAPI

Cromwell now allows a user to specify BigQuery jobs when using the PAPIv2 backend

### Configuration Changes

#### StatsD Instrumentation

There is a small change in StatsD's configuration path. Originally, the path to the config was `services.Instrumentation.config.statsd`
which now has been updated to `services.Instrumentation.config`. More info on its configuration can be found
[here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

#### cached-copy

A new experimental feature, the `cached-copy` localization strategy is available for the shared filesystem. 
More information can be found in the [documentation on localization](https://cromwell.readthedocs.io/en/stable/backends/HPC).

#### Yaml node limits

Yaml parsing now checks for cycles, and limits the maximum number of parsed nodes to a configurable value. It also
limits the nesting depth of sequences and mappings. See [the documentation on configuring
YAML](https://cromwell.readthedocs.io/en/stable/Configuring/#yaml) for more information.

### API Changes

#### Workflow Metadata

* It is now possible to use `includeKey` and `excludeKey` at the same time. If so, the metadata key must match the `includeKey` **and not** match the `excludeKey` to be included.
* It is now possible to use "`calls`" as one of your `excludeKey`s, to request that only workflow metadata gets returned.

### PostgreSQL support

Cromwell now supports PostgreSQL (version 9.6 or higher, with the Large Object
extension installed) as a database backend.
See [here](https://cromwell.readthedocs.io/en/stable/Configuring/#database) for
instructions for configuring the database connection.

## 42 Release Notes

### Womtool endpoint

The `/describe` endpoint now differentiates between an invalid workflow and a valid workflow with invalid inputs.

Specifically, the new `validWorkflow` key indicates whether the workflow file is valid by itself. If inputs are provided, they are not considered when calculating this field; if inputs are not provided, the value is identical to `valid`.

### Configuration Changes

 *  Virtual private networks can now be configured. See the section below for details.
 
#### Batch Request Timeouts

The timeout on Cromwell's requests to PAPIv2 can now be configured. See the sample PAPIv2.conf for more documentation:

```conf
backend {
  providers {
    PAPIv2 {
      config { 
        batch-requests {
          timeouts {
            read = 10 seconds
            connect = 10 seconds
          }
        }
      }
    }
  }
}
```

### Virtual Private Networks

Cromwell now allows PAPIV2 jobs to run on a private network by adding the network name inside `virtual-private-cloud` in backend configuration.
More info [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### AWS Backend

Now includes background job status polling to hopefully reduce the incidence of 'HTTP 429' errors for large workflows.

## 41 Release Notes

### Workflow Options

* It is now possible to supply custom `google-labels` in [workflow options](https://cromwell.readthedocs.io/en/stable/wf_options/Google/).

### AWS backend

It is now possible to use WDL disk attributes with the following formats on AWS.
```
disks: "local-disk 20 SSD"
```
```
disks: "/some/mnt 20 SSD"
```
Because Cromwell's AWS backend auto-sizes disks, the size specification is simply discarded.

### Time Formatting

In previous versions of Cromwell, times were converted to strings using
[the default Java formatter](https://docs.oracle.com/javase/8/docs/api/java/time/OffsetDateTime.html#toString--) which
generates a variety of ISO-8601 formats. String conversions also retained whatever server time zone generated that
specific time instance.

Going forward, times stored in Cromwell metadata, and later returned via the HTTP endpoint, are now converted to UTC
then formatted with exactly three digits of milliseconds.

For example:
- `2017-01-19T12:34:56-04:00` will now be formatted as
- `2017-01-19T16:34:56.000Z`

This change only affects newly formatted dates. Older dates already formatted and stored by previous versions of
Cromwell will not be updated however they will still return a
[valid ISO-8601 format](https://en.wikipedia.org/wiki/ISO_8601). The older format may be in various non-UTC time zones,
and may or may not include microseconds or even nanoseconds, for example `2017-01-19T12:34:56.123456789-04:00`.

### Config Changes

#### Heartbeat failure shutdown

When a Cromwell instance is unable to write heartbeats for some period of time it will automatically shut down. For more
information see the docs on [configuring Workflow Hearbeats](https://cromwell.readthedocs.io/en/stable/Configuring/).

NOTE: In the remote chance that the `system.workflow-heartbeats.ttl` has been configured to be less than `5 minutes`
then the new configuration value `system.workflow-heartbeats.write-failure-shutdown-duration` must also be explicitly
set less than the `ttl`.

#### nVidia Driver Attribute Change

The runtime attribute `nvidia-driver-version` was previously allowed only as a default runtime attribute in configuration.
Because WDL does not allow attribute names to contain `-` characters, this has been changed to `nvidiaDriverVersion`.
This field is now accepted within WDL files as well as within the configuration file.

#### Logging long running jobs

All backends can now emit slow job warnings after a configurable time running. 
NB This example shows how to configure this setting for the PAPIv2 backend:
```conf
# Emit a warning if jobs last longer than this amount of time. This might indicate that something got stuck.
backend {
  providers {
    PAPIv2 {
      config { 
        slow-job-warning-time: 24 hours
      }
    }
  }
}
```

### Runtime Attributes

#### GPU Attributes

* The `gpuType` attribute is no longer validated against a whitelist at workflow submission time. Instead, validation now happens at runtime. This allows any valid accelerator to be used.
* The `nvidiaDriverVersion` attribute is now available in WDL `runtime` sections. The default continues to be `390.46` which applies if and only if GPUs are being used.
* A default `gpuType` ("nvidia-tesla-k80") will now be applied if `gpuCount` is specified but `gpuType` is not.
* Similarly, a default `gpuCount` (1) will be applied if `gpuType` is specified but `cpuCount` is not. 

### Bug fixes

#### Better validation of workflow heartbeats

An error will be thrown on startup when the `system.workflow-heartbeats.heartbeat-interval` is not less than the
`system.workflow-heartbeats.ttl`.


## 40 Release Notes

### Config Changes

#### Cromwell ID in instrumentation path

When set, the configuration value of `system.cromwell_id` will be prepended to StatsD metrics. More info [here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

#### HealthMonitor Configuration

The HealthMonitor configuration has been refactored to provide a simpler interface:
* You no longer need to specify a monitor class in your `cromwell.conf` as this will now be inherited from the `reference.conf` value.
* You can now opt-in and opt-out of any combination of status monitors.
* The PAPI backends to monitor can now be listed in a single field.

##### Upgrading

You are no longer tied to the previous preset combinations of health checks. However if you just want to carry forward
the exact same set of health checks, you can use one of the following standard recipes:

###### From default, or `NoopHealthMonitorActor`:
If you're currently using the (default) NoopHealthMonitorActor, no action is required.

###### From `StandardHealthMonitorServiceActor`:
If you're currently using the `StandardHealthMonitorServiceActor`, replace this stanza:
```
services {
    HealthMonitor {
        class = "cromwell.services.healthmonitor.impl.standard.StandardHealthMonitorServiceActor"
    }
}
``` 
With this one:
```
services {
    HealthMonitor {
        config {
            check-dockerhub: true
            check-engine-database: true
        }
    }
}
``` 
###### From `WorkbenchHealthMonitorServiceActor`:
Replace this stanza:
```
services {
    HealthMonitor {
        class = "cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor"

        config {
            papi-backend-name = PAPIv1
            papi-v2-backend-name = PAPIv2

            google-auth-name = service-account
            gcs-bucket-to-check = "cromwell-ping-me-bucket"
        }
    }
}
``` 
With this one:
```
services {
    HealthMonitor {
        config {
            check-dockerhub: true
            check-engine-database: true
            check-gcs: true
            check-papi-backends: [PAPIv1, PAPIv2]

            google-auth-name = service-account
            gcs-bucket-to-check = "cromwell-ping-me-bucket"
    }
  }
}
``` 
### Workflow options changes

A new workflow option is added. If the `final_workflow_outputs_dir` is set 
`use_relative_output_paths` can be used. When set to `true` this will copy 
all the outputs relative to their execution directory. 
my_final_workflow_outputs_dir/~~MyWorkflow/af76876d8-6e8768fa/call-MyTask/execution/~~output_of_interest.
More information can be found in [the workflow options documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/#output-copying).

### Bug fixes

#### WDL 1.0 strings can contain escaped quotes

For example, the statement `String s = "\""` is now supported, whereas previously it produced a syntax error.

#### Empty call blocks in WDL 1.0

Cromwell's WDL 1.0 implementation now allows empty call blocks, e.g. `call task_with_no_inputs {}`. This brings 1.0 in line with draft-2, which has always supported this syntax.

#### Packed CWL bugfix

Fixed a bug that caused an error like `Custom type was referred to but not found` to be issued when using an imported type as a `SchemaDefRequirement` in packed CWL.

## 39 Release Notes

### Cromwell ID changes

When set, the configuration value of `system.cromwell_id` will now have a random suffix appended, unless the
configuration key `system.cromwell_id_random_suffix` is set to `false`.

The generated id also appears more places in the logs, including when picking up workflows from the database and during
shutdown.

### Bug fixes

#### Format fix for `write_map()` 

Fixed an issue that caused the `write_map()` function in Cromwell's WDL 1.0 implementation to produce output in the wrong format. Specifically, the output's rows and columns were swapped. WDL draft-2 was not affected.
  
Incorrect `write_map()` output in Cromwell 38 and earlier:
```
key1    key2    key3
value1  value2  value3
```
Corrected `write_map()` output in Cromwell 39 and later:
```
key1  value1
key2  value2
key3  value3
```

## 38 Release Notes

### HPC paths with Docker

The `ConfigBackendLifecycleActorFactory` path variables `script`, `out` and `err` are now consistent when running with
and without docker. Similarly, when killing a docker task the `kill-docker` configuration key is now used instead of
`kill`. For more information see the [online documentation](https://cromwell.readthedocs.io/en/stable/backends/SGE/).

### No-op Health Monitor is now the default health monitor

Previous versions of Cromwell defaulted to using a health monitor service that checked Docker Hub and engine database status.
Neither check was useful if the `status` endpoint was never consulted as is likely the case in most deployments. Cromwell 38
now defaults to a `NoopHealthMonitorServiceActor` which does nothing. The previous health service implementation is still
available as `StandardHealthMonitorServiceActor`.

### Bug fixes
- Fixed an issue that could cause Cromwell to consume disk space unnecessarily when using zipped dependencies

#### HTTP responses

- When returning errors as json the `Content-Type` header is set to `application/json`.

## 37 Release Notes

### Docker

- Adds support for retrieving docker digests of asia.gcr.io images
- Adds configuration settings for docker digest lookups. See the `docker` section of the `reference.conf` for more information 
- Attempt to automatically adjust the boot disk size on the Google Cloud Backend (version 2) if the size of the image is greater than the default disk size or the required disk size in the runtime attributes.
Only works for registries that support the version 2 of the manifest schema (https://docs.docker.com/registry/spec/manifest-v2-2/)
At this date (12/09/18) this includes GCR and Dockerhub.

### Added new call cache path+modtime hashing strategy.

Call caching hashes with this new strategy are based on the path and the last modified time of the file.

### Instance independent abort

For multi-instance Cromwell deployments sharing a single database, earlier versions of Cromwell required abort
requests to be sent specifically to the instance that was running the targeted workflow. Cromwell 37 now
allows abort commands to be sent to any Cromwell instance in the shared-database deployment. Configuration details
[here](https://cromwell.readthedocs.io/en/develop/Configuring/abort-configuration).


### Call cache blacklisting

The Google Pipelines API (PAPI) version 1 and 2 backends now offer the option of call cache blacklisting on a per-bucket basis.
More info [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#call-cache-copy-authorization-failure-prefix-blacklisting).

### WDL

- All memory units in WDL are now treated as base-2.
For instance `1 KB == 1 KiB == 1024 Bytes`.

### Backend name for call caching purposes

Previous versions of Cromwell incorporated the name of the backend on which a call was run into the call cache hashes generated for that call.
Unfortunately this made it impossible to change the name of a backend without losing all previously run calls as potential cache hits.
Cromwell 37 introduces the `name-for-call-caching-purposes` backend configuration option as a means of decoupling the backend name from the
value used for the backend name for call caching purposes.

### CWL

Support `InputResourceRequirement` hint

### Changing configuration options

#### Logging Token Distribution

In cases where its not obvious why jobs are queued in Cromwell, you can enable logging for the Job Execution Token Dispenser, using
the `system.hog-safety.token-log-interval-seconds` configuration value.

The default, `0`, means that no logging will occur. 

#### HTTP Filesystem

- The HTTP filesystem is now enabled for engine use by default. To continue without an HTTP filesystem, you can add the 
following content into the appropriate stanza of your configuration file:
```
engine {
  filesystems {
    http { 
      enabled: false 
    }
  }
}
``` 
- When the value `exit-code-timeout-seconds` is set, `check-alive` command is now only called once every timeout interval instead of each poll.

### Beta preview of new Womtool `/describe` endpoint

This new endpoint brings the functionality of Womtool to the world of web services. Submit workflows for validation and receive a JSON description in response.

The endpoint is still undergoing heavy development and should not be used in production. The final version will ship in a future release of Cromwell; watch this space.   

### Bug fixes

- Fixed a regression in Cromwell 36 that could cause operations on empty arrays to fail with a spurious type error (closes [#4318](https://github.com/broadinstitute/cromwell/issues/4318))

#### Abort On Hold Workflows

On Hold workflows may now be aborted.

#### Command fixes for AWS and TES

The AWS and TES backends can now handle calls that generate longer command lines. Like the other
backends, commands scripts are first written to a file, the file is downloaded to the execution
host, and then the localized script is run.

Also fixed are AWS `command {}` blocks that use `|` at the start of a line. For example:

```
command {
  echo hello world \
  | cat
}
```

## 36 Release Notes

### Extra configuration options

The value `exit-code-timeout-seconds` can now set in a backend configuration.
Details [here](https://cromwell.readthedocs.io/en/develop/backends/HPC/#exit-code-timeout)

### [AWS S3 file transfers are now encrypted](https://github.com/broadinstitute/cromwell/pull/4264)

### Bug fixes

#### Metadata Request Coalescing

Coalesce metadata requests to eliminate expensive and redundant queries and metadata construction.

#### Eliminate redundant SFS logging and metadata 

Eliminate superfluous logging and metadata publishing in the shared filesystem backend on poll intervals where there was not a state change.

#### AWS region configuration respected throughout

Previously US-EAST-1 was hardcoded in places.

## 35 Release Notes

### Submit workflow using URL

Cromwell now allows for a user to submit the URL pointing to workflow file to run a workflow.
More details on how to use it in: 
- `Server` mode can be found [here](https://cromwell.readthedocs.io/en/develop/api/RESTAPI/).
- `Run` mode can be found [here](https://cromwell.readthedocs.io/en/develop/CommandLine/#run).

### Languages

- Added an opt-in namespace cache for the WDL Draft 2 language factory. Please see the Cromwell example configuration for details. NOTE: if upgrading from a hotfix version of Cromwell
that relied upon this cache, the cache is now opt-in and must be turned on explicitly in config.
- To maintain conformance with the OpenWDL spec, Cromwell drops support for the `version draft-3` identifier in this release. In the rare case where end users may have been using `version draft-3`, `version 1.0` is a drop-in replacement with no effect on functionality.

### HTTP Workflow Inputs for Shared File System and Google Pipelines API Version 2 Backends

`http` and `https` workflow inputs are now supported for shared filesystem and Google Pipelines API (PAPI) version 2
backends. Configuration details are described [here](http://cromwell.readthedocs.io/en/develop/filesystems/HTTP).

### Call cache hint support

More efficient cache hit copying in multi-user environments is now supported through the `call_cache_hit_path_prefixes` workflow option.
Details [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#call-cache-hit-path-prefixes)

### Root workflow level file hash caching support

Cromwell now offers the ability to cache file hashes on a root workflow level basis, details [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#file-hash-caching).

### Extra configuration options

The value `dockerRoot` can now be set in a backend configuration. 
This will set the execution folder in the container (default: `/cromwell-executions`).

### Bug Fixes

#### API
- The `releaseHold` endpoint will now return `404 Not Found` for an unrecognized workflow ID and `400 Bad Request` for a malformed or invalid workflow ID.

#### Languages

- Fixed a bug that allowed values to be "auto-boxed" into a single-element `Array` of that type, which is not allowed in the WDL spec (Closes [#3478](https://github.com/broadinstitute/cromwell/issues/3478)).

#### PAPI version 1

- Restored standard output and error streaming for jobs.

## 34 Release Notes

### Query API

* Fixes a bug which stopped `includeSubworkflow=false` from paging correctly and subworkflows from being discounted correctly from `totalResultsCount`.
* Query results will now be returned in reverse chronological order, with the most-recently submitted workflows returned first.

### Requester Pays on GCS

Access of Google Cloud Storage buckets with Requester Pays enabled is now supported.
Please read the [relevant documentation](http://cromwell.readthedocs.io/en/develop/filesystems/GoogleCloudStorage#requester-pays) for information on how to enable it and the consequences.

### Private Docker Support on Pipelines API v2

Support for private Docker Hub images is now included in the Google Pipelines API v2 backend. PAPI v2 private Docker support is
equivalent to that in PAPI v1 but the configuration differs, please see
[Docker configuration](http://cromwell.readthedocs.io/en/develop/filesystems/Google#Docker) for more details.

### Updated MySQL client with 8.0 support

Updated the MySQL connector client from `5.1.42` to `5.1.46` which adds support for connecting to MySQL 8.0. See the
documentation on [Changes in MySQL Connector/J](https://dev.mysql.com/doc/relnotes/connector-j/5.1/en/news-5-1.html) for
more information.

## 33 Release Notes

### Query endpoint

#### Exclude workflows based on Labels

This gives the ability to **filter out** workflows based on labels. Two new parameters called `excludeLabelAnd` and `excludeLabelOr` can be used for this purpose.
More details on how to use them can be found [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Include/Exclude subworkflows

Cromwell now supports excluding subworkflows from workflow query results using the `includeSubworkflows` parameter. By default they are included in the results.
More information can be found at [REST API](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Query workflows by Submission time

Cromwell now supports querying workflows by submission time. This will help find workflows that are submitted but not started yet (i.e. workflows which are
in On Hold state). More information can be found [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Submission time in Workflow Query Response

Submission time of a workflow is now included in WorkflowQueryResult, which is part of the response for workflow query.

### File Localization (NIO) Hint

Cromwell now allows tasks in WDL 1.0 can now specify an optimization in their `parameter_meta` that some `File` inputs do not need to be localized for the task to run successfully.
Full details are available in the [documentation page for this optimization](http://cromwell.readthedocs.io/en/develop/optimizations/FileLocalization).

### Bug Fixes

Workflows which are in 'On Hold' state can now be fetched using the query endpoint.

## 32 Release Notes

### Backends

#### Pipelines API V2
Initial support for Google [Pipelines API version 2](https://cloud.google.com/genomics/reference/rest/).
Expect feature parity except for private dockerhub images which are not supported at the moment, but will be in the near future.
Additionally, the "refresh token" authentication mode is **NOT** supported on PAPI V2.

In addition, the following changes are to be expected:
* Error messages for failed jobs might differ from V1
* The Pipelines API log file content might differ from V1

**Important (If you're running Cromwell with a Google backend, read this)**:
The `actor-factory` value for the google backend (`cromwell.backend.impl.jes.JesBackendLifecycleActorFactory`) is being deprecated.
Please update your configuration accordingly.

| PAPI Version  |                                 actor-factory                                |
|---------------|:----------------------------------------------------------------------------:|
|      V1       | cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory |
|      V2alpha1 | cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory |
|      V2beta   | cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory   |

If you don't update the `actor-factory` value, you'll get a deprecation warning in the logs, and Cromwell will default back to **PAPI V1**

### Task Retries
Cromwell now supports retrying failed tasks up to a specified count by declaring a value for the [maxRetries](RuntimeAttributes.md#maxRetries) key through the WDL runtime attributes.

### Labels
* Cromwell has removed most of the formatting restrictions from custom labels. Please check the [README](README.md#label-format) for more detailed documentation.
* Custom labels won't be submitted to Google backend as they are now decoupled from Google's default labels.
* Cromwell now publishes the labels as soon as the workflow is submitted (whether started or on hold). If the labels are invalid, the workflow will not be submitted and request will fail.

### Scala 2.11 Removed
From version 32 onwards we will no longer be publishing build artifacts compatible with Scala 2.11. 

* If you don't import the classes into your own scala project then this should have no impact on you.
* If you **are** importing the classes into your own scala project, make sure you are using Scala 2.12.

### Input Validation
Cromwell can now validate that your inputs files do not supply inputs with no impact on the workflow. Strict validation will be disabled by default in WDL draft 2 and CWL but enabled in WDL draft 3. See the 'Language Factory Config' below for details.

### Language Factory Config
All language factories can now be configured on a per-language-version basis. All languages and versions will support the following options:
* `enabled`: Defaults to `true`. Set to `false` to disallow workflows of this language and version.
* `strict-validation`: Defaults to `true` for WDL draft 3 and `false` for WDL draft 2 and CWL. Specifies whether workflows fail if the inputs JSON (or YAML) file contains values which the workflow did not ask for (and will therefore have no effect). Additional strict checks may be added in the future.

### API

* More accurately returns 503 instead of 500 when Cromwell can not respond in a timely manner
* Cromwell now allows a user to submit a workflow but in a state where it will not automatically be picked up for execution. This new state is called 'On Hold'. To do this you need to set the parameter workflowOnHold to true while submitting the workflow.
* API end point 'releaseHold' will allow the user to send a signal to Cromwell to allow a workflow to be startable, at which point it will be picked up by normal execution schemes.

### GPU

The PAPI backend now supports specifying GPU through WDL runtime attributes:

```wdl
runtime {
    gpuType: "nvidia-tesla-k80"
    gpuCount: 2
    zones: ["us-central1-c"]
}
```

The two types of GPU supported are `nvidia-tesla-k80` and `nvidia-tesla-p100`

**Important**: Before adding a GPU, make sure it is available in the zone the job is running in: https://cloud.google.com/compute/docs/gpus/

### Job Shell

Cromwell now allows for system-wide or per-backend job shell configuration for running user commands rather than always
using the default `/bin/bash`. To set the job shell on a system-wide basis use the configuration key `system.job-shell` or on a
per-backend basis with `<config-key-for-backend>.job-shell`. For example:

```
# system-wide setting, all backends get this
-Dsystem.job-shell=/bin/sh
```

```
# override for just the Local backend
-Dbackend.providers.Local.config.job-shell=/bin/sh
```

For the Config backend the value of the job shell will be available in the `${job_shell}` variable. See Cromwell's `reference.conf` for an example
of how this is used for the default configuration of the `Local` backend.

### Bug Fixes

The imports zip no longer unpacks a single (arbitrary) internal directory if it finds one (or more). Instead, import statements should now be made relative to the base of the import zip root.

#### Reverting Custom Labels

Reverting to a prior custom label value now works.

["Retrieves the current labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow)
will return the most recently summarized custom label value.

The above endpoint may still return the prior value for a short period of time after using
["Updated labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#update-labels-for-a-workflow)
until the background metadata summary process completes.

#### Deleting Duplicate Custom Label Rows

If you never used the REST API to revert a custom label back to a prior value you will not be affected. This only applies to workflows previously updated using
["Updated labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#update-labels-for-a-workflow).

The database table storing custom labels will delete duplicate rows for any workflow label key. For efficiency purposes
the values are not regenerated automatically from the potentially large metadata table.

In rare cases where one tried to revert to a prior custom label value you may continue to see different results
depending on the REST API used. After the database update
["Retrieves the current labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow)
will return the most-recent-unique value while
["Get workflow and call-level metadata for a specified workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#get-workflow-and-call-level-metadata-for-a-specified-workflow)
will return the up-to-date value. For example, if one previously updated a value from `"value-1"` > `"value-2"` >
`"value-3"` > `"value-2"` then the former REST API will return `value-3` while the latter will return `value-2`.

#### Workflow options `google_project` output in metadata

Workflow metadata for jobs run on a Google Pipelines API backend will report the `google_project` specified via a
[workflow options json](http://cromwell.readthedocs.io/en/develop/wf_options/Google/#google-pipelines-api-workflow-options).

## 31 Release Notes

* **Cromwell server**  
The Cromwell server source code is now located under `server/src`. `sbt assembly` will build the runnable Cromwell JAR in 
`server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

* **Robustness**
    + The rate at which jobs are being started can now be controlled using the `system.job-rate-control` configuration stanza.  
    + A load controller service has been added to allow Cromwell to self-monitor and adjust its load accordingly.
The load controller is currently a simple on/off switch controlling the job start rate. It gathers metrics from different parts of the system
to inform its decision to stop the creation of jobs.
You can find relevant configuration in the `services.LoadController` section of the `cromwell.examples.conf` file,
as well as in the `load-control` section in `reference.conf`.
The load level of the monitored sub-systems are instrumented and can be found under the `cromwell.load` statsD path.
    + The statsD metrics have been re-shuffled a bit. If you had a dashboard you might find that you need to update it.
Changes include: 
        + Removed artificially inserted "count" and "timing" the path
        + Added a `load` section
        + Metrics were prefixed twice with `cromwell` (`cromwell.cromwell.my_metric`), now they're only prefixed once
        + Added `processed` and `queue` metrics under various metrics monitoring the throughput and amount of queued work respectively
        + Added a memory metric representing an estimation of the free memory Cromwell thinks it has left

* Added a configuration option under `docker.hash-lookup.enabled` to disable docker hash lookup.
 Disabling it will also disable call caching for jobs with floating docker tags.
 
* **API**    
    + Updated the `/query` response to include the total number of query results returned. See [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#workflowqueryresponse) for more information.

## 30.1 Release Notes

* A set of bug fixes following the migration of Cromwell to WOM (the Workflow Object Model) in version 30.

## 30 Release Notes

### Breaking changes

* The `customLabels` form field for workflow submission has been renamed to `labels`.

### Other changes

* **New Cromwell documentation**  
Our documentation has moved from our [README](https://github.com/broadinstitute/cromwell/blob/29_hotfix/README.md) to a new website: [Cromwell Documentation](http://cromwell.readthedocs.io/en/develop/). There are new [Tutorials](http://cromwell.readthedocs.io/en/develop/tutorials/FiveMinuteIntro/) and much of the documentation has been re-written. The source files are in the [/docs](https://github.com/broadinstitute/cromwell/tree/develop/docs) directory.

* **API**  
    + Cromwell now supports input files in the yaml format (JSON format is still supported).
    + Added a [GET version for the `labels` endpoint](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow) which will return current labels for a workflow.

* **Database**  
You have the option of storing the metadata in a separate SQL database than the database containing the internal engine
data. When switching connection information for an existing database containing historical data, the tables
should be manually replicated from one database instance to another using the tools appropriate for your specific
database types. Cromwell will not move any existing data automatically. This feature should be considered experimental
and likely to change in the future. See the [Database Documentation](https://cromwell.readthedocs.io/en/develop/Configuring/#database) or the `database` section in
[cromwell.examples.conf](https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf) for more
information.

* **StatsD**  
Added initial support for StatsD instrumentation. See the [Instrumentation Documentation](https://cromwell.readthedocs.io/en/develop/Instrumentation) for details on how to use it.

* **User Service Account auth mode for Google**  
Added a new authentication mode for [Google Cloud Platform](https://cromwell.readthedocs.io/en/develop/backends/Google) which will allow a user to supply the JSON key file in their workflow options to allow for per-workflow authentication via service account. This is analogous to the previously existing refresh token authentication scheme. As with the refresh token scheme it is encouraged that the **user_service_account_json** workflow option field is added to the **encrypted-fields** list in the configuration.

* **Bugfixes**  
Abort of Dockerized tasks on the Local backend should now work as expected. Cromwell uses `docker kill` to kill the Docker container.

## 29 Release Notes

### Breaking Changes

* **Command line**  
In preparation for supporting CWL scripts (yes, you read that right!), we have extensively revised the Command Line in Cromwell 29. For more details about the usage changes please see the [README](https://github.com/broadinstitute/cromwell#command-line-usage). And stay tuned to the [WDL/Cromwell blog](https://software.broadinstitute.org/wdl/blog) over the next couple of months for more news about CWL.

* **Request timeouts**   
Cromwell now returns more specific `503 Service Unavailable` error codes on request timeouts, rather than the more generic `500 Internal Server Error`. The response for a request timeout will now be plain text, rather than a JSON format.

* **Metadata endpoint**  
The response from the metadata endpoint can be quite large depending on your workflow. You can now opt-in to have Cromwell gzip your metadata file, in order to reduce file size, by sending the `Accept-Encoding: gzip` header. The default behavior now does not gzip encode responses.

* **Engine endpoints**  
Previously the engine endpoints were available under `/api/engine` but now the endpoints are under `/engine` so they don't require authentication. Workflow endpoints are still available under `/api/workflows`. We also deprecated the setting `api.routeUnwrapped` as a part of this internal consistency effort.

* **Call caching diff**  
We updated the response format of the [callcaching/diff](https://github.com/broadinstitute/cromwell#get-apiworkflowsversioncallcachingdiff) endpoint.

### Other changes

* **Cromwell server**  
When running in server mode, Cromwell now attempts to gracefully shutdown after receiving a `SIGINT` (`Ctrl-C`) or `SIGTERM` (`kill`) signal. This means that Cromwell waits for all pending database writes before exiting, as long as you include `application.conf` at the top of your config file. You can find detailed information about how to configure this feature in the [Cromwell Wiki](https://github.com/broadinstitute/cromwell/wiki/DevZone#graceful-server-shutdown).

* **Concurrent jobs**  
You can now limit the number of concurrent jobs for any backend. Previously this was only possible in some backend implementations. Please see the [README](https://github.com/broadinstitute/cromwell#backend-job-limits) for details.

### WDL

* **Optional WDL variables**  
Empty optional WDL values are now rendered as the `null` JSON value instead of the JSON string `"null"` in the metadata and output endpoints. You do not need to migrate previous workflows. Workflows run on Cromwell 28 and prior will still render empty values as `"null"`.

* **Empty WDL variables**  
Cromwell now accepts `null` JSON values in the input file and coerces them as an empty WDL value. WDL variables must be declared optional in order to be supplied with a `null` JSON value.

input.json
```json
{
    "null_input_values.maybeString": null,
    "null_input_values.arrayOfMaybeInts": [1, 2, null, 4]
}
```

workflow.wdl
```
workflow null_input_values {
    String? maybeString
    Array[Int?] arrayOfMaybeInts
}
```

## 28

### Bug Fixes

#### WDL write_* functions add a final newline

The following WDL functions now add a newline after the final line of output (the previous behavior of not adding this
newline was inadvertent):
- `write_lines`
- `write_map`
- `write_object`
- `write_objects`
- `write_tsv`

For example:

```
task writer {
  Array[String] a = ["foo", "bar"]
  command {
    # used to output: "foo\nbar"
    # now outputs: "foo\nbar\n"
    cat write_lines(a)
  }
}
```

#### `ContinueWhilePossible`

A workflow utilizing the WorkflowFailureMode Workflow Option `ContinueWhilePossible` will now successfully reach a terminal state once all runnable jobs have completed.
#### `FailOnStderr` 
When `FailOnStderr` is set to false, Cromwell no longer checks for the existence of a stderr file for that task. 

### WDL Functions

#### New functions: floor, ceil and round:

Enables the `floor`, `ceil` and `round` functions in WDL to convert floating point numbers to integers.

For example we can now use the size of an input file to influence the amount of memory the task is given. In the example below a 500MB input file will result in a request for a VM with 2GB of memory:

```
task foo {
    File in_file
    command { ... }
    runtime {
      docker: "..."
      memory: ceil(size(in_file)) * 4 
    }
}
```

### Call Caching

* Hash values calculated by Cromwell for a call when call caching is enabled are now published to the metadata.
It is published even if the call failed. However if the call is attempted multiple times (because it has been preempted for example),
since hash values are strictly identical for all attempts, they will only be published in the last attempt section of the metadata for this call.
If the hashes fail to be calculated, the reason is indicated in a `hashFailures` field in the `callCaching` section of the call metadata.
*Important*: Hashes are not retroactively published to the metadata. Which means only workflows run on Cromwell 28+ will have hashes in their metadata.

See the [README](https://github.com/broadinstitute/cromwell#get-apiworkflowsversionidmetadata) for an example metadata response.

* New endpoint returning the hash differential for 2 calls. 

`GET /api/workflows/:version/callcaching/diff`

See the [README](https://github.com/broadinstitute/cromwell#get-apiworkflowsversioncallcachingdiff) for more details.

### Workflow Submission

* The workflow submission parameters `wdlSource` and `wdlDependencies` have been deprecated in favor of `workflowSource` and
`workflowDependencies` respectively.  The older names are still supported in Cromwell 28 with deprecation warnings but will
be removed in a future version of Cromwell.

### Labels
* A new `/labels` endpoint has been added to update labels for an existing workflow. See the [README](README.md#patch-apiworkflowsversionidlabels) for more information.
* Label formatting requirements have been updated, please check the [README](README.md#label-format) for more detailed documentation.


### JES Backend

The JES backend now supports a `filesystems.gcs.caching.duplication-strategy` configuration entry.
It can be set to specify the desired behavior of Cromwell regarding call outputs when a call finds a hit in the cache.
The default value is `copy` which will copy all output files to the new call directory.
A second value is allowed, `reference`, that will instead point to the original output files, without copying them.


```hocon
filesystems {
  gcs {
    auth = "application-default"
    
    caching {
      duplication-strategy = "reference"
    }
  }
}
```

A placeholder file will be placed in the execution folder of the cached call to explain the absence of output files and point to the location of the original ones.


### Metadata Write Batching

Metadata write batching works the same as in previous versions of Cromwell, but the default batch size has been changed from 1 to 200.  It's possible that 200 is too high in some environments, but 200 is more likely to be an appropriate value
than the previous default.


## 27

### Migration

* Call Caching has been improved in this version of Cromwell, specifically the time needed to determine whether or not a job can be cached
 has drastically decreased. To achieve that the database schema has been modified and a migration is required in order to preserve the pre-existing cached jobs.
 This migration is relatively fast compared to previous migrations. To get an idea of the time needed, look at the size of your `CALL_CACHING_HASH_ENTRY` table.
 As a benchmark, it takes 1 minute for a table with 6 million rows.
 The migration will only be executed on MySQL. Other databases will lose their previous cached jobs.
 In order to run properly on MySQL, **the following flag needs to be adjusted**: https://dev.mysql.com/doc/refman/5.5/en/server-system-variables.html#sysvar_group_concat_max_len
 The following query will give you a minimum to set the group_concat_max_len value to:
 
 ```sql
SELECT MAX(aggregated) as group_concat_max_len FROM
      (
            SELECT cche.CALL_CACHING_ENTRY_ID, SUM(LENGTH(CONCAT(cche.HASH_KEY, cche.HASH_VALUE))) AS aggregated
            FROM CALL_CACHING_HASH_ENTRY cche
            GROUP BY cche.CALL_CACHING_ENTRY_ID
      ) aggregation
 ```

 Here is the SQL command to run to set the group_concat_max_len flag to the proper value:
 
 ```sql
SET GLOBAL group_concat_max_len = value
 ```
 
 Where `value` is replaced with the value you want to set it to.
 
 Note that the migration will fail if the flag is not set properly.
 
### Breaking Changes

* The update to Slick 3.2 requires a database stanza to
[switch](http://slick.lightbend.com/doc/3.2.0/upgrade.html#profiles-vs-drivers) from using `driver` to `profile`.

```hocon
database {
  #driver = "slick.driver.MySQLDriver$" #old
  profile = "slick.jdbc.MySQLProfile$"  #new
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

### Call Caching

Cromwell now supports call caching with floating Docker tags (e.g. `docker: "ubuntu:latest"`). Note it is still considered
a best practice to specify Docker images as hashes where possible, especially for production usages.

Within a single workflow Cromwell will attempt to resolve all floating tags to the same Docker hash, even if Cromwell is restarted
during the execution of a workflow. In call metadata the `docker` runtime attribute is now the same as the
value that actually appeared in the WDL:

```
   "runtimeAttributes": {
     "docker": "ubuntu:latest",
     "failOnStderr": "false",
     "continueOnReturnCode": "0"
   }
```

Previous versions of Cromwell rewrote the `docker` value to the hash of the Docker image.

There is a new call-level metadata value `dockerImageUsed` which captures the hash of the Docker image actually used to
run the call:

```
   "dockerImageUsed": "library/ubuntu@sha256:382452f82a8bbd34443b2c727650af46aced0f94a44463c62a9848133ecb1aa8"
```

### Docker

* The Docker section of the configuration has been slightly reworked 
An option to specify how a Docker hash should be looked up has been added. Two methods are available.
    "local" will try to look for the image on the machine where cromwell is running. If it can't be found, Cromwell will try to `pull` the image and use the hash from the retrieved image.
    "remote" will try to look up the image hash directly on the remote repository where the image is located (Docker Hub and GCR are supported)
Note that the "local" option will require docker to be installed on the machine running cromwell, in order for it to call the docker CLI.
* Adds hash lookup support for public [quay.io](https://quay.io/) images.

### WDL Feature Support
* Added support for the new WDL `basename` function. Allows WDL authors to get just the file name from a File (i.e. removing the directory path)
* Allows coercion of `Map` objects into `Array`s of `Pair`s. This also allows WDL authors to directly scatter over WDL `Map`s.

### Miscellaneous
* Adds support for JSON file format for google service account credentials. As of Cromwell 27, PEM credentials for PAPI are deprecated and support might be removed in a future version.

```
google {

  application-name = "cromwell"

  auths = [
    {
      name = "service-account"
      scheme = "service_account"
      json-file = "/path/to/file.json"
    }
  ]
}
```

### General Changes

* The `/query` endpoint now supports querying by `label`. See the [README](README.md#get-apiworkflowsversionquery) for more information.
* The `read_X` standard library functions limit accepted filesizes.  These differ by type, e.g. read_bool has a smaller limit than read_string.  See reference.conf for default settings.

## 26

### Breaking Changes

* Failure metadata for calls and workflows was being displayed inconsistently, with different formats depending on the originating Cromwell version. Failures will now always present as an array of JSON objects each representing a failure. Each failure will have a message and a causedBy field. The causedBy field will be an array of similar failure objects. An example is given below:

```
failures: [{
  message: "failure1",
  causedBy: [{
    message: "cause1",
    causedBy: []
   }, {
    message: "cause2",
    causedBy: []
  }]
 }, {
  message: "failure2",
  causedBy: []
}]
```

### Additional Upgrade Time

* Upgrading to Cromwell 26 will take additional time due to the migration of failure metadata. Cromwell will automatically run a database query during the upgrade which appears to be roughly linear to the number of rows in the METADATA_ENTRY table. You can estimate upgrade time using the following equation: `time to migrate (in seconds) ~= (rows in METADATA_ENTRY) / 65000` Note that due to differences in hardware and database speed, this is only a rough estimate.

### Config Changes

* Added a configuration option under `system.io` to throttle the number of I/O queries that Cromwell makes, as well as configure retry parameters.
 This is mostly useful for the JES backend and should be updated to match the GCS quota available for the project.
 
```
system.io {
  # Global Throttling - This is mostly useful for GCS and can be adjusted to match
  # the quota availble on the GCS API
  number-of-requests = 100000
  per = 100 seconds
  
  # Number of times an I/O operation should be attempted before giving up and failing it.
  number-of-attempts = 5
}
```

## 25

### External Contributors
* A special thank you to @adamstruck, @antonkulaga and @delocalizer for their contributions to Cromwell.
### Breaking Changes

* Metadata keys for call caching are changed. All call caching keys are now in a `callCaching` stanza. `Call cache read result` has moved here and is now `result`. The `allowResultReuse` and `effectiveCallCachingMode` have moved here. The `hit` boolean is a simple indication of whether or not it was a hit, with no additional information. An example using the new format is:
```
"callCaching": {
  "hit": false,
  "effectiveCallCachingMode": "ReadAndWriteCache",
  "result": "Cache Miss",
  "allowResultReuse": true
}
```

### Config Changes

* Added a field `insert-batch-size` to the `database` stanza which defines how many values from a batch insert will be processed at a time. This value defaults to 2000. 
* Moved the config value `services.MetadataService.metadata-summary-refresh-interval` to `services.MetadataService.config.metadata-summary-refresh-interval`
* Added ability to override the default zone(s) used by JES via the config structure by setting `genomics.default-zones` in the JES configuration
* The cromwell server TCP binding timeout is now configurable via the config key `webservice.binding-timeout`, defaulted
  to the previous value `5s` (five seconds) via the reference.conf.
* For MySQL users, a massive scalability improvement via batched DB writing of internal metadata events. Note that one must add `rewriteBatchedStatements=true` to their JDBC URL in their config in order to take advantage of this

### General Changes

* Cromwell's WDL parser now recognizes empty array literals correctly, e.g. `Array[String] emptyArray = []`.
* Cromwell now applies default labels automatically to JES pipeline runs.
* Added support for new WDL functions:
  * `length: (Array[X]) => Integer` - report the length of the specified array
  * `prefix: (String, Array[X]) => Array[String]` - generate an array consisting of each element of the input array prefixed
     by a specified `String`.  The input array can have elements of any primitive type, the return array will always have
     type `Array[String]`.
  * `defined: (Any) => Boolean` - Will return false if the provided value is an optional that is not defined. Returns true in all other cases.
* Cromwell's Config (Shared Filesystem) backend now supports invocation of commands which run in a Docker image as a non-root user.
  The non-root user could either be the default user for a given Docker image (e.g. specified in a Dockerfile via a `USER` directive),
  or the Config backend could pass an optional `"-u username"` as part of the `submit-docker` command.
* In some cases the SFS backend, used for Local, SGE, etc., coerced `WdlFile` to `WdlString` by using `.toUri`. This
resulted in strings prepended with `file:///path/to/file`. Now absolute file paths will not contain the uri scheme.
* Launch jobs on servers that support the GA4GH Task Execution Schema using the TES backend.
* **Call caching: Cromwell will no longer try to use the cache for WDL tasks that contain a floating docker tag.** 
  Call caching will still behave the same for tasks having a docker image with a specific hash.
  See https://github.com/broadinstitute/cromwell#call-caching-docker-tags for more details. 
* Added docker hash lookup. Cromwell will try to lookup the hash for a docker image with a floating tag, and use that hash when executing the job.
  This will be reflected in the metadata where the docker runtime attribute will contains the hash that was used.
  If Cromwell is unable to lookup the docker hash, the job will be run with the original user defined floating tag.
  Cromwell is currently able to lookup public and private docker hashes for images on Docker Hub and Google Container Engine for job running on the JES backend.
  For other backends, cromwell is able to lookup public docker hashes for Docker Hub and Google Container Engine.
  See https://github.com/broadinstitute/cromwell#call-caching-docker-tags for more details. 

### Database schema changes
* Added CUSTOM_LABELS as a field of WORKFLOW_STORE_ENTRY, to store workflow store entries.

## 24

* When emitting workflow outputs to the Cromwell log only the first 1000 characters per output will be printed
* Added support for conditional (`if`) statements.
* Globs for Shared File System (SFS) backends, such as local or SGE, now use bash globbing instead of Java globbing, consistent with the JES backend.

## 23

* The `meta` and `parameter_meta` blocks are now valid within `workflow` blocks, not just `task`
* The JES backend configuration now has an option `genomics-api-queries-per-100-seconds` to help tune the rate of batch polling against the JES servers. Users with quotas larger than default should make sure to set this value.
* Added an option `call-caching.invalidate-bad-cache-results` (default: `true`). If true, Cromwell will invalidate cached results which have failed to copy as part of a cache hit.
* Timing diagrams and metadata now receive more fine grained workflow states between submission and Running.
* Support for the Pair WDL type (e.g. `Pair[Int, File] floo = (3, "gs://blar/blaz/qlux.txt")`)
* Added support for new WDL functions:
  * `zip: (Array[X], Array[Y]) => Array[Pair[X, Y]]` - align items in the two arrays by index and return them as WDL pairs 
  * `cross: (Array[X], Array[Y]) => Array[Pair[X, Y]]` - create every possible pair from the two input arrays and return them all as WDL pairs
  * `transpose: (Array[Array[X]]) => Array[Array[X]]` compute the matrix transpose for a 2D array. Assumes each inner array has the same length.
* By default, `system.abort-jobs-on-terminate` is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <wdl> <inputs>`.
* Enable WDL imports when running in Single Workflow Runner Mode.
* Both batch and non-batch REST workflow submissions now require a multipart/form-data encoded body.
* Support for sub workflows (see [Annex A](#annex-a---workflow-outputs))
* Enable WDL imports when running in Single Workflow Runner Mode as well as Server Mode
* Support for WDL imports through an additional imports.zip parameter
* Support for sub workflows
* Corrected file globbing in JES to correctly report all generated files. Additionally, file globbing in JES now uses bash-style glob syntax instead of python style glob syntax
* Support declarations as graph nodes
* Added the ability to override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`. More details can be found in the README.md
* Fix bugs related to the behavior of Cromwell in Single Workflow Runner Mode. Cromwell will now exit once a workflow completes in Single Workflow Runner Mode. Additionally, when restarting Cromwell in Single Workflow Runner Mode, Cromwell will no longer restart incomplete workflows from a previous session.

### Annex A - Workflow outputs
    
The WDL specification has changed regarding [workflow outputs](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#outputs) to accommodate sub workflows.
This change is backward compatible in terms of runnable WDLs (WDL files using the deprecated workflow outputs syntax will still run the same). 
The only visible change lies in the metadata (as well as the console output in single workflow mode, when workflow outputs are printed out at the end of a successful workflow).

TL;DR Unless you are parsing or manipulating the "key" by which workflow outputs are referenced in the metadata (and/or the console output for single workflow mode), you can skip the following explanation.

*Metadata Response*
```
{
  ...
  outputs {
    "task_output_1": "hello",
    "task_output_2": "world"
            ^
       If you don't manipulate this part of the metadata, then skip this section
  }
}
```

In order to maintain backward compatibility, workflow outputs expressed with the deprecated syntax are "expanded" to the new syntax. Here is an example:

```
task t {
    command {
        #do something
    }
    output {
        String out1 = "hello"
        String out2 = "world"
    }
}
```

```
    workflow old_syntax {
        call t
        output {
            t.*
        }
    }
```

```
    workflow new_syntax {
        call t
        output {
            String wf_out1 = t.out1
            String wf_out2 = t.out2
        }
    }
```

The new syntax allows for type checking of the outputs as well as expressions. It also allows for explicitly naming to the outputs.
The old syntax doesn't give the ability to name workflow outputs. For consistency reasons, Cromwell will generate a "new syntax" workflow output for each task output, and name them.
Their name will be generated using their FQN, which would give 

```
output {
   String w.t.out1 = t.out1
   String w.t.out2 = t.out2
}
```
        
However as the FQN separator is `.`, the name itself cannot contain any `.`. 
For that reason, `.` are replaced with `_` :

*Old syntax expanded to new syntax*
```
output {
   String w_t_out1 = t.out1
   String w_t_out2 = t.out2
}
```

The consequence is that the workflow outputs section of the metadata for `old_syntax` would previously look like 
 
 ```
    outputs {
        "w.t.out1": "hello",
        "w.t.out2": "hello"
    }
 ```
 
but it will now look like 

```
    outputs {
        "w_t_out1": "hello",
        "w_t_out2": "hello"
    }
```

The same applies for the console output of a workflow run in single workflow mode.


## 0.22

* Improved retries for Call Caching and general bug fixes.
* Users will experience better scalability of status polling for Google JES.
* Now there are configurable caching strategies for a SharedFileSystem backend (i.e. Local, SFS) in the backend's stanza:
  See below for detailed descriptions of each configurable key.

```
backend {
  ...
  providers {
    SFS_BackendName {
      actor-factory = ...
      config {
        ...
        filesystems {
          local {
            localization: [
               ...
            ]
            caching {
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]
              # Possible values: file, path
              # "file" will compute an md5 hash of the file content.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              hashing-strategy: "file"

              # When true, will check if a sibling file with the same name and the .md5 extension exists, and if it does, use the content of this file as a hash.
              # If false or the md5 does not exist, will proceed with the above-defined hashing strategy.
              check-sibling-md5: false
            }
```
* Multiple Input JSON files can now be submitted in server mode through the existing submission endpoint: /api/workflows/:version.
    This endpoint accepts a POST request with a multipart/form-data encoded body. You can now include multiple keys for workflow inputs.

        Each key below can contain an optional JSON file of the workflow inputs. A skeleton file can be generated from wdltool using the "inputs" subcommand.
        NOTE: In case of key conflicts between multiple JSON files, higher values of x in workflowInputs_x override lower values. For example, an input
        specified in workflowInputs_3 will override an input with the same name that was given in workflowInputs or workflowInputs_2. Similarly, an input
        specified in workflowInputs_5 will override an input with the same name in any other input file.

        workflowInputs
        workflowInputs_2
        workflowInputs_3
        workflowInputs_4
        workflowInputs_5

* You can now limit the number of concurrent jobs for a backend by specifying the following option in the backend's config stanza:
```
backend {
  ...
  providers {
    BackendName {
      actor-factory = ...
      config {
        concurrent-job-limit = 5
```


## 0.21

* Warning: Significant database updates when you switch from version 0.19 to 0.21 of Cromwell.
  There may be a long wait period for the migration to finish for large databases.
  Please refer to MIGRATION.md for more details.

* There are significant architectural changes related to increases in performance and scaling.

* The biggest user-facing changes from 0.19 to 0.21 are related to the application.conf file, which has been restructured significantly.
The configuration for backends now is all contained within a `backend` stanza, which specifies 1 stanza per name per backend and a default backend, as follows:

```
backend {
    default=Local
    providers {
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }
        }
        JES {
            actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }
        }
        SGE {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }r
        }
    }
}
```
* A new `/stats` endpoint has been added to get workflow and job count for a Cromwell running in server mode.

* Renamed Workflow Options:
   “workflow_log_dir” -> “final_workflow_log_dir”
    “call_logs_dir” -> “final_call_logs_dir”
    “outputs_path” -> “final_workflow_outputs_dir”
    “defaultRuntimeOptions” -> “default_runtime_attributes”

* Timing diagrams endpoint has been updated to include additional state information about jobs.

* Add support for Google Private IPs through `noAddress` runtime attribute. If set to true, the VM will NOT be provided with a public IP address.
*Important*: Your project must be whitelisted in "Google Access for Private IPs Early Access Program". If it's not whitelisted and you set this attribute to true, the task will hang.
  Defaults to `false`.
  e.g:
```
task {
    command {
        echo "I'm private !"
    }

    runtime {
        docker: "ubuntu:latest"
        noAddress: true
    }
}
```

* The Local and the SGE backend have been merged into a generic
Shared File System (SFS) backend. This updated backend can be configured
to work with various other command line dispatchers such as LSF. See the
[README](README.md#sun-gridengine-backend) for more info.

* On the JES and SFS backends, task `command` blocks are now always
passed absolute paths for input `File`s.

* On the SFS backends, the call directory now contains two sub-directories:
    * `inputs` contains all the input files that have been localized for this task (see next below for more details)
    * `execution` contains all other files (script, logs, rc, potential outputs etc...)

* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc.
* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc.

For example:
```
# use a mysql database
database {
  driver = "slick.driver.MySQLDriver$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

## 0.20

* The default per-upload bytes size for GCS is now the minimum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.

* Updated Docker Hub hash retriever to parse json with [custom media
types](https://github.com/docker/distribution/blob/05b0ab0/docs/spec/manifest-v2-1.md).

* Added a `/batch` submit endpoint that accepts a single wdl with
multiple input files.

* The `/query` endpoint now supports querying by `id`, and submitting
parameters as a HTTP POST.
# Citing `Cromwell`

If you use `Cromwell` in your work we would prefer it if you would use the following reference in your work.

## BibTeX

```bibtex
@misc{https://doi.org/10.7490/f1000research.1114634.1,
  doi = {10.7490/f1000research.1114634.1},
  url = {https://f1000research.com/slides/6-1381},
  author = {Voss,  Kate and Auwera,  Geraldine Van Der and Gentry,  Jeff},
  title = {Full-stack genomics pipelining with GATK4 + WDL + Cromwell [version 1; not peer reviewed]},
  journal = {ISCB Comm J},
  publisher = {F1000Research},
  type = {slides},
  volume = {6},
  number = {1381},
  year = {2017}
}
```

## Textual

Voss K, Van der Auwera G and Gentry J. Full-stack genomics pipelining with GATK4 + WDL + Cromwell
[version 1; not peer reviewed]. _F1000Research_ 2017, **6**(ISCB Comm J):1381 (slides)
(https://doi.org/10.7490/f1000research.1114634.1) 
[![Build Status](https://travis-ci.com/broadinstitute/cromwell.svg?branch=develop)](https://travis-ci.com/broadinstitute/cromwell?branch=develop)
[![codecov](https://codecov.io/gh/broadinstitute/cromwell/branch/develop/graph/badge.svg)](https://codecov.io/gh/broadinstitute/cromwell)

## Welcome to Cromwell

Cromwell is an open-source Workflow Management System for bioinformatics. Licensing is [BSD 3-Clause](LICENSE.txt).

The [Cromwell documentation has a dedicated site](https://cromwell.readthedocs.io/en/stable).

First time to Cromwell? Get started with [Tutorials](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

### Community

Thinking about contributing to Cromwell? Get started by reading our [Contributor Guide](CONTRIBUTING.md).

Cromwell has a growing ecosystem of community-backed projects to make your experience even better! Check out our [Ecosystem](https://cromwell.readthedocs.io/en/stable/Ecosystem/) page to learn more.

Talk to us:
- [Join the Cromwell Slack workspace](https://join.slack.com/t/cromwellhq/shared_invite/zt-dxmmrtye-JHxwKE53rfKE_ZWdOHIB4g) to discuss the Cromwell workflow engine.
- [Join the OpenWDL Slack workspace](https://join.slack.com/t/openwdl/shared_invite/zt-ctmj4mhf-cFBNxIiZYs6SY9HgM9UAVw) to discuss the evolution of the WDL language itself.
    - More information about WDL is available in [that project's repository](https://github.com/openwdl/wdl).  

### Capabilities and roadmap

A majority of Cromwell users today run their workflows in [Terra](https://app.terra.bio/), a fully-managed cloud-native bioinformatics computing platform. See [here](https://support.terra.bio/hc/en-us/articles/360036379771-Get-started-running-workflows) for a quick-start guide.

Users with specialized needs who wish to install and maintain their own Cromwell instances can [download](https://github.com/broadinstitute/cromwell/releases) a JAR or Docker image. The development team accepts reproducible bug reports from self-managed instances, but cannot feasibly provide direct support.

[Cromwell's backends](https://cromwell.readthedocs.io/en/stable/backends/Backends/) receive development resources proportional to customer demand. The team is actively developing for Google Cloud and AWS. Maintenance of other backends is primarily community-based.

Cromwell [supports](https://cromwell.readthedocs.io/en/stable/LanguageSupport/) the WDL and CWL workflow languages. The Cromwell team is actively developing WDL, while maintenance for CWL is primarily community-based.  

### Security reports

If you believe you have found a security issue please contact `infosec@broadinstitute.org`.

### Issue tracking in JIRA

<!--
AEN external issue filing tested 2020-12-08 with `oednichols@gmail.com` / `https://broadworkbench.atlassian.net/browse/CROM-6681`
-->

Need to file an issue? Head over to [our JIRA](https://broadworkbench.atlassian.net/jira/software/c/projects/CROM/issues). You must create a free profile to view or create.

[Issues in Github](https://github.com/broadinstitute/cromwell/issues) remain available for discussion among community members but are not actively monitored by the development team.

![Cromwell JIRA](docs/img/cromwell_jira.png)

![Jamie, the Cromwell pig](docs/jamie_the_cromwell_pig.png)
### Contributing to Cromwell

Thank you for contributing to Cromwell. The project is sponsored by Broad Institute and has participants all over the world, who have collectively added tremendous value over the years. This page describes how we handle external contributions to maximize the impact of your work and reduce delays in getting your PRs accepted.

#### Run your idea by us first
If you're thinking of writing a non-trivial amount of code to solve a problem, we encourage you to reach out via an [issue](https://github.com/broadinstitute/cromwell/issues/new) to get feedback. It is likely we will have suggestions about how to proceed. It's also possible, though hopefully rare, that there is a hidden impediment that would prevent your solution from working. If we spot it at the idea stage, we can give you feedback much earlier than if we have to reject your pull request!

#### Maintenance considerations
The Cromwell team at Broad maintains and enhances the application to serve the needs of both Broad Institute and external users. Sometimes, we may identify a feature idea or pull request that works and is a good idea, but we may be unable to commit to maintaining it indefinitely. This may be because it does not align with the strategic direction of the project, or simply due to time constraints on the maintainers. Once again, it always helps to solicit early feedback.

#### Reviewing pull requests
Because pull requests require a substantial amount of time to review carefully, we prioritize and schedule them into our sprints alongside all of our other work. At present, the team operates on three-week sprints so if you happen to submit a PR early on in the sprint it may be a while before a team member has a chance to look at it. We realize this may be frustrating and strive to provide timely updates about PR status.
AWS Batch Backend Architecture
==============================

Overview
--------

The architecture of the code base follows very closely to the Google version.
Probably a little too closely, and lots of code was lifted from the Google
backend originally, then modified to work with AWS.

Fundamentally, Google Pipelines API (a.k.a. PAPI) works pretty differently from
AWS Batch. In Pipelines, all the infrastructure is completely managed by Google,
while AWS Batch exposes that infrastructure to a large degree so that customers
can fine tune it as necessary. An implementation that uses Fargate might be an 
alternative that is closer or an implementation that uses Step Functions although
that would be a separate backend.

From a Cromwell perspective, this means that unlike Pipelines, where
infrastructure details are defined in the configuration or the WDL, in AWS
Batch, these configuration details are handled outside. All the AWS Batch
backend needs to know is "what is the ARN for the job Queue"?

A good example of the difference can be seen in the 'disks' configuration. In
Pipelines, you need to specify the type of disk and size. In AWS, this will
be defined instead when you setup your environment (more on that later), so
all the AWS backend really needs to know is what mount points you need
defined.

This infrastructure and all the associated configuration still exists; however,
it is moved out of the Cromwell configuration.

AWS Batch
---------

Because AWS Batch is so different from PAPI, those familiar only with PAPI
would be best off with an overview of AWS Batch. If you are familiar with
the workings of Batch, feel free to skip this section, and move on.

[AWS Batch](https://aws.amazon.com/batch/) fundamentally is a service to allow batch jobs to run easily and
efficiently. To use it effectively, however, you need to understand its own
technical stack. To create a job, you need a "Job Queue". That job queue allows
jobs to be scheduled onto one or more "Compute Environments". This can
be managed through AWS Batch, but when AWS Batch sets up a compute environment,
it's simply setting up an Elastic Container Service (ECS) Cluster. The ECS
cluster, in turn is just a few managed CloudFormation templates, that is
controlling an AutoScaling group of EC2 instances.

What really makes an ECS instance an ECS instance is the presence of a configured
[Amazon ECS agent](https://github.com/aws/amazon-ecs-agent). This agent polls
the ECS service to determine if there are any tasks to run. An AWS Batch Job
will be turned into an ECS task. From there, the agent will pick it up and
manage it, sending updates back to ECS (and from there AWS Batch) on a regular
basis.

There are some limits that will impact the design of the AWS Batch backend, and
will be discussed later. These are:

* [AWS Batch Limits](https://docs.aws.amazon.com/batch/latest/userguide/service_limits.html)
* [8k container overrides limit.](https://docs.aws.amazon.com/cli/latest/reference/ecs/run-task.html)

The ECS workers used by the AWS Batch backend can be any instance type and should
be based on an AMI running the ECS agent and docker. An ECS optimized AMI is recommended.
An EC2 LaunchTemplate is used to provide some additional "on first boot" configuration that:
1. Installs AWS CLI v2,
1. Installs a script to mount an EBS as a `btrfs` file system that will auto-expand,
1. Configures docker to use that file system so that the "filesystem" of the container
will auto-expand,
1. Installs a `fetch_and_run.sh` script that allows the container to download 
generated shell scripts from S3 that contain the instructions of the workflow
task 

```text
                  +-------------+
                  |             |
                  |  AWS Batch  |
                  |             |
                  +------+------+
                         |
                         |
                         |
                         |
                         |
        +----------------v------------------+
        |                                   |
        |  Elastic Container Service (ECS)  |
        |                                   |
        +----------------+------------------+
                         |
                         |
                         |
                         |
                         |
+------------------------v-------------------------+
|                                                  |
|  AutoScaling Group                               |
|                                                  |
| +---------------------------------+              |
| |                                 |              |
| |  EC2 Instance                   |              |
| |                                 |              |
| |  +--------------------+         |              |
| |  |                    |         |              |
| |  |  Docker Container  |         |              |
| |  |                    |         |              |
| |  +--------------------+  ...    |              |
| |                                 |              |
| +---------------------------------+     ...      |
|                                                  |
+--------------------------------------------------+

```

Cromwell AWS Batch Backend
--------------------------

There are several scala classes as part of the AWS Batch Backend, but
the primary classes involved in running the backend are shown below. The
arrows represent the flow of job submission.

```text
    +----------------------------------------+
    |                                        |
    |  AwsBatchBackendLifecycleActorFactory  |
    |                                        |
    +------------------+---------------------+
                       |
                       |
                       |
                       |
                       |
    +------------------v----------------------+
    |                                         |
    |  AwsBatchAsyncBackendJobExecutionActor  |
    |                                         |
    +------------------+----------------------+
                       |
                       |
                       |
                       |
                       |
               +-------v-------+                 +-------------------------+
               |               |                 |                         |
               |  AwsBatchJob  +----------------->  AwsBatchJobDefinition  |
               |               |                 |                         |
               +---------------+                 +-------------------------+
```

1. The `AwsBatchBackendLifecycleActorFactory` class is configured by the user
   as the Cromwell backend. This factory provides an object from the
   `AwsBatchAsyncBackendJobExecutionActor` class to create and manage the job.
2. The `AwsBatchAsyncBackendJobExecutionActor` creates and manages the job.
   The job itself is encapsulated by the functionality in `AwsBatchJob`.
3. `AwsBatchJob` is the primary interface to AWS Batch. It creates the
   necessary `AwsBatchJobDefinition`, then submits the job using the SubmitJob
   API.
4. `AwsBatchJobDefinition` is responsible for the creation of the job definition.
   In AWS Batch, every job must have a definition. Note that the job definition
   can be overridden by the `SubmitJob`, so the `JobDefinition` contains core information such
   as the docker image type while the `SubmitJob` contains details that are more related to
   the actual task.

AWS Batch Job Instantiation
---------------------------
```text
             +--------------------+
             |                    |
             |  Cromwell Backend  |
             |                    |
             +---------+----------+
                       |
                       |
                   SubmitJob
                       |
                       |
                +------v------+
                |             |
                |  AWS Batch  |
                |             |
                +------^------+
                       |
                       |
                     Polls
                       |
                       |
                +------+------+
                |             |
                |  ECS Agent  |
                |             |
                +------+------+
                       |
           Creates, Launches and Monitors
                       |
              +--------v---------+ 
              |                  |
              |  Task Container  |
              |                  |
              +------------------+

```

When a Cromwell task begins, the Cromwell backend will call the SubmitJob
API of AWS Batch. From there, the backend will call the AWS Batch `DescribeJobs`
API to provide status to the Cromwell engine as requested.

Once the job is Submitted in AWS Batch, one of the EC2 instances assigned
to the compute environment (a.k.a. ECS Cluster) with a running agent will
pick up the Cromwell Job/AWS Batch Job/ECS Task and run it. Importantly,
AWS Batch calls ECS' `RunTask` API when submitting the job. It uses the
task definition, and overrides both the command text and the environment
variables.

Input files are read into the container from S3 and output files are copied back to
S3. Three additional files are also written to the S3 bucket using the names of these
environment variables:

* AWS_CROMWELL_RC_FILE (the return code of the task)
* AWS_CROMWELL_STDOUT_FILE (STDOUT of the task)
* AWS_CROMWELL_STDERR_FILE (STDERR of the task)

These files are placed in the correct location in S3 after task execution. In addition
STDOUT and STDERR are fed to the tasks cloudwatch log.

Input and Command Compression
-----------------------------

NOTE: All limits in this section are subject to change

In testing, specifically with large fan-in operations such as the MergeVCFs
task of the Haplotype caller test, that the container overrides length limit
of 8k was being exceeded. There are several limits described on AWS Batch,
and a limit for container overrides on ECS, all of which should be considered.

* Maximum payload size for RegisterJobDefinition calls: 24KiB
* Maximum payload size for SubmitJob calls: 30KiB
* Maximum JSON payload for ECS RunTask containerOverrides values: 8KiB

Effective limits, however, are much, much smaller. While both AWS Batch and
ECS have command and environment as part of their Job/Task definitions
respectively, AWS Batch passes both command and environment through to ECS
based solely on RunTask. While a lot of effort was initially placed on
balancing payloads between RegisterJobDefinition and SubmitJob, because
everything is passed as an override to RunTask, we're gated by the ECS
RunTask 8KiB limit.

Dependencies
------------

Two dependencies were added to Cromwell as part of this project, though
existing Cromwell dependencies were also leveraged. These two are:

* AWS Java SDK v2
* elerch/S3 Filesystem

The Java SDK version two carries with it a significant amount of additional
dependencies. These had a significant effect on the ability of Cromwell to
compile, but ultimately the overlapping dependencies came down to:

* com.fasterxml.jackson.core (jackson-annotations):
    Jackson version was bumped to accomodate throughput
* org.slf4j (jcl-over-slf4j):
    Ignored by aws sdk in favor of the version bundled in Cromwell
* nettyHandler:
    Version bumped to accomodate AWS SDK
* sttpV:
    Version bumped to accomodate AWS SDK

While the AWS SDK v2 has nio capabilities, it does not include
a FileSystemProvider to allow Cromwell to operate seamlessly with AWS S3.
This was found later in the project, and as such, the filesystem provider
is not integrated with the Cromwell configuration system. The S3 FS was
forked from Udacity's provider, which was based on AWS SDK for Java version 1.
Of particular note is the fact that all API configuration for the provider
is performed through environment variables. As such, configuration is possible,
but currently disjoint from Cromwell proper.

Authentication
--------------

Authentication is required for both the Cromwell AWS Backend and the S3
Filesystem provider. By default, in both cases, the default credential provider
chain is followed. AWS tools all follow the same prioritized set of
checks for access key/secret key (and optionally token), with the exception
of #2 below.

[Default credential provider chain](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html)
1. Configured permissions (environment)
2. (NOTE BELOW) Java properties
3. Access key/secret key as defined in $HOME/.aws/credentials and config
4. Container role
5. EC2 Instance role

NOTE: The Java properties check is specific and unique to the Java SDK. This
      does not apply to other SDKs or tools (e.g. the CLI).

Normally customers will be using an EC2 Instance role (recommended) or file configuration
as described in #3 (not recommended in production).

Permissions
-----------

Within AWS, everything must be authorized. This is a consistent rule, and as
such, AWS Services themselves are not immune to the rule. Therefore, customers
of AWS are responsible for granting services access to APIs within their account.
The flow described below represents the permissions needed by each stage, from 
Cromwell server through the task running. This includes the permissions needed for
the AWS Services involved in the processing of the work.

```text
+----------------------------+
|                            |  s3:GetObject on bucket for workflow and script bucket
|                            |  s3:ListObjects on script bucket
|                            |  s3:PutObject on script bucket
|          Cromwell          |  batch:RegisterTaskDefinition
|                            |  batch:SubmitJob
|                            |  batch:DescribeJobs
|                            |  batch:DescribeJobDefinitions
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |  AWSBatchServiceRole managed policy - described at:
|          AWS Batch         |
|                            |     https://docs.aws.amazon.com/batch/latest/userguide/service_IAM_role.html
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |  AWSServiceRoleForECS Service-linked role, documented at:
|                            |
| Elastic Container Service  |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/using-service-linked-roles.html
|                            |
|  (See discussion #1 below) |  AmazonEC2ContainerServiceAutoscaleRole managed policy - described at:
|                            |
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/autoscale_IAM_role.html
+-------------+--------------+
              |
              |
              |
+-------------v--------------+
|                            |
|                            |  AmazonEC2ContainerServiceforEC2Role managed policy, described at:
| ECS Agent (running on EC2) |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html (EC2)
|                            |    OR
|                            |  AmazonECSTaskExecutionRolePolicy managed policy, described at:  
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_execution_IAM_role.html (Fargate)
+-------------+--------------+ 
              |
              |
              |
+-------------v--------------+
|                            |  Task Role permissions. These are user defined, but ecs-tasks.amazon.com must have sts:AssumeRole trust relationship defined. Documentation:
|       Task Container       |
|                            |     https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_IAM_role.html
|                            |  s3:GetObject, s3:PutObject, s3:ListObjects
+----------------------------+
```


1. ECS has several sets of permissions for various items. AWS Batch, however,
   does not take advantage of certain features of ECS, most importantly
   ECS Services are out of scope of AWS Batch. ECS services require things
   like load balancing registration and DNS updates. While there is
   documentation regarding roles related to ECS services, these are
   irrelevant to the Cromwell use case.
2. Other than access to the main Cromwell bucket, the task container itself 
   does not need additional permissions unless
   the task in the WDL has been defined with a command that interfaces
   with AWS directly. This may include access to additional s3 buckets containing
   things like reference genome files.

   Task container permissions are currently supported through ECS and AWS
   Batch, but there is no configuration currently wired for the Cromwell
   AWS Backend to pass these settings to AWS. As such, the task container
   permissions must be managed by attaching a role to the EC2 Instance
   with permissions necessary for both the ECS Agent and the task container.

NOTE: ECS Agent permissions currently must use the permissions as outlined
      in the AmazonEC2ContainerServiceForEC2Role managed policy.

Future considerations
---------------------

AWS Batch Backend

* Should the 'disks' configuration be renamed to mount points or maybe
  ignored? This might make the wdl configuration more portable between backends
* The intent is to be able to override the queueArn between the
  default-runtime-attributes and the runtime attributes in the WDL. This
  appears broken at the moment.
* Job caching appears to be broken at the moment. Identical tasks need not be repeated
  if the results of a previous run of the task are still available.
* Retrying failed jobs is not currently attempted. Adding this would be beneficial 
  especially in conjunction with result caching
* Some S3 FS stuff can be removed: It looks like there’s a bunch of leftover
  unused code here from before having the nio implementation
  [here](https://github.com/broadinstitute/cromwell/tree/develop/filesystems/s3/src/main/scala/cromwell/filesystems/s3/batch)
  Also, I recommend adding another case for S3
  [here](https://github.com/broadinstitute/cromwell/blob/7a830a7fceaab9e7eaaa6802e58fe3cfd5d411a8/engine/src/main/scala/cromwell/engine/io/nio/NioFlow.scala#L95)
  otherwise files will be pulled down to be hashed instead of looking up the
  hash from S3 metadata. You might need to add the s3 Cromwell fs as a dependency
  to the engine to be able to do that.
* S3 Filesystem should be an official AWS thing (e.g. upplication or awslabs account)
* 8k container overrides limit should be bigger (ECS)
* We should understand why AWS Batch is putting everything in container overrides
* Authentication configuration consistency between S3FS and AWS Batch backend
* Full configuration of jobs in AWS Batch

Cromwell

* The style of integration with backends requires a lot of boilerplate code
  that I believe can be significantly reduced via heavier use of traits and
  supporting libraries
* There is a lot of data available within the backend if you know where
  to look, but if you don't, there is a lot of looking at inherited classes
  etc. to try to find them. Often, lack of code comments, etc and very similarly
  named variables can be confusing for folks confronting the code base for the
  first time. Workflow root vs call root is a great example. Adding additional
  comments and potentially context objects to encapsulate state would make
  the codebase more approachable.
* There is a significant amount of dependent libraries (with more added by
  the introduction of the S3 filesystem and the AWS SDK). Dependency management
  is challenging as a result. Adding significant new functionality is relatively
  painful when new dependencies are needed.
# Stackdriver task monitor*

This folder contains code for monitoring resource utilization in PAPIv2 tasks
through [Stackdriver Monitoring](https://cloud.google.com/monitoring).

[monitor.py](monitor.py) script
is intended to be used as a Docker image, via a background "monitoring action" in PAPIv2.
The image can be specified through `monitoring_image` workflow option.

It uses [psutil](https://psutil.readthedocs.io) to
continuously measure CPU, memory and disk space utilization
and disk IOPS, and periodically report them
as distinct metrics to Stackdriver Monitoring API.

The labels for each time point contain
- Cromwell-specific values, such as workflow ID, task call name, index and attempt.
- GCP instance values such as instance name, zone, number of CPU cores, total memory and disk size.

This approach enables:

1)  Users to easily plot real-time resource usage statistics across all tasks in
    a workflow, or for a single task call across many workflow runs,
    etc.

    This can be very powerful to quickly determine the outlier tasks
    that could use optimization, without the need for any configuration
    or code.

2)  Scripts to easily get aggregate statistics
    on resource utilization and to produce suggestions
    based on those.

[*] Detailed discussion: [PR 4510](https://github.com/broadinstitute/cromwell/pull/4510).
For information on Cromwell's backends, check out the [Cromwell Backend Documentation](http://cromwell.readthedocs.io/en/develop/backends/Backends/).# Cromwell Example Backends

This is a folder of example backend providers for Cromwell. You can read about
the providers here, and then copy paste one or more of the providers you want
to use to your Cromwell configuration file, represented here as the
[cromwell.examples.conf](cromwell.examples.conf) file in the base of the 
repository.

## What are the backend providers?

### Cloud Providers

 - [AWS](AWS.conf): Amazon Web Services ([documentation](https://cromwell.readthedocs.io/en/stable/tutorials/AwsBatch101/))
 - [BCS](BCS.conf) Alibaba Cloud Batch Compute (BCS) backend ([documentation](https://cromwell.readthedocs.io/en/stable/backends/BCS/))
 - [TES](TES.conf) is a backend that submits jobs to a server with protocol defined by GA4GH ([documentation](https://cromwell.readthedocs.io/en/stable/backends/TES/))
 - [PAPIv2](PAPIv2.conf): Google Pipelines API backend (version 2!) ([documentation](https://cromwell.readthedocs.io/en/stable/backends/Google/))

### Containers

 - [Docker](Docker.conf): an example backend that only runs workflows with docker in *every* command
 - [Singularity](singularity.conf): run Singularity containers locally ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#local-environments))
 - [Singularity+Slurm](singularity.slurm.conf): An example using Singularity with SLURM ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#job-schedulers))
 - [TESK](TESK.conf) is the same, but intended for Kubernetes. See the [TES docs](https://cromwell.readthedocs.io/en/stable/backends/TES/) at the bottom.
 - [udocker](udocker.conf): to interact with udocker locally [documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#udocker)
 - [udocker+Slurm](udocker.slurm.conf): to interact with udocker on SLURM ([documentation](https://cromwell.readthedocs.io/en/develop/tutorials/Containers/#udocker))

### Workflow Managers

 - [HtCondor](HtCondor.conf): a workload manager at UW-Madison ([documentation](https://cromwell.readthedocs.io/en/stable/backends/HTcondor/))
 - [LSF](LSF.conf): the Platform Load Sharing Facility backend ([documentation](https://cromwell.readthedocs.io/en/stable/backends/LSF/))
 - [SGE](SGE.conf): a backend for Sungrid Engine ([documentation](https://cromwell.readthedocs.io/en/stable/backends/SGE))
 - [slurm](slurm.conf): SLURM workload manager ([documentation](https://cromwell.readthedocs.io/en/stable/backends/SLURM/))

### Custom

 - [LocalExample](LocalExample.conf): What you should use if you want to define a new backend provider ([documentation](https://cromwell.readthedocs.io/en/stable/backends/Local/))


## How do I add a backend provider?

The section in the file called "backends" has a key, "providers" that looks like
this:

```

backend {

  # Override the default backend.
  #default = "LocalExample"

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
        ....
  }

}
```

The examples here also have this section. You would want to copy paste the content
of the file, specifically the section for the provider under backend -> providers,
into the backend -> providers section in the [cromwell.examples.conf](cromwell.examples.conf).
Here is what it would look like to add the [slurm](slurm.conf) backend
provider example. 

```
backend {

  # Override the default backend.
  #default = "LocalExample"

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
    slurm {
         ...
      }
    }

    # Second backend provider would be copy pasted here!

  }
}
```

This isn't json, so you don't need to add commas between the providers - just
copy paste them one after the other in the backend -> providers section.
Let's say we wanted slurm to be our default! We would do this:

```
backend {

  # Override the default backend.
  default = slurm

  # The list of providers. Copy paste the contents of a backend provider in this section
  providers {
    slurm {
         ...
      }
    }
  }
}
```

Don't forget to customize the sections for your purposes! If anything is
not explained clearly, please [open an issue](https://github.com/broadinstitute/cromwell/issues).

## What if a provider is missing?

If a provider is missing and you don't want to use the [LocalExample](LocalExample.conf)
to write a custom provider, please [let us know](https://github.com/broadinstitute/cromwell/issues)
and we can start discussion about how to define your backend.
###
### IMPORTANT: Please file new issues over in our Jira issue tracker!
###
### https://broadworkbench.atlassian.net/projects/BA/issues
###
### You may need to create an account before you can view/create issues.
###

<!--
Hi!  Thanks for taking the time to report feedback.

Before posting an issue over in Jira tracker, please check whether your question is already answered in our:
  forum  https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team
  documentation http://cromwell.readthedocs.io/en/develop/

Other forums:
FireCloud https://gatkforums.broadinstitute.org/firecloud/categories/ask-the-firecloud-team
WDL https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team
CWL https://www.biostars.org/
-->

<!-- Are you seeing something that looks like a bug?  Then great!  You're almost in the right place.  -->

<!-- You'll want to go to https://broadworkbench.atlassian.net/projects/BA/issues and then tell us: -->

<!-- Which backend are you running?  -->

<!-- Paste/Attach your workflow if possible: -->

<!-- Paste your configuration if possible, MAKE SURE TO OMIT PASSWORDS, TOKENS AND OTHER SENSITIVE MATERIAL: -->
An **Amazon AWS S3** FileSystem Provider **JSR-203** for Java 8 (NIO2) using the AWS SDK v2.

This version uses the released version of the AWS SDK v2 and is forked from
[Emil Lerch's S3 FileSystem](https://github.com/elerch/Amazon-S3-FileSystem-NIO2/tree/e8283b5).
That version uses a preview version of the AWS SDK v2 and is forked from
[Upplication's S3 FileSystem](https://github.com/upplication/Amazon-S3-FileSystem-NIO2).

The previous fork did not update the tests for AWS SDK v2, so the tests targeting AWS SDK v1 have been left out of this
fork.
For information on Cromwell's Integration Testing Suite, see the [Cromwell documentation on Centaur](https://cromwell.readthedocs.io/en/develop/developers/Centaur/).# Cromwell Docker Development

This is a container intended to help with development of Cromwell, in the
case that you don't want to install the dependencies on your host. It
includes the required software mentioned in [the developer docs](http://cromwell.readthedocs.io/en/develop/Building/) but with additional instructions to interact with the repository.

As with the build instructions, first start by cloning the Cromwell repository from GitHub:

```bash
$ git clone git@github.com:broadinstitute/cromwell.git
```

And change into this directory:

```bash
$ cd cromwell/scripts/docker-develop
```

This is where this README.md sits with a [Dockerfile](Dockerfile)

## Step 1. Build the Container
From this folder with the Dockerfile, build the container. In the command below
we are calling it `cromwell-dev` and you can choose to change this name if you want.

```bash
$ docker build -t cromwell-dev .
```

and change back to the root of the repository

```bash
cd ../../
```

## Step 2. Shell into Working Environment
If you require a specific version of Cromwell as a starting point, do the appropriate `git checkout` now. 

You are next going to want to bind the cromwell
source code to be somewhere in the container. Actually, we have a `/code` directory
ready for you to make this easy. Run this command from the base of the repository:

```bash
$ docker run -v $PWD/:/code -it cromwell-dev bash
```

## Step 4. The Development Steps

You are going to use the command `sbt assembly` to create a runnable Cromwell JAR. The output
will be `server/target`. Remember, here we are inside the container with scala and sbt:


```bash
$ sbt assembly
```

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

You can then interact with it in the container, or on your host if you like. Remember that
you have Java already in the container, so it makes sense to develop there.

Looking for the MySQL docs? Checkout the [Cromwell documentation](https://cromwell.readthedocs.io/en/develop/Configuring/#database)!For information on the WOMtool, see the [Cromwell documentation on WOMtool](https://cromwell.readthedocs.io/en/develop/WOMtool).# Customize tasks

Runtime attributes can be specified in one of two ways:

 1. Within a task you can specify runtime attributes to customize the environment for the call.  
 2. [Default runtime attributes](#default-values) for all tasks can be specified in [Workflow Options](wf_options/Overview.md).

_Task Example_

```
task jes_task {
  command {
    echo "Hello JES!"
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "4G"
    cpu: "3"
    zones: "us-central1-c us-central1-b"
    disks: "/mnt/mnt1 3 SSD, /mnt/mnt2 500 HDD"
  }
}
workflow jes_workflow {
  call jes_task
}
```


## Recognized Runtime attributes and Backends

Cromwell recognizes certain runtime attributes and has the ability to format these for some [Backends](/backends/Backends). See the table below for common attributes that apply to _most_ backends.

| Runtime Attribute    | LOCAL |  Google Cloud  | AWS Batch |  HPC  |
| -------------------- |:-----:|:-----:|:-----:|:------:|
| [cpu](#cpu)                  |       |   x   |   x   |  `cpu`  |
| [memory](#memory)                              |       |   x   |   x   |  `memory_mb` / `memory_gb`  |
| [disks](#disks)                                |       |   x   |       |  *  |
| [docker](#docker)                              |   x   |   x   |   x   |  `docker` (see below)  |
| [maxRetries](#maxretries)                      |   x   |   x   |   x   | * |
| [continueOnReturnCode](#continueonreturncode) |   x   |   x   |   x   | * |
| [failOnStderr](#failonstderr)                  |   x   |   x   |   x   |  *  |


> `*` The HPC [Shared Filesystem backend](/backends/HPC#shared-filesystem) (SFS) is fully configurable and any number of attributes can be exposed. Cromwell recognizes some of these attributes (`cpu`, `memory` and `docker`) and parses them into the attribute listed in the table which can be used within the HPC backend configuration.


### Google Cloud Specific Attributes
There are a number of additional runtime attributes that apply to the Google Cloud Platform:

- [zones](#zones)
- [preemptible](#preemptible)
- [bootDiskSizeGb](#bootdisksizegb)
- [noAddress](#noaddress)
- [gpuCount, gpuType, and nvidiaDriverVersion](#gpucount-gputype-and-nvidiadriverversion)
- [cpuPlatform](#cpuplatform)
- [useDockerImageCache](#usedockerimagecache)



## Expression support

Runtime attribute values are interpreted as expressions.  This means that it has the ability to express the value of a runtime attribute as a function of one of the task's inputs.  
_For example:_

```
task runtime_test {
  String ubuntu_tag
  Int memory_gb

  command {
    ./my_binary
  }

  runtime {
    docker: "ubuntu:" + ubuntu_tag
    memory: memory_gb + "GB"
  }
}
```

HPC backends may define other configurable runtime attributes beyond the five listed, to find out more visit the [SunGridEngine](/backends/SGE) tutorial.

## Default Values

Default values for runtime attributes can be specified via [Workflow Options](wf_options/overview).  
For example, consider this WDL file:

```wdl
task first {
  command { ... }
}

task second {
  command {...}
  runtime {
    docker: "my_docker_image"
  }
}

workflow w {
  call first
  call second
}
```

And this set of workflow options:

```json
{
  "default_runtime_attributes": {
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
}
```

Then, these values for `docker` and `zones` will be used for any task that does not explicitly override them in the WDL file. In return, the effective runtime for `task first` is:

```
{
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
```

And the effective runtime for `task second` is:

```
{
    "docker": "my_docker_image",
    "zones": "us-central1-c us-central1-b"
  }
```

Note how for `task second` the WDL value for `docker` is used instead of the default provided in the workflow options.


## Runtime Attribute Descriptions

### `cpu`

*Default: _1_*

The `cpu` runtime attribute represents the number of cores that a job requires, however each backend may interpret this differently:

- In Google Cloud: this is interpreted as "the minimum number of cores to use."
- In HPCs (SFS): this is configurable, but usually a reservation and/or limit of number of cores.

Example
```
runtime {
  cpu: 2
}
```


#### CWL

CWL splits the `cpu` requirement into `cpuMin` and `cpuMax`. If one of them is provided, `cpu` will inherit this value. If both of them are provided, `cpu` will take the value of `cpuMin`.
If none is provided, `cpu` will default to its default value.

Note: If provided, `cpuMin` and/or `cpuMax` will be available to the [HPC runtime attribute configuration](/tutorials/HPCIntro.md#specifying-the-runtime-attributes-for-your-hpc-tasks).


### `memory`
*Default: "2G"*

Memory is the amount of RAM that should be allocated to a task, however each backend may interpret this differently:

- Google Cloud: The minimum amount of RAM to use.
- SFS: Configurable, but usually a reservation and/or limit of memory.

The memory size is specified as an amount and units of memory, for example "4G":

```
runtime {
  memory: "4G"
}
```

Within the SFS backend, you can additionally specify `memory_mb` or `memory_gb` as runtime attributes within the configuration. More information can be found [here](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/#specifying-the-runtime-attributes-for-your-hpc-tasks).

#### CWL

CWL splits the `memory` requirement into `ramMin` and `ramMax`. If one of them is provided, `memory` will inherit this value. If both of them are provided, `memory` will take the value of `ramMin`. If none is provided, `memory` will default to its default value.

Note: If provided, `ramMin` and/or `ramMax` will be available to the [HPC runtime attribute configuration](/tutorials/HPCIntro.md#specifying-the-runtime-attributes-for-your-hpc-tasks).


### `disks`

This attribute specifies volumes that will be mounted to the VM for your job. These volumes are where you can read and write files that will be used by the commands within your workflow. 


They are specified as a comma separated list of disks. Each disk is further separated as a space separated triplet (e.g. `local-disk 10 SSD`) consisting of:

1. Mount point (absolute path), or `local-disk` to reference the mount point where Google Cloud will localize files and the task's current working directory will be
2. Disk size in GB (rounded to the next 375 GB for LOCAL)
3. Disk type.  One of: "LOCAL", "SSD", or "HDD" ([documentation](https://cloud.google.com/compute/docs/disks/#overview))

All tasks launched on Google Cloud *must* have a `local-disk`.  If one is not specified in the runtime section of the task, then a default of `local-disk 10 SSD` will be used.  The `local-disk` will be mounted to `/cromwell_root`.

For the AWS Batch backend, the disk volume is managed by AWS EBS with autoscaling capabilities.  As such, the Disk size and disk type will be ignored. If provided, the mount point will be verified at runtime.


The Disk type must be one of "LOCAL", "SSD", or "HDD". When set to "LOCAL", the size of the drive is constrained to 375 GB intervals so intermediate values will be rounded up to the next 375 GB. All disks are set to auto-delete after the job completes.

*Example 1: Changing the Localization Disk*

```
runtime {
  disks: "local-disk 100 SSD"
}
```

*Example 2: Mounting an Additional Two Disks*

```
runtime {
  disks: "/mnt/my_mnt 3 SSD, /mnt/my_mnt2 500 HDD"
}
```

### `docker`

When specified, Cromwell will run your task within the specified Docker image. 

```
runtime {
  docker: "ubuntu:latest"
}
```

- Local: Cromwell will automatically run the docker container.
- SFS: When a docker container exists within a task, the `submit-docker` method is called. See the [Getting started with containers](/tutorials/Containers/) guide for more information.
- GCP: This attribute is mandatory when submitting tasks to Google Cloud.
- AWS Batch: This attribute is mandatory when submitting tasks to AWS Batch.


### `maxRetries`

*Default: _0_*

This retry option is introduced to provide a method for tackling transient job failures. For example, if a task fails due to a timeout from accessing an external service, then this option helps re-run the failed the task without having to re-run the entire workflow. It takes an Int as a value that indicates the maximum number of times Cromwell should retry a failed task. This retry is applied towards jobs that fail while executing the task command. This method only applies to transient job failures and is a feeble attempt to retry a job, that is it cannot be used to increase memory in out-of-memory situations.

If using the Google backend, it's important to note that The `maxRetries` count is independent from the [preemptible](#preemptible) count. For example, the task below can be retried up to 6 times if it's preempted 3 times AND the command execution fails 3 times.

```
runtime {
  preemptible: 3
  maxRetries: 3
}
```

### `continueOnReturnCode`
*Default: _0_*

When each task finishes it returns a code. Normally, a non-zero return code indicates a failure. However you can override this behavior by specifying the `continueOnReturnCode` attribute.

When set to false, any non-zero return code will be considered a failure. When set to true, all return codes will be considered successful.

```
runtime {
  continueOnReturnCode: true
}
```

When set to an integer, or an array of integers, only those integers will be considered as successful return codes.

```
runtime {
  continueOnReturnCode: 1
}
```

```
runtime {
  continueOnReturnCode: [0, 1]
}
```

### `failOnStderr`

*Default: _false_*

Some programs write to the standard error stream when there is an error, but still return a zero exit code. Set `failOnStderr` to true for these tasks, and it will be considered a failure if anything is written to the standard error stream.

```
runtime {
  failOnStderr: true
}
```



### `zones`

The ordered list of zone preference (see [Region and Zones](https://cloud.google.com/compute/docs/zones) documentation for specifics).

*The zones are specified as a space separated list, with no commas:*

```
runtime {
  zones: "us-central1-c us-central1-b"
}
```

Defaults to the configuration setting `genomics.default-zones` in the Google Cloud configuration block, which in turn defaults to using `us-central1-b`.

### `preemptible`

*Default: _0_*

Passed to Google Cloud: "If applicable, preemptible machines may be used for the run."

Take an Int as a value that indicates the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one.  
*eg. With a value of 1, Cromwell will request a preemptible VM, if the VM is preempted, the task will be retried with a non-preemptible VM.*

```
runtime {
  preemptible: 1
}
```




### `bootDiskSizeGb`

In addition to working disks, Google Cloud allows specification of a boot disk size. This is the disk where the docker image itself is booted (**not the working directory of your task on the VM**).
Its primary purpose is to ensure that larger docker images can fit on the boot disk.
```
runtime {
  # Yikes, we have a big OS in this docker image! Allow 50GB to hold it:
  bootDiskSizeGb: 50
}
```

Since no `local-disk` entry is specified, Cromwell will automatically add `local-disk 10 SSD` to this list.


### `noAddress`

This runtime attribute adds support to disable assigning external IP addresses to VMs provisioned by the Google backend. If set to true, the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs.

Note well!  You must enable "Private Google Access" for this feature to work. See "How To Setup" below.

For example, the task below will succeed:
```
command {
  echo "hello!"
  
}

runtime {
  docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
  noAddress: true
}
```

The task below will fail for two reasons:
 1. The command is accessing an external service, in this case GitHub.
 2. The docker image is available in DockerHub and not the Google Container Registry. 
```
command {
  git clone https://github.com/broadinstitute/cromwell.git
  
}

runtime {
  docker: "docker.io/alpine/git:latest"
  noAddress: true
}
```

#### How to Setup

Configure your Google network to use "Private Google Access". This will allow your VMs to access Google Services including Google Container Registry, as well as Dockerhub images.

1. Using `gcloud compute networks subnets list`, identify the subnet and region you will be using with Cromwell. If multiple, run the next step for each region and subnet you wish to use.
1. `gcloud compute networks subnets update [SUBNET-NAME] --region [REGION]  --enable-private-ip-google-access`

That's it!  You can now run with `noAddress` runtime attribute and it will work as expected.

### `gpuCount`, `gpuType`, and `nvidiaDriverVersion`

Attach GPUs to the instance when running on the Pipelines API([GPU documentation](https://cloud.google.com/compute/docs/gpus/)).
Make sure to choose a zone for which the type of GPU you want to attach is available.

The types of compute GPU supported are:

* `nvidia-tesla-k80` 
* `nvidia-tesla-v100`
* `nvidia-tesla-p100`
* `nvidia-tesla-p4`
* `nvidia-tesla-t4`

For the latest list of supported GPU's, please visit [Google's GPU documentation](nvidia-drivers-us-public).

The default driver is `418.87.00`, you may specify your own via the `nvidiaDriverVersion` key.  Make sure that driver exists in the `nvidia-drivers-us-public` beforehand, per the [Google Pipelines API documentation](https://cloud.google.com/genomics/reference/rest/Shared.Types/Metadata#VirtualMachine). 

```
runtime {
    gpuType: "nvidia-tesla-k80"
    gpuCount: 2
    nvidiaDriverVersion: "418.87.00"
    zones: ["us-central1-c"]
}
```

### `cpuPlatform`

This option is specific to the Google Cloud backend, specifically [this](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform) feature when a certain minimum CPU platform is desired.

A usage example:

```
runtime {
    cpu: 2
    cpuPlatform: "Intel Cascade Lake"
}
```
Note that when this options is specified, make sure the requested CPU platform is [available](https://cloud.google.com/compute/docs/regions-zones/#available) in the `zones` you selected.

The following CPU platforms are currently supported by the Google Cloud backend:
- `Intel Cascade Lake`
- `Intel Skylake`     
- `Intel Broadwell`   
- `Intel Haswell`     
- `Intel Ivy Bridge`  
- `Intel Sandy Bridge`
- `AMD Rome`

### 'useDockerImageCache'

This option is specific to the Google Cloud backend, moreover it is only supported by Google Life Sciences API starting from version v2 beta.
In order to use this feature Cromwell has to have PAPI v2 backend configured with this feature enabled.  
More information about this feature and it's configuration can be found [in the Google backend section of documentation](backends/Google.md).
Command line utilities for interacting with the Workflow Object Model (WOM). You can download the latest WOMtool from the [Cromwell releases page on Github](https://github.com/broadinstitute/cromwell/releases/latest).

## Requirements

The following is the toolchain used for development of womtool.  Other versions may work, but these are recommended.

* [Scala 2.12](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/)
* [Git](https://git-scm.com/)

## Building

`sbt assembly` will build a runnable JAR in `womtool/target/scala-2.12/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

## Command Line Usage

Run the JAR file with no arguments to get the usage message:

`$ java -jar womtool.jar`

```bash
java -jar /path/to/womtool.jar <action> <parameters>

Actions:
validate [--list-dependencies] <WDL file>

  Performs full validation of the WDL file including syntax
  and semantic checking. -l or --list-dependencies is an optional flag to 
  list files referenced in import statements.

inputs <WDL file>

  Print a JSON skeleton file of the inputs needed for this
  workflow.  Fill in the values in this JSON document and
  pass it in to the 'run' subcommand.

highlight <WDL file> <html|console>

  Reformats and colorizes/tags a WDL file. The second
  parameter is the output type.  "html" will output the WDL
  file with <span> tags around elements.  "console" mode
  will output colorized text to the terminal
  
parse <WDL file>

  Compares a WDL file against the grammar and writes out an
  abstract syntax tree if it is valid, and a syntax error
  otherwise.  Note that higher-level AST checks are not done
  via this sub-command and the 'validate' subcommand should
  be used for full validation

graph <WDL file>

  Reads a WDL file against the grammar and prints out a
  .dot of the DAG if it is valid, and a syntax error
  otherwise. Note that graph currently DOES NOT WORK on
  version 1.0 workflows.

womgraph <WDL or CWL file> [ancillary files]

  Reads a WDL or CWL file from the first argument and
  converts it to a WOM representation then prints out a graph
  of the WOM produced.
  Any imported files can be supplied as subsequent arguments.
```

### `validate`

Given a WDL file, this runs the full syntax checker over the file and resolves imports in the process.  If any syntax errors are found, they are written out.  Otherwise the program exits.

Error if a `call` references a task that doesn't exist:

`$ java -jar womtool.jar validate 2.wdl`

```
ERROR: Call references a task (BADps) that doesn't exist (line 22, col 8)

  call BADps
       ^
```

Error if namespace and task have the same name:

`$ java -jar womtool.jar validate 5.wdl`

```
ERROR: Task and namespace have the same name:

Task defined here (line 3, col 6):

task ps {
     ^

Import statement defined here (line 1, col 20):

import "ps.wdl" as ps
                   ^
```

##### --list-dependencies or -l flag

For a successful validation, this will output the list of files referenced in import statements in workflows and their subworkflows.

`$ java -jar womtool.jar validate -l myWdl.wdl`

```hocon
Success!
List of Workflow dependencies are:
/path/to/my/import/myImport.wdl
/path/to/another/import/anotherImport.wdl
https://path-to-http-import/httpImport.wdl
```

### `inputs`

Examine a WDL file with one workflow in it, compute all the inputs needed for that workflow and output a JSON template that the user can fill in with values.  The keys in this document should remain unchanged.  The values tell you what type the parameter is expecting.  For example, if the value were `Array[String]`, then it's expecting a JSON array of JSON strings, like this: `["string1", "string2", "string3"]`

`$ java -jar womtool.jar inputs 3step.wdl`


```
{
  "three_step.cgrep.pattern": "String"
}
```

This inputs document is used as input to the `run` subcommand.

### `highlight`

Formats a WDL file and semantically tags it.  This takes a second parameter (`html` or `console`) which determines what the output format will be.

test.wdl
```
task abc {
  String in
  command {
    echo ${in}
  }
  output {
    String out = read_string(stdout())
  }
}

workflow wf {
  call abc
}
```

### `parse`

Given a WDL file input, this does grammar level syntax checks and writes out the resulting abstract syntax tree.

`$ echo "workflow wf {}" | java -jar womtool.jar parse /dev/stdin`


```
(Document:
  imports=[],
  definitions=[
    (Workflow:
      name=<stdin:1:10 identifier "d2Y=">,
      body=[]
    )
  ]
)
```

This WDL file can be formatted in HTML as follows:

`$ java -jar womtool.jar highlight test.wdl html`

```
<span class="keyword">task</span> <span class="name">abc</span> {
  <span class="type">String</span> <span class="variable">in</span>
  <span class="section">command</span> {
    <span class="command">echo ${in}</span>
  }
  <span class="section">output</span> {
    <span class="type">String</span> <span class="variable">out</span> = <span class="function">read_string</span>(<span class="function">stdout</span>())
  }
}

<span class="keyword">workflow</span> <span class="name">wf</span> {
  <span class="keyword">call</span> <span class="name">abc</span>
}
```

### `graph`
 
The syntax of the graph command is:
```
womtool graph [--all] wdlFile.wdl
```

Given a WDL file input, command generates the data-flow graph through the system in `.dot` format.

For example the fork-join WDL:
```
task mkFile {
  command {
    for i in `seq 1 1000`
    do
      echo $i
    done
  }
  output {
    File numbers = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow forkjoin {
  call mkFile
  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }
  call join { input: wcCount = wc.count, grepCount = grep.count }
  output {
    join.proportion
  }
}
```

Produces the DAG:
```
digraph forkjoin {
  "call forkjoin.mkFile" -> "call forkjoin.wc"
  "call forkjoin.mkFile" -> "call forkjoin.grep"
  "call forkjoin.wc" -> "call forkjoin.join"
  "call forkjoin.grep" -> "call forkjoin.join"
}
```

#### The `--all` flag

If this flag is set, all WDL graph nodes become nodes in the generated DAG, even if they are not "executed". Typically this will mean task declarations and call outputs. 
For example in the above example, with `--all` you would get:

```
digraph forkjoin {
  "call forkjoin.grep" -> "String forkjoin.grep.pattern"
  "call forkjoin.grep" -> "output { forkjoin.grep.count = read_int(stdout()) }"
  "call forkjoin.grep" -> "File forkjoin.grep.in_file"
  "call forkjoin.wc" -> "output { forkjoin.wc.count = read_int(stdout()) }"
  "call forkjoin.grep" -> "call forkjoin.join"
  "call forkjoin.wc" -> "File forkjoin.wc.in_file"
  "call forkjoin.mkFile" -> "call forkjoin.grep"
  "call forkjoin.join" -> "output { forkjoin.join.proportion = read_int(stdout()) }"
  "call forkjoin.join" -> "Int forkjoin.join.wcCount"
  "call forkjoin.wc" -> "call forkjoin.join"
  "call forkjoin.mkFile" -> "output { forkjoin.mkFile.numbers = stdout() }"
  "call forkjoin.mkFile" -> "call forkjoin.wc"
  "call forkjoin.join" -> "Int forkjoin.join.grepCount"
}
```
## Overview

You can configure Cromwell settings either through configuration files or the Java command line.

Check out the tutorial on [How to Configure Cromwell](tutorials/ConfigurationFiles) for more information.

### Configuration examples

You can find a description of options and example stanzas in the [Cromwell Example Configuration][cromwell-examples-conf],
along with backend provider examples in the [Example Providers Folder][cromwell-examples-folder].

### Custom configuration files

You write configuration files in
[HOCON](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation).

To run using your configuration file, you should copy relevant stanzas from `cromwell.examples.conf` into a new
file, modify it as appropriate, then pass it to Cromwell via:

```
$ java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar ...
``` 

To create your own configuration file, start by creating a new text file, for example `my.conf`.

At the start of your file, include the file `application.conf` at the top before your custom configurations.

```hocon
# include the application.conf at the top
include required(classpath("application"))
```

From there, copy or add other configuration values and/or stanzas with your customizations.

```hocon
# include the application.conf at the top
include required(classpath("application"))

# Add customizations
webservice.port = 58000
```

Your configuration file can specify configuration as JSON-like stanzas or as dot-separated values. These next two examples are are equivalent.

_JSON-like stanza:_

```hocon
include required(classpath("application"))
webservice {
  port = 8000
  interface = 0.0.0.0
}
```

_Dot-separated values:_

```hocon
include required(classpath("application"))
webservice.port = 8000
webservice.interface = 0.0.0.0
```

## Configuration via command line

In addition to using configuration files, you can use dot-separated configuration names to specify values directly on the Java command line:

```
$ java -Dwebservice.port=8080 cromwell.jar ...
```

## Advanced configuration

**WARNING:** These advanced configuration values can significantly affect the performance of Cromwell. 

### Server

By default the Cromwell server will bind to `0.0.0.0` on port `8000`.  
You can then access it through a browser at `http://localhost:8000`.  
To change these settings, simply edit the following values in your configuration file:

```
webservice {
  port = 9000
  interface = 0.0.0.0
}
```

The above configuration will use port `9000`.

Cromwell uses `akka-http` to serve requests. For more advanced configuration settings, refer to the [akka-http](https://doc.akka.io/docs/akka-http/current/scala/http/configuration.html) documentation.

For example, to increase the request timeout to 30 seconds you can add this stanza to your configuration file:

```
akka.http.server.request-timeout = 30s
```

### I/O

**I/O Throttling**

Certain [backends](backends/Backends) impose I/O limits. For example the Pipelines API imposes a quota on the number of queries that can be made per second.

You can effectively control and throttle the number of requests and resources allocated to those operations in the `system.io` configuration:

```
system.io {
  number-of-requests = 100000
  per = 100 seconds
}
```

**I/O Resilience**

I/O operations can fail for a number of reason from network failures to server errors. Some of those errors are not fatal and can be retried.

Cromwell will retry I/O operations on such retryable errors, for a limited number of times before giving up and failing. This number (more precisely the number of attempts that will be made) can be set using the following configuration option:

```
system.io {
  number-of-attempts = 5
}
```

### Workflows

**Max Concurrent Workflows**

Cromwell has a configurable cap on the number of workflows running at a time. You can adjust the limit from the default `5000` by setting:

```hocon
system.max-concurrent-workflows = 5000
```

**New Workflow Poll Rate**

Cromwell will look for new workflows to start on a regular interval, configured as a number of seconds. You can change the polling rate from the default `2` seconds by editing the value:

```hocon
system.new-workflow-poll-rate = 2
```

**Max Workflow Launch Count**

On every poll, Cromwell will take at limited number of new submissions, provided there are new workflows to launch and the `system.max-concurrent-workflows` number has not been reached. While the default is to launch up to `50` workflows, you can override this by setting:

```hocon
system.max-workflow-launch-count = 50
```

***Abort configuration***

Cromwell will scan for abort requests using default configuration values equivalent to those below. In most circumstances
there shouldn't be a need to override these defaults.

```hocon
system {
  abort {
    # How frequently Cromwell should scan for aborts.
    scan-frequency: 30 seconds

    # The cache of in-progress aborts. Cromwell will add entries to this cache once a WorkflowActor has been messaged to abort.
    # If on the next scan an 'Aborting' status is found for a workflow that has an entry in this cache, Cromwell will not ask
    # the associated WorkflowActor to abort again.
    cache {
      # Guava cache concurrency.
      concurrency: 1
      # How long entries in the cache should live from the time they are added to the cache.
      ttl: 20 minutes
      # Maximum number of entries in the cache.
      size: 100000
    }
  }
}
```

### Database

**Using a MySQL Database**

Cromwell tracks the execution of workflows and stores outputs of task invocations in a SQL database. Cromwell supports either an external MySQL database, or a temporary in-memory database.

By default, Cromwell uses an in-memory database which will only live for the duration of the JVM.  This provides a quick way to run workflows locally without having to set up MySQL, though it also makes workflow executions somewhat transient.

To configure Cromwell to instead point to a MySQL database, first create the empty database.  In the example below, the database name is `cromwell`.

Then, edit your configuration file `database` stanza, as follows:

```hocon
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

To see the full list of possible parameters and values for the `db` stanza see [the slick documentation](http://slick.lightbend.com/doc/3.2.0/api/index.html#slick.jdbc.JdbcBackend$DatabaseFactoryDef@forConfig(String,Config,Driver):Database).

**Cromwell server on MySQL Database**

You can use [docker-compose](https://github.com/broadinstitute/cromwell/tree/develop/scripts) to link together a Cromwell docker image (built locally with `sbt docker` or available on [Dockerhub](https://hub.docker.com/r/broadinstitute/cromwell/)) with a MySQL docker image.

To change the version of Cromwell used, [change the tag in `compose/cromwell/Dockerfile`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/compose/cromwell/Dockerfile).

**Local**

`docker-compose up` from this directory will start a Cromwell server running on a MySQL instance with local backend.

The default configuration file used can be [found at `compose/cromwell/app-config/application.conf`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/compose/cromwell/app-config/application.conf).
To override it, simply mount a volume containing your custom `application.conf` to `/app-config` ([see `jes-cromwell/docker-compose.yml` for an example](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/jes-cromwell/docker-compose.yml)).

**Google Cloud**

The [`jes-cromwell` directory](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-compose-mysql/jes-cromwell) is an example of how to customize the original compose file with a configuration file and environment variables.

It uses the application default credentials of the host machine. To use it make sure your gcloud is up to date and that your [application-default credentials](https://developers.google.com/identity/protocols/application-default-credentials) are set up.

Then run `docker-compose -f docker-compose.yml -f jes-cromwell/docker-compose.yml up` to start a Cromwell server with a Google Cloud backend on MySQL.

**MySQL**

The data directory in the MySQL container is [mounted to `compose/mysql/data`](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-compose-mysql/compose/mysql/init), which allows the data to survive a `docker-compose down`.

To disable this feature, simply remove the `./compose/mysql/data:/var/lib/mysql` line in the [volume section of `docker-compose.yml`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/docker-compose.yml).

Note that in such case, the data will still be preserved by a `docker-compose stop` that stops the container but doesn't delete it.

**Notes**

To run Cromwell in the background, add `-d` at the end of the command:
`docker-compose up -d`.

To then see the logs for a specific service, run `docker-compose logs -f <service>`. 
For example `docker-compose logs -f cromwell`.

For more information about docker compose: [Docker compose doc](https://docs.docker.com/compose/).

**Insert Batch Size**

Cromwell queues up and then inserts batches of records into the database for increased performance. You can adjust the
number of database rows batch inserted by Cromwell as follows:

```hocon
database {
  insert-batch-size = 2000
}
```

**Separate Metadata Database**

This feature should be considered _experimental_ and likely to change in the future.

Cromwell stores metadata about each job and workflow intended. This metadata is intended for end users, and includes paths to job results, start and end times, etc. The metadata grows at a significantly faster rate than the rest of the internal engine data.

To use a separate database for metadata, under the `database` config section, configure a sub-path for `metadata` with custom settings.

```hocon
database {
  # Store metadata in a file on disk that can grow much larger than RAM limits.
  metadata {
    profile = "slick.jdbc.HsqldbProfile$"
    db {
      driver = "org.hsqldb.jdbcDriver"
      url = "jdbc:hsqldb:file:metadata-db-file-path;shutdown=false;hsqldb.tx=mvcc"
      connectionTimeout = 3000
    }
  }
}
```

If no override is found for `metadata`, Cromwell falls back to using the settings under the root `database` configuration.

**Database Time Zones**

Cromwell's default configuration assumes that its MySQL database is set to UTC.

The following MySQL configurations typically default to UTC and work with Cromwell out of the box:
- Google CloudSQL
- An official MySQL image running in Docker

These configurations may use the system, or local, time zone instead:
- MySQL installed natively on a workstation or server

If Cromwell fails to start with a message like
```
The server time zone value 'XXX' is unrecognized or represents more than one time zone.
```
you can resolve the problem by adding the option `&serverTimezone=UTC` to your database connection URL:
```hocon
url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true&serverTimezone=UTC"
```

Using this option does not alter your database's underlying timezone; rather, it causes Cromwell to "speak UTC" when communicating with the DB, and the DB server performs the conversion for you. 

**Using Cromwell with Postgresql**

To use Postgresql as the database, you will need to install and enable the
Large Object extension.  If the extension is present, setting up the database
requires just these commands:

```
$ createdb cromwell
$ psql -d cromwell -c "create extension lo;"
```

Postgresql configuration in Cromwell is very similar to MySQL.  An example:

```hocon
database {
  profile = "slick.jdbc.PostgresProfile$"
  db {
    driver = "org.postgresql.Driver"
    url = "jdbc:postgresql://localhost:5432/cromwell"
    user = "user"
    password = "pass"
    port = 5432
    connectionTimeout = 5000
  }
}
```

**Using Cromwell with file-based database (No server required)**

SQLite is currently not supported. However, HSQLDB does support running with a persistence file.
To set this up the following configuration can be used:
```hocon
database {
  profile = "slick.jdbc.HsqldbProfile$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = """
    jdbc:hsqldb:file:cromwell-executions/cromwell-db/cromwell-db;
    shutdown=false;
    hsqldb.default_table_type=cached;hsqldb.tx=mvcc;
    hsqldb.result_max_memory_rows=10000;
    hsqldb.large_data=true;
    hsqldb.applog=1;
    hsqldb.lob_compressed=true;
    hsqldb.script_format=3
    """
    connectionTimeout = 120000
    numThreads = 1
   }
}
```

Explanation of the options (see also http://hsqldb.org/doc/2.0/guide/dbproperties-chapt.html):

* `jdbc:hsqldb:file:cromwell-executions/cromwell-db/cromwell-db;` This will make sure
   all persistence files will end up in a folder `cromwell-db` inside `cromwell-executions`.
* `shutdown=false`. This makes sure the database will not be shutdown unless Cromwell explicitly does so.
* `hsqlldb.default_table_type=cached`. 
   By default hsqldb uses in memory tables, this will ensure data is written to disk and 
   decrease memory usage.
* `hsqldb.result_max_memory_rows=10000` . Limits the amount of rows in memory for temp tables. 
* `hsqldb.tx=mvcc` this is a  cromwell default for running with hsqldb.
* `hsqldb.large_data=true`. Cromwell creates huge DBs that need to be opened.
* `hsqldb.applog=1`. Log errors relating to the database.
* `hsqldb.lob_compressed=true`. Compress lobs. This saves some space. Do note that lobs are 
  compressed individually. The total database will still contain a lot of redundancy because a
  lot of lobs will be similar.
* `hsqldb.script_format=3`. Compress script. (uses gzip internally). 
   The script can still be opened normally after decompressing with gzip.
* `connectionTimeout = 120000` opening the large database files again when running cromwell will 
  take some time. The default timeout of 3000 ms (3s) is not enough. So it is set to 120000ms (120s).
* `numThreads = 1`. This will limit the CPU usage of Cromwell, which can be useful in HPC environments.

Comparison to MySQL (or PostgreSQL) server:
Advantages:

* No need to set up a server
* No worries about database users, passwords and permissions. This will be handled by filesystem permissions.

Disadvantages:

* Cromwell requires more memory
* The database files will consume a lot of disk space (multiple gigabytes are not uncommon)
* Cromwell's interaction with the database is slower.

Comparison to the default in-memory database:
Advantages:

* Much less memory needed.
* Call-caching enabled

Disadvantages:

* Slower.

### Abort

**Control-C (SIGINT) abort handler**

For backends that support aborting jobs, Cromwell can be configured to automatically try to abort all calls when it receives a Control-C, also known as SIGINT. All currently running calls will also set their status to `Aborted`.

To explicitly turn this feature on or off, set the configuration option:

```hocon
system {
  abort-jobs-on-terminate=true
}
```

Or, via `-Dsystem.abort-jobs-on-terminate=true` command line option.

By default, this value is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <workflow source> <inputs>`.

Read the [Abort](execution/ExecutionTwists/#abort) section to learn more about how abort works.

### Call caching

Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  
To learn more see [Call Caching](cromwell_features/CallCaching).

To enable Call Caching, add the following to your Cromwell configuration:

```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

When `call-caching.enabled=true` (default: `false`), Cromwell will be able to to reference or copy results from previously run jobs (when appropriate).
When `invalidate-bad-cache-results=true` (default: `true`), Cromwell will invalidate any cache results which contain files that cannot be accessed within a cache-hit. This is usually desired, but might be unwanted if this failure occurs for external reasons, such as a difference in user authentication.

Cromwell also accepts [Workflow Options](wf_options/Overview#call-caching-options) to override the cache read/write behavior.  

### Local filesystem options

When running a job on the Config (Shared Filesystem) backend, Cromwell provides some additional options in the backend's 
config section:

```HOCON
      config {
        filesystems {
          local {
            # When localizing a file, what type of file duplication should occur. 
            # possible values: "hard-link", "soft-link", "copy", "cached-copy".
            # For more information check: https://cromwell.readthedocs.io/en/stable/backends/HPC/#shared-filesystem
            localization: [
              "hard-link", "soft-link", "copy"
            ]

            caching {
              # When copying a cached result, what type of file duplication should occur. 
              # possible values: "hard-link", "soft-link", "copy", "cached-copy".
              # For more information check: https://cromwell.readthedocs.io/en/stable/backends/HPC/#shared-filesystem
              # Attempted in the order listed below:
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]

              # Possible values: md5, xxh64, fingerprint, path, path+modtime
              # For extended explanation check: https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching
              # "md5" will compute an md5 hash of the file content.
              # "xxh64" will compute an xxh64 hash of the file content. Much faster than md5
              # "fingerprint" will take last modified time, size and hash the first 10 mb with xxh64 to create a file fingerprint.
              # This strategy will only be effective if the duplication-strategy (above) is set to "hard-link", as copying changes the last modified time.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              # "path+modtime" will compute an md5 hash of the file path and the last modified time. The same conditions as for "path" apply here.
              # Default: "md5"
              hashing-strategy: "md5"
              
              # When the 'fingerprint' strategy is used set how much of the beginning of the file is read as fingerprint. 
              # If the file is smaller than this size the entire file will be read.
              # Default: 10485760 (10MB). 
              fingerprint-size: 10485760

              # When true, will check if a sibling file with the same name and the .md5 extension exists, and if it does, use the content of this file as a hash.
              # If false or the md5 does not exist, will proceed with the above-defined hashing strategy.
              # Default: false
              check-sibling-md5: false
            }
          }
        }
      }
```

#### Call cache strategy options for local filesystem

* hash based options. These read the entire file. These strategies work with containers.
    * `xxh64` (community-supported*). This uses the 64-bit implementation of the [xxHash](https://www.xxhash.com)
             algorithm. This algorithm is optimized for file integrity hashing and provides a more than 10x speed improvement over
             md5.
    * `md5`. The well-known md5sum algorithm
* Path based options. These are based on filepath. Extremely lightweight, but only work with the `soft-link` file 
caching strategy and can therefore never work with containers.
    * `path` creates a md5 hash of the path.
    * `path+modtime` creates a md5 hash of the path and its modification time.
* Fingerprinting. This strategy works with containers.
    * `fingerprint` (community-supported*) tries to create a fingerprint for each file by taking its last modified time (milliseconds since
       epoch in hexadecimal) + size (bytes in hexadecimal) + the xxh64 sum of the first 10 MB** of the file. 
       It is much more lightweight than the hash based options while still unique enough that collisions are unlikely. This 
       strategy works well for workflows that generate multi-gigabyte files and where hashing these files on the 
       cromwell instance provides CPU or I/O problems. 
       NOTE: This strategy requires hard-linking as a dupliation strategy, as copying changes the last modified time.

(*) The `fingerprint` and `xxh64` strategies are features that are community supported by Cromwell's HPC community. There
is no official support from the core Cromwell team.

(**) This value is configurable.
 
### Workflow log directory

To change the directory where Cromwell writes workflow logs, change the directory location via the setting:

```hocon
workflow-options {
    workflow-log-dir = "cromwell-workflow-logs"
}
```


**Preserving Workflow Logs**

By default Cromwell erases the per workflow logs when the workflow completes to reduce disk usage. You can change this behavior by setting the following value to `false`:

```hocon
workflow-options {
    workflow-log-temporary = true
}
```


**Exception monitoring via Sentry**

Cromwell supports [Sentry](https://docs.sentry.io) which is a service that can be used to monitor exceptions reported in an application’s logs.

To enable Sentry monitoring in Cromwell, enter your DSN URL using the system property:

```properties
sentry.dsn=DSN_URL
```


### Job shell configuration

Cromwell allows for system-wide or per-backend job shell configuration for running user commands rather than always
using the default `/bin/bash`. To set the job shell on a system-wide basis use the configuration key `system.job-shell` or on a
per-backend basis with `<config-key-for-backend>.job-shell`. For example:

```
# system-wide setting, all backends get this
-Dsystem.job-shell=/bin/sh
```

```
# override for just the Local backend
-Dbackend.providers.Local.config.job-shell=/bin/sh
```

For the Config backend the value of the job shell will be available in the `${job_shell}` variable. See Cromwell's `reference.conf` for an example
of how this is used for the default configuration of the `Local` backend.

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends

### Workflow Heartbeats

**Cromwell ID**

Each Cromwell instance is assigned a `cromwell_id`. By default, the Cromwell ID is `cromid-<7_digit_random_hex>`.
A custom identifier may replace the "cromid" portion of the string. For example:

```hocon
system {
  cromwell_id = "main"
}
```

This would generates a `cromwell_id` of `main-<7_digit_random_hex>`. Each time Cromwell restarts the random part of the
ID will change, however the `main` prefix would remain the same.

If the random part of the Cromwell ID should not be generated, set the configuration value:

```hocon
system {
  cromwell_id_random_suffix = false
}
```

**Heartbeat TTL**

When a Cromwell instance begins running or resuming a workflow it stores the above `cromwell_id` within the database row
for the workflow, along with a timestamp called the "heartbeat". As the workflow continues to run the Cromwell instance
will intermittently update the heartbeat for the running workflow.  If the Cromwell dies, after some time-to-live (TTL),
the workflow has been abandoned, and will be resumed by another available Cromwell instance.

Adjust the heartbeat TTL via the configuration value:

```hocon
system.workflow-heartbeats {
  ttl = 10 minutes
}
```

The default TTL is 10 minutes. The shortest allowable value for the TTL option is 10 seconds.

**Heartbeat Interval**

The interval for writing heartbeats may be adjusted via:

```hocon
system.workflow-heartbeats {
  heartbeat-interval = 2 minutes
}
```

The default interval is 2 minutes. The shortest interval option is 3.333 seconds. The interval may not be greater than
the TTL.

**Heartbeat Failure Shutdown**

Cromwell will automatically shutdown when unable to write heartbeats for a period of time. This period of time may be
adjusted via:

```hocon
system.workflow-heartbeats {
  write-failure-shutdown-duration = 5 minutes
}
```

The default shutdown duration is 5 minutes. The maximum allowed shutdown duration is the TTL.

**Heartbeat Batch Size**

Workflow heartbeats are internally queued by Cromwell and written in batches. When the configurable batch size is
reached, all of the heartbeats within the batch will be written at the same time, even if the heartbeat interval has not
elapsed.

This batch threshold may be adjusted via:

```hocon
system.workflow-heartbeats {
  write-batch-size = 100
}
```

The default batch size is 100.

**Heartbeat Threshold**

Cromwell writes one batch of workflow heartbeats at a time. While the internal queue of heartbeats-to-write passes above
a configurable threshold then [instrumentation](developers/Instrumentation.md) may send a metric signal that the
heartbeat load is above normal.

This threshold may be configured via the configuration value:

```hocon
system.workflow-heartbeats {
  write-threshold = 100
}
```

The default threshold value is 100, just like the default for the heartbeat batch size.

### YAML

**Maximum number of nodes**

Cromwell will throw an error when detecting cyclic loops in Yaml inputs. However one can craft small acyclic YAML
documents that consume significant amounts of memory or cpu. To limit the amount of processing during parsing, there is
a limit on the number of nodes parsed per YAML document.

This limit may be configured via the configuration value:

```hocon
yaml {
  max-nodes = 1000000
}
```

The default limit is 1,000,000 nodes.

**Maximum nesting depth**

There is a limit on the maximum depth of nested YAML. If you decide to increase this value, you will likely need to also
increase the Java Virtual Machine's thread stack size as well using
[either `-Xss` or `-XX:ThreadStackSize`](https://docs.oracle.com/javase/8/docs/technotes/tools/unix/java.html).

This limit may be configured via the configuration value:

```hocon
yaml {
  max-depth = 1000
}
```

The default limit is a maximum nesting depth of 1,000.
**Welcome to Cromwell**

Cromwell is a Workflow Management System geared towards scientific workflows. Cromwell is open sourced under the [BSD 3-Clause license](https://github.com/broadinstitute/cromwell/blob/develop/LICENSE.txt).

![Jamie, the Cromwell pig](jamie_the_cromwell_pig.png)
# Language Support

Below are the Domain Specific Languages (DSL) that Cromwell currently supports and will soon support for describing your workflow.

## Current Language Support

### WDL Draft 2
Cromwell started life as a WDL engine and WDL draft2 was our first language!
For many examples on how to use WDL and some great getting-started resources you can view [the OpenWDL site](https://github.com/openwdl/wdl#getting-started-with-wdl).

Cromwell supports the majority of [Draft-2 of the WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md).

> *Known Issues:*
>
> - Be careful when using `Object`. They are superceded by 'struct' in WDL 1.0 and are being removed outright in WDL 2.0.
> - Cromwell does not support nested `scatter`s in draft-2.


### WDL 1.0

Cromwell also supports WDL version 1.0.

As well as the changes to the WDL spec between draft-2 and 1.0, Cromwell also supports nested scatters and the [localization_optional](optimizations/FileLocalization.md) optimization in WDL 1.0.  


### CWL 1.0

Cromwell provides support for Common Workflow Language (CWL), beginning with the core spec, and most heavily used requirements.
If you spot a CWL feature that Cromwell doesn't support, please notify us using an issue on our github page!


## Future Language Support

### WDL 'development'

As the SPEC is being improved and honed, Cromwell continues to support the current `development` version of WDL. That
means that when (or shortly after) new versions are published, Cromwell will be ready to support them.
# Ecosystem

Many community projects have been created to support the use of the core Cromwell engine. Check them out to make your Cromwell experience even better!

* [Cromshell](https://github.com/broadinstitute/cromshell): Shell script for interacting with cromwell servers created at the [Broad Institute](https://github.com/broadinstitute).
* [Oliver](https://github.com/stjudecloud/oliver): An opinionated Cromwell orchestration manager created by the [St. Jude Cloud team](https://github.com/stjudecloud).# Cromwell Releases

Cromwell releases are available at the [GitHub Releases](https://github.com/broadinstitute/cromwell/releases/latest) page. 
You are strongly encouraged to use the latest release of Cromwell whenever possible.

Cromwell is distributed as a conda package on [conda-forge](https://conda-forge.org/).
These instructions need to be followed for [installing the miniconda distribution](https://docs.conda.io/en/latest/miniconda.html) and 
[activating the conda-forge channel](https://conda-forge.org/#about). After this Cromwell can be installed in the 
base environment with `conda install cromwell` or a separate environment for Cromwell can be created with 
`conda create -n cromwell cromwell`. If you are using Cromwell for bioinformatics workflows, you might like to take
a look at [bioconda](http://bioconda.github.io)  as well. 
The conda installation of Cromwell comes with a wrapper that locates the jar for you and allows for running Cromwell or Womtool with a 
`cromwell run`, `womtool validate` or other command. Conda also installs the required Java dependency 
in the environment automatically.

Mac users with Homebrew can also get Cromwell with the command `brew install cromwell`.

This documentation frequently refers to a "Cromwell jar" with a name like `cromwell-<version>.jar`. 
This is the main artifact in Cromwell releases that contains all executable Cromwell code and default configuration.   

A distribution of Java 11 is required to run Cromwell. Cromwell is developed, tested, and containerized using
[AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/).

For users running a Cromwell server [a docker image](https://hub.docker.com/r/broadinstitute/cromwell) has been made available.

### Apple Silicon support statement (updated 2020-11-17)

#### Cromwell JAR works out of the box

The Cromwell JAR works on any standard Java installation. A user can install an x86 Java runtime on an Apple Silicon Mac and the Rosetta 2 translation layer runs Cromwell at near-native speed.

Once natively-compiled Java runtimes become available, performance will increase with no change in functionality. 

#### Docker Desktop support is in progress

The Cromwell Docker image will not run on M1 Macs until Docker Desktop ships the appropriate update. For more details, please see [their official announcement](https://www.docker.com/blog/apple-silicon-m1-chips-and-docker/).

By extension, the absence of Docker means that Cromwell's local Docker backend is not yet supported.

Even when Docker Desktop goes native on Apple Silicon, any tool images running on the local backend will need to cross-compile for the x86 and Arm architectures. This is because the Rosetta 2 translation layer [does not support virtualization](https://developer.apple.com/documentation/apple_silicon/about_the_rosetta_translation_environment). Please contact the tool maintainers for more information. 

For built-in documentation of Cromwell command line usage, run the Cromwell JAR file with no arguments:

`$ java -jar cromwell-<versionNumber>.jar`

You will get a usage message like the following:

```bash
cromwell 29
Usage: java -jar /path/to/cromwell.jar [server|run] [options] <args>...

  --help                   Cromwell - Workflow Execution Engine
  --version                
Command: server
Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.
Command: run [options] workflow-source
Run the workflow and print out the outputs in JSON format.
  workflow-source          Workflow source file or workflow url .
  --workflow-root <value>  Workflow root
  -i, --inputs <value>     Workflow inputs file.
  -o, --options <value>    Workflow options file.
  -t, --type <value>       Workflow type.
  -v, --type-version <value>
                           Workflow type version.
  -l, --labels <value>     Workflow labels file.
  -p, --imports <value>    A zip file to search for workflow imports.
  -m, --metadata-output <value>
                           An optional JSON file path to output metadata.
```

Cromwell's Server and Run modes can be invoked with the `server` and `run` arguments respectively. More information on these Cromwell modes can be found in [Modes](Modes).

The Cromwell jar file can be built as described in [Building](Building). 

## `server`

`server` mode accepts no arguments and runs Cromwell as a web server that accepts REST requests. The default mode for most applications of Cromwell, suitable for production use. See the documentation for [Cromwell's REST endpoints](/api/RESTAPI) for how to interact with Cromwell in `server` mode.

## `run`

`run` mode executes a single workflow in Cromwell and then exits. It is designed for local prototyping or demos and has limited features compared to `server`.

* **`workflow-source`**  
The single required argument. It can be either a local path or a remote URL pointing to the workflow source file.
 
* **`--inputs`**  
An optional file of workflow inputs.  Although optional, it is a best practice to use an inputs file to satisfy workflow
requirements rather than hardcoding inputs directly into a workflow source file.

* **`--options`**  
An optional file of workflow options.  Some options are global (supported by all backends), while others are backend-specific. See the [Workflow Options](wf_options/Overview) for more details.

* **`--type`**  
An optional parameter to specify the language for the workflow source. As of Cromwell 29 any value specified for this parameter is currently ignored and internally the value `WDL` is used.

* **`--type-version`**  
An optional parameter to specify the version of the language for the workflow source. Currently any specified value is ignored.

* **`--labels`**  
An optional parameter to specify a file of JSON key-value label pairs to associate with the workflow.

* **`--imports`**  
You have the option of importing WDL workflows or tasks to use within your workflow, known as sub-workflows.
If you use sub-workflows within your primary workflow then you must include a ZIP file with the WDL import files.
See the documentation on [Imports](Imports) for more information.

* **`--metadata-output`**  
You can specify a filename where Cromwell will write workflow metadata JSON such as start/end timestamps, status, inputs and outputs. By default Cromwell does not write workflow metadata. The metadata format in the `--metadata-output` file is the same as described for the [REST API](api/RESTAPI#get-workflow-and-call-level-metadata-for-a-specified-workflow).

* **`--version`**  
The `--version` option prints the version of Cromwell and exits.

* **`--help`**  
The `--help` option prints the full help text above and exits.
## Logging Properties

### Setting Logging Properties

Cromwell accepts two properties for controlling logging. You can set these properties via a Java system property on the command line using `-D`:

```bash
$ java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

Alternatively, you can also set the log level via an environment variable:

```bash
export LOG_LEVEL=DEBUG
java -jar cromwell.jar server
```

*If you set same property via a system property, and an environment variable, the system property overrides the environment variable.*

### Log Format

Cromwell outputs log in one of two formats, either `pretty` or `standard`. You can change the format of the logs by setting the property to `LOG_MODE`.

* In `standard` mode, your logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs.

* In `pretty` mode, your logs are output in a colorful, easier to read format, more appropriate for a single workflow run.

The default mode for server is `standard`, while the default when running a single worklow is `pretty`. You can explicitly specify the format by running cromwell with:

```bash
java -DLOG_MODE=pretty -jar cromwell.jar server
```

### Log Level

By default, Cromwell outputs messages at a `LOG_LEVEL` of `INFO`. Sometimes, you may want more or less information logged. For example, while debugging an issue you may want to increase the amount information in the logs temporarily. Alternatively, the standard level may be too verbose, and you may only want Cromwell to log warnings and errors.

You can set the level via the property `LOG_LEVEL` to any one of the values: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, or `OFF`. The default log level is `INFO`.

```bash
java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

## Workflow Logs

While a workflow is running, Cromwell generates a log file specifically for the workflow. After the workflow completes, to clear up local disk space, Cromwell deletes the local copy of this log file. See the [Configuration](Configuring#workflow-log-directory) section on logs for more information on preventing cromwell from deleting each workflow log.

Before Cromwell deletes the files and before the workflow completes, you can configure Cromwell to copy the workflow
logs to various locations. Normally, you'll want to copy the log to a remote bucket or directory. To specify the remote
directory to copy the logs to, use the separate [Workflow Option](wf_options/Overview#output-copying)
`final_workflow_log_dir`.

## Call Logs

As each call in a workflow runs, it generates output to the standard output and standard error. This output is stored per call in call log files. Additionally, depending on the backend, specific per call backand logs may be generated.

All of these call logs may be copied at the end of a workflow to a remote directory. Configure this directory by setting the [Workflow Option](wf_options/Overview#output-copying) `final_call_logs_dir`.
In order to support the composition and reuse of workflows, WDL allows the execution of an entire workflow as a step in a larger workflow.  When a workflow calls another workflow, that second workflow is referred to as a sub-workflow.  Note that sub-workflows can themselves contain sub-workflows and so on, and there is no explicit limit as to how deeply workflows can be nested.  Cromwell supports execution of such workflows.

However, a single WDL file can contain only a single workflow definition.  In order to reference a sub-workflow, the import directive can be used to bring that sub-workflow into existence and referenced by it's alias and name.  See the [documentation on Imports](Imports.md) for more details of how to declare and reference tasks and workflows via imports.

**Execution**

A sub-workflows is executed exactly as a task would be.
*This means that if another call depends on an output of a sub-workflow, this call will run when the whole sub-workflow completes (successfully).*
For example, in the following case :

`main.wdl`
```
import "sub_wdl.wdl" as sub

workflow main_workflow {

    call sub.hello_and_goodbye { input: hello_and_goodbye_input = "sub world" }
    
    # call myTask { input: hello_and_goodbye.hello_output }
    
    output {
        String main_output = hello_and_goodbye.hello_output
    }
}
```

`sub_wdl.wdl`
```
task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
}

task goodbye {
  String addressee
  command {
    echo "Goodbye ${addressee}!"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow hello_and_goodbye {
  String hello_and_goodbye_input
  
  call hello {input: addressee = hello_and_goodbye_input }
  call goodbye {input: addressee = hello_and_goodbye_input }
  
  output {
    String hello_output = hello.salutation
    String goodbye_output = goodbye.salutation
  }
}
```

`myTask` will start only when hello_and_goodbye completes (which means all of its calls are done), even though `myTask` only needs the output of hello in the hello_and_goodbye sub-workflow. 
If hello_and_goodbye fails, then `myTask` won't be executed.
Only workflow outputs are visible outside a workflow, which means that references to outputs produced by a sub-workflow will only be valid if those outputs are exposed in the workflow output section.

Sub-workflows are executed in the context of a main workflow, which means that operations that are normally executed once per workflow (set up, clean up, outputs copying, log copying, etc...)
will NOT be re-executed for each sub-workflow. For instance if a resource is created during workflow initialization, sub-workflows will need to share this same resource.
Workflow outputs will be copied for the main root workflow but not for intermediate sub-workflows.

Restarts, aborts, and call-caching work exactly as they would with tasks. 
All tasks run by a sub-workflow are eligible for call caching under the same rules as any other task.
However, workflows themselves are not cached as such. Which means that running the exact same workflow twice with call caching on will trigger each task to cache individually,
but not the workflow itself.

The root path for sub-workflow execution files (scripts, output files, logs) will be under the parent workflow call directory.
For example, the execution directory for the above main workflow would look like the following:

```
cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/ <- main workflow id
└── call-hello_and_goodbye <- call directory for call hello_and_goodbye in the main workflow
    └── hello_and_goodbye <- name of the sub-workflow 
        └── a6365f91-c807-465a-9186-a5d3da98fe11 <- sub-workflow id
            ├── call-goodbye
            │   └── execution
            │       ├── rc
            │       ├── script
            │       ├── script.background
            │       ├── script.submit
            │       ├── stderr
            │       ├── stderr.background
            │       ├── stdout
            │       └── stdout.background
            └── call-hello
                └── execution
                    ├── rc
                    ├── script
                    ├── script.background
                    ├── script.submit
                    ├── stderr
                    ├── stderr.background
                    ├── stdout
                    └── stdout.background

```

**Metadata**

Each sub-workflow will have its own workflow ID. This ID will appear in the metadata of the parent workflow, in the call section corresponding to the sub-workflow, under the "subWorkflowId" attribute.
For example, querying the `main_workflow` metadata above (minus the `myTask` call) , could result in something like this:

`GET /api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [
      {
        "executionStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "end": "2016-11-17T14:13:41.117-05:00",
        "attempt": 1,
        "start": "2016-11-17T14:13:39.236-05:00",
        "subWorkflowId": "a6365f91-c807-465a-9186-a5d3da98fe11"
      }
    ]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {},
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```

The sub-workflow ID can be queried separately:

`GET /api/workflows/v2/a6365f91-c807-465a-9186-a5d3da98fe11/metadata`

```
{
  "workflowName": "hello_and_goodbye",
  "calls": {
    "sub.hello_and_goodbye.hello": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Hello sub world!"
        },
        "runtimeAttributes": {
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49830",
        "backend": "Local",
        "end": "2016-11-17T14:13:40.712-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
        "callRoot": "/cromwell/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:40.704-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:40.704-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:40.712-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ],
    "sub.hello_and_goodbye.goodbye": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Goodbye sub world!"
        },
        "runtimeAttributes": {
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49831",
        "backend": "Local",
        "end": "2016-11-17T14:13:41.115-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
        "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:41.112-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:41.112-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:41.115-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ]
  },
  "outputs": {
    "goodbye_output": "Goodbye sub world!",
    "hello_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
  "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
  "inputs": {
    "hello_and_goodbye_input": "sub world"
  },
  "status": "Succeeded",
  "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
  "end": "2016-11-17T14:13:41.116-05:00",
  "start": "2016-11-17T14:13:39.236-05:00"
}
```

It's also possible to set the URL query parameter `expandSubWorkflows` to `true` to automatically include sub-workflows metadata (`false` by default).

`GET api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata?expandSubWorkflows=true`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [{
      "executionStatus": "Done",
      "subWorkflowMetadata": {
        "workflowName": "hello_and_goodbye",
        "calls": {
          "sub.hello_and_goodbye.hello": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Hello sub world!"
            },
            "runtimeAttributes": {
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49830",
            "backend": "Local",
            "end": "2016-11-17T14:13:40.712-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
            "callRoot": "cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:40.704-05:00"
            }, {
              "startTime": "2016-11-17T14:13:40.704-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:40.712-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }],
          "sub.hello_and_goodbye.goodbye": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Goodbye sub world!"
            },
            "runtimeAttributes": {
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49831",
            "backend": "Local",
            "end": "2016-11-17T14:13:41.115-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
            "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:41.112-05:00"
            }, {
              "startTime": "2016-11-17T14:13:41.112-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:41.115-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }]
        },
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
        "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "status": "Succeeded",
        "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
        "end": "2016-11-17T14:13:41.116-05:00",
        "start": "2016-11-17T14:13:39.236-05:00"
      },
      "shardIndex": -1,
      "outputs": {
        "goodbye_output": "Goodbye sub world!",
        "hello_output": "Hello sub world!"
      },
      "inputs": {
        "hello_and_goodbye_input": "sub world"
      },
      "end": "2016-11-17T14:13:41.117-05:00",
      "attempt": 1,
      "start": "2016-11-17T14:13:39.236-05:00"
    }]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {

  },
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```**Support Forum**

If you have questions that aren't covered by this documentation you can ask them in the
[Support Forum](http://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team).
Please note that despite the `wdl` currently in its name, this is a forum that supports Cromwell users too.

**Website and User Guide**

The [WDL website](https://software.broadinstitute.org/wdl/) is the best place to go for more information on WDL.
In particular new users should check out the [user guide](https://software.broadinstitute.org/wdl/userguide/)
which has many tutorials, examples and other bits to get you started.
Sometimes you might want to break up a long WDL file into smaller components for easier maintenance. Sometimes you'll want to do it to reuse components in multiple workflows.  Have no fear, imports are here!  Imports allow you to reference other WDL files that contain entire workflows or even just raw tasks.

To import a WDL, you can use the `import` WDL construct at the top of your workflow:

```
import "<resource>" as <alias>
```

There are two types of resources that are supported in imports: *http(s)* and *file-path based*.  Any public http(s) based URL can be used as the resource for an import, such as a website, github or a GA4GH compliant TES endpoint.  For example:

```wdl
import "http://mywdlrepository/my.wdl" as http_import1
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import2
```
To use a file-based import resource, provide a ZIP bundle of your resources and then use a path relative to that ZIP in your import statement. For example:

```wdl
import "my-wdl-in-the-root-directory.wdl" as file_import1
import "my/wdl/sub/directory/example.wdl" as file_import2
```

Imports from your submitted workflow are evaluated relative to the base of the zip file. In other cases, import paths are relative to the file you are currently importing from, for example:

>If there exists a `my/wdl/sub/directory/imports/importing_an_import.wdl`, `my/wdl/sub/directory/example.wdl` could import it relatively like so:
>```wdl
>import "imports/importing_an_import.wdl" as file_import3
>```

Here's a complete example showing both http(s) and file-based imports workflow in WDL:

_workflow.wdl_
```wdl
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import
import "imports/imported.wdl" as provided_import

workflow my_workflow {
    ...
}
```

Note that we said "`import "imports/imported.wdl`" in the workflow so we must have a `imports/imported.wdl` structure in the imports file:
_imports.zip_
```
imports
└── imported.wdl
```

A more common scenario might have the imports at the root of the imports.zip:

_workflow.wdl_
```wdl
import "my_wdl_1.wdl"
import "my_wdl_2.wdl"
```
_imports.zip_
```
my_wdl_1.wdl
my_wdl_2.wdl
```

---
The mechanism to provide the ZIP file of resources to be imported differs between [Run](Modes#run) and [Server](Modes#server) mode.

In [Run](Modes#run) mode, a sample command to run _workflow.wdl_ would be:  
```
$ java -jar cromwell.jar run workflow.wdl --imports imports.zip
```

In [Server](Modes#server) mode, pass in a ZIP file using the parameter `workflowDependencies` via the [Submit](api/RESTAPI#submit-a-workflow-for-execution) endpoint.


Cromwell can be thought of as performing three fundamental roles:

1. Front end - handling all REST requests
1. Runner - picking up and running workflows
1. Summarizer - summarizing metadata

In the simplest Cromwell deployment, a single Cromwell instance performs all three of these roles.
But it is also possible to create Cromwell deployments with many Cromwell instances, gaining advantages in
resiliency and scalability. There is only one restriction on roles: there must be exactly one Cromwell
instance performing the role of summarizer. Apart from this restriction, any instance in a multi-Cromwell
deployment can perform one or more of these roles. As a real world example, the current production Terra
deployment uses 7 Cromwell instances, each dedicated to a single role: 3 front ends, 3 runners, and one summarizer.

These roles are not explicit in Cromwell configuration, rather they are implied by configuration settings as
described below.

** Summarizer configuration **

The frequency of metadata summarization is determined by the value for the configuration key
`services.MetadataService.config.metadata-summary-refresh-interval`, which has a default value of `1 second`.
As stated above, there must be exactly one Cromwell instance performing the role of summarizer, so
all Cromwell instances which are not performing the summarizer role should specify `Inf` for this value.

** Runner configuration **

Cromwell instances in the runner role should periodically scan the workflow store to pick up and run
unclaimed workflows. The relevant configuration parameters are described below along with their default values.
The default values may be adequate for instances in the runner role, but will need to be overridden for
non-runner instances to effectively turn running off.

```hocon
system {
...
  # Number of seconds between polls of the workflow store.
  # Set this to a very large value for non-runners (e.g. 999999)
  new-workflow-poll-rate = 20

  # Cromwell will launch up to N submitted workflows at a time, regardless of how many open workflow slots exist
  # Set this to 0 for a non-runner.
  max-workflow-launch-count = 50

  # The maximum number of workflows to run concurrently.
  # Set this to 0 for a non-runner.
  max-concurrent-workflows = 5000
...
}
```

The documentation on [workflow heartbeats](https://cromwell.readthedocs.io/en/stable/Configuring/#workflow-heartbeats) describes how multiple Cromwell
runners collaborate to run workflows from a single workflow store.

** Front end configuration **

Cromwell instances should not require any configuration changes to operate in the front end role. If a particular
Cromwell instance is not intended to operate in the front end role then requests should not be directed to
that instance. If there are multiple front end instances then it may be desirable to configure a load balancer
in front of these instances to direct requests to the front end instances only.

## Server

The default mode for most applications of Cromwell, suitable for production use. Server mode starts Cromwell as a web server that exposes REST endpoints. All features and APIs are available.

By default the server will be accessible at `http://localhost:8000`. See the [Server Section](Configuring#server) of the configuration for more information on how to configure it. A description of the endpoints can be found in the [API Section](api/RESTAPI).

Follow the [Server Tutorial](tutorials/ServerMode) to get your Cromwell server up and running in a few steps.

## Run

A good way to get started with Cromwell and experiment quickly. Run mode launches a single workflow from the command line and exits `0` or `1` to indicate the result. Appropriate for prototyping or demo use on a user's local machine. Features are limited and the web API is not supported.

Sending a `SIGINT` signal (via `CTRL-C` for example) will by default abort all running jobs and then exit.
This behavior can be configured, and is explained in more details in the [Abort](Configuring#abort) section of the configuration.
# Retry with More Memory

With this feature one can specify an array of strings which when encountered in the `stderr` file by Cromwell, 
allows the task to be retried with more memory. The retry will be counted against the `maxRetries` count mentioned in 
the `runtimeAtrributes` in the task. There are 2 settings for this feature:

* `system.memory-retry-error-keys` : the error keys that need to be set in Cromwell config
* `memory_retry_multiplier` : [optional] the factor by which the memory should be multiplied while retrying. This needs 
to be passed in through workflow options and should be in the range `1.0 ≤ multiplier ≤ 99.0` (note: if set to `1.0` the task
will retry with same amount of memory). If this is not specified, Cromwell will retry the task (if applicable) but not 
change the memory amount.

For example, if the error keys set in Cromwell config are as below, and the multiplier passed through workflow options is 
`"memory_retry_multiplier": 1.1` 
```hocon
system {
  memory-retry-error-keys = ["OutOfMemory", "Killed"]
}
```  
this tells Cromwell to retry the task with 1.1x memory when it sees either `OutOfMemoryError` or `Killed` in the `stderr` 
file. 

If the task has runtime attributes as below 
```hocon
runtimeAtrributes {
  memory: "1 GB"
  maxRetries: 1
}
``` 
the task will be retried 1 more time if it runs out of memory, and this time with "1.1 GB".

If the task return code is 0, the task will not be retried with more memory, even if the `stderr` file contains a
string present in `system.memory-retry-error-keys`. Similarly, if the runtime attribute `continueOnReturnCode` is
specified as a true, or the return code of the task matches a value specified by `continueOnReturnCode`, the task
will be considered successful and will not be retried with more memory.

Please note that this feature currently only works in Google Cloud backend. Also, Pipelines API might adjust the 
memory value based on their standards for memory for a VM. So it's possible that even though the request says 1.1 GB 
memory, it actually allocated a bit more memory to the VM.

Two environment variables called `${MEM_UNIT}` and `${MEM_SIZE}` are also available inside the command block of a task,
making it easy to retrieve the new value of memory on the machine.
# Hog Factors

## Introduction

Cromwell has only a finite amount of resources at its disposal. WDL and CWL workflows allow scattered tasks to be
run a huge number of times with very simple syntax.
This makes it easy for a very small number of workflows to hog all of the resources of Cromwell, forcing all other 
workflows (even a simple 'hello_world') to wait in line behind them. 
Sometimes that's exactly what you want, but often you would like Cromwell to remain responsive to new users' small 
workflows even while continuing to process production workflows from established users.

Cromwell 35 provides new ways of stopping one workflow or group of workflows from locking out everyone else by 
introducing hog factors, hog groups and hog limits. This page describes what they are and how they work.

## Concepts

### Hog Group 

The Hog Group is a way of grouping workflow from different submissions together to restrain their overall resource 
usage as a whole.

- Every top-level workflow is assigned to a hog group when Cromwell receives it. 
    + Exactly how this happens is [configurable](#configuration).
- Every sub-workflow or call started by a workflow is associated with the same hog group as its parent workflow.
- Multiple top-level workflows can be assigned to the same hog group to allow them to be grouped together.

Thus:

- Every workflow is assigned to a hog group when submitted.
- Every hog group may have many workflows assigned to it.

### Hog Factor

The hog factor is an integer greater than or equal to 1. It represents a trade-off between: 

- Fully utilizing all resources available to Cromwell to complete jobs, for as long as there are jobs to be processed.
- Reserving resources for requests from other hog groups - even if we have jobs waiting to run that could be using them.

Here are a few mental models which might be helpful to thinking about the hog factor:

- A hog factor of 2 means that "2 greedy users would be able to hog the entire resources of Cromwell" 
- A hog factor of 100 means "any 1 group is only ever allowed to use 1/100th of the resources of the total Cromwell server"

### Hog Limit

A Hog Limit is how much of a given resource a hog group is allowed use. Hog limits are not set directly; they are 
values that Cromwell calculates internally.
For example, a single hog group may be limited by Cromwell to a hog limit of 200 jobs per group. Therefore no matter how 
many workflows, sub-workflows, and jobs are queued in a greedy hog group, the whole group is limited to 200 concurrent 
running jobs.


## Configuration

Cromwell accepts the following configuration values for hog factors in the `hog-safety` stanza of `system` in the configuration 
file:
```conf
system {
  hog-safety {
    hog-factor = 1
    workflow-option = "hogGroup"
    token-log-interval-seconds = 0
  }
}
```

Additionally, you can override system-level hog-factor on a backend level, by setting it in the particular backend configuration like this:
```conf
backend {
  providers {
    PAPIv2 {
      config {
        hog-factor: 2
      }
    }
  }
}
```

### Setting a hog-factor

The hog factor option sets the integer described in the [Hog Factor](#hog-factor) section above.
The default value is `1` (which is equivalent to not limiting by hog group). 

### Assignment of hog groups

Within the configuration file, you can specify the workflow option that will determine the hog group. The default is
`hogGroup`. So if a workflow arrives with the following workflow options file, Cromwell will assign `hogGroupA` as the 
workflow's hog group:

```json
{
  "hogGroup": "hogGroupA"
}
```

- Any workflow option value can be used so long as it is a simple `String` value:
    + You can come up with a new field and set it specifically for assigning hog groups.
    + You can choose a field that is already being used for other reasons
- If a workflow is submitted without a value for the designated field in its workflow options, the workflow ID is used 
as the hog group identifier. 

### Logging

Because the system is not a simple first-in-first-out, it can be valuable to see the status of all 
the existing queues inside Cromwell. 

To have this information logged on a regular basis, you can enable periodic queue logging. This will
also alert you at the same frequency when events are happening, such as hog groups being at their individual limits.

* You can enable logging for the Job Execution Token Dispenser, using
the `system.hog-safety.token-log-interval-seconds` configuration value.
* The default, `0`, means that no logging will occur.

## Effects

### Job Execution

#### Reserving Space

- Cromwell allows administrators to designate an overall maximum concurrent job limit per
 [backend](../backends/Backends.md#backend-job-limits). 
- Within that limit, a hog factor allows us to limit the maximum concurrent jobs started *per hog group*.
    + This is what allows new jobs to run immediately even if many jobs from bigger workflows are already queued up.

#### Round robin allocation

Rather than starting jobs on a strict first-come first served basis, Cromwell now assigns in a round-robin
fashion between hog groups and then on a first-come-first-served *within* a hog group.

In other words if the hog groups had the following entries queued up:
```
 A: jobA1, jobA2, jobA3, ..., jobA1000000
 B: jobB1, jobB2
 C: jobC1
 D: jobD1, jobD2
```

Then Cromwell would start the jobs in the following order, even though `jobA1000000` was added before `jobD1`:
```
jobA1, jobB1, jobC1, jobD1, jobA2, jobB2, jobD2, jobA3, ..., jobA1000000
```

#### Example: How job execution is affected by hog factors

##### An administrator sets up a Cromwell server

- A Cromwell administrator sets the overall maximum concurrent job limit to 100,000 PAPIv2 jobs.
- The administrator also sets the hog factor to be 25.
- Cromwell will therefore calculate a per-hog-group concurrent job limit of 4,000 PAPIv2 jobs.

##### Our first hog group hits its limit

- 100 workflows are running in hog group "A" and between them have generated 20,000 jobs for PAPIv2. 
    + Cromwell initially starts 4,000 jobs.
    + Cromwell then starts the remaining 16,000 new jobs as existing jobs from this group finish.
    + New workflows in this group will not be able to start jobs either
        + Their jobs are queued behind the existing jobs from this hog group.
    + Note that Cromwell is currently only using 1/25th of its overall limit because its hog factor is 25.

##### Another hog group appears

- Now hog group B submits 1,000 workflows and between them they generate 200,000 jobs for PAPIv2.
- Even though 16,000 of group A's jobs are still queued, Cromwell starts 4,000 of group B's jobs immediately,
- The remaining 196,000 of group B's jobs are queued up waiting for group B's existing jobs to complete.

##### Where do we stand?

- Cromwell knows about 220,000 jobs that could be started
- Cromwell has an overall limit of 100,000
- Cromwell is running 8,000 jobs in two hog groups.

*In other words, not so great - perhaps we should have set the hog factor lower...?*

#####But wait, more workflows appear...

- Now another 23 hog groups ("C" through "Y") submit workflows of a similar scale to hog group A.
- One by one, the workflows of each hog group fill up their share of the overall concurrent job limit.
- So Cromwell is now running 100,000 jobs and each hog group has been allocated 4,000 of those.

##### What about poor hog group "Z"?

- A final group submits workflows under hog group "Z".
- Alas, even though hog group "Z" is not running anything yet, we cannot start their workflows because we're 
at the global maximum of 100,000.


*In other words, perhaps we should have set the hog factor higher...?*

##### So what now?

- As jobs in other hog group complete, we will begin to see hog group "Z" jobs started alongside new jobs from the 
other hog groups.
- Going forward Cromwell will start jobs from all groups at the same rate, even though hog group Z's jobs arrived
later than those from hog group A. Thus, over time, each group will approach approximately 1/26th of the total pool.
  
## FAQs

#### Can I opt out of using hog groups?

Yes, to various degrees:

- No matter what, your workflows will be assigned to a hog group. 
- To opt out of reserving Cromwell's resources for new hog groups, leave the hog factor set to 1.
- To opt out of round-robin allocation between workflows, and preserve a strict first-in-first-out allocation of jobs,
assign all workflows to the same hog-group in their workflow options.
    + To set this as the default, you can add a value to the default workflow options. 
    + For an example see the `workflow-options` / `default` stanza of [cromwell.examples.conf][cromwell-examples-conf].

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf
Labels in Cromwell are a way to group together workflows that are related or associated to each other.
For example, if you ran workflows to analyze data from a diabetes study, you can assign the label `project:diabetes-cohort`.  

**Custom Labels JSON**

In order to assign labels to a workflow, the first step is to create a JSON file with key-value pairs that define a label. For the example above, the labels JSON should look like:

```
{
  "project":"diabetes-cohort"
}
```

When choosing key-value pairs, it's important to make sure you're adhering to Cromwell supported label syntax below.  

There are two ways to add labels to a workflow:  
1. Upon workflow submission set the `labels` parameter of the [Submit endpoint](api/RESTAPI#submit-a-workflow-for-execution), or  
2. Setting the `-l` argument when running in [Command Line](/CommandLine) mode.

Labels can be added to existing workflows by using the [Labels patch endpoint](api/RESTAPI#update-labels-for-a-workflow).

After adding labels to your workflows, you can take advantage of features like [Query](api/RESTAPI#get-workflows-matching-some-criteria) to filter tagged workflows. The Google backend supports labelling cloud resources and you can learn more about that [here](backends/Google#google-labels).

#### Label Format

When labels are supplied to Cromwell, it will fail any request containing invalid label strings. Below are the requirements for a valid label key/value pair in Cromwell:

* Label keys may not be empty but label values may be empty.
* Label key and values have a max char limit of 255.

For [default labels](backends/Google#google-labels) applied by the Google backend, Cromwell will modify workflow/task/call names to fit the schema, according to the following rules:

* Any capital letters are converted to lowercase.
* Any character which is not one of `[a-z]`, `[0-9]` or `-` will be replaced with `-`.
* If the start character does not match `[a-z]` then prefix with `x--`
* If the final character does not match `[a-z0-9]` then suffix with `--x`
* If the string is too long, only take the first 30 and last 30 characters and add `---` between them.
Call Caching allows Cromwell to detect when a job has been run in the past so that it doesn't have to re-compute results, saving both time and money.  Cromwell searches the cache of previously run jobs for one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will use the results of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  In order for call caching to be used on any previously run jobs,
it is best to configure Cromwell to [point to a MySQL database](../Configuring.md#database) instead of the default
in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to
utilize results from all calls that are in that database.

**Configuring Call Caching**

*Call Caching is disabled by default.*  Call Caching can be enabled in your Cromwell
[Configuration](../Configuring.md#call-caching) and the behavior can be modified via
[Workflow Options](../wf_options/Overview.md). If you are adding Workflow options, do not set
[`read_from_cache` or `write_to_cache`](../wf_options/Overview.md#call-caching-options) = false, as it will impact the
following process.

Once enabled, Cromwell by default will search the call cache for every `call` statement invocation.

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache.
* If there was a cache hit, outputs are either **copied from the original cached job to the new job's output directory**
or **referenced from the original cached job** depending on the Cromwell
[Configuration](../Configuring.md#call-caching) settings.

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.  Additionally, if any files from a previous job directory are removed, call caching will fail due to missing files.

***File hash caching***

Cromwell offers the option to cache file hashes within the scope of a root workflow to prevent repeatedly requesting the hashes of the
same files multiple times. File hash caching is off by default and can be turned on with the configuration option `system.file-hash-cache=true`.

***Call cache blacklisting***
Cromwell offers the ability to filter cache hits based on copying failures. 

Call cache blacklisting configuration looks like:

```
call-caching {

  enabled = true

  # In a multi-user environment this should be false so unauthorized users don't invalidate results for authorized users. 
  invalidate-bad-cache-results = false

  blacklist-cache {
     # The call caching blacklist cache is off by default. This cache is used to blacklist cache hits based on cache
     # hit ids or buckets of cache hit paths that Cromwell has previously failed to copy for permissions reasons.
     enabled: true

     # All blacklisting values below are optional. In order to use groupings (blacklist caches shared among root
     # workflows) a value must be specified for `groupings.workflow-option` in configuration and the workflows to
     # be grouped must be submitted with workflow options specifying the same group.
     groupings {
       workflow-option: call-cache-blacklist-group
       concurrency: 10000
       ttl: 2 hours
       size: 1000
     }

     buckets {
       # Guava cache concurrency.
       concurrency: 10000
       # How long entries in the cache should live from the time of their last access.
       ttl: 1 hour
       # Maximum number of entries in the cache.
       size: 1000
     }

     hits {
       # Guava cache concurrency.
       concurrency: 10000
       # How long entries in the cache should live from the time of their last access.
       ttl: 1 hour
       # Maximum number of entries in the cache.
       size: 20000
     }
  }
}
```

**** Blacklist cache grouping ****

By default Cromwell's blacklist caches work at the granularity of root workflows, but Cromwell can also be configured to
share a blacklist cache among a group of workflows. 
If a value is specified for `call-caching.blacklisting.groupings.workflow-option` and a workflow option is specified
having a matching key, all workflows specifying the same value will share a blacklist cache. 

For example, if Cromwell configuration contains `call-caching.blacklisting.groupings.workflow-option = "project"` and
a workflow is submitted with the options

```json
{
  "project": "Mary"
}
```

then this workflow will share a blacklist cache with any other workflows whose workflow options contain `"project": "Mary"`.

Grouping of blacklist caches can significantly improve blacklisting effectiveness and overall call caching performance.
Workflows should be grouped by their effective authorization to ensure the same filesystem/object store permissions
exist for every workflow in the group.

**** Hit blacklisting ****

If a cache hit fails copying for any reason, Cromwell will record that failure in the blacklist cache and will not use
the hit again. Hit blacklisting is particularly effective at improving call caching performance in conjunction with the 
grouping feature described above.

**** Path prefix (GCS bucket) blacklisting on 403 Forbidden errors ****

In a multi-user environment user A might cache hit to one of user B's results
but that doesn't necessarily mean user A is authorized to read user B's outputs from the filesystem. Call cache blacklisting
allows Cromwell to record which file path prefixes were involved in cache result copy authorization failures.
If Cromwell sees that the file paths for a candidate cache hit have a blacklisted prefix, Cromwell will quickly 
fail the copy attempt without doing any potentially expensive I/O.

Path prefix blacklisting could be supported by any backend type though it is currently implemented only for Google
(PAPI) backends. For Google backends the GCS bucket is considered the prefix for blacklisting purposes.


***Call cache whitelisting***
 
In a multi-user environment where access to job outputs may be restricted among different users, it can be useful to limit
cache hits to those that are more likely to actually be readable for cache hit copies.
Cromwell now supports a `call_cache_hit_path_prefixes` workflow option for this purpose. This is particularly useful in the PAPI backend where the workflow
root can be specified in workflow options via `jes_gcs_root`. The value of `call_cache_hit_path_prefixes` should be an array of strings representing  
prefixes that call cache hit output files should have in order to be considered as a cache hit. Using PAPI as an example and assuming Alice and Bob have
made their data accessible to each other, Alice could submit a workflow with these options:

```
{
  "call_cache_hit_path_prefixes": [ "gs://alice_bucket", "gs://bob_bucket" ]
}
```

With these workflow options Cromwell would only look for cache hits for Alice's jobs in Alice's or Bob's buckets.

As a further optimization the PAPI backend has the concept of "this" bucket on a per-workflow basis, where "this" bucket is
the bucket that contains the current workflow root.
If `call_cache_hit_path_prefixes` is specified in 
workflow options on the PAPI backend, Cromwell will automatically prepend "this" bucket to the call cache hit path prefixes to search.
For example, if Charles specified the same workflow options as in the example above and his workflow root was under `gs://charles_bucket`,
Cromwell would search cache hits in all of the `gs://alice_bucket`, `gs://bob_bucket` and `gs://charles_bucket` buckets without having to specify
`gs://charles_bucket` bucket explicitly in `call_cache_hit_path_prefixes`.

If no `call_cache_hit_path_prefixes` are specified then all matching cache hits will be considered.

***Call cache failure logging***

When Cromwell fails to cache a job from a previous result the reason will be logged. To reduce the verbosity of the logs
only the first three failure reasons will be logged per shard of each job. Cromwell will continue to try copying
previous results for the call, and when no candidates are left Cromwell will run the job on the backend.

**Docker Tags**

Certain Docker tags can impact call caching performance. 
Docker tags are a convenient way to point to a version of an image (`ubuntu:14.04`), or even the latest version (`ubuntu:latest`).
For that purpose, tags are mutable, meaning that the image they point to can change, while the tag name stays the same.
While this is very convenient in some cases, using mutable, or "floating" tags in tasks affects the reproducibility of a workflow. 
If you were to run the same workflow using `ubuntu:latest` now, and again in a year (or even in a month) may run with different docker images.
This has an even bigger impact when Call Caching is turned on in Cromwell, and could lead to unpredictable behaviors if a tag is updated in the middle of a workflow or even a scatter for example.

In order to ensure perfect reproducibility, Docker provides another way of identifying an image version by using the specific digest of the image, which is an immutable identifier. The digest is guaranteed to be different if 2 images have different byte content. For more information see [Docker's api specs](https://docs.docker.com/registry/spec/api/#/content-digests).
A docker image can be referenced using the digest (e.g. `ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950`).
This image refers to a specific image of ubuntu that does not depend on a floating tag.
A workflow containing this Docker image run now and a year from now will run in the exact same container.

However, in order to remove unpredictable behaviors, when Cromwell finds a job ready to be run, it will first look at its docker runtime attribute, and apply the following logic:

* If the job doesn't specify a docker image it will be dispatched and all call caching settings (read/write) will apply normally.
* If the job does specify a docker runtime attribute, then:
    * if the docker image uses a hash, all call caching settings apply normally
    * if the docker image uses a floating tag:
        * Cromwell will attempt to look up the immutable digest of the image with this floating tag. Upon success it will pass both the floating tag and this digest value to the backend.
        * All backends currently included with Cromwell will utilize this digest value to run the job.
        * Within a single workflow, all floating tags within a given workflow will resolve to the same digest value even if Cromwell is restarted when the workflow is running.
        * If Cromwell fails to lookup the digest (for instance an unsupported docker registry, wrong credentials, ...) it will run the job with the user provided floating tag.
        * The actual docker image (floating tag or digest) used for the job will be reported in the `dockerImageUsed` attribute of the call metadata.

**Docker Lookup**

Cromwell provides two methods to lookup a Docker hash from a Docker tag:

* _Local_  
    In this mode, Cromwell will first attempt to find the image on the local machine where it's running using the `docker` CLI. If the image is locally present, then its digest will be used.
    If the image is not present locally, Cromwell will execute a `docker pull` to try and retrieve it. If this succeeds, the newly retrieved digest will be used. Otherwise the lookup will be considered failed.
    Note that Cromwell runs the `docker` CLI the same way a human would. This means two things:
     * The machine Cromwell is running on needs to have Docker installed and a Docker daemon running.
     * The current `docker` CLI credentials on that machine will be used to pull the image.
    
* _Remote_  
    In this mode, Cromwell will attempt to retrieve the hash by contacting the remote docker registry where the image is stored. This currently supports Docker Hub and GCR.
    
    Docker registry and access levels supported by Cromwell for docker digest lookup in "remote" mode:
    
    <!-- Pasted into then regenerated at https://www.tablesgenerator.com/markdown_tables -->

    |               | DockerHub | DockerHub |   GCR  |   GCR   |   ECR  |   ECR   |   ACR  |   ACR   |
    |:-------------:|:---------:|:---------:|:------:|:-------:|:------:|:-------:|:------:|:-------:|
    |               |   Public  |  Private  | Public | Private | Public | Private | Public | Private |
    | Pipelines API |     X     |     X     |    X   |    X    |        |         |        |         |
    |   AWS Batch   |     X     |           |    X   |         |        |         |        |         |
    |      BCS      |           |           |        |         |        |         |        |    X    |
    |     Other     |     X     |           |    X   |         |        |         |        |         |

    <!-- Pasted then regenerated at https://www.tablesgenerator.com/markdown_tables -->

**Runtime Attributes**

As well as call inputs and the command to run, call caching considers the following [runtime
attributes](../RuntimeAttributes.md) of a given task when determining whether to call cache:

* [`ContinueOnReturnCode`](../RuntimeAttributes.md#continueonreturncode)
* [`Docker`](../RuntimeAttributes.md#docker)
* [`FailOnStderr`](../RuntimeAttributes.md#failonstderr)

If any of these attributes have changed from a previous instance of the same task, that instance will not be call-cached
from. Other runtime attributes, including [`memory`](../RuntimeAttributes.md#memory),
[`cpu`](../RuntimeAttributes.md#cpu), and [`disks`](../RuntimeAttributes.md#disks), are not considered by call caching
and therefore may be changed without preventing a cached result from being used.
# Workflow Options Overview

Workflow options can affect the execution of a single workflow without having to change configuration options or restart Cromwell. 

You provide workflow options to Cromwell in a JSON format. This can be supplied at workflow-submit time either via the [CLI](../CommandLine.md) or the [REST endpoint](../api/RESTAPI.md):

```json
{
	"option_name_1": "option value 1",
	"option_name_2": "option value 2"
}
```

Unless otherwise specified you can expect workflow options to override any hard-coded defaults in Cromwell or defaults provided in the [configuration file](../Configuring.md), but to be overridden by any values provided in the workflow definition file itself (WDL or CWL).

Some workflow options apply only to tasks running on the [Google Pipelines API backend](Google).

# Global Workflow Options 

The following workflow options apply to all workflows and their calls regardless of the backend being used.

## Runtime Attributes

Some options allow you to override or set defaults for runtime attributes.

### Setting Default Runtime Attributes

You can supply a default for any [Runtime Attributes](../RuntimeAttributes.md) by adding a `default_runtime_attributes` map to your workflow options file. Use the key to provide the attribute name and the value to supply the default. 

These defaults replace any defaults in the Cromwell configuration file but are themselves replaced by any values explicitly provided by the task in the WDL or CWL file.

Example `options.json`:
```json
{
    "default_runtime_attributes": {
        "docker": "ubuntu:latest",
        "continueOnReturnCode": [4, 8, 15, 16, 23, 42]
    }
}
```

In this example, if a task in a workflow specifies a `docker:` attribute, the task will get what it specifies. However if any task does not provide a value then it will be treated as though it had specified `ubuntu:latest`.

### Specific Runtime Attributes

|Option|Value|Description|
|---|---|---|
|`continueOnReturnCode`|`true` or `false` or integer array|Globally overrides the `continueOnReturnCode` [runtime attribute](../RuntimeAttributes.md) for all tasks| 
|`backend`|An [available](../Configuring.md) backend|Set the **default** backend specified in the Cromwell configuration for this workflow only.|

Example `options.json`:
```json
{
    "continueOnReturnCode": false,
    "backend": "Local"
}
```

In this example, all tasks will be given to the `Local` backend unless they provide a value explicitly in their `runtime { ... }` block. In addition, the `continueOnReturnCode` value for all tasks is hard-coded to `false`, regardless of what the tasks put in their `runtime` block. **TODO or is just a default ala `default_runtime_attributes`?**

## Workflow Failure

The `workflow_failure_mode` option can be given the following values. This overrides any default set by the `workflow-options.workflow-failure-mode` [configuration](../Configuring.md) options.

|Value|Description|
|---|---|
|`ContinueWhilePossible`|Continues to start and process calls in the workflow, as long as they did not depend on the failing call.|
|`NoNewCalls`|No *new* calls are started but existing calls are allowed to finish.|

Example `options.json`:
```json
{
    "workflow_failure_mode": "ContinueWhilePossible"
}
```

## Output Copying
|Option|Value|Description|
|---|---|---|
|`final_workflow_outputs_dir`|A directory available to Cromwell|Specifies a path where final workflow outputs will be written. If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.|
|`use_relative_output_paths`| A boolean | When set to `true` this will copy all the outputs relative to their execution directory. my_final_workflow_outputs_dir/~~MyWorkflow/af76876d8-6e8768fa/call-MyTask/execution/~~output_of_interest . Cromwell will throw an exception when this leads to collisions. When the option is not set it will default to `false`.|
|`final_workflow_log_dir`|A directory available to Cromwell|Specifies a path where per-workflow logs will be written. If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted.|
|`final_call_logs_dir`|A directory available to Cromwell|Specifies a path where final call logs will be written.  If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.|

Note that these directories should be using the same filesystem as the workflow. Eg if you run on Google's PAPI, you should provide `gs://...` paths.

Example `options.json`:
```json
{
    "final_workflow_outputs_dir": "/Users/michael_scott/cromwell/outputs",
    "use_relative_output_paths": true,
    "final_workflow_log_dir": "/Users/michael_scott/cromwell/wf_logs",
    "final_call_logs_dir": "/Users/michael_scott/cromwell/call_logs"
}
```

With `"use_relative_output_paths": false` (the default) the outputs will look like this

```
final_workflow_outputs_dir/my_workflow/ade68a6d876e8d-8a98d7e9-ad98e9ae8d/call-my_one_task/execution/my_output_picture.jpg
final_workflow_outputs_dir/my_workflow/ade68a6d876e8d-8a98d7e9-ad98e9ae8d/call-my_other_task/execution/created_subdir/submarine.txt
```

The above result will look like this when `"use_relative_output_paths": true`:
```
final_workflow_outputs_dir/my_output_picture.jpg
final_workflow_outputs_dir/created_subdir/submarine.txt
```

This will create file collisions in `final_workflow_outputs_dir` when a workflow is run twice. When cromwell
detects file collisions it will throw an error and report the workflow as failed.

## Call Caching Options

These options can override Cromwell's configured call caching behavior for a single workflow. See the [Call Caching](../cromwell_features/CallCaching.md) section for more details and how to set defaults. The call caching section will also explain how these options interact when, for example, one is set `true` and the other is `false`.

**Note:** If call caching is disabled, these options will be ignored and the options will be treated as though they were `false`.

|Option|Values|Description|
|---|---|---|
|`write_to_cache`|`true` or `false`|If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](../cromwell_features/CallCaching.md) section for more details.|
|`read_from_cache`|`true` or `false`|If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](../cromwell_features/CallCaching.md) section for more details.|

Example `options.json`:
```json
{
    "write_to_cache": true,
    "read_from_cache": true
}
```

## Retry with More Memory Multiplier

The `memory_retry_multiplier` workflow option sets the factor by which the memory should be multiplied while retrying 
when Cromwell encounters one of the error keys (specified in Cromwell config using `system.memory-retry-error-keys`) in 
the `stderr` file. The factor should be in the range `1.0 ≤ multipler ≤ 99.0` (note: if set to `1.0` the task will retry 
with same amount of memory). If this is not specified, Cromwell will retry the task (if applicable) but not change the 
memory amount. See the [Retry with More Memory](../cromwell_features/RetryWithMoreMemory.md) section for 
more details.

Example `options.json`:
```json
{
    "memory_retry_multiplier" : 1.1
}
```
# Google Pipelines API Workflow Options

These workflow options provide Google-specific information for workflows running tasks on the Google PAPI backend.

<!-- Pasted into then regenerated at https://www.tablesgenerator.com/markdown_tables -->

| Keys                               | Possible Values | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|------------------------------------|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `jes_gcs_root`                     | `string`        | Where outputs of the workflow will be written.  Expects this to be a GCS URL (e.g. `gs://my-bucket/workflows`).  If this is not set, this defaults to the value within `backend.jes.config.root` in the [Configuration](../Configuring).                                                                                                                                                                                                                                                                      |
| `google_compute_service_account`   | `string`        | Alternate service account to use on the compute instance (e.g. `my-new-svcacct@my-google-project.iam.gserviceaccount.com`).  If this is not set, this defaults to the value within `backend.jes.config.genomics.compute-service-account` in the [Configuration](../Configuring) if specified or `default` otherwise.                                                                                                                                                                                          |
| `google_project`                   | `string`        | Google project used to execute this workflow.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `auth_bucket`                      | `string`        | A GCS URL that only Cromwell can write to.  The Cromwell account is determined by the `google.authScheme` (and the corresponding `google.userAuth` and `google.serviceAuth`). Defaults to the the value in [jes_gcs_root](#jes_gcs_root).                                                                                                                                                                                                                                                                     |
| `monitoring_script`                | `string`        | Specifies a GCS URL to a script that will be invoked prior to the user command being run.  For example, if the value for monitoring_script is `"gs://bucket/script.sh"`, it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.                                                                                                                                                                                                         |
| `monitoring_image`                 | `string`        | Specifies a Docker image to monitor the task. This image will run concurrently with the task container, and provides an alternative mechanism to `monitoring_script` (the latter runs *inside* the task container). For example, one can use `quay.io/broadinstitute/cromwell-monitor`, which reports cpu/memory/disk utilization metrics to [Stackdriver](https://cloud.google.com/monitoring/).                                                                                                             |
| `monitoring_image_script`          | `string`        | Specifies a GCS URL to a script that will be invoked on the container running the `monitoring_image`.  This script will be invoked instead of the ENTRYPOINT defined in the `monitoring_image`.  Unlike the `monitoring_script` no files are automatically de-localized.                                                                                                                                                                                                                                      |
| `google_labels`                    | `object`        | An object containing only string values. Represent custom labels to send with PAPI job requests. Per the PAPI specification, each key and value must conform to the regex `[a-z]([-a-z0-9]*[a-z0-9])?`.                                                                                                                                                                                                                                                                                                       |
| `enable_ssh_access`                | `boolean`       | If set to true, will enable SSH access to the Google Genomics worker machines. Please note that this is a community contribution and is not officially supported by the Cromwell development team.                                                                                                                                                                                                                                                                                                            |
| `delete_intermediate_output_files` | `boolean`       | **Experimental:** Any `File` variables referenced in call `output` sections that are not found in the workflow `output` section will be considered an intermediate `File`. When the workflow finishes and this option is set to `true`, all intermediate `File` objects will be deleted from GCS. Cromwell must be run with the configuration value `system.delete-workflow-files` set to `true`. The default for both values is `false`. NOTE: The behavior of this option on other backends is unspecified. |
| `enable_fuse`                      | `boolean`       | Specifies if workflow tasks should be submitted to Google Pipelines with an additional `ENABLE_FUSE` flag. It causes container to be executed with `CAP_SYS_ADMIN`. Use it only for trusted containers. Please note that this is a community contribution and is not officially supported by the Cromwell development team.                                                                                                                                                                                   |

<!-- Pasted into then regenerated at https://www.tablesgenerator.com/markdown_tables -->

# Example
```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "google_compute_service_account": "my-new-svcacct@my-google-project.iam.gserviceaccount.com",
  "auth_bucket": "gs://my-auth-bucket/private",
  "monitoring_script": "gs://bucket/monitoring_script.sh",
  "monitoring_image": "quay.io/broadinstitute/cromwell-monitor",
  "monitoring_image_script": "gs://bucket/monitoring_image_script.sh",
  "enable_ssh_access": false,
  "google_labels": {
    "custom-label": "custom-value"
  },
  "delete_intermediate_output_files": false,
  "enable_fuse": false
}
```
**Warning!**

 - Cromwell is NOT on its own a security appliance!
 - Only YOU are responsible for your own security! 
 - Please be sure to check with your security team before setting up your Cromwell server
 - Some recommendations and suggestions on security can be found below

__This is intended as helpful guidance only, and is not endorsed by the Broad Institute.__

Cromwell running in server mode accepts all connections on the configured webservice port.  Without taking additional measures to protect your Cromwell server, this can leave your Cromwell server and therefore all information stored by and accessible by your Cromwell server vulnerable to anyone who is able to access the server.  For instance, an unprotected Cromwell server running against the [Google Cloud backend](backends/Google) would leave you vulnerable to outside users running workflows that will cost you money, and potentially access the data files you use within your workflows!

The simplest way to restrict access is by putting an authenticating proxy server in between users and the Cromwell server:  

1. Configure a firewall rule on the Cromwell server host to deny access to the webservice port (e.g. 8000) from all addresses except a secure proxy host.  
2. Configure `<YourFavoriteWebProxy>` on the proxy host with `<YourFavoriteAuthMechanism>`, to proxy authenticated traffic from the world to the Cromwell server. Using Apache `httpd` web server for example with basic `htpassword` file-based authentication, the configuration might look something like:

```Apache
<Location /cromwell>
    Order deny,allow
    Allow from all
    AuthType Basic
    AuthName "Password Required"
    AuthUserFile /path/to/my/htpasswdfile
    Require user someone someoneelse
    ProxyPass http://101.101.234.567:8000   # address of cromwell server web service
</Location>
```

Users now hit `http://my.proxy.org/cromwell` with authenticated requests, and they're forwarded to port 8000 on the Cromwell server host. 

**Multiple Servers on one Host**

The above scheme extends easily to multiple Cromwell instances, for use by different groups within an organization, for example. If multiple instances are running on the same host, then each instance should be run as its own dedicated user (such as service account when using a [Google Cloud backend](backends/Google)), e.g. `cromwell1`, `cromwell2` etc so that processes running under one Cromwell instance cannot access the files of another. Different webservice ports must also be configured separately in order to not clash. If persistent database storage is being used, then each instance must be configured with its own database. The proxy configuration above is extended simply by adding another `Location`:

```Apache
<Location /cromwell1>
    Order deny,allow
    Allow from all
    AuthType Basic
    AuthName "Password Required"
    AuthUserFile /path/to/my/htpasswdfile1
    Require user stillanotherperson andanother
    ProxyPass http://101.101.234.567:8001
</Location>
```

**Multiple Tenants in one Cromwell Server**

Even with a proxy in place, a single Cromwell server does not provide _authorization_ for individual users hitting the endpoints.  This allows, for instance, one _authenticated_ user to view the metadata for a workflow run by another _authenticated_ user.  Due to this limitation, it is important that all users of this Cromwell server be trusted.  

With the exception of the use of the `google_compute_service_account` scheme in the [Configuration](Configuring#authentication) for the [Google Cloud backend](backends/Google)), Cromwell servers are setup to access files for use in workflows using a single set of credentials.  This means that all users of these Cromwell servers will have access to the same files that any other user has access to.  The `refresh_token` or `google_compute_service_account` scheme is the only way to ensure data is protected among multiple users of a Cromwell server, however the aforementioned caveats of authorization for endpoints still applies.  

**Protecting Secrets**

Various parts of the Cromwell server [Configuration](Configuring) contain sensitive information, e.g. username and password for your persistent database.  It is strongly recommended to protect the configuration files from any untrusted users, for instance by limiting who can access your Cromwell server host or using a technology such as [HashiCorp Vault](https://www.vaultproject.io/).  

Additionally, the contents of the Cromwell database can contain sensitive information, so it is recommended to limit the access of the database only to trusted users.# How CWL is Parsed

The ultimate entry point to parsing CWL is in `CwlV1_0LanguageFactory.validateNamespace`. 

Within that method, the process of parsing CWL can be thought of as performing these 4 steps:

1. Produce a canonical one-file representation of the CWL as a JSON object in memory
    - For the entrypoint into this process, see `CwlPreprocessor.preProcessCwl`
2. Use Circe to turn that representation into scala case classes
3. Convert the case classes into an appropriate WOM representations
4. Mix in inputs to produce a WOM Executable.

## Getting a canonical flat file

### From a source file and dependencies zip:

Note that the recursive implementation of `SALAD and flatten` in code can be found in: `CwlPreprocessor.preProcessCwl`.

The process looks like:

1. Write the workflow source file to disk
2. Make sure the dependencies are unzipped next to it
3. `SALAD and flatten` (see `CwlCanonicalizer.getCanonicalCwl`):
    - Run `cwltool --print-pre` against the source file to get a canonicalize JSON representation (aka `SALAD`) it.
        - NB: `cwltool --print-pre` must be able to resolve dependencies but it does not flatten them into the file.
    - Recursively `SALAD and flatten` any references into their own JSON representations of the contents.
    - Replace every reference in-place with the expanded JSON content.
4. At the end of this process we will have:
    - A single flat JSON 
    - Every imported step expanded to contain an in-line description
    - No more external references - everything required is provided in-line

### From a workflow URL:

The process looks similar to "source file and dependencies" except we never write anything locally:

1. `SALAD and flatten`:
    - Run `cwltool --print-pre` against the URL to get a canonicalize JSON representation (aka `SALAD`) it.
        - NB: `cwltool --print-pre` must be able to resolve dependencies but it does not flatten them into the file.
    - Recursively `SALAD and flatten` any references into their own JSON representations of the contents.
    - Replace every reference in-place with the expanded JSON content.
2. At the end of this process we will have:
    - A single flat JSON 
    - Every imported step expanded to contain an in-line description
    - No more external references - everything required is provided in-line


## Using Circe to get case classes

This process is more or less automatic - the case classes define the fields which they anticipate existing in the JSON.
and Circe does the rest.
    - Optional fields in case classes are allowed to not exist in the JSON
    - All fields in the JSON must be represented by fields in the appropriate case classes
    - In cases where the JSON may be structured in one of several ways, we use Shapeless coproducts to specify "one of" many options.

If a file is failing to parse the first stop should always be to double-check that the case classes accurately represent
the range of possible CWL JSON.

## Convert the case classes into WOM representations

- This is a recursive process of examining the CWL case classes and creating WOM values instead.
- Note that the CWL Language Factory (at least as of November 2018) does not implement the provision of WOM Bundles for 
imports. Instead, CWL Language Factory implements only the all-in-one method of converting CWL source and inputs files 
together into a WOM Executable.
    - Pragmatically, this just means that using the CWL language factory to process an import statement in a WDL file 
    is currently impossible (ie we cannot import CWL from WDL... yet).
# Database Reference


## All tables

|Table|Description|
|--|--|
|Call Caching Aggregation| |
|Call Caching Detritus| |
|Call Caching| |
|Call Caching Hash| |
|Call Caching Simpleton| |
|Custom Label| |
|Docker Hash Store| |
|Job Key Value| |
|Job Store| |
|Job Store Simpleton| |
|Metadata Entry| |
|Sub Workflow Store| |
|Summary Status| |
|Workflow Metadata Summary| |
|Workflow Store| |

## Call Caching Aggregation

| Field                             | Type         | Null | Key | Default | Extra          |
|-----------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_AGGREGATION_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| CALL_CACHING_ENTRY_ID             | int(11)      | NO   | MUL | NULL    |                |
| BASE_AGGREGATION                  | varchar(255) | NO   | MUL | NULL    |                |
| INPUT_FILES_AGGREGATION           | varchar(255) | YES  |     | NULL    |                |

## Call Caching Detritus

| Field                          | Type         | Null | Key | Default | Extra          |
|--------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_DETRITUS_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| DETRITUS_KEY                   | varchar(255) | YES  |     | NULL    |                |
| DETRITUS_VALUE                 | longtext     | YES  |     | NULL    |                |
| CALL_CACHING_ENTRY_ID          | int(11)      | YES  | MUL | NULL    |                |

## Call Caching Entry

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_ENTRY_ID     | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | YES  | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| RETURN_CODE               | int(11)      | YES  |     | NULL    |                |
| ALLOW_RESULT_REUSE        | tinyint(4)   | YES  |     | 1       |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |

## Call Caching Hash Entry

| Field                      | Type         | Null | Key | Default | Extra          |
|----------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_HASH_ENTRY_ID | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| HASH_KEY                   | varchar(255) | NO   |     | NULL    |                |
| HASH_VALUE                 | varchar(255) | NO   |     | NULL    |                |
| CALL_CACHING_ENTRY_ID      | int(11)      | YES  | MUL | NULL    |                |

## Call Caching Simpleton Entry

| Field                           | Type         | Null | Key | Default | Extra          |
|---------------------------------|--------------|------|-----|---------|----------------|
| CALL_CACHING_SIMPLETON_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| HASH_KEY                        | varchar(255) | NO   |     | NULL    |                |
| HASH_VALUE                      | varchar(255) | NO   |     | NULL    |                |
| CALL_CACHING_ENTRY_ID           | int(11)      | YES  | MUL | NULL    |                |

## Custom Labels

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| CUSTOM_LABEL_ENTRY_ID   | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| CUSTOM_LABEL_KEY        | varchar(255) | YES  | MUL | NULL    |                |
| CUSTOM_LABEL_VALUE      | varchar(255) | YES  |     | NULL    |                |
| WORKFLOW_EXECUTION_UUID | varchar(100) | NO   | MUL | NULL    |                |

## Docker hash store

| Field                      | Type         | Null | Key | Default | Extra          |
|----------------------------|--------------|------|-----|---------|----------------|
| DOCKER_HASH_STORE_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID    | varchar(255) | NO   | MUL | NULL    |                |
| DOCKER_TAG                 | varchar(255) | NO   |     | NULL    |                |
| DOCKER_HASH                | varchar(255) | NO   |     | NULL    |                |


## Job Key Value

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| JOB_KEY_VALUE_ENTRY_ID    | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | NO   | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |
| STORE_KEY                 | varchar(255) | NO   |     | NULL    |                |
| STORE_VALUE               | varchar(255) | NO   |     | NULL    |                |

## Sub Workflow Store

| Field                          | Type         | Null | Key | Default | Extra          |
|--------------------------------|--------------|------|-----|---------|----------------|
| SUB_WORKFLOW_STORE_ENTRY_ID    | int(11)      | NO   | PRI | NULL    | auto_increment |
| ROOT_WORKFLOW_ID               | int(11)      | NO   | MUL | NULL    |                |
| PARENT_WORKFLOW_EXECUTION_UUID | varchar(255) | NO   | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME      | varchar(255) | NO   |     | NULL    |                |
| CALL_INDEX                     | int(11)      | NO   |     | NULL    |                |
| CALL_ATTEMPT                   | int(11)      | NO   |     | NULL    |                |
| SUB_WORKFLOW_EXECUTION_UUID    | varchar(255) | NO   |     | NULL    |                |

## Summary Status

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| SUMMARY_STATUS_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| SUMMARY_TABLE_NAME      | varchar(255) | NO   | MUL | NULL    |                |
| SUMMARIZED_TABLE_NAME   | varchar(255) | NO   |     | NULL    |                |
| MAXIMUM_ID              | bigint(20)   | NO   |     | NULL    |                |

## Job Store

| Field                     | Type         | Null | Key | Default | Extra          |
|---------------------------|--------------|------|-----|---------|----------------|
| JOB_STORE_ENTRY_ID        | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID   | varchar(255) | YES  | MUL | NULL    |                |
| CALL_FULLY_QUALIFIED_NAME | varchar(255) | YES  |     | NULL    |                |
| JOB_INDEX                 | int(11)      | YES  |     | NULL    |                |
| JOB_ATTEMPT               | int(11)      | YES  |     | NULL    |                |
| JOB_SUCCESSFUL            | tinyint(4)   | NO   |     | NULL    |                |
| RETURN_CODE               | int(11)      | YES  |     | NULL    |                |
| EXCEPTION_MESSAGE         | longtext     | YES  |     | NULL    |                |
| RETRYABLE_FAILURE         | tinyint(4)   | YES  |     | NULL    |                |


## Job Store Simpleton

| Field                        | Type         | Null | Key | Default | Extra          |
|------------------------------|--------------|------|-----|---------|----------------|
| JOB_STORE_SIMPLETON_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| SIMPLETON_KEY                | varchar(255) | NO   |     | NULL    |                |
| SIMPLETON_VALUE              | longtext     | YES  |     | NULL    |                |
| WDL_TYPE                     | varchar(255) | NO   |     | NULL    |                |
| JOB_STORE_ENTRY_ID           | int(11)      | YES  | MUL | NULL    |                |


## Metadata

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| METADATA_JOURNAL_ID     | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID | varchar(255) | NO   | MUL | NULL    |                |
| METADATA_KEY            | varchar(255) | NO   |     | NULL    |                |
| CALL_FQN                | varchar(255) | YES  |     | NULL    |                |
| JOB_SCATTER_INDEX       | int(11)      | YES  |     | NULL    |                |
| JOB_RETRY_ATTEMPT       | int(11)      | YES  |     | NULL    |                |
| METADATA_VALUE          | longtext     | YES  |     | NULL    |                |
| METADATA_TIMESTAMP      | datetime     | NO   |     | NULL    |                |
| METADATA_VALUE_TYPE     | varchar(10)  | YES  |     | NULL    |                |


## Workflow Metadata Summary

| Field                              | Type         | Null | Key | Default | Extra          |
|------------------------------------|--------------|------|-----|---------|----------------|
| WORKFLOW_METADATA_SUMMARY_ENTRY_ID | bigint(20)   | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID            | varchar(100) | NO   | UNI | NULL    |                |
| WORKFLOW_NAME                      | varchar(100) | YES  | MUL | NULL    |                |
| WORKFLOW_STATUS                    | varchar(50)  | YES  | MUL | NULL    |                |
| START_TIMESTAMP                    | datetime     | YES  |     | NULL    |                |
| END_TIMESTAMP                      | datetime     | YES  |     | NULL    |                |
| SUBMISSION_TIMESTAMP               | datetime     | YES  |     | NULL    |                |

## Workflow Store

| Field                   | Type         | Null | Key | Default | Extra          |
|-------------------------|--------------|------|-----|---------|----------------|
| WORKFLOW_STORE_ENTRY_ID | int(11)      | NO   | PRI | NULL    | auto_increment |
| WORKFLOW_EXECUTION_UUID | varchar(255) | YES  | UNI | NULL    |                |
| WORKFLOW_DEFINITION     | longtext     | YES  |     | NULL    |                |
| WORKFLOW_INPUTS         | longtext     | YES  |     | NULL    |                |
| WORKFLOW_OPTIONS        | longtext     | YES  |     | NULL    |                |
| WORKFLOW_STATE          | varchar(20)  | YES  | MUL | NULL    |                |
| SUBMISSION_TIME         | datetime     | NO   |     | NULL    |                |
| IMPORTS_ZIP             | longblob     | YES  |     | NULL    |                |
| CUSTOM_LABELS           | longtext     | NO   |     | NULL    |                |
| WORKFLOW_TYPE           | varchar(30)  | YES  |     | NULL    |                |
| WORKFLOW_TYPE_VERSION   | varchar(255) | YES  |     | NULL    |                |
| WORKFLOW_ROOT           | varchar(100) | YES  |     | NULL    |                |
| CROMWELL_ID             | varchar(100) | YES  |     | NULL    |                |
| HEARTBEAT_TIMESTAMP     | timestamp    | YES  |     | NULL    |                |

## Cromwell Backend Development

### What’s a Backend?

Cromwell is a workflow execution service, which means it takes a representation of work to do (tasks) as well as the
dependencies and code flow between them. When it comes to actually executing the specific task, Cromwell delegates that
work to a “backend”. A backend is therefore responsible for the actual execution of a single task on some underlying
computing platform, such as Google Cloud Platform, AWS, TES, Local, etc. Backends are implemented as a software layer
between the Cromwell engine and that underlying platform.

The underlying platform services requests for a job to run while the backend shim provides the interface layer between
Cromwell and the platform. In general, the more sophisticated the platform the thinner the shim needs to be and vice
versa.

A job from Cromwell’s perspective is a collection of the following information:

* A unix command line to run
* A mapping of where its input files currently live to where the command line expects them to be
* A mapping of where the command line will write its outputs to where they should eventually wind up
* An optional Docker image which if supplied will be the environment in which the command will be run
* A collection of arbitrary key/value pairs which can be used by the execution platform to tune the request, e.g. amount
  of memory or the number of CPUs.

The Docker image is not required for a backend but is highly recommended. The bioinformatics workflow field is rapidly
moving towards a Docker model so one will find better support for their backend if it is built around using Docker
containers.

The input/output file mappings will depend on the needs of the platform. In the simplest example these values would be
identical, living on a shared filesystem. However an example of where this would be useful would be a model where an
input file lives in a cloud bucket and is directly copied onto the machine running the command line followed by copying
the program’s outputs back to a cloud bucket.

The underlying platform can be extremely simple or have arbitrary levels of complexity as long as one can map from the
above concepts to running a command. That mapping could happen completely in the platform, the Cromwell backend, or a
mixture of the two. The implementer of a backend will need to strike the balance which works best for both their needs
and the particulars of the platform itself.

Throughout the rest of this document the following terms are used, the difference between them may be subtle but
important:

* Workflow: a container of one or more calls, possibly expressing execution dependencies on each other. May contain
  scatters, conditional logic, or subworkflow invocations.
* Task: The abstract definition of a thing to run. Think of this like a function in a programming language.
* Call: An instantiated request in Cromwell of a thing to run. To further the function analogy this would be an
  invocation of that function.
* Job: The physical manifestation of a call by the backend. Examples of this would be an SGE job, a unix process, or a
  Google Life Sciences operation ID.

### Backend Lifecycle

Within a running workflow, the backend is used in the following manner:

* Initialization: Initialization routine called the first time a workflow uses a backend
* Execute: The workflow requests that the backend run a job
* Recover: Attempt to reconnect to a previously started job. An example of this would be if Cromwell was restarted and
  wanted to reattach to currently running jobs.
* Abort: Request that a running job be halted
* Finalization: When a workflow is complete, allows for any workflow level cleanup to take place.

The initialization and finalization steps are optional and are provided for cases where a backend needs that behavior;
not all backends will need these.

The implementations of both the recover and abort steps are up to the backend developer. For instance the recover
function could be implemented to actually perform an execution and/or the abort function could be written to do nothing
at all. Implementing these in a more robust manner is recommended for most cases but neither are universally
appropriate.

### How do I create a Backend?

This section is assuming that you both know the underlying execution platform you wish to use and that you know how to
programmatically submit work to it.

The Cromwell engine uses backend-specific types extending Akka `Actor` while processing a workflow:

* [`BackendLifecycleActorFactory`](https://github.com/broadinstitute/cromwell/blob/9bf1622ca8988365477b77b9f26ce388b54fc58c/backend/src/main/scala/cromwell/backend/BackendLifecycleActorFactory.scala#L17)
  The entry point into the backend implementation specified in Cromwell configuration. e.g.
  a [sample configuration for the Local backend](https://github.com/broadinstitute/cromwell/blob/2b19f00976ee258142185917083460d724f7fe3d/cromwell.example.backends/cromwell.examples.conf#L370)
* [`BackendWorkflowInitializationActor`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L27)
  Handles the initialization phase
* Usually a derivative
  of [`StandardAsyncExecutionActor`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L75)
  to run an individual job within the workflow
* [`BackendWorkflowFinalizationActor`](https://github.com/broadinstitute/cromwell/blob/a40de672c565c4bbd40f57ff96d4ee520dc2b4fc/backend/src/main/scala/cromwell/backend/BackendWorkflowFinalizationActor.scala#L10)
  Handles the finalization phase

These three actors are all represented by a trait which a backend developer would need to implement. Both the
initialization and finalization actors are optional.

The following explanations of the traits mention multiple types defined in the Cromwell codebase. In particular
`BackendJobDescriptor` wraps all the information necessary for a backend to instantiate a job and `JobKey` provides the
information to uniquely identify a job. It is recommended that one look in the Cromwell codebase for more information on
these and other types.

#### BackendWorkflowInitializationActor

If a backend developer wishes to take advantage of the initialization phase of the backend lifecycle they must implement
this trait. There are three functions which must be implemented:

* [`abortInitialization: Unit`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L168)
   specifies what to do, if anything, when a workflow is requested to abort while a backend initialization is in
   progress.
* [`validate: Future[Unit]`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L178)
   is provided so that the backend can ensure that all of the calls it will handle conform to the rules of that backend.
   For instance, if a backend requires particular runtime attributes to exist.
* [`beforeAll: Future[Option[BackendInitializationData]]`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L173)
   is the actual initialization functionality. If a backend requires any work to be done prior to handling a call, that
   code must be called from here.

#### StandardAsyncExecutionActor

Nearly all production backend implementations in Cromwell extend
the [StandardAsyncExecutionActor](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L75)
trait. The minimum overrides for implementations of this trait are enumerated below, but overrides of other methods will
also be required. There are several backend implementations in the Cromwell codebase that can serve as references, for
example:

* [Google Life Sciences / Pipelines API](https://github.com/broadinstitute/cromwell/blob/0aff35336b4e2ba19b18530a68e622df1462d9b7/supportedBackends/google/pipelines/common/src/main/scala/cromwell/backend/google/pipelines/common/PipelinesApiAsyncBackendJobExecutionActor.scala#L95)
  common layer
    * [Life Sciences API (beta)](https://github.com/broadinstitute/cromwell/blob/a49e1fc65703ccfda2840d1d9266fad2bdbb7339/supportedBackends/google/pipelines/v2beta/src/main/scala/cromwell/backend/google/pipelines/v2beta/PipelinesApiAsyncBackendJobExecutionActor.scala#L27)
    * [Pipelines API (alpha)](https://github.com/broadinstitute/cromwell/blob/a49e1fc65703ccfda2840d1d9266fad2bdbb7339/supportedBackends/google/pipelines/v2alpha1/src/main/scala/cromwell/backend/google/pipelines/v2alpha1/PipelinesApiAsyncBackendJobExecutionActor.scala#L27)
* [GA4GH Task Execution Service (TES)](https://github.com/broadinstitute/cromwell/blob/6bf7af3c12a411db26786ac34646238fc053ec97/supportedBackends/tes/src/main/scala/cromwell/backend/impl/tes/TesAsyncBackendJobExecutionActor.scala#L55)
* [AWS Batch](https://github.com/broadinstitute/cromwell/blob/470d482e8ba2a9e2bc544896a4e6ceea57d55bb2/supportedBackends/aws/src/main/scala/cromwell/backend/impl/aws/AwsBatchAsyncBackendJobExecutionActor.scala#L83)

Overrides required for compilation and basic execution of `StandardAsyncExecutionActor` implementations:

* [`type StandardAsyncRunInfo`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L89)
  encapsulates the type of the run info when a job is started.
* [`type StandardAsyncRunState`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L92)
  encapsulates the type of the run status returned during each poll.
* [`def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L97)
  should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be
  carried around in the state type
* [`def standardParams: StandardAsyncExecutionActorParams`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L103)
  a standard set of parameters passed to the backend.
* [`def isTerminal(runStatus: StandardAsyncRunState): Boolean`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L829)
  Returns true when a job is complete, either successfully or unsuccessfully.
* At least one of:
  * [`def execute(): ExecutionHandle`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L738) executes the job specified in the params
  * [`def executeAsync(): Future[ExecutionHandle]`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L748) asynchronously executes the job specified in the params
* At least one of:
  * [`def pollStatus(handle: StandardAsyncPendingExecutionHandle): StandardAsyncRunState`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L798) returns the run status for the job
  * [`def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[StandardAsyncRunState]`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L808) asynchronously returns the run status for the job

#### BackendWorkflowFinalizationActor

There is only one function to override for this trait if the backend developer chooses to use the finalization
functionality:

* [`afterAll: Future[Unit]`](https://github.com/broadinstitute/cromwell/blob/a40de672c565c4bbd40f57ff96d4ee520dc2b4fc/backend/src/main/scala/cromwell/backend/BackendWorkflowFinalizationActor.scala#L37)
  is the hook to perform any desired functionality, and is called when the workflow is completed.
**Contribution Guidelines**

Cromwell is an open-source project and we heartily welcome community contributions to both our code and our documentation. Here are some guidelines for adding documentation and recommendations on where we could use help the most. 

First off, here are useful links:

* Cromwell documentation: [cromwell.readthedocs.io](http://cromwell.readthedocs.io)
* Source on Github: [github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell/tree/develop/docs)
* Builds on ReadTheDocs: [readthedocs.org/projects/cromwell](https://readthedocs.org/projects/cromwell/builds/)
* How to build and view the documentation locally: [mkdocs.readthedocs.io](https://mkdocs.readthedocs.io/en/stable/#installation)
* Public Google Docs - Table of Contents:  [drive.google.com/open?id=1myTVzWx5HG720nPBUAF9vE8XTacrzCW70RpgsjMnchM](https://drive.google.com/open?id=1myTVzWx5HG720nPBUAF9vE8XTacrzCW70RpgsjMnchM)

### "needs docs"

There are plenty of areas of the Cromwell documentation that need to be updated, improved, or added. Within the Cromwell repo on Github there are many issues labeled as ["needs docs"](https://github.com/broadinstitute/cromwell/issues?q=is%3Aopen+is%3Aissue+label%3A%22needs+docs%22), so feel free to start there. 

If you would like to request additional documentation, you can create a [Github issue in the Cromwell repo](https://github.com/broadinstitute/cromwell/issues/new).

### Writing Tips

1. Keep it clear, accurate, and concise.
2. Put the most important information first.
3. Use the second person, use “you” instead of “the user”.
4. No passive verbs (everything is done by something).
5. Link to the original source, don't repeat documentation. 

### Formatting

The documentation is written in Markdown. Click here for a [Github Guide on Markdown](https://guides.github.com/features/mastering-markdown/), and click here for more [tips from MkDocs](http://www.mkdocs.org/user-guide/writing-your-docs/).

### Styling

**Links:**

* Absolute: `[link text](www.destinationURL.com)`
_Example:_ `[Broad Institute](www.broadinstitute.org)` _produces this link_ [Broad Institute](https://www.broadinstitute.org).
* Relative: `[link text](Destination_Page)`, where `Destination_Page` is the file name without the `.md` extension  
_Example:_ `[How to use the Cromwell CLI](CommandLine)` _produces this link_ [How to use the Cromwell CLI](CommandLine).
* Anchor link: `[anchor text](../Path/To/Page#Anchor)`  
_Example:_ `[HPC filesystems](backends/HPC#filesystems)` _produces this link_ [HPC filesystems](backends/HPC#filesystems).

**Code:**

* To style a word of code, use a backtick (\`) before and after the word.  
_Example:_ \`file.json\` _produces_ `file.json`. 
* To style a block of code, use three backticks (\`\`\`) before and after the block of code.  
_Example:_  
	\`\`\`  
		workflow myWorkflow {  
		call myTask  
		}  
	\`\`\`  
_produces this block_  
```  
	workflow myWorkflow {  
	call myTask  
	}  
```  
* To use syntax highlighting, include the language after the first three backticks (\`\`\`).  
_Example:_   
	\`\`\`json  
		\{  
		    "MyWorkflow.MyTask.VariableTwo": "Variable2"  
		\}  
	\`\`\`	
_produces this block_  
```json  
	{
	    "MyWorkflow.MyTask.VariableTwo": "Variable2"
	}
```	

**Images**

* Relative: `![](ImgName.png)`
	* _Example:_ `![](../jamie_the_cromwell_pig.png)` _produces this image_  
	![](../jamie_the_cromwell_pig.png)  
* Absolute: `![](URLofImg.png)`
	* _Example:_ `![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)` _produces this image_  
	![](https://www.broadinstitute.org/sites/all/themes/custom/at_broad/logo.png)  

### REST API & Menu

**REST API:**

1. Edit the [`cromwell.yaml`](https://github.com/broadinstitute/cromwell/blob/develop/engine/src/main/resources/swagger/cromwell.yaml) to make any changes to the REST API content.  
2. Regenerate the REST API markdown file by running `sbt generateRestApiDocs` from the main Cromwell directory.
3. Commit both the changes to the `cromwell.yaml` and the (re)generated [RESTAPI.md](https://github.com/broadinstitute/cromwell/blob/develop/docs/api/RESTAPI.md).
4. Once your branch is merged to the [`develop` branch](https://github.com/broadinstitute/cromwell/tree/develop), you will see your changes on the [REST API page](api/RESTAPI/).

**Left-side menu:**

To add or remove items from the menu, edit [mkdocs.yml](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml) in Cromwell.

### FAQs

**_Why isn't my documentation showing up?_** 

* **Is your PR merged?**  
If not, [kindly ask the team to merge it](https://github.com/broadinstitute/cromwell/pulls). Once your PR is merged to develop, it will trigger an automatic build.  

* **Has the build finished?**   
[Check build status here](https://readthedocs.org/projects/cromwell/builds/).

* **Did you add the file(s) to the YAML file?**  
If not, [add it here](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml).

**_How do I add to the REST API documentation?_**

* **Don't forget to regenerate**  
After you edit the [`cromwell.yaml`](https://github.com/broadinstitute/cromwell/blob/develop/engine/src/main/resources/swagger/cromwell.yaml), run `sbt generateRestApiDocs` and commit all changes.  
_Hint:_ Once you have regenerated the docs correctly, the hidden timestamp at the top of the [`RESTAPI.md` file](https://raw.githubusercontent.com/broadinstitute/cromwell/develop/docs/api/RESTAPI.md) will show the current time.
Most users should not need to build Cromwell and can use pre-built Cromwell [releases](Getting).

If for some reason you require a non-release version of Cromwell or are developing new Cromwell
features or fixes, the following are required to build Cromwell from source:

* [Scala 2.12](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/)
* [Git](https://git-scm.com/)

You can also use the [development image](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-develop), and build a development container to work inside:

```bash
$ docker build -t cromwell-dev .
$ docker run -it cromwell-dev bash
```

First start by cloning the Cromwell repository from GitHub:

```bash
$ git clone git@github.com:broadinstitute/cromwell.git
```

Next change into the `cromwell` directory:

```bash
$ cd cromwell
```

If you require a specific version of Cromwell as a starting point, do the appropriate `git checkout` now. 

Finally build the Cromwell jar:

```bash
$ sbt assembly
```

NOTE: This command will run for a long time the first time.  
NOTE: Compiling will not succeed on directories encrypted with ecryptfs (ubuntu encrypted home dirs for example), due to long file paths.

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`. It will also build a runnable Womtool JAR in `womtool/target/scala-2.12/` with a name like `womtool-<VERSION>.jar`.

To build a [Docker](https://www.docker.com/) image, run:

```bash
$ sbt server/docker
```

This will build and tag a Docker image with a name like `broadinstitute/cromwell:<VERSION>-SNAP`.
Centaur is an integration testing suite for the [Cromwell](http://github.com/broadinstitute/cromwell) execution engine.  Its purpose is to exercise the functionality of a specific deployment of Cromwell, to ensure that it is functioning properly 'in the wild'.  

## Prerequisites

Centaur expects to find a Cromwell server properly configured and running in server mode, listening on port 8000.  
This can be configured by modifying the `cromwellUrl` parameter in `application.conf`.

You can get a build of your current cromwell code with [these instructions](Building.md).
The server can be run with `java -jar <Cromwell JAR> server`, checkout [this page](../CommandLine.md) 
for more detailed instructions. 
You can now run the tests from another terminal.

## Running

There are two ways to invoke the integration tests:

* `sbt "centaur / IntegrationTest / test"` - compiles and run via sbt directly, simple but also has the problem of running 2x cores tests in parallel which can overwhelm your Cromwell server if running in a development environment

* `src/ci/bin/testCentaurLocal.sh` - runs the same tests using the continuous integration pipeline configuration

### Tags

All tests are tagged with their name and their TESTFORMAT, and also any custom tags specified in the `.test` file.

Tag names are all lower case, so a test named "tagFoo" has a tag "tagfoo".

To run only those tests which have been tagged with a specified tag `tagFoo`:
```
sbt "centaur / IntegrationTest / testOnly * -- -n tagfoo"
```

Or to instead exclude all tests which have been tagged with a specified tag `tagFoo`:
```
sbt "centaur / IntegrationTest / testOnly * -- -l tagfoo"
```

## Adding custom tests

You can add your own tests to the test suite by adding `-Dcentaur.optionalTestPath=DIR` on your sbt invocation, 
e.g. `sbt -Dcentaur.optionalTestPath=/some/path/to/tests test`. The value of `DIR` is expected to be a directory
which contains one or more test case files.
 
The same result can be achieved more permanently by adding the custom directory into the `application.conf` file directly: 
```
centaur {
  optionalTestPath = "/some/path/to/tests"
}
```

## Defining test cases

Each test case file is a HOCON file with the following structure:
```
name: NAME  // Required: Name of the test
testFormat: TESTFORMAT // Required: One of WorkflowSuccessTest, WorkflowFailureTest, runtwiceexpectingcallcaching
backends: [BACKENDNAME1, BACKENDNAME2, ...] // Optional list of backends. If supplied, this test will be ignored if these backends are not supported by the Cromwell server
basePath: /an/optional/field  // Optional, location for the files {} entries to be found relative to
tags: [ "any", "custom", "tags" ]  // Optional, a set of custom tags to apply to this test
ignore: false  // Optional, whether centaur will ignore this test when running

files {
  wdl: path/to/wdl  // Required: path to the WDL file to submit
  inputs: optional/path/to/inputs  // Optional, a path to an inputs JSON to include in the submission
  options: optional/path/to/options  // Optional, a path to an options JSON to include in the submission
}

// Optional, some metadata to verify on workflow completion:
metadata {
  fully.qualified.key.name1: VALUE1
  fully.qualified.key.name2: VALUE2
  // Examples:
  // failures is a list, the first entry (0) might be the error you are looking for. If multiple errors are expected the entire list can be checked. 
  // It has a "message" and a "causedBy" field.
  "failures.0.message": "Cromwell senses you did not use WomTool validate."
  "failures.0.causedBy": "BetweenKeyboardAndChairException"
}

filesystemcheck: "local" // possible values: "local", "gcs". Used in conjunction with outputExpectations to define files we expect to exist after running this workflow.
outputExpectations: {
    "/path/to/my/output/file1": 1
    "/path/to/file/that/should/not/exist": 0
}
```

The tags are optional. If supplied they will allow people to turn on or off this test case by including or excluding tags when running (see above).

The `basePath` field is optional, but if supplied all paths will be resolved from that directory. If it is not supplied, all paths will be resolved from the directory the test case file is in.

The `testFormat` field can be one of the following, case insensitive:
* `workflowsuccess`: The workflow being supplied is expected to successfully complete
* `workflowfailure`: The workflow being supplied is expected to fail

The `metadata` is optional. If supplied, Centaur will retrieve the metadata from the successfully completed workflow and compare the values retrieved to those supplied. At the moment the only fields supported are strings, numbers and booleans.

You can find which metadata is recorded by running a workflow ```java -jar <Cromwell JAR> run -m metadata.json my_workflow.wdl```.
This will save the metadata in `metadata.json`.

For any metadata values or outputExpectations which require workflow ID (i.e, file paths), use `<<UUID>>` as a placeholder instead. For example:
* `"calls.hello.hello.stdout": "gs://google-project/jes/root/wdl/<<UUID>>/call-task/task-stdout.log"`

In case the absolute path the cromwell root is used (for example: `/home/my_user/projects/cromwell/cromwell-executions`)
 you can use `<<WORKFLOW_ROOT>>` as a replacement. 
* `"calls.hello.hello.exit_code": "<<WORKFLOW_ROOT>>/call-hello/execution/exit_code"`

In case testing of the caching is required `<<CACHE_HIT_UUID>>` can be used. 
The testFormat should be `runtwiceexpectingcallcaching`.
# Overview

Cromwell's instrumentation support can be useful to collect utilization data in long-running, high-volume
production environments. The default implementation of this ignores these metrics, but Cromwell includes alternate implementations that can forward metrics to a
specific server.

### StatsD

While this instrumentation support can be used in smaller environments it will still require setting up a
[StatsD](https://github.com/etsy/statsd) server outside of Cromwell and it's possible not enough data would be produced to be useful. 
Cromwell collects metrics while running and sends them to an internal service. 

Make sure to configure your StatsD service:

```hocon
services.Instrumentation {
	class = "cromwell.services.instrumentation.impl.statsd.StatsDInstrumentationServiceActor"

	config {
	    hostname = "localhost" # Replace with your host
	    port = 8125 # Replace with your port
	    # prefix = "my_prefix" # All metrics will be prefixed by this value if present.
	    flush-rate = 1 second # Rate at which metrics are sent to the StatsD server
  	}
}
```

There is also an additional configuration value that can be set: 

```hocon
# Rate at which Cromwell updates its gauge values (number of workflows running, queued, etc...)
system.instrumentation-rate = 5 seconds
```

If you have multiple Cromwell instances, and would like to separate the instrumentation path for each instance, set the `system.cromwell_id` with the unique identifier for each Cromwell instance. For example,
```hocon
system.cromwell_id = "cromwell-instance-1"
```
will prepend all the metrics with path `cromwell.cromwell-instance-1...` for that instance.


##### Metrics

The current StatsD implementation uses metrics-statsd to report instrumentation values.
metrics-statsd reports all metrics with a gauge type.
This means all metrics will be under the gauge section. We might add or remove metrics in the future depending on need and usage.
These are the current high level categories:

* `backend`
* `rest-api`
* `job`
* `workflow`
* `io`


### Stackdriver

Cromwell now supports sending metrics to [Google's Stackdriver API](https://cloud.google.com/monitoring/api/v3/). To use the Stackdriver instrumentation
specify this in your config:
```hocon
services.Instrumentation {
    class = "cromwell.services.instrumentation.impl.stackdriver.StackdriverInstrumentationServiceActor"

    config {
        # auth scheme can be `application_default` or `service_account`
        auth = "service-account"
        google-project = "my-project"
        # rate at which aggregated metrics will be sent to Stackdriver. It needs to be equal or greater than 1 minute.
        # Google's Stackdriver API needs each metric to be sent not more than once per minute.
        flush-rate = 1 minute
        # below 3 keys are attached as labels to each metric. `cromwell-perf-test-case` is specifically meant for perf env.
        cromwell-instance-role = "role"
        cromwell-perf-test-case = "perf-test-1"
    }
 }
```
The 2 label keys are optional. If specified, each metric will have label(s) added in the form of a (key, value) pair.
So for example, if `cromwell-instance-role = "backend"` is mentioned in config, each metric data point sent to Stackdriver
will have a label (cromwell_instance_role, backend) added to it.

There is another optional label that can be added to each metric. `cromwell_id` represents the identifier for different Cromwell instances.
```hocon
# Unique Cromwell instance identifier
system.cromwell_id = "cromwell-instance-1"
```


##### Metric type and Label keys naming convention
More details on the this can be found [here](https://cloud.google.com/monitoring/api/v3/metrics-details#metric-kinds).

You must adhere to the following spelling rules for metric type names:
- You can use upper and lower-case letters, digits, and underscores (_) in the names.
- You can use periods (.) in the domain part of the names.
- You can use forward slashes (/) to separate path elements.
- You can start each path element with a letter or digit.
- The maximum length of a metric type name is 200 characters.

You must adhere to the following spelling rules for metric label names:
- You can use upper and lower-case letters, digits, underscores (_) in the names.
- You can start names with a letter or digit.
- The maximum length of a metric label name is 100 characters.



Coming soon
# IoActor: basic concepts

* **Word count:** 225

## Actor Hierarchy

IO subsystem consists of the following actors.

![IoActor hierarchy](IoActor_basic_hierarchy.png)

All messages are sent to `IoActorProxy`, which then routes them either to `IoActor` or `IoPromiseProxyActor` based on 
whether the message contains a Promise. `IoPromiseProxyActor` receives a message with Promise, forwards message to the 
`IoActor`, and once response from `IoActor` comes back, completes the promise.

`IoActor` has back-pressure (always enabled) and throttling (configurable) features implemented in order to prevent 
overflowing the incoming messages queue.

## Message types and processing logic

`IoActor` accepts messages of base type `IoCommand[T]` with or without clientContext (which is a Promise, which will be 
passed through, back to the `IoPromiseProxyActor` in the end)

Message processing logic is implemented using Akka Stream API and can be represented with the following graph:
![IoActor message processing graph](IoActor_message_processing_graph.png)

There are 2 message processing flows:  
1. "GCS batch flow" to process messages, which are subtypes of `GcsBatchIoCommand`  
1. "Default flow" for all other subtypes of `IoCommand`
  
Both "GCS batch flow" and "Default flow" allow to configure desired parallelism level, which can be done during 
`IoActor's` creation using `gcsParallelism` and `nioParallelism` parameters respectively. 
The particular implementation logic for each base IO command is encapsulated inside `ParallelGcsBatchFlow` and 
`NioFlow` classes.

The final `IoResult` is sent both to the "reply sink" which is responsible for sending result to the original sender, and 
"instrumentation sink", which is responsible for sending `InstrumentationServiceMessage` to the `ServiceRegistryActor`  
# Workflow Execution: Execution Store and Value Store Examples

## Introduction

This page provides run-throughs to give insight into how the
[Execution Store](executionStore.md) and [Value Store](valueStore.md) work in
practice in some example situations. 

## Handling a single task call

To begin, consider this simple workflow. It has a single task call whose
result is exposed as an output String:

```wdl
version 1.0

workflow single_task_workflow {
  call single_task
  
  output {
    String string_out = single_task.string_out
  }
}

task single_task {
  command {
    echo hello
  }
  output {
    String string_out = "hello"
  }
}
```

The **Execution Store** will keep track of statuses as the workflow runs:

| | `single_task` | `string_out` |
|---|---|---|
|1|`NotStarted`|`NotStarted`|
|2|`QueuedInCromwell`|`NotStarted`|
|3|`Starting`|`NotStarted`|
|4|`Running`|`NotStarted`|
|5|`Done`|`NotStarted`|
|6|`Done`|`Running`|
|7|`Done`|`Done`|

In step 1, the workflow has just started and the ExecutionStore is created in its initial
state. The Value Store doesn't track statuses and so begins empty: `{ }`.

In steps 2-4, the Execution Store tracks the `single_task` job as the engine is executing it.

As the Execution Store is updated to indicate task completion is step 5, the Value Store is also updated to
include the output value of the task:
```json
{
  "single_task.string_out": "hello"
}
```

By step 6, Cromwell can use the fact that the task is complete to decide that the output node is ready to be
evaluated. And the input to the output expression is available for lookup in the Value Store.

In step 7, all workflow nodes have run and the workflow is complete. The Value Store is updated once again to
additionally contain the output node value:
```json
{
  "single_task.string_out": "hello",
  "string_out": "hello"
}
``` 

Cromwell can use this information to trigger the "workflow complete" logic.

## Handling scatters

When Cromwell runs scattered tasks, the Execution Store cannot tell ahead of time how many
`JobKey`s it will need to represent all of the shards in the scatter. It can get around
this problem by putting a placeholder `JobKey` for the scatter node in the Execution Store. When
the scatter key is evaluated, it expands the Execution Store to include new `JobKey`s representing
every shard in the scatter.

As with the single task example, the Value Store starts empty, and is updated with the results of each
shard only as and when they are generated.

To see that in action, Consider this workflow:

```wdl
version 1.0

workflow scattered_task_workflow {
  scatter (x in range(2)) {
    call scattered_task
  }
  output {
    Int results_count = length(scattered_task.string_out)
  }
}

task scattered_task {
  command {
    echo hello
  }
  output {
    String string_out = "hello"
  }
}
```

#### Scatter Expansion

As the workflow starts, the execution store has three entries. An `x` represents the array-input for the scatter,
a`ScatterNode` represents the placeholder for expanding the scatter, and a `results_count` represents the workflow output. 

The start of workflow execution looks like this:

| | `x` | `ScatterNode` | `results_count` |
|---|---|---|---|
|1|`NotStarted`|`NotStarted`| `NotStarted` |
|2|`Running`|`NotStarted`| `NotStarted` | 
|3|`Done`|`NotStarted`| `NotStarted` | 

Once `x` is evaluated the value store gains an entry:
```json
{
  "x": "[0, 1]"
}
```

The scatter node now becomes runnable because its upstream dependency (`x`) is `Done` in the Execution Store.

The evaluation of `ScatterNode` updates the execution store in a number of ways:

* One call key for each index of `scattered_task` is added.
* The `scattered_task` gets an un-indexed key too. This key is used to mark when all of the shards of the call are complete.
* The gathered value `scattered_task.string_out` represents the "gathered" results of the task's output. It only runs 
once the un-indexed `scattered_task` key is Done and gathers output values into an array.
This gather key also acts as the upstream dependency of the `results_count` output expression.
* The `ScatterNode` is marked as `Done` so that it doesn't get triggered to run again.

Following the scatter-expansion evaluation of `ScatterNode`, the Execution Store looks like this:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|4|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`

The Value Store is not changed at this time because no new values have been generated.

#### Parallel Shard Execution

The two scattered shards are now immediately runnable because they have no upsteam dependencies.
As the two jobs are run, the Execution Store map updates to track their statuses:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|4|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`
|5|`Done`|`Done`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`
|6|`Done`|`Done`|`QueuedInCromwell`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`
|7|`Done`|`Done`|`Starting`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`
|8|`Done`|`Done`|`Starting`|`Starting`|`NotStarted`|`NotStarted`|`NotStarted`
|9|`Done`|`Done`|`Running`|`Starting`|`NotStarted`|`NotStarted`|`NotStarted`
|10|`Done`|`Done`|`Running`|`Running`|`NotStarted`|`NotStarted`|`NotStarted`
|11|`Done`|`Done`|`Running`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`
|12|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`

As the results for each shard come in, the value store is also updated to include them:

At step 11 (shard 1 has finished but shard 0 has not):
```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:1": "hello"
}
```

At step 12:
```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello"
}
```

#### Scatter Completion and Gathering

Once all of the sharded keys for `scattered_task` are complete, the un-indexed marker key for that call becomes
runnable. And once the marker is complete, the gather key for the output also becomes runnable.

The progression in the Execution Store goes like:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|12|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`
|13|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`
|14|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`

As the gather node completes in step 14, the value store is also updated to contain the unindexed, gathered result of 
the `scattered_task.string_out` output:

```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello",
  "scattered_task.string_out": ["hello", "hello"]
}
```

When the `scattered_task.string_out` gather node completes, the upstream dependencies of the `results_count` output are
finally satisfied and it becomes runnable too. It runs to produce the workflow outputs:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|14|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`
|15|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`Running`
|16|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`

For step 16, completion of the output evaluation creates an entry in the Value Store which can be exposed as a workflow output as the 
workflow completes:

```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello",
  "scattered_task.string_out": ["hello", "hello"],
  "results_count": 2
}
```### WDL Expression Evaluation

#### Expressions in WOM

Expressions in WOM expose the following methods:

* List the names of inputs which the expression will need in order to evaluate.
* List the files which would need to be available in order to evaluate.
* Evaluate the type of value which the expression will evaluate to.
* Evaluate the expression.  

#### How WDL expressions become WomExpressions

Relating back to the [WDL parsing](../workflowParsing/wdlParsingOverview.md) process:

* WDL 1.0 engine functions become contextless WDLOM `ExpressionElement`s during the transliteration phase
* Reference resolution between WDLOM elements occurs during the linking phase. 
* WDLOM `ExpressionElement`s are wrapped into `WdlomWomExpression`s during the graph building phase.

During the graph construction phase static expression elements are mixed together with evaluation functions
to produce the final WOM expressions, and any differences between language versions are baked in. 

#### Where evaluation functions are coded

The four evaluators types described in "Expressions in WOM" are coded in four package objects for each version. They live in:

* WDL `version 1.0`: [`wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/linking/expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/linking/expression)
* WDL `version development`: [`wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/linking/expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/biscayne/src/test/scala/wdl/transforms/biscayne/linking/expression)

These evaluators are really just long pattern matches from various WDLOM element types to the relevant evaluator for that element.

Because the evaluations are largely the same for now between WDL versions, you'll notice that the imported evaluators mostly come from files 
in [`wdl.transforms.base.linking.expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/linking/expression). However, in WDL development, there are a number of imported functions from biscayne-specific directories.
If future WDL versions diverge more starkly from the 1.0 base, it is likely that the imports into these pattern match evaluators will come from a
wider variety of origins.

#### How Evaluation Functions Build WOM Expressions

The top-level Graph construction functions for WDL 1.0 and the development version are:

* WDL `version 1.0`: [`wdl.draft3.transforms.wdlom2wom.workflowDefinitionElementToWomWorkflowDefinition`](https://github.com/broadinstitute/cromwell/blob/develop/wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/wdlom2wom/package.scala)
* WDL `version development`: [`wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/wdlom2wom/package`](https://github.com/broadinstitute/cromwell/blob/develop/wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/wdlom2wom/package.scala)

In each case:

* The package objects import the evaluation functions (from "Where evaluation functions are coded")
    * These imports fill in the implicit parameters required by `WorkflowDefinitionElementToWomWorkflowDefinition.convert`.
    * The imported conversion functions are used by the `convert` function every time it needs to make a WomExpression from WDLOM.
* The package objects create a language-specific conversion, `workflowDefinitionElementToWomWorkflowDefinition`.
    * This is used in the graph construction phase mentioned in "How WDL expressions become WomExpressions".
* The language-specific conversion is used by the appropriate LanguageFactory when it gets asked to construct WOM from WDL. 
# Workflow Execution: Major Actors

* **Word Count:** 245

## Major Actor Hierarchy

At the highest level, these are the main actors involved in workflow execution.

![high level overview diagram](WorkflowExecutionHighLevelOverview.png)

## Actors and their Purposes

### WorkflowManagerActor

The `WorkflowManagerActor` is responsible for:

* Polling the `WorkflowStore` at pre-configured intervals.
* Starting new workflows
* Tracking, supervising and aborting running workflows
* Parent actor for all `WorkflowActor`s

### WorkflowActor(s)

The `WorkflowActor` is responsible for:
 
* Co-ordinating the stages of a workflow lifecycle from parsing through to finalization.
* Parent actor of the `WorkflowExecutionActor` which runs the workflow's jobs.

### WorkflowExecutionActor(s)

The `WorkflowExecutionActor` is responsible for:

* Starting jobs and sub-workflows as soon as they are able to run.
    * Based on values in the (in-memory) ValueStore and ExecutionStore objects.
* Parent actor for all `EngineJobExecutionActor`s and `SubWorkflowExecutionActor`s.

### EngineJobExecutionActor(s)

Each `EngineJobExecutionActor` (EJEA) is responsible for:

* Running a single job.
    * A "job" is a command line instruction to run on a backend.
    * Multiple shards for a single call each get their own EJEA.
    * Multiple attempts to run the same job operate within the same EJEA
* Respects hog-limiting
* Checks the call cache and job store to avoid running the job if it doesn't have to.
* Triggers job initialization, execution and finalization at appropriate times.

### SubWorkflowExecutionActor(s)

Each `SubWorkflowExecutionActor` is responsible for:

* Running a single sub-workflow.
* Parent actor for a new `WorkflowExecutionActor` (see above) created to run the sub-workflow.

## Major Actor Hierarchy (in context)

The above diagram omitted a lot of details. This diagram attempts to show a little more of the
context:

![high level overview in context diagram](WorkflowExecutionHighLevelOverviewInContext.png)

## See Also 

* EngineJobExecutionActor (**TODO**)
* Backend Execution Actors (**TODO**)
Coming soon
# Workflow Execution: The Execution Store

## Purpose

The Execution Store is a data structure owned and maintained inside each 
`WorkflowExecutionActor` (see [Major Actors](majorActors.md)).

Its purpose is to hold the current _status_ of what the workflow is doing, 
what it has done, and what it has left to do. If the WOM graph holds a static
representation of the workflow which doesn't change as the workflow is run, the
Execution Store holds the dynamic status of each node in the graph at any given moment
in time.

The Execution Store is also used to answer questions like "which nodes in the
workflow graph are ready to run because their dependencies are satisfied?".

**Note:** The Execution Store does **not** hold the values which are generated by
running the various nodes in the WOM graph. That is the domain of the [Value Store](valueStore.md). 

## Data Structure

The Execution Store is a mapping from WOM nodes (and a shard index, if necessary)
to the execution status of those nodes in the actual workflow.

## Examples

Some worked through examples of how the Execution Store and Value Store change as workflows progress
are given on the [Execution and Value Store Examples](executionAndValueStoreExamples.md) page.# Workflow Execution: Database Tables: Workflow, Subworkflow, and Job Stores

## Database Tables: Workflow, Subworkflow, and Job Stores

Cromwell uses the workflow, subworkflow and job store tables to hold data related to submitted or running workflows.
Once a workflow reaches a terminal state all data for that workflow should be deleted from these tables.

### Workflow Store / `WORKFLOW_STORE_ENTRY`

`WORKFLOW_STORE_ENTRY` holds data received in a workflow submission (workflow sources, inputs, options etc.)
and workflow-scoped execution data (e.g. submission time, status, fields to support
running [Horizontal Cromwell](../horicromtal.md) etc).

### Job Store / `JOB_STORE_ENTRY`

`JOB_STORE_ENTRY` holds data for *completed* jobs within a workflow. Jobs that are still running or have not yet been
started will not have rows in this table. The main purpose of the job store table is to support resuming execution of
a workflow when Cromwell is restarted by recovering the outputs of completed jobs. This table is closely related to
`JOB_STORE_SIMPLETON_ENTRY` which holds the [simpleton](../general/simpletons.md) values comprising a job's outputs,
and loosely related to the [job key/value store (`JOB_KEY_VALUE_ENTRY`)](jobKeyValueStore.md) which holds other
job-scoped data important in recovering jobs on Cromwell restart.

### Subworkflow Store / `SUB_WORKFLOW_STORE_ENTRY`

`SUB_WORKFLOW_STORE_ENTRY` holds data for subworkflows that have begun execution. The rows in this table persist the fact
that particular subworkflows corresponding to a call FQN and index were started and assigned a workflow ID.
The completed jobs within these subworkflows will be recorded in the job store described above, linking to the
subworkflows in this table by the subworkflow's ID.
# Workflow Execution: The Value Store

## Purpose

The Value Store is a data structure owned and maintained inside each 
`WorkflowExecutionActor` (see [Major Actors](majorActors.md)).

Its purpose is to hold the set of _values_ produced by the workflow so far. If
the WOM graph holds a static representation of the workflow which doesn't change
as the workflow is run, the Value Store records the values assigned to every 
task output and value definition evaluated so far during workflow execution.

**Note:** The Value Store does **not** hold the execution status of the various 
nodes in the WOM graph. Nor does it determine when downstream nodes are ready
to run. That is the domain of the [Execution Store](executionStore.md). 

## Data Structure

The Value Store data structure is a mapping from output ports on WOM `GraphNode`s 
(and shard index if necessary) to the appropriate `WomValue`. 

## Examples

Some worked through examples of how the Execution Store and Value Store change as workflows progress
are given on the [Execution and Value Store Examples](executionAndValueStoreExamples.md) page.Coming soon
# 5000ft view of Cromwell Actor System 

Below diagram shows Cromwell's main actors and their interactions on a very high level.

![5000ft actor system overview diagram](5000_ft_Cromwell_Actors.png)
Coming soon.
### Intro to Workflow Object Model

"WOM" is the acronym for the Workflow Object Model, living in the `wom` package. WOM is a directed acyclic graph that captures workflow inputs, outputs, calls, and the dependencies between them.

Examine the workflow below, making note of how the outputs of early calls become inputs for later calls. The first example is the `numbers` output of the `mkFile` call serving as the `in_file` input to the `grep` call. 
```
version 1.0

##
# Checks a simple branch and join operation.
# We start with a task, branch into two parallel executions, and then rejoin to calculate the result.
##

workflow forkjoin {
  call mkFile

  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }

  call join { input: wcCount = wc.count, grepCount = grep.count }

  output {
    Int joined_proportion = join.proportion
  }
}

task mkFile {
  command <<<
    for i in `seq 1 1000`
    do
      echo $i
    done
  >>>
  output {
    File numbers = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  input {
    String pattern
    File in_file
  }
  command <<<
    [ -f "~{in_file}" ] && grep '~{pattern}' ~{in_file} | wc -l
  >>>
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  input {
    File in_file
  }
  command <<<
    [ -f "~{in_file}" ] && cat ~{in_file} | wc -l
  >>>
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task join {
  input {
    Int grepCount
    Int wcCount
  }
  command <<<
    expr ~{wcCount} / ~{grepCount}
  >>>
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}
```  

Now, compare the workflow source to its WOM graph (generated with the `womtool womgraph` command).

![Graph of forkjoin](forkjoin_graph.svg)

Input values are ovals, outputs are hexagons. You can see that a hexagon in one node (call) becomes an oval in the next. 

The `grep` call is special because one of its two inputs, `pattern`, is not specified by any previous call within the bounds of the graph. This causes Cromwell to generate a blue "external graph input node" (`wom.graph.ExternalGraphInputNode` if you're looking at the code). Its value must be specified by the user in the inputs file of the workflow under key `grep.pattern`, and Cromwell will pass it into the `grep` call.   

Finally, the `proportion` output of the `join` call is piped out of the bounds of the workflow graph by becoming the green `joined_proportion` graph output node (`wom.graph.ExpressionBasedGraphOutputNode` in the code). 

Graph nodes inherit from trait `wom.graph.GraphNode`. The graph is constructed in class `wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition`. The `convertGraphElements` function is especially interesting. It accepts a `Set` of `WorkflowGraphElement` objects, which represent individual pieces of the parsed workflow, and converts them to WOM nodes. Then it links the WOM nodes together with edges and emits the finished graph as a `wom.graph.Graph`.
### WDL Source to WOM Conversion

#### Parsing Flowchart

For the various versions of WDL supported in Cromwell, the conversion to WOM follows
these paths: 

![Parsing Flowchart](wdlmap.svg)

#### Process Description

You can think of WDL parsing in Cromwell in terms of the following major steps:

1. Lexing: Converting the raw WDL string into a one dimensional stream of "tokens".
2. Parsing: Converting the stream of tokens into an abstract syntax tree (AST).
3. Transliteration: Transforming the language-specific AST into a standard set of Scala objects
4. Import Resolution: Recursively processing any import statements into WOM bundles.
5. Linking: Discovering, resolving and recording all references within the AST and imports.
6. WOM Building: Creating a set of WOM objects
7. Input Validation: Link any provided inputs to inputs on the WOM objects.


#### Intermediate Data Formats

* **WDL Object Model (WDLOM)**:
    * A Scala case class representation of WDL grammar ASTs.
* **Linked inputs**:
    * The original WDL source's WDLOM representation
    * And WOM bundles imported
    * Links from any references to their definitions
        * Including custom type references, variable references, task calls, subworkflow calls 
* **WOM Bundle**:
    * A set of tasks, workflows and custom types, and the fully qualified names by which they can be referenced.
    * In Cromwell's WOM format (the WOM format is the ultimate destination for _all_ languages, including WDL and CWL)
* **Validated WOM Namespace**:
    * The conjunction of a WOM bundle with an input set.
    * The entry point workflow (or sometimes task, in CWL) is known.
### How To: Add Engine Functions 

Adding a new engine functions requires:

1. Allowing the engine function to be parsed into a WDLOM `ExpressionElement`.
2. Expression evaluators for processing the new element. 
### Converting WDL into WDLOM

To provide a concrete example, we will see how Cromwell parses the following line of WDL into WDLOM:

```wdl
Int foo = min(100, max(1,2))
``` 

#### How Hermes interprets the WDL

Hermes can be asked to show its parse tree for a valid WDL file by running:

```bash
$ hermes analyze grammar.hgr <WDL FILE>
```

This allows us to see how Hermes interprets our line of WDL:

```
(Declaration:
  type=<string:5:5 type "SW50">,
  name=<string:5:9 identifier "Zl9yb3VuZA==">,
  expression=(FunctionCall:
    name=<string:5:19 identifier "bWlu">,
    params=[
      <string:5:23 integer "MTAw">,
      (FunctionCall:
        name=<string:5:28 identifier "bWF4">,
        params=[
          <string:5:32 integer "MQ==">,
          <string:5:34 integer "Mg==">
        ]
      )
    ]
  )
)
```

In graphical form this is (with string hashes replaced by values, for convenience):

![Hermes AST Graph](wdlToWdlom_hermes.svg)

#### How WDLOM represents WDL

WDLOM tries to be a programmer-friendlier, WDL version agnostic data model to hold WDL syntax trees.

It would use the following data structure to represent this declaration:

```scala
InputDeclarationElement(
    typeElement = PrimitiveTypeElement(WomIntegerType),
    name = "foo",
    expression = Min(
        arg1 = PrimitiveLiteralExpressionElement(WomInteger(100)),
        arg2 = Max(
            arg1 = PrimitiveLiteralExpressionElement(WomInteger(1)),
            arg2 = PrimitiveLiteralExpressionElement(WomInteger(2))
        )
    )
)
```

Again, attempting to show this graphically:

![WDLOM AST Graph](wdlToWdlom_wdlom.svg)

#### Transliteration functions from Hermes ASTs to WDLOM

The various classes in the `wdl.transforms.base.ast2wdlom` package implement conversions from various types of AST element into various types of WDLOM.

Let's look at [`AstToInputDeclarationElement`](https://github.com/broadinstitute/cromwell/blob/master/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/ast2wdlom/AstToInputDeclarationElement.scala):

* The conversion from AST to WDLOM relies on two other conversions, `astNodeToTypeElement` and `astNodeToExpressionElement`, 
to validate the attributes on the AST it is given, and construct WDLOM from Hermes ASTs.
* The resulting function can itself be used as a building block to construct higher-level WDLOM types. 

#### How AST-to-WDLOM building blocks are chained together

To see how these building blocks are pieced together we can look at WDL 1.0's [`ast2wdlom`](https://github.com/broadinstitute/cromwell/blob/master/wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/ast2wdlom/ast2wdlom.scala)  package object.

Notice how - and where - the `astToInputDeclarationElement` value is declared:

* It follows the `implicit val`s `astNodeToTypeElement` and `astNodeToExpressionElement` and so can use those in its processing.
* It is, as an `implicit val`, used later on by `astNodeToInputsSectionElement`, which is then used by `astNodeToWorkflowBodyElement` and `astNodeToTaskSectionElement`, and so on. 

#### How different WDL versions convert to WDLOM differently

These package objects ultimately control how different WDL versions are parsed differently. These differences can be seen most easily by comparing package objects for different WDL versions:

```bash
diff \
  wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/ast2wdlom/ast2wdlom.scala \
  wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/ast2wdlom/ast2wdlom.scala   
``` 

In an extreme case, these could contain fundamentally different conversion logic from AST to WDLOM for different WDL versions. Luckily, between
WDL 1.0 and the WDL `development` version, there are only a few minor changes, amongst which are:

* Importing a new Hermes parser.
* Specifying that additional functions are available when constructing the `astNodeToExpressionElement` conversion.
* Specifying that a new `Directory` type is available when constructing the `astNodeToTypeElement` conversion.
# Workflow jobs: retry decision logic

* **Word count:** 132

## Key concepts

* Workflow may consist of multiple jobs.
* Job "retry" actually means creating a new job.
* Job may have attributes which are being stored in the `JOB_KEY_VALUE_ENTRY` table in the database, identified by a 
`ScopedKey`, which comprises workflow id, call fully qualified name, job index, job attempt number, and attribute name.
* There are 2 types of retries:
  * backend-specific retries (e.g., VM preemption in PAPI)
  * general retries 
* Backend-specific and general retries have separate retry counters, which are being stored in `JOB_KEY_VALUE_ENTRY` in
the end of the job execution attempt, and pre-fetched from the table in the beginning of the next attempt.

The retry logic is shown on the sequence diagram below with the example of PAPIv2 backend and VM preemption as an example 
of backend-specific retry reason.

![Job retry logic (example with PAPIv2 and VM preemption)](Workflow_job_retry_logic_(example_with_PAPIv2_and_VM_preemption).png)
Cromwell uses [Continuous Integration](https://en.wikipedia.org/wiki/Continuous_integration) (CI) testing, along with [Continuous Delivery](https://en.wikipedia.org/wiki/Continuous_delivery) (CD) to the Cromwell-as-a-Service (CaaS) `DEV` environment. [Continuous Deployment](https://en.wikipedia.org/wiki/Continuous_deployment) is not implemented.

## CI testing in Travis and Jenkins

Any suite of tests running under 2.5 hours, using 2 cpus, and 6gb of memory executes on [Travis](https://travis-ci.com/broadinstitute/cromwell/). Travis tests every pull request by a trusted contributor. Larger test suites run on [Jenkins instances](https://github.com/broadinstitute/dsp-jenkins#readme). Examples include the DSP Workbench CI testing ([swatomation](https://fc-jenkins.dsp-techops.broadinstitute.org/job/swatomation-pipeline/)), a nightly test of `develop` on (a snapshot of) the [$5 Genome WDL](https://fc-jenkins.dsp-techops.broadinstitute.org/job/cromwell-cron-aws/), and the [Cromwell-Perf tests](https://fc-jenkins.dsp-techops.broadinstitute.org/job/cromwell-perf-cron/) that call-cache thousands of jobs.

## CD to CaaS DEV

<a href="../CaaS_DEV_CD.svg" title="Cromwell automatically deploys to CaaS DEV using a series of hooks on GitHub, Travis, Jenkins, and Slack" target="_blank">![CaaS DEV CD](CaaS_DEV_CD.svg)</a>

One instance utilizing continuous delivery is the `develop` branch to CaaS `DEV`. While manual testing could occur on `DEV`, users primarily test on `PROD`.

## Manual Deployments

<a href="../Cromwell_Deployment_Strategies.svg" title="Each hosting platform deploys Cromwell differently" target="_blank">![Cromwell Deployment Strategies](Cromwell_Deployment_Strategies.svg)</a>
 
The Cromwell developers stage Terra and CaaS `PROD` deployments. All other deployments are performed by respective groups, who decide when and how to redeploy. Some upgrade Cromwell, while others deploy completely new instances, including a new database schema.

## Code Coverage

Only Travis Pull Requests generate [maximal code coverage reports](https://codecov.io/gh/broadinstitute/cromwell/pulls). All other CI either doesn't report coverage, or under-reports due to skipped tests.

## Vulnerability Scanning

In collaboration with the DSP Information Security team, scans include but are not limited to:

- Committed Git secrets
- Vulnerable Java dependencies
- Penetration testing
## Travis builds by user access

For infrastructures that require secured credentials, cloud backend tests only run for developers with write access to the broadinstitute/cromwell GitHub. Secure tests are skipped for all other users.

Other backends run tests for any user.

| Backend       | Read-only users | Write/Admin users |
|---------------|:---------------:|:-----------------:|
| AWS           |                 |        ✅         |
| BCS           |                 |        ✅         |
| Local         |       ✅        |        ✅         |
| PAPI V2alpha1 |                 |        ✅         |
| PAPI V2beta   |                 |        ✅         |
| SLURM         |       ✅        |        ✅         |
| TES           |       ✅        |        ✅         |

## Upgrade / Horicromtal / etc.

| CI Test Type                  | Cromwell Config                                                  | Centaur Config                                         |
|-------------------------------|------------------------------------------------------------------|--------------------------------------------------------|
| Engine Upgrade                | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| Horicromtal                   | `papi_[v2beta or v2alpha1]_horicromtal_application.conf`**       | `centaur_application_`<br>`horicromtal.conf`           |
| Horicromtal<br>Engine Upgrade | `papi_v2beta_application.conf`**                                 | `centaur_application_`<br>`horicromtal_no_assert.conf` |
| Papi Upgrade                  | `papi_v1_v2alpha1_upgrade_application.conf`**                    | `centaur_application.conf`*                            |
| Papi Upgrade<br>New Workflows | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| WDL Upgrade                   | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| (other)                       | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |

| CI Test Type                  | ScalaTest Spec              | Test Directory                      |
|-------------------------------|-----------------------------|-------------------------------------|
| Engine Upgrade                | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`            |
| Horicromtal                   | `CentaurTestSuite`          | `standardTestCases`***              |
| Horicromtal<br>Engine Upgrade | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`***         |
| PAPI Upgrade                  | `PapiUpgradeTestCaseSpec`   | `papiUpgradeTestCases`              |
| PAPI Upgrade<br>New Workflows | `CentaurTestSuite`          | `papiUpgradeNewWorkflowsTestCases`  |
| WDL Upgrade                   | `WdlUpgradeTestCaseSpec`    | `standardTestCases`****             |
| (other)                       | `CentaurTestSuite`          | `standardTestCases`                 |

<small>
\* Centaur Config always uses `centaur_application.conf` except when overridden with `papi_v2alpha1_centaur_application.conf`
or `papi_v2beta_centaur_application.conf`
  ([48 preview link](https://github.com/broadinstitute/cromwell/blob/a7d0601/src/ci/bin/test.inc.sh#L455-L457))  
\*\* Cromwell Config overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L213-L221))  
\*\*\* Test Directory overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L440-L449))  
\*\*\*\* Test Directory only tests tagged with `wdl_upgrade`
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/centaur/src/main/resources/standardTestCases/write_lines.test#L3))  
</small>

- Engine Upgrade: Retrieves the [Cromwell Version](https://github.com/broadinstitute/cromwell/blob/47/project/Version.scala#L8) then retrieves the previous jar/docker-image from DockerHub. Centaur starts with the prior version, then restarts with the compiled source code.
- Horicromtal: Runs a [docker-compose](https://github.com/broadinstitute/cromwell/blob/47/src/ci/docker-compose/docker-compose-horicromtal.yml) with:
    1. db-mstr: started first
    2. sum-back: runs summarizer
    3. front-back: exposes HTTP
- Horicromtal Engine Upgrade: Combination of Horicromtal and Engine Upgrade
- PAPI Upgrade: Tests run with an older version of Papi and upon restart use a newer version of Papi
- PAPI Upgrade New Workflows: Test definition [does not run any tests](https://travis-ci.org/broadinstitute/cromwell/jobs/475378412)
- WDL Upgrade: Upgrades WDL from draft-2 to 1.0 before testing
- (other): Runs `*.test` files listing the configured backend names

## RDBMS

| Backend | MySQL  | PostgreSQL  | MariaDB  |
|---------|:------:|:-----------:|:--------:|
| AWS     |   ✅   |             |          |
| BCS     |   ✅   |             |          |
| Local   |   ✅   |      ✅     |          |
| PAPI V2 |   ✅   |             |    ⭕    |
| SLURM   |   ✅   |             |          |
| TES     |   ✅   |             |          |

<small>
⭕ Tests Horicromtal Engine Upgrade versus standard Centaur suite
</small>

All backends run against MySQL. The Local backend also test PostgreSQL, allowing contributors ensure WDLs work with PostgreSQL. MariaDB is tested on a specialized upgrade, where the MySQL connector client is used first, and the MariaDB client is used after restart.
# The 'localization_optional' Optimization

Available in Cromwell version 33 and higher.

## Scope

The 'localization_optional' optimization can be applied to a task's individual input declarations containing files, specifically `File` and `File?` values and any complex types containing them. 
It allows you to save time and money by identifying files which do not need to be localized for the task to succeed.

## Condition

The optimization signals to Cromwell that a task has been written in such a way that:
 
 * The task **will work** if Cromwell does localize the specified file inputs
   * For example if a file is localized for a local dockerized execution environment.

**And**:

 * The task will **also** work if Cromwell **does not** localize the same file input
   * For example the file remains in a cloud object store and the command is constructed using its URL rather than a local path.

## Effect on File Localization

If the [backend](#backend-support) has been set up to respect `localization_optional`, Cromwell will 
choose not to localize the appropriate file input.

### Effect on Call Caching:

None! 

Files marked for optional localization are still treated in exactly the same way as other `File` inputs for call caching.

## Language Support

### WDL 1.0 (or later)

In a WDL 1.0 `task`, this optimization is specified by adding a `localization_optional` field to 
an input's entry in the task's `parameter_meta` section. Here's an example:

```wdl
task nio_task {
  input {
    File foo_file
    File bar_file
  }
  
  parameter_meta {
    foo_file: {
      description: "a foo file",
      localization_optional: true
    }
    bar_file: {
      description: "a bar file"
    }
  }
  
  command <<<
    # This tool must work for **BOTH** local file paths **AND** object store URL values:
    java -jar my_tool_1.jar ~{foo_file}
    
    # Because the optimization is not applied to 'bar_file' in parameter_meta, this file **WILL** be localized:
    java -jar my_tool_2.jar ~{bar_file}
  >>>
}
```

## Backend Support

This optimization is currently only applied to localization in the Pipelines API (GCE) backends.
# The 'volatile' Optimization

Available in Cromwell version 49 and higher.

### Effect on Call Caching:

The 'volatile' optimization is applied to tasks in their `meta` section.
Call caching will be disabled for any call to that task during the execution of the workflow. 

This is particularly useful if:

* One task can produce stochastic results but you still want to use call caching in the rest of the workflow.
* You want to guarantee that a task is never call cached for any other reason.

## Language Support

### WDL

In a WDL `task`, this optimization is specified by adding a `volatile` field to 
the task's `meta` section. Here's an example:

```wdl
version 1.0
 
task make_random_int {
  
  meta {
    volatile: true
  }
  
  command <<<
    echo $RANDOM
  >>>

  output {
    Int random = read_string(stdout())
  }
}
```

## Backend Support

The volatile keyword applies equally to all backends.
# Optimization

These optimizations are *Cromwell-specific* functionality which can be triggered from within your workflow descriptions.

You can think of them as ways of telling Cromwell that the task or workflow has a certain property which allows Cromwell to do something clever. 
 

## Portability Warning

Note that these optimizations are *outside* of the language specifications and so not all workflow engines will respect them.
In order to maintain portability of workflows, write defensively with respect to these optimizations: 

  * Remember that a Cromwell instance might have your optimization turned off.
  * Remember that your workflow might need to run on a version of Cromwell which predates the optimization.
  * Remember that to share your WDL most widely, it will need to be able to run on engines other than Cromwell - and those engines won't necessarily respect these optimizations.
# The 'Checkpoint File' Optimization

## Overview

Available in Cromwell 55 and higher.

This optimization provides a way to mitigate the problem of long-running tasks getting preempted partway through execution. It allows the user to save intermediates of a task at regularly scheduled checkpoints. After an interruption, the task can be restarted from the last saved checkpoint, thereby avoiding having to re-run the entire computation again.

### Description

Specifying a `checkpointFile` value in a task's `runtime` section designates a checkpoint file which will periodically be
copied to cloud storage every 10 minutes. This checkpoint file will then be restored automatically on subsequent attempts if the job is interrupted. Once the final output has been successfully generated, the checkpoint file will be deleted.

To use this feature effectively, the WDL task must be written intentionally to use the checkpoint file. See example below. 

### Effects on cloud charges

Charges will accrue from storing the checkpoint file during the running of the task, and additional charges may apply to the transfer between the VM and the cloud storage bucket depending on their locations. These costs should be minor, especially balanced against the performance and cost benefits of being able to restore from the checkpoint when a worker VM gets preempted.

Since the checkpoint file is deleted after successful completion of the task, no further charges will accrue after completion. However, if the task is aborted or otherwise stopped externally, ie through interruption of Cromwell's operation, the checkpoint file will NOT be deleted and storage charges will continue to accrue indefinitely, or until the file is deleted manually. 

### Effect on Call Caching

The presence or absence of the `checkpointFile` attribute is not considered when determining whether to call cache.  

### Example

The following WDL demonstrates the use of the `checkpointFile` optimization. It has a command that is checkpoint-aware:

* It starts by attempting to restore state from the `my_checkpoint` file (or starts at `1` if the checkpoint is empty)
* Then it counts up to 100, printing out the current counter value and a date timestamp at each value.

To make the checkpointing work, the `runtime` section specifies `checkpointFile: "my_checkpoint"`.

```wdl
version 1.0

workflow count_wf {
  # Count to 2100 at 1/second => 35 minutes to complete, but
  # luckily the state can be checkpointed every 10 minutes in
  # case of preemption: 
  call count { input: count_to = 2100 }
}

task count {
  input {
    Int count_to
  }

  command <<<
    # Note: Cromwell will stage the checkpoint file on recovery attempts.
    # This task checks the 'my_checkpoint' file for a counter value, or else
    # initializes the counter at '1':
    FROM_CKPT=$(cat my_checkpoint | tail -n1 | awk '{ print $1 }')
    FROM_CKPT=${FROM_CKPT:-1}

    echo '--' >> my_checkpoint
    for i in $(seq $FROM_CKPT ~{count_to})
    do
      echo $i $(date) >> my_checkpoint
      sleep 1
    done
  >>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
    # Note: This checkpointFile attribute is what signals to Cromwell to save
    # the designated checkpoint file:
    checkpointFile: "my_checkpoint"
  }

  output {
    # Note: This task also uses the checkpoint as its output. This is not
    # required for checkpointing to work:
    Array[String] out = read_lines("my_checkpoint")
  }
}
```

## Backend Support

Cromwell supports the `checkpointFile` attribute on the following backends:

* The Google PAPIv2 (alpha1) backend
* The Google Life Sciences (beta) backend
# Filesystems

Most workflows represent their inputs and outputs in the form of files. Those files are stored in filesystems. There exists many filesystems. This section describes which filesystems Cromwell supports.

## Overview

Filesystems are configurable. The `reference.conf`, which is the configuration inherited by any Cromwell instance, contains the following:

```hocon
# Filesystems available in this Crowmell instance
# They can be enabled individually in the engine.filesystems stanza and in the config.filesystems stanza of backends
# There is a default built-in local filesytem that can also be referenced as "local" as well.
filesystems {
  drs {
      class = "cromwell.filesystems.drs.DrsPathBuilderFactory"
      # Use to share a unique global object across all instances of the factory
      global {
        # Class to instantiate and propagate to all factories. Takes a single typesafe config argument
        class = "cromwell.filesystems.drs.DrsFileSystemConfig"
        config {
          martha {
            url = "https://martha-url-here"
            # The number of times to retry failures connecting or HTTP 429 or HTTP 5XX responses, default 3.
            num-retries = 3
            # How long to wait between retrying HTTP 429 or HTTP 5XX responses, default 10 seconds.
            wait-initial = 10 seconds
            # The maximum amount of time to wait between retrying HTTP 429 or HTTP 5XX responses, default 30 seconds.
            wait-maximum = 30 seconds
            # The amount to multiply the amount of time to wait between retrying HTTP or 429 or HTTP 5XX responses.
            # Default 2.0, and will never multiply the wait time more than wait-maximum.
            wait-mulitiplier = 2.0
            # The randomization factor to use for creating a range around the wait interval.
            # A randomization factor of 0.5 results in a random period ranging between 50% below and 50% above the wait
            # interval. Default 0.1.
            wait-randomization-factor = 0.1
          }
        }
      }
   }
  gcs {
    class = "cromwell.filesystems.gcs.GcsPathBuilderFactory"
  }
  oss {
    class = "cromwell.filesystems.oss.OssPathBuilderFactory"
  }
  s3 {
    class = "cromwell.filesystems.s3.S3PathBuilderFactory"
  }
  http {
    class = "cromwell.filesystems.http.HttpPathBuilderFactory"
  }
}
```

It defines the filesystems that can be accessed by Cromwell.
Those filesystems can be referenced by their name (`drs`, `gcs`, `oss`, `s3`, `http` and `local`) in other parts of the configuration.

**Note:**
- **OSS and S3 filesystems are experimental.** 
- **DRS filesystem has initial support only. Also, currently it works only with [GCS filesystem](../GoogleCloudStorage) in [PapiV2 backend](http://cromwell.readthedocs.io/en/develop/backends/Google).**


Also note that the local filesystem (the one on which Cromwell runs on) is implicitly accessible but can be disabled. 
To do so, add the following to any `filesystems` stanza in which the local filesystem should be disabled: `local.enabled: false`.

### Engine Filesystems

Cromwell is conceptually divided in an engine part and a backend part. One Cromwell instance corresponds to an "engine" but can have multiple backends configured.
The `engine.filesystems` section configures filesystems that Cromwell can use when it needs to interact with files outside of the context of a backend.

For instance, consider the following WDL:

```wdl
version 1.0

workflow my_workflow {
    String s = read_string("/Users/me/my_file.txt")
    output {
        String out = s
    }
}
```

This workflow is valid WDL and does not involve any backend, or even a task. However it does involve interacting with a filesystem to retrieve the content of `my_file.txt`
With a default configuration Cromwell will be able to run this workflow because the local filesystem is enabled by default.
If the file is located on a different filesystem (a cloud filesystem for instance), we would need to modify the configuration to tell Cromwell how to interact with this filesystem:

```hocon
engine {
  filesystems {
    gcs {
      auth = "application-default"
    }
  }
}
```

(See the [Google section](../backends/Google.md) for information about the `auth` field.)

We can now run this workflow

```wdl
version 1.0

workflow my_workflow {
    String s = read_string("gs://mybucket/my_file.txt")
    output {
        String out = s
    }
}
```

#### Default "engine" Filesystems

If you don't change anything in your own configuration file, the following default is inherited from `reference.conf`:
```
engine {
  filesystems {
    local {
      enabled: true
    }
    http {
      enabled: true
    }
  }
}
```

**Note**: since our configuration files are HOCON, to disable filesystems you *must* add `enabled: false` into your 
overriding configuration file. It is **not** sufficient to simply omit a filesystem from your stanza. 

For example: adding this to your configuration file will remove the `http` filesystem and leave `local` for use in the 
engine:
```
engine {
  filesystems {
    http {
      enabled: false
    }
  }
}
```

Whereas this example will leave `http` unchanged and merely re-assert the default enabling of `local`. In other
words, **this will do nothing**:
```
engine {
  filesystems {
    local {
      enabled: true
    }
  }
}
```

#### Engine filesystems and CWL

Note that CWL *always* needs to access file attributes from within the engine, so if you are using CWL, please make sure
that **every filesystem you might use** is added to the `engine.filesystems` stanza.

### Backend Filesystems

Similarly to the engine, you can also configure backend filesystems individually. Some backends might require the use of a specific filesystem.
For example, the [Pipelines API](../tutorials/PipelinesApi101.md) backend requires Google Cloud Storage.
Let's take another example:

```wdl
version 1.0

task my_pipelines_task {
    input {
        File input_file
    }
    String content = read_string(input_file)
    
    command {
        echo ~{content}
    }
    
    runtime {
        docker: "ubuntu"
    }
}
workflow my_workflow {
    call my_pipelines_task { input: input_file = "gs://mybucket/my_file.txt" }
}
```

Suppose this workflow is submitted to a Cromwell running a Pipelines API backend. This time the `read_string` function is in the context of a task run by the backend.
The filesystem configuration used will be the one in the `config` section of the Pipelines API backend.

### Supported Filesystems

-  Shared File System (SFS)

-  Google Cloud Storage (GCS) - [Cromwell Doc](GoogleCloudStorage.md) / [Google Doc](https://cloud.google.com/storage/)

-  Simple Storage Service (S3) - [Amazon Doc](https://aws.amazon.com/documentation/s3/)

-  Object Storage Service (OSS) - [Alibaba Cloud Doc](https://www.alibabacloud.com/product/oss)

-  HTTP - support for `http` or `https` URLs for [workflow inputs only](http://cromwell.readthedocs.io/en/develop/filesystems/HTTP)

- File Transfer Protocol (FTP) - [Cromwell Doc](FileTransferProtocol.md)
# Google Cloud Storage (GCS)

## Overview 

Cromwell supports workflows referencing objects stored in [Google Cloud Storage](https://cloud.google.com/storage/).
The Cromwell configuration for GCS is as follow:

```hocon
filesystems {
  gcs {
    # A reference to a potentially different auth for manipulating files via engine functions.
    auth = "application-default"

    # Google project which will be billed for requests on buckets with requester pays enabled
    project = "google-billing-project"

    caching {
      # When a cache hit is found, the following duplication strategy will be followed to use the cached outputs
      # Possible values: "copy", "reference". Defaults to "copy"
      # "copy": Copy the output files
      # "reference": DO NOT copy the output files but point to the original output files instead.
      #              Will still make sure than all the original output files exist and are accessible before
      #              going forward with the cache hit.
      duplication-strategy = "copy"
    }
  }
}
```

- The `auth` field refers to the authentication schema that should be used to authenticate requests. See [here](../backends/Google.md) for more info.
- The `project` field has to do with the Requester Pays feature (see below).
- The `caching.duplication-strategy` field determines how Cromwell should behave w.r.t output files when call is being cached. The default strategy `copy` is to copy the file to its new call location. As mentioned, `reference` will not copy the file and simply point the results to the existing location.
See the [Call Caching documentation](../cromwell_features/CallCaching.md) for more information.

## Requester Pays

GCS has a feature called Requester Pays (RP). This section describes how Cromwell supports it and the consequences on cost. Please first read the [official documentation](https://cloud.google.com/storage/docs/requester-pays) if you're not already familiar with it.

The billing project Cromwell uses to access a bucket with requester pays is determined as follows:

- If a `google_project` was set in the [workflow options](../wf_options/Google.md) when the workflow was submitted, this value is used
- Otherwise, the value of the `project` field in the `gcs` filesystem configuration is used
- Otherwise, if the machine Cromwell runs on is authenticated using gcloud and a default project is set, this value will be used

**Important Note #1**: In order for a project to be billable to access a bucket with requester pays, the credentials used need to have the `serviceusage.services.use` permission on this project. 

**Important Note #2**: Pipelines API version 1 does **not** support buckets with requester pays, so while Cromwell itself might be able to access bucket with RP, jobs running on Pipelines API V1 with file inputs and / or outputs will **not** work.
For full requester pays support, use the [Pipelines API v2 Cromwell backend](https://github.com/broadinstitute/cromwell/blob/develop/CHANGELOG.md#pipelines-api-v2). 

**Important Note #3**: Access to requester pays buckets from Cromwell is seamless, this also means that Cromwell will not report in the logs or metadata when it access a bucket with requester pays. It is the user's responsibility to be aware of the extra cost of running workflows access requester pays buckets.
# HTTP Inputs

## Overview

For shared filesystem and Google Pipelines API (PAPI) version 2 backends Cromwell can support workflow inputs specified by `http` and `https` URLs.
Please note this is not true "filesystem" support for HTTP URLs;
if inputs to a workflow are specified by HTTP URLs the outputs of steps will nevertheless appear at local or GCS paths and not HTTP
URLs.

### Configuration

Cromwell's default configuration defines an instance of the HTTP filesystem named `http`. There is no additional configuration
required for the HTTP filesystem itself so adding HTTP filesystem support to a backend is a simple as
adding a reference to this filesystem within the backend's `filesystems` stanza. e.g. Cromwell's default `Local` shared filesystem
backend is configured like this (a PAPI version 2 backend would be configured in a similar way):

```
backend {
  default = "Local"
  providers {
    Local {
      ...
      config {
        filesystems {
          local {
            ...
          }
          http { }
        }
      }
      ...
    }
    ...
  }
}
```

If there is a need to turn off this `http` filesystem in the default `Local` backend the following Java property
allows for this: `-Dbackend.providers.Local.config.filesystems.http.enabled=false`.

### Caveats

Using HTTP inputs in Cromwell can produce some unexpected behavior:
- Files specified by HTTP URIs will be renamed locally, so programs that rely on file extensions or other filenaming conventions may not function properly.
- Files located in the same remote HTTP-defined directory will not be colocated locally. This can cause problems if a program is expecting an index file (e.g. `.fai`) to appear in the same directory as the associated data file (e.g. `.fa`) without specifying the index location.
# Data Repository Service (DRS)

The Cromwell configuration for DRS is as follows:

**Filesystem Configuration**

```hocon
drs {
    # A reference to a potentially different auth required to contact Martha service.
    auth = "application-default"
}
```

The `auth` field refers to the authentication schema that should be used to authenticate requests to Martha service.

The `drs` section needs to be added to
- `engine.filesystems` block
- [PapiV2](http://cromwell.readthedocs.io/en/develop/backends/Google) backend's `filesystems` block


**Localization Configuration**

DRS localization must be configured with the docker image to use.

```hocon
drs {
    localization {
        docker-image = "broadinstitute/drs-localizer:latest"
    }
}
```


**Example**

A sample configuration for DRS filesystem might look like:

```hocon
engine {
  filesystems {
    # ... other filesystems here, probably gcs, and then ...
    drs {
      auth = "application-default"
    }
  }
}

backend {
    # ... other global backend config here, probably just setting the default ...
    providers {
        # ... other providers here ...
        Papi {
            actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
            config {
                # ... other config here ...
                filesystems {
                    # ... other filesystems here, probably gcs, and then ...
                    drs {
                        auth = "application-default"
                    }
                }
            }
        }
    }
}

drs {
    localization {
        docker-image = "broadinstitute/drs-localizer:latest"
    }
}
```
# File Transfer Protocol (FTP)

Cromwell supports communication with basic FTP servers (not FTPS, and not SFTP either).
Please read the [#Filesystems](Filesystems.md) section before to get an explanation of how filesystems can be configured in general.
For a full example of FTP configuration see [here](#Example)

**Note**: Be aware that FTP is in many cases very inefficient, not secure and generally doesn't scale well. If possible, consider using an alternate file system.

## Overview

Cromwell handles FTP connections as follows:

- Cromwell maintains one FTP FileSystem per FTP server per user
- An FTP FileSystem maintains a pool of connections with a fixed size (see `max-connection-per-server-per-user` configuration below) to that server, for that user
- When an FTP FileSystem hasn't been used in a certain amount of time (see `cache-ttl` configuration below), the associated connections are closed and it is destroyed
- In a given pool, when a connection has been idle for a certain amount of time (see `idle-connection-timeout` configuration below), it is closed

```
 Cromwell
+---------------------------------------------------------------------------------------------------+
|                                                                                                   |
|  FileSystem 1                     FileSystem 2                     FileSystem 3                   |
| +-----------------------------+  +-----------------------------+  +-----------------------------+ |
| |                             |  |                             |  |                             | |
| | ftp://user1@server1.com     |  | ftp://user2@server1.com     |  | ftp://user1@server2.com     | |
| |                             |  |                             |  |                             | |
| |  Connection Pool            |  |  Connection Pool            |  |  Connection Pool            | |
| | +-------------------------+ |  | +-------------------------+ |  | +-------------------------+ | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 1         | |  | |  - Connection 1         | |  | |  - Connection 1         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 2         | |  | |  - Connection 2         | |  | |  - Connection 2         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 3         | |  | |  - Connection 3         | |  | |  - Connection 3         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | +-------------------------+ |  | +-------------------------+ |  | +-------------------------+ | |
| |                             |  |                             |  |                             | |
| +-----------------------------+  +-----------------------------+  +-----------------------------+ |
+---------------------------------------------------------------------------------------------------+
```


## Global Configuration

The FTP filesystem supports a global configuration. The goal is to configure values that will apply Cromwell-wide, as opposed to other configuration that could be backend (or engine) specific (see below).

Default configuration:

```hocon
filesystems.ftp.global.config = {
    # This value should be sufficiently high to cover for the duration of the longest expected I/O operation.
    # It is a time to live value after which a filesystem (unique per user per server) will be closed and evicted from the cache if unused for the specified duration.
    cache-ttl = 1 day
    
    # How long to wait trying to obtain a connection from the pool before giving up. Don't specify for no timeout
    # obtain-connection-timeout = 1 hour
    
    # Maximum number of connections that will be established per user per ftp server. This is across the entire Cromwell instance.
    # Setting this number allows to workaround FTP server restrictions for the number of connections for a single user per IP address.
    # Has to be >= 2 to allow copying (Copying and FTP file requires downloading and uploading it)
    max-connection-per-server-per-user = 30
    
    # Time after which a connection will be closed if idle. This is to try to free connections from a filesystem when it's not heavily used.
    idle-connection-timeout = 1 hour
    
    # FTP connection port to use
    connection-port: 21
    
    # FTP connection mode
    connection-mode = "passive"
}

```

Note that this configuration being global, it applies the same way regardless of the FTP server. There is currently no good way to make this configuration dependent on the server.

## Instance Configuration

Configuration that can be applied to a backend or engine filesystem stanza:

```hocon
ftp {
  # optional
  auth {
    username = "username"
    password = "password"
    # Optional
    account = "account"
  }
}
```

A default authentication can be provided in the configuration. It can be overridden using the following workflow options: `ftp-username`, `ftp-password`, `ftp-account`
If authentication is omitted the connection is established as an anonymous user.

## Example

Example of configuration mixing the above two sections:

```hocon
filesystems.ftp.global.config {
    # We only override what changes from the default
    max-connection-per-server-per-user = 20
}

# Configure the default auth for the engine
engine.filesystems.ftp {
  auth {
    username = "me"
    password = "my_password"
  }
}

# Configure the default auth for the local backend
backend.providers.Local.config.filesystems.ftp {
  auth {
    username = "me"
    password = "my_password"
  }
}
```
**Configuration the Local Backend**

<!--
### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)


### Goals

At the end of this tutorial you'll have seen how to add a configuration section for the Local backend to your configuration file, and seen what changing some of the various options does.

### Let's get started!


### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
-->

_Drop us a line in the [Forum](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team) if you have a question._

\*\*\* **UNDER CONSTRUCTION** \*\*\*  
[![Pennywell pig in red wellies - Richard Austin Images](http://www.richardaustinimages.com/wp-content/uploads/2015/04/fluffyAustin_Pigets_Wellies-500x395.jpg)](http://www.richardaustinimages.com/product/pennywell-pigs-under-umbrella-2/)
# Five minute Introduction to Cromwell

### Prerequisites:

* A Unix-based operating system (yes, that includes Mac!)
* A Java 8 runtime environment 
	* You can see what you have by running `$ java -version` on a terminal. You're looking for a version that's at least `1.8` or higher.
	* If not, you can download Java [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
* A sense of adventure!

### Goals

At the end of this five minute introduction you will have:

- Downloaded Cromwell!
- Written your first workflow
- Run it through Cromwell

### Step 1: Downloading Cromwell

We host our Cromwell releases on GitHub! You can find the latest version on our [Releases](https://github.com/broadinstitute/cromwell/releases/latest) page.

* Look for the latest version at the top of the page, and find the link to download the jar. It'll have a name like `cromwell-XY.jar` where `XY` is the version. Download the jar file.
* WARNING! If you're on a Mac, the security settings might try to stop you from running Cromwell! Don't worry, if this happens just go to `System Preferences > Security & Privacy > General` and find the `cromwell` jar listed on the page. Click `Open anyway`. The `cromwell-XY.jar` will now automatically download to your `Downloads` directory.
* Put your downloaded Cromwell somewhere you can find it later, like in a Cromwell directory in your home directory.

For example, in a terminal:
```sh
cd ~
mkdir cromwell
cp ~/Downloads/cromwell-XY.jar cromwell/
cd cromwell/
```
_(if you're not using a Mac, the final command might be different for you)_


### Step 2: Writing your first workflow description

This bit is easy, you're just going to copy and paste something from the internet.

Open your favorite editor. Paste in the following content and save it as `myWorkflow.wdl` in your new `cromwell` directory:

```wdl
workflow myWorkflow {
	call myTask
}

task myTask {
	command {
		echo "hello world"
	}
	output {
		String out = read_string(stdout())
	}
}
```

Don't worry, **you don't need to understand too much about the workflow contents to continue for now**. In brief, it tells Cromwell to run a task to run `echo "hello world"`, and then return the output as a String. If you'd like to learn more about how to author WDL, you can find all the WDL resources you could ever want [here](https://github.com/openwdl/wdl).

### Step 3: Running the workflow

Ok, we have Cromwell, we have a workflow, let's put it all together! 

Make sure you're in the cromwell directory with the `.jar` file and the `.wdl` file. Now type in:
```sh
java -jar cromwell-XY.jar run myWorkflow.wdl
```

Cromwell will print out a fair old chunk of logging information, which can be configured (once you've completed this tutorial and [Configuration Files](ConfigurationFiles), you might want to investigate the [Logging](../Logging) page)

Ultimately, the workflow should succeed and you'll end up with the following output printed out when Cromwell finishes:
```json
{
	"myWorkflow.myTask.out": "hello world"
}
```

Ok, you can stop your timer! You just installed and ran your first workflow in Cromwell, congratulations!

### Next Steps

Pat yourself on the back for completing this tutorial, bravo! Then continue on to one of the follow pages:

* [Server Mode](ServerMode)
* [Configuration Files](ConfigurationFiles)
## Getting started on Alibaba Cloud with the Batch Compute Service

### Prerequisites

This tutorial page relies on completing the previous tutorials:

- [Configuration Files](ConfigurationFiles.md)

### Goals

In this tutorial you'll learn to run the first workflow against the Batch Compute service on Alibaba Cloud.

### Let's get started!

####

#### Configuring Alibaba Cloud

- Go to <a href="https://www.aliyun.com/" target="_blank">Alibaba Cloud</a> and activate <a href="https://www.aliyun.com/product/oss">Alibaba Cloud OSS</a> and <a href="https://www.aliyun.com/product/batchcompute">Alibaba Cloud BatchCompute</a> services. 
- Follow <a href="https://help.aliyun.com/document_detail/63724.html" target="_blank">AccessKey Guide</a> to retrieve an access-id and access-key pair. We will refer to this pair as `<test-access-id>` and `<test-access-key>`, respectively.
- Log on to the <a href="https://oss.console.aliyun.com/" target="_blank">OSS console</a> and choose a region to create a new bucket. We will use `<test-region>` and `<test-bucket>` to refer the chosen region and bucket.
- Find the corresponding OSS endpoint in <a href="https://help.aliyun.com/document_detail/31837.html" target="_blank">OSS region and endpoint</a>. We will refer to it as `<test-oss-endpoint>`.

#### Preparing workflow source files

Copy over the sample `echo.wdl` and `echo.inputs` files to the same directory as the Cromwell jar. 
This workflow takes a string value as an output file name and writes "Hello World!" to the file. 

***echo.wdl***

```
task echo {
  String out

  command {
    echo Hello World! > ${out}
  }

  output {
    File outFile = "${out}"
    Array[String] content = read_lines(outFile)
  }
}

workflow wf_echo {
  call echo
  output {
    echo.outFile
    echo.content
  }
}
```

***echo.inputs***

```
{
  "wf_echo.echo.out": "output"
}
```

#### Configuration file for Alibaba Cloud

Copy over the sample `bcs.conf` file to the same directory that contains your sample WDL, inputs and the Cromwell jar. Replace `<test-bucket>`, `<test-region>`, `<test-access-id>`, `<test-access-key>`, `<test-oss-endpoint>` in the configuration file with actual values.  

***bcs.conf***

```
include required(classpath("application"))

backend {
  default = "BCS"
  
  providers {
    BCS {
      actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
      config {
        root = "oss://<test-bucket>/cromwell-dir"
        region = "<test-region>"
        access-id = "<test-access-id>"
        access-key = "<test-access-key>"
        
        filesystems {
          oss {
            auth {
              endpoint = "<test-oss-endpoint>"
              access-id = "<test-access-id>"
              access-key = "<test-access-key>"
            }
          }
        }
        
        default-runtime-attributes {
          failOnStderr: false
          continueOnReturnCode: 0
          cluster: "OnDemand ecs.sn1ne.large img-ubuntu"
          vpc: "192.168.0.0/16"
        } 
      }
    }
  }
}
```

#### Run workflow

`java -Dconfig.file=bcs.conf -jar cromwell.jar run echo.wdl --inputs echo.inputs`

#### Outputs

The end of your workflow logs should report the workflow outputs. 

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_echo.echo.outFile": "oss://<test-bucket>/cromwell-dir/wf_echo/38b088b2-5131-4ea0-a161-4cf2ca8d15ac/call-echo/output",
    "wf_echo.echo.content": ["Hello World!"]
  },
  "id": "38b088b2-5131-4ea0-a161-4cf2ca8d15ac"
}
```
# Tutorial Map

This map shows the tutorials in this section and the orders in which they can be read:

![Map](UserGuideMap.png)
**Viewing Metadata**

<!--
### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Server Mode](ServerMode)


### Goals

At the end of this tutorial you'll have seen how to query the metadata for a Cromwell workflow and the information that it provides.

### Let's get started

### Next Steps

After completing this tutorial you might find the following page interesting:

* [Timing Diagrams](TimingDiagrams)
* [Configuration Files](ConfigurationFiles)
-->

_Drop us a line in the [Forum](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team) if you have a question._

\*\*\* **UNDER CONSTRUCTION** \*\*\*  
[![Pennywell pig in red wellies - Richard Austin Images](http://www.richardaustinimages.com/wp-content/uploads/2015/04/fluffyAustin_Pigets_Wellies-500x395.jpg)](http://www.richardaustinimages.com/product/pennywell-pigs-under-umbrella-2/)
## Getting started on Google Cloud with the Genomics Pipelines API

## Pipelines API v2

### Basic Information

Initial support for Google [Pipelines API version 2](https://cloud.google.com/genomics/reference/rest/) was added in Cromwell 32.
Expect feature parity with v1 except:

* PAPI v2 private Docker support is equivalent to that in PAPI v1 but the configuration differs, please see [Docker configuration](http://cromwell.readthedocs.io/en/develop/filesystems/Google#Docker) for more details.
* The "refresh token" authentication mode is **NOT** supported on PAPI V2.

In addition, the following changes are to be expected:

* Error messages for failed jobs might differ from V1
* The Pipelines API log file content might differ from V1

### Setting up PAPIv2

For now the easiest way to try PAPIv2 is to migrate an existing set up from PAPIv1 (see below). After that, copy the PAPIv2 sample configuration in [cromwell.examples.conf](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/PAPIv2.conf) in place of the PAPIv1 backend.

#### Permissions:

Google recommends using a service account to authenticate to GCP.  

You may create a service account using the `gcloud` command, consider running the following script and replace MY-GOOGLE-PROJECT:

```
#!/bin/bash
export LC_ALL=C 
RANDOM_BUCKET_NAME=$(head /dev/urandom | tr -dc a-z | head -c 32 ; echo '')

#Create a new service account called "my-service-account", and from the output of the command, take the email address that was generated
EMAIL=$(gcloud beta iam service-accounts create my-service-account --description "to run cromwell"  --display-name "cromwell service account" --format json | jq '.email' | sed -e 's/\"//g')

# add all the roles to the service account
for i in storage.objectCreator storage.objectViewer lifesciences.workflowsRunner lifesciences.admin iam.serviceAccountUser storage.objects.create
do
    gcloud projects add-iam-policy-binding MY-GOOGLE-PROJECT --member serviceAccount:"$EMAIL" --role roles/$i
done

# create a bucket to keep the execution directory
gsutil mb gs://"$RANDOM_BUCKET_NAME"

# give the service account write access to the new bucket
gsutil acl ch -u "$EMAIL":W gs://"$RANDOM_BUCKET_NAME"

# create a file that represents your service account.  KEEP THIS A SECRET.
gcloud iam service-accounts keys create sa.json --iam-account "$EMAIL"
```

## Pipelines API v1

### Deprecation warning

Please note that Google intends to deprecate PAPIv1 in the near future (circa mid 2019 or perhaps earlier). 

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have run your first workflow against the Google Pipelines API.

### Let's get started!


**Configuring a Google Project**

Install the <a href="https://cloud.google.com/sdk/downloads" target="_blank">Google Cloud SDK</a>. 
Create a <a href="https://cloud.google.com/resource-manager/docs/creating-managing-projects" target="_blank">Google Cloud Project</a> and give it a project id (e.g. sample-project).  We’ll refer to this as `<google-project-id>` and your user login (e.g. username@gmail.com) as `<google-user-id>`.  

On your Google project, open up the <a href="https://console.developers.google.com/apis/library" target="_blank">API Manager</a> and enable the following APIs:
        
* Google Compute Engine API
* Cloud Storage
* Google Cloud Life Sciences API

Authenticate to Google Cloud Platform  
`gcloud auth login <google-user-id>`

Set your default account (will require to login again)  
`gcloud auth application-default login`

Set your default project  
`gcloud config set project <google-project-id>`

Create a Google Cloud Storage (GCS) bucket to hold Cromwell execution directories.
We will refer to this bucket as `google-bucket-name`, and the full identifier as `gs://google-bucket-name`.  
`gsutil mb gs://<google-bucket-name>`  


**Workflow Source Files**

Copy over the sample `hello.wdl` and `hello.inputs` files to the same directory as the Cromwell jar. 
This workflow takes a string value as specified in the inputs file and writes it to stdout. 


***hello.wdl***
```
task hello {
  String addressee  
  command {
    echo "Hello ${addressee}! Welcome to Cromwell . . . on Google Cloud!"  
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello

  output {
     hello.message
  }
}
```

***hello.inputs***
```
{
  "wf_hello.hello.addressee": "World"
}
```

**Google Configuration File**

Copy over the sample `google.conf` file utilizing <a href="https://developers.google.com/identity/protocols/application-default-credentials" target="_blank">Application Default credentials</a> to the same directory that contains your sample WDL, inputs and Cromwell jar.
Replace `<google-project-id>` and `<google-bucket-name>`in the configuration file with the project id and bucket name. Replace `<google-billing-project-id>` with the project id that has to be billed for the request (more information for Requester Pays can be found at:
<a href="https://cloud.google.com/storage/docs/requester-pays" target="_blank">Requester Pays</a>)

***google.conf***
```
include required(classpath("application"))

google {

  application-name = "cromwell"

  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    }
  ]
}

engine {
  filesystems {
    gcs {
      auth = "application-default"
      project = "<google-billing-project-id>"
    }
  }
}

backend {
  default = “PAPIv2”
  providers {
    PAPIv2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        // Google project
        project = "<google-project-id>"

        // Base bucket for workflow executions
        root = "gs://<google-bucket-name>/cromwell-execution"

        // Polling for completion backs-off gradually for slower-running jobs.
        // This is the maximum polling interval (in seconds):
        maximum-polling-interval = 600

        // Optional Dockerhub Credentials. Can be used to access private docker images.
        dockerhub {
          // account = ""
          // token = ""
        }

        genomics {
          // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
          // Pipelines and manipulate auth JSONs.
          auth = "application-default"
          
          // Endpoint for APIs, which defaults to us-central1. To run with a location different from us-central1,
          // change the endpoint-url to start with the location, such as https://europe-west2-lifesciences.googleapis.com/
          endpoint-url = "https://lifesciences.googleapis.com/"
          
          // This allows you to use an alternative service account to launch jobs, by default uses default service account
          compute-service-account = "default"
		   
          // Cloud Life Sciences API is limited to certain locations. See https://cloud.google.com/life-sciences/docs/concepts/locations
          // and note that changing the location also requires changing the endpoint-url.
          location = "us-central1"	

          // Pipelines v2 only: specify the number of times localization and delocalization operations should be attempted
          // There is no logic to determine if the error was transient or not, everything is retried upon failure
          // Defaults to 3
          localization-attempts = 3
        }

        filesystems {
          gcs {
            // A reference to a potentially different auth for manipulating files via engine functions.
            auth = "application-default"
            project = "<google-billing-project-id>"
          }
        }
      }
    }
  }
}
```

**Run Workflow**

`java -Dconfig.file=google.conf -jar cromwell-67.jar run hello.wdl -i hello.inputs`

**Outputs**

The end of your workflow logs should report the workflow outputs.  

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_hello.hello.message": "Hello World! Welcome to Cromwell . . . on Google Cloud!"
  },
  "id": "08213b40-bcf5-470d-b8b7-1d1a9dccb10e"
}
```

Success!  

### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
## Getting started on AWS with AWS Batch (beta)

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have configured your local environment to run workflows using Cromwell on AWS Batch.

### Let's get started!

To create all the resources for running a Cromwell server on AWS using CloudFormation, launch the [Cromwell Full Stack Deployment](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/).  Alternatively, this page will walk through the specific steps to configure and run a local Cromwell server using AWS Batch.

1. [Authenticating a local Cromwell server with AWS](#authenticating-a-local-cromwell-server-with-aws)
2. [Configuring the AWS environment](#configuring-the-aws-environment)
3. [Configuring Cromwell](#configuring-cromwell)
4. [Workflow Source Files](#workflow-source-files)
5. [Running Cromwell and AWS](#running-cromwell-and-aws)
6. [Outputs](#outputs)

#### Authenticating a local Cromwell server with AWS

The easiest way to allow a local Cromwell server to talk to AWS is to:

1. Install the AWS CLI through Amazon's [user guide](https://docs.aws.amazon.com/cli/latest/userguide/installing.html).
2. [Configure the AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html) by calling `aws configure` (provide your `Access Key` and `Secret Access Key` when prompted).

Cromwell can access these credentials through the default authentication provider. For more options, see the [Configuring authentication of Cromwell with AWS](/backends/AWSBatch#configuring-authentication) section below.


#### Configuring the AWS environment

Next you'll need the following setup in your AWS account:
- The core set of resources (S3 Bucket, IAM Roles, AWS Batch)
- Custom Compute Resource (Launch Template or AMI) with Cromwell Additions

Information and instructions to setup an AWS environment to work properly with Cromwell can be found on [AWS for Genomics Workflow](https://docs.opendata.aws/genomics-workflows/core-env/introduction/). By deploying the CloudFormation templates provided by AWS, the stack will output the S3 bucket name and two AWS Batch queue ARNs (default and high-priority) used in the Cromwell configuration.




#### Configuring Cromwell

Now we're going to configure Cromwell to use the AWS resources we just created by updating a `*.conf` file to use the `AWSBackend` at runtime. This requires three pieces of information:

- The [AWS Region](https://docs.aws.amazon.com/general/latest/gr/rande.html) where your resources are deployed.
- S3 bucket name where Cromwell will store its execution files.
- The ARN of the AWS Batch queue you want to use for your tasks.

You can replace the placeholders (`<your region>`, `<your-s3-bucket-name>` and `<your-queue-arn>`) in the following config:

##### `aws.conf`

```hocon
include required(classpath("application"))

aws {

  application-name = "cromwell"
  auths = [
    {
      name = "default"
      scheme = "default"
    }
  ]
  region = "<your-region>"
}

engine {
  filesystems {
    s3.auth = "default"
  }
}

backend {
  default = "AWSBatch"
  providers {
    AWSBatch {
      actor-factory = "cromwell.backend.impl.aws.AwsBatchBackendLifecycleActorFactory"
      config {
        
        numSubmitAttempts = 6
        numCreateDefinitionAttempts = 6

        // Base bucket for workflow executions
        root = "s3://<your-s3-bucket-name>/cromwell-execution"

        // A reference to an auth defined in the `aws` stanza at the top.  This auth is used to create
        // Jobs and manipulate auth JSONs.
        auth = "default"

        default-runtime-attributes {
          queueArn: "<your arn here>"
        }

        filesystems {
          s3 {
            // A reference to a potentially different auth for manipulating files via engine functions.
            auth = "default"
          }
        }
      }
    }
  }
}

```

For more information about this configuration or how to change the behaviour of AWS Batch, visit the [AWS Backend](/backends/AWSBatch) page.


#### Workflow Source Files 

Lastly, create an example workflow to run. We're going to define a simple workflow that will `echo` a string to the console and return the result to Cromwell. Within AWS Batch (like other cloud providers), we're required to specify a Docker container for every task.

##### `hello.wdl`

```wdl
task hello {
  String addressee = "Cromwell"
  command {
    echo "Hello ${addressee}! Welcome to Cromwell . . . on AWS!"
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello
  output { hello.message }
}
```

#### Running Cromwell and AWS 

Provided all of the files are within the same directory, we can run our workflow with the following command:

> **Note**: You might have a different Cromwell version number here

```bash
java -Dconfig.file=aws.conf -jar cromwell-36.jar run hello.wdl
```

This will:
1. Start Cromwell in `run` mode,
2. Prepare `hello.wdl` as a job and submit this to your AWS Batch queue. You can monitor the job within your [AWS Batch dashboard](https://console.aws.amazon.com/batch/home).
3. Run the job, write execution files back to S3, and report progress back to Cromwell.

#### Outputs

The end of your workflow logs should report the workflow outputs.

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_hello.hello.message": "Hello World! Welcome to Cromwell . . . on AWS!"
  },
  "id": "08213b40-bcf5-470d-b8b7-1d1a9dccb10e"
}
```

Success!

### Next steps

You might find the following tutorials and guides interesting to tackle next:

* [Server Mode](/tutorials/ServerMode)
* [AWS Batch Backend](/backends/AWSBatch)
* [Persisting Data Between Restarts](/tutorials/PersistentServer)
## Timing Diagrams

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Server Mode](ServerMode)

### Goals

At the end of this tutorial you'll have seen how to view timing diagram for a workflow as it's running, and interpreting the information on the timing diagram for a completed workflow.

### Let's get started

Want to see the tasks that ran in your workflow laid out in a handy timing diagram? Good news! That's exactly what's about to happen to you!

#### Preparing files

To get a good example of these timing diagrams, we're going to make a workflow that scatters a task across a few indices, and see how that gets represented. Open your favorite text editor, copy the following text in and save it in your cromwell directory as `timingWorkflow.wdl`:
```wdl
workflow timingWorkflow {
	scatter(i in range(15)) {
		call sleep { input: sleep_time = i }
	}
}

task sleep {
	Int sleep_time
	command {
		echo "I slept for ${sleep_time}"
		sleep ${sleep_time}
	}
	output {
		String out = read_string(stdout())
	}
}
```

In brief, this workflow will scatter 15 tasks, each one will sleep for a time proportional to their scatter index. Is that hard to imagine? Never mind, we'll see it in diagramatic form soon!

#### Submit to cromwell:

If it's not running already, start the cromwell server:
```sh
java -jar cromwell-29.jar server
```

Submit the workflow:
```sh
curl -X POST --header "Accept: application/json"\
	-v "localhost:8000/api/workflows/v1" \
	-F workflowSource=@timingWorkflow.wdl
```

Amongst the curl output you should see the workflow ID come back, eg:
```
[...] Workflow 8d18b845-7143-4f35-9543-1977383b7d2f submitted.
```

I can now enter the following address into my web browser (i.e. Chrome) and see the timing diagram for the workflow. You'll need to swap out my workflow ID for the one that you received (they're all randomly generated) 
```
http://localhost:8000/api/workflows/v1/8d18b845-7143-4f35-9543-1977383b7d2f/timing
```

Once your tasks complete you should see output like mine, displayed below:

![Timing diagram](timingDiagram.png)

### Next Steps

Nice work! Now you know how to investigate which tasks in your workflow are spending the most time, and which ones are holding up the rest of the workflow.
Perhaps after completing such an amazing feat you might find the following pages interesting:

* [Viewing Metadata](MetadataEndpoint)
* [Configuration Files](ConfigurationFiles)
**Frequent Errors**

<!--
### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntroduction)

### Goals

At the end of this tutorial you'll have seen some a few common error cases in Cromwell, and learnt how easy ways to address them.

### Let's get started

### Next Steps

TBD (you're already looking pretty good...!)

-->

_Drop us a line in the [Forum](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team) if you have a question._

\*\*\* **UNDER CONSTRUCTION** \*\*\*  
[![Pennywell pig in red wellies - Richard Austin Images](http://www.richardaustinimages.com/wp-content/uploads/2015/04/fluffyAustin_Pigets_Wellies-500x395.jpg)](http://www.richardaustinimages.com/product/pennywell-pigs-under-umbrella-2/)
## Persisting Data Between Restarts

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)
* [Docker](https://docs.docker.com/engine/installation/)

### Goals

Cromwell remembers everything it knows!

### Let's get started!

- Start the MySQL docker container with the following line:

```bash
docker run -p 3306:3306 --name NameOfTheContainer -e MYSQL_ROOT_PASSWORD=YourPassword -e MYSQL_DATABASE=DatabaseName -e MYSQL_USER=ChooseAName -e MYSQL_PASSWORD=YourOtherPassword -d mysql/mysql-server:5.5
```

- Update your `application.conf` file.

```hocon
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://localhost/DatabaseName?rewriteBatchedStatements=true&useSSL=false"
    user = "ChooseAName"
    password = "YourOtherPassword"
    connectionTimeout = 5000
  }
}
```
Add the line above, below the all other lines in your `application.conf`. Replace `"DatabaseName"`, `"ChooseAName"` and `"YourOtherPassword"` with the values you choose in step 2, preserving the double quotes.

Test it by running your server with the updated `application.conf`:
```bash
java -Dconfig.file=/path/to/application.conf/ -jar cromwell-[version].jar ...
```

### Next steps

You might find the following tutorials interesting to tackle next:

* [Server Mode](ServerMode)
**Introduction to Call Caching**

<!--
### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Server Mode](ServerMode)
* [Configuration Files](ConfigurationFiles)

### Goals

At the end of this tutorial you'll have seen how to store the results for the workflows that you've run, and seen how those results are re-used the next time that you run the same job.

### Let's get started

### Next Steps

TBD (you're already looking pretty good...!)

-->

_Drop us a line in the [Forum](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team) if you have a question._

\*\*\* **UNDER CONSTRUCTION** \*\*\*  
[![Pennywell pig in red wellies - Richard Austin Images](http://www.richardaustinimages.com/wp-content/uploads/2015/04/fluffyAustin_Pigets_Wellies-500x395.jpg)](http://www.richardaustinimages.com/product/pennywell-pigs-under-umbrella-2/)
## Containers

Containers are encapsulated environments that include an operating system, libraries, and software. For example, if you have a host machine running Centos, you can run an isolated container with Ubuntu 18.04. At a high level, it's useful to think of a container as a program or binary.
 
To promote reproducibility and portability, it's considered best practice to define containers for a WDL and CWL task to run in - this ensures that running the same task on a different system will run the exact same software. 

Docker images are the most common container format, but it is not advisable for certain systems to run Docker itself, and for this reason Cromwell can be configured to support a number of alternatives.

* [Prerequisites](#prerequisites)
* [Goals](#goals)
* [Specifying Containers in your Workflow](#specifying-containers-in-your-workflow)
* [Docker](#docker)
    * [Docker on a Local Backend](#docker-on-a-local-backend)
    * [Docker on Cloud](#docker-on-cloud)
    * [Docker on HPC](#docker-on-hpc)
* [Singularity](#singularity)
    * [Installation](#installation)
    * [Configuring Cromwell for Singularity](#configuring-cromwell-for-singularity)
        * [Local environments](#local-environments)
        * [Job schedulers](#job-schedulers)
    * [Without Setuid](#without-setuid)
    * [Singularity Cache](#singularity-cache)
* [udocker](#udocker)
    * [Installation](#installation-1)
    * [Configuration](#configuration)
    * [Caching](#caching)
* [Configuration in Detail](#configuration-in-detail)
    * [Enforcing container requirements](#enforcing-container-requirements)
    * [Docker Digests](#docker-digests)
    * [Docker Root](#docker-root)
    * [Docker Config Block](#docker-config-block)
* [Best Practices](#best-practices)
    * [Image Versions](#image-versions)
* [Notes](#notes)
    * [How does Cromwell know when a job or container has completed?](#how-does-cromwell-know-when-a-job-or-container-has-completed)
    * [Cromwell: Run-in-background](#cromwell-run-in-background)
* [Next Steps](#next-steps)

### Prerequisites
This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)
* [Configuration Files](ConfigurationFiles.md)
* Recommended: [Getting started on HPC clusters](HPCIntro.md)

### Goals

At the end of this tutorial, you'll become familiar with container technologies and how to configure Cromwell to use these independently, or with job schedulers.


### Specifying Containers in your Workflow

Containers are specified on a per-task level, this can be achieved in WDL by specifying a [`docker`](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#docker) tag in the `runtime` section. For example, the following script should run in the `ubuntu:latest` container:

```wdl
task hello_world {
    String name = "World"
    command {
        echo 'Hello, ${name}'
    }
    output {
        File out = stdout()
    }
    runtime {
        docker: 'ubuntu:latest'
    }
}

workflow hello {
    call hello_world
}
```

Similarly in CWL, you can specify a [`DockerRequirement`](https://www.commonwl.org/v1.0/CommandLineTool.html#DockerRequirement) inside the requirements section:

```cwl
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
    name:
        type: string
        default: "World"
        inputBinding:
          prefix: "Hello, "
outputs:
    out: stdout

requirements:
    DockerRequirement:
        dockerPull: "ubuntu:latest"
``` 
___
### Docker

[Docker](https://www.docker.com) is a popular container technology that is natively supported by Cromwell and WDL.


#### Docker on a Local Backend

On a single machine (laptop or server), no extra configuration is needed to allow docker to run, provided Docker is installed.

You can install Docker for Linux, Mac or Windows from [Docker Hub](https://hub.docker.com/search/?type=edition&offering=community)

#### Docker on Cloud

It is strongly advised that you provide a Docker image to tasks that will run on Cloud backends, and in fact most Cloud providers require it.

It might be possible to use an alternative container engine, but this is not recommended if Docker is supported.

#### Docker on HPC

Docker can allow running users to gain superuser privileges, called the [Docker daemon attack surface](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface). In HPC and multi-user environments, Docker recommends that "only trusted users should be allowed to control your Docker Daemon".

For this reason, this tutorial will also explore other technologies that support the reproducibility and simplicity of running a workflow that use docker containers; Singularity and udocker.

___

### Singularity

Singularity is a container technology designed for use on HPC systems in particular, while ensuring an appropriate level of security that Docker cannot provide.

#### Installation
Before you can configure Cromwell on your HPC system, you will have to install Singularity, which is documented [here](https://www.sylabs.io/guides/3.0/admin-guide/admin_quickstart.html#installation).
In order to gain access to the full set of features in Singularity, it is strongly recommended that Singularity is installed by root, with the `setuid` bit enabled, as is ([documented here](https://www.sylabs.io/guides/2.6/admin-guide/security.html#how-does-singularity-do-it)).
This likely means that you will have to ask your sysadmin to install it for you.
Because `singularity` ideally needs `setuid`, your admins may have some qualms about giving Singularity this privilege.
If that is the case, you might consider forwarding [this letter](https://www.sylabs.io/guides/3.0/user-guide/installation.html#singularity-on-a-shared-resource) to your admins.

If you are not able to get Singularity installed with these privileges, you can attempt a user install.
If this is the case, you will have to alter your Cromwell configuration to work in "sandbox" mode, which is explained in [this part](#without-setuid) of the documentation. 

#### Configuring Cromwell for Singularity

Once Singularity is installed, you'll need to modify the `config` block inside `backend.providers` in your Cromwell configuration. In particular, this block contains a key called `submit-docker`, which will contain a script that is run whenever a job needs to run that uses a Docker image. If the job does not specify a Docker image, the regular `submit` block will be used.

As the configuration will require more knowledge about your execution environment, see the local and job scheduler sections below for example configurations.

##### Local environments

On local backends, you have to configure Cromwell to use a different 
`submit-docker` script that would start Singularity instead of docker. 
Singularity requires docker images to be prefixed with the prefix `docker://`.

Using containers isolates the filesystem that the script is allowed to interact with, for that reason we'll bind in the current working directory as `${docker_cwd}`, and we'll use the container-specific script path `${docker_script}`. 

An example submit script for Singularity is:
```bash
singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
```

As the `Singularity exec` command does not emit a job-id, we must include the `run-in-background` tag within the the provider section in addition to the docker-submit script. As Cromwell watches for the existence of the `rc` file, the `run-in-background` option has the caveat that we require the Singularity container to successfully complete, otherwise the workflow might hang indefinitely.

To ensure reproducibility and an isolated environment inside the container, 
`--containall` is an **important** function. By default, Singularity will mount
the user's home directory and import the user's environment as well as some 
other things that make Singularity easier to use in an interactive shell. 
Unfortunately settings in the home directory and the user's environment may 
affect the outcome of the tools that are used. This means different users may
get different results. Therefore, to ensure reproducibility while using 
Singularity, the `--containall` flag should be used. This will make sure the 
environment is cleaned and the HOME directory is not mounted.

Putting this together, we have an example base configuration for a local environment:
```hocon
include required(classpath("application"))

backend {
    default: singularity
    providers: {
        singularity {
            # The backend custom configuration.
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

            config {
                run-in-background = true
                runtime-attributes = """
                  String? docker
                """
                submit-docker = """
                  singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
                """
            }
        }
    }
}
```

##### Job schedulers

To run Singularity on a job scheduler, the singularity command needs to be passed to the scheduler as a wrapped command.

For example, in SLURM, we can use the normal SLURM configuration as explained in the [SLURM documentation](../backends/SLURM), however we'll add a `submit-docker` block to execute when a task is tagged with a docker container. 

When constructing this block, there are a few things to keep in mind:
- Make sure Singularity is loaded (and in PATH). If `module` is installed for 
  example you can call `module load Singularity`. If the cluster admin has made
  a Singularity module available. Alternatively you can alter the `PATH` 
  variable directly or simply use `/path/to/singularity` 
  directly in the config.
- We should treat worker nodes as if they do not have stable access to the 
  internet or build access, so we will pull the container before the task is 
  submit to the cluster.
- It's a good idea to use a Singularity cache so that same images should only
  have to be pulled once. Make sure you set the `SINGULARITY_CACHEDIR` 
  environment variable to a location on the filesystem that is reachable by the
  worker nodes!
- If we are using a cache we need to ensure that submit processes started by
  Cromwell do not pull to the same cache at the same time. This may corrupt the
  cache. We can prevent this by implementing a filelock with `flock` and 
  pulling the image before the job is submitted. The flock and pull command 
  needs to be placed *before* the submit command so all pull commands are 
  executed on the same node. This is necessary for the filelock to work.
- As mentioned above the `--containall` flag is **important** for 
  reproducibility.

```
submit-docker = """
    # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
    # based on the users home.
    if [ -z $SINGULARITY_CACHEDIR ]; 
        then CACHE_DIR=$HOME/.singularity/cache
        else CACHE_DIR=$SINGULARITY_CACHEDIR
    fi
    # Make sure cache dir exists so lock file can be created by flock
    mkdir -p $CACHE_DIR  
    LOCK_FILE=$CACHE_DIR/singularity_pull_flock
    # Create an exclusive filelock with flock. --verbose is useful for 
    # for debugging, as is the echo command. These show up in `stdout.submit`.
    flock --verbose --exclusive --timeout 900 $LOCK_FILE \
    singularity exec --containall docker://${docker} \
    echo "successfully pulled ${docker}!"

    # Submit the script to SLURM
    sbatch \
      [...]
      --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
  """
```

Putting this all together, a complete SLURM + Singularity config might look like this: 

```
backend {
  default = slurm

  providers {
    slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"                                                                                     
      config {
        runtime-attributes = """
        Int runtime_minutes = 600
        Int cpus = 2
        Int requested_memory_mb_per_core = 8000
        String? docker
        """

        submit = """
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${out} \
              -e ${err} \
              -t ${runtime_minutes} \
              ${"-c " + cpus} \
              --mem-per-cpu=${requested_memory_mb_per_core} \
              --wrap "/bin/bash ${script}"
        """

        submit-docker = """
            # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
            # based on the users home.
            if [ -z $SINGULARITY_CACHEDIR ]; 
                then CACHE_DIR=$HOME/.singularity/cache
                else CACHE_DIR=$SINGULARITY_CACHEDIR
            fi
            # Make sure cache dir exists so lock file can be created by flock
            mkdir -p $CACHE_DIR  
            LOCK_FILE=$CACHE_DIR/singularity_pull_flock
            # Create an exclusive filelock with flock. --verbose is useful for 
            # for debugging, as is the echo command. These show up in `stdout.submit`.
            flock --verbose --exclusive --timeout 900 $LOCK_FILE \
            singularity exec --containall docker://${docker} \
            echo "successfully pulled ${docker}!"

            # Submit the script to SLURM
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${cwd}/execution/stdout \
              -e ${cwd}/execution/stderr \
              -t ${runtime_minutes} \
              ${"-c " + cpus} \
              --mem-per-cpu=${requested_memory_mb_per_core} \
              --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}

```

#### Without Setuid
In addition, if you or your sysadmins were not able to give `setuid` permissions to `singularity`, you'll have to modify the config further to ensure the use of sandbox images:

```
submit-docker = """
    [...]

    # Build the Docker image into a singularity image
    # We don't add the .sif file extension because sandbox images are directories, not files
    DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
    IMAGE=${cwd}/$DOCKER_NAME
    singularity build --sandbox $IMAGE docker://${docker}

    # Now submit the job
    # Note the use of --userns here
    sbatch \
      [...]
      --wrap "singularity exec --userns --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
"""
```

#### Singularity Cache
By default, Singularity will cache the Docker images you pull in `~/.singularity`, your home directory.

However, if you are sharing your Docker images with other users or have limited space in your user directory, you can redirect this caching location by exporting the `SINGULARITY_CACHEDIR` variable in your `.bashrc` or at the start of the `submit-docker` block.
```
export SINGULARITY_CACHEDIR=/path/to/shared/cache
```

For further information on the Singularity Cache, refer to the [Singularity 2 caching documentation](https://www.sylabs.io/guides/2.6/user-guide/build_environment.html#cache-folders) (this hasn't yet been updated for Singularity 3).

___

### udocker

[udocker](https://github.com/indigo-dc/udocker) is a tool designed to "execute simple docker containers in user space without requiring root privileges".

In essence, udocker provides a command line interface that mimics `docker`, and implements the commands using one of four different container backends:

* PRoot
* Fakechroot
* runC
* Singularity

#### Installation
udocker can be installed without any kind of root permissions. Refer to udocker's installation documentation [here](https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md) for more information.

#### Configuration

(As of [2019-02-18](https://github.com/indigo-dc/udocker/issues/112)) udocker does not support looking up docker container by digests, hence you'll have to make ensure `hash-lookup` is disabled. Refer to [this section](#docker-digests) for more detail.

To configure `udocker` to work in a local environment, you must tag the provider's configuration to `run-in-background` and update the `submit-docker` to use udocker:
```
run-in-background = true
submit-docker = """
    udocker run -v ${cwd}:${docker_cwd} ${docker} ${job_shell} ${docker_script}
"""
```

With a job queue like SLURM, you just need to wrap this script in an `sbatch` submission like we did with Singularity:

```
submit-docker = """
    # Pull the image using the head node, in case our workers don't have network access
    udocker pull ${docker}
    
    sbatch \
      -J ${job_name} \
      -D ${cwd} \
      -o ${cwd}/execution/stdout \
      -e ${cwd}/execution/stderr \
      -t ${runtime_minutes} \
      ${"-c " + cpus} \
      --mem-per-cpu=${requested_memory_mb_per_core} \
      --wrap "udocker run -v ${cwd}:${docker_cwd} ${docker} ${job_shell} ${docker_script}"
"""
```

#### Caching
udocker caches images in a single directory, which defaults to [`~/.udocker`](https://github.com/indigo-dc/udocker/blob/master/udocker.py#L137), meaning that caching is done on a per-user basis. 
However, like Singularity, if you want to share a cache with other users in your project,you you can override the location of the udocker cache directory either using:
* A config file [described here](https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md#9-configuration), containing a line such as `topdir = "/path/to/cache"`.
* Using the environment variable `$UDOCKER_DIR`

___

### Configuration in Detail
The behaviour of Cromwell with containers can be modified using a few other options.

#### Enforcing container requirements
You can enforce the use of a container by not including the `submit` block in the provider section.

However note that some interpolated variables (`${stdout}`, `${stderr}`) are different between these two blocks.

#### Docker Digests

Each Docker repository has a number of tags that can be used to refer to the latest image of a particular type. 
For instance, when you run a normal Docker image with `docker run image`, it will actually run `image:latest`, the `latest` tag of that image.

However, by default Cromwell requests and runs images using their `sha` hash, rather than using tags.
This strategy is actually preferable, because it ensures every execution of the task or workflow will use the exact same version of the image, but some engines such as `udocker` don't support this feature.

If you are using `udocker` or want to disable the use of hash-based image references, you can set the following config option:
```
docker.hash-lookup.enabled = false
```

Nb: By disabling hash-lookup, call caching will not work for any container using a floating tag.

#### Docker Root
If you want to change the root directory inside your containers, where the task places input and output files, you can edit the following option:

```
backend {
  providers {
    LocalExample {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
      
        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above.
        dockerRoot = "/cromwell-executions"
      }
    }
  }
}
```


#### Docker Config Block
Further docker configuration options available to be put into your config file are as follows. 
For the latest list of parameters, refer to the [example configuration file][cromwell-examples-conf],
and [specific backend provider examples][cromwell-examples-folder].

```
docker {
  hash-lookup {
    # Set this to match your available quota against the Google Container Engine API
    #gcr-api-queries-per-100-seconds = 1000

    # Time in minutes before an entry expires from the docker hashes cache and needs to be fetched again
    #cache-entry-ttl = "20 minutes"

    # Maximum number of elements to be kept in the cache. If the limit is reached, old elements will be removed from the cache
    #cache-size = 200

    # How should docker hashes be looked up. Possible values are "local" and "remote"
    # "local": Lookup hashes on the local docker daemon using the cli
    # "remote": Lookup hashes on docker hub, gcr, gar, quay
    #method = "remote"
  }
}
```


### Best Practices

#### Image Versions

When choosing the image version for your pipeline stages, it is highly recommended that you use a hash rather than a tag, for the sake of reproducibility
For example, in WDL, you could do this:
```wdl
runtime {
    docker: 'ubuntu:latest'
}
```

But what you should do is this:
```wdl
runtime {
    docker: 'ubuntu@sha256:7a47ccc3bbe8a451b500d2b53104868b46d60ee8f5b35a24b41a86077c650210'
}
```

You can find the `sha256` of an image using `docker images --digests`
 
 
### Notes

#### How does Cromwell know when a job or container has completed?
Cromwell uses the presence of the `rc` (returncode) file to determine whether a task has succeeded or failed. This `rc` file is generated as part of the `script` within the execution directory, where the script is assembled at runtime. This is important as if the script executes successfully but the container doesn't terminate, Cromwell will continue the execution of the workflow and the container will persist hogging system resources.

Within the configurations above:
- `singularity`: The exec mode does not run a container on the background

#### Cromwell: Run-in-background

By enabling Cromwell's run-in-background mode, you remove the necessity for the `kill`, `check-alive` and `job-id-regex` blocks, which disables some safety checks when running workflows:

- If there is an error starting the container or executing the script, Cromwell may not recognise this error and hang. For example, this may occur if the container attempts to exceed its allocated resources (runs out of memory); the container daemon may terminate the container without completing the script.
- If you abort the workflow (by attempting to close Cromwell or issuing an abort command), Cromwell does not have a reference to the container execution and will not be able to terminate the container.

This is only necessary in local environments where there is no job manager to control this, however if your container technology can emit an identifier to stdout, then you are able to remove the run-in-background flag. 

### Next Steps

Congratulations for improving the reproducibility of your workflows! You might find the following cloud-based tutorials interesting to test your workflows (and ensure the same results) in a completely different environment:

- [Getting started with AWS Batch](AwsBatch101.md)
- [Getting started on Google Pipelines API](PipelinesApi101.md)
- [Getting started on Alibaba Cloud](BCSIntro/)


[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends
# Server Mode

Cromwell is best experienced in "server" mode, as discussed in the [Modes section of the docs](../Modes).

# Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)

# Goals

At the end of this tutorial you'll have run Cromwell in Server Mode, allowing you to submit more than one workflow in parallel, and query workflows even after they have completed.

# Prepare Files

Paste the following into  a file called `hello.wdl`:
```wdl
task hello {
  String name

  command {
    echo 'Hello ${name}!'
  }
  output {
    File response = stdout()
  }
}

workflow test {
  call hello
}
```

Paste the following into a file called `inputs.json`:
```json
{
  "test.hello.name": "World"
}
```

# Run the server

1. Run `java -jar cromwell-[version].jar server` (replace [version] with actual version).  Note that there is a `server` argument, this is the special sauce!
2. Visit <a href="http://localhost:8000">localhost:8000</a> in your browser

# Start the job

1. Navigate to Workflows section and click "Show/Hide"  
![](workflows.png)
2. Navigate to `/workflows/{version}` which has a green "POST" on the left.  
![](submit.png)
3. Find workflowSource file, "Choose File" and navigate to `hello.wdl`.  
![](workflowSource.png)  
4. Find inputs file and navigate to `inputs.json`.  
![](inputs.png)  
5. Navigate to the bottom of this section and click "Try it out!"  
![](try.png)
6. Observe output from the server process.

# What happened?

* [Did it work?  Check the status endpoint](../api/RESTAPI#api-workflows-version-id-status-get)
* Holy logs!  [How can I just see my outputs?](../api/RESTAPI#api-workflows-version-id-outputs-get)
* [Check out metadata related to your workflow.](../api/RESTAPI#api-workflows-version-id-metadata-get)
* [All kinds of other interesting info.](../api/RESTAPI)

### Next Steps

Nice job! Now that you have cromwell running in server mode you've reached the upper echilons of Cromwell prowess! After reaching these dizzy heights, you might also find the following pages interesting:

* [Viewing Metadata](MetadataEndpoint)
* [Timing Diagrams](TimingDiagrams)
* [Configuration Files](ConfigurationFiles)
**General Debugging Tips**

<!--
### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [FrequentErrors](FrequentErrors)

### Goals

At the end of this tutorial you'll have learnt some general methods to debug problems that you see in your Cromwell running in server mode.

### Let's get started

### Next Steps

TBD (you're already looking pretty good...!)
-->

_Drop us a line in the [Forum](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team) if you have a question._

\*\*\* **UNDER CONSTRUCTION** \*\*\*  
[![Pennywell pig in red wellies - Richard Austin Images](http://www.richardaustinimages.com/wp-content/uploads/2015/04/fluffyAustin_Pigets_Wellies-500x395.jpg)](http://www.richardaustinimages.com/product/pennywell-pigs-under-umbrella-2/)
## Getting started on HPC clusters

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)

### Goals

At the end of this tutorial you'll have set up Cromwell to run against your HPC cluster. We'll use SGE as an example but this applies equally to LSF and others.

### Let's get started!

####

#### Telling Cromwell the type of backend

Start by defining your new backend configuration under the section `backend`. For now, we'll give your backend the name `SGE`, but you can use any name you would like.

```hocon
backend {
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        # to be filled in
      }
    }
  }
}
```

The `actor-factory` above tells cromwell that you will be using the `config` section to tell cromwell how to submit jobs, abort jobs, etc.

You'll likely also want to change the default backend to your new backend, by setting this configuration value:

```hocon
backend.default = SGE
```

#### Specifying the runtime attributes for your HPC tasks

In the config section for your backend, you can define the different [runtime attributes](../RuntimeAttributes) that your HPC tasks will support. Any runtime attribute configured here will be read from the WDL tasks, and then passed into the command line used to submit jobs to the HPC cluster.

All runtime attributes must be defined in a single multi-line block. The syntax of this block is the same as defining the inputs for a WDL task.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Int cpu = 1
  Float? memory_gb
  String? sge_queue
  String? sge_project
  """
}
```

In the example above, we have defined four different WDL variables defined, `cpu`, `memory_gb`, `sge_queue`, and `sge_project`. Below you will find more information on `cpu` and `memory`, and the ability to add custom runtime attributes like the `sge_queue` and `sge_project`.

**cpu**

When you declare a runtime attribute with the name `cpu`, it must be an `Int`. This integer will validated to always be `>= 1`.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Int cpu = 1
  # ...
  """
}
```

**memory**

When running a workflow, the memory runtime attribute in the task will specify the units of memory. For example, this jobs specifies that it only needs 512 megabytes of memory when running.

```wdl
task hello {
  command { echo hello }
  runtime { memory: "512 MB" }
}
```

However, it's possible that when submitting jobs to your HPC cluster you want to specify the units in gigabytes.

To specify the memory units that the submit command should use, append the units to the memory runtime attribute. For example:

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Float? memory_gb
  # ...
  """
}
```

Now, no matter what unit of memory is used within the task, the value will be converted into gigabytes before it is passed to your submit command.

**custom attributes**

You can also declare other runtime attributes that a WDL task may use. For example, suppose you would like to allow the WDL to specify an sge queue in a task, like:

```wdl
task hello {
  command { echo hello }
  runtime { sqe_queue: "short" }
}
```

You declare your runtime attribute in your config by adding any other custom value to the `runtime-attributes` section:

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  String? sge_queue
  # ...
  """
}
```

In this case, we've stated that the `sge_queue` is optional. This allows us to reuse WDLs from other pipeline authors who may not have set an `sge_queue`.

Alternatively, you can also set a default for the declared runtime attributes.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  String sge_queue = "short"
  # ...
  """
}
```

##### Call Caching based on runtime attributes

The rules for call caching in HPC backends are:
* `docker`: Will be considered when call caching.
* Memory options: Will *not* be considered when call caching.
* CPU options: Will *not* be considered when call caching.
* Custom Attributes: Will *not* be considered when call caching (by default).
 
Although custom attributes will not be considered when call caching by default, you can override this in a `runtime-attributes-for-caching` section. Eg:
```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  String sge_queue = "short"
  String singularity_image
  # ...
  """
  runtime-attributes-for-caching {
    sge_queue: false
    singularity_image: true
  }
}
```

* Note: Only *custom* attributes can be altered like this. Memory, CPU and docker will always have their default cache-consideration behavior.
* Note: Unlike memory, cpu and docker attributes which inherit validation and hash-lookup behavior, any custom attributes will be simple primitive comparisons.
    * For example, a `docker` attribute will be cached by looking up docker hashes against a docker repository, but a custom `singularity` attribute would be a primitive string match.

#### How Cromwell should start an HPC job

When Cromwell runs a task, it will fill in a template for the job using the declared runtime attributes. This specific template will vary depending on the requirements of your HPC cluster. For example, say you normally submit jobs to SGE using:

```bash
qsub -terse -V -b y -N my_job_name \
  -wd /path/to/working_directory \
  -o /path/to/stdout.qsub \
  -e /path/to/stderr.qsub \
  -pe smp 1 -l mem_free=0.5g -q short \
  /usr/bin/env bash myScript.bash
```

For this particular SGE cluster, the above sets the working directory, stdout and stderr paths, the number of cpus to 1, the memory to half a gigabyte, and runs on the short queue.

Converting this into a template using our runtime attributes requires defining `submit` as one would a WDL task `command`:

```hocon
backend.providers.SGE.config {
  submit = """
  qsub \
  -terse \
  -V \
  -b y \
  -N ${job_name} \
  -wd ${cwd} \
  -o ${out}.qsub \
  -e ${err}.qsub \
  -pe smp ${cpu} \
  ${"-l mem_free=" + memory_gb + "g"} \
  ${"-q " + sge_queue} \
  ${"-P " + sge_project} \
  /usr/bin/env bash ${script}
  """
}
```

When the job finishes submitting, Cromwell will need to retrieve the job id, so that it can abort the job if necessary. This job should be written to the stdout after submission, where Cromwell will then read the job id. Because the job id may be surrounded by other text, a custom regular expression should capture the actual job id. Because the submit above uses `-terse`, the job id will be the entire contents of the stdout, but should be all digits:

```hocon
backend.providers.SGE.config {
  job-id-regex = "(\\d+)"
}
```

#### How Cromwell should abort an HPC job

When aborting an HPC job, Cromwell will run a command confifured under the key `kill`, passing in the WDL variable `job_id`:

```hocon
backend.providers.SGE.config {
  kill = "qdel ${job_id}"
}
```

#### How Cromwell checks if an HPC job is alive

Whenever Cromwell restarts it checks to see if a job has completed by searching for return code in a file called `rc`. If this file isn't available, in this case Cromwell runs an extra check to make sure the job is still alive. You can configure the command used for this check via:

```hocon
backend.providers.SGE.config {
  check-alive = "qstat -j ${job_id}"
}
```

#### Other backend settings

On some systems, the administrators may limit the number of HPC jobs a user may run at a time. To configure this limit, you can use the value `concurrent-job-limit` to limit the number of jobs.

```hocon
backend.providers.SGE.config {
  concurrent-job-limit = 100
}
```

#### Putting the config section all together

With the above sections, we can combine them all together to create a completly working HPC backend.

```hocon
backend {
  default = SGE
  
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        concurrent-job-limit = 100
    
        runtime-attributes = """
        Int cpu = 1
        Float? memory_gb
        String? sge_queue
        String? sge_project
        """
    
        submit = """
        qsub \
        -terse \
        -V \
        -b y \
        -N ${job_name} \
        -wd ${cwd} \
        -o ${out} \
        -e ${err} \
        -pe smp ${cpu} \
        ${"-l mem_free=" + memory_gb + "g"} \
        ${"-q " + sge_queue} \
        ${"-P " + sge_project} \
        /usr/bin/env bash ${script}
        """

        job-id-regex = "(\\d+)"

        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
      }
    }
  }
}
```

Running Cromwell with this in our configuration file will now submit jobs to SGE!

### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
## Configuration Files

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)


### Goals

At the end of this tutorial you'll have set up a configuration file for Cromwell and used it to modify Cromwell's behavior.

### Let's get started

#### Customizing Cromwell with Configuration Files

When Cromwell runs, it contains a large number of default options useful for getting started. For example, by default Cromwell doesn't require an external database while running all workflow jobs on your local machine.

Soon you may want to start storing the results of your Cromwell runs in an external MySQL database. Or, you may want to run jobs on your organizations compute farm, or even run jobs in the cloud via the Pipelines API. All of these changes to the defaults will be done by setting configuration values.

When you have many configuration settings you would like to set, you specify them in a custom configuration file. See the [configuration](../Configuring) page for more specific information on the configuration file, and for links to the example configuration file.

#### Configuration file syntax

Cromwell configuration files are written in a syntax called HOCON. See the [HOCON documentation](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation) for more information on all the ways one can create a valid configuration file.

#### Creating your first configuration file

To get started customizing Cromwell via a configuration file, create a new empty text file, say `your.conf`. Then add this include at the top:

```hocon
include required(classpath("application"))
```

The default Cromwell configuration values are set via Cromwell's `application.conf`. To ensure that you always have the defaults from the `application.conf`, you must include it at the top of your new configuration file.

#### Running Cromwell with your configuration file

Once you have created a new configuration file, you can pass the path to Cromwell by setting the system property `config.file`:

```bash
java -Dconfig.file=/path/to/your.conf -jar cromwell-[VERSION].jar server
```

Cromwell should start up as normal. As you haven't actually overridden any values yet, Cromwell should be running with the same settings.

#### Setting a configuration value

To override a configuration value, you can specify new values in your configuration file. For example, say you want to change the default port that cromwell listens from `8000` to `8080`. In your config file you can set:

```hocon
# below the include line from before
webservice {
  port = 8080
}
```

When you then run Cromwell updated config file, cromwell will now be listening on 8080 or 8000.

#### Finding more configuration properties

In addition to the common configuration properties listed on the [configuration](../Configuring) page, there are also a large number of example configuration stanzas commented in [cromwell.examples.conf][cromwell-examples-conf], and
backend provider examples in [cromwell.example.backends][cromwell-examples-folder].


[cromwell.examples.conf](https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf).

### Next Steps

After completing this tutorial you might find the following pages interesting:

* [Configuring the Local Backend](LocalBackendIntro)
* [Server Mode](ServerMode.md)

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends
In order to run a workflow, Cromwell uses the backends available to it to create jobs and monitor them until they are complete. When no new jobs can start and all jobs have finished (Success, Failure, Aborted), the workflow terminates.

In an ideal situation, all jobs succeed and the workflow succeeds.
Unfortunately things don't always go as planned. Here is what to expect from Cromwell when things go off the rails.

## Failure Modes

Cromwell supports two failure modes, which specify how Cromwell behaves when a job fails during the execution of a workflow.

* `NoNewCalls` **(default)**  
	* Cromwell does not start any new call as soon as a job fails. Cromwell will still monitor the rest of the jobs until they complete (successfully or not).  
* `ContinueWhilePossible`  
	* Cromwell attempts to run as many jobs as possible until no more can be started. When all running jobs are complete, the workflow fails.

The failure mode can be set in the [Configuration](../Configuring/) or [Workflow Options](../wf_options/Overview#workflow-failure).

_For example:_

![](ABdependency.png)

This simple diagram represents 4 jobs:

* Job A and Job B are independent; 
* Job A1 depends on A; 
* Job B1 depends on B.

Let's look at the case where A and B are both running, and for some reason B fails (shaded red).

**`NoNewCalls`**  
If the failure mode is `NoNewCalls` Cromwell waits for A to complete (shaded green). Regardless of whether A is successful or not, Cromwell then fails the workflow, without starting A1 or B1 (shaded grey).

![](NNC_B_fail.png)

**`ContinueWhilePossible`**  
If the failure mode is `ContinueWhilePossible` and A succeeds (green), then Cromwell starts A1 (green) and waits for it to complete. At this point all jobs that can run have completed. B1 (grey) cannot run since B failed (red), therefore Cromwell fails the workflow without starting B1.

![](CWP_B_fail.png)

### Retryable failures

Retryable failures are **not** failures that can trigger a workflow failure. An example of a retryable failure is when a [preemptible VM](../RuntimeAttributes/#preemptible) is preempted. 

In the example above, if B's failure is retryable then B will be retried (shaded yellow and green). The workflow will keep running normally, regardless of which failure mode is enabled.

![](CWP_B_retryable_fail_then_success.png)

Using the previous example, let's imagine **B** failed from a **non-retryable** failure (shaded red). After B failed, **A** fails from a **retryable** failure. Now Cromwell's behavior will depend on the failure mode.

**`NoNewCalls`**  
If the failure mode is `NoNewCalls`, then A **will not be** retried (yellow). A1 (grey) will not start, because A did not complete.

![](NCC_B_fail_A_retryable.png)

**`ContinueWhilePossible`**  
If the failure mode is **`ContinueWhilePossible`**, then A **will be** retried (yellow and green). If A is successful then A1 will start (green).

![](CWP_B_fail_A_retryable.png)

## Abort

In both [Run](../Modes/#run) and [Server](../Modes/#server) mode, you can abort a running workflow. This section explains what that entails.

When aborting a workflow, either through the [abort endpoint](../api/RESTAPI#abort-a-running-workflow) or by terminating the [Cromwell run process](../Modes) (if [configured](../Configuring#abort) to do so), Cromwell does the following:

1. Changes the status of the workflow to `Aborting`,
2. Does not start any new jobs,
3. Asks every running job to abort,
4. Waits for all running jobs to complete,
5. Finalizes the workflow,
6. Changes the status of the workflow to `Aborted`.

The action of aborting a job is backend specific. Cromwell can only ask a backend to abort a job and wait for the backend to notify it when it is aborted.  

For example, if you are running Cromwell on the Google backend and abort a job, Google will send an abort request to the Google Pipelines API. When Pipelines API indicates that the status of the job is aborted, Cromwell will mark it as such.
Remember that abort action is entirely dependent on the backend. In this particular case, Pipelines API does not guarantee the success of an abort request (see [Pipelines API documentation on abort](https://cloud.google.com/genomics/reference/rest/v1alpha2/operations/cancel)).

You'll also notice that the workflow is finalized even though being aborted.
Finalization is the last step in the execution of a workflow and a chance for each backend to do some work before the workflow is terminated.
Backends won't be denied the chance to finalize the workflow even if it's being aborted.

_Note that by the time the backend is asked to abort a job, the job may have succeeded or failed already. In this case Cromwell will report the job's status (successful or failed)._

_If a job fails with a retryable failure (e.g is preempted), it will **not** be attempted again when the workflow is aborting._

## Restart

When Cromwell restarts (for example to upgrade to a new version) it will reconnect to all workflows that were in progress. On the Google and HPC backends only, Cromwell will additionally attempt to reconnect to all running jobs. Note that a workflow
does not "belong" to any one Cromwell instance (it belongs to the cluster), so a different instance in a horizontal cluster might reconnect to the workflow instead of the original.


If the workflow was in state `Aborting`, Cromwell will ask all running jobs to abort again. No new jobs will be started.

Once all jobs have been reconnected to, the workflow will keep running normally.

During the reconnection process Cromwell might ask backends to reconnect to jobs that were never started before the restart. In that case, the job will be mark as failed with an explanation message. This failure is benign and only an artifact of the fact that Cromwell was restarted.  
If the backend does not support reconnection to an existing job, jobs will be marked as failed with an explanation message as well. The backend status of the jobs will be "Unknown".

## Graceful Shutdown

When Cromwell is run as a server, it will by default attempt to gracefully shutdown, stopping its different services in a specific order to avoid losing critical data.
This behavior, documented below, can be turned off in the configuration via `system.graceful-server-shutdown = false`.

Upon receiving a `SIGINT` or `SIGTERM` signal, the JVM will initiate its shutdown process. Prior to this Cromwell will attempt to shutdown its own services in the following way:

1. Workflows in `Submitted` state are no longer started
2. Cromwell unbinds from the address/port it was listening on. From this point the Cromwell server is unreachable via the endpoints.
3. All actors generating data that needs to be persisted receive a message asking them to gracefully stop.
This means that they are given some time (see below for how much and how to change it) to return to a known "consistent" state.
For example, an actor waiting for a response from the database before sending information to the metadata will wait for that response before shutting itself down.
4. All active connections from the REST endpoints are completed and closed. At this point any client that made a request before the shutdown process started should have received a response.
5. All actors responsible for data persistence are in turn being asked to gracefully shutdown. 
For example, all queued up metadata writes are executed.
6. Database connection pools are shutdown.
7. Actor system shuts down.
8. JVM exits.
    
This multi-stage process is designed to minimize the risk of data loss during shutdown. However in order to prevent this process from lasting forever, each stage (called phase) has its own timeout.
If the phase does not complete within the given timeout, actors will be forcefully stopped and the next phase will start.

This logic is implemented using [Akka Coordinated Shutdown Extension](http://doc.akka.io/docs/akka/current/scala/actors.html#coordinated-shutdown). Currently Cromwell is running [version 2.5.4](https://doc.akka.io/docs/akka/2.5.4/scala/actors.html#coordinated-shutdown).
It comes with a set of pre-defined phases, that can be added on and modified. Those phases can be linked together to form a Graph. Cromwell shutdown graphs looks as such:

![Scaladoc](CromwellShutdownProcess.png)

Pre-defined but unused phases have been omitted (cluster related phases for example that are irrelevant in Cromwell).

You'll notice the presence of a `PhaseAbortAllWorkflows` phase. This phase is at the same level as the `PhaseServiceRequestsDone` phase which corresponds to our step #3 above.
The reason for a specific abort phase is so that its timeout can be configured differently than the normal shutdown phase.

Indeed, stopping all workflows and aborting them is very similar from an outside perspective. We send a message (resp. "Stop" and "Abort") and wait for a response.

In the case where you want to be able to give more time to abort, as it will likely involve more work, you can edit the value of `coordinated-shutdown.phases.abort-all-workflows.timeout` which defaults to 1 hour.
Phases timeouts default to 5 seconds, except the stop-io-activity phase which defaults to 30 minutes. This is because depending on the Database load at the time of the shutdown, it might take a significant amount of time to flush all pending writes.

All of the timeouts are configurable in the `akka.coordinated-shutdown.phases` section ([see the latest `reference.conf`](https://raw.githubusercontent.com/akka/akka/master/akka-actor/src/main/resources/reference.conf)).
To change the default timeout, change the value of `akka.coordinated-shutdown.default-phase-timeout`.
Requests that Cromwell can't process return a failure in the form of a JSON response respecting the following JSON schema:

```
{
  "$schema": "http://json-schema.org/draft-04/schema#",
  "description": "Error response schema",
  "type": "object",
  "properties": {
    "status": {
      "enum": [ "fail", "error"]
    },
    "message": {
      "type": "string"
    },
    "errors": {
      "type": "array",
      "minItems": 1,
      "items": { "type": "string" },
      "uniqueItems": true
    }
  },
  "required": ["status", "message"]
}
```

The `status` field can be `"fail"` or `"error"`.  

`"fail"` means that the request was invalid and/or data validation failed. `"fail"` status is most likely returned with a 4xx HTTP Status code.  

*For example,*

```
{
  "status": "fail",
  "message": "Workflow input processing failed.",
  "errors": [
    "Required workflow input 'helloworld.input' not specified."
  ]
}
```

`"error"` means that an error occurred while processing the request. `"error"` status is most likely returned with a 5xx HTTP Status code.  

*For example,*

```
{
  "status": "error",
  "message": "Connection to the database failed."
}
```

The `message` field contains a short description of the error.

The `errors` field is optional and may contain additional information about why the request failed.<!--
This file was generated by `sbt generateRestApiDocs` on Thu, 25 Mar 2021 19:28:57 -0400

!!! DO NOT CHANGE THIS FILE DIRECTLY !!!

If you wish to change something in this file, either change cromwell.yaml or GenerateRestApiDocs.scala then
regenerate.
-->
# Cromwell Server REST API


<a name="overview"></a>
**Overview**  
Describes the REST API provided by a Cromwell server


**Version information**  
*Version* : 30


**License information**  
*License* : BSD  
*License URL* : https://github.com/broadinstitute/cromwell/blob/develop/LICENSE.txt  
*Terms of service* : null


**Produces**  

* `application/json`





<a name="describe"></a>
## Machine-readable description of a workflow, including inputs and outputs
```
POST /api/womtool/{version}/describe
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**FormData**|**workflowInputs**  <br>*optional*|JSON or YAML file containing the inputs as an object.|file||
|**FormData**|**workflowSource**  <br>*optional*|The workflow source file to submit for execution. Either workflow source or workflow url is required.|file||
|**FormData**|**workflowType**  <br>*optional*|The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.|enum (WDL, CWL)||
|**FormData**|**workflowTypeVersion**  <br>*optional*|The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.|enum (draft-2, 1.0, v1.0)||
|**FormData**|**workflowUrl**  <br>*optional*|URL which points to the workflow. Either workflow source or workflow url is required.|string||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Workflow description.|[WorkflowDescription](#workflowdescription)|


#### Consumes

* `multipart/form-data`


#### Tags

* Womtool


<a name="submit"></a>
## Submit a workflow for execution
```
POST /api/workflows/{version}
```


#### Description
Submits a workflow to Cromwell. Note that this endpoint can accept an unlimited number of input files via workflowInputs_N but swagger needs them to be explicitly defined so we have provided 5 as an example.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**FormData**|**labels**  <br>*optional*|JSON object of labels to apply to this workflow.|file||
|**FormData**|**workflowDependencies**  <br>*optional*|ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to this workflow.|file||
|**FormData**|**workflowInputs**  <br>*optional*|JSON or YAML file containing the inputs as an object. For WDL workflows a skeleton file can be generated from WOMtool using the "inputs" subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.|file||
|**FormData**|**workflowInputs_2**  <br>*optional*|A second JSON or YAML file containing inputs.|file||
|**FormData**|**workflowInputs_3**  <br>*optional*|A third JSON or YAML file containing inputs.|file||
|**FormData**|**workflowInputs_4**  <br>*optional*|A fourth JSON or YAML file containing inputs.|file||
|**FormData**|**workflowInputs_5**  <br>*optional*|A fifth JSON or YAML file containing inputs.|file||
|**FormData**|**workflowOnHold**  <br>*optional*|Put workflow on hold upon submission. By default, it is taken as false.|boolean||
|**FormData**|**workflowOptions**  <br>*optional*|JSON file containing configuration options for the execution of this workflow.|file||
|**FormData**|**workflowRoot**  <br>*optional*|The root object to be run. Only necessary for CWL submissions containing multiple objects (in an array).|string||
|**FormData**|**workflowSource**  <br>*optional*|The workflow source file to submit for execution. Either workflow source or workflow url is required.|file||
|**FormData**|**workflowType**  <br>*optional*|The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.|enum (WDL, CWL)||
|**FormData**|**workflowTypeVersion**  <br>*optional*|The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.|enum (draft-2, 1.0, v1.0)||
|**FormData**|**workflowUrl**  <br>*optional*|URL which points to the workflow. Either workflow source or workflow url is required.|string||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**201**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Invalid submission request|No Content|
|**500**|Internal Error|No Content|


#### Consumes

* `multipart/form-data`


#### Tags

* Workflows


<a name="backends"></a>
## List the supported backends
```
GET /api/workflows/{version}/backends
```


#### Description
Returns the backends supported by this Cromwell server, as well as the default backend.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[BackendResponse](#backendresponse)|


#### Tags

* Workflows


<a name="submitbatch"></a>
## Submit a batch of workflows for execution
```
POST /api/workflows/{version}/batch
```


#### Description
In instances where you want to run the same workflow multiple times with varying inputs you may submit a workflow batch. This endpoint is fundamentally the same as the standard submission endpoint with the exception that the inputs JSON will be an array of objects instead of a single object.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**FormData**|**labels**  <br>*optional*|JSON object of labels to apply to this workflow.|file||
|**FormData**|**workflowDependencies**  <br>*optional*|ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to these workflows.|file||
|**FormData**|**workflowInputs**  <br>*required*|JSON file containing the inputs as an array of objects. Every element of the array will correspond to a single workflow. For WDL workflows a skeleton file can be generated from WOMtool using the "inputs" subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.|file||
|**FormData**|**workflowOnHold**  <br>*optional*|Put workflow on hold upon submission. By default, it is taken as false.|boolean||
|**FormData**|**workflowOptions**  <br>*optional*|JSON file containing configuration options for the execution of this workflow.|file||
|**FormData**|**workflowSource**  <br>*optional*|The workflow source file to submit for execution. Either workflow source or workflow url is required.|file||
|**FormData**|**workflowType**  <br>*optional*|The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.|enum (WDL, CWL)||
|**FormData**|**workflowTypeVersion**  <br>*optional*|The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.|enum (draft-2, 1.0, v1.0)||
|**FormData**|**workflowUrl**  <br>*optional*|URL which points to the workflow. Either workflow source or workflow url is required.|string||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|< [WorkflowIdAndStatus](#workflowidandstatus) > array|
|**400**|Malformed Workflow ID|No Content|
|**500**|Internal Error|No Content|


#### Consumes

* `multipart/form-data`


#### Tags

* Workflows


<a name="callcachediff"></a>
## Explain hashing differences for 2 calls
```
GET /api/workflows/{version}/callcaching/diff
```


#### Description
This endpoint returns the hash differences between 2 completed (successfully or not) calls.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**Query**|**callA**  <br>*required*|Fully qualified name, including workflow name, of the first call.|string||
|**Query**|**callB**  <br>*required*|Fully qualified name, including workflow name, of the second call|string||
|**Query**|**indexA**  <br>*optional*|Shard index for the first call for cases where the requested call was part of a scatter.|integer||
|**Query**|**indexB**  <br>*optional*|Shard index for the second call for cases where the requested call was part of a scatter.|integer||
|**Query**|**workflowA**  <br>*required*|Workflow Id of the first workflow|string||
|**Query**|**workflowB**  <br>*required*|Workflow Id of the second workflow|string||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**404**|No matching cache entry. Cromwell versions prior to 28 will not have recorded information necessary for this endpoint and thus will also appear to not exist.|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="querypost"></a>
## Get workflows matching some criteria
```
POST /api/workflows/{version}/query
```


#### Description
Query workflows by start dates, end dates, names, ids, labels, or statuses.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**Body**|**parameters**  <br>*required*|Same query parameters as GET /query endpoint, submitted as a json list. Example: [{"status":"Success"},{"status":"Failed"}]|< [WorkflowQueryParameter](#workflowqueryparameter) > array||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowQueryResponse](#workflowqueryresponse)|
|**400**|Malformed Workflow ID|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="queryget"></a>
## Get workflows matching some criteria
```
GET /api/workflows/{version}/query
```


#### Description
Query for workflows which match various criteria. When a combination of criteria are applied the endpoint will return


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**Query**|**additionalQueryResultFields**  <br>*optional*|Currently only 'labels' is a valid value here. Use it to include a list of labels with each result.|< string > array(multi)||
|**Query**|**end**  <br>*optional*|Returns only workflows with an equal or earlier end datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.|string (date-time)||
|**Query**|**excludeLabelAnd**  <br>*optional*|Excludes workflows with the specified label.  If specified multiple times, excludes workflows with all of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"|< string > array(multi)||
|**Query**|**excludeLabelOr**  <br>*optional*|Excludes workflows with the specified label.  If specified multiple times, excludes workflows with any of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"|< string > array(multi)||
|**Query**|**id**  <br>*optional*|Returns only workflows with the specified workflow id.  If specified multiple times, returns workflows with any of the specified workflow ids.|< string > array(multi)||
|**Query**|**includeSubworkflows**  <br>*optional*|Include subworkflows in results. By default, it is taken as true.|boolean (boolean)||
|**Query**|**label**  <br>*optional*|Returns workflows with the specified label keys.  If specified multiple times, returns workflows with all of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"|< string > array(multi)||
|**Query**|**labelor**  <br>*optional*|Returns workflows with the specified label keys.  If specified multiple times, returns workflows with any of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"|< string > array(multi)||
|**Query**|**name**  <br>*optional*|Returns only workflows with the specified name.  If specified multiple times, returns workflows with any of the specified names.|< string > array(multi)||
|**Query**|**start**  <br>*optional*|Returns only workflows with an equal or later start datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.|string (date-time)||
|**Query**|**status**  <br>*optional*|Returns only workflows with the specified status.  If specified multiple times, returns workflows in any of the specified statuses.|< string > array(multi)||
|**Query**|**submission**  <br>*optional*|Returns only workflows with an equal or later submission time. Can be specified at most once. If both submission time and start date are specified, submission time should be before or equal to start date.|string (date-time)||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowQueryResponse](#workflowqueryresponse)|
|**403**|Workflow in terminal status|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="abort"></a>
## Abort a running workflow
```
POST /api/workflows/{version}/{id}/abort
```


#### Description
Request Cromwell to abort a running workflow. For instance this might be necessary in cases where you have submitted a workflow with incorrect inputs or no longer need the results. Cromwell will schedule a halt of all currently running jobs from this workflow.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**403**|Workflow in terminal status|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="labels"></a>
## Retrieves the current labels for a workflow
```
GET /api/workflows/{version}/{id}/labels
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[LabelsResponse](#labelsresponse)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="updatelabels"></a>
## Update labels for a workflow
```
PATCH /api/workflows/{version}/{id}/labels
```


#### Description
Update multiple labels for an existing workflow. When supplying a label with a key unique to the workflow submission, a new label key/value entry is appended to that workflow's metadata. When supplying a label with a key that is already associated to the workflow submission, the original label value is updated with the new value for that workflow's metadata.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|Workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**Body**|**labels**  <br>*required*|Custom labels submitted as JSON. Example: {"key-1":"value-1","key-2":"value-2"}|object||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[LabelsResponse](#labelsresponse)|
|**400**|Malformed Workflow ID|No Content|
|**403**|Workflow in terminal status|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="logs"></a>
## Get the logs for a workflow
```
GET /api/workflows/{version}/{id}/logs
```


#### Description
Returns paths to the standard out and standard error files that were generated during the execution of all calls in a workflow. A call has one or more standard out and standard error logs, depending on if the call was scattered or not. In the latter case, one log is provided for each instance of the call that has been run.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="metadata"></a>
## Get workflow and call-level metadata for a specified workflow
```
GET /api/workflows/{version}/{id}/metadata
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|
|**Query**|**excludeKey**  <br>*optional*|When specified, filters metadata to not return any field with a name which begins with this value. This key is used relative to the root of the response *and* relative to each call's metadata fields. Use 'calls' to filter out all call level metadata.|< string > array(multi)||
|**Query**|**expandSubWorkflows**  <br>*optional*|When true, metadata for sub workflows will be fetched and inserted automatically in the metadata response.|boolean||
|**Query**|**includeKey**  <br>*optional*|When specified, filters metadata to only return fields with names which begins with this value. This key is used relative to the root of the response *and* relative to each call's metadata fields.|< string > array(multi)||


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowMetadataResponse](#workflowmetadataresponse)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="outputs"></a>
## Get the outputs for a workflow
```
GET /api/workflows/{version}/{id}/outputs
```


#### Description
Retrieve the outputs for the specified workflow. Cromwell will return any outputs which currently exist even if a workflow has not successfully completed.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="releasehold"></a>
## Switch a workflow from 'On Hold' to 'Submitted' status
```
POST /api/workflows/{version}/{id}/releaseHold
```


#### Description
Request Cromwell to release the hold on a workflow. It will switch the status of a workflow from 'On Hold' to 'Submitted' so it can be picked for running. For instance this might be necessary in cases where you have submitted a workflow with workflowOnHold = true.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**403**|Workflow in terminal status|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="status"></a>
## Retrieves the current state for a workflow
```
GET /api/workflows/{version}/{id}/status
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="timing"></a>
## Get a visual diagram of a running workflow
```
GET /api/workflows/{version}/{id}/timing
```


#### Description
Returns a javascript file which will render a Gantt chart for the requested workflow. The bars in the chart represent start and end times for individual task invocations. This javascript is intended to be embedded into another web page.


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**id**  <br>*required*|A workflow ID|string||
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[WorkflowIdAndStatus](#workflowidandstatus)|
|**400**|Malformed Workflow ID|No Content|
|**404**|Workflow ID Not Found|No Content|
|**500**|Internal Error|No Content|


#### Tags

* Workflows


<a name="enginestatus"></a>
## Return the current health status of any monitored subsystems
```
GET /engine/{version}/status
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|All subsystems report an "ok" status|[StatusResponse](#statusresponse)|
|**500**|At least one subsystem does not have an "ok" status|[StatusResponse](#statusresponse)|


#### Tags

* Engine


<a name="engineversion"></a>
## Return the version of this Cromwell server
```
GET /engine/{version}/version
```


#### Parameters

|Type|Name|Description|Schema|Default|
|---|---|---|---|---|
|**Path**|**version**  <br>*required*|Cromwell API Version|string|`"v1"`|


#### Responses

|HTTP Code|Description|Schema|
|---|---|---|
|**200**|Successful Request|[VersionResponse](#versionresponse)|


#### Tags

* Engine


<a name="definitions"></a>
## Definitions

<a name="backendresponse"></a>
### BackendResponse

|Name|Description|Schema|
|---|---|---|
|**defaultBackend**  <br>*required*|The default backend of this server|string|
|**supportedBackends**  <br>*required*|The backends supported by this server|< string > array|


<a name="callmetadata"></a>
### CallMetadata
Call level metadata


|Name|Description|Schema|
|---|---|---|
|**backend**  <br>*optional*|The type of backend on which the call executed (e.g. JES, SGE, Local)|string|
|**backendLogs**  <br>*optional*|Paths to backend specific logs for this call|object|
|**backendStatus**  <br>*optional*|Status in backend-specific terms.  Currently this will only be defined for the JES backend.|string|
|**end**  <br>*optional*|End datetime of the call execution in ISO8601 format with milliseconds|string (date-time)|
|**executionStatus**  <br>*required*|Status in Cromwell execution terms.|string|
|**failures**  <br>*optional*||[FailureMessage](#failuremessage)|
|**inputs**  <br>*required*|Mapping of input fully qualified names to stringified values|object|
|**jobId**  <br>*optional*|Backend-specific job ID|string|
|**returnCode**  <br>*optional*|Call execution return code|integer|
|**start**  <br>*optional*|Start datetime of the call execution in ISO8601 format with milliseconds|string (date-time)|
|**stderr**  <br>*optional*|Path to the standard error file for this call|string|
|**stdout**  <br>*optional*|Path to the standard output file for this call|string|


<a name="descriptortype"></a>
### DescriptorType
One from a list of descriptor type strings (e.g. CWL, WDL). Note that these files can also include associated Docker/container files and test parameters that further describe a version of a tool

*Type* : enum (CWL, WDL)


<a name="descriptortypeandversion"></a>
### DescriptorTypeAndVersion
A workflow descriptor file type and version.


|Name|Description|Schema|
|---|---|---|
|**descriptorType**  <br>*required*||[DescriptorType](#descriptortype)|
|**descriptorTypeVersion**  <br>*required*|**Example** : `"1.0"`|string|


<a name="failuremessage"></a>
### FailureMessage
Failure messages


|Name|Description|Schema|
|---|---|---|
|**failure**  <br>*required*|The failure message|string|
|**timestamp**  <br>*required*|The time at which this failure occurred|string (date-time)|


<a name="labelsresponse"></a>
### LabelsResponse

|Name|Description|Schema|
|---|---|---|
|**id**  <br>*required*|The identifier of the workflow  <br>**Example** : `"label-key-1"`|string|
|**labels**  <br>*required*|The labels which have been updated  <br>**Example** : `"label-value-1"`|string|


<a name="mapvaluetype"></a>
### MapValueType
A type representing a map from one type to another.


|Name|Schema|
|---|---|
|**keyType**  <br>*required*|[ValueType](#valuetype)|
|**valueType**  <br>*required*|[ValueType](#valuetype)|


<a name="statusresponse"></a>
### StatusResponse
Returns the status of monitored subsystems.


|Name|Schema|
|---|---|
|**serviceName**  <br>*optional*|[serviceName](#statusresponse-servicename)|

<a name="statusresponse-servicename"></a>
**serviceName**

|Name|Schema|
|---|---|
|**messages**  <br>*optional*|< string > array|
|**ok**  <br>*optional*|boolean|


<a name="toolinputparameter"></a>
### ToolInputParameter
An input parameter for a tool or workflow.


|Name|Description|Schema|
|---|---|---|
|**default**  <br>*required*|The in-language expression used to evaluate a default value for this parameter, if none is supplied.|string|
|**name**  <br>*required*|The name of this input value (formatted as expected by the tool)|string|
|**optional**  <br>*required*|Whether the tool allows this value to not be specified|boolean|
|**typeDisplayName**  <br>*required*|An easy-to-read display name for the type of the input|string|
|**valueType**  <br>*required*||[ValueType](#valuetype)|


<a name="tooloutputparameter"></a>
### ToolOutputParameter
An output parameter for a tool or workflow.


|Name|Description|Schema|
|---|---|---|
|**name**  <br>*required*|The name of this input value (formatted as expected by the tool)|string|
|**typeDisplayName**  <br>*required*|An easy-to-read display name for the type of the output|string|
|**valueType**  <br>*required*||[ValueType](#valuetype)|


<a name="valuetype"></a>
### ValueType
The type expected for a given value.


|Name|Description|Schema|
|---|---|---|
|**arrayType**  <br>*optional*||[ValueType](#valuetype)|
|**mapType**  <br>*optional*||[MapValueType](#mapvaluetype)|
|**objectFieldTypes**  <br>*optional*||< [objectFieldTypes](#valuetype-objectfieldtypes) > array|
|**optionalType**  <br>*optional*||[ValueType](#valuetype)|
|**tupleTypes**  <br>*optional*||< [ValueType](#valuetype) > array|
|**typeName**  <br>*optional*|The type of this value|enum (String, File, Directory, Float, Int, Boolean, Optional, Array, Tuple, Map, Object, Pair)|

<a name="valuetype-objectfieldtypes"></a>
**objectFieldTypes**

|Name|Schema|
|---|---|
|**fieldName**  <br>*optional*|string|
|**fieldType**  <br>*optional*|[ValueType](#valuetype)|


<a name="versionresponse"></a>
### VersionResponse
Returns the version of Cromwell


|Name|Description|Schema|
|---|---|---|
|**cromwell**  <br>*optional*|The version of the Cromwell Engine  <br>**Example** : `"30"`|string|


<a name="workflowdescription"></a>
### WorkflowDescription

|Name|Description|Schema|
|---|---|---|
|**errors**  <br>*required*|The set of validation failure messages  <br>**Example** : `[ "The 'errors' field will be filled if 'valid' is false", "We might also provide warnings to a 'valid' workflow here", "Otherwise, 'errors' will be the empty array" ]`|< string > array|
|**inputs**  <br>*required*|A list of inputs for this tool  <br>**Example** : `[ {<br>  "name" : "my_wf.string_input",<br>  "valueType" : {<br>    "typeName" : "String"<br>  },<br>  "optional" : false,<br>  "default" : null,<br>  "typeDisplayName" : "String"<br>}, {<br>  "name" : "my_wf.array_input",<br>  "valueType" : {<br>    "typeName" : "Array",<br>    "arrayType" : {<br>      "typeName" : "String"<br>    }<br>  },<br>  "optional" : false,<br>  "default" : null,<br>  "typeDisplayName" : "Array[String]"<br>}, {<br>  "name" : "my_wf.optional_input",<br>  "valueType" : {<br>    "typeName" : "Optional",<br>    "optionalType" : {<br>      "typeName" : "String"<br>    }<br>  },<br>  "optional" : true,<br>  "default" : "hello",<br>  "typeDisplayName" : "String?"<br>}, {<br>  "name" : "my_wf.map_input",<br>  "valueType" : {<br>    "typeName" : "Map",<br>    "mapType" : {<br>      "keyType" : {<br>        "typeName" : "String"<br>      },<br>      "valueType" : {<br>        "typeName" : "Int"<br>      }<br>    }<br>  },<br>  "optional" : false,<br>  "default" : null,<br>  "typeDisplayName" : "Map[String, Int]"<br>}, {<br>  "name" : "my_wf.object_input",<br>  "valueType" : {<br>    "typeName" : "Object",<br>    "objectFieldTypes" : [ {<br>      "fieldName" : "int_field",<br>      "fieldType" : {<br>        "typeName" : "Int"<br>      }<br>    }, {<br>      "fieldName" : "int_array_field",<br>      "fieldType" : {<br>        "typeName" : "Array",<br>        "arrayType" : {<br>          "typeName" : "Int"<br>        }<br>      }<br>    } ]<br>  },<br>  "optional" : false,<br>  "default" : null,<br>  "typeDisplayName" : "Object"<br>}, {<br>  "name" : "my_wf.int_string_pair_input",<br>  "valueType" : {<br>    "typeName" : "Pair",<br>    "pairTypes" : [ {<br>      "leftType" : [ {<br>        "typeName" : "Int"<br>      } ]<br>    }, {<br>      "rightType" : [ {<br>        "typeName" : "String"<br>      } ]<br>    } ]<br>  },<br>  "optional" : false,<br>  "default" : null,<br>  "typeDisplayName" : "Pair[Int, String]"<br>} ]`|< [ToolInputParameter](#toolinputparameter) > array|
|**isRunnableWorkflow**  <br>*required*|Indicates whether this file can be run on its own (e.g. a WDL workflow)|boolean|
|**name**  <br>*required*|For a source file with one workflow and zero or more tasks, the name of the workflow. For a single task, the name of the task. For a source file with multiple tasks but no workflows, the empty string.|string|
|**outputs**  <br>*required*|A list of outputs for this tool  <br>**Example** : `[ {<br>  "name" : "my_wf.string_output",<br>  "valueType" : {<br>    "typeName" : "String"<br>  },<br>  "typeDisplayName" : "String"<br>}, {<br>  "name" : "my_wf.array_output",<br>  "valueType" : {<br>    "typeName" : "Array",<br>    "arrayType" : {<br>      "typeName" : "String"<br>    }<br>  },<br>  "typeDisplayName" : "Array[String]"<br>}, {<br>  "name" : "my_wf.map_output",<br>  "valueType" : {<br>    "typeName" : "Map",<br>    "mapType" : {<br>      "keyType" : {<br>        "typeName" : "String"<br>      },<br>      "valueType" : {<br>        "typeName" : "Int"<br>      }<br>    }<br>  },<br>  "typeDisplayName" : "Map[String, Int]"<br>}, {<br>  "name" : "my_wf.object_output",<br>  "valueType" : {<br>    "typeName" : "Object",<br>    "objectFieldTypes" : [ {<br>      "fieldName" : "int_field",<br>      "fieldType" : {<br>        "typeName" : "Int"<br>      }<br>    }, {<br>      "fieldName" : "int_array_field",<br>      "fieldType" : {<br>        "typeName" : "Array",<br>        "arrayType" : {<br>          "typeName" : "Int"<br>        }<br>      }<br>    } ]<br>  },<br>  "typeDisplayName" : "Object"<br>}, {<br>  "name" : "my_wf.int_string_pair_output",<br>  "valueType" : {<br>    "typeName" : "Pair",<br>    "tupleTypes" : [ {<br>      "typeName" : "Int"<br>    }, {<br>      "typeName" : "String"<br>    } ]<br>  },<br>  "typeDisplayName" : "Pair[Int, String]"<br>} ]`|< [ToolOutputParameter](#tooloutputparameter) > array|
|**submittedDescriptorType**  <br>*required*||[DescriptorTypeAndVersion](#descriptortypeandversion)|
|**valid**  <br>*required*|Indicates that the workflow is valid and that the inputs, if provided, are compatible with the workflow.|boolean|
|**validWorkflow**  <br>*required*|Indicates whether the workflow file is valid by itself. If inputs are provided, they are not considered when calculating this field; if inputs are not provided, the value is identical to `valid`.|boolean|


<a name="workflowidandstatus"></a>
### WorkflowIdAndStatus

|Name|Description|Schema|
|---|---|---|
|**id**  <br>*required*|The identifier of the workflow  <br>**Example** : `"00001111-2222-3333-aaaa-bbbbccccdddd"`|string|
|**status**  <br>*required*|The status of the workflow  <br>**Example** : `"Submitted"`|string|


<a name="workflowmetadataresponse"></a>
### WorkflowMetadataResponse
Workflow and call level metadata


|Name|Description|Schema|
|---|---|---|
|**calls**  <br>*optional*||[CallMetadata](#callmetadata)|
|**end**  <br>*optional*|End datetime of the workflow in ISO8601 format with milliseconds|string (date-time)|
|**failures**  <br>*optional*||[FailureMessage](#failuremessage)|
|**id**  <br>*required*|The identifier of the workflow|string|
|**inputs**  <br>*optional*|Map of input keys to input values|object|
|**outputs**  <br>*optional*|Map of output keys to output values|object|
|**start**  <br>*optional*|Start datetime of the workflow in ISO8601 format with milliseconds|string (date-time)|
|**status**  <br>*required*|The status of the workflow|string|
|**submission**  <br>*required*|Submission datetime of the workflow in ISO8601 format with milliseconds|string (date-time)|


<a name="workflowqueryparameter"></a>
### WorkflowQueryParameter
Workflow query parameters


|Name|Description|Schema|
|---|---|---|
|**end**  <br>*optional*|Returns only workflows with an equal or earlier end datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.|string (date-time)|
|**excludeLabelAnd**  <br>*optional*|Excludes workflows with the specified label.  If specified multiple times, excludes workflows with all of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"  <br>**Pattern** : `"^([a-z][-a-z0-9]*[a-z0-9])?[:]([a-z][-a-z0-9]*[a-z0-9])?$"`|string (array)|
|**excludeLabelOr**  <br>*optional*|Excludes workflows with the specified label.  If specified multiple times, excludes workflows with any of the specified label keys. Specify the label key and label value pair as separated with "label-key:label-value"  <br>**Pattern** : `"^([a-z][-a-z0-9]*[a-z0-9])?[:]([a-z][-a-z0-9]*[a-z0-9])?$"`|string (array)|
|**id**  <br>*optional*|Returns only workflows with the specified workflow id.  If specified multiple times, returns workflows with any of the specified workflow ids.  <br>**Pattern** : `"^[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}$"`|string|
|**includeSubworkflows**  <br>*optional*|Include subworkflows in results. By default, it is taken as true.|string (boolean)|
|**name**  <br>*optional*|Returns only workflows with the specified name.  If specified multiple times, returns workflows with any of the specified names.  <br>**Pattern** : `"^[a-zA-Z][a-zA-Z0-9_]*$"`|string|
|**page**  <br>*optional*|When pageSize is set, which page of results to return. If not set, the first page of 'pageSize' results will be returned.|integer|
|**pageSize**  <br>*optional*|The number of results to return at a time|integer|
|**start**  <br>*optional*|Returns only workflows with an equal or later start datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.|string (date-time)|
|**status**  <br>*optional*|Returns only workflows with the specified status.  If specified multiple times, returns workflows in any of the specified statuses.|enum (Submitted, Running, Aborting, Failed, Succeeded, Aborted)|
|**submission**  <br>*optional*|Returns only workflows with an equal or later submission time. Can be specified at most once. If both submission time and start date are specified, submission time should be before or equal to start date.|string (date-time)|


<a name="workflowqueryresponse"></a>
### WorkflowQueryResponse
Response to a workflow query


|Name|Schema|
|---|---|
|**results**  <br>*required*|< [WorkflowQueryResult](#workflowqueryresult) > array|
|**totalResultsCount**  <br>*required*|integer|


<a name="workflowqueryresult"></a>
### WorkflowQueryResult
Result for an individual workflow returned by a workflow query


|Name|Description|Schema|
|---|---|---|
|**end**  <br>*optional*|Workflow end datetime|string (date-time)|
|**id**  <br>*required*|Workflow ID|string|
|**name**  <br>*required*|Workflow name|string|
|**start**  <br>*optional*|Workflow start datetime|string (date-time)|
|**status**  <br>*required*|Workflow status|string|
|**submission**  <br>*optional*|Workflow submission datetime|string (date-time)|


<a name="workflowsubmitresponse"></a>
### WorkflowSubmitResponse

|Name|Description|Schema|
|---|---|---|
|**id**  <br>*required*|The identifier of the workflow  <br>**Example** : `"00001111-2222-3333-aaaa-bbbbccccdddd"`|string|
|**status**  <br>*required*|The status of the workflow  <br>**Example** : `"Submitted"`|string|

Cromwell provides a generic way to configure a backend relying on most High Performance Computing (HPC) frameworks, and with access to a shared filesystem.

The two main features that are needed for this backend to be used are a way to submit a job to the compute cluster and to get its status through the command line.
You can find example configurations for a variety of those backends here:

* [SGE](SGE)
* [LSF](LSF)
* [SLURM](SLURM)
* [HTCondor](HTcondor)

## FileSystems

### Shared FileSystem
HPC backends rely on being able to access and use a shared filesystem to store workflow results.

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.
Any input files to a call need to be localized into the `<call_dir>/inputs` directory. There are different localization strategies that Cromwell will try until one works:

* `hard-link` - This will create a hard link to the file
* `soft-link` - Create a symbolic link to the file. This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file
* `cached-copy` An experimental feature. This copies files to a file cache in 
`<workflow_root>/cached-inputs` and then hard links them in the `<call_dir>/inputs` directory. 

`cached-copy` is intended for a shared filesystem that runs on multiple physical disks, where docker containers are used. 
Hard-links don't work between different physical disks and soft-links don't work with docker. Copying uses a lot of
space if a multitude of tasks use the same input. `cached-copy` copies the file only once to the physical disk containing
the `<workflow_root>` and then uses hard links for every task that needs the input file. This can save a lot of space.

The default order in `reference.conf` is `hard-link`, `soft-link`, `copy`

Shared filesystem localization is defined in the `config` section of each backend. The default stanza for the Local and HPC backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```

### Additional FileSystems

HPC backends (as well as the Local backend) can be configured to be able to interact with other type of filesystems, where the input files can be located for example.
Currently the only other filesystem supported is Google Cloud Storage (GCS). See the [Google section](Google) of the documentation for information on how to configure GCS in Cromwell.
Once you have a google authentication configured, you can simply add a `gcs` stanza in your configuration file to enable GCS:

```
backend.providers.MyHPCBackend {
  filesystems {
    gcs {
      # A reference to a potentially different auth for manipulating files via engine functions.
      auth = "application-default"
    }
  }
}
```

### Exit code timeout

If the cluster forcefully kills a job, it is unable to write its exit code anymore.
To address this the option `exit-code-timeout-seconds` can be used.
Cromwell will check the aliveness of the job with the `check-alive` script, every `exit-code-timeout-seconds` (polling).
When a job is no longer alive and another `exit-code-timeout-seconds` seconds have passed without an RC file being made, Cromwell can mark the job as failed.
If retries are enabled the job is submitted again.
This option will enable polling with the `check-alive` option, this could cause high load on whatever system `check-alive` calls.

When the option `exit-code-timeout-seconds` is **not** set cromwell will only execute the `check-alive` option after a restart of a cromwell server.

```
backend {
  providers {
    <backend name> {
      config {
        exit-code-timeout-seconds = 120
        # other config options
      }
    }
  }
}
```
**Alibaba Cloud BCS Backend**

This backend adds support for execution jobs on Alibaba Cloud's BatchCompute service in a workflow.

### Configuring Backend

The backend is specified via the actor factory `BcsBackendLifecycleActorFactory`:

```hocon
backend {
  providers {
    BCS {
      config {
        actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
        # ... other configuration
      }
    }
  }
}
```

You'll likely also want to change the default backend to your new backend, by setting this configuration value:

```hocon
backend {
  providers {
    default = BCS
  }
}
```

Before reading further in this section please see the [Getting started on Alibaba Cloud](../tutorials/BCSIntro.md) for instructions on configuring to Alibaba Cloud services.

The configuration file for Alibaba Cloud will look like the following.

```hocon
backend {
  providers {
    BCS {
      config {
        actor-factory = "cromwell.backend.impl.bcs.BcsBackendLifecycleActorFactory"
        root = "oss://<test-bucket>/cromwell-dir"
        region = "<test-region>"
        access-id = "<test-access-id>"
        access-key = "<test-access-key>"
       
        filesystems {
        # ... to be filled in
        }
        
        default-runtime-attributes {
        # ... to be filled in
        }
      }
    }
  }
}
```

- `<test-bucket>` : OSS bucket name.
- `<test-region>` : Region in Alibaba Cloud chosen to deploy cromwell, it must be the same as the region of `<test-bucket>`.
- `<test-access-id>` : Access ID to access Alibaba Cloud services through restful API.
- `<test-access-key>` : Access key to access Alibaba Cloud services through restful API.

The values above are necessary for Cromwell to submit and poll status of workflow jobs to and from Alibaba Cloud BatchCompute service.
The `filesystems` stanza in the backend config defines how to configure a filesystem in Alibaba Cloud. Details of filesystem related configurations will be explained in the next section. 

### File Systems

Currently, this backend only works with objects on an Alibaba Cloud OSS filesystem. It's necessary to supply all values in 
the configuration key `backend.providers.BCS.config.filesystems.auth` in order to read/write OSS file system objects in Alibaba Backend jobs. A typical config looks like this:

- `<test-oss-endpoint>` - API endpoint to access OSS bucket `<test-bucket>`.
- `<test-access-id>` - Access ID to access Alibaba Cloud services through restful API. 
- `<test-access-key>` - Access key to access Alibaba Cloud services through restful API. 
- `<refresh-interval>` - The interval of auth refreshing if you are using an STS(Alibaba Cloud Security Token Service) way to access the OSS filesystem.

```hocon
backend {
  providers {
    BCS {
      config {
        # BCS related configurations mentioned above
       
        filesystems {
          oss {
            auth {
              endpoint = "<test-oss-endpoint>"
              access-id = "<test-access-id>"
              access-key = "<test-access-key>"
            }
            refresh-interval = 1800
          }
        }
        
        default-runtime-attributes {
        # ... to be filled in
        }
      }
    }
  }
}
```

### Runtime Attributes

This backend supports additional runtime attributes that are specified in the configuration key `backend.providers.BCS.config.runtime-attributes`. 
It uses the same syntax as specifying runtime attributes in a task in WDL. A typical runtime attributes example for BCS backend looks like this:

```hocon
backend {
  providers {
    BCS {
      config {
        # BCS and OSS related configurations mentioned above
       
        default-runtime-attributes {
          cluster: "OnDemand ecs.sn1ne.large img-ubuntu-vpc"
          imageId: "img-ubuntu-vpc"
          mounts: "oss://<test-bucket>/inputs/ /home/inputs/ false"
          dockerTag: "ubuntu/latest oss://<test-bucket>/registry/ubuntu/"
          docker: "registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2"
          userData: "key value"
          reserveOnFail: true
          autoReleaseJob: true
          verbose: false
          systemDisk: "cloud 50"
          dataDisk: "cloud 250 /home/data/"
          timeout: 3000
          isv: "abc"
        }
      }
    }
  }
}
```

#### cluster

There are two different ways of specifying an Alibaba Cloud BatchCompute cluster in which workflow jobs run.

- Reserved cluster - A pre-created cluster ID in BatchCompute service like this:

```hocon
      default-runtime-attributes {
        cluster: "cls-your-cluster-id"
      }
```

- Auto cluster - Cluster configuration to create a new runtime cluster bound to the workflow job:

  - `<resource-type>` - Type of resource, can only support `OnDemand` and `Spot` currently.
  - `<instance-type>` - Type of VM instance. Go to <a href="https://help.aliyun.com/document_detail/25378.html" target="_blank">Alibaba Cloud BatchCompute Instance Type</a> to choose a suitable type for you.
  - `<image-id>` - Image ID of Alibaba Cloud BatchCompute service to create a VM.

```hocon
      default-runtime-attributes {
        cluster: "<resource-type> <instance-type> <image-id>"
        # Maybe like cluster: "OnDemand ecs.sn1ne.large img-ubuntu"
      }
```

#### imageId

The BCS job image ID can be specified by the runtime `cluster`, while if you are going to use Call Caching, another optional runtime 
`imageId` should be specified. The change of image ID will lead to a cache miss and the call will be executed as normal.

```hocon
      default-runtime-attributes {
        imageId: "img-ubuntu-vpc"
      }
```

#### mounts

BCS jobs can mount both OSS and [Alibaba Cloud NAS](https://www.aliyun.com/product/nas) to local filesystem as a file or a directory in VM.
It uses distribute-caching and lazy-load techniques to optimize concurrently read requests of the OSS file system. 
You can mount your OSS objects to VM like this:

- `<mount-src>` - An OSS object path or OSS prefix or NAS address to mount from, such as
  `oss://<test-bucket>/inputs/ /home/inputs/ false` for OSS
  and `nas://0266e49fea-yio75.cn-beijing.nas.aliyuncs.com:/ /home/nas/ true` for NAS. See the [NAS mount](https://www.alibabacloud.com/help/doc-detail/50494.htm) for more details of NAS mount.
- `<mount-destination>` - An unix file path or directory path to mount to in VM.
- `<write-support>` - Writable for mount destination, only works for directory.

```hocon
default-runtime-attributes {
  mounts: "<mount-src> <mount-destination> <write-support>"
}
```



#### docker

This backend supports docker images pulled from OSS registry or Alibaba Cloud Container Registry.

##### OSS registry
```hocon
default-runtime-attributes {
  dockerTag: "<docker-image> <oss-registry-path>"
}
```

- `<docker-image>` - Docker image name such as: ubuntu:latest.
- `<oss-registry-path>` - Image path in OSS filesyetem where you pushed your docker image.


##### Alibaba Cloud Container Registry

```hocon
default-runtime-attributes {
  docker: "<docker-image-with-tag>"
}
```
- `docker-image-with-tag` - Docker image stored in Alibaba Cloud Container Registry, such as `registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2`.
#### userData

If a runtime cluster is specified, it's possible to pass some environment variables to VM when running BCS jobs.
It looks like this:

```hocon
 default-runtime-attributes {
   userData: "key1 value1, key2, value2"
 }
```

#### autoReleaseJob

The Alibaba Cloud BatchCompute service limits the number of simultaneous jobs per user. Jobs created by the backend are
deleted when the job finishes. However it is possible to tell the BCS backend to not delete the job when the call
finishes by setting `autoReleaseJob` to `false`:

```hocon
 default-runtime-attributes {
   autoReleaseJob: false
 }
```

#### systemDisk

If it's necessary to run a job with a particular system disk type or disk size, a runtime attribute named `systemDisk` can be used to
specify disk type and size.

- `<disk-type>` - Disk type to be used, can only support `cloud` or `cloud_efficiency` currently.
- `<disk-size-in-GB>` - Disk size to be used.

```hocon
 default-runtime-attributes {
   systemDisk: "<disk-type> <disk-size-in-GB>"
 }
```

#### dataDisk

The system disk size can support up to 500GB. One can mount another data disk in VM if needed.

- `<disk-type>` - Disk type to be used, can only support `cloud` or `cloud_efficiency` currently.
- `<disk-size-in-GB>` - Disk size to be used.
- `<mount-point>` - Destination the data disk mounted to in VM.

```hocon
 default-runtime-attributes {
   dataDisk: "<disk-type> <disk-size-in-GB> <mount-point>"
 }
```
#### isv

If a BCS ISV is used in the task, a runtime attribute name `isv` can be used to specify it.

```hocon
 default-runtime-attributes {
   isv: "abc"
 }
```


###CallCaching
BCS supports CallCaching feature when the docker image is from Alibaba Cloud Container Registry.
The configuration file will look like the following:
```hocon
call-caching {
  enabled = true
  invalidate-bad-cache-results = true

}

docker {
  hash-lookup {
    enabled = true
    method = "remote"
    alibabacloudcr {
      num-threads = 5
      auth {
        access-id = xxxx
        access-key = yyyy
        security-token = zzzz
      }
    }
  }
}

backend {
  providers {
    BCS {
      config {
        # BCS and OSS related configurations mentioned above
        filesystems {
          oss {
            caching {
               duplication-strategy = "reference"
               invalidate-bad-cache-results = true
            }
            # ... to be filled in
          }
        }
        default-runtime-attributes {
          docker: "registry.cn-shanghai.aliyuncs.com/batchcompute/myubuntu:0.2"
          # ... to be filled in
        }
      }
    }
  }
}
```

- `docker.hash-lookup.method` - BCS only supports `remote` method for hash-lookup
- `filesystems.oss.caching.duplication-strategy` - BCS only supports `reference` for duplication strategy.**TES Backend**

The TES backend submits jobs to a server that complies with the protocol described by the [GA4GH schema](https://github.com/ga4gh/task-execution-schemas).

This backend creates three files in the `<call_dir>`:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process

The `script` file contains:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

`<container_call_root>` would be equal to the runtime attribute `dockerWorkingDir`  or `/cromwell-executions/<workflow_uuid>/call-<call_name>/execution` if this attribute is not supplied.

**Configuring**

Configuring the TES backend is straightforward; one must only provide the TES API endpoint for the service. 

```hocon
backend {
  default = "TES"
  providers {
    TES {
      actor-factory = "cromwell.backend.impl.tes.TesBackendLifecycleActorFactory"
      config {
        endpoint = "https://<some-url>/v1/tasks"
        root = "cromwell-executions"
        dockerRoot = "/cromwell-executions"
        concurrent-job-limit = 1000
      }
    }
  }
}
```

**Supported File Systems**  

Currently this backend only works with files on a Local or Shared File System. 

**Docker**

This backend supports the following optional [Runtime Attributes](../RuntimeAttributes) and [Workflow Options](../wf_options/Overview/) for working with Docker:

* `docker`: Docker image to use such as "Ubuntu".
* `dockerWorkingDir`: defines the working directory in the container.

**CPU, Memory and Disk** 

This backend supports CPU, memory and disk size configuration through the use of the following [Runtime Attributes](../RuntimeAttributes) and [Workflow Options](../wf_options/Overview/):  

* `cpu` defines the amount of CPU to use. 
    * Type: Integer (ex: 4)
* `memory` defines the amount of memory to use. 
    * Type: String (ex: "4 GB" or "4096 MB")
* `disk` defines the amount of disk to use. 
    * Type: String (ex: "1 GB" or "1024 MB")
* `preemptible` defines whether or not to use preemptible VMs. 
    * Type: Boolean (ex: "true" or "false")

If they are not set, the TES backend may use default values.

**TESK**

[TESK](https://github.com/EMBL-EBI-TSI/TESK) is an implementation of the TES interface that uses Kubernetes and FTP.
When running Cromwell with a TESK backend, you will want to customize the way Cromwell process globs, as kubernetes will not work well with hard links in a lot of cases which is the default behavior in Cromwell.
By adding this to the `config` section of the TES backend in Cromwell, Cromwell will use symlinks instead.  

`glob-link-command = "ls -L GLOB_PATTERN 2> /dev/null | xargs -I ? ln -s ? GLOB_DIRECTORY"`
**HTCondor Backend**

Allows to execute jobs using HTCondor which is a specialized workload management system for compute-intensive jobs created by the Center for High Throughput Computing in the Department of Computer Sciences at the University of Wisconsin-Madison (UW-Madison).

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    HtCondor {
      config {
        actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
        # ... other configuration
      }
    }
  }
}
```

This backend makes the same assumption about the filesystem that the local backend does: the Cromwell process and the jobs both have read/write access to the CWD of the job.

The CWD will contain a `script.sh` file which will contain the same contents as the Local backend:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

The job is launched with a configurable script command such as:

```
chmod 755 ${script}
cat > ${cwd}/execution/submitFile <<EOF
Iwd=${cwd}/execution
requirements=${nativeSpecs}
leave_in_queue=true
request_memory=${memory_mb}
request_disk=${disk_kb}
error=${err}
output=${out}
log_xml=true
request_cpus=${cpu}
executable=${script}
log=${cwd}/execution/execution.log
queue
EOF
condor_submit ${cwd}/execution/submitFile
```

The HtCondor backend gets the job ID from parsing the `submit.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

The command used to submit the job is specified under the configuration key `backend.providers.HtCondor.config.submit`. It uses the same syntax as a command in WDL, and will be provided the variables:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `cwd` - The path where the script should be run.
* `out` - The path to the stdout.
* `err` - The path to the stderr.
* `job_name` - A unique name for the job.

This backend also supports docker as optional feature. Configuration key `backend.providers.HtCondor.config.submit-docker` is specified for this end. When the WDL contains a docker runtime attribute, this command will be provided with two additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.
* `docker_script` - The path of the `script` within the docker container.
* `docker_out` - The path of the `out` within the docker container.
* `docker_err` - The path of the `err` within the docker container.

```
chmod 755 ${script}
cat > ${cwd}/execution/dockerScript <<EOF
#!/bin/bash
docker run --rm -i -v ${cwd}:${docker_cwd} ${docker} /bin/bash ${docker_script}
EOF
chmod 755 ${cwd}/execution/dockerScript
cat > ${cwd}/execution/submitFile <<EOF
Iwd=${cwd}/execution
requirements=${nativeSpecs}
leave_in_queue=true
request_memory=${memory_mb}
request_disk=${disk_kb}
error=${cwd}/execution/stderr
output=${cwd}/execution/stdout
log_xml=true
request_cpus=${cpu}
executable=${cwd}/execution/dockerScript
log=${cwd}/execution/execution.log
queue
EOF
condor_submit ${cwd}/execution/submitFile
```

This backend support additional runtime attributes that are specified in the configuration key `backend.providers.HtCondor.config.runtime-attributes`. It uses the same syntax as specifying runtime attributes in a task in WDL.

There are five special runtime attribute configurations, `cpu`, `memory_mb`, `disk_kb`, `nativeSpecs`, `docker`.
Optional values are defined with the prefix `?` attached to the type.

```
backend {
  providers {
    HtCondor {
      config {
        # ... other configuration
	    runtime-attributes = """
	       Int cpu = 1
	       Float memory_mb = 512.0
	       Float disk_kb = 256000.0
	       String? nativeSpecs
	       String? docker
	    """
      }
    }
  }
}
```

**Native Specifications**

The use of runtime attribute 'nativeSpecs' allows to the user to attach custom HtCondor configuration to tasks.
An example of this is when there is a need to work with 'requirements' or 'rank' configuration.

```
"runtimeAttributes": {
    cpu = 2
    memory = "1GB"
    disk = "1GB"
    nativeSpecs: "TARGET.Arch == \"INTEL\" && TARGET.Memory >= 64"
}
```

`nativeSpecs` attribute needs to be specified as String.

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
*This sample configuration is a community contribution and therefore not officially supported.*

The following configuration can be used as a base to allow Cromwell to interact with a [volcano](https://www.github.com/volcano-sh/volcano) cluster and dispatch jobs to it:

```hocon
Volcano {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    runtime-attributes = """
    Int runtime_minutes = 600
    Int cpus = 2
    Int requested_memory_mb_per_core = 8000
    String queue = "short"
    """

    submit = """
        vcctl job run -f ${script}
    """
    kill = "vcctl job delete -N ${job_id}"
    check-alive = "vcctl job view -N ${job_id}"
    job-id-regex = "(\\d+)"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)

If you have any questions about Volcano, please open an issue at
https://www.github.com/volcano-sh/volcano/issues
or contact us at
Slack Channel : https://volcano-sh.slack.com
Mailing List : https://groups.google.com/forum/#!forum/volcano-sh
otherwise, feel free to mail to : klaus1982.cn@gmail.com
The following configuration can be used as a base to allow Cromwell to interact with an [LSF](https://en.wikipedia.org/wiki/Platform_LSF) cluster and dispatch jobs to it:

```hocon
LSF {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    submit = "bsub -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} /usr/bin/env bash ${script}"
    kill = "bkill ${job_id}"
    check-alive = "bjobs ${job_id}"
    job-id-regex = "Job <(\\d+)>.*"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
**Google Cloud Backend**

Google Genomics Pipelines API is a Docker-as-a-service from Google. It was formerly called JES (Job Execution Service);
you may see outdated references to the older JES terminology in Cromwell configuration files and code.

This section offers detailed configuration instructions for using Cromwell with the Pipelines API in all supported
authentication modes. Before reading futher in this section please see the
[Getting started on Google Pipelines API](../tutorials/PipelinesApi101) for instructions common to all authentication modes
and detailed instructions for the application default authentication scheme in particular.
The instructions below assume you have created a Google Cloud Storage bucket and a Google project enabled for the appropriate APIs.

**Configuring Authentication**

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are four different
authentication schemes that might be used:

* `application_default` (default, recommended) - Use [application default](https://developers.google.com/identity/protocols/application-default-credentials) credentials.
* `service_account` - Use a specific service account and key file (in PEM format) to authenticate.
* `user_account` - Authenticate as a user.
* `user_service_account` - Authenticate each individual workflow using service account credentials supplied in the workflow options.

The `auths` block in the `google` stanza defines the authentication schemes within a Cromwell deployment:

```hocon
google {
  application-name = "cromwell"
  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    },
    {
      name = "service-account"
      scheme = "service_account"
      service-account-id = "my-service-account"
      pem-file = "/path/to/file.pem"
    },
    {
      name = "user-service-account"
      scheme = "user_service_account"
    }
  ]
}
```

These authentication schemes can be referenced by name within other portions of the configuration file.  For example, both
the `genomics` and `filesystems.gcs` sections within a Google configuration block must reference an auth defined in this block.
The auth for the `genomics` section governs the interactions with Google itself, while `filesystems.gcs` governs the localization
of data into and out of GCE VMs.

**Application Default Credentials**

By default, application default credentials will be used.  Only `name` and `scheme` are required for application default credentials.

To authenticate, run the following commands from your command line (requires [gcloud](https://cloud.google.com/sdk/gcloud/)):

```
$ gcloud auth login
$ gcloud config set project my-project
```

**Service Account**

First create a new service account through the [API Credentials](https://console.developers.google.com/apis/credentials) page.  Go to **Create credentials -> Service account key**.  Then in the **Service account** dropdown select **New service account**.  Fill in a name (e.g. `my-account`), and select key type of JSON.

Creating the account will cause the JSON file to be downloaded.  The structure of this file is roughly like this (account name is `my-account`):

```
{
  "type": "service_account",
  "project_id": "my-project",
  "private_key_id": "OMITTED",
  "private_key": "-----BEGIN PRIVATE KEY-----\nBASE64 ENCODED KEY WITH \n TO REPRESENT NEWLINES\n-----END PRIVATE KEY-----\n",
  "client_email": "my-account@my-project.iam.gserviceaccount.com",
  "client_id": "22377410244549202395",
  "auth_uri": "https://accounts.google.com/o/oauth2/auth",
  "token_uri": "https://accounts.google.com/o/oauth2/token",
  "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
  "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/my-account%40my-project.iam.gserviceaccount.com"
}
```

Most importantly, the value of the `client_email` field should go into the `service-account-id` field in the configuration (see below).  The
`private_key` portion needs to be pulled into its own file (e.g. `my-key.pem`).  The `\n`s in the string need to be converted to newline characters.

While technically not part of Service Account authentication mode, one can also override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`.  The service account you provide must have been granted Service Account Actor role to Cromwell's primary service account. As this only affects Google Pipelines API and not GCS, it's important that this service account, and the service account specified in `JES.config.genomics.auth` can both read/write the location specified by `JES.config.root`

**User Service Account**

A [JSON key file for the service account](../wf_options/Google.md) must be passed in via the `user_service_account_json` field in the [Workflow Options](../wf_options/Google.md) when submitting the job. Omitting this field will cause the workflow to fail. The JSON should be passed as a string and will need to have no newlines and all instances of `"` and `\n` escaped. 

In the likely event that this service account does not have access to Cromwell's default google project the `google_project` workflow option must be set. In the similarly likely case that this service account can not access Cromwell's default google bucket, the `jes_gcs_root` workflow option should be set appropriately.

For information on the interaction of `user_service_account_json` with private Docker images please see the `Docker` section below.  

**Docker**

It's possible to reference private Docker images to which only particular Docker Hub accounts have access:

```
task mytask {
  command {
    ...
  }
  runtime {
    docker: "private_repo/image"
    memory: "8 GB"
    cpu: "1"
  }
  ...
}
```

In order for a private image to be used the appropriate Docker configuration must be provided. If the Docker images being used
are public there is no need to add this configuration.

For Pipelines API (PAPI) version 1:
```
backend {
  default = "PAPIv1"
  providers {
    PAPIv1 {
      actor-factory = "cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory"
      config {
        dockerhub {
          token = "base64-encoded-docker-hub-username:password"
        }
      }
    }
  }
}
```

`token` is the standard base64-encoded username:password for the appropriate Docker Hub account.

For PAPI version 2 alpha 1:

```
backend {
  default = "PAPIv2"
  providers {
    PAPIv2 {
      actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
      config {
        dockerhub {
          token = "base64-encoded-docker-hub-username:password"
          key-name = "name/of/the/kms/key/used/for/encrypting/and/decrypting/the/docker/hub/token"
          auth = "reference-to-the-auth-cromwell-should-use-for-kms-encryption"
        }
      }
    }
  }
}
```

For PAPI version 2 beta:

```
backend {
  default = "PAPIv2"
  providers {
    PAPIv2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        dockerhub {
          token = "base64-encoded-docker-hub-username:password"
          key-name = "name/of/the/kms/key/used/for/encrypting/and/decrypting/the/docker/hub/token"
          auth = "reference-to-the-auth-cromwell-should-use-for-kms-encryption"
        }
      }
    }
  }
}
```

`key-name` is the name of the Google KMS key Cromwell should use for encrypting the Docker `token` before including it
in the PAPI job execution request. This `key-name` will also be included in the PAPI job execution
request and will be used by PAPI to decrypt the Docker token used by `docker login` to enable access to the private Docker image.
 
`auth` is a reference to the name of an authorization in the `auths` block of Cromwell's `google` config.
Cromwell will use this authorization for encrypting the Google KMS key.

The equivalents of `key-name`, `token` and `auth` can also be specified in workflow options which take
precedence over values specified in configuration. The corresponding workflow options are named `docker_credentials_key_name`,
`docker_credentials_token`, and `user_service_account_json`. While the config value `auth` refers to an auth defined in the 
`google.auths` stanza elsewhere in Cromwell's
configuration, `user_service_account_json` is expected to be a literal escaped Google service account auth JSON.
See the `User Service Account` section above for more information on using user service accounts.
If the key, token or auth value is provided in workflow options then the corresponding private Docker configuration value
is not required, and vice versa. Also note that for the `user_service_account_json` workflow option to work an auth of type `user_service_account`
must be defined in Cromwell's `google.auths` stanza; more details in the `User Service Account` section above.

Example PAPI v2 workflow options for private Docker configuration:

```
{
  "docker_credentials_key_name": "name/of/the/kms/key/used/for/encrypting/and/decrypting/the/docker/hub/token",
  "docker_credentials_token": "base64_username:password",
  "user_service_account_json": "<properly escaped user service account JSON file>"
}
```

Important

If any of the three private Docker configuration values of key name, auth, or Docker token are missing, PAPI v2 will not perform a `docker login`.
If the Docker image to be pulled is not public the `docker pull` will fail which will cause the overall job to fail.

If using any of these private Docker workflow options it is advisable to add
them to the `workflow-options.encrypted-fields` list in Cromwell configuration.


**Monitoring**

In order to monitor metrics (CPU, Memory, Disk usage...) about the VM during Call Runtime, a workflow option can be used to specify the path to a script that will run in the background and write its output to a log file.

```
{
  "monitoring_script": "gs://cromwell/monitoring/script.sh"
}
```

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.  This feature is meant to run a script in the background during long-running processes.  It's possible that if the task is very short that the log file does not flush before de-localization happens and you will end up with a zero byte file.

**Google Cloud Storage Filesystem**

On the Google Pipelines backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the Google backend the `jes_gcs_root` [Workflow Option](../wf_options/Google) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the Google backend.

**Pipeline timeout**

Google sets a default pipeline timeout of 7 days, after which the pipeline will abort. Setting `pipeline-timeout` overrides this limit to a maximum of 30 days.

```hocon
backend.providers.PAPIv2.config {
    pipeline-timeout: 14 days
}
```

**Enabling FUSE capabilities**

*This is a community contribution and not officially supported by the Cromwell team.*
By default Cromwell task containers doesn't allow to mount any FUSE filesystems. It happens because containers are launched without specific linux capabilities being enabled. 
Google pipelines backend supports running containers with the enabled capabilities and so does Cromwell. 

If you need to use fuses within task containers then you can set `enable_fuse` workflow option. 

```
{
    "enable_fuse": true
}
```

Differently you can enable support for fuses right in your backend configuration.

```
backend.providers.Papiv2.config {
    genomics {
        enable-fuse = true
    }
}
```

There is a list of limitations regarding the usage of FUSE filesystems:

+ Any inputs brought in via a FUSE filesystem will not be considered for call caching.
+ Any outputs stored via a FUSE filesystem will not be recreated if a task is replayed from a call-cache hit.
+ If the filesystem is writable, your job is potentially no longer idempotent - Cromwell may decide to retry your job for you, and you might get unforeseen file collisions or even incorrect results if that happens.

#### Google Labels

Every call run on the Pipelines API backend is given certain labels by default, so that Google resources can be queried by these labels later. 
The current default label set automatically applied is:

| Key | Value | Example | Notes |
|-----|-------|---------|-------|
| cromwell-workflow-id | The Cromwell ID given to the root workflow (i.e. the ID returned by Cromwell on submission) | cromwell-d4b412c5-bf3d-4169-91b0-1b635ce47a26 | To fit the required [format](#label-format), we prefix with 'cromwell-' |
| cromwell-sub-workflow-name | The name of this job's sub-workflow | my-sub-workflow | Only present if the task is called in a subworkflow. |
| wdl-task-name | The name of the WDL task | my-task | |
| wdl-call-alias | The alias of the WDL call that created this job | my-task-1 | Only present if the task was called with an alias. |

Any custom labels provided as '`google_labels`' in the [workflow options](../wf_options/Google) are also applied to Google resources by the Pipelines API.

## Using NCBI Sequence Read Archive (SRA) Data

The v2alpha1 and v2beta backends support accessing [NCBI
SRA](https://www.ncbi.nlm.nih.gov/sra) accessions natively.  To configure this
support you'll need to enable it in your config file like so:

```hocon
filesystems {
  sra {
    class = "cromwell.filesystems.sra.SraPathBuilderFactory"
    docker-image = "fusera/fusera:alpine"
    ngc = "bmNiaV9nYXAfiwgAAAAAAAADBcHBDYQgEADAv1XQAGYXcfErUe5x0diCiFESA0Y8/VD8zTzrlXwMDEsoII9usPT5znZSmTqUohaSg5Gay14TbxsluMGOSBuqDEKefvbwCzv3BAAKoexb5uIbjjg7dq/p9mH7A5VTImxjAAAA"
  }
}
```

This filesystem has two required configuration options:
* `docker-image`: The [fusera](https://github.com/mitre/fusera) docker image to
  use to provide access.  This can be a custom image, but using the public
  [fusera/fusera:alpine](https://hub.docker.com/r/fusera/fusera/) image is
  recommended.
* `ngc`: A base-64 encoded NGC file.  This is provided through the NCBI
  interface.  Please see [the
  documentation](https://www.ncbi.nlm.nih.gov/books/NBK63512/#Download.are_downloaded_files_encrypted)
  for more information on obtaining your NGC.  The `ngc` value provided above
  is the sample credential file.

### Virtual Private Network

Cromwell can arrange for jobs to run in specific GCP private networks via the `config.virtual-private-cloud` stanza of a PAPI v2 backend.
There are two ways of specifying private networks:

* [Literal network and subnetwork values](#virtual-private-network-via-literals) that will apply to all projects
* [Google project labels](#virtual-private-network-via-labels) whose values in a particular Google project will specify the network and subnetwork

#### Virtual Private Network via Literals

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        virtual-private-cloud {
          network-name = "vpc-network"
          subnetwork-name = "vpc-subnetwork"
        }
        ...
      }
    }
  }
}
```

The `network-name` and `subnetwork-name` should reference the name of your private network and subnetwork within that
network respectively. The `subnetwork-name` is an optional config.

For example, if your `virtual-private-cloud` config looks like the one above, then Cromwell will use the value of the
configuration key, which is `vpc-network` here, as the name of private network and run the jobs on this network.
If the network name is not present in the config Cromwell will fall back to trying to run jobs on the default network.

If the `network-name` or `subnetwork-name` values contain the string `${projectId}` then that value will be replaced
by Cromwell with the name of the project running the Pipelines API.

If the `network-name` does not contain a `/` then it will be prefixed with `projects/${projectId}/global/networks/`.

Cromwell will then pass the network and subnetwork values to the Pipelines API. See the documentation for the
[Cloud Life Sciences API](https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.pipelines/run#Network)
for more information on the various formats accepted for `network` and `subnetwork`.

#### Virtual Private Network via Labels

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        virtual-private-cloud {
          network-label-key = "my-private-network"
          subnetwork-label-key = "my-private-subnetwork"
          auth = "reference-to-auth-scheme"
        }
        ...
      }
    }
  }
}
```


The `network-label-key` and `subnetwork-label-key` should reference the keys in your project's labels whose value is the name of your private network
and subnetwork within that network respectively. `auth` should reference an auth scheme in the `google` stanza which will be used to get the project metadata from Google Cloud.
The `subnetwork-label-key` is an optional config.

For example, if your `virtual-private-cloud` config looks like the one above, and one of the labels in your project is

```
"my-private-network" = "vpc-network"
```

Cromwell will get labels from the project's metadata and look for a label whose key is `my-private-network`.
Then it will use the value of the label, which is `vpc-network` here, as the name of private network and run the jobs on this network.
If the network key is not present in the project's metadata Cromwell will fall back to trying to run jobs using literal
network labels, and then fall back to running on the default network.

### Custom Google Cloud SDK container
Cromwell can't use Google's container registry if VPC Perimeter is used in project.
Own repository can be used by adding `cloud-sdk-image-url` reference to used container:

```
google {
  ...
  cloud-sdk-image-url = "eu.gcr.io/your-project-id/cloudsdktool/cloud-sdk:354.0.0-alpine"
  cloud-sdk-image-size-gb = 1
}
```

### Parallel Composite Uploads

Cromwell can be configured to use GCS parallel composite uploads which can greatly improve delocalization performance. This feature
is turned off by default but can be enabled backend-wide by specifying a `gsutil`-compatible memory specification for the key
`genomics.parallel-composite-upload-threshold` in backend configuration. This memory value represents the minimum size an output file
must have to be a candidate for `gsutil` parallel composite uploading:

```
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        genomics {
          ...
          parallel-composite-upload-threshold = 150M
          ...
        }
        ...
      }
    }
  }
}
```

Alternatively this threshold can be specified in workflow options using the key `parallel-composite-upload-threshold`,
which takes precedence over a setting in configuration. The default setting for this threshold is `0` which turns off
parallel composite uploads; a value of `0` can also be used in workflow options to turn off parallel composite uploads
in a Cromwell deployment where they are turned on in config.

#### Issues with composite files

Please see the [Google documentation](https://cloud.google.com/storage/docs/gsutil/commands/cp#parallel-composite-uploads)
describing the benefits and drawbacks of parallel composite uploads.

The actual error message observed when attempting to download a composite file on a system without a compiled `crcmod`
looks like the following:

```
/ # gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp gs://my-bucket/composite.bam .
Copying gs://my-bucket/composite.bam...
==> NOTE: You are downloading one or more large file(s), which would
run significantly faster if you enabled sliced object downloads. This
feature is enabled by default but requires that compiled crcmod be
installed (see "gsutil help crcmod").

CommandException:
Downloading this composite object requires integrity checking with CRC32c,
but your crcmod installation isn't using the module's C extension, so the
hash computation will likely throttle download performance. For help
installing the extension, please see "gsutil help crcmod".

To download regardless of crcmod performance or to skip slow integrity
checks, see the "check_hashes" option in your boto config file.

NOTE: It is strongly recommended that you not disable integrity checks. Doing so
could allow data corruption to go undetected during uploading/downloading.
/ #
```

As the message states, the best option would be to have a compiled `crcmod` installed on the system.
Turning off integrity checks on downloads does get around this issue but really isn't a great idea.

#### Parallel composite uploads and call caching

Because the parallel composite upload threshold is not considered part of the hash used for call caching purposes, calls
which would be expected to generate non-composite outputs may call cache to results that did generate composite
outputs. Calls which are executed and not cached will always honor the parallel composite upload setting at the time of
their execution.

### Migration from Google Cloud Genomics v2alpha1 to Google Cloud Life Sciences v2beta

1. If you currently run your workflows using Cloud Genomics v2alpha1 and would like to switch to Google Cloud Life 
Sciences v2beta, you will need to do a few changes to your configuration file: `actor-factory` value should be changed 
from `cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory` to `cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory`.
2. Parameter `genomics.endpoint-url` value should be changed from `https://genomics.googleapis.com/` to 
`https://lifesciences.googleapis.com/`.
3. Also you should add a new mandatory parameter `genomics.location` to your backend configuration. Currently Google Cloud 
Life Sciences API is available only in `us-central1` and `europe-west2` locations.

### Alpha support for WDL optional outputs on PAPI v2

Cromwell 53 adds alpha-quality support for WDL optional outputs on PAPI v2 backends. Constructs such as: 

```
  struct MyStruct {
    String name
    File? file
  }
  .
  .
  .
  output {
    File? file_does_not_exist = "does_not_exist"
    Pair[String, File?] pair_file_does_not_exist = ("this", "does_not_exist")
    Map[String, File?] map_file_does_not_exist = { "does_not_exist": "does_not_exist" }
    Array[File?] array_file_does_not_exist = ["does_not_exist"]
    MyStruct struct_file_does_not_exist = object { name: "this", file: "does_not_exist" } 
  }
```

will not produce errors if the file `does_not_exist` does not exist. This support for optional files is considered alpha
quality for two reasons:

1. As seen in the example above, support for optional files extends to complex WDL types but there is a restriction that
all `File` components of non-primitive types must be optional. e.g. Cromwell would not allow the assignment of a 
missing file to the right side of a pair of type `Pair[File, File?]` since the left member of the pair is a non-optional
file. This restriction exists solely due to technical limitations in how type evaluation works in Cromwell today and
may be removed in a future Cromwell release.

2. Call caching does not work for calls with empty optional outputs. Cromwell currently does not recognize
that it is okay for optional output files to be missing, will incorrectly claim that any cache hits with missing 
optional output files are unusable, and will proceed to search for more cache hits which if found will also be unusable,
before eventually giving up and running the job. This behavior may be corrected in a future Cromwell release.

### Reference Disk Support

Cromwell 55 and later support mounting reference disks from prebuilt GCP disk images as an alternative to localizing large
input reference files on PAPI v2. Please note the configuration of reference disk manifests has changed starting with
Cromwell 57 and now uses the format documented below. 

Within the `config` stanza of a PAPI v2 backend the `reference-disk-localization-manifests`
key specifies an array of reference disk manifests:  

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        reference-disk-localization-manifests = [
          {
            "imageIdentifier" : "projects/broad-dsde-cromwell-dev/global/images/broad-references-disk-image",
            "diskSizeGb" : 500,
            "files" : [ {
              "path" : "gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.nhr",
              "crc32c" : 407769621
            }, {
              "path" : "gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.sa",
              "crc32c" : 1902048083
            },
            ...
          },
          ...
        ]
        ...
      }
    }
  }
}
```

Reference disk usage is an opt-in feature, so workflow submissions must specify this workflow option:

```json
{
  ...
  "use_reference_disks": true,
  ...
}
```

Using the first file in the manifest above as an example, assume a PAPI v2 backend is configured to use this manifest and the appropriate
`use_reference_disks` workflow option is set to `true` in the workflow submission. If a call in that workflow 
specifies the input `gs://my-references/enormous_reference.bam` and because that input matches the path of a file on the
reference image without the leading `gs://`, Cromwell would
arrange for a reference disk based on this image to be mounted and for the call's input to refer to the 
copy of the file on the reference disk, bypassing localization of the input.     

The Cromwell git repository includes a Java-based tool to facilitate the creation of manifests called
[CromwellRefdiskManifestCreatorApp](https://github.com/broadinstitute/cromwell/tree/develop/CromwellRefdiskManifestCreator).
Please see the help command of that tool for more details.

Alternatively for public data stored under `gs://gcp-public-data--broad-references` there exists a shell script to
extract reference data to a new disk and then convert that disk to a public image. For more information see
[create_images.sh](https://github.com/broadinstitute/cromwell/tree/develop/scripts/reference_disks/create_images.sh).

### Docker Image Cache Support

To optimize job execution time, Cromwell 55 and later support the use of Docker image caches on the PAPI v2 lifesciences beta backend. Docker image caches are not available on the PAPI v2 genomics alpha backend.
Configuration looks like:

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        docker-image-cache-manifest-file = "gs://path/to/a/docker/image/cache/manifest.json"
        ...
      }
    }
  }
}
```

Docker image cache manifest JSONs have a format like:

```json
{
  "biocontainers/samtools:1.3.1": "projects/broad-dsde-cromwell-dev/global/images/v1-docker-biocontainers-samtools-1-3-1",
  "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest": "projects/broad-dsde-cromwell-dev/global/images/v1-docker-gcr-io-gcp-runtimes-ubuntu-16-0-4-latest",
  ...
}
```

Docker image cache usage is an opt-in feature, so workflow submissions must specify this workflow option:

```json
{
  ...
  "use_docker_image_cache": true,
  ...
}
```

Individual tasks within a workflow can turn off Docker image caching through the use of a runtime attribute:

```wdl
task my_task {
  ...
  runtime {
    ...
    useDockerImageCache: false
  }
}
```

If Cromwell is running a workflow on PAPI v2 beta with Docker image caching enabled and a task specifies a
Docker image which corresponds to a configured Docker image cache JSON, Cromwell will arrange for the
job's VM to mount a disk built from the corresponding disk image. In the event that multiple
manifests describe disk images containing the specified Docker image, Cromwell will choose the disk image with the
smallest `diskSizeGb` value.

Conversely, Docker image caching can be turned off at the workflow level (either turned off explicitly or left at the
default setting of `false`) but turned on at the individual task level:

```wdl
task my_task {
  ...
  runtime {
    ...
    useDockerImageCache: true
  }
}
```

These settings could be useful for cost reasons: mounting Docker image caches adds nonzero cost
which might not be offset by eliminating Docker image pull times for long-running jobs.
# AWS Batch backend (beta)

Check out the [getting started guide](../tutorials/AwsBatch101.md) for the bulk of our documentation.

AWS support is fairly new to Cromwell and this reference section will expand as features are added and documented.

### Disks

Cromwell performs automatic disk sizing on your behalf when running with the AWS backend, so attributes like 
```
disks: "local-disk"
```
or
```
disks: "/some/mnt"
```
are adequate to specify a disk that will suffice to complete your task.

To facilitate the running of workflows originally authored for Pipelines API on Google Cloud Platform, Cromwell's AWS backend can also interpret attributes like
```
disks: "local-disk 20 SSD"
```
and
```
disks: "/some/mnt 20 SSD"
```
The size information and HDD/SSD have no effect on this backend and Cromwell simply drops them.
#Backends

A backend is a way to run the commands of your workflow. Cromwell allows for backends conforming to
the Cromwell backend specification to be plugged into the Cromwell engine. Additionally, backends are included with the
Cromwell distribution:

* **[Local](Local)**
* **[HPC](HPC)**, including **[Sun Grid Engine](SGE), [LSF](LSF), [HTCondor](HTcondor) & [SLURM](SLURM)** 
    * Run jobs as subprocesses or via a dispatcher.
    * Supports launching in Docker containers.
    * Use `bash`, `qsub`, and `bsub` to run scripts.
* **[Google Cloud](Google)** 
    * Launch jobs on Google Compute Engine through the Google Genomics Pipelines API.
* **[GA4GH TES](TES)** 
    * Launch jobs on servers that support the GA4GH Task Execution Schema (TES).
* **[Alibaba Cloud](BCS)** 
    * Launch jobs on Alibaba Cloud BatchCompute service.
* **[AWS Batch (beta)](AWS.md)**
    * Use Job Queues on AWS Batch

HPC backends are put under the same umbrella because they all use the same generic configuration that can be specialized to fit the need of a particular technology.

Backends are specified in the `backend.providers` configuration. Each backend has a configuration that looks like:

```hocon
BackendName {
  actor-factory = "FQN of BackendLifecycleActorFactory class"
  config {
    ...
  }
}
```

The structure within the `config` block will vary from one backend to another; it is the backend implementation's responsibility
to be able to interpret its configuration.

The providers section can contain multiple backends which will all be available to Cromwell.

## Backend Job Limits

All backends support limiting the number of concurrent jobs by specifying the following option in the backend's configuration
stanza:

```
backend {
  ...
  providers {
    BackendName {
      actor-factory = ...
      config {
        concurrent-job-limit = 5
```

## Backend Filesystems

Each backend will utilize a filesystem to store the directory structure and results of an executed workflow.
The backend/filesystem pairings are as follows:

* Local and HPC backend use the [Shared Local Filesystem](HPC/#filesystems).
* Google backend uses the [Google Cloud Storage Filesystem](Google/#google-cloud-storage-filesystem).
* Alibaba Cloud backend uses the OSS Storage FileSystem.

Additional filesystems capabilities can be added depending on the backend.
For instance, an HPC backend can be configured to work with files on Google Cloud Storage. See the [HPC documentation](HPC) for more details.
The following configuration can be used as a base to allow Cromwell to interact with a [SLURM](https://slurm.schedmd.com/) cluster and dispatch jobs to it:

```hocon
SLURM {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    runtime-attributes = """
    Int runtime_minutes = 600
    Int cpus = 2
    Int requested_memory_mb_per_core = 8000
    String queue = "short"
    """

    submit = """
        sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${runtime_minutes} -p ${queue} \
        ${"-c " + cpus} \
        --mem-per-cpu ${requested_memory_mb_per_core} \
        --wrap "/bin/bash ${script}"
    """
    kill = "scancel ${job_id}"
    check-alive = "squeue -j ${job_id}"
    job-id-regex = "Submitted batch job (\\d+).*"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
**Local Backend**

The local backend will simply launch a subprocess for each job invocation and wait for it to produce a Return Code file (rc file) which will contain the exit code of the job's command.
It is pre-enabled by default and there is no further configuration needed to start using it.

It uses the local filesystem on which Cromwell is running to store the workflow directory structure.

You can find the complete set of configurable settings with explanations in the 
[Cromwell Example Configuration File][cromwell-examples-conf],
along with backend provider examples in the [Example Providers Folder][cromwell-examples-folder].

The Local backend makes use of the same generic configuration as HPC backends. The same [filesystem considerations](HPC#filesystems) apply.

**Note to OSX users**: Docker on Mac restricts the directories that can be mounted. Only some directories are allowed by default.
If you try to mount a volume from a disallowed directory, jobs can fail in an odd manner. Before mounting a directory make sure it is in the list
of allowed directories. See the [Docker documentation](https://docs.docker.com/docker-for-mac/osxfs/#namespaces) for how to configure those directories.

[cromwell-examples-conf]: https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends
# AWS Batch Backend

AWS Batch is a set of batch management capabilities that dynamically provision the optimal quantity and type of compute resources (e.g., CPU or memory optimized instances) based on the volume and specific resource requirements of the batch jobs submitted.

This section provides details on how to configure the AWS Batch backend with Cromwell. For instructions on common configuration and deployment tutorial, see [Getting started with AWS Batch](https://cromwell.readthedocs.io/en/develop/tutorials/AwsBatch101/). 



## Resources and Runtime Attributes

Cromwell and AWS Batch recognizes number of runtime attributes, more information can be found in the [customize tasks](/RuntimeAttributes#recognized-runtime-attributes-and-backends) page.


## Running Cromwell on an EC2 instance

Cromwell can be run on an EC2 instance and submit jobs to AWS Batch, AWS provide [CloudFormation stacks and guides](https://docs.opendata.aws/genomics-workflows/orchestration/cromwell/cromwell-overview/) to building the correct IAM permissions.



## Scaling Requirements
For a Cromwell server that will run multiple workflows, or workflows with many steps (e.g. ones with large scatter steps), it is recommended to setup a database to store workflow metadata.  The application config file will expect a SQL database location. Follow [these instructions](https://docs.aws.amazon.com/AmazonRDS/latest/AuroraUserGuide/aurora-serverless.create.html) on how to create a serverless Amazon Aurora database. 

## Configuring Cromwell for AWS Batch

Within the `*.conf` file, you have a number of options to change the Cromwell's interaction with AWS Batch.

### Filesystems
> More information about filesystems can be found on the [Filesystems page](/filesystems/Filesystems/).
> 

Amazon's S3 storage is a supported filesystem in both the engine and backend, this means that S3 files can be referenced at a workflow level, and as input files, provided they are prefixed by `'s3://'`.

* filesystems
* filesystems.s3.auth
* filesystems.s3.caching.duplication-strategy

### Configuring Authentication

To allow Cromwell to talk to AWS, the `default` authentication scheme uses the [default authentication provider](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html) with the following AWS search paths:
- Environment Variables - `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
- Java system properties - `aws.accessKeyId` and `aws.secretKey`
- Default credential profiles file - Created by the AWS CLI, typically located at `~/.aws/credentials`
- _Instance profile credentials_ - Only relevant on EC2 instances

### Allowing private Docker containers

AWS Batch allows the use of private Docker containers by providing `dockerhub` credentials. Under the specific backend's configuration, you can provide the following object:

```hocon
(backend.providers.AWSBatch.config.)dockerhub = {
  // account = ""
  // token = ""
}
```

### More configuration options

* `(backend.providers.AWSBatch.config.)concurrent-job-limit` specifies the number of jobs that Cromwell will allow to be running in AWS at the same time. Tune this parameter based on how many nodes are in the compute environment.
* `(backend.providers.AWSBatch.config.)root` points to the S3 bucket where workflow outputs are stored. This becomes a path on the root instance, and by default is cromwell_root. This is monitored by preinstalled daemon that expands drive space on the host, ie AWS EBS autoscale.  This path is used as the 'local-disk' for containers.
**Sun GridEngine Backend**

The GridEngine and similar backends use programs such as `qsub` to launch a job and will poll the filesystem to determine if a job is completed.

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        # ... other configuration
      }
    }
  }
}
```

This backend makes the same assumption about the filesystem that the local backend does: the Cromwell process and the jobs both have read/write access to the CWD of the job.

The CWD will contain a `script.sh` file which will contain the same contents as the Local backend:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

The job is launched with a configurable command such as:

```bash
qsub \
    -terse \
    -V \
    -b n \
    -N ${job_name} \
    -wd ${cwd} \
    -o ${out}.qsub \
    -e ${err}.qsub \
    -pe smp ${cpu} \
    ${"-l m_mem_free=" + memory_gb + "gb"} \
    ${"-q " + sge_queue} \
    ${"-P " + sge_project} \
    ${script}
```

The SGE backend gets the job ID from parsing the `submit.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

The command used to submit the job is specified under the configuration key `backend.providers.SGE.config.submit`. It uses the same syntax as a command in WDL, and will be provided the variables:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `cwd` - The path where the script should be run.
* `out` - The path to the stdout.
* `err` - The path to the stderr.
* `job_name` - A unique name for the job.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out}.qsub \
            -e ${err}.qsub \
            ${script}
        """
      }
    }
  }
}
```

If the backend supports docker, the optional configuration keys `backend.providers.<backend>.config.submit-docker`
and  `backend.providers.<backend>.config.kill-docker` may be specified. When the WDL contains a docker runtime
attribute, this command will be provided three additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.
* `docker_cid` - The host path to which the [container ID file](https://docs.docker.com/engine/reference/run/#pid-equivalent) should be written.
* `docker_script` - The path of the `script` inside the docker container.
* `docker_out` - The path of the `out` inside the docker container.
* `docker_err` - The path of the `err` inside the docker container.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit-docker = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out}.qsub \
            -e ${err}.qsub \
            -l docker,docker_images="${docker}"
            -xdv ${cwd}:${docker_cwd}
            ${script}
        """
      }
    }
  }
}
```

If the backend would like to support additional runtime attributes they may be specified in the configuration key `backend.providers.<backend>.config.runtime-attributes`. It uses the same syntax as specifying runtime attributes in a task in WDL.

There are two special runtime attribute configurations, `cpu`, and `memory_<unit>`.

When the runtime attribute configuration `Int cpu` is specified, it is always validated as a positive integer.

When the runtime attribute configuration `Int memory_<unit>` or `Float memory_<unit>` is specified, it is provided to submit by the runtime attribute in WDL `memory`.

For example, if the backend specifies the configuration for `backend.providers.<backend>.config.runtime-attributes` as:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = "Float memory_mb"
      }
    }
  }
}
```

And the WDL specifies a task with:

```
task hello_gigabyte {
  command { echo "hello world" }
  runtime { memory: "1 GB" }
}
```

Then for this call, the backend will be provided an additional variable `memory_mb` set to `1000.0`.

Other runtime attributes may be defined by specifying them in under the runtime attributes configuration.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb
        String sge_project
        """
      }
    }
  }
}
```

These variables will then be passed from the WDL into the submit configuration. If one would like to have a default value, just like in WDL, the configuration may specify that the value have a default. The default must match the defined type or an error will be produced.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb = 2.0
        String sge_project = "default"
        """
      }
    }
  }
}
```

Optional values may also be used by appending `?` to the type:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float? memory_mb
        String? sge_project
        """
      }
    }
  }
}
```

The value will be passed to the submit configuration if provided, and omitted otherwise.

There are also configuration values related to how jobs are rechecked on startup and aborted.

The option is `backend.providers.<backend>.config.run-in-background`. When `true` the backend runs the submit configuration and records the unix process id (PID). To abort the job, the PID is stopped with the unix command `kill`. Upon a cromwell restart, the PID is checked via the unix command `ps` to see if it is still alive, before cromwell goes back to polling for the `rc` file.

When `backend.providers.<backend>.config.run-in-background` is `false`, the default, the backend must specify how read the job identifier from the stdout of the submit, how to kill the job, and how to check if the job is still running during a cromwell restart. These three configuration values are `job-id-regex`, `kill`, and `check-alive`, respectively:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        job-id-regex = "(\\d+)"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        """
      }
    }
  }
}
```

The `job-id-regex` should contain one capture group while matching against the whole line or stdout file. The `check-alive` should return zero if the job is still alive.

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
