# loinc2hpo Changelog



## current
* v 1.1.8

## v1.1.8
* Update to phenol 1.5.0

### loinc2hpo-core
* Adding enums for FHIR observation codes

## v0.0.2

### loinc2hpogui

* Added FHIR to HPO conversion dialog

### loinc2hpo-core

* added FHIR parsing

## v1.0.1
This is the first official release. Many features are added.

### loinc2hpogui
Multiple features are added. This will be the baseline for future tracking.

### loinc2hpo-core

* The core change is completely switching to FHIR to parse `observation` and `patient` resources. 

* Redesigned the annotation class. 
  - Use Code (system/namespace, code) to ensure uniqueness
  - An HPO term for a coded result is wrapped in a class that also indicate whether the term should be negated.
  - A complete annotation for a Loinc code contains many `Code` - `(HPO term, isNegated)` list. 
  
* Redesigned the logic process from observation to HPO term. 
  - The app first looks at whether the observation has an interpretation field. 
    - If it does, it will first try to find an annotation for the interpretation code directly; 
    - if it fails, it will try to convert convert the interpretation code to the internal code and then find the corresponding HPO term. 
  
  - If the app fails the last step, it will try to use the raw value and interpret it with the reference ranges.  
  
## v1.0.2

* Build Jar with all dependencies with maven-assembly-plugin

* Add META-INF/services to Core module because it appears that is what Jar requires

## v1.0.3

* Refactor pom files. 

* Refactor gitignore file. 

## v1.0.4

* Allow adding multiple labels to Github issues

* Disable the function to clear annotation fields during manual query

* Loinc entries change color if they have been annotated

* Add tooltips to HPO listview and treeview

* Allow user to switch to previously selected Loinc list

* Allow user to categorize Loinc entries

## v1.1.0 

* New develop version

Additional changes for this version

* Change menu `Edit` to `Configuration`
- [ ] update tutorial

* Create new features that allow user to manipulate a session
  
* Automatically retrieve information from auto-saved data for last session

## v1.1.1

* Session data now only saves terms for low, intermediate, and high value, instead for all 6 internal codes

* Basic and Advanced data are stored separately

## v1.1.2

* Session data now saves to a universal TSV file

* Restrict internal mappings

  Qn will not be mapped to "Presence" or "Absence" and Ord (of "Presence" type) will not be mapped to "high", "low", "normal"
  
* Internal codes changed match FHIR

  "system" is renamed to "FHIR";
  Code for "presence" changed from "P" to "POS", code for "not presence" changed from "NP" to "NEG" to be consistent with FHIR

* Show version in "About" message

* Refactor to remove deprecated classes

* Allow user to specify the path to hp OWL and Obo

## v1.1.3

* Add function to simulate patient data (patient and observation resources)

* Add function to allow uploading simulated data to hapi-fhir server

* Add function to allow downloading patient data from fhir server

## v1.1.4

* Add function to restart the app when necessary

* Prevent app from crashing when user's local HPO is outdated

* Added classes to appTempData patient phenopacket (phenotype only)
* Refactored to use phenol 1.0.0

## v1.1.5

* Add feature to copy annotation for one LOINC and paste to multiple selections of similar LOINC

* Add function to annotate LOINC panels

* Add function to convert FHIR messages for LOINC panels

* Add a LOINC list for tests with "unspecified specimen". Messages on those should be not converted to HPO. 

## v1.1.6

* Add feature to allow easy addition of LOINC lists: allow user to change colors of LOINC lists

* Added HPO parser for owl format.

* Refactored to use phenol 1.2.6

## v1.1.7

* Add algorithms to appTempData patients

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/709c959bb0024403a667affaf2b9f476)](https://www.codacy.com/app/peter.robinson/loinc2hpo?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=monarch-initiative/loinc2hpo&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/loinc2hpo/badge/?version=latest)](https://loinc2hpo.readthedocs.io/en/latest/?badge=latest)


# loinc2hpo
A Java library to map tests results from [LOINC](https://loinc.org/) codes to  
[Human Phenotype Ontology.](https://hpo.jax.org/app/) terms.
For details, please see [Zhang et al. (2012)](https://pubmed.ncbi.nlm.nih.gov/31119199/) Semantic integration of clinical laboratory tests from electronic 
health records for deep phenotyping and biomarker discovery. *NPJ Digit Med*. 2019;2:32.


## LOINC2HPO annotation
The LOINC to HPO mapping file is available at the
[loinc2hpoAnnotation](https://github.com/TheJacksonLaboratory/loinc2hpoAnnotation) repository.

## Using loinc2hpo
loinc2hpo is intended to be used as a software library. Selected functions are demonstrated in the CLI app.

## documentation
Please refer to http://loinc2hpo.readthedocs.io/en/latest/.

## spring framework app
We are developing a separate app that will specialize in one functionality of this app - converting FHIR observations into HPO terms. The new app will be coded with the Spring framework and we strive to achieve enterprise-level quality. Please refer to the app with the following link: https://github.com/OCTRI/fhir2hpo

## funding
We gratefully acknowledge funding by NCATS (CD2H project, A NATIONAL CENTER FOR DIGITAL HEALTH INFORMATICS INNOVATION), 1U24TR002306

============
FHIR Mapping
============

loinc2hpo can be used to convert LOINC-encoded laboratory results from FHIR.


Mapping LOINC to candidate HPO terms
====================================

LOINC observations have four main categories, quantitative(Qn), ordinal(Ord),
nominal (Nom) and narrative(Nar).
The majority of them are ``Qn`` (~50%) and ``Ord`` (~26%) and a smaller percentage
are ``Nom`` and ``Nar`` type (<10%). Currently, loinc2hpo maps QN, Ord, and Nom codes
but is not able to map Nar codes.

FHIR Interpretation codes
=========================

Test outcomes can represented with a code from the
`FHIR interpretation code valueset <https://www.hl7.org/fhir/valueset-observation-interpretation.html>`_.
HPO terms are qualitative. For instance,
`Thrombocytopenia HP:0001873 <https://hpo.jax.org/app/browse/term/HP:0001873>`_ does not
indicate whether there is a mild, moderate or severe degree of low blood platelets. Therefore
we map FHIR codes for different degrees of severity to the same internal code (these codes
are shown on the index page of this documentation).


Table 1: FHIR interpretation code set Mapping to internal code system

+------------------------------------+---------------------------+
|FHIR interpretation code value set  |internal code              |
+-------+----------------------------+---------------------------+
|Code   | Meaning                    |Code    | Meaning          |
+=======+============================+========+==================+
|<      |Off scale low               |L       |low               |
+-------+----------------------------+--------+------------------+
|>      |Off scale high              |H       |high              |
+-------+----------------------------+--------+------------------+
|A      |Abnormal                    |A       |abnormal          |
+-------+----------------------------+--------+------------------+
|AA     |Critically abnormal         |A       |abnormal          |
+-------+----------------------------+--------+------------------+
|AC     |Anti-complementary          |POS     |present           |
|       |substances present          |        |                  |
+-------+----------------------------+--------+------------------+
|B      |Better                      |N       |normal            |
+-------+----------------------------+--------+------------------+
|D      |Significant change down     |L       |low               |
+-------+----------------------------+--------+------------------+
|DET    |Detected                    |POS     |present           |
+-------+----------------------------+--------+------------------+
|H      |High                        |H       |high              |
+-------+----------------------------+--------+------------------+
|HH     |Critically high             |H       |high              |
+-------+----------------------------+--------+------------------+
|HM     |Hold for Medical Review     |U       |unknown           |
+-------+----------------------------+--------+------------------+
|HU     |Very high                   |H       |high              |
+-------+----------------------------+--------+------------------+
|I      |Intermediate                |N       |normal            |
+-------+----------------------------+--------+------------------+
|IE     |Insufficient evidence       |U       |unknown           |
+-------+----------------------------+--------+------------------+
|IND    |Indeterminate               |U       |unknown           |
+-------+----------------------------+--------+------------------+
|L      |Low                         |L       |low               |
+-------+----------------------------+--------+------------------+
|LL     |Critically low              |L       |low               |
+-------+----------------------------+--------+------------------+
|LU     |Very low                    |L       |low               |
+-------+----------------------------+--------+------------------+
|MS     |Moderately susceptible.     |U       |unknown           |
|       |(microbiology)              |        |                  |
+-------+----------------------------+--------+------------------+
|N      |Normal                      |N       |normal            |
+-------+----------------------------+--------+------------------+
|ND     |Not Detected                |NEG     |not present       |
+-------+----------------------------+--------+------------------+
|NEG    |Negative                    |NEG     |not present       |
+-------+----------------------------+--------+------------------+
|NR     |Non-reactive                |NEG     |not present       |
+-------+----------------------------+--------+------------------+
|NS     |Non-susceptible             |U       |unknown           |
+-------+----------------------------+--------+------------------+
|null   |No range defined or normal  |U       |unknown           |
|       |ranges don't apply          |        |                  |
+-------+----------------------------+--------+------------------+
|OBX    |Interpretation qualifiers   |U       |unknown           |
|       |in separate OBX segments    |        |                  |
+-------+----------------------------+--------+------------------+
|POS    |Positive                    |POS     |positive          |
+-------+----------------------------+--------+------------------+
|QCF    |Quality Control Failure     |U       |unknown           |
+-------+----------------------------+--------+------------------+
|R      |Resistant                   |U       |unknown           |
+-------+----------------------------+--------+------------------+
|RR     |Reactive                    |POS     |present           |
+-------+----------------------------+--------+------------------+
|S      |Susceptible                 |U       |unknown           |
+-------+----------------------------+--------+------------------+
|SDD    |Susceptible-dose dependent  |U       |unknown           |
+-------+----------------------------+--------+------------------+
|SYN-R  |Synergy - resistant	     |U       |unknown           |
+-------+----------------------------+--------+------------------+
|SYN-S  |Synergy - susceptible	     |U       |unknown           |
+-------+----------------------------+--------+------------------+
|TOX    |Cytotoxic substance present |POS     |present           |
+-------+----------------------------+--------+------------------+
|U      |Significant change up       |H       |high              |
+-------+----------------------------+--------+------------------+
|VS     |Very susceptible.           |U       |unknown           |
|       |(microbiology)              |        |                  |
+-------+----------------------------+--------+------------------+
|W      |Worse                       |A       |abnormal          |
+-------+----------------------------+--------+------------------+
|WR     |Weakly reactive             |POS     |present           |
+-------+----------------------------+--------+------------------+


The following graph summarizes the mapping.

  .. image:: images/annotation_scheme.png
    :width: 400
    :alt: Annotation scheme

Nominal outcomes
================

Nominal observations can use coding systems that are more difficult to handle.
For example, `Loinc 600-7 <https://loinc.org/600-7/>`_ (Bacteria identified in Blood by Culture)
may use a SNOMED concept to represent the
finding that *Staphylococcus aureus* is detected::

  "coding":[{
    "system": "http://snomed.info/sct",
    "code": "3092008",
    "display": "Staphylococcus aureus"
  }]

We currently map this to  `Bacteremia HP:0031864 <https://hpo.jax.org/app/browse/term/HP:0031864>`_,
but this term does not contain information about which bacterium was idenfitied in the blood.

We are currently extending the annotations to enable one to indicate what kind of bacteremia. In the above mentioned case,
one would use the NCBI `Taxonomy ID: 1280 <https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=1280>`_
for *Staphylococcus aureus*.






Getting started
===============

Overview
--------

loinc2hpo is a Java library. It requires at least Java 11. It is intended to be used by other applications
for specific purposes. This tutorial shows how to install loinc2hpo and how to use it in a typical Java program.

Installation
------------
First clone the library from GitHub. ::

    git clone https://github.com/monarch-initiative/loinc2hpo.git

Now use maven to install the library. ::

    cd loinc2hpo
    mvn install

We plan to place loinc2hpo on maven central in the future, which will make this step unnecessary.

Using loinc2hpo in maven projects
---------------------------------

To use the loinc2hpo in your own Java project, add the following to your pom file.

.. code-block:: XML

  <properties>
    <loinc2hpo.version>1.7.0</loinc2hpo.version>
  </properties>
  (...)
  <dependencies>
    <dependency>
      <groupId>org.monarchinitiative</groupId>
      <artifactId>loinc2hpo-core</artifactId>
      <version>${loinc2hpo.version}</version>
    </dependency>
    <dependency>
      <groupId>org.monarchinitiative</groupId>
      <artifactId>loinc2hpo-fhir</artifactId>
      <version>${loinc2hpo.version}</version>
    </dependency>
  </dependencies>

The ``loinc2hpo-fhir`` module is only required for working with FHIR of course.


Core module
~~~~~~~~~~~

The loinc2hpo annotation file referenced in the following code is available from the
`loinc2hpoAnnotation repository <https://github.com/TheJacksonLaboratory/loinc2hpoAnnotation>`_.
To run the code, ingest the annotation file and then pass



.. code-block:: java

  import org.monarchinitiative.loinc2hpocore.codesystems.ShortCode;
  import org.monarchinitiative.phenol.ontology.data.TermId;
  import org.monarchinitiative.loinc2hpocore.Loinc2Hpo;
  import org.monarchinitiative.loinc2hpocore.loinc.LoincId;
  import org.monarchinitiative.loinc2hpocore.codesystems.Outcome;

  String annot_file = "loinc2hpo-annotations.tsv";
  Loinc2Hpo loinc2Hpo = new Loinc2Hpo(annot_file);
  LoincId loincId = new LoincId("26515-7");
  Optional<Hpo2Outcome> opt = loinc2Hpo.query(loincId, Outcome.LOW());
  if (opt.isPresent()) {
    Hpo2Outcome hpo2outcome = opt.get();
    TermId hpoId = hpo2outcome.getHpoId();
    Outcome outcome = hpo2outcome.outcome();
    // do something with the HPO term and the outcome (low in this example).
  }


FHIR module
~~~~~~~~~~~

The FHIR module is intended to be used with `HAPI FHIR <https://hapifhir.io/>`_.
It can be used with the FHIR specifications DSTU3, R4, or R5.

.. code-block:: java

    import org.monarchinitiative.loinc2hpofhir.Loinc2HpoFhir;
    import org.hl7.fhir.r5.model.*;
    import org.monarchinitiative.loinc2hpocore.codesystems.Outcome;
    import org.monarchinitiative.loinc2hpocore.loinc.LoincId;

    String annot_file = "loinc2hpo-annotations.tsv";
    Loinc2HpoFhir loinc2hpoFHIR = new Loinc2HpoFhir(String path);
    // The following is a R5 Observation
    Observation observation = getObservationFromSomewhere(); // your code does this
    Optional<Hpo2Outcome> opt = loinc2hpoFHIR.r5(observation);
    if (opt.isPresent()) {
      Hpo2Outcome hpo2outcome = opt.get();
      TermId hpoId = hpo2outcome.getHpoId();
      Outcome outcome = hpo2outcome.outcome();
      // do something with the HPO term and the outcome.
    }


The ``Loinc2HpoFhir`` has analogous methods called ``dstu3`` and ``r4`` for the other
FHIR versions.Introduction to LOINC
=====================

`LOINC <https://loinc.org/>`_ (Logical Observation Identifiers Names and Codes)
provides a set of universal names and ID codes for identifying laboratory and clinical
test results.

LOINC currently provides ~92,000 entries that define the names and IDs of laboratory tests.
The following shows three examples of LOINC codes:

  .. image:: images/loinc_examples.png

Each LOINC entry represents a laboratory test.

Parts of LOINC entry
--------------------

- ``LOINC``: unique identifier
- ``Name``: structured term name
- ``Component``: defines the analyte in the test
- ``Property``: defines "kinds of quantities". It can be divided into five several categories, mass, substance, catalytic activity, and number or counts. Each category is further divided into subclasses, for example MCnc or "mass concentration" is a subclass of "mass" property, while ``NCnc`` or "number of concentration (count/vol)" and ``Naric`` or "number aeric (number per area)" are subclasses of counts.
- ``Time``: defines whether a measurement was made at a moment, or aggregated from a series of physiologic states. The three examples are all ``PT`` ("points"), meaning that they are measurements at a single time point. As an example, a test on "daily urine amount" will be labeled as ``24H`` ("24 hours").
- ``Aspect``: defines a modifier for a measurement over a duration. For example, "8H^max heart rate" means the "max" heart rate measured during an 8-hour period. Min, max, first, last, mean typically appear here.
- ``System``: can be considered as the specimen for the test, such as "serum", "blood", "urine", "cerebrospinal fluid" etc.
- ``Scale``" defines the scale of the measurement. Scale is the most important information for our application. The following table summarizes possible values of ``scale``.
- ``Method``: defines the method used for the measurement.


Table 1: LOINC Scale Types

+----------------+------+-------------------------------------------------------------------------------------+
| Scale Type     | Abbr.| Description                                                                         |
+================+======+=====================================================================================+
| Quantitative   | Qn   | The result of the test is a numeric value that relates to a continuous numeric      |
|                |      | scale.                                                                              |
+----------------+------+-------------------------------------------------------------------------------------+
| Ordinal        | Ord  | Ordered categorical responses, e.g., positive, negative;                            |
+----------------+------+-------------------------------------------------------------------------------------+
| Quantitative   | OrdQn| Test can be reported as either Ord or Qn.                                           |
+----------------+------+-------------------------------------------------------------------------------------+
| Nominal        | Nom  | Nominal or categorical responses that do not have a natural ordering.               |
+----------------+------+-------------------------------------------------------------------------------------+
| Narrative      | Nar  | Text narrative.                                                                     |
+----------------+------+-------------------------------------------------------------------------------------+
| “Multi”        | Multi| Many separate results structured as one text “glob”                                 |
+----------------+------+-------------------------------------------------------------------------------------+
| Document       | Doc  | A document that could be in many formats (XML, narrative, etc.)                     |
+----------------+------+-------------------------------------------------------------------------------------+
| Set            | Set  | Used for clinical attachments                                                       |
+----------------+------+-------------------------------------------------------------------------------------+



``Qn``, ``Ord`` and ``Nom`` are the three most frequently used LOINC codes,
accounting for about 99% of data. ``Qn`` typically makes up about 80% of cases.














Introduction to FHIR
====================

Fast Healthcare Interoperability Resources (FHIR) is a set of standards created by HL7 that aims to improve the exchange and interoperability of healthcare information. Consider the following case:

  .. image::images/fhir_intro.png

The patient visits clinics, hospitals and receive medical tests in clinical labs, or visits drug stores. It would be nice if the physicians can retrieve his test results from other providers, or the drug store can automatically receive the prescriptions for the patient. The above is not possible without a standardization of information exchange. This is how FHIR come into play.

FHIR is also extremely helpful in reusing Electronic Health Record for scientific research. EHRs provide a wealth of medical information, but their lack of standardization presents a huge hurdle for their reusability. With FHIR, one can easily process patient record stored in different EHR systems.


Resources
---------

FHIR uses RESTful API technology for data exchange. It takes all medical information as individual pieces of resources. For example, ``Administration``-related resources are mainly about identifications (e.g. about patient, about practioner, about device ...). Diagnostics resources are those about observations, specimen etc. There are also resources about medications and health insurances. Resources are typically stored and exchanged in JSon format. There are unique identifiers in each resource that can cross reference to other another so that it is possible to find, say, all resources related to one patient.

  .. image:: images/FHIR_resource1.png
  .. image:: images/FHIR_resource2.png
  .. image:: images/FHIR_resource3.png

``Observation``

Among all the resources, ``Observation`` is the most central one in EHR. ``Observation`` describes a medical test result for one patient. The following JSon snippet is an exerpt from a real world (ref: `hl7.org <https://www.hl7.org/fhir/observation-example-f001-glucose.json.html>`_).

Information about the ID of the resource. ::

  {
  "resourceType": "Observation",
  "id": "f001",
  "text": {
    "status": "generated",
    "div": "<div xmlns=\"http:\/\/www.w3.org\/1999\/xhtml\"><p><b>Generated Narrative with Details<\/b><\/p><p><b>id<\/b>: f001<\/p><p><b>identifier<\/b>: 6323 (OFFICIAL)<\/p><p><b>status<\/b>: final<\/p><p><b>code<\/b>: Glucose [Moles\/volume] in Blood <span>(Details : {LOINC code '15074-8' = 'Glucose [Moles\/volume] in Blood', given as 'Glucose [Moles\/volume] in Blood'})<\/span><\/p><p><b>subject<\/b>: <a>P. van de Heuvel<\/a><\/p><p><b>effective<\/b>: 02\/04\/2013 9:30:10 AM --&gt; (ongoing)<\/p><p><b>issued<\/b>: 03\/04\/2013 3:30:10 PM<\/p><p><b>performer<\/b>: <a>A. Langeveld<\/a><\/p><p><b>value<\/b>: 6.3 mmol\/l<span> (Details: UCUM code mmol\/L = 'mmol\/L')<\/span><\/p><p><b>interpretation<\/b>: High <span>(Details : {http:\/\/hl7.org\/fhir\/v2\/0078 code 'H' = 'High', given as 'High'})<\/span><\/p><h3>ReferenceRanges<\/h3><table><tr><td>-<\/td><td><b>Low<\/b><\/td><td><b>High<\/b><\/td><\/tr><tr><td>*<\/td><td>3.1 mmol\/l<span> (Details: UCUM code mmol\/L = 'mmol\/L')<\/span><\/td><td>6.2 mmol\/l<span> (Details: UCUM code mmol\/L = 'mmol\/L')<\/span><\/td><\/tr><\/table><\/div>"
  },
  "identifier": [
    {
      "use": "official",
      "system": "http:\/\/www.bmc.nl\/zorgportal\/identifiers\/observations",
      "value": "6323"
    }
  ],
  "status": "final",

Information about the test, identified with a LOINC code.  ::

  "code": {
    "coding": [
      {
        "system": "http:\/\/loinc.org",
        "code": "15074-8",
        "display": "Glucose [Moles\/volume] in Blood"
      }
    ]
  },

Information about the subject/patient of this resource. ::

  "subject": {
    "reference": "Patient\/f001",
    "display": "P. van de Heuvel"
  },

Some meta data about this test: ::

  "effectivePeriod": {
    "start": "2013-04-02T09:30:10+01:00"
  },
  "issued": "2013-04-03T15:30:10+01:00",
  "performer": [
    {
      "reference": "Practitioner\/f005",
      "display": "A. Langeveld"
    }
  ],

Measured value and reference range: ::

  "valueQuantity": {
    "value": 6.3,
    "unit": "mmol\/l",
    "system": "http:\/\/unitsofmeasure.org",
    "code": "mmol\/L"
  },

  "referenceRange": [
    {
      "low": {
        "value": 3.1,
        "unit": "mmol\/l",
        "system": "http:\/\/unitsofmeasure.org",
        "code": "mmol\/L"
      },
      "high": {
        "value": 6.2,
        "unit": "mmol\/l",
        "system": "http:\/\/unitsofmeasure.org",
        "code": "mmol\/L"
      }
    }
  ]

Interpretation from physicians: ::

  "interpretation": {
    "coding": [
      {
        "system": "http:\/\/hl7.org\/fhir\/v2\/0078",
        "code": "H",
        "display": "High"
      }
    ]
  },
  }

``Patient``

``Patient`` manages all relevant information about the patient, such as name, address, sex, etc. The following snippet is an example in JSon (ref: `hl7.org <https://www.hl7.org/fhir/patient-example-f001-pieter.json>`_). The code is pretty self-explanatory.

::

  {
  "resourceType": "Patient",
  "id": "f001",
  "text": {
    "status": "generated",
    "div": "<div xmlns=\"http:\/\/www.w3.org\/1999\/xhtml\"><p><b>Generated Narrative with Details<\/b><\/p><p><b>id<\/b>: f001<\/p><p><b>identifier<\/b>: 738472983 (USUAL), ?? (USUAL)<\/p><p><b>active<\/b>: true<\/p><p><b>name<\/b>: Pieter van de Heuvel <\/p><p><b>telecom<\/b>: ph: 0648352638(MOBILE), p.heuvel@gmail.com(HOME)<\/p><p><b>gender<\/b>: male<\/p><p><b>birthDate<\/b>: 17\/11\/1944<\/p><p><b>deceased<\/b>: false<\/p><p><b>address<\/b>: Van Egmondkade 23 Amsterdam 1024 RJ NLD (HOME)<\/p><p><b>maritalStatus<\/b>: Getrouwd <span>(Details : {http:\/\/hl7.org\/fhir\/v3\/MaritalStatus code 'M' = 'Married', given as 'Married'})<\/span><\/p><p><b>multipleBirth<\/b>: true<\/p><h3>Contacts<\/h3><table><tr><td>-<\/td><td><b>Relationship<\/b><\/td><td><b>Name<\/b><\/td><td><b>Telecom<\/b><\/td><\/tr><tr><td>*<\/td><td>Emergency Contact <span>(Details : {http:\/\/hl7.org\/fhir\/v2\/0131 code 'C' = 'Emergency Contact)<\/span><\/td><td>Sarah Abels <\/td><td>ph: 0690383372(MOBILE)<\/td><\/tr><\/table><h3>Communications<\/h3><table><tr><td>-<\/td><td><b>Language<\/b><\/td><td><b>Preferred<\/b><\/td><\/tr><tr><td>*<\/td><td>Nederlands <span>(Details : {urn:ietf:bcp:47 code 'nl' = 'Dutch', given as 'Dutch'})<\/span><\/td><td>true<\/td><\/tr><\/table><p><b>managingOrganization<\/b>: <a>Burgers University Medical Centre<\/a><\/p><\/div>"
  },
  "identifier": [
    {
      "use": "usual",
      "system": "urn:oid:2.16.840.1.113883.2.4.6.3",
      "value": "738472983"
    },
    {
      "use": "usual",
      "system": "urn:oid:2.16.840.1.113883.2.4.6.3"
    }
  ],
  "active": true,
  "name": [
    {
      "use": "usual",
      "family": "van de Heuvel",
      "given": [
        "Pieter"
      ],
      "suffix": [
        "MSc"
      ]
    }
  ],
  "telecom": [
    {
      "system": "phone",
      "value": "0648352638",
      "use": "mobile"
    },
    {
      "system": "email",
      "value": "p.heuvel@gmail.com",
      "use": "home"
    }
  ],
  "gender": "male",
  "birthDate": "1944-11-17",
  "deceasedBoolean": false,
  "address": [
    {
      "use": "home",
      "line": [
        "Van Egmondkade 23"
      ],
      "city": "Amsterdam",
      "postalCode": "1024 RJ",
      "country": "NLD"
    }
  ],
  "maritalStatus": {
    "coding": [
      {
        "system": "http:\/\/hl7.org\/fhir\/v3\/MaritalStatus",
        "code": "M",
        "display": "Married"
      }
    ],
    "text": "Getrouwd"
  },
  "multipleBirthBoolean": true,
  "contact": [
    {
      "relationship": [
        {
          "coding": [
            {
              "system": "http:\/\/hl7.org\/fhir\/v2\/0131",
              "code": "C"
            }
          ]
        }
      ],
      "name": {
        "use": "usual",
        "family": "Abels",
        "given": [
          "Sarah"
        ]
      },
      "telecom": [
        {
          "system": "phone",
          "value": "0690383372",
          "use": "mobile"
        }
      ]
    }
  ],
  "communication": [
    {
      "language": {
        "coding": [
          {
            "system": "urn:ietf:bcp:47",
            "code": "nl",
            "display": "Dutch"
          }
        ],
        "text": "Nederlands"
      },
      "preferred": true
    }
  ],
  "managingOrganization": {
    "reference": "Organization\/f001",
    "display": "Burgers University Medical Centre"
  }
  }

``Patient`` is essential in cases when the interpretation of an observation, e.g. height or weight, dependents on the sex, age or other relevant information of the patient. In this case, one can retrieve ``Patient`` with the link in the "subject" field of ``Observation``.
loinc2hpo Documentation
=======================

Loinc2hpo is a Java library designed to convert the findings of laboratory tests to HPO codes.
For instance, if the test is `LOINC 26515-7 Platelets [#/volume] in Blood <https://loinc.org/26515-7/>`_
and the outcome of the test is an abnormally low value, then we can infer the
`Human Phenotype Ontology (HPO) <https://hpo.jax.org/app/>`_ term
`Thrombocytopenia HP:0001873 <https://hpo.jax.org/app/browse/term/HP:0001873>`_.

The goal of this library is to encode EHR (Electronic Health Record) laboratory
data using HPO terms to extend the kinds of analysis that can be performed.
Laboratory results can be leveraged as phenotypic features for analysis.

  .. image:: images/mission.png
     :align: center
     :scale: 60 %

The library currently has three modules.

loinc2hpo-core
==============

This library contains the core functionality. It imports the annotation file from
the `loinc2hpoAnnotation <https://github.com/TheJacksonLaboratory/loinc2hpoAnnotation>`_ repository
(loinc2hpo-annotations.tsv), and for any combination of LOINC Id (laboratory test) and
outcome, it finds the appropriate HPO term if one exists. We use set of internal codes
to represent lab outcomes.

+---------+------------------------------------------------------------------------------+
| Code    | Explanation                                                                  |
+=========+==============================================================================+
| L       | Low (below normal range). Used for quantitative tests (Qn).                  |
+---------+------------------------------------------------------------------------------+
| H       | High (above normal range). Used for quantitative tests (Qn).                 |
+---------+------------------------------------------------------------------------------+
| N       | Normal (within normal range). Used for quantitative tests (Qn).              |
+---------+------------------------------------------------------------------------------+
| NEG     | Negative (not present, a normal result). Used for ordinal tests (Ord)        |
+---------+------------------------------------------------------------------------------+
| POS     | Positive (present, an abnormal result). Used for ordinal tests (Ord)         |
+---------+------------------------------------------------------------------------------+
| NOM     | Nominal (an abnormal result). Used for nominal tests (Nom)                   |
+---------+------------------------------------------------------------------------------+

loinc2hpo-fhir
==============

This library provides an interface that extracts LOINC-encoded data from
`FHIR - Fast Healthcare Interoperability Resources <hl7.org/fhir>`_ data. Specifically,
it provides an interface that takes a FHIR `Observation <https://www.hl7.org/fhir/observation.html>`_
and attempts to extract a LOINC code and an outcome; if successful, these are passed to the
core loinc2hpo module to get the corresponding HPO term.

loinc2hpo-cli
=============

This is a command-line interface tool that can be used to obtain descriptive statistics or
perform quality control of the input files.

Contents
========

.. toctree::
   :maxdepth: 1

   intro_to_LOINC
   intro_to_FHIR
   FHIR_mapping.rst
   getting_started

GitHub repo
-----------
The source code of Loinc2hpo can be found at GitHub:
https://github.com/monarch-initiative/loinc2hpo

Contact
-------

Peter Robinson
peter.robinson@jax.org

`The Jackson Laboratory <https://www.jax.org/>`_
10 Discovery Drive
Farmington, CT
USA

Xingmin Aaron Zhang
kingmanzhang@gmail.com

Curation
--------

We have developed a JavaFX application to curate loinc2hpo data. This app is not needed to
use the library, but may be of interest to potential contributors: https://github.com/pnrobinson/loinc2hpoMiner.
