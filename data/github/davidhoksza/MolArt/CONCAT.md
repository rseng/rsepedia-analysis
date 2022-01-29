# Changelog

The project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

##  [Unreleased]

### Added

### Removed

### Fixed


##  [1.12.0]

### Added

- Ability to specify PyMOL export file name

### Removed

### Fixed

- Fixed unpSeqNumber passed in HighlightCustomElements
- Fixed not showing structure viewer error messages

##  [1.11.1]

### Added

### Removed

### Fixed

- Fixed unpSeqNumber passed in HighlightCustomElements


##  [1.11.0]

### Added

- Added `variantFilterDefaultBehavior` making ProtVista's |variant filter dialog to behave
as standard checkboxes. By default, when consequence clicked in ProtVista and that consequence
is already checked, all other checkboxes become unchecked which is not intuitive for some users.

### Removed

### Fixed

- Mapping in HighlightCustomElements when the structure starts before UniProt sequence (e.g. 2I4I)

##  [1.10.0] - 2021-10-27

### Added

- AlphaFold export to PyMOL
- Ability to show additional information about hovered element.
- Choose which structure to load by default.

### Removed

- Removed min-width of structure and sequence view

### Fixed

- User-defined visuals do not create default visuals (causing problems once removed)
- Failing when AlphaFold or SMR not available

##  [1.9.1] - 2021-10-15

### Added

### Removed

### Fixed

- User visuals stay colored after feature selection


##  [1.9.0] - 2021-10-15

### Added

- Controlling user visuals via MolArt API (add visual, remove visual)

### Removed

### Fixed


##  [1.8.0] - 2021-09-29

### Added

- AlphaFold DB structures available.
- Ability to specify default structure transparency level using lmInitSurfaceTransparency.
- Download all annotation data in CSV.
- Highlithing a region in the sequence view also highlights corresponding part of the structure in the
structure view. 

### Removed

### Fixed

- Taking into account inserted segments in structure with respect to sequence (e.g.6I53:A).
- Structure IDs show up as uppercase. 



##  [1.7.0] - 2020-10-22

### Added

- Get information about which regions in sequence have a structure mapped.
- Emitting an event every time a structure has been loaded.

### Removed

### Fixed

- Possible conflict of jQuery with Semantics UI jQuery object
- Returned "main" field in package.json

##  [1.6.0] - 2020-10-11

### Added

- Integrated new version of ProtVista which correctly handles changes in Proteins API 
- Ability to focus at a residue in the structure view.
- Convenience function to find out whether MolArt is loaded
- Ability to provide custom 3D structure
- Ability to provide custom sequence and sequence-structure mapping
- PredictProtein is passed sequence when UniProt ID is not available
- Possibility to pass also chain specification when restricting the list of structures
- Event listeners for sequence and structure mouse hover events
- Ability to control highliting of position in both sequence and structure views outside of MolArt
- Input options sanitization (e.g. Uniprot ID can be passed with trailing spaces)
- Destroy method of the MolArt object not only destroys LiteMol object, but also removes everything from the root container.
- SMR records sorted by coverage

### Removed

### Fixed

- SwissProt structures showing twice in ProtVista when PDBe mapping not availble 
- Beginning of region when shifting observed ranges
- Included updated ProtVista which handles problems with shifted ruler
- Issue with shifted arrow icon in categories headers in some scenarios (removed verctial bottom alignment)
- MouseMove emmiter events generated numbers outside of the range of sequence numbers.
- Dehighlight in structure view when mouse moved outside of sequence view.
- Using local jQuery for container object to avoid clashing with plugins which claim global jQuery instance
- Improved error message when LiteMol cannot process the structure file

##  [1.5.0] - 2019-09-30

### Added

- Ability to show unobserved structure regions in the sequence view

### Removed

- Showing TaxId in the structure description in the sequence view

### Fixed

- Keeping highlighted first or last residue in the structure when mouse moves past the end of sequence 
in the sequence view
- When clicking an overlay arrow icon in the sequence view and a feature was selected, i.e. the yellow bar was 
acive, the bar did not get hidden

## [1.4.0] - 2019-08-22

### Added

- PredictProtein as a default category 
- Possibility to specify different categories and features labels for PyMol export

### Fixed

- Fix of issue when a user selection is empty or non-empty but does not intersect with current structure
- Fixed issue with features being out of categories in PyMol export when special characters were present
- Fixing of PyMOL export issue when beginning or end of selection falls into an unobserved region
- Export to PyMOL was giving wrong selections if sequence and structure numbering did not match.
- Dependencies update

## [1.3.3] - 2019-05-09

### Fixed

- Fixed issue with positioning of the components when ProteinAPI takes too long to load.

## [1.3.2] - 2019-04-08

### Added

- Export to PyMOL restricted to current chain

### Fixed

- Fixed moving MolArt logo when scrolling 

## [1.3.1] - 2019-01-26

### Fixed

- Fix PyMol export selections when sequence and structure numbering does not match


## [1.3.0] - 2019-01-13

### Added
- PyMOL export now contains CA selection.
- PyMOL export now has hierarchical categories the same way the users sees them in sequence viewer.
- List of mapped structures in the sequence view can be sorted by name (by default it is sorted by coverage).
This can be set in the constructor.
- User defined highlights. The user can define custom selections to be shown on the structure. Can be useful, 
e.g. for showing pathogenic mutations on every structure. The format follows LiteMol's selection and 
 highlight "language".

### Fixed
- Fixed issue with logo showing out of plugin boundaries in some cases.
- Handled situation when category and feature type has the same name when exporting to PyMOL (e.g. ANTIGEN).


# MolArt (MOLeculAR structure annoTator)

MolArt is a responsive, easy-to-use JavaScript plugin which enables users to
view annotated protein sequence (including variation data from large scale
studies) and overlay the annotations over a corresponding experimental
or predicted protein structure.
It couples protein sequence annotation capabilities provided by [ProtVista](https://github.com/ebi-uniprot/ProtVista)
 (or more precisely its [modified responsive version](https://github.com/davidhoksza/protvista) 
 implemented when developing MolArt) with structure visualization 
 capabilities provided by [LiteMol](https://github.com/dsehnal/LiteMol). 
 Since it does not have any software dependencies and all the data are obtained on the fly,
  it is easy to integrate it to any web page.

Examples of MolArt's use can be found at https://cusbg.github.io/MolArt/.

[Documentation](https://github.com/davidhoksza/MolArt/tree/master/docs) shows usage of the plugin.

The plugin is being developed at the Luxembourg Center for Systems Biomedicine, University of Luxembourg.

To cite MolArt use: **MolArt: A molecular structure annotation and visualization tool** (doi: [10.1093/bioinformatics/bty489](https://doi.org/10.1093/bioinformatics/bty489))


<div style="text-align:center;">
    <img src="gitweb/teaser.png" style="width:80%"/>
</div>


## Features overview

- Visualization of protein structure as provided by LiteMol
- Annotation of protein sequence as provided by ProtVista
- Annotations of protein sequence from [PredictProtein](https://en.wikipedia.org/wiki/Predictprotein) 
- Mapping of structure on corresponding substring in the sequence
- Automatic retrieval of sequence data based on UniProt ID and corresponding experimental structures from PDB
- Retrieval of predicted models from [AlphaFold Protein Structure DB](https://alphafold.ebi.ac.uk/) (ADB) and [SWISS-MODEL Repository](https://swissmodel.expasy.org/repository) (SMR) if 
no PDB structure is available
- Controlling transparency of the structure to see both cartoon and surface view of the structure
- Hovering over position in sequence to highlight it in structure and vice versa
- Color overlay any sequence feature over the structure
- Color overlay all sequence features of given type over the structure
- Color overlay individual variation over the structure
- Color overlay all mutations to given amino acid over the structure
- Color overlay mutation frequency of residues over the structure
- Exports of the structure and annotations to [PyMol](https://pymol.org/2/) for advanced inspection 
(the export does not include variants)
- Upload of custom annotations
- Ability to control which structures will be shown
- Ability to select and highlight specific residues (not necessarily corresponding to annotations) and specific atoms .
The selection can be visualized either as a surface, a balls and stick representation or simply van der 
Waals-based spheres with given color and transparency.
- Ability to provide custom sequence and mapping



## Data sources

- Sequence and annotation data
  - Sequence information comes from [UniProt website REST API](https://www.uniprot.org/help/api) or can be provided by the user
  - Sequence annotations are provided by 
    - ProtVista plugin which utilizes the EBI's 
  [Proteins REST API](https://www.ebi.ac.uk/proteins/api/doc/). 
  Proteins API includes access to variation, proteomics and antigen services containing 
  "annotations imported and mapped from large scale data sources, such as 1000 Genomes, 
  ExAC (Exome Aggregation Consortium), COSMIC (Catalogue Of Somatic Mutations In Cancer), 
  PeptideAtlas, MaxQB (MaxQuant DataBase), EPD (Encyclopedia of Proteome Dynamics) and HPA, 
  along with UniProtKB annotations for these feature types".
    - [PredictProtein API](https://en.wikipedia.org/wiki/Predictprotein) which gives access to a range
    of integrated sequence-based [prediction methods](https://www.predictprotein.org/about) such as
    conservation or prediction of functional changes. If the user prefers only experimental data, 
    the PredictProtein annotations can be easily turned off with the exclusion parameter as 
    described in the [documentation](https://github.com/davidhoksza/MolArt/tree/master/docs) 
  
- Structure mapping
    - Automatic
        - To obtain the mapping between UniProt and PDB, MolArt is using the [SIFTS API](https://www.ebi.ac.uk/pdbe/api/doc/sifts.html), part of the [PDBe REST API](http://www.ebi.ac.uk/pdbe/pdbe-rest-api).
        - In case the SIFTS mapping yields no PDB structures, ADB and SMR are queried using their respective APIs for available models.
    - Custom
        - When working with sequence which is not in UniProt or when the mapping is not known, 
        sequence-structure mapping can be provided. This includes
            - Molecule-level mapping, i.e. which PDB structures correspond to the sequence
            - Residue-level mapping, i.e. which regions in sequence correspond to which regions in the structure

- Structure data
  - In case an experimental structure is available in PDB for given UniProt ID, this structure is downloaded by LiteMol. In this case, MolArt instructs LiteMol to use the mmCIF format.
  - In case there is no experimental structure in PDB, but a model exists in ADB 
  or SMR, MolArt instructs LiteMol to use the mmCIF-format in case of ADB and PDB-format in case of SMR.

## How to use MolArt

- Obtain the JavaScript file with MolArt and link it from your web page
- Create a container DIV (or SPAN) element which will hold the viewer
- Create a JavaScript object and pass it reference to the DIV

The detail description of how to incorporate MolArt into your project can be found in the [developer documentation](https://github.com/davidhoksza/MolArt/tree/master/docs).

### Examples of use
Exmpales of how to use MolArt are located in the ``examples`` folder.
- See the ``bare.html`` file for a bare bones example of how to use the plugin.
- The ``plugin-page.html`` is slightly styled plugin usage.
- For an advanced example, see the ``web`` directory. It contains a simple web application which enables querying Uniprot (only top 10 matches are retrieved) and for every found record one can click the UniProt ID which creates a new tab with new instance of MolArt for that UniProt ID.

## Contributing

We would be happy to hear about your use cases, experiences and ideas/feature requests. Either raise an issue [here](https://github.com/davidhoksza/MolArt/issues) or get in touch by mail (david.hoksza at gmail).

## Support

Please submit your issues through the MolArt's repository issue tracker available [here](https://github.com/davidhoksza/MolArt/issues).

## License

This project is licensed under the Apache 2.0 license, quoted below.

Copyright (c) 2018 David Hoksza

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Support

<p align="center">
  <img src="img/logo-elixir.png" />
</p>

PrankWeb is a part of services provided by ELIXIR – European research infrastructure for biological information.
See [services](https://www.elixir-czech.cz/services) provided [ELIXIR's Czech Republic Node](https://www.elixir-czech.cz/).

# Developers documentation

## Obtaining MOLART

#### Using precompiled distribution file

The easiest way is to simply download the precompiled distribution file
``molart.js`` from the [dist](https://github.com/davidhoksza/MolArt/tree/master/dist) directory.

#### Building MOLART from source code

Other option, especially if you want to do some modification to the code before using it, is to build MOLART
directly from the source codes.

###### Obtain the source

Download the Github repository with the following command:

```
git clone https://github.com/davidhoksza/MolArt
```

###### Building the source codes

MOLART uses [Gulp](https://gulpjs.com/) as its build system. So if Gulp is not installed on your system yet, you will first need to download it and only then run its default task to build MOLART.

```
npm install -g gulp
npm install
gulp
```

This process will result in a single file in the ``dist`` directory with single file
named ``molart.js``, the only file needed to use the plugin.

#### Using MOLART as a NPM package

If your application uses NPM as the packaging system, MOLART can be
also easily obtain with NPM using:

```
npm install git://github.com/davidhoksza/MolArt
```

## Using the plugin

#### Include the code into a web page

All javascript files, style sheets, SVG and font files used by ProtVista and LiteMol are included in the MOLART distribution js file, so that is the only thing you need to embed into your web page. You can embed it into your HTML using the script tag:

```
<script type="text/javascript" src="molart.js"></script>
```

or, if you are using Browserify or Webpack to bundle your Javascript the above script tag is of course not required and the require` statement discussed bellow is sufficient.

#### Create a container

Create a DIV or SPAN element.

```
<div id="pluginContainer"></div>
```

If this is the only element in the page then the plugin will span the whole width and height of the window. However, you can limit the width, height or position of the plugin which will then take the provided space and resize accordingly. For example, you can use the following styles to center the plugin in the middle of the window and make it 50% width and 80% height of the window:

```
<div id="pluginContainer" style="position: absolute;  width: 50%; height: 80vh; left:25%; top:10%"></div>
```

However, please bear in mind that there exists a minimal width (800px) the plugin needs to comfortably accommodate all its components. So if you set the width of the container below this threshold, horizontal scrollbar will appear.

#### Create instance of MOLART

Finally, you need to create an instance of the plugin and specify a UniProt ID and the ID of the container.

```javascript
var MolArt = require('MolArt');
molart = new MolArt({
    uniprotId: 'P37840',
    containerId: 'pluginContainer'
});
```


If you are using a bundling system to build your application, not using NPM and not using the `script` tag to embed MolArt, you need to change the `require` location to point to the plugin script. The above syntax should work just fine if you are using NPM.

#### Destroy MOLART

If you need to create an instance repeatedly , e.g. every time a user clicks on UniprotID you open/create a tab and spawn a new instance of MOLART and when she closes the tab, you remove all its content. In such situations LiteMol (or rather THREE.js) does not release all the WebGL related structures and still tries to draw something into HTML elements which are not available any more. To get rid of the warning messages, you can call the `destroy` method. However, if there are still some callback functions active, which try to access LiteMol, you will get the

## Options, parameters and events

All parameters for ProtVista are also available in MOLART, which simply
takes them and passes them to ProtVista. These include ability to exclude
cateogires, customization of their order and also ability to specify
custom categories/data sources. The full description with examples can be
found in [ProtVista's documentation](http://ebi-uniprot.github.io/ProtVista/developerGuide.html#starting-protvista).

#### Define visibility and order of categories

In order to exclude categories or customize their order, simply pass additional
parameters to the `MolArt` object as you would when using ProtVista only. The list 
of all categories available in ProtVista can be found [here](https://github.com/davidhoksza/ProtVista/blob/90000af6e11131d138faab050d89e30f27e03e19/src/config.json#L2).

Additionally, MolArt adds [PredictProtein](https://en.wikipedia.org/wiki/Predictprotein) annotations category with ID 
`PREDICT_PROTEIN` and this category is treated the same way as the ProtVista annotations. Therefore, it can be 
turned off using the ``exclusion`` parameter. 

```javascript
molart = new MolArt({
    uniprotId: 'P37840',
    containerId: 'pluginContainer',
    categoryOrder: ['PTM'],
    exclusions: ['PREDICT_PROTEIN', 'SEQUENCE_INFORMATION', 'STRUCTURAL', 'TOPOLOGY', 'MUTAGENESIS', 'MOLECULE_PROCESSING']
});
```

Since the core of MOLART lies in mapping between the sequence and structure,
the customization does not impact the experimental/predicted structures category
which always appears first.

#### Custom data sources

You can also provide custom annotations which will be automatically
mapped over the structures. The format of the annotation data
is defined in ProtVista documentation in [Adding your own sources](http://ebi-uniprot.github.io/ProtVista/developerGuide.html#adding-your-own-sources) section.
However, unlike in ProtVista, in MOLART you can pass the data directly
in the constructor and the annotations can thus be
generated on the fly, if needed. Moreover, MOLART lets you to
use multiple data sources.

The following example shows how to mix "local"
with "external" data sources. In the following example, three categories
are created two of which consists of randomly generated data while the third
downloads data from http://localhost/externalFeatures_P37840.json (if available).

```javascript
function initializeTestDataSet(sequence, catName){

    const ix1 = Math.floor(Math.random() * sequence.length);
    const ix2 = ix1 + Math.floor(Math.random() * (sequence.length - ix1));

    return {
        sequence: sequence,
        features: [
            {
                type: "ACT_SITE",
                category: catName,
                begin: String(ix1),
                end: String(ix1),
                color: "#00F5B8"
            },
            {
                type: "MY_REGION",
                category: catName,
                begin: String(ix1),
                end: String(ix2),
                color: "#FF7094"
            }
        ]
    };
}

const sequence = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA';
const customDataSources = [
    {
        source: 'RANDOM',
        useExtension: false,
        data: initializeTestDataSet(sequence, 'MY_CATEGORY1')
    },
    {
        source: 'RANDOM',
        useExtension: false,
        data: initializeTestDataSet(sequence, 'MY_CATEGORY2')
    }
    ,
    {
        source: 'RANDOM',
        useExtension: true,
        url: 'http://localhost/externalFeatures_'
    }
];

const MolArt = require('MolArt');
molart = new MolArt({
    uniprotId: 'P37840',
    containerId: 'pluginContainer',
    customDataSources: customDataSources
});
```

#### Custom sequence and sequence-structure mapping

MolArt is able to handle situations when a user wants to visualize 
- a sequence which is not available in UniProt 
- PDB structures where the sequence-structure mapping is not available in Protein API
- structures which are not in PDB

###### User-provided sequence-structure mapping

To provide a custom sequence-structure mapping, one needs to pass the MolArt constructor the information about
which structures map to given sequence and which regions in the sequence map to which regions in the structure. This
is basically the information which MolArt automatically gets from the `https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/` 
and `https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/` PDBe REST API endpoints. Also the format is 
somewhat similar to what these endpoints provide.

The format also allows to define custom source of data which can be either a valid URI (`structure.uri`) or simply
a string with the structure data (`structure.data`). In any case, the user needs to pass in information about the
format of the data (`structure.format`) which can be either `mmCIF` or `PDB`.

```javascript

const MolArt = require('MolArt');
molart = new MolArt({
    uniprotId: 'P29274',
    containerId: 'pluginContainer',
    sequenceStructureMapping: [
        {
          id: '5uig',
          chainId: 'A',
          structure: {
            format: 'mmCIF', //valid parameters are PDB and mmCIF
            data: window.uig_a, //the structure in the PDB or mmCIF format
            uri: 'https://www.ebi.ac.uk/pdbe/static/entry/5uig_updated.cif' //data and uri parameters are mutually exclusive
          },
          start: 27, // where the structure begins with respect to the full molecule sequence
          end: 438, // where the structure ends with respect to the full molecule sequence
          seqStart: 1, // where the structure begins with respect to the sequence (the sequence does not have to covert the full molecule, therefore seqStart does not have to be 1)
          seqEnd: 316, // where the structure ends with respect to the sequence
          coverage: [
              {
                  start: {
                      residue_number: 30, // position of the region with respect to the full molecule sequence
                      author_residue_number: 4, // position with respect to the beginning of the structure (in this case first three residues are not observed, i.e. residues 27, 28, 29 with respect to the full molecule)
                      author_insertion_code: undefined,
                  },
                  end: {
                      residue_number: 148,
                      author_residue_number: 174,
                      author_insertion_code: undefined,        
                  }
              },
              {
                  start: {
                      author_residue_number: 159,
                      author_insertion_code: undefined,
                      residue_number: 185
                  },
                  end: {
                      author_residue_number: 1048,
                      author_insertion_code: undefined,
                      residue_number: 282
                  }
              },
              {
                  start: {
                      author_residue_number: 1056,
                      author_insertion_code: undefined,
                      residue_number: 290
                  },
                  end: {
                      author_residue_number: 311,
                      author_insertion_code: undefined,
                      residue_number: 433
                  }
              }
            ]
        }
    ]
});
```

###### User-provided sequence

To provide custom sequence, one simply needs to do so in the MolArt constructor. Obviously, in such a case
the sequence must be accompanied by a sequence-structure mapping as there is not UniProtID which could be 
used to obtain the mapping. Moreover, it is not possible to provide both sequence and UniProt ID at the
same time.

```javascript
const molstar = new MolArt({
      containerId: 'pluginContainer',
      sequence:
          'MPIMGSSVYITVELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAI\n' +
          'PFAITISTGFCAACHGCLFIACFVLVLTQSSIFSLLAIAIDRYIAIRIPLRYNGLVTGTR\n' +
          'AKGIIAICWVLSFAIGLTPMLGWNNCGQPKEGKNHSQGCGEGQVACLFEDVVPMNYMVYF\n' +
          'NFFACVLVPLLLMLGVYLRIFLAARRQLKQMESQPLPGERARSTLQKEVHAAKSLAIIVG\n' +
          'LFALCWLPLHIINCFTFFCPDCSHAPLWLMYLAIVLSHTNSVVNPFIYAYRIREFRQTFR\n' +
          'KIIRSHVLRQQEPFKAAGTSARVLAAHGSDGEQVSLRLNGHPPGVWANGSAPHPERRPNG\n' +
          'YALGLVSGGSAQESQGNTGLPDVELLSHELKGVCPEPPGLDDPLAQDGAGVS',
      sequenceStructureMapping: [
          {
              id: '5uig',
              chainId: 'A',
              structure: {
                    format: 'mmCIF', //valid parameters are PDB and mmCIF
                    data: window.uig_a, //the structure in the PDB or mmCIF format
                    uri: 'https://www.ebi.ac.uk/pdbe/static/entry/5uig_updated.cif' //data and uri parameters are mutually exclusive
              },
              start: 27, // where the structure begins with respect to the full molecule sequence
              end: 438, // where the structure ends with respect to the full molecule sequence
              seqStart: 1, // where the structure begins with respect to the sequence (the sequence does not have to covert the full molecule, therefore seqStart does not have to be 1)
              seqEnd: 316, // where the structure ends with respect to the sequence
              coverage: [
                  {
                      start: {
                          residue_number: 30, // position of the region with respect to the full molecule sequence
                          author_residue_number: 4, // position with respect to the beginning of the structure (in this case first three residues are not observed, i.e. residues 27, 28, 29 with respect to the full molecule)
                          author_insertion_code: undefined,
                      },
                      end: {
                          residue_number: 148,
                          author_residue_number: 174,
                          author_insertion_code: undefined,
                      }
                  },
                  {
                      start: {
                          author_residue_number: 159,
                          author_insertion_code: undefined,
                          residue_number: 185
                      },
                      end: {
                          author_residue_number: 1048,
                          author_insertion_code: undefined,
                          residue_number: 282
                      }
                  },
                  {
                      start: {
                          author_residue_number: 1056,
                          author_insertion_code: undefined,
                          residue_number: 290
                      },
                      end: {
                          author_residue_number: 311,
                          author_insertion_code: undefined,
                          residue_number: 433
                      }
                  }
              ]
          }
      ]    
});
```

#### Other options

- ```sortStructures``` - when set to ```id``` the lists of experimental and predicted
structures are sorted by their name. This changes the default behavior when the lists
are sorted by coverage, i.e. how much of the sequence is covered by a structure. 

- ```categoriesTooltips``` - array of arrays of size two containing category code and
its tooltip text. This option allows to set user-defined titles for categories. An example can be found in the plugin page example (```categoriesTooltips: [['DOMAINS_AND_SITES', 'Describes various domains and sites'], ['PTM', 'Post-translational modifications']]```).
The tooltips for tracks can be set using the ProtVista [custom config file](http://ebi-uniprot.github.io/ProtVista/developerGuide.html#further-customization) option.
Since all options passed to MolArt's constructor are further passed to ProtVista,
this option is available in MolArt as well.

- ```enableTooltips``` (default ```true```) - when set to ```false``` hovering over a category
does not show a tooltip for that category as ProtVista does by default. This option
was added because default category tooltips in ProtVista are not very illustrative.

- ```highlightByHovering``` (default ```false```) - when set to ```true```
hovering over an annotation in the sequence view highlights the corresponding part
of the structure in the structure view. By default, an annotation is highlighted
only when clicked on in the sequence view.

- ```alwaysLoadPredicted``` (default ```false```) - when set to ```true```
MolArt always connects to ADB (AlphaFold DB) and SMR (SwissProt Model Repository) and downloads available models; by default, 
it will query ADB and SMR only when no experimental structure is available.

- ```pdbIds``` (default ```undefined```) - list of PDB IDs and possibly chain IDs (such as ```['1ryp:b', '4r17']```) 
which are
supposed to be shown in the Experimental structures category. If given PDB ID is presented without a chain, 
all chains of given protein will be listed.
Structures in the mapping outside of this list will not be shown. If
not set or the list is empty, no restriction takes place.

- ```smrIds``` (default ```undefined```) - list of SMR IDs (such as ```['3q27.1']```) which are
supposed to be shown in the Predicted structures category. 
Structures in the mapping outside of this list will not be shown. If
not set or the list is empty, no restriction takes place.

- ```defaultStructureId``` - ID of the structure to be displayed when MolArt starts. This can be either
PDB ID, SMR ID or AlphaFold ID. For example, if we want to load by default the AlphaFold structure,
we can call MolArt with:

```javascript
new MolArt({
    uniprotId: 'O00571',        
    containerId: 'pluginContainer',
    alwaysLoadPredicted: true,
    defaultStructureId: 'AF-O00571-F1',    
})
```

- ```lmInitSurfaceTransparency``` - specifies the default transparency level (0-100) of the surface 
representation of the displayed structure.

- ```extraHighlights``` (default ```undefined```) - allows to highlight a list of residues and
even restrict atoms of the residues. Moreover, one can specify the type of highlight. Specifically, 
on needs to pass an object containing 
an array where each elements defines a selection and how that selection should be viusalized.
The LiteMol documentation specifies available parameters for the [balls and sticks](https://webchemdev.ncbr.muni.cz/LiteMol/SourceDocs/interfaces/litemol.bootstrap.visualization.molecule.ballsandsticksparams.html)
and [surface](https://webchemdev.ncbr.muni.cz/LiteMol/SourceDocs/interfaces/litemol.bootstrap.visualization.molecule.surfaceparams.html)
visualizations. If the object has controlVisibility set to `true`, a dropdown will be shown in the
header where the user can turn on and of the defined highlights. 

    ```
    extraHighlights: {
        controlVisibility: false, //whether the list of custom highlights will be shown as a dropdown through which the can control visibility of the individual highlights
        label: 'Extra selections',
        content: [
            {
                label: 'Extra 1 - CA, CE',
                showOnStart: true,
                sequenceNumbers: [58, 50],
                atomNames: ['CA', 'CE'],
                visual: {
                    type: 'BallsAndSticks',
                    //https://webchemdev.ncbr.muni.cz/LiteMol/SourceDocs/interfaces/litemol.bootstrap.visualization.molecule.ballsandsticksparams.html
                    params: { useVDW: true, vdwScaling: 1, bondRadius: 0.13, detail: 'Automatic' },
                    color: {r:1, g: 0, b: 0},
                    alpha: 1
                }

            }, {
                label: 'Extra 2',
                showOnStart: false,
                sequenceNumbers: [60, 61, 62],
                //atomNames: ['C'],
                visual: {
                    type: 'Surface',
                    //https://webchemdev.ncbr.muni.cz/LiteMol/SourceDocs/interfaces/litemol.bootstrap.visualization.molecule.surfaceparams.html
                    params: {isWireframe: true, probeRadius: 0, density: 1.25, smoothing: 3},
                    color: {r: 0, g: 0, b: 1},
                    alpha: 0.8
                }
            }
        ]
    }
    
    ```
- ```labelCallback``` - A function which is called when a position in the structured is hovered. It should
return a string which will be displayed in the structure view together with the default information about the 
hovered atom. This can be useful in situations where you, for example, highlight (e.g. by a sphere)
residues which you know are often mutated and then hovering over that highlighted residue you can show additional
information about the mutation in the structure view. The function is passed one argument containing information about the hovered atom and corresponding
residue. See the actual fields by running: 
```javascript
    new MolArt({
        uniprotId: 'P14867',        
        containerId: 'pluginContainer',
        labelCallback: (info) => {console.log(info); return info.resInfo.unpSeqNumber}
});
```

 ### Events
 
 MolArt object generates several events. The app which owns the MolArt object can register listeners 
 for those events and thus react on what happens inside the plugin.
 
 * ``lmReady``
    Emitted every time when any of the ProtVista's track is loaded. This means that once
    *lmReady* is emitted, the sequence view is ready to accept user events.
 * ``pvReady``
    Emitted once LiteMol (the structure view) loads a structure. This happens also when
    a user selects a structure in the sequence view. 
 * ``structureLoaded``
    Emitted every time the active structure changes.
 * ``structureMouseOn``
    * ```molart.on("structureMouseOn", residue => console.log(residue))```
        * Triggers when a residue is hovered over in the structure view.
         The listener is passed an argument with a ``residue`` object containing information about 
         the corresponding residue. The ``residue`` object
        is a dictionary having the available information as returned from LiteMol:
             * ```javascript
               authName: "SER"
               authSeqNumber: 13
               chain:
                   asymId: "A"
                   authAsymId: "A"
                   entity:
                       entityId: "1"
                       index: 0
                       __proto__: Object
                   index: 0
                   __proto__: Object
               index: 12
               insCode: null
               isHet: 0
               name: "SER"
               seqNumber: 281
               __proto__: Object 
               ```
               
 * ``structureMouseOff``
     * Triggers when the mouse stops being over a residue.
 * ``sequenceMouseOn``
    * Triggers when a mouse hovers over a position in the sequence view. 
    The listener is passed an argument which corresponds to the actual position in the sequence. 
 * ``sequenceMouseOff``
    * Triggers when the mouse leaves the space of the sequence view. It is also triggered
    when it leaves a category DIV.
* ``pyMOLExportFileName``
    * Allows user to change the PyMOL export filename:
    ```javascript
        pyMOLExportFileName: {
            PDB: 'PDB-{id}_{chain}.py',
            SMR: 'SMR-{id}_{chain}.py',
            AF: 'AF-{id}_{chain}.py',
            USER: 'user-{id}_{chain}.py'
    }
    ```
 
 ### MolArt API

 Molart can be partially controlled from by accessing methods available in the MolArt object.

 * ``getPvController``
    * Returns the ProtVista (sequence view) controller.

 * ``getLmController``
     * Returns the LiteMol (structure view) controller.

 * ``getGlobals``
     * Returns the globals object which includes various objects (including the Lm and Pv controllers).

 * ``highlightInSequence(pos)``
     * Highlights the specified residue in the sequence view.

 * ``deHighlightInSequence()``
     * Removes highlighting in the sequence view.

 * ``highlightInStructure(pos)``
     * Highlights the specified residue in the structure view. The position is the position in sequence (which might differ from the position in the structure).

 * ``deHighlightInStructure()``
     * Removes highlighting in the structure view.
    
 * ``focusInStructure(seqPos, neighborhoodSize)``
     * Focuses at a given residue. If *neighborhoodSize* provided, all residues having atoms
     in a sphere of given radius (in Angstroms) will be dispalayed as well.
    
 * ``structureLoaded()``
     * True if a structure is loaded. Can be false either if there does not exist anystructure
     for the molecule or when the stucture it not loaded yet.
   
 * ``getLmController().setVisualInteractive(visualId, selectionDef, visualDev)``
   * Enables users to add custom highlights into the visualization. For example, one can draw red balls 
   on the position of C-alpha atoms for residues which are have some property. Another option si to
   "draw" a mesh over a range of residues. The method has three parameters:
     * ``visualId`` identifies the visual, so that it can be later removed
     * ``selectionDef`` specifies residues and possibly also atoms over which the visual definition should
     be applied, e.g. ``{sequenceNumbers: [63,64], atomNames: ['CA', 'CB']}`` (``atomNames`` is optional)
     * ``visualDef`` specifies the type of visual to be applied over the selection, e.g.
     ``{type: 'Surface', params: {isWireframe: true, probeRadius: 0, density: 1.25, smoothing: 3}, color: {r: 0, g: 0, b: 1}, alpha: 0.8 }``.
     The visual definition is passed directly to LiteMol and thus any definition valid in LiteMol
     is also valid for MolArt. Check the [interactive_highlight.html](../examples/interactive_highlight.html)  
     for some examples.
     
 * ``getLmController().clearVisualInteractive(visualId)``
   * Removes previously created custom highlight.

 * ``getSeqStrRange()``
     * Get ranges in the sequence which have a structure mapped. The return value
     is an array of arrays of size 2 as there can be undetermined parts of the structure 
     resulting in multiple observed regions. 
     If there is not structure set, the method throws an error (can be checked in advance
     using the ``structureLoaded`` method). 